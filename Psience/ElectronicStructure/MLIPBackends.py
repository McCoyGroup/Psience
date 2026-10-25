"""
Machine-learned interatomic potential backends: AIMNet2, HIPNN, and anything exposed as an ASE
calculator (MACE, MACE-POLAR, UMA, UPET, Orb)
"""

import contextlib
import functools
import os
import warnings

import numpy as np

import McUtils.Numputils as nput
import McUtils.Devutils as dev

from .Backends import (
    ElectronicStructureBackend, quiet_mode, resolve_torch_device, hf_authenticate, environment_override
)
from .ChargeResponse import FluctuatingChargeResponse

__all__ = [
    "AIMNet2Backend",
    "HIPNNBackend",
    "ASECalculatorBackend",
    "MACEBackend",
    "MACEPolarBackend",
    "UMABackend",
    "UPETBackend",
    "OrbBackend",
]


def _stack_structures(results):
    """concatenate per-structure `[value, d1, ...]` lists along the leading structure axis"""
    return [np.concatenate(terms, axis=0) for terms in zip(*results)]


@ElectronicStructureBackend.register('aimnet2')
class AIMNet2Backend(ElectronicStructureBackend):
    """
    AIMNet2 (and AIMNet2-NSE) energies, charges, and fluctuating-charge polarizabilities, with
    derivatives by autograd. Borrows structure from `AIMNet2ASE` to call the model directly.
    """
    distance_units = 'Angstroms'
    property_units = {
        'energy': 'ElectronVolts',
        'charges': 'ElementaryCharge',
        'polarizability': None,  # model-native (charges in e, distances in Angstrom, energies in eV)
    }
    analytic_derivative_orders = {
        'energy': 2,
        'charges': 1,
        'polarizability': 1,
    }
    property_shapes = {
        'charges': lambda n: (n,),
        'polarizability': (3, 3),
    }

    def __init__(self, atoms, model='aimnet2', charge=None, device=None, model_dir=None, multiplicity=None, quiet=True):
        """
        :param model: an AIMNet2 model name (`'aimnet2'`, `'aimnet2nse'`, ...) or a loaded `AIMNet2Calculator`
        :param device: the torch device (default: CUDA if available)
        :param model_dir: a cache directory for downloaded weights (default `$MODEL_CACHE_DIR`)
        """
        super().__init__(atoms, charge=charge, multiplicity=multiplicity, quiet=quiet)
        self.model_name = model
        self.model = self.setup_model(model, device=device, model_dir=model_dir)

    @classmethod
    def handle_specialization(cls, tag):
        return {'model': tag}

    # region model setup

    @classmethod
    @contextlib.contextmanager
    def _overload_aimnet_modeldir(cls, root_dir):
        if root_dir is None:
            root_dir = os.environ.get("MODEL_CACHE_DIR")

        if root_dir is not None:
            import aimnet.calculators.model_registry as reg
            cur = reg._maybe_download_asset
            try:
                def _maybe_download_asset(file: str, url: str, sha256: str = None, **etc) -> str:
                    # recent `aimnet` versions also pass the expected `sha256`
                    import requests
                    filename = os.path.join(root_dir, 'aimnet', file)
                    if not os.path.exists(filename):
                        os.makedirs(os.path.dirname(filename), exist_ok=True)
                        print(f"Downloading {url} -> {filename}")
                        with open(filename, "wb") as f:
                            response = requests.get(url, timeout=60)
                            f.write(response.content)
                    if sha256 is not None:
                        import hashlib
                        with open(filename, 'rb') as f:
                            digest = hashlib.sha256(f.read()).hexdigest()
                        if digest != sha256:
                            raise ValueError(f"checksum mismatch for {filename}: expected {sha256}, got {digest}")
                    return filename
                reg._maybe_download_asset = _maybe_download_asset
                yield _maybe_download_asset
            finally:
                reg._maybe_download_asset = cur
        else:
            yield None

    @staticmethod
    def _patch_aimnet_neighbor_list():
        """Retry native neighbor construction when reported counts exceed storage.

        Some AIMNet/nvalchemi versions return the full neighbor counts without
        raising an overflow error when the adaptive buffer is too small. Install
        this check once, retaining the native implementation and its arguments.
        Legacy aimnet2calc installations do not need this patch.
        """
        try:
            from aimnet.calculators import calculator
        except ImportError:
            return False

        neighbor_type = getattr(calculator, "AdaptiveNeighborList", None)
        if neighbor_type is None:
            return False
        original = neighbor_type.__call__
        if getattr(original, "__psience_neighbor_completeness__", False):
            return True

        @functools.wraps(original)
        def checked_neighbors(self, *args, **kwargs):
            while True:
                allocated = self.max_neighbors
                nbmat, counts, shifts = original(self, *args, **kwargs)
                required = int(counts.max().item()) if counts.numel() else 0
                if required <= nbmat.shape[1]:
                    return nbmat, counts, shifts
                if required <= allocated:
                    raise RuntimeError(
                        "AIMNet returned an incomplete neighbor matrix despite "
                        "sufficient allocation: "
                        f"{required} neighbors, {nbmat.shape[1]} columns, "
                        f"{allocated} allocated"
                    )
                self.max_neighbors = self._round_to_16(
                    max(required, int(self.max_neighbors * 1.5))
                )

        checked_neighbors.__psience_neighbor_completeness__ = True
        neighbor_type.__call__ = checked_neighbors
        return True

    @classmethod
    def setup_model(cls, model, device=None, model_dir=None):
        """
        Load an `AIMNet2Calculator` (or return `model` if it already is one).
        """
        cls._patch_aimnet_neighbor_list()
        if isinstance(model, str):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                with dev.OutputRedirect():
                    try:
                        from aimnet.calculators import AIMNet2Calculator
                        modeldir_override = cls._overload_aimnet_modeldir(model_dir)
                    except ImportError:
                        from aimnet2calc import AIMNet2Calculator
                        modeldir_override = contextlib.nullcontext()

            # the weights are downloaded when the calculator is built, so the cache-directory
            # override has to be active here (not just around the import)
            with quiet_mode(), modeldir_override:
                device = resolve_torch_device(device)

                try:
                    calc = AIMNet2Calculator(model, device=device)
                except TypeError:
                    calc = AIMNet2Calculator(model)  # can't set device in older aimnet

        else:
            calc = model

        return calc

    # endregion

    # region autodiff

    @staticmethod
    def autodiff(expr, coord, pad_dim, base_shape=None, create_graph=False, retain_graph=False):
        import torch
        # here forces have shape (N, 3) and coord has shape (N+1, 3)
        # return hessian with shape (N, 3, N, 3)
        coord_ndim = len(coord.shape)
        if base_shape is not None:
            coord = coord.reshape(base_shape + (-1, 3))
        shape = coord.shape + expr.shape
        grads = []
        npad = len(pad_dim)
        rav_shape = tuple(
            np.prod([expr.shape[coord_ndim*i+k] for k in range(coord_ndim)], dtype=int)
            for i in range((expr.ndim - npad) // coord_ndim)
        ) + pad_dim
        if len(rav_shape) > 0:
            cache = {}
            inds = np.array(np.unravel_index(np.arange(np.prod(expr.shape, dtype=int)), rav_shape)).T
            if npad == 0:
                shinds = np.sort(inds, axis=1)
            else:
                shinds = np.concatenate([np.sort(inds[:, :-(npad)], axis=1), inds[:, -npad:]], axis=-1)
            for i,_f in enumerate(expr.flatten()):
                ravioli = tuple(inds[i])
                key = tuple(shinds[i])
                if ravioli == key:
                    g = torch.autograd.grad(_f,
                                            coord,
                                            retain_graph=True,
                                            create_graph=create_graph)
                    g = g[0]
                    grads.append(g)
                    cache[key] = g
                else:
                    grads.append(cache[key])
        else:
            grads = []
            if len(expr.shape) == 0:
                grads.append(
                    torch.autograd.grad(expr,
                                        coord,
                                        retain_graph=True,
                                        create_graph=create_graph)[0]
                )
            else:
                for n, _f in enumerate(expr):
                    g = torch.autograd.grad(_f,
                                            coord,
                                            retain_graph=True,
                                            create_graph=create_graph)[0]
                    grads.append(g)

        # Autograd stacks output components first. Move coordinate axes in
        # front of them; reshaping alone interleaves vector-property Jacobians.
        deriv = torch.stack(grads).reshape(expr.shape + coord.shape)
        axes = tuple(range(expr.ndim, expr.ndim + coord.ndim)) + tuple(range(expr.ndim))
        deriv = deriv.permute(axes)
        shape_tuple = ((slice(None, -1),) + (slice(None),) * (coord_ndim-1)) * ((len(shape) - npad) // coord_ndim)
        return deriv, shape_tuple

    def prep_eval(self, coords, order, charge=None, multiplicity=None, keep_graph=None, forces=None, hessian=None, **opts):
        """
        Run the model on `coords` (`(..., natoms, 3)`, Angstrom), setting up autograd for
        `order` energy derivatives.

        :return: `(base_shape, flat_coords, data)`
        """
        import torch

        charge, multiplicity = self.resolve_state(charge, multiplicity)
        base_shape = coords.shape[:-2]
        coords = coords.reshape((-1,) + coords.shape[-2:])

        if forces is None:
            forces = order > 0
        if hessian is None:
            hessian = order > 1

        arg_dict = {
            'coord': torch.tensor(coords.reshape(-1, 3), dtype=torch.float32, device=self.model.device),
            'numbers': torch.tensor(
                np.repeat(np.array([self.numbers]), coords.shape[0], axis=0).flatten(),
                dtype=torch.int64, device=self.model.device
            ),
            'charge': torch.tensor(charge, dtype=torch.float32, device=self.model.device),
            'cell': None,
            'mol_idx': torch.tensor(
                np.repeat(np.arange(coords.shape[0]), coords.shape[1]).flatten(),
                dtype=torch.int64,
                device=self.model.device
            )
        }

        if multiplicity is not None:
            arg_dict['mult'] = torch.tensor(multiplicity, dtype=torch.float32, device=self.model.device)
        elif 'nse' in str(self.model_name).lower():
            # NSE (open-shell, two-channel) models need an explicit spin key; don't guess one
            raise ValueError(
                f"AIMNet2 model {self.model_name!r} is a spin-aware (NSE) model but no spin multiplicity was supplied; "
                "set `multiplicity` on the molecule (`mol.spin`) or pass `multiplicity=` explicitly"
            )

        data = self.model.prepare_input(arg_dict)
        if hessian and data['mol_idx'][-1] > 0:
            raise NotImplementedError('higher-derivative calculation is not supported for multiple molecules')
        data = self.model.set_grad_tensors(data,
                                           forces=keep_graph or forces,
                                           stress=False,
                                           hessian=keep_graph or hessian
                                           )
        if keep_graph and not (forces or hessian):
            data['coord'].requires_grad_(True)
            self.model._saved_for_grad['coord'] = data['coord']
        with torch.jit.optimized_execution(False):
            data = self.model.model(data)
        # Recent AIMNet2 models leave Coulomb and DFTD3 corrections in
        # calculator-owned modules rather than the Torch model itself.
        # Match AIMNet2Calculator.eval before taking energy derivatives.
        if hasattr(self.model, '_run_external_modules'):
            data, external_terms = self.model._run_external_modules(
                data, forces=forces or hessian, stress=False, hessian=hessian
            )
            opts['coulomb_terms'] = external_terms
        data = self.model.get_derivatives(
            data,
            forces=forces,
            hessian=hessian,
            **opts
        )

        return base_shape, coords, data

    def _apply_autodiff(self, data, root_key, deriv_key, order, property_shape,
                        diff_var=None,
                        key_tag='deriv',
                        pad_coord=1,
                        base_shape=None,
                        create_graph=None,
                        retain_graph=None):
        if pad_coord != 1:
            raise NotImplementedError("all AIMNet coordinates are padded by 1, need to relax this condition")
        keys = []
        if order > 0:
            if diff_var is None:
                diff_var = self.model._saved_for_grad['coord']
            if deriv_key is not None:
                keys.append(f'{key_tag}_1')
                data[f'{key_tag}_1'] = (-data[deriv_key], (slice(None), slice(None)))
                for o in range(2, order + 1):
                    keys.append(f'{key_tag}_{o}')
                    data[f'{key_tag}_{o}'] = self.autodiff(data[f'{key_tag}_{o - 1}'][0],
                                                           diff_var,
                                                           property_shape,
                                                           base_shape=base_shape,
                                                           create_graph=create_graph if create_graph is not None else (
                                                                       (o == 1) or (o < order)),
                                                           retain_graph=retain_graph if retain_graph is not None else (
                                                                       o < order)
                                                           )
            else:
                if 'forces' in data: data['forces'].unbind()
                data[f'{key_tag}_0'] = (data[root_key], (slice(None),))
                for o in range(1, order + 1):
                    keys.append(f'{key_tag}_{o}')
                    data[f'{key_tag}_{o}'] = self.autodiff(data[f'{key_tag}_{o - 1}'][0],
                                                           diff_var,
                                                           property_shape,
                                                           base_shape=base_shape,
                                                           create_graph=create_graph if create_graph is not None else (
                                                                       (o == 1) or (o < order)),
                                                           retain_graph=retain_graph if retain_graph is not None else (
                                                                       o < order)
                                                           )
        return keys

    @staticmethod
    def _replace_nans(e, replacement):
        if replacement is True: replacement = 1e8
        if replacement is False or not nput.is_numeric(replacement):
            return e
        else:
            e = np.asanyarray(e)
            if e.ndim == 0:
                if np.isnan(e):
                    return replacement
                else:
                    return e
            else:
                e[np.isnan(e)] = replacement
            return e

    def process_aimnet_derivs(self, base_shape, coords, data, root_key, order,
                              extra_deriv_coord=None,
                              deriv_key=None,
                              property_shape=(),
                              pad_coord=1,
                              replace_nans=True,
                              output_shape=None):
        """
        Collect `data[root_key]` and its coordinate derivatives through `order` (and, with
        `extra_deriv_coord=(var, n)`, mixed derivatives through order `n` in `var`) as numpy arrays.
        """
        if extra_deriv_coord is not None:
            deriv_keys = []
            root_keys = []
            if deriv_key is not None: raise ValueError("mixed derivatives with a deriv_key not supported")
            var, suborder = extra_deriv_coord
            self._apply_autodiff(data, root_key, deriv_key, suborder, property_shape,
                                 diff_var=var,
                                 pad_coord=pad_coord,
                                 base_shape=base_shape,
                                 retain_graph=True,
                                 create_graph=True
                                 )
            for n in range(suborder+1):
                deriv, subshape = data[f"deriv_{n}"]
                if len(subshape) > 0 and len(deriv.shape) > 0:
                    deriv = deriv.__getitem__(subshape)
                data[f"subderiv_{n}"] = deriv
                dkeys = self._apply_autodiff(data, f"subderiv_{n}", None, suborder,
                                             data[f"subderiv_{n}"].shape,
                                             key_tag=f"deriv_{n}",
                                             pad_coord=pad_coord,
                                             base_shape=base_shape
                                             )
                deriv_keys += dkeys
                root_keys.append(f"subderiv_{n}")
        else:
            deriv_keys = self._apply_autodiff(data, root_key, deriv_key, order, property_shape,
                                              diff_var=None, pad_coord=pad_coord)
            root_keys = [root_key]

        if output_shape is None:
            output_shape = property_shape
        else:
            output_shape = tuple(output_shape)
        slice_shape = tuple(slice(None, o) for o in output_shape)
        for key in deriv_keys:
            deriv,subshape = data[key]
            if len(subshape + slice_shape) > 0:
                deriv = deriv.__getitem__(subshape + slice_shape)
            data[key] = deriv

        ko = self.model.keys_out
        afk = self.model.atom_feature_keys
        try:
            self.model.keys_out = self.model.keys_out + root_keys + deriv_keys
            if root_key not in self.model.keys_out:
                self.model.keys_out.append(root_key)
            self.model.atom_feature_keys = self.model.atom_feature_keys + deriv_keys
            expansion = self.model.process_output(data)
        finally:
            self.model.keys_out = ko
            self.model.atom_feature_keys = afk

        ndim = np.prod(coords.shape[-2:], dtype=int)

        if extra_deriv_coord is None:
            expansion = [
                self._replace_nans(
                    expansion[k].detach().cpu().numpy().reshape(
                        base_shape + (ndim,) * o + output_shape
                    ),
                    replace_nans
                )
                for o,k in enumerate(root_keys + deriv_keys)
            ]
        else:
            var, nsub = extra_deriv_coord
            targ_shape = list(var.shape[len(base_shape):])
            targ_shape[0] -= pad_coord
            vdim = np.prod(targ_shape, dtype=int)
            _ = []
            for n in range(extra_deriv_coord[1]+1):
                subexpansion = []
                for i in range(order + 1):
                    esub = expansion[f"subderiv_{n}" if i == 0 else f"deriv_{n}_{i}"].detach().cpu().numpy()
                    e_new = esub.reshape(
                        base_shape + (ndim,) * i + (vdim,) * n + output_shape
                    )
                    e_new = self._replace_nans(e_new, replace_nans)
                    subexpansion.append(e_new)
                _.append(subexpansion)
            expansion = _

        return expansion

    # endregion

    # region evaluation

    def evaluate_energy(self, coords, order, charge=None, multiplicity=None, **opts):
        """
        `[E, dE/dx, ...]` on a flat batch `(nstruct, natoms, 3)` (Hessians are taken one structure at a time)
        """
        if order < 2 or len(coords) == 1:
            base_shape, coords, data = self.prep_eval(coords, order, charge=charge, multiplicity=multiplicity, **opts)
            return self.process_aimnet_derivs(base_shape, coords, data, 'energy', order, deriv_key='forces')
        return _stack_structures(
            self.evaluate_energy(c[np.newaxis], order, charge=charge, multiplicity=multiplicity, **opts)
            for c in coords
        )

    def evaluate_charges(self, coords, order, charge=None, multiplicity=None, **opts):
        """
        `[q, dq/dx, ...]` on a flat batch `(nstruct, natoms, 3)` (derivatives are taken one
        structure at a time: the padded per-atom autograd doesn't separate batched molecules)
        """
        if order < 1 or len(coords) == 1:
            base_shape, coords, data = self.prep_eval(coords, order, charge=charge, multiplicity=multiplicity,
                                                      keep_graph=order > 0, forces=False, hessian=False, **opts)
            # AIMNet2 calculator charges are already in elementary-charge units.
            return self.process_aimnet_derivs(
                base_shape, coords, data,
                'charges',
                order,
                property_shape=(coords.shape[-2] + 1,),
                output_shape=(coords.shape[-2],)
            )
        return _stack_structures(
            self.evaluate_charges(c[np.newaxis], order, charge=charge, multiplicity=multiplicity, **opts)
            for c in coords
        )

    def evaluate_electrostatic_potential_derivatives(self, coords, coord_order, electrostatic_order,
                                                     charge=None, multiplicity=None, **opts):
        """
        The model's Coulomb energy and its derivatives with respect to the atomic charges (through
        `electrostatic_order`) and coordinates (through `coord_order`), as
        `[[E, dE/dx, ...], [dE/dq, ...], [d2E/dq2, ...]]`
        """
        self.model.model.outputs.lrcoulomb.key_out = 'coul_energy'
        base_shape, coords, data = self.prep_eval(coords, 0, charge=charge, multiplicity=multiplicity,
                                                  keep_graph=electrostatic_order > 0 or coord_order > 0,
                                                  **opts)
        data['coul_energy'] = data['coul_energy'].sum()

        return self.process_aimnet_derivs(
            base_shape, coords, data,
            'coul_energy',
            coord_order,
            property_shape=data['coul_energy'].shape[len(base_shape):],
            extra_deriv_coord=(data['charges'], electrostatic_order)
        )

    def evaluate_field_response(self, coords, coord_order=0, field_order=2, charge=None, multiplicity=None, **opts):
        """
        The fluctuating-charge dipole and polarizability expansions,
        `[[mu, dmu/dx, ...], [alpha, dalpha/dx, ...]][:field_order]`.

        :param coords: `(..., natoms, 3)`, Angstrom
        """
        charge, multiplicity = self.resolve_state(charge, multiplicity)
        coords = np.asanyarray(coords)
        potential_by_charges_by_coords = self.evaluate_electrostatic_potential_derivatives(
            coords, coord_order, field_order, charge=charge, multiplicity=multiplicity, **opts
        )
        return FluctuatingChargeResponse.expand(charge, coords, potential_by_charges_by_coords)

    def evaluate_properties(self, coords, orders, charge, multiplicity, **opts):
        out = {}
        for prop, order in orders.items():
            if prop == 'energy':
                out[prop] = self.evaluate_energy(coords, order, charge=charge, multiplicity=multiplicity, **opts)
            elif prop == 'charges':
                out[prop] = self.evaluate_charges(coords, order, charge=charge, multiplicity=multiplicity, **opts)
            elif prop == 'polarizability':
                out[prop] = self.evaluate_field_response(coords, order, 2, charge=charge, multiplicity=multiplicity, **opts)[1]
            else:
                raise ValueError(f"unknown property {prop!r}")
        return out

    # endregion


@ElectronicStructureBackend.register('hipnn')
class HIPNNBackend(ElectronicStructureBackend):
    """
    A `hippynn` model loaded from a training directory, with energies (`energy_key`, gradients
    from `gradient_key`, higher derivatives by autograd) and dipoles (`dipole_key`).
    """
    distance_units = 'Angstroms'
    property_units = {
        'energy': 'Kilocalories/Mole',
        'dipole': ('ElementaryCharge', 'Angstroms'),
    }
    analytic_derivative_orders = {
        'energy': 1,
        'dipole': 0,
        'polarizability': 1,  # charge-equilibration models; see `evaluate_field_response`
    }
    property_shapes = {'dipole': (3,), 'polarizability': (3, 3)}
    energy_key = 'energies'
    gradient_key = 'Grad'
    dipole_key = 'dipole'
    dipole_gradient_key = 'dipole_grad'

    def __init__(self, atoms, model_dir,
                 energy_key=dev.default,
                 gradient_key=dev.default,
                 dipole_key=dev.default,
                 charge=None, multiplicity=None, device=None, quiet=True):
        """
        :param model_dir: the directory `hippynn` saved the experiment to
        :param energy_key: the graph node holding the energy (default `'energies'`)
        :param gradient_key: the graph node holding its gradient (default `'Grad'`; `None` to use autograd)
        :param dipole_key: the graph node holding the dipole (default `'dipole'`)
        """
        super().__init__(atoms, charge=charge, multiplicity=multiplicity, quiet=quiet)
        self.model_dir = model_dir
        self.model = self.setup_model(model_dir, device=device)
        if not dev.is_default(energy_key, allow_None=False):
            self.energy_key = energy_key
        if not dev.is_default(gradient_key, allow_None=False):
            self.gradient_key = gradient_key
        if not dev.is_default(dipole_key, allow_None=False):
            self.dipole_key = dipole_key

    @classmethod
    def handle_specialization(cls, tag):
        return {'model_dir': tag}

    @classmethod
    def setup_model(cls, model_dir, *, device=None, **opts):
        import torch
        import hippynn.graphs
        import hippynn.experiment.serialization
        cur_dir = os.getcwd()

        device = resolve_torch_device(device)
        try:
            os.chdir(model_dir)
            with torch.no_grad():
                model = hippynn.experiment.serialization.load_model_from_cwd(model_device=device)
        finally:
            os.chdir(cur_dir)

        predictor = hippynn.graphs.Predictor.from_graph(model)
        return predictor

    def print_model_properties(self):
        self.model.graph.print_structure()

    def enable_charge_equilibration_outputs(self):
        """
        Rebuild the predictor so it also exposes the charge-equilibration (`ChEQ`) Coulomb energy
        and its inputs, for fluctuating-charge responses.
        """
        import hippynn.graphs
        predictor = hippynn.graphs.Predictor.from_graph(self.model.graph,
                                                        additional_outputs=[
                                                            'ChEQ.coul_energy'
                                                        ]
                                                        )
        cheq = predictor.graph.node_from_name('ChEQ')
        predictor.outputs.extend(cheq.parents[-2:])
        self.model = predictor

    @staticmethod
    def autodiff(expr, coord_list, pad_dim, create_graph=False, retain_graph=False):
        import torch
        shape = coord_list.shape + expr.shape[1:]
        grads = []
        npad = len(pad_dim)
        rav_shape = tuple(
            (expr.shape[1+2*i]*expr.shape[2+2*i])
            for i in range((expr.ndim - 1 - npad) // 2)
        )
        cache = {}
        if len(rav_shape) > 0:
            inds = np.array(np.unravel_index(np.arange(np.prod(expr.shape[1:-npad], dtype=int)), rav_shape)).T
            shinds = np.sort(inds, axis=1)
            for n,subexpr in enumerate(expr):
                for i,_f in enumerate(subexpr.flatten()):
                    ravioli = tuple(inds[i])
                    key = tuple(shinds[i])
                    if ravioli == key:
                        g = torch.autograd.grad(_f,
                                                coord_list,
                                                retain_graph=True,
                                                create_graph=create_graph)
                        g = g[0][n]
                        grads.append(g)
                        cache[key] = g
                    else:
                        grads.append(cache[key])
        else:
            for n,_f in enumerate(expr):
                g = torch.autograd.grad(_f,
                                        coord_list,
                                        retain_graph=True,
                                        create_graph=create_graph)[0]
                grads.append(g[n])

        deriv = torch.stack(grads).view(shape)
        shape_tuple = (slice(None),)
        return deriv, shape_tuple

    def process_derivatives(self, expr, coords, order):
        terms = [(expr, slice(None,))]
        pad_dim = expr.shape[1:]
        for o in range(1, order + 1):
            terms.append(
                self.autodiff(terms[-1][0],
                              coords,
                              pad_dim,
                              create_graph=(o < order),
                              retain_graph=(o < order)
                              )
            )
        return [
            deriv.__getitem__(subshape)
            for deriv, subshape in terms
        ]

    def process_results(self, base_shape, coords, outputs, order,
                        property_key=None,
                        gradient_key=None,
                        property_shape=(),
                        output_shape=None,
                        **etc):
        if property_key is None:
            property_key = self.energy_key
        property_shape = tuple(property_shape)
        terms = [outputs[property_key]]
        if order > 0:
            if gradient_key is not None:
                g = outputs[gradient_key]
                terms.append(g)
                order = order - 1
            if order > 0:
                terms.extend(
                    self.process_derivatives(
                        terms[-1],
                        coords,
                        order
                    )[1:]
                )

        n = coords.shape[-1] * coords.shape[-2]
        if output_shape is None:
            output_shape = property_shape
        else:
            output_shape = tuple(output_shape)
        return [
            e.detach().numpy().reshape(base_shape + (n,)*i + output_shape)
            for i,e in enumerate(terms)
        ]

    def prep_eval(self, coords, order, property_key=dev.default, gradient_key=dev.default, **opts):
        """
        Run the predictor on `coords` (`(..., natoms, 3)`, Angstrom), adding `property_key` and
        `gradient_key` to its outputs if needed.

        :return: `(base_shape, torch_coords, outputs)`
        """
        import torch

        if dev.is_default(property_key, allow_None=False):
            property_key = self.energy_key
        if dev.is_default(gradient_key, allow_None=False):
            gradient_key = self.gradient_key

        base_shape = coords.shape[:-2]
        coords = coords.reshape((-1,) + coords.shape[-2:])

        torch_coords = torch.tensor(coords, dtype=torch.float32, device=self.model.model_device)
        arg_dict = {
            'coordinates': torch_coords,
            'species': torch.tensor(
                np.repeat(np.array([self.numbers]), coords.shape[0], axis=0),
                dtype=torch.int64, device=self.model.model_device
            )
        }

        req_grad = self.model.requires_grad
        outputs = None
        output_names = None
        output_dbnames = None
        try:
            self.model.requires_grad = order > 0
            if property_key is not None and property_key not in self.model.out_names:
                prop = self.model.graph.node_from_name(property_key)
                prop.bname = prop.db_name
                if outputs is None:
                    outputs = list(self.model.outputs)
                    output_names = list(self.model.out_names)
                    output_dbnames = list(self.model.out_dbnames)
                self.model.add_output(prop)
            if gradient_key is not None and gradient_key not in self.model.out_names:
                prop = self.model.graph.node_from_name(gradient_key)
                prop.bname = prop.db_name
                if outputs is None:
                    outputs = list(self.model.outputs)
                    output_names = list(self.model.out_names)
                    output_dbnames = list(self.model.out_dbnames)
                self.model.add_output(gradient_key)

            with quiet_mode(self.quiet):
                with torch.jit.optimized_execution(False):
                    data = self.model(**arg_dict, **opts)
        finally:
            self.model.requires_grad = req_grad
            if outputs is not None:
                self.model.outputs[:] = outputs
                self.model.out_names[:] = output_names
                self.model.out_dbnames[:] = output_dbnames

        return base_shape, torch_coords, data

    def evaluate_energy(self, coords, order, **opts):
        if order < 2 or len(coords) == 1:
            base_shape, tcoords, data = self.prep_eval(coords, order, **opts)
            return self.process_results(base_shape, tcoords, data, order,
                                        property_key=self.energy_key, gradient_key=self.gradient_key)
        return _stack_structures(self.evaluate_energy(c[np.newaxis], order, **opts) for c in coords)

    def evaluate_dipole(self, coords, order, **opts):
        base_shape, tcoords, data = self.prep_eval(coords, order, property_key=self.dipole_key, gradient_key=None, **opts)
        return self.process_results(base_shape, tcoords, data, order,
                                    property_key=self.dipole_key,
                                    gradient_key=self.dipole_gradient_key,
                                    property_shape=(tcoords.shape[-1],))

    def evaluate_field_response(self, coords, coord_order=0, field_order=2, charge=None, multiplicity=None, **opts):
        # The pre-migration `HIPNNDipolePolarizabilityEvaluator` called AIMNet2-only machinery
        # (`process_aimnet_derivs`, `model.outputs.lrcoulomb`) on a HIPNN predictor, so it never ran.
        # Getting the ChEQ hardness/electronegativity derivatives out of a hippynn graph still needs
        # to be written.
        raise NotImplementedError("fluctuating-charge responses aren't implemented for HIPNN models yet")

    def evaluate_properties(self, coords, orders, charge, multiplicity, **opts):
        out = {}
        for prop, order in orders.items():
            if prop == 'energy':
                out[prop] = self.evaluate_energy(coords, order, **opts)
            elif prop == 'dipole':
                out[prop] = self.evaluate_dipole(coords, order, **opts)
            elif prop == 'polarizability':
                out[prop] = self.evaluate_field_response(coords, order, 2, charge=charge, multiplicity=multiplicity, **opts)[1]
            else:
                raise ValueError(f"unknown property {prop!r}")
        return out


@ElectronicStructureBackend.register('ase')
class ASECalculatorBackend(ElectronicStructureBackend):
    """
    Any ASE calculator. Energies come with forces and (if the calculator has `get_hessian`)
    Hessians; any other result the calculator produces (`'charges'`, `'dipole'`, model-specific
    keys) is available without derivatives, and is computed in the same calculator pass as the
    energy and forces.
    """
    distance_units = 'Angstroms'
    default_multiplicity = 1
    property_units = {
        'energy': 'ElectronVolts',
        'charges': 'ElementaryCharge',
        'dipole': ('ElementaryCharge', 'Angstroms'),
    }
    analytic_derivative_orders = {
        'energy': 2,
    }
    property_shapes = {
        'charges': lambda n: (n,),
        'dipole': (3,),
    }
    hessian_function_attribute = 'get_hessian'

    def __init__(self, atoms, charge=None, multiplicity=None, quiet=True, calculator=None, info=None, calc=None,
                 **calculator_settings):
        """
        :param calculator: an ASE calculator (built by `setup_calculator(**calculator_settings)` if not given)
        :param info: extra entries for the ASE `Atoms.info` dict (charge and spin are set from the state)
        :type info: dict | None
        :param calc: an alias for `calculator`
        """
        super().__init__(atoms, charge=charge, multiplicity=multiplicity, quiet=quiet)
        if calculator is None:
            calculator = calc
        if calculator is None:
            calculator = self.setup_calculator(**calculator_settings)
        elif len(calculator_settings) > 0:
            raise ValueError(f"got calculator settings {list(calculator_settings)} along with a calculator")
        self.model = calculator
        self.info = {} if info is None else dict(info)

    @classmethod
    def setup_calculator(cls, **settings):
        raise ValueError(f"{cls.__name__} needs an ASE `calculator`")

    @property
    def calculator(self):
        return self.model

    def get_analytic_derivative_order(self, prop):
        # anything else the calculator returns is available, without derivatives
        return self.analytic_derivative_orders.get(prop, 0)

    def prep_ase(self, coords, charge=None, multiplicity=None, info=None):
        """
        An `ASEMolecule` at `coords` with this backend's calculator attached.
        """
        from McUtils.ExternalPrograms import ASEMolecule
        charge, multiplicity = self.resolve_state(charge, multiplicity)
        return ASEMolecule.from_coords(
            self.atoms,
            coords.reshape((-1,) + coords.shape[-2:])[0],
            calculator=self.model,
            charge=charge,
            spin=multiplicity,
            info=dict(self.info, **({} if info is None else info))
        )

    def calculate_results(self, coords, keys, charge=None, multiplicity=None, info=None):
        """
        The calculator's raw `results` for `keys` at each structure of `coords`.

        :param coords: `(..., natoms, 3)`, Angstrom
        :return: arrays keyed by result name, with the leading batch shape of `coords`
        :rtype: dict
        """
        coords = np.asanyarray(coords, dtype=float)
        base_shape = coords.shape[:-2]
        flat = coords.reshape((-1,) + coords.shape[-2:])
        mol = self.prep_ase(flat, charge=charge, multiplicity=multiplicity, info=info)
        res = mol.calculate_props(list(keys), geoms=flat)
        return {k: v.reshape(base_shape + v.shape[1:]) for k, v in res.items()}

    def evaluate_properties(self, coords, orders, charge, multiplicity, info=None, **opts):
        e_order = orders.get('energy')
        others = [p for p in orders if p != 'energy']
        out = {}
        nstruct = coords.shape[0]
        ncoord = 3 * self.natoms
        if e_order is not None and e_order > 1:
            mol = self.prep_ase(coords, charge=charge, multiplicity=multiplicity, info=info)
            out['energy'] = list(mol.calculate_energy(coords, order=e_order, hessian_func_attr=self.hessian_function_attribute))
            keys = others
        else:
            keys = list(others)
            if e_order is not None:
                keys = ['energy'] + keys
                if e_order > 0:
                    keys.append('forces')
        if len(keys) > 0:
            res = self.calculate_results(coords, keys, charge=charge, multiplicity=multiplicity, info=info)
            if 'energy' not in out and e_order is not None:
                out['energy'] = [np.asanyarray(res['energy']).reshape(nstruct)]
                if e_order > 0:
                    out['energy'].append(-res['forces'].reshape(nstruct, ncoord))
            for p in others:
                out[p] = [res[p]]
        return out


@ElectronicStructureBackend.register('mace')
class MACEBackend(ASECalculatorBackend):
    """
    MACE foundation models through `mace.calculators` (`'extra_large'` = MACE-OMOL by default)
    """
    model_types = {
        'extra_large':'omol/extra_large',
        "small":'off/small',
        "medium":'off/medium',
        "large":'off/large',
        'anicc':'anicc/none',
        'off':'off/large',
        'mp':'mp/medium-mpa-0',
        'polar':'polar/polar-1-l'
    }
    _mace_cache = {}  # weird OOM errors...

    @classmethod
    def handle_specialization(cls, tag):
        return {'model': tag}

    @classmethod
    def _overload_mace_modeldir(cls, root_dir):
        if root_dir is None:
            root_dir = os.environ.get("MODEL_CACHE_DIR")
        if root_dir is not None:
            root_dir = os.fspath(root_dir)
        return environment_override(XDG_CACHE_HOME=root_dir)

    @classmethod
    def setup_calculator(cls, model='extra_large', model_type=None, device=None, model_dir=None, **settings):
        model = os.fspath(model)
        if model_type is None:
            model = model.lower()
            model = cls.model_types.get(model, model)
            if '/' not in model or os.path.isabs(model):
                raise ValueError("set model_type when supplying a MACE checkpoint path or an unqualified model")
            model_type, model = model.split('/', 1)
        if isinstance(model_type, str):
            with quiet_mode():
                if model_type == 'omol':
                    from mace.calculators import mace_omol as model_type
                elif model_type == 'mp':
                    from mace.calculators import mace_mp as model_type
                elif model_type == 'anicc':
                    from mace.calculators import mace_anicc as model_type
                elif model_type == 'polar':
                    from mace.calculators import mace_polar as model_type
                elif model_type == 'off':
                    from mace.calculators import mace_off as model_type
                else:
                    import mace.calculators
                    model_type = getattr(mace.calculators, model_type)

        device = resolve_torch_device(device)
        with quiet_mode(), cls._overload_mace_modeldir(model_dir):
            cache_root = os.path.abspath(os.environ.get(
                'XDG_CACHE_HOME', os.path.join(os.path.expanduser('~'), '.cache')
            ))
            cache_key = (model_type, model, str(device), cache_root)
            try:
                hash(cache_key)
            except TypeError:
                cacheable = False
            else:
                cacheable = len(settings) == 0
            if not cacheable or cache_key not in cls._mace_cache:
                calc = model_type(model=model, device=device, **settings)
                if cacheable:
                    cls._mace_cache[cache_key] = calc
            else:
                calc = cls._mace_cache[cache_key]

        return calc


@ElectronicStructureBackend.register('mace-polar')
@ElectronicStructureBackend.register('mace_polar')
class MACEPolarBackend(MACEBackend):
    """
    MACE-POLAR: energies plus the Gaussian-smeared atomic multipoles it predicts
    (`density_coefficients`, real `m = -l..l` moments in e Angstrom^l in the `graph_longrange`
    GTO normalization) and their spin channels (`spin_charge_density`), with smearing width
    `smearing` (Angstrom).
    """
    model_types = {
        'small': 'polar-1-s',
        'medium': 'polar-1-m',
        'large': 'polar-1-l',
    }
    default_model = 'polar-1-l'
    property_units = dict(
        ASECalculatorBackend.property_units,
        density_coefficients=None,  # e Angstrom^l per l
        spin_charge_density=None,
        electrostatic_energy='ElectronVolts',
    )
    property_shapes = dict(ASECalculatorBackend.property_shapes)

    def __init__(self, atoms, model=None, charge=None, multiplicity=None, external_field=None,
                 calculator=None, device=None, model_dir=None, default_dtype='float64', quiet=True,
                 smearing=None, info=None, **settings):
        """
        :param model: `'polar-1-s'`/`'polar-1-m'`/`'polar-1-l'` (or `'small'`/`'medium'`/`'large'`)
        :param external_field: uniform external field (`atoms.info["external_field"]`), in the model's native units
        :param calculator: an existing MACE-POLAR ASE calculator to reuse
        :param smearing: override for the model's smearing width (Angstrom)
        """
        if calculator is None:
            if model is None:
                model = self.default_model
            model = self.model_types.get(model, model)
            settings = dict(settings, model=f'polar/{model}', device=device, model_dir=model_dir, default_dtype=default_dtype)
        super().__init__(atoms, calculator=calculator, charge=charge, multiplicity=multiplicity, quiet=quiet,
                         info=info, **settings)
        self.external_field = None if external_field is None else np.asanyarray(external_field, dtype=float)
        if self.external_field is not None:
            self.info['external_field'] = list(self.external_field)
        self._smearing = smearing

    @classmethod
    def handle_specialization(cls, tag):
        return {'model': tag}

    @property
    def model_smearing(self):
        """the smearing width the loaded model uses (Angstrom)"""
        models = getattr(self.model, 'models', None)
        if models:
            return getattr(models[0], 'atomic_multipoles_smearing_width', None)
        return None

    @property
    def smearing(self):
        """the smearing width in use (Angstrom): the override, else the model's"""
        return self._smearing if self._smearing is not None else self.model_smearing

    polar_properties = ('energy', 'density_coefficients')

    def evaluate_polar_properties(self, coords, properties=None, charge=None, multiplicity=None):
        """
        Runs the MACE-POLAR calculator and returns its raw results.

        :param coords: nuclear coordinates, `(..., natoms, 3)`, in Angstrom
        :param properties: calculator result keys (e.g. `'electrostatic_energy'`, `'spin_charge_density'`, `'dipole'`)
        :param charge: total charge override for this call
        :param multiplicity: spin multiplicity override for this call
        :return: arrays keyed by property, each with the leading batch shape of `coords`
        :rtype: dict
        """
        if properties is None:
            properties = self.polar_properties
        return self.calculate_results(coords, properties, charge=charge, multiplicity=multiplicity)


@ElectronicStructureBackend.register('uma')
class UMABackend(ASECalculatorBackend):
    """
    FAIRChem UMA models (`'uma-s-1p2'` on the `omol` task by default)
    """
    model_types = {
        'medium':'uma-m-1p1/omol',
        'small':'uma-s-1p2/omol',
        'crystals':'uma-s-1p2/omc',
        'crystals/small':'uma-s-1p2/omc',
        'crystals/medium':'uma-m-1p1/omc',
        'inorganics/small':'uma-s-1p2/omat',
        'inorganics/medium':'uma-m-1p1/omat',
    }

    @classmethod
    def handle_specialization(cls, tag):
        return {'model': tag}

    @classmethod
    def setup_calculator(cls, model="uma-s-1p2",
                         device=None, task_name=None, model_dir=None,
                         login=None,
                         hf_token=None,
                         predict_unit_options=None, **settings):
        from fairchem.core import pretrained_mlip, FAIRChemCalculator

        if login is None:
            if hf_token is None:
                hf_token = os.environ.get("HF_TOKEN")
            login = hf_token is not None
        if login:
            hf_authenticate(hf_token)

        model = cls.model_types.get(model, model)
        if task_name is None:
            if '/' in model:
                model, task_name = model.split('/', 1)
            else:
                task_name = 'omol'

        if predict_unit_options is None:
            predict_unit_options = {}
        predict_unit_options = predict_unit_options.copy()
        if model_dir is None:
            model_dir = os.environ.get("MODEL_CACHE_DIR")
        if model_dir is not None and 'cache_dir' not in predict_unit_options:
            predict_unit_options['cache_dir'] = os.path.join(model_dir, 'fairchem')

        model = model.lower()
        with quiet_mode():
            if 'device' not in predict_unit_options:
                device = resolve_torch_device(device)
                predict_unit_options['device'] = device
            predictor = pretrained_mlip.get_predict_unit(model, **predict_unit_options)
            calc = FAIRChemCalculator(predictor, task_name=task_name, **settings)

        return calc


@ElectronicStructureBackend.register('upet')
class UPETBackend(ASECalculatorBackend):
    """
    UPET (PET-MAD and relatives) models through `upet.calculator`
    """
    model_types = {
        # --- Current / recommended ---
        "mad": "pet-mad-s",
        "mad-xs": "pet-mad-xs",

        "oam": "pet-oam-xl",
        "oam-l": "pet-oam-l",

        "omat-xs": "pet-omat-xs",
        "omat-s": "pet-omat-s",
        "omat-m": "pet-omat-m",
        "omat-l": "pet-omat-l",
        "omat-xl": "pet-omat-xl",

        "omatpes": "pet-omatpes-l",

        "spice": "pet-spice-s",
        "spice-l": "pet-spice-l",

        # --- Legacy ---
        "omad-xs": "pet-omad-xs",
        "omad-s": "pet-omad-s",
        "omad-l": "pet-omad-l",
    }

    @classmethod
    def handle_specialization(cls, tag):
        return {'model': tag}

    @classmethod
    def setup_calculator(cls, model="pet-mad-s", device=None, version=None, login=None,
                         hf_token=None,
                         model_dir=None,
                         **settings):
        from upet.calculator import UPETCalculator

        if model_dir is None:
            model_dir = os.environ.get("MODEL_CACHE_DIR")

        model = cls.model_types.get(model, model)
        if version is None:
            if '/' in model:
                model, version = model.split('/', 1)
            else:
                version = 'lastest'

        with environment_override(HF_HUB_CACHE=model_dir):
            if login is None:
                if hf_token is None:
                    hf_token = os.environ.get("HF_TOKEN")
                login = hf_token is not None
            if login:
                hf_authenticate(hf_token)

            model = model.lower()
            with quiet_mode():
                if device is None:
                    device = resolve_torch_device(device)
                calc = UPETCalculator(model=model, version=version, device=device, **settings)

        return calc


@ElectronicStructureBackend.register('orb')
class OrbBackend(ASECalculatorBackend):
    """
    Orb models through `orb_models` (`'omol'` = `orbmol-v3-conservative-omol` by default)
    """
    model_types = {
        'v1':"orbmol-v1-conservative",
        'v2':"orbmol-v2",
        'omol':"orbmol-v3-conservative-omol",
        'v3':"orbmol-v3-conservative-omol",
        'v3-direct':"orbmol-v3-direct-omol",
        'v1-direct':"orbmol-v1-direct",
        'omat':"orb-v3-conservative-inf-omat",
        'v3-omat':"orb-v3-conservative-inf-omat",
        'v3-omat-direct':"orb-v3-direct-inf-omat",
        'v3-omat20':"orb-v3-conservative-20-omat",
        'v3-omat20-direct':"orb-v3-direct-20-omat",
        'd3':"separate-d3-5layer",
        'd3-3layer':"separate-d3-3layer",
        'd4':"separate-d3-5layer",
        'd4-3layer':"separate-d4-3layer"
    }

    @classmethod
    def handle_specialization(cls, tag):
        return {'model': tag}

    @classmethod
    def setup_calculator(cls, model='omol', device=None, version=None, login=None, **model_opts):
        from orb_models.forcefield import pretrained, ORB_PRETRAINED_MODELS
        from orb_models.forcefield.inference.calculator import ORBCalculator

        model = cls.model_types.get(model, model)
        model = ORB_PRETRAINED_MODELS.get(cls.model_types.get(model, model), model)

        device = "cpu"  # or device="cuda"
        orbff, atoms_adapter = ORB_PRETRAINED_MODELS[model](
            device=device,
            **model_opts,
        )

        calc = ORBCalculator(orbff, atoms_adapter=atoms_adapter, device=device)
        return calc
