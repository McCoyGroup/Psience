import abc
import itertools
import math
import sys, os
import numpy as np
import warnings

from McUtils.Data import AtomData, UnitsData
from McUtils.Zachary import TensorDerivativeConverter, FiniteDifferenceDerivative, CoordinateFunction
import McUtils.Numputils as nput
import McUtils.Devutils as dev
import McUtils.Coordinerds as coordops
from McUtils.Scaffolding import Logger

from ..Data import PotentialSurface, DipoleSurface
from .CoordinateSystems import MolecularEmbedding
from .Properties import NormalModesManager

__all__ = [
    "MolecularEvaluator",
    "EnergyEvaluator",
    "RDKitEnergyEvaluator",
    "AIMNet2EnergyEvaluator",
    "XTBEnergyEvaluator",
    "PotentialFunctionEnergyEvaluator",
    "DipoleEvaluator",
    "DipoleFunctionDipoleEvaluator"
]

class MolecularEvaluator:
    def __init__(self, embedding:MolecularEmbedding, normal_modes:NormalModesManager):
        self.embedding = embedding
        self.normal_modes = normal_modes
    @property
    def coords(self):
        return self.embedding.coords

    @property
    def masses(self):
        return self.embedding.masses

    def evaluate(self,
                 func,
                 use_internals=None,
                 order=None,
                 strip_embedding=False
                 ):
        if use_internals is None:
            use_internals = self.embedding.internals is not None
        if use_internals:
            coords = self.embedding.internal_coordinates
            if strip_embedding:
                embedding_coords = self.embedding.embedding_coords
                if embedding_coords is not None:
                    good_coords = np.setdiff1d(np.arange(3 * len(self.embedding.masses)), embedding_coords)
                    coords = coords.reshape(coords.shape[:-2] + (3 * len(self.embedding.masses),))
                    coords = coords[..., good_coords]
            if order is None:
                return func(coords).view(np.ndarray)

            terms = func(coords, order=order)
            # raise Exception([t.shape for t in terms], coords.shape)
            if strip_embedding:
                embedding_coords = self.embedding.embedding_coords
                if embedding_coords is not None:
                    good_coords = np.setdiff1d(np.arange(3 * len(self.embedding.masses)), embedding_coords)

                    const = terms[0]
                    terms = terms[1:]
                    new = []
                    ncs = 3 * len(self.embedding.masses)
                    if self.embedding.coords.ndim > 2:
                        npts = len(coords)
                    else:
                        npts = 1
                    for n, ders in enumerate(terms):
                        dt = np.zeros((npts,) + (ncs,)*(n+1))
                        idx_pos = (...,) + np.ix_(*[good_coords]*(n+1))
                        dt[idx_pos] = ders.view(np.ndarray)
                        if self.embedding.coords.ndim == 2:
                            dt = dt[0]
                        new.append(dt)
                    terms = [const] + new

            const = terms[0]
            jacs = self.embedding.get_internals_by_cartesians(order)

            terms = TensorDerivativeConverter(
                jacs,
                terms[1:],
                jacobians_name='dXdR',
                values_name='f'
            ).convert()  # , check_arrays=True)

            return [const.view(np.ndarray)] + [t.view(np.ndarray) for t in terms]
        else:
            if order is None:
                return func(self.embedding.coords).view(np.ndarray)
            else:
                return [x.view(np.ndarray) for x in func(self.embedding.coords, order=order)]

    def evaluate_at(self,
                    func,
                    coords,
                    use_internals=None,
                    order=None,
                    strip_embedding=False
                    ):
        return type(self)(
            MolecularEmbedding(self.embedding.masses, coords, self.embedding.internals),
            self.normal_modes
        ).evaluate(func, use_internals=use_internals, order=order, strip_embedding=strip_embedding)

    def get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None,
                                  use_internals=False,
                                  coordinate_expansion=None,
                                  strip_embedding=False,
                                  shift=True
                                  ):
        displacements = np.asanyarray(displacements)

        if which is not None:
            which = tuple(
                np.ravel_multi_index(idx, (len(self.embedding.masses), 3))
                    if not isinstance(idx, (int, np.integer)) else
                idx
                for idx in which
            )

        if use_internals and self.embedding.internals is None:
            raise ValueError("can't displace in internals without internal coordinate spec")
        base_coords = self.embedding.coords if not use_internals else self.embedding.internal_coordinates
        if strip_embedding:
            ecs = self.embedding.embedding_coords
            all_coords = np.arange(len(self.embedding.masses) * 3)
            which = np.setdiff1d(all_coords, ecs)[which,]

        if coordinate_expansion is not None:
            if which is not None:
                base_disps = np.zeros(displacements.shape[:-1] + coordinate_expansion[0].shape[:1])
                base_disps[..., which] = displacements
                displacements = base_disps
            _ = []
            new_disps = 0
            shared = displacements.ndim-1
            for n,disp in enumerate(coordinate_expansion):
                disp = disp / math.factorial(n+1)
                if nput.is_numeric(disp) and disp == 0: continue
                for i in range(n+1):
                    if i == 0:
                        disp = np.tensordot(displacements, disp, axes=[-1, 0])
                    else:
                        disp = nput.vec_tensordot(displacements, disp, axes=[-1, shared], shared=shared)
                new_disps = new_disps + disp
            if use_internals:
                displacements = new_disps.reshape(new_disps.shape[:-1] + (-1,))
            else:
                displacements = new_disps.reshape(new_disps.shape[:-1] + (-1, 3))
            which = None

        if which is not None:
            if displacements.shape[-1] != len(which):  # displacements provided in atom coordinates
                displacements = displacements.reshape(
                    displacements.shape[:-2] +
                    (np.prod(displacements.shape[-2:], dtype=int),)
                )
            if displacements.ndim > 1:
                for _ in range(displacements.ndim - 1):
                    base_coords = np.expand_dims(base_coords, 0)
                base_coords = np.broadcast_to(base_coords, displacements.shape[:-1] + base_coords.shape[-2:])
            base_coords = base_coords.copy()
            flat_coords = base_coords.reshape(displacements.shape[:-1] + (-1,))

            if shift:
                flat_coords[..., which] += displacements
            else:
                flat_coords[..., which] = displacements
            base_coords = flat_coords.reshape(base_coords.shape)

        elif sel is not None or axes is not None:
            if displacements.ndim > 2:
                for _ in range(displacements.ndim - 2):
                    base_coords = np.expand_dims(base_coords, 0)
                base_coords = np.broadcast_to(base_coords, displacements.shape[:-2] + base_coords.shape[-2:])
            base_coords = base_coords.copy()

            if sel is None:
                sel = np.arange(len(self.embedding.masses))
            if axes is None:
                axes = np.arange(3)

            sel = np.asanyarray(sel)[:, np.newaxis]
            axes = np.asanyarray(axes)[np.newaxis, :]

            if shift:
                base_coords[..., sel, axes] += displacements
            else:
                base_coords[..., sel, axes] = displacements

        else:
            if displacements.shape[-2:] != base_coords.shape:
                raise ValueError("displacements with shape {} passed but coordinates have shape {}".format(
                    displacements.shape,
                    base_coords.shape
                ))
            if shift:
                base_coords = np.expand_dims(base_coords, list(range(displacements.ndim - 2)))
                base_coords = base_coords + displacements
            else:
                base_coords = displacements

        if use_internals:
            # track the embedding info...
            base_coords = self.embedding.internal_coordinates.system(base_coords,
                                                                     **self.embedding.internal_coordinates.converter_options
                                                                     )
            if isinstance(use_internals, str):
                if use_internals == 'convert':
                    base_coords = base_coords.convert(self.embedding.coords.system)
                elif use_internals == 'reembed':
                    base_coords = self.embedding.embed_coords(base_coords.convert(self.embedding.coords.system))
        else:
            base_coords = self.embedding.coords.system(base_coords)
        return base_coords

    def get_scan_coordinates(self,
                             domains,
                             internals=False,
                             which=None, sel=None, axes=None,
                             coordinate_expansion=None,
                             strip_embedding=False,
                             shift=True,
                             return_displacements=False
                             ):

        displacement_mesh = np.moveaxis(
            np.array(
                np.meshgrid(*[np.linspace(*d) for d in domains], indexing='ij')
            ),
            0, -1
        )
        disps = self.get_displaced_coordinates(displacement_mesh, shift=shift,
                                              use_internals=internals, which=which, sel=sel, axes=axes,
                                              coordinate_expansion=coordinate_expansion,
                                              strip_embedding=strip_embedding
                                              )
        if return_displacements:
            return displacement_mesh, disps
        else:
            return disps

    def get_nearest_displacement_atoms(self,
                                       points,
                                       sel=None, axes=None, weighting_function=None,
                                       return_distances=False
                                       ):

        pts = np.asanyarray(points)
        smol = pts.ndim == 1
        if smol: pts = pts[np.newaxis]
        pts = pts.reshape(-1, pts.shape[-1])

        if axes is None:
            axes = [0, 1, 2]
        axes = np.asanyarray(axes)

        if sel is None:
            sel = np.arange(len(self.masses))
        sel = np.asanyarray(sel)

        ref = self.coords
        masses = self.masses
        dists = np.linalg.norm(
            pts[:, np.newaxis, :] - ref[np.newaxis, sel[:, np.newaxis], axes[np.newaxis, :]],
            axis=-1
        )
        if weighting_function is None:
            weighting_function = np.sqrt
        dists = dists * weighting_function(masses)[np.newaxis, sel]
        nearest = np.argmin(dists, axis=-1)
        atom_idx = sel[nearest]

        if return_distances:
            return atom_idx, dists[np.arange(len(nearest)), nearest]
        else:
            return atom_idx

    def get_nearest_displacement_coordinates(self,
                                             points,
                                             sel=None, axes=None, weighting_function=None,
                                             modes_nearest=False,
                                             return_distances=False
                                             ):
        """
        Displaces the _nearest_ atom (in a mass-weighted sense) to the given point
        This allows for the development of functions / potentials that are in the small-displacement limit
        where everything _except_ the nearest atom to the given Cartesian position is at its equilibrium value

        :param points:
        :param sel:
        :param axes:
        :param weighting_function:
        :return:
        """

        pts = np.asanyarray(points)
        smol = pts.ndim == 1
        if smol: pts = pts[np.newaxis]
        base_shape = pts.shape[:-1]
        pts = pts.reshape(-1, pts.shape[-1])

        if axes is None:
            axes = [0, 1, 2]
        axes = np.asanyarray(axes)

        if sel is None:
            sel = np.arange(len(self.masses))
        sel = np.asanyarray(sel)

        ref = self.coords
        atom_idx = self.get_nearest_displacement_atoms(
            points,
            sel=sel, axes=axes, weighting_function=weighting_function,
            return_distances=return_distances
        )
        if return_distances:
            atom_idx, dists = atom_idx

        if modes_nearest:
            npts = len(points)
            diag = np.arange(npts)
            mode_origin = self.normal_modes.modes.basis.origin
            mode_origin = mode_origin[:, axes]
            origin_coords = np.broadcast_to(mode_origin[np.newaxis], (npts,) + mode_origin.shape)

            origin_coords = origin_coords[diag, atom_idx]
            nearest_disps = pts - origin_coords

            modes = self.normal_modes.modes.basis.matrix
            mat = np.reshape(modes, (1, len(self.masses), -1, modes.shape[-1]))
            mat = mat[:, :, axes, :]
            mat = np.broadcast_to(mat, (npts,) + mat.shape[1:])
            mat_bits = mat[diag, atom_idx, :, :]  # N x 3 x N_modes
            pseudo_inverses = np.linalg.inv(mat_bits @ mat_bits.transpose(0, 2, 1))
            # print(mat_bits.shape, pseudo_inverses.shape, nearest_disps.shape)
            ls_coords = mat_bits.transpose(0, 2, 1) @ pseudo_inverses @ nearest_disps[:, :, np.newaxis]
            # ls_coords = np.reshape(ls_coords, ls_coords.shape[:2])
            # print(modes.shape, ls_coords.shape, pseudo_inverses.shape, modes.shape)
            test = modes[np.newaxis, :, :] @ ls_coords
            ls_coords = np.reshape(ls_coords, ls_coords.shape[:2])
            # print(test.shape, nearest_disps.shape)
            # print(test.reshape(test.shape[0], -1, 3))
            # print(nearest_disps)
            # print(pts)
            # print(ls_coords)
            norms = np.linalg.norm(ls_coords, axis=-1)
            sorting = np.argsort(norms)
            # print(norms[sorting[:5]])
            # print(ls_coords[sorting[:5]])
            # print(pts[sorting[:5]])


            # raise Exception(ls_coords)
            coords = self.normal_modes.modes.basis.to_new_modes().unembed_coords(
                np.reshape(ls_coords, ls_coords.shape[:2])
            )

            # print(coords)

        else:
            coords = np.broadcast_to(ref[np.newaxis], (len(pts),) + ref.shape).copy()
            coords[np.arange(len(atom_idx))[:, np.newaxis], atom_idx[:, np.newaxis], axes[np.newaxis, :]] = pts
            coords = coords.reshape(base_shape + coords.shape[-2:])
        if smol: coords = coords[0]

        if return_distances:
            return coords, dists
        else:
            return coords

    def get_nearest_scan_coordinates(self, domains, sel=None, axes=None):
        displacement_mesh = np.moveaxis(
            np.array(
                np.meshgrid(*[np.linspace(*d) for d in domains], indexing='ij')
            ),
            0, -1
        )
        return self.get_nearest_displacement_coordinates(displacement_mesh, axes=axes, sel=sel)

class PropertyEvaluator(metaclass=abc.ABCMeta):

    def __init__(self,
                 embedding=None,
                 permutation=None,
                 use_internals=True,
                 strip_embedding=True,
                 flatten_internals=True,
                 reembed_cartesians=False,
                 **defaults):
        self.defaults = defaults
        self.embedding = embedding
        if use_internals and not self.supports_internals:
            use_internals = 'reembed'
        self.use_internals = use_internals
        self.strip_embedding = strip_embedding
        self.flatten_internals = flatten_internals
        self.reembed_cartesians = reembed_cartesians
        self.permutation = permutation

    def use_internal_coordinate_handlers(self):
        return (self.use_internals and self.embedding.internals is not None)

    def embed_coords(self, coords, embed_reembedded=True):
        coords = np.asanyarray(coords)
        if self.embedding is not None:
            if (
                    self.use_internals
                    and self.embedding.internals is not None
                    and (embed_reembedded or not dev.str_is(self.use_internals, "reembed"))
            ):
                base_coords = coords
                coords = self.embedding.get_internals(coords=coords, strip_embedding=self.strip_embedding)
                if not self.strip_embedding and self.flatten_internals:
                    if base_coords.ndim == coords.ndim:
                        coords = coords.reshape(base_coords.shape[:-2] + (-1,))
            elif self.reembed_cartesians:
                coords = self.embedding.embed_coords(coords)
        if self.permutation is not None:
            if self.use_internals and self.embedding.internals is not None:
                if not self.strip_embedding and not self.flatten_internals:
                    raise NotImplementedError("don't know how to permute not-flat internals...")
                coords = coords[..., self.permutation]
            else:
                coords = coords[..., self.permutation, :]
        return coords

    def unembed_coords(self, coords, embed_reembedded=True):
        # coords = np.asanyarray(coords)
        if self.embedding is not None:
            if (
                    self.use_internals
                    and self.embedding.internals is not None
                    and (embed_reembedded or not dev.str_is(self.use_internals, "reembed"))
            ):
                coords = self.embedding.get_cartesians(coords=coords, strip_embedding=self.strip_embedding)

        if self.permutation is not None:
            inv = np.argsort(self.permutation)
            if self.use_internals and self.embedding.internals is not None:
                if not self.strip_embedding and not self.flatten_internals:
                    raise NotImplementedError("don't know how to permute not-flat internals...")
                coords = coords[..., inv]
            else:
                coords = coords[..., inv, :]
        return coords

    def embed_derivs(self, coords, derivs, embed_reembedded=True):
        if len(derivs) == 0: return derivs

        if self.embedding is not None:
            if (
                    self.use_internals
                    and self.embedding.internals is not None
                    and (embed_reembedded or not dev.str_is(self.use_internals, "reembed"))
            ):
                expansion = self.embedding.get_cartesians_by_internals(
                    order=len(derivs),
                    coords=coords,
                    strip_embedding=self.strip_embedding
                )
                derivs = nput.tensor_reexpand(expansion, derivs, len(derivs))

        return derivs

    def unembed_derivs(self, base_coords, coords, derivs, embed_reembedded=True):
        if len(derivs) == 0: return derivs

        if self.permutation is not None:
            inv = np.argsort(self.permutation)
            new_derivs = []
            base_shape = base_coords.shape[:-2]
            ndim = len(base_shape)
            f_shape = derivs[0].shape[ndim + 1:]
            fdim = len(f_shape)
            ncoord = derivs[0].shape[ndim]
            nat = ncoord // 3  # not always used, but never an error
            for i, d in enumerate(derivs):
                # TODO: fix assumption that we always request all orders
                if self.use_internals and self.embedding.internals is not None:
                    for ax in range(-(i + 1), 0):
                        d = np.take(d, inv, axis=(ax - fdim))
                else:
                    for ax in range(i + 1):
                        coord_shape = tuple(ncoord for _ in range(ax)) + (nat, 3) + tuple(ncoord for _ in range(i - ax))
                        d_shape = d.shape
                        d = d.reshape(base_shape + coord_shape + f_shape)
                        d = np.take(d, inv, axis=ndim + ax).reshape(d_shape)
                new_derivs.append(d)
            derivs = new_derivs

        if self.embedding is not None:
            if (
                    self.use_internals
                    and self.embedding.internals is not None
                    and (embed_reembedded or not dev.str_is(self.use_internals, "reembed"))
            ):
                expansion = self.embedding.get_internals_by_cartesians(
                    order=len(derivs),
                    coords=base_coords,
                    strip_embedding=self.strip_embedding
                )
                derivs = nput.tensor_reexpand(expansion, derivs, len(derivs))
            elif self.reembed_cartesians:
                derivs = self.embedding.unembed_derivs(coords, derivs)
        return derivs

    @abc.abstractmethod
    def evaluate_term(self, coords, order, **opts):
        ...

    @classmethod
    @abc.abstractmethod
    def from_mol(cls, mol, **opts):
        ...

    fd_defaults=dict(
        stencil=None,
        mesh_spacing=.005,
        displacement_function=None,
        prep=None,
        lazy=False,
        cache_evaluations=True,
        parallelizer=None,
    )
    def get_fd_opts(self, **opts):
        return dict(self.fd_defaults, **opts)

    def finite_difference_derivs(self, coords, order, batched_orders=None, **opts):
        opts = dev.OptionsSet(opts)
        fd_opts = opts.filter(FiniteDifferenceDerivative)
        opts = opts.exclude(FiniteDifferenceDerivative)

        if batched_orders is None:
            batched_orders = self.batched_orders

        base_shape = coords.shape[:-2]
        coord_shape = coords.shape[-2:]
        flat_coords = coords.reshape(-1, np.prod(coord_shape, dtype=int))

        def derivs(structs):
            reconst = structs.reshape(structs.shape[:-1] + coord_shape)
            ders = self.evaluate_term(reconst, self.analytic_derivative_order, **opts)
            if batched_orders:
                return ders[-1]
            else:
                return ders

        der = FiniteDifferenceDerivative(derivs,
                                         function_shape=((0,), (0,) * self.analytic_derivative_order),
                                         **self.get_fd_opts(**fd_opts)
                                         )
        tensors = der.derivatives(flat_coords).derivative_tensor(
            list(range(1, (order - self.analytic_derivative_order) + 1))
        )
        if flat_coords.shape[0] > 1:
            tensors = [np.moveaxis(d, n + 1, 0) for n, d in enumerate(tensors)]

        res = [
            t.reshape(base_shape + t.shape[1:])
                if t.ndim > 0 and flat_coords.shape[0] != 1 else
            t
            for t in tensors
        ]

        return res

    def internal_finite_difference_derivs(self, coords, order, batched_orders=None, **opts):
        opts = dev.OptionsSet(opts)
        fd_opts = opts.filter(FiniteDifferenceDerivative)
        opts = opts.exclude(FiniteDifferenceDerivative)

        if batched_orders is None:
            batched_orders = self.batched_orders

        base_shape = coords.shape[:-1]
        coord_shape = coords.shape[-1:]
        flat_coords = coords.reshape(-1, np.prod(coord_shape, dtype=int))
        reembed = dev.str_is(self.use_internals, 'reembed')
        def derivs(structs):
            reconst = structs.reshape(structs.shape[:-1] + coord_shape)
            if reembed:
                reconst = self.unembed_coords(reconst)
                if not batched_orders:
                    res = [
                        self.evaluate_term(reconst, o, **opts)
                        for o in range(self.analytic_derivative_order + 1)
                    ]
                else:
                    res = self.evaluate_term(reconst, self.analytic_derivative_order, **opts)
                return res[0] + self.embed_derivs(reconst, res[1:])
            else:
                ders = self.evaluate_term(reconst, self.analytic_derivative_order, **opts)
                if batched_orders:
                    return ders[-1]
                else:
                    return ders

        der = FiniteDifferenceDerivative(derivs,
                                         function_shape=((0,), (0,) * self.analytic_derivative_order),
                                         **self.get_fd_opts(**fd_opts)
                                         )
        tensors = der.derivatives(flat_coords).derivative_tensor(
            list(range(1, (order - self.analytic_derivative_order) + 1))
        )
        if flat_coords.shape[0] > 1:
            tensors = [np.moveaxis(d, n+1, 0) for n,d in enumerate(tensors)]
        res = [
            t.reshape(base_shape + t.shape[1:])
                if t.ndim > 0 and flat_coords.shape[0] != 1 else
            t
            for t in tensors
        ]

        return res

    supports_internals = False
    property_units = None
    target_property_units = None
    distance_units = 'Angstroms'
    batched_orders = False
    analytic_derivative_order = 1
    def evaluate_expansion(self,
                 coords,
                 order=0,
                 analytic_derivative_order=None,
                 batched_orders=None,
                 logger=None,
                 fd_handler=None,
                 **opts):

        opts = dev.OptionsSet(dict(self.defaults, **opts))
        fd_opts = opts.filter(FiniteDifferenceDerivative)
        opts = opts.exclude(FiniteDifferenceDerivative)

        if nput.is_numeric(order):
            order = list(range(order+1))
        # TODO: we assume order sorted, handle when it's not
        expansion = []
        if analytic_derivative_order is None:
            analytic_derivative_order = self.analytic_derivative_order
        if batched_orders is None:
            batched_orders = self.batched_orders
        if batched_orders:
            terms = self.evaluate_term(coords, min(analytic_derivative_order, order[-1]), **opts)
            expansion.extend([terms[o] for o in order if o <= analytic_derivative_order])
        else:
            for i in order:
                if i > analytic_derivative_order: break
                expansion.append(self.evaluate_term(coords, i, **opts))
        if order[-1] > analytic_derivative_order:
            if fd_handler is None: fd_handler = self.finite_difference_derivs

            # self.fd_derivs(**opts)
            #TODO: handle disjoin list of orders...
            expansion.extend(
                fd_handler(coords, order[-1], batched_orders=batched_orders, **fd_opts)
            )

        if self.property_units is None:
            eng_conv = 1
        else:
            eng_conv = UnitsData.convert(self.property_units, self.target_property_units)

        if self.distance_units is None:
            bohr_conv = 1
        else:
            bohr_conv = UnitsData.convert(self.distance_units, "BohrRadius")

        return [
            e * eng_conv / (bohr_conv ** o)
            for o,e in zip(order, expansion)
        ]

    def evaluate(self,
                 coords,
                 order=0,
                 logger=None,
                 **opts):
        use_internals = self.use_internals
        reembed = self.reembed_cartesians
        perm = self.permutation
        base_coords = coords
        coords = self.embed_coords(coords, embed_reembedded=False)

        try:
            self.use_internals = False
            self.reembed_cartesians = False
            self.permutation = None

            # TODO: make this less hacky...
            expansion = self.evaluate_expansion(
                coords,
                order=order,
                logger=logger,
                **opts
            )
        finally:
            self.use_internals = use_internals
            self.reembed_cartesians = reembed
            self.permutation = perm

        expansion = [expansion[0]] + self.unembed_derivs(base_coords, coords, expansion[1:], embed_reembedded=False)
        return expansion

    # evaluator_types = {}
    @classmethod
    @abc.abstractmethod
    def get_evaluators(cls):
        ...

    @classmethod
    def get_evaluators_by_attributes(cls):
        return None

    _profile_dispatch = dev.uninitialized
    default_evaluator_type = 'expansion'
    @classmethod
    def profile_generator_dispatch(cls):
        cls._profile_dispatch = dev.handle_uninitialized(
            cls._profile_dispatch,
            dev.OptionsMethodDispatch,
            args=(cls.get_evaluators,),
            kwargs=dict(
                default_method=cls.default_evaluator_type,
                attributes_map=cls.get_evaluators_by_attributes()
            )
        )
        return cls._profile_dispatch

    @classmethod
    @abc.abstractmethod
    def get_default_function_evaluator_type(cls):
        ...

    @classmethod
    def resolve_evaluator(cls, name):
        if name is not None and not isinstance(name, (str, dict, tuple)) and callable(name):
            if hasattr(name, 'evaluate_term'):
                return name, {}
            else:
                eval = cls.get_default_function_evaluator_type()
                eval.bind_default(name)
                return eval, {}
        return cls.profile_generator_dispatch().resolve(name)

    class quiet_mode:
        def __init__(self, quiet=True):
            self.quiet = quiet
            self._stdout = None
            self._devnull = None
        def __enter__(self):
            if self.quiet:
                self._stdout = sys.stdout
                self._devnull = open(os.devnull, 'w+').__enter__()
                sys.stdout = self._devnull
        def __exit__(self, exc_type, exc_val, exc_tb):
            if self.quiet:
                self._devnull.__exit__(exc_type, exc_val, exc_tb)
                sys.stdout = self._stdout

class PropertyFunctionEvaluator(PropertyEvaluator):
    # Slightly strange use case...we mostly use this as a component to
    # implement individual methods

    def __init__(self,
                 potential_function,
                 property_units=None,
                 distance_units='Angstroms',
                 batched_orders=False,
                 supports_internals=True,
                 analytic_derivative_order=0,
                 **defaults
                 ):
        super(type(self).__bases__[0], self).__init__(**defaults)
        self.property_function = potential_function
        self.batched_orders = batched_orders
        self.analytic_derivative_order = analytic_derivative_order
        self.property_units = property_units
        self.distance_units = distance_units
        self.supports_internals = supports_internals

    def evaluate_term(self, coords, order, **opts):
        return self.property_function(coords, order=order, **opts)

    default_property_function = None
    @classmethod
    def bind_default(cls, potential):
        cls.default_property_function = potential
    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        return prop_func
    @classmethod
    def initializer(cls, mol, **opts):
        return cls(None, **opts)
    @classmethod
    def initialize_from_mol(root,
                            cls,
                            mol,
                            property_function=None,
                            *,
                            embedding=None,
                            use_internals=None,
                            reembed_cartesians=None,
                            property_units=None,
                            distance_units=None,
                            batched_orders=None,
                            analytic_derivative_order=None,
                            **opts):
        if property_function is None:
            property_function = cls.default_property_function
            cls.default_property_function = None
        property_function, opts = cls.get_property_function(property_function, mol, **opts)
        if property_function is None:
            raise ValueError(f"can't evaluate with property function {property_function}")
        init_opts = {
            k: v for k, v in dict(
                embedding=mol.embedding if embedding is None else embedding,
                use_internals=use_internals,
                reembed_cartesians=reembed_cartesians,
                property_units=property_units,
                distance_units=distance_units,
                batched_orders=batched_orders,
                analytic_derivative_order=analytic_derivative_order,
            ).items()
            if v is not None
        }
        return cls(
            property_function,
            **opts,
            **init_opts
        )

class EnergyEvaluator(PropertyEvaluator):

    target_property_units = 'Hartrees'
    @classmethod
    def get_evaluators(cls):
        return {
            'rdkit': RDKitEnergyEvaluator,
            'aimnet2': AIMNet2EnergyEvaluator,
            'xtb': XTBEnergyEvaluator,
            'expansion': PotentialExpansionEnergyEvaluator
        }
    @classmethod
    def get_evaluators_by_attributes(cls):
        return {
            ('potential_function',):PotentialFunctionEnergyEvaluator
        }
    @classmethod
    def get_default_function_evaluator_type(cls):
        return PotentialFunctionEnergyEvaluator

    def minimizer_function_by_order(self, order, allow_fd=False, **opts):
        conv = UnitsData.convert(self.property_units, "Hartrees")
        if self.analytic_derivative_order >= order:
            def func(crd, _):
                res = self.evaluate_term(crd.reshape(1, -1, 3), order, **opts)
                if self.batched_orders:
                    return res[-1] * conv
                else:
                    return res * conv
        elif allow_fd:
            def func(crd, _):
                res = self.finite_difference_derivs(crd.reshape(1, -1, 3), order, **opts)[-1]
                return res[np.newaxis] * conv # it obliterates the 1
        else:
            func = None
        return func
    def minimizer_func(self, **opts):
        return self.minimizer_function_by_order(0, **opts)
    def minimizer_jacobian(self, **opts):
        return self.minimizer_function_by_order(1, allow_fd=True, **opts)
    def minimizer_hessian(self, **opts):
        return self.minimizer_function_by_order(2, **opts)

    def minimizer_internal_function_by_order(self, order, allow_fd=False, **opts):
        conv = UnitsData.convert(self.property_units, "Hartrees")
        if self.analytic_derivative_order >= order:
            opts = dev.OptionsSet(opts)
            # fd_opts = opts.filter(FiniteDifferenceDerivative)
            opts = opts.exclude(FiniteDifferenceDerivative)

            reembed = dev.str_is(self.use_internals, 'reembed')
            def func(crd, _):
                # res = self.evaluate_term(crd, order, **opts)
                if reembed:
                    crd = self.unembed_coords(crd)
                    if not self.batched_orders:
                        res = [
                            self.evaluate_term(crd, o, **opts)
                            for o in range(order + 1)
                        ]
                    else:
                        res = self.evaluate_term(crd, order, **opts)
                    dist_conv = UnitsData.convert("BohrRadius", self.distance_units)
                    res = [res[0]] + self.embed_derivs(crd, [r*dist_conv**(n+1) for n,r in enumerate(res[1:])])
                    if not self.batched_orders:
                        res = res[-1]
                else:
                    res = self.evaluate_term(crd, order, **opts)
                if self.batched_orders:
                    return res[-1] * conv
                else:
                    return res * conv
        elif allow_fd:
            def func(crd, _):
                res = self.internal_finite_difference_derivs(crd, order, **opts)[-1]
                return res[np.newaxis] * conv  # it obliterates the 1
        else:
            func = None
        return func

    def get_internal_coordinate_projector(self, coord_spec, mask=None):

        base_internals = self.embedding.internals
        if base_internals.get('specs') is not None:
            specs = base_internals.get('specs')
            sidx = [tuple(s) for s in specs]
            cpos = [tuple(c) for c in coord_spec]
            inds = [sidx.index(c) for c in cpos]
            num_coords = len(inds)
            redundant_transformation = base_internals.get('redundant_transformation')
            base_tensor = np.zeros((len(specs), num_coords))
            for n, i in enumerate(inds):
                base_tensor[..., i, n] = 1

            if redundant_transformation is not None:
                base_tensor = nput.maximum_similarity_transformation(
                    redundant_transformation,
                    base_tensor,
                    apply_transformation=True
                )

            proj = nput.orthogonal_projection_matrix(base_tensor)
            def constraint(coords):
                nstruct = len(coords.shape[:-1])
                return np.broadcast_to(
                    np.expand_dims(proj, list(range(nstruct))),
                    coords.shape[:-1] + proj.shape
                )

        elif base_internals.get('zmatrix') is not None:
            zmat = base_internals.get('zmatrix')
            num_specs = coordops.num_zmatrix_coords(
                zmat,
                coord_spec
            )
            inds = coordops.zmatrix_indices(
                zmat,
                coord_spec
            )
            num_coords = len(inds)

            base_tensor = np.zeros((num_specs, num_coords))
            for n, i in enumerate(inds):
                base_tensor[..., i, n] = 1
            proj = nput.orthogonal_projection_matrix(base_tensor)

            def constraint(coords):
                nstruct = len(coords.shape[:-1])
                return np.broadcast_to(
                    np.expand_dims(proj, list(range(nstruct))),
                    coords.shape[:-1] + proj.shape
                )
        else:
            raise ValueError(f"can't get {coord_spec} for {base_internals}")

        return constraint

    def get_internal_coordinate_constraints(self, constraints):
        if dev.is_dict_like(constraints):
            coord_list = list(constraints.keys())
            funs = list(constraints.values())

            def mask_func(internal_vals):
                return (f(i) for i, f in zip(internal_vals, funs))

            return self.get_internal_coordinate_projector(coord_list, mask=mask_func)
        elif coordops.is_coordinate_list_like(constraints):
            return self.get_internal_coordinate_projector(constraints)
        elif nput.is_numeric_array_like(constraints):
            return nput.orthogonal_projection_matrix(constraints)
        else:
            return constraints
    internal_optimizer_defaults = dict(
        max_displacement=.01
    )

    optimizer_defaults = dict(
        method='quasi-newton',
        unitary=False,
        # generate_rotation=False,
        # dtype='float64',
        orthogonal_directions=None,
        convergence_metric=None,
        tol=1e-8,
        max_iterations=25,
        damping_parameter=None,
        damping_exponent=None,
        restart_interval=None,
        max_displacement=.2,
        line_search=False,
        optimizer_settings=None,
        mode=None,
        func=None,
        jacobian=None,
        hessian=None,
        logger=None,

        generate_rotation=False,
        dtype='float64',
        orthogonal_projection_generator=None,
        region_constraints=None,
        max_displacement_norm=None,
        oscillation_damping_factor=None,
        termination_function=None,
        prevent_oscillations=None,
        use_max_for_error=True,
        track_best=True
    )
    @classmethod
    def get_optimizer_options(self):
        return (
                tuple(self.optimizer_defaults.keys())
                + ('coordinate_constraints',)
                + FiniteDifferenceDerivative.__props__
        )

    def get_coordinate_projector(self, coord_spec, mask=None):
        # tf_fun = nput.internal_conversion_function(coord_spec, order=1)
        def constraint(coords):
            coords = coords.reshape((-1, len(self.embedding.masses), 3))
            bases, _, _ = nput.internal_basis(coords, coord_spec)
            # base_tensors = tf_fun(coords)
            # if mask is not None:
            #     checks = mask(base_tensors[0])
            #     #TODO: handle the mixed-shape case
            #     base_tensors = base_tensors[1][..., :, checks]
            base_tensors = np.concatenate(bases, axis=1)

            # print(bases[0])
            proj = nput.orthogonal_projection_matrix(base_tensors, orthonormal=True)
            # print(proj)
            # raise Exception(...)

            return proj

        return constraint

    def get_coordinate_constraints(self, constraints):
        if dev.is_dict_like(constraints):
            coord_list = list(constraints.keys())
            funs = list(constraints.values())
            def mask_func(internal_vals):
                return (f(i) for i,f in zip(internal_vals, funs))

            return self.get_coordinate_projector(coord_list, mask=mask_func)
        elif coordops.is_coordinate_list_like(constraints):
            return self.get_coordinate_projector(constraints)
        elif nput.is_numeric_array_like(constraints):
            return nput.orthogonal_projection_matrix(constraints)
        else:
            return constraints

    scipy_no_hessian_methods = {'bfgs'}
    def optimize_iterative(self,
                 coords,
                 coordinate_constraints=None,
                 **opts
                 ):
        from McUtils.Numputils import iterative_step_minimize
        opts = dict(self.optimizer_defaults, **{o:k for o,k in opts.items() if k is not None})

        (
            method,
            unitary,
            # generate_rotation=False,
            # dtype='float64',
            orthogonal_directions,
            orthogonal_projection_generator,
            convergence_metric,
            tol,
            max_iterations,
            damping_parameter,
            damping_exponent,
            restart_interval,
            max_displacement,
            line_search,
            optimizer_settings,
            mode,
            func,
            jacobian,
            hessian,
            logger
        ) = (
            opts.pop(k) for k in
            [
                "method",
                "unitary",
                # generate_rotation=False,
                # dtype='float64',
                "orthogonal_directions",
                "orthogonal_projection_generator",
                "convergence_metric",
                "tol",
                "max_iterations",
                "damping_parameter",
                "damping_exponent",
                "restart_interval",
                "max_displacement",
                "line_search",
                "optimizer_settings",
                "mode",
                "func",
                "jacobian",
                "hessian",
                "logger"
            ]
        )

        fopts = dev.OptionsSet(opts).exclude(None, props=self.optimizer_defaults.keys())
        opt_opts = dev.OptionsSet(opts).filter(None, props=self.optimizer_defaults.keys())

        if optimizer_settings is None:
            optimizer_settings = {}

        #TODO: let evaluations cache results when batched
        if func is None:
            func = self.minimizer_func(**fopts)
        if jacobian is None:
            jacobian = self.minimizer_jacobian(**fopts)
        if hessian is None:
            hessian = self.minimizer_hessian(**fopts)

        logger = Logger.lookup(logger, construct=True)

        if mode == 'scipy' or method == 'scipy':
            from scipy.optimize import minimize

            def sfunc(crd):
                return func(crd, None)[0]
            def sjacobian(crd):
                return jacobian(crd, None)[0]

            if self.analytic_derivative_order > 1:
                def shessian(crd):
                    return hessian(crd, None)[0]
            else:
                shessian = None

            callback = optimizer_settings.pop('callback', None)
            if logger.active:
                if callback is None:
                    callback = lambda intermediate_result: logger.log_print(
                        "Step: {intermediate_result}",
                        intermediate_result=intermediate_result
                    )
                else:
                    callback = lambda intermediate_result, cb=callback: [
                        logger.log_print(
                            "Step: {intermediate_result}",
                            intermediate_result=intermediate_result
                        ),
                        cb(intermediate_result)
                    ][-1]

            method = (
                method
                    if dev.str_is(mode, 'scipy') else
                opts.pop('scipy_method', self.optimizer_defaults.get('method', 'bfgs'))
            )
            method = 'bfgs' if method=='quasi-newton' else method
            if method in self.scipy_no_hessian_methods:
                shessian = None
            min = minimize(sfunc,
                           coords.flatten(),
                           method='bfgs' if method=='quasi-newton' else method,
                           tol=tol,
                           jac=sjacobian,
                           hess=shessian,
                           callback=callback,
                           options=dict({'maxiter':max_iterations}, **optimizer_settings)
                           )
            return True, min.x.reshape(coords.shape), min
        else:
            if method is None or isinstance(method, str):
                method = dict(
                    {
                        'method': method,
                        'func': func,
                        'jacobian': jacobian,
                        'hessian': hessian,
                        'damping_parameter': damping_parameter,
                        'damping_exponent': damping_exponent,
                        'restart_interval': restart_interval,
                        'line_search': line_search,
                    },
                    **optimizer_settings
                )
                method = {k:v for k,v in method.items() if v is not None}
            if orthogonal_projection_generator is None:
                orthogonal_projection_generator = self.get_coordinate_constraints(coordinate_constraints)
            opt_coords, converged, (errs, its) = iterative_step_minimize(
                coords.flatten(),
                method,
                unitary=unitary,
                orthogonal_directions=orthogonal_directions,
                convergence_metric=convergence_metric,
                tol=tol,
                max_iterations=max_iterations,
                max_displacement=max_displacement,
                logger=logger,
                orthogonal_projection_generator=orthogonal_projection_generator,
                **opt_opts
            )
            return converged, opt_coords.reshape(coords.shape), {'errors':errs, 'iterations':its}


    def optimize(self,
                 coords,
                 use_internals=None,
                 func=None,
                 jacobian=None,
                 hessian=None,
                 logger=None,
                 coordinate_constraints=None,
                 orthogonal_projection_generator=None,
                 **opts
                 ):
        if use_internals is None: use_internals = self.use_internals
        if use_internals and self.embedding.internals is not None:
            opts = dict(self.internal_optimizer_defaults, **opts)
            fopts = dev.OptionsSet(opts).exclude(None, props=self.optimizer_defaults.keys())
            if func is None:
                func = self.minimizer_internal_function_by_order(0, **fopts)
            if jacobian is None:
                jacobian = self.minimizer_internal_function_by_order(1, allow_fd=True, **fopts)
            if hessian is None:
                hessian = self.minimizer_internal_function_by_order(2, **fopts)
            try:
                self.use_internals = use_internals
                coords = self.embed_coords(coords)
                if orthogonal_projection_generator is None:
                    orthogonal_projection_generator = self.get_internal_coordinate_constraints(coordinate_constraints)
                convergence, opt_coords, opt_settings = self.optimize_iterative(coords,
                                                                         func=func,
                                                                         jacobian=jacobian,
                                                                         hessian=hessian,
                                                                         logger=logger,
                                                                         orthogonal_projection_generator=orthogonal_projection_generator,
                                                                         **opts
                                                                         )
                coords = self.unembed_coords(opt_coords)
                return convergence, coords, opt_settings
            finally:
                self.use_internals = use_internals

        else:
            return self.optimize_iterative(coords,
                                    logger=logger,
                                    coordinate_constraints=coordinate_constraints,
                                    orthogonal_projection_generator=orthogonal_projection_generator,
                                    **opts
                                    )

class RDKitEnergyEvaluator(EnergyEvaluator):
    def __init__(self, rdmol, force_field='mmff', **defaults):
        super().__init__(**defaults)
        # if hasattr(mol, 'rdmol'):
        #     mol = mol.rdmol
        self.rdmol = rdmol
        self.force_field = force_field

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.rdmol, embedding=mol.embedding, **opts)

    property_units = 'Kilocalories/Mole'
    analytic_derivative_order = 1
    def evaluate_term(self, coords, order, force_field_type=None, **opts):
        if force_field_type is None:
            force_field_type = self.force_field
        if order == 0:
            return self.rdmol.calculate_energy(coords, force_field_type=force_field_type, **opts)
        elif order == 1:
            return self.rdmol.calculate_gradient(coords, force_field_type=force_field_type, **opts)
        else:
            raise ValueError(f"order {order} not supported")

    def optimize(self,
                 coords,
                 method=None,
                 force_field_type=None,
                 unitary=False,
                 # generate_rotation=False,
                 # dtype='float64',
                 orthogonal_directions=None,
                 convergence_metric=None,
                 tol=1e-8,
                 max_iterations=100,
                 damping_parameter=None,
                 damping_exponent=None,
                 reset_interval=None,
                 **opts
                 ):
        if force_field_type is None:
            force_field_type = self.force_field

        if method == 'rdkit':
            maxIters = opts.pop('maxIters', max_iterations)
            optimizer = opts.pop('optimizer', method)
            return self.rdmol.optimize_structure(
                geoms=coords,
                force_field_type=force_field_type,
                optimizer=optimizer,
                maxIters=maxIters,
                **opts
            )
        else:
            return super().optimize(
                coords,
                method=method,
                tol=tol,
                max_iterations=max_iterations,
                **opts
            )

class AIMNet2EnergyEvaluator(EnergyEvaluator):
    """
    Borrows structure from AIMNet2ASE to call appropriately
    """
    def __init__(self, atoms, model='aimnet2', charge=0, multiplicity=None, quiet=True, **defaults):
        super().__init__(**defaults)
        self.eval = self.setup_aimnet(model)
        self.model = model
        self.atoms = atoms
        self.numbers = [AtomData[atom, "Number"] for atom in atoms]
        self.charge = charge
        self.multiplicity = multiplicity
        self.quiet = quiet

    @classmethod
    def from_mol(cls, mol, charge=None, multiplicity=None, **opts):
        return cls(mol.atoms, embedding=mol.embedding, multiplicity=multiplicity, **opts)

    @classmethod
    def setup_aimnet(cls, model):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            from aimnet2calc import AIMNet2Calculator

        with cls.quiet_mode():
            calc = AIMNet2Calculator(model)

        return calc

    @staticmethod
    def autodiff(expr, coord, create_graph=False, retain_graph=False):
        import torch
        # here forces have shape (N, 3) and coord has shape (N+1, 3)
        # return hessian with shape (N, 3, N, 3)
        shape = coord.shape + expr.shape
        grads = []
        # N = coord.shape[0] - 1
        # for p in itertools.combinations(coord.shape, expr.ndim):
        rav_shape = tuple(
            (expr.shape[2*i]*expr.shape[2*i+1])
            for i in range(expr.ndim // 2)
        )
        cache = {}
        inds = np.array(np.unravel_index(np.arange(np.prod(expr.shape, dtype=int)), rav_shape)).T
        shinds = np.sort(inds, axis=1)
        for i,_f in enumerate(expr.flatten()):
            ravioli = tuple(inds[i])
            key = tuple(shinds[i])
            if ravioli == key:
                g = torch.autograd.grad(_f,
                                        coord,
                                        # allow_unused=True,
                                        retain_graph=True,
                                        create_graph=create_graph)
                g = g[0]
                grads.append(g)
                cache[key] = g
            else:
                grads.append(cache[key])

        deriv = torch.stack(grads).view(shape)
        shape_tuple = (slice(None, -1), slice(None)) * (len(shape) // 2)
        return deriv, shape_tuple

    # def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
    #     super().calculate(atoms, properties, system_changes)
    #     self.update_tensors()
    #
    #     if self.atoms.cell is not None and self.atoms.pbc.any():
    #         # assert self.base_calc.cutoff_lr < float('inf'), 'Long-range cutoff must be finite for PBC'
    #         cell = self.atoms.cell.array
    #     else:
    #         cell = None
    #
    #     results = self.base_calc({
    #         'coord': torch.tensor(self.atoms.positions, dtype=torch.float32, device=self.base_calc.device),
    #         'numbers': self._t_numbers,
    #         'cell': cell,
    #         'mol_idx': self._t_mol_idx,
    #         'charge': self._t_charge,
    #         'mult': self._t_mult,
    #     }, forces='forces' in properties, stress='stress' in properties)

    property_units = 'ElectronVolts'
    batched_orders = True
    analytic_derivative_order = 6
    def prep_eval(self, coords, order, **opts):
        import torch

        base_shape = coords.shape[:-2]
        coords = coords.reshape((-1,) + coords.shape[-2:])
        # numbers = np.broadcast_to(
        #     np.array(self.numbers)[np.newaxis],
        #     (len(coords), len(self.numbers))
        # )

        forces = order > 0
        hessian = order > 1

        arg_dict = {
            'coord': torch.tensor(coords.reshape(-1, 3), dtype=torch.float32, device=self.eval.device),
            'numbers': torch.tensor(
                np.repeat(np.array([self.numbers]), coords.shape[0], axis=0).flatten(),
                dtype=torch.int64, device=self.eval.device
            ),
            'charge': torch.tensor(self.charge, dtype=torch.float32, device=self.eval.device),
            'cell': None,
            'mol_idx': torch.tensor(
                np.repeat(np.arange(coords.shape[0]), coords.shape[1]).flatten(),
                dtype=torch.int64,
                device=self.eval.device
            )
        }
        if self.multiplicity is not None:
            arg_dict['mult'] = torch.tensor(self.multiplicity, dtype=torch.float32, device=self.eval.device)

        data = self.eval.prepare_input(arg_dict)
        if hessian and data['mol_idx'][-1] > 0:
            raise NotImplementedError('higher-derivative calculation is not supported for multiple molecules')
        data = self.eval.set_grad_tensors(data,
                                          forces=forces,
                                          stress=False,
                                          hessian=hessian
                                          )
        with torch.jit.optimized_execution(False):
            data = self.eval.model(data)
        data = self.eval.get_derivatives(
            data,
            forces=forces,
            hessian=hessian,
            **opts
        )

        return base_shape, coords, data

    def process_aimnet_derivs(self, base_shape, coords, data, root_key, order, deriv_key=None, property_shape=()):
        if order > 0:
            if deriv_key is not None:
                data['deriv_1'] = (-data[deriv_key], (slice(None,), slice(None,)))
                for o in range(2, order+1):
                    data[f'deriv_{o}'] = self.autodiff(data[f'deriv_{o-1}'][0],
                                                       self.eval._saved_for_grad['coord'],
                                                       # create_graph=(o==1),
                                                       create_graph=(o < order),
                                                       retain_graph=(o < order)
                                                       )
            else:
                data['deriv_0'] = (data[root_key], slice(None,))
                for o in range(1, order+1):
                    data[f'deriv_{o}'] = self.autodiff(data[f'deriv_{o - 1}'][0],
                                                       self.eval._saved_for_grad['coord'],
                                                       create_graph=(o==1),
                                                       retain_graph=(o < order)
                                                       )

        deriv_keys = [
            f'deriv_{o}'
            for o in range(1, order + 1)
        ]
        for key in deriv_keys:
            deriv,subshape = data[key]
            data[key] = deriv.__getitem__(subshape)
        ko = self.eval.keys_out
        afk = self.eval.atom_feature_keys
        try:
            self.eval.keys_out = self.eval.keys_out + deriv_keys
            self.eval.atom_feature_keys = self.eval.atom_feature_keys + deriv_keys[:1]
            expansion = self.eval.process_output(data)
        finally:
            self.eval.keys_out = ko
            self.eval.atom_feature_keys = afk

        ndim = np.prod(coords.shape[-2:], dtype=int)

        expansion = [
            expansion[k].detach().cpu().numpy().reshape(
                base_shape + (ndim,) * o + property_shape
            )
            for o,k in enumerate([root_key] + deriv_keys)
        ]

        return expansion
    def evaluate_term(self, coords, order, **opts):

        base_shape, coords, data = self.prep_eval(coords, order, **opts)

        return self.process_aimnet_derivs(
            base_shape, coords, data,
            'energy',
            order,
            deriv_key='forces'
        )

    def optimize(self,
                 coords,
                 method=None,
                 **opts
                 ):
        tol, max_iterations, max_displacement = [
            opts.pop(k, self.optimizer_defaults[k])
            for k in ["tol", "max_iterations", "max_displacement"]
        ]
        if method == 'ase':
            from McUtils.ExternalPrograms import ASEMolecule
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                from aimnet2calc import AIMNet2ASE

            with self.quiet_mode():
                calc = AIMNet2ASE(self.model,
                                  charge=self.charge,
                                  mult=self.multiplicity
                                  )

            mol = ASEMolecule.from_coords(
                self.atoms,
                coords.reshape((-1,) + coords.shape[-2:])[0],
                calculator=calc
            )
            return mol.optimize_structure(coords,
                                          fmax=tol, steps=max_iterations,
                                          maxstep=max_displacement
                                          )
        else:
            return super().optimize(
                coords,
                method=method,
                tol=tol,
                max_iterations=max_iterations,
                max_displacement=max_displacement,
                **opts
            )

class XTBEnergyEvaluator(EnergyEvaluator):
    """
    Uses XTB to calculate, not tested since `xtb` is terrible to install without `conda`
    """
    def __init__(self, atoms, method="GFN2-xTB", charge=0, quiet=True, **defaults):
        super().__init__(**defaults)

        # self.eval = self.setup_aimnet(model)
        # self.model = model
        self.atoms = atoms
        self.numbers = np.array([AtomData[atom, "Number"] for atom in atoms])
        # self.charge = charge
        # self.multiplicity = multiplicity
        self.quiet = quiet
        self.method = method
        self.charge = charge

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.atoms, embedding=mol.embedding, charge=mol.charge, **opts)

    def get_single_point(self, coords):
        from xtb.interface import Calculator
        from xtb.utils import get_method
        calc = Calculator(get_method(self.method), self.numbers, coords,
                          charge=self.charge)
        if self.quiet:
            calc.set_verbosity('muted')
        return calc.singlepoint()


    distance_units = "BohrRadius"
    property_units = 'Hartrees'
    batched_orders = True
    analytic_derivative_order = 1
    def evaluate_term(self, coords, order, **opts):
        # import torch

        base_shape = coords.shape[:-2]
        coords = coords.reshape((-1,) + coords.shape[-2:])

        if order > 1:
            raise NotImplementedError('`xtb` only implements up through gradients')

        energies = np.empty(coords.shape[0], dtype=float)
        if order > 0:
            grads = np.empty(coords.shape, dtype=float)
        else:
            grads = None
        for i,coord in enumerate(coords):
            res = self.get_single_point(coord)
            energies[i] = res.get_energy()
            if order > 0:
                grads[i] = res.get_gradient()

        energies = energies.reshape(base_shape)
        if order > 0:
            grads = grads.reshape(base_shape +(-1,))
            return [energies, grads]
        else:
            return [energies]

        # if order > 2:
        # if order > (forces=False, stress=False, hessian=False)
        # if order == 0:
        #     return self.rdmol.calculate_energy(coords, **opts)
        # elif order == 2:
        #     return self.rdmol.calculate_gradient(coords, **opts)
        # else:
        #     raise ValueError(f"order {order} not supported")

        # @staticmethod
        # def calculate_hessian(forces, coord):
        #     # here forces have shape (N, 3) and coord has shape (N+1, 3)
        #     # return hessian with shape (N, 3, N, 3)
        #     hessian = - torch.stack([
        #         torch.autograd.grad(_f, coord, retain_graph=True)[0]
        #         for _f in forces.flatten().unbind()
        #     ]).view(-1, 3, coord.shape[0], 3)[:-1, :, :-1, :]
        #     return hessian

    # def optimize(self,
    #              coords,
    #              method='quasi-newton',
    #              tol=1e-8,
    #              max_iterations=25,
    #              damping_parameter=None,
    #              damping_exponent=None,
    #              restart_interval=None,
    #              max_displacement=None,
    #              **opts
    #              ):
    #     if method == 'xtb-base':
    #         ...
    #     else:
    #         return super().optimize(
    #             coords,
    #             method=method,
    #             tol=tol,
    #             max_iterations=max_iterations,
    #             damping_parameter=damping_parameter,
    #             damping_exponent=damping_exponent,
    #             restart_interval=restart_interval,
    #             **opts
    #         )

class PotentialFunctionEnergyEvaluator(EnergyEvaluator):
    # slightly weird, want single inheritance, but using `PropertyFunctionEvaluator` to define the
    # interface--I'm sure there's a better way to do this, but multiple inheritance always feels
    # fragile and python doesn't implement traits...

    def __init__(self,
                 potential_function,
                 property_units=None,
                 energy_units='Hartrees',
                 distance_units='BohrRadius',
                 **opts
                 ):
        if property_units is None:
            property_units = energy_units
        PropertyFunctionEvaluator.__init__(
            self,
            potential_function,
            property_units=property_units,
            distance_units=distance_units,
            **opts
        )

    def evaluate_term(self, coords, order, **opts):
        return PropertyFunctionEvaluator.evaluate_term(
            self,
            coords,
            order,
            **opts
        )

    def use_internal_coordinate_handlers(self):
        return PropertyFunctionEvaluator.use_internal_coordinate_handlers(self)

    def evaluate(self,
                 coords,
                 order=0,
                 logger=None,
                 fd_handler=None,
                 **opts):
        if fd_handler is None and self.use_internal_coordinate_handlers():
            fd_handler = self.internal_finite_difference_derivs
        return PropertyFunctionEvaluator.evaluate(
            self,
            coords,
            order,
            fd_handler=fd_handler,
            **opts
        )

    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        return prop_func, opts

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 potential_function=None,
                 **opts
                 ):
        if property_function is None:
            property_function = potential_function
        return PropertyFunctionEvaluator.initialize_from_mol(
            cls,
            mol,
            property_function=property_function,
            **opts
        )

    default_property_function = None
    @classmethod
    def bind_default(cls, potential):
        cls.default_property_function = potential


class PotentialExpansionEnergyEvaluator(PotentialFunctionEnergyEvaluator):

    def __init__(self,
                 expansion,
                 center=None,
                 ref=None,
                 batched_orders=None,
                 transforms=None,
                 transformed_derivatives=False,
                 **opts
                 ):
        if not callable(expansion):
            expansion = PotentialSurface.from_derivatives(expansion,
                                                          center=center, ref=ref,
                                                          transforms=transforms,
                                                          transformed_derivatives=transformed_derivatives
                                                          )
        super().__init__(
            expansion,
            batched_orders=True,
            **opts
        )

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 expansion=None,
                 **opts
                 ):
        return super().from_mol(
            mol,
            property_function=expansion if property_function is None else property_function,
            **opts
        )

    @classmethod
    def get_property_function(cls, expansion, mol, transforms=None, **ignored):
        if not callable(expansion):
            transforms = transforms
            transformed_derivatives = False
            if expansion is None:
                expansion = mol.potential_derivatives
                # handle partial quatics
            if len(expansion) > 2 and expansion[2].shape[0] < expansion[1].shape[0]:
                if transforms is None:
                    modes = mol.get_normal_modes(use_internals=False, project_transrot=False).remove_mass_weighting()
                    transforms = [[modes.coords_by_modes], [modes.modes_by_coords]]
                if nput.is_numeric_array_like(transforms[0]):
                    tf = np.asanyarray(transforms[0])
                    if tf.ndim > 2:
                        tf = tf[0]
                else:
                    tf = np.asanyarray(transforms[0][0])
                _ = []
                for i,e in enumerate(expansion):
                    for j in range(min([i+1, 2])):
                        e = np.tensordot(tf, e, axes=[1, -1])
                    _.append(e)
                expansion = _
                transformed_derivatives = True
            expansion = PotentialSurface.from_mol(mol,
                                                  expansion=expansion,
                                                  transforms=transforms,
                                                  transformed_derivatives=transformed_derivatives
                                                  )
        return expansion, ignored

class DipoleEvaluator(PropertyEvaluator):
    target_property_units = ("ElementaryCharge", "BohrRadius")
    @classmethod
    def get_evaluators(cls):
        return {
                'expansion':DipoleExpansionEnergyEvaluator
            }
    @classmethod
    def get_default_function_evaluator_type(cls):
        return DipoleFunctionDipoleEvaluator

class DipoleFunctionDipoleEvaluator(DipoleEvaluator):

    def __init__(self,
                 potential_function,
                 property_units=('ElementaryCharge', 'BohrRadius'),
                 distance_units='BohrRadius',
                 **opts
                 ):
        PropertyFunctionEvaluator.__init__(
            self,
            potential_function,
            property_units=property_units,
            distance_units=distance_units,
            **opts
        )

    def embed_coords(self, coords):
        return PropertyFunctionEvaluator.embed_coords(
            self,
            coords
        )

    def unembed_coords(self, coords):
        return PropertyFunctionEvaluator.unembed_coords(
            self,
            coords
        )

    def unembed_derivs(self, base_coords, coords, derivs):
        return PropertyFunctionEvaluator.unembed_derivs(
            self,
            base_coords,
            coords,
            derivs
        )

    def evaluate_term(self, coords, order, **opts):
        return PropertyFunctionEvaluator.evaluate_term(
            self,
            coords,
            order,
            **opts
        )

    def evaluate(self,
                 coords,
                 order=0,
                 logger=None,
                 **opts):
        return PropertyFunctionEvaluator.evaluate(
            self,
            coords,
            order,
            **opts
        )

    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        return prop_func, opts

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 dipole_function=None,
                 **opts
                 ):
        if property_function is None:
            property_function = dipole_function
        return PropertyFunctionEvaluator.initialize_from_mol(
            cls,
            mol,
            property_function=property_function,
            **opts
        )

    default_property_function = None
    @classmethod
    def bind_default(cls, potential):
        cls.default_property_function = potential

class DipoleExpansionEnergyEvaluator(DipoleFunctionDipoleEvaluator):

    def __init__(self,
                 expansion,
                 center=None,
                 ref=None,
                 property_units=('ElementaryCharge', "BohrRadius"),
                 distance_units='BohrRadius',
                 analytic_derivative_order=None,
                 batched_orders=None,
                 **opts
                 ):
        if not callable(expansion):
            expansion = DipoleSurface.from_derivatives(expansion, center=center, ref=ref)
        if analytic_derivative_order is None:
            analytic_derivative_order = len(expansion.expansion_tensors)
        super().__init__(
            expansion,
            property_units=property_units,
            distance_units=distance_units,
            analytic_derivative_order=analytic_derivative_order,
            batched_orders=True,
            **opts
        )

    @classmethod
    def expansion_from_mol_charges(cls, mol, charges=None):
        if charges is None:
            charges = mol.calculate_charges()
        charges = np.asanyarray(charges)
        com = mol.center_of_mass
        coords = mol.coords
        disps = coords - com[np.newaxis]
        dip_contribs = charges[:, np.newaxis] * disps
        deriv = np.zeros(dip_contribs.shape + (3,), dtype=dip_contribs.dtype)
        for i in range(3):
            deriv[:, i, i] = dip_contribs[:, i]
        return [np.sum(dip_contribs, axis=0), deriv.reshape(-1, 3)]

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 expansion=None,
                 **opts
                 ):
        return super().from_mol(
            mol,
            property_function=expansion
            if property_function is None else property_function,
            **opts
        )
    @classmethod
    def get_property_function(cls, expansion, mol, center=None, transforms=None, **ignored):
        if not callable(expansion):
            transformed_derivatives = False
            if expansion is None:
                expansion = mol.dipole_derivatives
                if expansion is None:
                    expansion = cls.expansion_from_mol_charges(mol)
                # handle partial quatics

            dip_exp = mol.hamiltonian.dipole_expansion(expansion=expansion)
            expansion = dip_exp.get_terms(transformation=transforms)
            transformed_derivatives = True
            if transforms is None:
                transforms = [
                    [dip_exp.embedding.get_coords_by_modes(mass_weighted=False)],
                    [dip_exp.embedding.get_modes_by_coords(mass_weighted=False)]
                ]
                # print("FWD:", [f.shape for f in transforms[0]])
                # print("REV:", [f.shape for f in transforms[1]])


            # if len(expansion) > 2 and expansion[2].shape[0] < expansion[1].shape[0]:
            #     raise Exception("?")
            #     if transforms is None:
            #         modes = mol.get_normal_modes(use_internals=False, project_transrot=False).remove_mass_weighting()
            #         transforms = [[modes.coords_by_modes], [modes.modes_by_coords]]
            #     if nput.is_numeric_array_like(transforms[0]):
            #         tf = np.asanyarray(transforms[0])
            #         if tf.ndim > 2:
            #             tf = tf[0]
            #     else:
            #         tf = np.asanyarray(transforms[0][0])
            #     _ = []
            #     for i, e in enumerate(expansion[1:]):
            #         e = np.tensordot(tf, e, axes=[1, -2])
            #         _.append(e)
            #     expansion = [expansion[0]] + _
            #     transformed_derivatives = True

            expansion = DipoleSurface.from_mol(mol,
                                               center=center,
                                               expansion=expansion,
                                               transforms=transforms,
                                               transformed_derivatives=transformed_derivatives
                                               )
        return expansion, ignored


class ChargeEvaluator(PropertyEvaluator):
    target_property_units = "ElementaryCharge"
    default_evaluator_type = 'rdkit'
    @classmethod
    def get_evaluators(cls):
        return {
                'rdkit':RDKitChargeEvaluator,
                'aimnet2':AIMNet2ChargeEvaluator
            }
    @classmethod
    def get_default_function_evaluator_type(cls):
        return ChargeFunctionChargeEvaluator

class ChargeFunctionChargeEvaluator(ChargeEvaluator):
    def __init__(self,
                 potential_function,
                 property_units='ElementaryCharge',
                 distance_units='BohrRadius',
                 **opts
                 ):
        PropertyFunctionEvaluator.__init__(
            self,
            potential_function,
            property_units=property_units,
            distance_units=distance_units,
            **opts
        )

    def evaluate_term(self, coords, order, **opts):
        return PropertyFunctionEvaluator.evaluate_term(
            self,
            coords,
            order,
            **opts
        )

    def evaluate(self,
                 coords,
                 order=0,
                 logger=None,
                 **opts):
        return PropertyFunctionEvaluator.evaluate(
            self,
            coords,
            order,
            **opts
        )

    @classmethod
    def from_mol(cls,
                 mol,
                 property_function=None,
                 **opts
                 ):
        return PropertyFunctionEvaluator.initialize_from_mol(
            cls,
            mol,
            property_function=property_function,
            **opts
        )

    default_property_function = None
    @classmethod
    def bind_default(cls, potential):
        cls.default_property_function = potential

class RDKitChargeEvaluator(ChargeEvaluator):
    def __init__(self, rdmol, model='gasteiger', **defaults):
        super().__init__(**defaults)
        self.rdmol = rdmol
        self.model = model

    @classmethod
    def from_mol(cls, mol, **opts):
        return cls(mol.rdmol)

    property_units = 'ElementaryCharge'
    analytic_derivative_order = 0
    def evaluate_term(self, coords, order, model=None, **opts):
        if model is None:
            model = self.model
        return np.array(
            self.rdmol.evaluate_charges(coords, model=model, **opts)
        ) / np.sqrt(2) # TODO: I think my base dipoles are just scaled incorrectly...


class AIMNet2ChargeEvaluator(ChargeEvaluator):

    def __init__(self, atoms, model='aimnet2', charge=0, multiplicity=None, quiet=True, **defaults):
        super().__init__(**defaults)
        self.base = AIMNet2EnergyEvaluator(atoms, model=model, charge=charge, multiplicity=multiplicity, quiet=quiet)

    @classmethod
    def from_mol(cls, mol, charge=None, multiplicity=None, **opts):
        return cls(mol.atoms, multiplicity=multiplicity, **opts)

    property_units = 'ElementaryCharge'
    distance_units = 'Angstroms'
    batched_orders = True
    analytic_derivative_order = 6
    scaling_factor = 0.59667 # necessary to make water work...but I don't know what the actual units are
    def evaluate_term(self, coords, order, **opts):
        base_shape, coords, data = self.base.prep_eval(coords, order, **opts)

        expansions = self.base.process_aimnet_derivs(
            base_shape, coords, data,
            'charges',
            order,
            property_shape=(coords.shape[-2],)
        )
        s = self.scaling_factor / np.sqrt(2)
        expansions = [e * s for e in expansions]

        return expansions

class ReducedDimensionalPotentialHandler:

    def __init__(self, mol):
        self.mol = mol

    # @property
    # def cartesian_derivatives(self):
    #     ...

    @classmethod
    def get_potential_params(cls, spec, re, g, local_derivs,
                             quartic_potential_cutoff=5e-2,
                             poly_expansion_order=2,
                             **opts):
        # if (use_morse and len(local_derivs) > 3) or (use_morse is None and len(spec) == 2):
        f2, f3, f4 = [local_derivs[i + 1] if i < len(local_derivs) else 0 for i in range(3)]
        if np.abs(f3 / f4) < quartic_potential_cutoff:
            opts['method'] = 'poly'
            opts['g'] = g
            sg = np.sqrt(g)
            opts['w_coeffs'] = [
                         sg * np.sign(d) * np.power(np.abs(d), 1 / (k + 2))
                         # * (.25 ** (k + 2))
                         for k, d in enumerate(local_derivs[1:4])
                     ][:poly_expansion_order - 1]
        else:
            opts['method'] = 'morse'
            opts['g'] = g
            sg = np.sqrt(g)
            w = sg * np.sqrt(f2)
            wx = -((g / (4 * w)) ** 2) * (f4 - 5 / 3 * (f3**2) / f2)
            params = [w, wx]
            # if len(local_derivs) > 4:
            #     params = params + [
            #         sg * np.sign(d) * np.power(np.abs(d) / math.factorial(k + 4), 1 / (k + 4))  # * (.25 ** (k + 2))
            #         for k, d in enumerate(local_derivs[4:])
            #     ]
            opts['w'] = w
            opts['wx'] = wx
        opts['re'] = re

        return opts

    @classmethod
    def get_morse_potential(cls, spec, *, w, wx, g, re):
        return CoordinateFunction.morse(
                spec,
                re=re,
                w=w,
                wx=wx,
                g=g
            )
    @classmethod
    def get_poly_potential(cls, spec, *, w_coeffs=None, coeffs=None, g, re):
        if w_coeffs is not None:
            sg = np.sqrt(g)
            coeffs = [0] + [np.sign(c) * (np.abs(c) / sg) ** (k + 2) / math.factorial(k + 2) for k, c in enumerate(w_coeffs)]
        return CoordinateFunction.polynomial(
            spec,
            center=re,
            coeffs=coeffs,
            ref=0
        )
    @classmethod
    def get_potential_types(cls):
        return {
            'morse':cls.get_morse_potential,
            'poly':cls.get_poly_potential
        }
    @classmethod
    def get_potentials_by_attributes(cls):
        return {
            ('w', 'wx'): 'morse',
            ('coeffs',): 'poly',
            ('w_coeffs',): 'poly'
        }

    _pot_dispatch = dev.uninitialized
    default_evaluator_type = 'poly'
    @classmethod
    def potential_generator_dispatch(cls) -> 'dev.OptionsMethodDispatch':
        cls._pot_dispatch = dev.handle_uninitialized(
            cls._pot_dispatch,
            dev.OptionsMethodDispatch,
            args=(cls.get_potential_types,),
            kwargs=dict(
                attributes_map=cls.get_potentials_by_attributes()
            )
        )
        return cls._pot_dispatch

    def get_g(self, which, return_tf=True, order=1):
        r, tf = nput.internal_coordinate_tensors(self.mol.coords, which,
                                                 return_inverse=True,
                                                 masses=self.mol.atomic_masses,
                                                 order=order)
        b = np.diag(np.repeat(1 / np.sqrt(self.mol.atomic_masses), 3)) @ r[1]
        G = b.T @ b
        if return_tf:
            return G, (r, tf)
        else:
            return G


    def get_anharmonic_parameters(self,
                                  which,
                                  evaluator=None,
                                  energy_expansion=None,
                                  params_handler=None,
                                  # allow_fd=False,
                                  # fd_options=None,
                                  **opts
                                  ):
        if params_handler is None:
            params_handler = self.get_potential_params

        if nput.is_numeric(which[0]):
            which = [which]

        G, (r, tf) = self.get_g(which, return_tf=True, order=4)
        if energy_expansion is None:
            energy_expansion = self.mol.get_cartesian_potential_derivatives(evaluator=evaluator, order=4)
        derivs = nput.tensor_reexpand(tf, energy_expansion, order=4)

        bond_params = []
        for n, spec in enumerate(which):
            diag_derivs = []
            for k, d in enumerate(derivs):
                idx = (n,) * (k + 1)
                diag_derivs.append(d[idx])
            bond_params.append(
                params_handler(spec, r[0][n], G[n, n], diag_derivs, **opts)
            )

        return bond_params

    def get_1d_potentials(self,
                          which,
                          evaluator=None,
                          energy_expansion=None,
                          params_handler=None,
                          potential_params=None,
                          coordinate_potential_handler=None,
                          **opts
                          ):
        if dev.is_dict_like(which):
            potential_params = which
            which = list(which.keys())

        if potential_params is not None:
            if dev.is_dict_like(potential_params):
                which = [
                    coordops.canonicalize_internal(w)
                    for w in which
                ]
                potential_params = {
                    coordops.canonicalize_internal(k): v
                    for k, v in potential_params.items()
                }
                potential_params = [
                    potential_params.get(w)
                    for w in which
                ]
            else:
                potential_params = list(potential_params)
            subwhich = [
                w for w, p in enumerate(potential_params)
                if p is None
            ]
            remwhich = [
                w for w, p in enumerate(potential_params)
                if p is not None
            ]
            gwhich = [
                w for w in remwhich
                if potential_params[w].get('g') is None or potential_params[w].get('re') is None
            ]
            if len(gwhich) > 0:
                G, (r, tf) = self.get_g([which[i] for i in gwhich], return_tf=True)
                for n,w in enumerate(gwhich):
                    potential_params[w] = potential_params[w].copy()
                    if potential_params[w].get('g') is None:
                        potential_params[w]['g'] = G[n, n]
                    if potential_params[w].get('re') is None:
                        potential_params[w]['re'] = r[0][n]
            if len(subwhich) > 0:
                subparams = self.get_anharmonic_parameters([which[k] for k in subwhich],
                                                           evaluator=evaluator,
                                                           energy_expansion=energy_expansion,
                                                           params_handler=params_handler,
                                                           **opts
                                                           )
                for k, p in zip(subwhich, subparams):
                    potential_params[k] = p

        else:
            potential_params = self.get_anharmonic_parameters(
                which,
                evaluator=evaluator,
                energy_expansion=energy_expansion,
                params_handler=params_handler,
                **opts
            )

        if coordinate_potential_handler is None:
            coordinate_potential_handler = self.potential_generator_dispatch()
        pots = []
        for spec, params in zip(which, potential_params):
            handler, props = coordinate_potential_handler.resolve(params)
            pots.append(
                handler(spec, **props)
            )

        return pots