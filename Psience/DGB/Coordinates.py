
import abc, numpy as np, functools
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput

from ..Molecools import StructuralProperties

from .Evaluators import *
from .Interpolation import *


__all__ = [
    "DGBCoords",
    "DGBCartesians",
    "DGBInternals",
    "DGBWatsonModes",
]

class DGBCoords(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def centers(self) -> 'np.ndarray':
        ...
    @property
    def shape(self):
        return self.centers.shape

    @property
    @abc.abstractmethod
    def kinetic_energy_evaluator(self) -> 'DGBKineticEnergyEvaluator':
        ...

    @property
    @abc.abstractmethod
    def pairwise_potential_evaluator_type(self) -> 'type[DGBPairwisePotentialEvaluator]':
        ...
    def pairwise_potential_evaluator(self, potential_functions) -> 'DGBPairwisePotentialEvaluator':
        if not isinstance(potential_functions, dict):
            raise ValueError("don't know what to do with pairwise functions {}".format(potential_functions))
        if 'functions' in potential_functions:
            opts = potential_functions.copy()
            potential_functions = potential_functions['functions']
            del opts['functions']
        else:
            opts = {}
        return self.pairwise_potential_evaluator_type(self, potential_functions, **opts)
    @abc.abstractmethod
    def __getitem__(self, item) -> 'DGBCoords':
        ...
    def take_indices(self, subinds) -> 'DGBCoords':
        return self[:, subinds] # default
    def drop_indices(self, subinds) -> 'DGBCoords':
        remaining = np.setdiff1d(np.arange(self.centers.shape[-1]), subinds)
        return self.take_indices(remaining) # default
    @abc.abstractmethod
    def gmatrix(self, coords:np.ndarray) -> np.ndarray:
        ...

    @classmethod
    def embedded_mode_function(cls, func, modes, masses=None):
        if isinstance(func, DGBInterpolator):
            centers = func.centers
            pot_vals = func.derivs
            opts = func.opts
            centers = DGBWatsonModes.embed_coords(centers, modes)
            pot_vals = DGBWatsonModes.embed_derivs(pot_vals, modes)
            return DGBWatsonInterpolator(centers, pot_vals, modes, **opts)
        else:
            @functools.wraps(func)
            def embedded_function(mode_coords, deriv_order=None):
                carts = DGBWatsonModes.unembed_coords(mode_coords, modes, masses)
                vals = func(carts, deriv_order=deriv_order)
                if deriv_order is not None:
                    vals = DGBWatsonModes.embed_derivs(vals, modes)
                return vals
            return embedded_function

    @classmethod
    def embedded_subcoordinate_function(cls, func, sel, ndim):
        if isinstance(func, DGBInterpolator):
            raise NotImplementedError("reembedding of interpolator not supported")
        @functools.wraps(func)
        def embedded_function(subcoords, deriv_order=None):
            full_shape = (subcoords.shape[0], ndim)
            if sel is None:
                coords = subcoords.reshape(full_shape)
            else:
                coords = np.zeros(full_shape)
                coords[:, sel] = subcoords
            vals = func(coords, deriv_order=deriv_order)
            if sel is not None and deriv_order is not None:
                fshape = vals[0].ndim
                _ = []
                for n, d in enumerate(vals):
                    for j in range(n):
                        d = np.take(d, sel, axis=j + fshape)
                    _.append(d)
                vals = _
            return vals

        return embedded_function

    @classmethod
    def embedded_cartesian_function(cls, func, atom_sel, xyz_sel, natoms, ndim):
        if (atom_sel is not None or xyz_sel is not None):
            full_sel = np.arange(natoms * ndim).reshape(natoms, ndim)
            if atom_sel is None:
                full_sel = full_sel[:, xyz_sel]
            elif xyz_sel is None:
                full_sel = full_sel[atom_sel, :]
            else:
                full_sel = full_sel[np.ix_(atom_sel, xyz_sel)]
            flat_sel = full_sel.flatten()
        else:
            flat_sel = None

        if flat_sel is None:
            return func

        if isinstance(func, DGBInterpolator):
            centers = func.centers
            centers = centers.reshape((centers.shape[0], natoms, ndim))
            pot_vals = func.derivs
            opts = func.opts

            if atom_sel is None:
                centers = centers[..., xyz_sel]
            elif xyz_sel is None:
                centers = centers[..., atom_sel, :]
            else:
                centers = centers[..., np.ix_(atom_sel, xyz_sel)]
            fshape = pot_vals[0].ndim
            _ = []
            for n, d in enumerate(pot_vals):
                for j in range(n):
                    d = np.take(d, flat_sel, axis=j + fshape)
                _.append(d)
            pot_vals = _
            if isinstance(func, DGBWatsonInterpolator):
                return DGBCartesianWatsonInterpolator(centers, pot_vals, func.modes, **opts)
            else:
                return DGBInterpolator(centers, pot_vals, **opts)
        else:
            if flat_sel is None:
                return func
            @functools.wraps(func)
            def embedded_function(subcart_coords, deriv_order=None, flat_sel=flat_sel):
                full_shape = (subcart_coords.shape[0], natoms, ndim)
                if subcart_coords.shape == full_shape:
                    carts = subcart_coords
                else:
                    subcart_coords = np.reshape(
                        subcart_coords,
                        (
                            subcart_coords.shape[0],
                            len(atom_sel) if atom_sel is not None else natoms,
                            len(xyz_sel) if xyz_sel is not None else ndim
                        )
                    )
                    if atom_sel is None and xyz_sel is None:
                        carts = subcart_coords.reshape(full_shape)
                    else:
                        carts = np.zeros(full_shape)
                        if atom_sel is None:
                            carts[..., xyz_sel] = subcart_coords
                        elif xyz_sel is None:
                            carts[..., atom_sel, :] = subcart_coords
                        else:
                            carts[..., np.ix_(atom_sel, xyz_sel)] = subcart_coords
                vals = func(carts, deriv_order=deriv_order)
                if deriv_order is not None:
                    fshape = vals[0].ndim
                    _ = []
                    for n, d in enumerate(vals):
                        for j in range(n):
                            d = np.take(d, flat_sel, axis=j + fshape)
                        _.append(d)
                    vals = _
                return vals

        return embedded_function

    class DGBEmbeddedFunction:
        def __init__(self, embedded_function, original_function, coords):
            self.og_fn = original_function
            self.base_coords = coords
            self.embed_fn = embedded_function
        def __call__(self, coords, deriv_order=None):
            return self.embed_fn(coords, deriv_order=deriv_order)
    @abc.abstractmethod
    def embed_function(self, fn) -> 'DGBEmbeddedFunction':
        ...

    def as_cartesians(self) -> 'tuple[DGBCartesians, tuple[np.ndarray, np.ndarray]]':
        raise NotImplementedError("{} can't be converted to Cartesians".format(
            type(self).__name__
        ))

class DGBCartesians(DGBCoords):
    def __init__(self, coords, masses, *, natoms=None, atom_sel=None, ndim=None, xyz_sel=None):
        coords = np.asanyarray(coords)
        masses = np.asanyarray(masses)
        if masses.ndim > 1:
            raise ValueError("expected a vector of masses...")
        coords = coords.reshape((len(coords), len(masses), -1))

        self.coords = coords
        self.masses = masses
        self.natoms = len(masses) if natoms is None else natoms# the original numebr of atoms
        self.ndim = self.coords.shape[-1] if ndim is None else ndim
        self.atom_sel = atom_sel
        self.xyz_sel = xyz_sel
    @property
    def centers(self):
        return self.coords.reshape((self.coords.shape[0], self.coords.shape[1]*self.coords.shape[2]))
    @property
    def cart_shape(self):
        return self.coords.shape
    @property
    def kinetic_energy_evaluator(self):
        return DGBCartesianEvaluator(
            np.broadcast_to(self.masses[:, np.newaxis], self.coords.shape[1:]).flatten()
        )
    @property
    def pairwise_potential_evaluator_type(self):
        return DGBCartesianPairwiseEvaluator
    @classmethod
    def resolve_masses(cls, coords, masses=None, atoms=None):
        if masses is None:
            if atoms is not None:
                atoms = [
                    AtomData[a, "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
                        if isinstance(a, str) else
                    a
                    for a in atoms
                ]
                masses = np.array(atoms)
            else:
                raise ValueError("either atoms or masses must not be `None`")

        if isinstance(masses, (int, np.integer, float, np.floating)):
            masses = [masses] * coords.shape[1]
        return masses
    @classmethod
    def from_cartesians(cls, centers, masses=None, atoms=None):
        centers = np.asanyarray(centers)
        masses = cls.resolve_masses(centers, masses=masses, atoms=atoms)
        return cls(centers, masses)

    def infer_shape_sel(self, selector):
        n, na, nd = self.coords.shape
        test = np.arange(na*nd).reshape((1, na, nd))
        inferrable = np.broadcast_to(test, self.coords.shape)[selector][0]

        test_atoms = test[0, :, 0] // nd
        sub_atoms = inferrable[:, 0] // nd
        # we map the atoms back to their positions in the original list
        atom_sel = np.searchsorted(test_atoms, sub_atoms)
        if len(atom_sel) == na and np.all(atom_sel == np.arange(na)):
            atom_sel = None

        test_xyz = test[0, 0] % nd
        sub_xyz = inferrable[0] % nd
        xyz_sel, _ = nput.find(test_xyz, sub_xyz)
        if len(xyz_sel) == nd and np.all(xyz_sel == np.arange(nd)):
            xyz_sel = None

        return atom_sel, xyz_sel
    @staticmethod
    def _merge_sel(new_sel, old_sel):
        if old_sel is not None:
            new_sel = np.asanyarray(old_sel)[new_sel,]
        return new_sel
    def __getitem__(self, item) -> 'DGBCartesians':
        c = self.coords[item]
        if not isinstance(c, np.ndarray) or c.ndim != 3:
            raise ValueError("bad slice {}".format(item))

        if c.ndim < 3:
            raise ValueError("zero-length slice?")

        if c.shape[1:] == self.coords.shape[1:]:
            return type(self)(
                c,
                self.masses,
                atom_sel=self.atom_sel,
                xyz_sel=self.xyz_sel,
                natoms=self.natoms,
                ndim=self.ndim
            )
        else:
            atom_sel, xyz_sel = self.infer_shape_sel(item)
            return type(self)(
                c,
                self.masses if atom_sel is None else self.masses[atom_sel,],
                atom_sel=self._merge_sel(atom_sel, self.atom_sel),
                xyz_sel=self._merge_sel(xyz_sel, self.xyz_sel),
                natoms=self.natoms,
                ndim=self.ndim
            )
    def take_indices(self, subinds):
        subinds = np.asanyarray(subinds)
        atom_inds = np.unique(subinds // self.coords.shape[2])
        xyz_inds = np.unique(subinds % self.coords.shape[2])
        return self[:, atom_inds, :][:, :, xyz_inds]
    def embed_function(self, function):
        """
        Embeds assuming we got a function in Cartesians _before_ any selections happened

        :param function:
        :return:
        """
        return self.DGBEmbeddedFunction(
            self.embedded_cartesian_function(
                function,
                self.atom_sel,
                self.xyz_sel,
                self.natoms,
                self.ndim
            ),
            function,
            self
        )

    def gmatrix(self, coords:np.ndarray) -> np.ndarray:
        mass_spec = np.broadcast_to(self.masses[:, np.newaxis], (len(self.masses), self.coords.shape[2])).flatten()
        mass_spec = np.diag(1 / mass_spec)
        return np.broadcast_to(mass_spec[np.newaxis], (coords.shape[0],) + mass_spec.shape)
class DGBInternals(DGBCoords):
    def __init__(self, coords, gmat_function=None, vprime_function=None):
        raise NotImplementedError('internals coming soon?')
class DGBWatsonModes(DGBCoords):
    def __init__(self, coords, modes,
                 *,
                 coriolis_inertia_function=None,
                 masses=None,
                 subselection=None
                 ):
        self.coords = coords
        self.masses = masses
        self.modes = modes
        self.subsel = subselection
        self.ci_func = coriolis_inertia_function

    @property
    def centers(self):
        return self.coords
    @property
    def kinetic_energy_evaluator(self):
        return DGBWatsonEvaluator(self.modes, self.ci_func)
    @property
    def pairwise_potential_evaluator_type(self):
        return DGBWatsonPairwiseEvaluator
    @staticmethod
    def zeta_momi(watson_coords, modes, masses):

        carts = DGBWatsonModes.unembed_coords(watson_coords, modes, masses=masses)
        if carts.shape[-1] < 3:
            # pad with zeros
            carts = np.concatenate(
                [
                    carts,
                    np.zeros(carts.shape[:2] + (3 - carts.shape[-1],))
                ],
                axis=-1
            )

        # scaled_modes = modes.matrix.T / np.sqrt(modes.freqs[:, np.newaxis]) # used in the VPT version?
        zeta, (B_e, eigs) = StructuralProperties.get_prop_coriolis_constants(carts,
                                                                             modes.matrix.T,
                                                                             masses
                                                                             )

        return zeta, B_e

    @classmethod
    def default_coriolis_inertia_function(cls, modes, masses):
        def coriolis_inertia_function(watson_coords, *, modes=modes, masses=masses, deriv_order=None):
            if deriv_order is not None:
                raise NotImplementedError("don't have support for Coriolis derivatives in DGB")
            zeta, mom_i = cls.zeta_momi(watson_coords, modes, masses)
            # freqs = modes.freqs
            # freq_term = np.sqrt(freqs[np.newaxis, :] / freqs[:, np.newaxis])
            # zeta = zeta * freq_term[np.newaxis, np.newaxis]
            bad_moms = np.where(mom_i <= 0)[0]
            mom_i[bad_moms] = 1
            B_e = 1 / (2 * mom_i)  # * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass"))
            B_e[bad_moms] = 0
            return np.sum(B_e, axis=-1), sum(
                B_e[:, a, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
                 * zeta[:, a, :, :, np. newaxis, np.newaxis]*zeta[:, a, np. newaxis, np.newaxis, :, :]
                for a in range(B_e.shape[1])
            )

        return coriolis_inertia_function

    @classmethod
    def embed_coords(cls, carts, modes, shift=True):
        if shift:
            flat_carts = (carts - modes.origin[np.newaxis]).reshape((len(carts), -1))
        else:
            flat_carts = (carts).reshape((len(carts), -1))
        return (flat_carts[:, np.newaxis, :] @ modes.inverse.T[np.newaxis]).reshape(
            flat_carts.shape[0],
            modes.matrix.shape[1]
        )
    @classmethod
    def unembed_coords(cls, mode_coords, modes, masses=None, shift=True):
        origin = modes.origin
        carts = (mode_coords[:, np.newaxis, :] @ modes.matrix.T[np.newaxis]).reshape(
            mode_coords.shape[:1] + origin.shape
        )
        if shift:
            carts = carts + origin[np.newaxis, :, :]
        # if masses is not None:
        #     carts = carts.reshape((carts.shape[0], len(masses), -1))
        #     carts = StructuralProperties.get_eckart_embedded_coords(
        #         masses,
        #         origin,
        #         carts
        #     )
        return carts
    @classmethod
    def embed_derivs(cls, derivs, modes):
        fshape = derivs[0].ndim
        _ = []
        for n, d in enumerate(derivs):
            for j in range(n):
                d = np.tensordot(d, modes.matrix, axes=[fshape, 0])
            _.append(d)
        return _
    @classmethod
    def from_cartesians(cls,
                        coords,
                        modes,
                        masses=None,
                        coriolis_inertia_function=None
                        ):
        coords = cls.embed_coords(coords, modes)
        if coriolis_inertia_function is None:
            coriolis_inertia_function = cls.default_coriolis_inertia_function(modes, masses)
        return cls(coords,
                   modes,
                   coriolis_inertia_function=coriolis_inertia_function,
                   masses=masses
                   )
    def as_cartesians(self, masses=None) -> 'tuple[DGBCartesians, tuple[np.ndarray, np.ndarray]]':
        mode_coords = self.coords
        modes = self.modes
        if masses is None:
            masses = self.masses

        carts = self.unembed_coords(mode_coords, modes, self.masses)
        if masses is None:
            raise ValueError("need `masses` to be able to convert back to Cartesians")
        carts = carts.reshape((carts.shape[0], len(masses), -1))

        return DGBCartesians(carts, masses), (modes.matrix, modes.inverse)

    def __getitem__(self, item):
        c = self.coords[item]
        if not isinstance(c, np.ndarray) or c.ndim != self.coords.ndim:
            raise ValueError("bad slice {}".format(item))
        if c.shape[1] != self.coords.shape[1]:
            # means we took some subselection of modes
            test_sel = np.broadcast_to(
                np.arange((self.coords.shape[1]))[np.newaxis],
                self.coords.shape
            ) # something that can be indexed the same way
            subsel = test_sel[item][0]
            modes = self.modes[subsel]
            ci_func = self.embedded_subcoordinate_function(
                self.ci_func,
                subsel,
                self.coords.shape[1]
            )
        else:
            modes = self.modes
            subsel = None
            ci_func = self.ci_func
        if self.subsel is not None:
            subsel = self.subsel[subsel]
        return type(self)(
            c,
            modes,
            coriolis_inertia_function=ci_func,
            masses=self.masses,
            subselection=subsel
        )
    def gmatrix(self, coords:np.ndarray) -> np.ndarray:
        base_spec = np.eye(coords.shape[1])
        return np.broadcast_to(base_spec[np.newaxis], (coords.shape[0],) + base_spec.shape)

    def embed_function(self, fn):
        return self.DGBEmbeddedFunction(
            self.embedded_mode_function(  # a function that takes in normal mode coordinates
                fn,
                self.modes,
                masses=self.masses
            ),
            fn,
            self
        )