"""
Provides support for handling modes that arise from
"""

import numpy as np, scipy
import McUtils.Numputils as nput
from McUtils.Coordinerds import CoordinateSystem, CartesianCoordinateSystem3D, InternalCoordinateSystem

# from .MoleculeInterface import AbstractMolecule
# from .Transformations import MolecularTransformation

__all__ = [
    "MixtureModes"
]
class MixtureModes(CoordinateSystem):
    """
    A `McUtils.Coordinerds.CoordinateSystem` object that expresses coordinates as
    a rotation on some base set of coordinates with some associated frequencies.
    """
    name="MixtureModes"
    def __init__(self,
                 basis,
                 coeffs,
                 freqs=None,
                 origin=None,
                 masses=None,
                 inverse=None,
                 mass_weighted=False,
                 frequency_scaled=False,
                 g_matrix=None,
                 name=None,
                 ):
        if (
                isinstance(coeffs, np.ndarray) or
                (not nput.is_numeric(coeffs[0]) and nput.is_numeric(coeffs[0][0]))
        ):
            coeffs = [coeffs]
        coeffs = [np.asanyarray(c) for c in coeffs]
        full_coeffs = coeffs
        coeffs = coeffs[0]

        super().__init__(
            matrix=coeffs,
            inverse=inverse,
            name=self.name if name is None else name,
            basis=basis,
            dimension=(coeffs.shape[1],),
            origin=origin
        )
        self.freqs = freqs
        self.masses = masses
        self.mass_weighted = mass_weighted
        self.frequency_scaled = frequency_scaled
        self.g_matrix = g_matrix
        self._extended_coeffs = full_coeffs
        self._inverse_coeffs = None

    @classmethod
    def prep_modes(cls, modes):
        if isinstance(modes, cls) or all(hasattr(modes, b) for b in ['basis', 'matrix', 'inverse', 'freqs']):
            return modes

        matrix = None
        inverse = None
        basis = None
        freqs = None
        if isinstance(modes, dict):
            opts = modes.copy()
            for k in ["matrix", "modes"]:
                if k in opts:
                    matrix = opts[k]
                    del opts[k]
        elif hasattr(modes, 'matrix'):
            matrix = modes.matrix
            opts = {}

            for k in ['inverse', 'basis', 'freqs']:
                if hasattr(modes, k): opts[k] = getattr(modes, k)
        else:
            matrix = np.asanyarray(modes)
            opts = {}

        if 'inverse' in opts:
            inverse = opts['inverse']
            del opts['inverse']
        if 'basis' in opts:
            basis = opts['basis']
            del opts['basis']
        if 'freqs' in opts:
            freqs = opts['freqs']
            del opts['freqs']

        if matrix is None and inverse is None:
            raise ValueError(f"can't prep {cls.__name__} without matrix or inverse")

        if basis is None:
            if matrix.shape[0] % 3 == 0:
                basis = CartesianCoordinateSystem3D
            else:
                basis = InternalCoordinateSystem(dimension=(None, matrix.shape[0]))

        return cls(
            basis,
            matrix,
            inverse=inverse,
            freqs=freqs,
            **opts
        )

    def __getitem__(self, item):
        """
        Takes a slice of the modes
        :param item:
        :type item:
        :return:
        :rtype:
        """

        if isinstance(item, (int, np.integer)):
            item = (item,)
        elif not isinstance(item[0], (int, np.integer)):
            item = tuple(item[0])

        sub_modes = self.matrix[:, item]
        inv = self._inv
        if inv is not None:
            inv = inv[item, :]
        freq = self.freqs[item,]
        return self.modify(
            matrix=sub_modes,
            freqs=freq,
            inverse=inv
        )

    def modify(self,
               matrix=None,
               *,
               freqs=None,
               origin=None,
               masses=None,
               inverse=None,
               name=None,
               mass_weighted=None,
               frequency_scaled=None,
               g_matrix=None
               ):
        return type(self)(
            self.basis,
            self.matrix if matrix is None else matrix,
            freqs=self.freqs if freqs is None else freqs,
            origin=self.origin if origin is None else origin,
            masses=self.masses if masses is None else masses,
            inverse=self.inverse if inverse is None else inverse,
            name=self.name if name is None else name,
            mass_weighted=self.mass_weighted if mass_weighted is None else mass_weighted,
            frequency_scaled=self.frequency_scaled if frequency_scaled is None else frequency_scaled,
            g_matrix=self.g_matrix if g_matrix is None else g_matrix,
        )

    def rotate(self, rot, in_place=False):
        raise NotImplementedError("too confusing...")

    def transform(self, tf, inv=None, origin=None):
        raise NotImplementedError("ambiguous...")
        #TODO: handle Cartesian variant where tf gets broadcasted

        if inv is None:
            if tf.shape[0] != tf.shape[1]:
                raise ValueError("if given only a transformation, need corresponding inverse")
            else:
                inv = np.linalg.inv(tf)
        base_inv = self.inverse
        base_mat = self.matrix
        raise Exception(base_inv @ base_mat, base_mat.shape, base_inv.shape)
        new_mat = tf@base_mat
        new_inv = base_inv@inv

        if origin is None:
            origin = self.origin
            if origin is not None:
                origin = tf@origin.flatten()

        return type(self)(
            self.basis,
            new_mat,
            freqs=self.freqs,
            masses=self.masses,
            origin=origin,
            inverse=new_inv
        )

    @property
    def cartesian_modes(self):
        return self.origin.ndim == 2

    def embed_coords(self, carts):
        flat_carts = (carts - self.origin[np.newaxis]).reshape((len(carts), -1))
        return (flat_carts[:, np.newaxis, :] @ self.inverse.T[np.newaxis]).reshape(
            flat_carts.shape[0],
            self.matrix.shape[1]
        )
    def unembed_coords(self, mode_coords):
        origin = self.origin
        carts = (mode_coords[:, np.newaxis, :] @ self.matrix.T[np.newaxis]).reshape(
            mode_coords.shape[:1] + origin.shape
        )
        carts = carts + origin[np.newaxis]
        return carts

    @property
    def total_transformation(self):
        return self._extended_coeffs
    @property
    def inverse_transformation(self):
        if self._inverse_coeffs is None:
            self._inverse_coeffs = nput.inverse_transformation(self.total_transformation,
                                                               len(self.total_transformation),
                                                               reverse_expansion=[self.inverse])
        return self._inverse_coeffs
    def embed_derivs(self, derivs):
        return nput.tensor_reexpand(self.total_transformation, derivs)
    def unembed_derivs(self, derivs):
        return nput.tensor_reexpand(self.inverse_transformation, derivs)

    @property
    def is_cartesian(self):
        if self.masses is not None:
            return self.matrix.shape[0] // 3 == len(self.masses)
        else:
            return 'Cartesian' in self.basis.name
    @property
    def coords_by_modes(self):
        return self.inverse
    @property
    def modes_by_coords(self):
        return self.matrix

    def _eval_G(self, masses):
        if not self.is_cartesian:  # a hack
            if self.mass_weighted:
                G = self.inverse.T @ self.inverse
            else:
                # if self.origin is None:
                #     raise ValueError("can't get mass weighting matrix (G^-1/2) without structure")
                # G = self.modes_by_coords @
                raise NotImplementedError("non-mass-weighted internal G-matrix not supported")
        else:
            G = np.diag(1 / np.repeat(masses, 3))
        return G

    def _get_gmatrix(self, masses=None):
        if masses is None:
            G = self.g_matrix
        else:
            G = None
        if G is None:
            if masses is None:
                masses = self.masses
                G = self.g_matrix = self._eval_G(self.masses)
            else:
                G = self._eval_G(self.masses)
        g12 = nput.fractional_power(G, 1/2)
        gi12 = nput.fractional_power(G, -1/2)
        return masses, g12, gi12

    @classmethod
    def compute_local_hessian(cls, f, g):
        a = np.diag(np.power(np.diag(g) / np.diag(f), 1 / 4))
        return a @ f @ a

    @classmethod
    def compute_local_gmatrix(cls, f, g):
        a = np.diag(np.power(np.diag(f) / np.diag(g), 1 / 4))
        return a @ g @ a

    def compute_hessian(self):
        pinv = self.coords_by_modes @ self.modes_by_coords
        return pinv @ np.diag(self.freqs ** 2) @ pinv.T
    def compute_gmatrix(self):
        if self.mass_weighted:
            return self.modes_by_coords.T @ self.modes_by_coords
        elif self.g_matrix is not None:
            return self.modes_by_coords.T @ self.g_matrix @ self.modes_by_coords
        elif self.is_cartesian:
            return self.modes_by_coords.T @ np.diag(np.repeat(1/self.masses, 3))  @ self.modes_by_coords
        else:
            return None
            # raise NotImplementedError("non-mass-weighted internal G-matrix not supported")
    def compute_freqs(self):
        # self = self.remove_frequency_scaling().remove_mass_weighting()
        f = self.compute_hessian()
        g = self.compute_gmatrix()
        freqs2, _ = scipy.linalg.eigh(f, g, type=3)
        return np.sign(freqs2) * np.sqrt(np.abs(freqs2))
    @property
    def local_hessian(self):
        # self = self.remove_frequency_scaling()
        f = self.compute_hessian()
        g = self.compute_gmatrix()
        return self.compute_local_hessian(f, g)

    @property
    def local_gmatrix(self):
        # self = self.remove_frequency_scaling()
        f = self.compute_hessian()
        g = self.compute_gmatrix()
        return self.compute_local_gmatrix(f, g)

    @property
    def local_freqs(self):
        return np.diag(self.local_hessian)

    def make_mass_weighted(self, masses=None):
        if self.mass_weighted: return self
        masses, g12, gi12 = self._get_gmatrix(masses=masses)
        L = g12 @ self.modes_by_coords
        Linv = self.coords_by_modes @ gi12
        origin = (self.origin.flatten()[np.newaxis, :] @ gi12).reshape(self.origin.shape)

        return self.modify(L,
                           inverse=Linv,
                           masses=masses,
                           origin=origin,
                           mass_weighted=True
                           )
    def remove_mass_weighting(self, masses=None):
        if not self.mass_weighted: return self
        masses, g12, gi12 = self._get_gmatrix(masses=masses)
        L = gi12 @ self.modes_by_coords
        Linv = self.coords_by_modes @ g12
        origin = (self.origin.flatten()[np.newaxis, :] @ g12).reshape(self.origin.shape)

        return self.modify(L,
                           inverse=Linv,
                           masses=masses,
                           origin=origin,
                           mass_weighted=False
                           )

    def _frequency_scaling(self, freqs=None):
        # L = self.matrix.shape.T
        if freqs is None:
            freqs = self.freqs
        conv = np.sqrt(freqs)
        return freqs, conv

    def make_frequency_scaled(self, freqs=None):
        if self.frequency_scaled: return self
        freqs, conv = self._frequency_scaling(freqs=freqs)
        L = self.matrix * conv[np.newaxis, :]
        Linv = self.inverse / conv[:, np.newaxis] # del_Q X

        return self.modify(L,
                           inverse=Linv,
                           freqs=freqs,
                           frequency_scaled=True
                           )

    def remove_frequency_scaling(self, freqs=None):
        if not self.frequency_scaled: return self
        freqs, conv = self._frequency_scaling(freqs=freqs)
        L = self.matrix / conv[np.newaxis, :]
        Linv = self.inverse * conv[:, np.newaxis] # del_Q X

        return self.modify(L,
                           inverse=Linv,
                           freqs=freqs,
                           frequency_scaled=False
                           )

    def make_dimensionless(self, freqs=None, masses=None):
        new = self
        if not new.mass_weighted: new = new.make_mass_weighted(masses=masses)
        if not new.frequency_scaled: new = new.make_frequency_scaled(freqs=freqs)

        return new

    def make_dimensioned(self, freqs=None, masses=None):
        new = self
        if not new.mass_weighted: new = new.remove_mass_weighting(masses=masses)
        if not new.frequency_scaled: new = new.remove_frequency_scaling(freqs=freqs)
        return new

    # def take_submodes(self, pos, axes=(0, 1, 2)):
    #     if self.is_cartesian:
    #         pos = np.asanyarray(pos)
    #         axes = np.asanyarray(axes)
    #         pos = (axes[np.newaxis, :] + pos[:, np.newaxis]*len(axes)).flatten()
    #     atom_disps = self.matrix[pos]
    #     f_base_sub = atom_disps.T @ np.diag(self.freqs**2) @ atom_disps
    #     return type(self).from_fg(
    #         self.basis.take_subbasis(pos),
    #         f_base_sub,
    #         self.g_matrix[np.ix_(pos, pos)],
    #         ...
    #     )

    def apply_projection(self, proj, project_transrot=True, masses=None, origin=None):
        if project_transrot:
            if masses is None:
                m = self.masses
            else:
                m = masses
            if origin is None:
                o = self.origin
            else:
                o = origin

            tr_proj = nput.translation_rotation_projector(
                np.asanyarray(o).reshape((-1, 3)),
                masses=m,
                mass_weighted=self.mass_weighted
            )
            proj = tr_proj @ proj @ tr_proj



        mbc = proj @ self.modes_by_coords
        cbm = self.coords_by_modes @ proj

        return self.modify(
            matrix=mbc,
            inverse=cbm,
            masses=masses,
            origin=origin
        )

    @classmethod
    def _atom_projector(cls, n, i, orthogonal_projection=False):
        if nput.is_numeric(i):
            i = [i]
        z = np.zeros((3 * n, 3 * n))
        for i in i:
            x = np.arange(3 * i, 3 * (i + 1))
            z[x, x] = 1
        if orthogonal_projection:
            z = np.eye(3 * n) - z
        return z
    def apply_constraints(self,
                          coordinate_constraints,
                          atoms=None,
                          masses=None,
                          origin=None,
                          orthogonal_projection=True,
                          ):

        if nput.is_numeric(coordinate_constraints[0]):
            coordinate_constraints = [coordinate_constraints]

        if masses is not None:
            m = masses
        else:
            m = self.masses
        if origin is not None:
            o = origin
        else:
            o = self.remove_mass_weighting().origin
        basis, _, _ = nput.internal_basis(
            np.asanyarray(o).reshape((-1, 3)),
            coordinate_constraints,
            masses=m,
            project_transrot=False
            )

        gi12 = np.diag(np.repeat(1 / np.sqrt(m), 3))
        if self.mass_weighted:
            projections = [
                nput.projection_matrix(gi12 @ b)
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(gi12 @ b)
                for b in basis
            ]
        else:
            projections = [
                nput.projection_matrix(b)
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(b)
                for b in basis
            ]

        if atoms is not None:
            if nput.is_numeric(atoms):
                atoms = [atoms]
            nats = len(self.masses)
            a_proj = self._atom_projector(nats, atoms)
            projections = [
                a_proj @ proj @ a_proj
                for proj in projections
            ]

        if orthogonal_projection:
            proj = projections[0]
            for p in projections[1:]:
                proj = p @ proj @ p
        else:
            proj = np.sum(projections, axis=0)

        new = self.apply_projection(proj, masses=m, origin=o, project_transrot=False)
        # if self.mass_weighted:
        #     new = new.make_mass_weighted(masses=m)

        return new

    def apply_transformation(self, tf, **opts):
        from .LocalizedModes import LocalizedModes

        return LocalizedModes(self, tf, **opts)
