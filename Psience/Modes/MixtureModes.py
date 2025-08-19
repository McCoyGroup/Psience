"""
Provides support for handling modes that arise from
"""

import enum
import numpy as np, scipy
import McUtils.Numputils as nput
from collections import namedtuple
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
        if inverse is not None:
            inverse = np.asanyarray(inverse)
        if masses is not None:
            masses = np.asanyarray(masses)
        if freqs is not None:
            freqs = np.asanyarray(freqs)
        if g_matrix is not None:
            g_matrix = np.asanyarray(g_matrix)

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


    def to_state(self, serializer=None):
        return {
            'basis':serializer.serialize(self.basis),
            'matrix':self.matrix,
            'inverse':self.inverse,
            'freqs':self.freqs,
            'masses':self.masses,
            'mass_weighted':self.mass_weighted,
            'frequency_scaled':self.frequency_scaled,
            'g_matrix':self.g_matrix
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(
            serializer.deserialize(data['basis']),
            data['matrix'],
            inverse=data['inverse'],
            freqs=data['freqs'],
            masses=data['masses'],
            mass_weighted=data['mass_weighted'],
            frequency_scaled=data['frequency_scaled'],
            g_matrix=data['g_matrix']
        )

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
        elif isinstance(item, slice):
            ...
            # item = np.arange(self.matrix.shape[1])[item]
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
    def compute_local_transformations(cls, f, g):
        return [
            np.diag(np.power(np.diag(g) / np.diag(f), 1 / 4)),
            np.diag(np.power(np.diag(f) / np.diag(g), 1 / 4))
        ]

    @classmethod
    def compute_local_hessian(cls, f, g):
        a = np.diag(np.power(np.diag(g) / np.diag(f), 1 / 4))
        return a @ f @ a

    @classmethod
    def compute_local_gmatrix(cls, f, g):
        a = np.diag(np.power(np.diag(f) / np.diag(g), 1 / 4))
        return a @ g @ a

    def compute_hessian(self, system='modes'):
        if system == 'modes':
            pinv = self.coords_by_modes @ self.modes_by_coords
            return pinv @ np.diag(np.sign(self.freqs) * (self.freqs ** 2)) @ pinv.T
        elif system == 'coords':
            pinv = self.modes_by_coords
            return pinv @ np.diag(np.sign(self.freqs) * (self.freqs ** 2)) @ pinv.T
        else:
            raise ValueError(f'unknown system for normal modes "{system}", valid are "modes", "coords"')

    def compute_gmatrix(self, system='modes', return_fractional=False):
        if system == 'modes':
            if self.mass_weighted:
                g = self.modes_by_coords.T @ self.modes_by_coords
            elif self.g_matrix is not None:
                g = self.modes_by_coords.T @ self.g_matrix @ self.modes_by_coords
            elif self.is_cartesian:
                g = self.modes_by_coords.T @ np.diag(np.repeat(1/self.masses, 3))  @ self.modes_by_coords
            else:
                return None
            if return_fractional:
                g12 = nput.fractional_power(g, 1/2)
                gi12 = nput.fractional_power(g, -1/2)
                return g, g12, gi12
            else:
                return g
                # raise NotImplementedError("non-mass-weighted internal G-matrix not supported")
        elif system == 'coords':
            g, g12, gi12 = self._get_gmatrix()
            if return_fractional:
                return g, g12, gi12
            else:
                return g
        else:
            raise ValueError(f'unknown system for normal modes "{system}", valid are "modes", "coords"')
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

    @property
    def local_mode_transformations(self):
        f = self.compute_hessian()
        g = self.compute_gmatrix()
        return self.compute_local_transformations(f, g)


    def get_nearest_mode_transform(self,
                                   alternate_modes:np.ndarray,
                                   mass_weighted=None,
                                   atoms=None,
                                   maximum_similarity=True,
                                   unitarize=True,
                                   masses=None
                                   ):
        if mass_weighted is not None:
            if not mass_weighted:
                if masses is None:
                    masses = self.masses
                    modes = self.remove_mass_weighting()
                else:
                    modes = self.remove_mass_weighting(masses)
            else:
                if masses is None:
                    masses = self.masses
                    modes = self.make_mass_weighted()
                else:
                    modes = self.make_mass_weighted(masses)
        else:
            modes = self

        inv = modes.coords_by_modes
        modes = modes.modes_by_coords
        if self.is_cartesian and atoms is not None:
            proj = self._atom_projector(modes.shape[0], atoms)
            # if mass_weighted:
            #     gi12 = np.diag(np.repeat(1/np.sqrt(masses), 3))
            #     proj = nput.find_basis(gi12 @ proj @ gi12)
            modes = proj @ modes
            alternate_modes = proj @ alternate_modes

        if maximum_similarity:
            tf = nput.maximum_similarity_transformation(modes, alternate_modes, apply_transformation=False)
            return tf, tf.T
        else:
            tf = inv @ alternate_modes
            if unitarize:
                tf = nput.unitarize_transformation(tf)
                inv = tf.T
            else:
                inv = np.linalg.pinv(tf)
            return tf, inv

    localization_type = 'ned'
    zero_freq_cutoff = 99 / 219474.56
    def get_projected_localized_mode_transformation(self,
                                                    projectors,
                                                    masses=None, origin=None,
                                                    localization_type=None,
                                                    allow_mode_mixing=False,
                                                    maximum_similarity=False,
                                                    unitarize=True,
                                                    zero_freq_cutoff=None,
                                                    orthogonal_projection=False,
                                                    atoms=None #ignored but tedious to have around
                                                    ):
        if zero_freq_cutoff is None:
            zero_freq_cutoff = self.zero_freq_cutoff
        if localization_type is None:
            localization_type = self.localization_type

        if masses is None:
            masses = self.masses
            mw = self.make_mass_weighted()
        else:
            mw = self.make_mass_weighted(masses=masses)

        f = mw.compute_hessian('coords')

        if self.is_cartesian and localization_type == 'direct':
            if origin is None:
                origin = mw.remove_mass_weighting().origin
            tr_tf, tr_inv = nput.translation_rotation_invariant_transformation(
                origin.reshape(-1, 3),
                masses,
                mass_weighted=True
            )
        else:
            tr_tf = tr_inv = np.eye(mw.matrix.shape[0])

        projectors = [
            (
                nput.orthogonal_projection_matrix(p)
                    if orthogonal_projection else
                nput.projection_matrix(p)
            )
                if p.shape[0] != p.shape[1] else
            np.asanyarray(p)
            for p in projectors
        ]

        g_matrix = tr_tf.T @ tr_tf
        if allow_mode_mixing:
            proj = np.sum(projectors, axis=0)
            # proj = proj @ g12
            # proj = self._atom_projector(len(masses), atoms)
            freqs, modes, inv = self.get_normal_modes(
                tr_inv @ proj @ f @ proj.T @ tr_inv.T,
                g_matrix,
                remove_transrot=True,
                mass_weighted=True,
                zero_freq_cutoff=zero_freq_cutoff
            )
            # freqs1, _, _ = self.get_normal_modes(
            #      f,
            #     g_matrix,
            #     remove_transrot=True,
            #     mass_weighted=True,
            #     zero_freq_cutoff=zero_freq_cutoff
            # )
            # raise Exception(freqs * 219474.63, freqs1 * 219474.63)
        else:
            # freqs = []
            modes = []
            # inv = []
            for proj in projectors:
                f_matrix = tr_inv @ proj @ f @ proj.T @ tr_inv.T
                # import McUtils.Formatters as mfmt
                # print(
                #     mfmt.TableFormatter("{:.3f}").format(
                #         f_matrix * 219474.56
                #     )
                # )
                sub_freqs, sub_modes, sub_inv = self.get_normal_modes(
                    f_matrix,
                    g_matrix,
                    remove_transrot=True,
                    mass_weighted=True,
                    zero_freq_cutoff=zero_freq_cutoff
                )
                # print(sub_freqs * 219474.56)

                # freqs.append(sub_freqs)
                modes.append(sub_modes)
                # inv.append(sub_inv)
            # freqs = np.concatenate(freqs)
            modes = np.concatenate(modes, axis=1)
            # inv = np.concatenate(inv, axis=0)

        if localization_type == 'direct':
            modes = tr_tf @ modes

        return mw.get_nearest_mode_transform(
            modes,
            mass_weighted=True,
            maximum_similarity=maximum_similarity,
            unitarize=unitarize
        )

    def get_atom_localized_mode_transformation(self,
                                               atoms,
                                               masses=None, origin=None,
                                               localization_type='ned',
                                               allow_mode_mixing=False,
                                               maximum_similarity=False,
                                               orthogonal_projection=False,
                                               unitarize=True,
                                               zero_freq_cutoff=None
                                               ):
        if nput.is_numeric(atoms):
            atoms = [atoms]
        if masses is None:
            m = self.masses
        else:
            m = masses
        nats = len(m)
        if orthogonal_projection:
            atoms = np.setdiff1d(np.arange(nats), atoms)

        proj = (
            [self._atom_projector(nats, atoms)]
                if allow_mode_mixing else
            [self._atom_projector(nats, a) for a in atoms]
        )
        # gi12 = np.diag(np.repeat(1 / np.sqrt(m), 3))
        # proj = nput.find_basis(gi12 @ proj @ gi12)

        return self.get_projected_localized_mode_transformation(
            proj,
            origin=origin,
            masses=masses,
            localization_type=localization_type,
            maximum_similarity=maximum_similarity,
            unitarize=unitarize,
            allow_mode_mixing=allow_mode_mixing,
            zero_freq_cutoff=zero_freq_cutoff
        )

    def get_coordinate_projected_localized_mode_transformation(self,
                                                               coordinate_constraints,
                                                               atoms=None,
                                                               masses=None, origin=None,
                                                               localization_type='ned',
                                                               allow_mode_mixing=False,
                                                               maximum_similarity=False,
                                                               orthogonal_projection=True,
                                                               unitarize=True
                                                               ):


        if self.is_cartesian:
            if nput.is_numeric(coordinate_constraints[0]):
                coordinate_constraints = [coordinate_constraints]
            nmw = self.remove_mass_weighting()
            if masses is None:
                m = self.masses
            else:
                m = masses
            basis, _, _ = nput.internal_basis(nmw.origin.reshape((-1, 3)), coordinate_constraints, masses=m)

            g12 = np.diag(np.repeat(1/np.sqrt(m), 3))
            projections = [
                nput.projection_matrix(g12 @ b)
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(g12 @ b)
                for b in basis
            ]

            if allow_mode_mixing and orthogonal_projection:
                proj = projections[0]
                for p in projections[1:]:
                    proj = p @ proj @ p
                projections = [proj]

            if atoms is not None:
                if nput.is_numeric(atoms):
                    atoms = [atoms]
                nats = len(m)
                a_proj = self._atom_projector(nats, atoms)
                projections = [
                    a_proj @ proj @ a_proj
                    for proj in projections
                ]
        else:
            if nput.is_numeric(coordinate_constraints):
                coordinate_constraints = [coordinate_constraints]
            if not (
                # Fast check
                isinstance(coordinate_constraints, np.ndarray)
                and np.issubdtype(coordinate_constraints.dtype, np.integer)
            ) and not all(nput.is_int(c) for c in coordinate_constraints):
                raise ValueError(f"internal modes can only be constrained by coordinate index (got {coordinate_constraints})")
            basis = np.zeros((len(coordinate_constraints),self.coords_by_modes.shape[1] ), dtype=float)
            coordinate_constraints = np.asanyarray(coordinate_constraints)
            for n,i in enumerate(coordinate_constraints):
                basis[n, i] = 1


            # g12 = nput.fractional_power(self.g_matrix, 1/2)
            # projections = [
            #     nput.projection_matrix(g12 @ b[:, np.newaxis])
            #         if not orthogonal_projection else
            #     nput.orthogonal_projection_matrix(g12 @ b[:, np.newaxis])
            #     for b in basis
            # ]

            projections = [
                nput.projection_matrix(b[:, np.newaxis])
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(b[:, np.newaxis])
                for b in basis
            ]

        # if orthogonal_projection:
        #     projections = [proj @ gi12 for proj in projections]
        # else:
        #     projections = [proj @ gi12 for proj in projections]


        return self.get_projected_localized_mode_transformation(
            projections,
            origin=origin,
            masses=masses,
            localization_type=localization_type,
            maximum_similarity=maximum_similarity,
            unitarize=unitarize,
            allow_mode_mixing=allow_mode_mixing
        )

    def get_internal_localized_mode_transformation(
            self,
            expansion_coordinates: "Iterable[Iterable[int]|dict]",
            fixed_atoms=None,
            mass_weighted=False,
            project_transrot=True,
            atoms=None,
            maximum_similarity=False,
            orthogonal_projection=False,
            projection=False,
            allow_mode_mixing=False,
            unitarize=True,
            origin=None,
            masses=None,
            localization_type='ned'
    ):
        # from .ObliqueModes import ObliqueModeGenerator

        nmw = self.remove_mass_weighting()
        base_derivs = nput.internal_coordinate_tensors(
            nmw.origin,
            expansion_coordinates,
            fixed_atoms=fixed_atoms
        )

        if projection or orthogonal_projection:
            if masses is None:
                m = self.masses
            else:
                m = masses
            base_derivs = np.diag(np.repeat(1 / np.sqrt(m), 3)) @ base_derivs

            if origin is None:
                origin = self.origin
            if masses is None:
                m = self.masses
            else:
                m = masses
            projector = nput.translation_rotation_projector(
                origin.reshape(-1, 3),
                m,
                mass_weighted=True
            )
            base_derivs = projector @ base_derivs

            # g12 = np.diag(np.repeat(np.sqrt(m), 3))
            gi12 = np.diag(np.repeat(1 / np.sqrt(m), 3))
            if allow_mode_mixing:
                projections = [
                    nput.projection_matrix(gi12 @ base_derivs)
                        if not orthogonal_projection else
                    nput.orthogonal_projection_matrix(gi12 @ base_derivs)
                ]
            else:
                projections = [
                    nput.projection_matrix(gi12 @ b[:, np.newaxis])
                        if not orthogonal_projection else
                    nput.orthogonal_projection_matrix(gi12 @ b[:, np.newaxis])
                    for b in base_derivs.T
                ]

            # if orthogonal_projection:
            #     projections = [proj @ gi12 for proj in projections]
            # else:
            #     projections = [proj @ gi12 for proj in projections]

            if atoms is not None:
                if nput.is_numeric(atoms):
                    atoms = [atoms]
                nats = len(m)
                a_proj = self._atom_projector(nats, atoms)
                projections = [
                    a_proj @ proj @ a_proj
                    for proj in projections
                ]

            return self.get_projected_localized_mode_transformation(
                projections,
                origin=origin,
                masses=masses,
                localization_type=localization_type,
                maximum_similarity=maximum_similarity,
                unitarize=unitarize,
                allow_mode_mixing=allow_mode_mixing
            )

        else:

            if mass_weighted:
                if masses is None:
                    m = self.masses
                else:
                    m = masses
                base_derivs = np.diag(np.repeat(1 / np.sqrt(m), 3)) @ base_derivs

            if project_transrot:
                if origin is None:
                    origin = self.origin
                if masses is None:
                    m = self.masses
                else:
                    m = masses
                projector = nput.translation_rotation_projector(
                        origin.reshape(-1, 3),
                        m,
                        mass_weighted=mass_weighted
                    )
                base_derivs = projector @ base_derivs

            return self.get_nearest_mode_transform(
                base_derivs,
                atoms=atoms,
                mass_weighted=mass_weighted,
                maximum_similarity=maximum_similarity,
                unitarize=unitarize
            )

    def get_displacement_localized_mode_transformation(self,
                                                       mode_blocks=None,
                                                       atoms=None,
                                                       mass_weighted=True,
                                                       unitarize=True,
                                                       **maximizer_opts
                                                       ):
        if not unitarize:
            raise ValueError("Jacobi iterations only give a unitary transformation")

        if mass_weighted:
            modes = self.make_mass_weighted().matrix
        else:
            modes = self.remove_mass_weighting().matrix

        if atoms is not None:
            modes = self._atom_projector(modes.shape[0], atoms) @ modes
        if mode_blocks is None:
            mode_blocks = [np.arange(modes.shape[1])]
        elif nput.is_numeric(mode_blocks[0]):
            mode_blocks = [mode_blocks]

        tf_inv = np.zeros((modes.shape[1], modes.shape[1]), dtype=float)
        for block in mode_blocks:
            _, sub_inv, _ = nput.jacobi_maximize(
                modes[:, block],
                nput.displacement_localizing_rotation_generator,
                **maximizer_opts
            )
            tf_inv[np.ix_(block, block)] = sub_inv

        return tf_inv, tf_inv.T

    def get_mass_scaled_mode_transformation(self,
                                            mass_scaling,
                                            *,
                                            atoms,
                                            localization_cutoff=.8,
                                            num_modes=None,
                                            project_transrot=False,
                                            unitarize=True,
                                            **diag_opts
                                            ):
        from .NormalModes import NormalModes
        self = self.make_mass_weighted()

        f = self.compute_hessian('coords')

        if not nput.is_numeric(mass_scaling):
            mass_scaling = np.asanyarray(mass_scaling)

        if nput.is_int(atoms):
            atoms = [atoms]
        atoms = tuple(atoms)

        scaled_masses = np.array(self.masses)
        scaled_masses[atoms,] *= mass_scaling

        if project_transrot:
            proj = nput.translation_rotation_projector(
                self.remove_mass_weighting().origin,
                masses=scaled_masses,
                mass_weighted=True
            )
        else:
            proj = None
        freqs0, q1, _ = NormalModes.get_normal_modes(f, scaled_masses,
                                                     projector=proj,
                                                     mass_weighted=True,
                                                     **diag_opts
                                                     )

        atom_inds = np.concatenate([
            n * 3 + np.arange(3)
            for n in atoms
        ])
        atom_proj = q1[atom_inds, :]
        max_local_modes = np.linalg.norm(atom_proj, axis=0)

        max_num_modes = len(atoms) * 3
        if localization_cutoff is None:
            num_modes = max_num_modes if num_modes is None else num_modes
            mode_pos = np.argsort(-max_local_modes)[:num_modes]
        else:
            mode_pos = np.where(max_local_modes > localization_cutoff)
            if len(mode_pos) == 0 or len(mode_pos[0]) == 0:
                return None, None

            mode_vals = max_local_modes[mode_pos]
            num_modes = max_num_modes if num_modes is None else num_modes
            ord = np.argsort(-mode_vals)[:num_modes]
            mode_pos = mode_pos[0][ord,]

        tf = self.coords_by_modes@q1[:, mode_pos]
        ord2 = np.argsort(np.argmax(np.abs(tf), axis=0))
        tf = tf[:, ord2]

        return tf, tf.T


    class LocalizationMethods(enum.Enum):
        MaximumSimilarity = 'target_modes'
        AtomLocalized = 'atoms'
        DisplacmentMinimized = 'displacements'
        Internals = 'coordinates'
        CoordinateConstraints = 'constraints'
        Projected = 'projections'
        MassScaled = 'mass_scaling'

    @property
    def localizer_dispatch(self):
        return {
            self.LocalizationMethods.MaximumSimilarity.value:(self.get_nearest_mode_transform, 'target_modes'),
            self.LocalizationMethods.AtomLocalized.value:(self.get_atom_localized_mode_transformation, 'atoms'),
            self.LocalizationMethods.DisplacmentMinimized.value:(self.get_displacement_localized_mode_transformation, 'mode_blocks'),
            self.LocalizationMethods.Internals.value:(self.get_internal_localized_mode_transformation, 'internals'),
            self.LocalizationMethods.CoordinateConstraints.value:(self.get_coordinate_projected_localized_mode_transformation, 'coordinate_constraints'),
            self.LocalizationMethods.Projected.value:(self.get_projected_localized_mode_transformation, 'projections'),
            self.LocalizationMethods.MassScaled.value:(self.get_mass_scaled_mode_transformation, 'mass_scaling'),
        }

    def localize(self,
                 method=None,
                 *,
                 atoms=None,
                 target_modes=None,
                 internals=None,
                 mode_blocks=None,
                 coordinate_constraints=None,
                 projections=None,
                 reorthogonalize=None,
                 mass_scaling=None,
                 unitarize=True,
                 **opts
                 ):

        from .LocalizedModes import LocalizedModes

        if method is None:
            if target_modes is not None:
                method = 'target_modes'
            elif internals is not None:
                method = 'coordinates'
            elif mode_blocks is not None:
                method = 'displacements'
            elif mass_scaling is not None:
                method = 'mass_scaling'
            elif atoms is not None:
                method = 'atoms'
            elif coordinate_constraints is not None:
                method = 'constraints'
            elif projections is not None:
                method = 'projections'
            else:
                method = 'displacements'

        if isinstance(method, str):
            method = self.LocalizationMethods(method)

        if hasattr(method, 'value'):
            method = self.localizer_dispatch.get(method.value, method)

        args = ()
        try:
            method, arg_names = method
        except TypeError:
            ...
        else:
            if isinstance(arg_names, str):
                arg_names = [arg_names]
            all_kw = dict(
                atoms=atoms,
                target_modes=target_modes,
                internals=internals,
                mode_blocks=mode_blocks,
                reorthogonalize=reorthogonalize,
                coordinate_constraints=coordinate_constraints,
                unitarize=unitarize,
                projections=projections,
                mass_scaling=mass_scaling
            )
            args = tuple(all_kw.get(k) for k in arg_names)
            if 'atoms' not in arg_names:
                opts['atoms'] = atoms # taken by all of the localizers

        tf, inv = method(*args, unitarize=unitarize, **opts)
        if tf is None: return None
        # inverse = tf.T @ self.inverse # assumes a unitary localization

        if reorthogonalize is None:
            reorthogonalize = not unitarize
        if reorthogonalize:
            modes = self.make_mass_weighted().matrix @ tf
            tf = tf @ nput.fractional_power(modes.T @ modes, -1 / 2)
            inv = nput.fractional_power(modes.T @ modes, 1 / 2) @ inv

        return LocalizedModes(self, tf, inverse=inv)

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
        conv = np.sign(freqs) * np.sqrt(np.abs(freqs))
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

    def make_oblique(self):
        from .ObliqueModes import ObliqueModeGenerator

        f = self.compute_hessian()
        g = self.compute_gmatrix()

        _, _, u, ui = ObliqueModeGenerator(f, g).run()
        return self.apply_transformation((u, ui))

