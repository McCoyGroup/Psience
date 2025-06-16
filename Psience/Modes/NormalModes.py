import enum
import functools

import numpy as np, scipy.linalg as slag, itertools, collections
import scipy.linalg

from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput

from .MixtureModes import MixtureModes

__all__ = [
    "NormalModes",
    "ReactionPathModes"
]

__reload_hook__ = [".MixtureModes"]

class NormalModes(MixtureModes):
    name="NormalModes"

    def __init__(self,
                 basis,
                 coeffs,
                 freqs=None,
                 origin=None,
                 masses=None,
                 inverse=None,
                 name=None,
                 mass_weighted=False,
                 frequency_scaled=False,
                 g_matrix=None
                 ):
        super().__init__(
            basis,
            coeffs,
            freqs=freqs,
            origin=origin,
            masses=masses,
            inverse=inverse,
            mass_weighted=mass_weighted,
            frequency_scaled=frequency_scaled,
            g_matrix=g_matrix,
            name=name
        )

    ModeData = collections.namedtuple("ModeData", ['freqs', 'modes', 'inverse'])
    default_zero_freq_cutoff = 1.0e-4 # 20 wavenumbers...
    @classmethod
    def get_normal_modes(cls,
                         f_matrix,
                         mass_spec,
                         # mass_units="AtomicMassUnits",
                         remove_transrot=True,
                         dimensionless=False,
                         mass_weighted=None,
                         zero_freq_cutoff=None,
                         return_gmatrix=False,
                         projector=None
                         ):

        f_matrix = np.asanyarray(f_matrix)
        base_shape = f_matrix.shape[:-2]
        f_matrix = f_matrix.reshape((-1,) + f_matrix.shape[-2:])
        if isinstance(mass_spec[0], str):  # atoms were supplied
            mass_spec = np.array([AtomData[a, "Mass"] for a in mass_spec]) * UnitsData.convert(
                "AtomicMassUnits",
                "AtomicUnitOfMass"
            )
        mass_spec = np.asanyarray(mass_spec)

        if mass_spec.ndim == 1:
            if mass_weighted is None: mass_weighted = dimensionless # generally do the right thing for Cartesians
            mass_spec = np.broadcast_to(mass_spec[:, np.newaxis], (len(mass_spec), 3)).flatten()
            mass_spec = np.diag(1 / mass_spec)

        mass_spec = np.reshape(mass_spec, (-1,) + mass_spec.shape[-2:])
        if mass_spec.shape[0] < f_matrix.shape[0]:
            mass_spec = np.broadcast_to(mass_spec, f_matrix.shape)

        gi12 = nput.fractional_power(mass_spec, -1 / 2)
        g12 = nput.fractional_power(mass_spec, 1 / 2)
        if projector is None:
            if f_matrix.shape[0] == 1:
                freq2, modes = slag.eigh(f_matrix[0], mass_spec[0], type=3)
                V = gi12[0] @ modes
                V = V[np.newaxis]
                freq2 = freq2[np.newaxis]
                modes = modes[np.newaxis]
            else:
                mw_F = g12 @ f_matrix @ np.moveaxis(g12, -1, -2)
                freq2, V = np.linalg.eigh(mw_F)
                modes = g12 @ V
        else:
            pG12 = projector @ g12
            mw_F = pG12 @ f_matrix @ np.moveaxis(pG12, -1, -2)
            freq2, V = np.linalg.eigh(mw_F)
            modes = g12 @ V

        # modes = g12 @ V
        # therefore, inv = V.T @ gi12
        # and modes[:, nz] = g12 @ V[:, nz]; inv[nz, :] = (V.T)[nz, :] @ gi12 (this is only the _left_ inverse)

        blocks = False
        if remove_transrot:
            if zero_freq_cutoff is None:
                zero_freq_cutoff = cls.default_zero_freq_cutoff
            if zero_freq_cutoff > 0:
                nonzero = np.abs(freq2) >= zero_freq_cutoff**2
            else:
                nonzero = np.full(freq2.shape, True)

            if gi12.shape[0] != modes.shape[0]:
                gi12 = np.broadcast_to(gi12, f_matrix.shape)

            # need to check if shape is good
            nz_counts = np.sum(nonzero, axis=1)
            blocks = len(freq2) == 1 or np.sum(np.abs(np.diff(nz_counts))) != 0
            freq2_ = []
            modes_ = []
            inv_ = []
            for nz2, f2, m2, v2, gi in zip(nonzero, freq2, modes, V, gi12):
                freq2_.append(f2[nz2])
                modes_.append(m2[:, nz2])
                inv_.append(v2.T[nz2, :] @ gi)
            freq2 = freq2_
            modes = modes_
            inv = inv_
            del freq2_
            del modes_
            del inv_
            if not blocks:
                freq2 = np.array(freq2)
                modes = np.array(modes)
                inv = np.array(inv)

            # else:
            #     if nz_counts[0] == 0:
            #         inv = np.moveaxis(V, -1, -2) @ gi12
            #     else:
            #         nz_spec = [np.where(nz)[0] for nz in nonzero]
            #         freq2 = freq2[nonzero].reshape((modes.shape[0], nz_counts[0]))
            #         print(nz_spec, modes.shape)
            #         modes = np.moveaxis(
            #             nput.vector_take(
            #                 np.moveaxis(modes, -1, -2),
            #                 nz_spec,
            #                 shared=1
            #             ),
            #             -1, -2
            #         )
            #         inv = np.moveaxis(
            #                 nput.vector_take(
            #                     np.moveaxis(V, 1, 0),
            #                     nz_spec,
            #                     shared=1
            #                 )
            #         print(modes.shape, inv.shape)
        else:
            inv = np.moveaxis(V, -1, -2) @ gi12
        # else:
        #     if np.linalg.det(modes) == 1 and np.allclose(modes.T @ modes, np.eyelen(modes)):
        #         inv = modes.T
        #     else:
        #         inv = np.linalg.inv(modes)

        if blocks:
            freqs = [np.sign(f2) * np.sqrt(np.abs(f2)) for f2 in freq2]
        else:
            freqs = np.sign(freq2) * np.sqrt(np.abs(freq2))

        if mass_weighted:
            if blocks:
                modes_ = []
                inv_ = []
                for m,i,g,gi in zip(modes, inv, g12, gi12):
                    modes_.append(gi @ m)
                    inv_.append(i @ g)
                modes = modes_
                inv = inv_
            else:
                # gi12 = slag.fractional_matrix_power(mass_spec, -1 / 2)
                # g12 = slag.fractional_matrix_power(mass_spec, 1 / 2)
                inv = inv @ g12
                modes = gi12 @ modes
        if dimensionless:
            if blocks:
                modes_ = []
                inv_ = []
                for f,m,i in zip(freqs, modes, inv):
                    weighting = np.sqrt(np.abs(f))
                    modes_.append(m / weighting[:, np.newaxis, :])
                    inv_.append(i * weighting[:, :, np.newaxis])
                modes = modes_
                inv = inv_
            else:
                weighting = np.sqrt(np.abs(freqs))
                modes = modes / weighting[np.newaxis, :]
                inv = inv * weighting[:, np.newaxis]
        # if mode == 'reasonable':
        # modes, inv = inv.T, modes.T

        if blocks:
            modes_ = []
            inv_ = []
            for f, m, i in zip(freqs, modes, inv):
                sorting = np.argsort(f)
                modes_.append(m[:, sorting])
                inv_.append(i[sorting, :])
            modes = modes_
            inv = inv_
            modes, inv = [np.moveaxis(i, -1, -2) for i in inv], [np.moveaxis(m, -1, -2) for m in modes]
            # TODO: iterative reshaping
        else:
            modes, inv = np.moveaxis(inv, -1, -2), np.moveaxis(modes, -1, -2)
            freqs = np.reshape(freqs, base_shape + freqs.shape[1:])
            modes = np.reshape(modes, base_shape + modes.shape[1:])
            inv = np.reshape(inv, base_shape + inv.shape[1:])

        if len(base_shape) == 0:
            freqs = freqs[0]
            modes = modes[0]
            inv = inv[0]

        mode_data = cls.ModeData(freqs, modes, inv)
        if return_gmatrix:
            mass_spec = mass_spec.reshape(base_shape + mass_spec.shape[-2:])
            return mode_data, mass_spec
        else:
            return mode_data

    @classmethod
    def from_fg(cls,
                basis,
                f_matrix,
                mass_spec,
                remove_transrot=True,
                dimensionless=False,
                zero_freq_cutoff=None,
                mass_weighted=None,
                origin=None,
                projector=None,
                **opts
                ):
        """
        Generates normal modes from the specified F and G matrices

        :param basis:
        :param f_matrix: second derivatives of the potential
        :param mass_spec:
        :param mass_units:
        :param remove_transrot:
        :param opts:
        :return:
        """

        (freqs, modes, inv), g_matrix = cls.get_normal_modes(f_matrix, mass_spec,
                                                             remove_transrot=remove_transrot,
                                                             dimensionless=dimensionless,
                                                             zero_freq_cutoff=zero_freq_cutoff,
                                                             mass_weighted=mass_weighted,
                                                             return_gmatrix=True,
                                                             projector=projector
                                                             )

        mass_weighted = mass_weighted or dimensionless

        if mass_weighted and origin is not None:
            origin = np.asanyarray(origin)
            origin = (origin.flatten() @ nput.fractional_power(g_matrix, -1/2)).reshape(origin.shape)

        return cls(basis, modes, inverse=inv, freqs=freqs,
                   mass_weighted=mass_weighted,
                   frequency_scaled=dimensionless,
                   g_matrix=g_matrix,
                   origin=origin,
                   **opts
                   )

    default_projected_zero_freq_cutoff = None
    @classmethod
    def from_molecule(cls, mol,
                      dimensionless=False,
                      use_internals=None,
                      potential_derivatives=None,
                      project_transrot=True,
                      zero_freq_cutoff=None,
                      masses=None,
                      **opts
                      ):
        from ..Molecools import Molecule

        mol = mol  # type:Molecule
        if use_internals is None:
            use_internals = mol.internal_coordinates is not None
        if potential_derivatives is None:
            potential_derivatives = mol.potential_derivatives
            if potential_derivatives is None and mol.energy_evaluator is not None:
                potential_derivatives = mol.calculate_energy(order=2)[1:]
        if use_internals:
            hess = nput.tensor_reexpand(
                mol.get_cartesians_by_internals(1, strip_embedding=True, reembed=True),
                [0, potential_derivatives[1]],
                order=2
            )[1]
            if zero_freq_cutoff is None:
                zero_freq_cutoff = cls.default_projected_zero_freq_cutoff
                if zero_freq_cutoff is None:
                    zero_freq_cutoff = 1 * UnitsData.convert("Wavenumbers", "Hartrees")
            return cls.from_fg(
                mol.internal_coordinates.system,
                hess,
                mol.get_gmatrix(masses=masses),
                dimensionless=dimensionless,
                origin=mol.get_internals(strip_embedding=True),
                zero_freq_cutoff=zero_freq_cutoff,
                **opts
            )
        else:
            if not project_transrot and mol.normal_modes.modes is not None:
                return mol.normal_modes.modes.basis.to_new_modes()

            if masses is None:
                masses = mol.atomic_masses
            hess = potential_derivatives[1]
            if project_transrot:
                proj = mol.get_translation_rotation_projector(mass_weighted=True)
                mw = np.diag(np.repeat(1 / np.sqrt(masses), 3))
                mw_inv = np.diag(np.repeat(np.sqrt(masses), 3))
                proj = mw_inv @ proj @ mw
                hess = proj @ hess @ proj.T
                if zero_freq_cutoff is None:
                    zero_freq_cutoff = cls.default_projected_zero_freq_cutoff
                    if zero_freq_cutoff is None:
                        zero_freq_cutoff = 1 * UnitsData.convert("Wavenumbers", "Hartrees")

            return cls.from_fg(
                mol.coords.system,
                hess,
                masses,
                dimensionless=dimensionless,
                masses=masses,
                origin=mol.coords,
                zero_freq_cutoff=zero_freq_cutoff,
                **opts
            )

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
                                                    zero_freq_cutoff=None
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

        f = mw.modes_by_coords @ np.diag(self.freqs ** 2) @ mw.modes_by_coords.T

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
        else:
            # freqs = []
            modes = []
            # inv = []
            for proj in projectors:
                f_matrix = tr_inv @ proj @ f @ proj.T @ tr_inv.T
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

            g12 = nput.fractional_power(self.g_matrix, 1/2)
            projections = [
                nput.projection_matrix(g12 @ b[:, np.newaxis])
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(g12 @ b[:, np.newaxis])
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

    class LocalizationMethods(enum.Enum):
        MaximumSimilarity = 'target_modes'
        AtomLocalized = 'atoms'
        DisplacmentMinimized = 'displacements'
        Internals = 'coordinates'
        CoordinateConstraints = 'constraints'
        Projected = 'projections'

    @property
    def localizer_dispatch(self):
        return {
            self.LocalizationMethods.MaximumSimilarity.value:(self.get_nearest_mode_transform, 'target_modes'),
            self.LocalizationMethods.AtomLocalized.value:(self.get_atom_localized_mode_transformation, 'atoms'),
            self.LocalizationMethods.DisplacmentMinimized.value:(self.get_displacement_localized_mode_transformation, 'mode_blocks'),
            self.LocalizationMethods.Internals.value:(self.get_internal_localized_mode_transformation, 'internals'),
            self.LocalizationMethods.CoordinateConstraints.value:(self.get_coordinate_projected_localized_mode_transformation, 'coordinate_constraints'),
            self.LocalizationMethods.Projected.value:(self.get_projected_localized_mode_transformation, 'projections'),
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
                unitarize=unitarize
            )
            args = tuple(all_kw.get(k) for k in arg_names)
            if 'atoms' not in arg_names:
                opts['atoms'] = atoms # taken by all of the localizers

        tf, inv = method(*args, unitarize=unitarize, **opts)
        # inverse = tf.T @ self.inverse # assumes a unitary localization

        if reorthogonalize is None:
            reorthogonalize = not unitarize
        if reorthogonalize:
            modes = self.make_mass_weighted().matrix @ tf
            tf = tf @ nput.fractional_power(modes.T @ modes, -1 / 2)
            inv = nput.fractional_power(modes.T @ modes, 1 / 2) @ inv

        return LocalizedModes(self, tf, inverse=inv)

class ReactionPathModes(NormalModes):

    @classmethod
    def get_rp_modes(cls,
                     gradient,
                     f_matrix,
                     mass_spec,
                     # mass_units="AtomicMassUnits",
                     remove_transrot=True,
                     dimensionless=False,
                     mass_weighted=None,
                     zero_freq_cutoff=None,
                     return_gmatrix=False,
                     projector=None,
                     zero_gradient_cutoff=2.5e-5, #5 cm-1
                     use_max_gradient_cutoff=True,
                     return_indices=False
                     ):

        gradient = np.asanyarray(gradient)
        f_matrix = np.asanyarray(f_matrix)
        base_shape = f_matrix.shape[:-2]
        if f_matrix.ndim > 2:
            if gradient.ndim == 1:
                gradient = gradient[np.newaxis, :, np.newaxis]
            elif gradient.shape[:-2] == base_shape:
                gradient = gradient.reshape((-1,) + gradient.shape[-2:])
            elif gradient.shape[:-1] == base_shape:
                gradient = gradient.reshape((-1, gradient.shape[-1], 1))
            else:
                raise ValueError("can't make gradient of shape {} work with F-matrix of shape {}".format(
                    gradient.shape,
                    f_matrix.shape
                ))
        else:
            if gradient.ndim == 1:
                gradient = gradient[np.newaxis, :, np.newaxis]
            else:
                gradient = gradient[np.newaxis]

        f_matrix = f_matrix.reshape((-1,) + f_matrix.shape[-2:])
        if projector is not None:
            projector = np.asanyarray(projector)
            projector = projector.reshape((-1,) + projector.shape[-2:])

        if isinstance(mass_spec[0], str):  # atoms were supplied
            mass_spec = np.array([AtomData[a, "Mass"] for a in mass_spec]) * UnitsData.convert(
                "AtomicMassUnits",
                "AtomicUnitOfMass"
            )
        mass_spec = np.asanyarray(mass_spec)

        og_mass = mass_spec

        if mass_spec.ndim == 1:
            if mass_weighted is None: mass_weighted = dimensionless # generally do the right thing for Cartesians
            mass_spec = np.broadcast_to(mass_spec[:, np.newaxis], (len(mass_spec), 3)).flatten()
            mass_spec = np.diag(1 / mass_spec)

        mass_spec = np.reshape(mass_spec, (-1,) + mass_spec.shape[-2:])
        if mass_spec.shape[0] < f_matrix.shape[0]:
            mass_spec = np.broadcast_to(mass_spec, f_matrix.shape)

        gi12 = nput.fractional_power(mass_spec, -1 / 2)
        g12 = nput.fractional_power(mass_spec, 1 / 2)
        raw_grad = gradient
        gradient = g12 @ gradient

        grad_norm = nput.vec_norms(gradient, axis=1)
        if zero_gradient_cutoff is not None:
            if use_max_gradient_cutoff:
                regular_mode_pos = np.where(np.max(np.abs(raw_grad), axis=1) < zero_gradient_cutoff)
            else:
                regular_mode_pos = np.where(grad_norm < zero_gradient_cutoff)
            if len(regular_mode_pos) > 0:
                regular_mode_pos = regular_mode_pos[0]
            rem_pos = np.delete(np.arange(len(grad_norm)), regular_mode_pos)
        else:
            regular_mode_pos = []
            rem_pos = np.arange(len(f_matrix))

        if len(regular_mode_pos) > 0:

            reg_modes = cls.get_normal_modes(
                f_matrix[regular_mode_pos,],
                og_mass,
                # mass_units="AtomicMassUnits",
                remove_transrot=remove_transrot,
                dimensionless=dimensionless,
                mass_weighted=mass_weighted,
                zero_freq_cutoff=zero_freq_cutoff,
                # return_gmatrix=True,
                projector=projector[regular_mode_pos,] if projector is not None else projector
            )
        else:
            reg_modes = None

        # if zero_gradient_cutoff is not None and grad_norm < zero_gradient_cutoff:
        #     return None


        gradient = nput.vec_normalize(gradient[rem_pos,], norms=grad_norm[rem_pos,], axis=1)
        grad_projector = nput.orthogonal_projection_matrix(gradient, orthonormal=True)
        if projector is not None:
            projector = grad_projector @ projector[rem_pos,] @ grad_projector
        else:
            projector = grad_projector

        if len(rem_pos) > 0:
            f_matrix = f_matrix[rem_pos,]
            (freqs, modes, inv), g_matrix = cls.get_normal_modes(
                f_matrix,
                og_mass,
                # mass_units="AtomicMassUnits",
                remove_transrot=remove_transrot,
                dimensionless=dimensionless,
                mass_weighted=mass_weighted,
                zero_freq_cutoff=zero_freq_cutoff,
                return_gmatrix=True,
                projector=projector
            )


            gi12 = gi12[rem_pos,]
            g12 = g12[rem_pos,]
            gradient_inv = np.moveaxis(gradient, -1, -2)
            f_grad = gradient_inv @ f_matrix @ np.moveaxis(gradient_inv, -1, -2)
            grad_g = np.moveaxis(gradient, -1, -2) @ g_matrix @ gradient
            if not mass_weighted:
                gradient = gi12 @ gradient
                gradient_inv = gradient_inv @ g12

            if not remove_transrot:
                # if isinstance(freqs, list):
                    freqs_ = []
                    modes_ = []
                    inv_ = []
                    for g,f,m,i in zip(gradient, freqs, modes, inv):
                        zero_pos = np.where(np.abs(f) < 1e-6)[0] # projected out
                        overlaps = (i[zero_pos, :] @ g)**2
                        max_sim_pos = np.argmax(overlaps)
                        # print(max_sim_pos)
                        del_idx = zero_pos[max_sim_pos]
                        freqs_.append(np.delete(f, del_idx))
                        modes_.append(np.delete(m, del_idx, axis=-1))
                        inv_.append(np.delete(i, del_idx, axis=-2))
                    if isinstance(freqs, list):
                        freqs, modes, inv = freqs_, modes_, inv_
                    else:
                        freqs, modes, inv = np.array(freqs_), np.array(modes_), np.array(inv_)

            rows, cols = np.diag_indices(grad_g.shape[-1])
            loc_2 = f_grad[..., rows, cols] * grad_g[..., rows, cols]
            freq = np.sign(loc_2) * np.sqrt(np.abs(loc_2))

            if isinstance(freqs, list):
                freqs = [
                    np.concatenate([f1, f2], axis=0)
                    for f1, f2 in zip(freq, freqs)
                ]

                modes = [
                    np.concatenate([g, m], axis=1)
                    for g,m in zip(gradient, modes)
                ]

                inv = [
                    np.concatenate([gi, i], axis=0)
                    for gi,i in zip(gradient_inv, inv)
                ]

                if reg_modes is not None:
                    reg_freqs, reg_modes, reg_inv = reg_modes
                    full_freqs:list[np.ndarray] = [None] * (len(reg_freqs) + len(freqs))
                    full_modes:list[np.ndarray] = [None] * (len(reg_freqs) + len(freqs))
                    full_inv:list[np.ndarray] = [None] * (len(reg_freqs) + len(freqs))
                    for r,f,m,i in zip(regular_mode_pos, reg_freqs, reg_modes, reg_inv):
                        full_freqs[r] = f
                        full_modes[r] = m
                        full_inv[r] = i
                    for r,f,m,i in zip(rem_pos, freqs, modes, inv):
                        full_freqs[r] = f
                        full_modes[r] = m
                        full_inv[r] = i

                    freqs, modes, inv = full_freqs, full_modes, full_inv

                if len(base_shape) == 0:
                    freqs = freqs[0]
                    modes = modes[0]
                    inv = inv[0]
                # elif len(base_shape) > 1:
                #     ...
            else:
                freqs = np.concatenate([freq, freqs], axis=1)

                modes = np.concatenate(
                    [
                        gradient,
                        modes
                    ],
                    axis=2
                )

                inv = np.concatenate(
                    [
                        gradient_inv,
                        inv
                    ],
                    axis=1
                )

                if reg_modes is not None:
                    reg_freqs, reg_modes, reg_inv = reg_modes
                    full_freqs:np.ndarray = np.empty(
                        (len(reg_freqs) + len(freqs),) + freqs.shape[1:],
                        dtype=freqs.dtype
                    )
                    full_freqs[regular_mode_pos,] = reg_freqs
                    full_freqs[rem_pos,] = freqs

                    full_modes:np.ndarray = np.empty(
                        (len(reg_modes) + len(modes),) + modes.shape[1:],
                        dtype=modes.dtype
                    )
                    full_modes[regular_mode_pos,] = reg_modes
                    full_modes[rem_pos,] = modes


                    full_inv:np.ndarray = np.empty(
                        (len(reg_inv) + len(inv),) + inv.shape[1:],
                        dtype=inv.dtype
                    )
                    full_inv[regular_mode_pos,] = reg_inv
                    full_inv[rem_pos,] = inv

                    freqs, modes, inv = full_freqs, full_modes, full_inv

                freqs = freqs.reshape(base_shape + freqs.shape[1:])
                modes = modes.reshape(base_shape + modes.shape[1:])
                inv = inv.reshape(base_shape + inv.shape[1:])

        else:
            freqs, modes, inv = reg_modes

        mode_data = cls.ModeData(freqs, modes, inv)
        if return_gmatrix or return_indices:
            res = (mode_data,)
            if return_gmatrix:
                res = res + (mass_spec,)
            if return_indices:
                res = res + ((regular_mode_pos, rem_pos),)
            return res
        else:
            return mode_data

    @classmethod
    def from_grad_fg(
            cls,
            basis,
            gradient,
            f_matrix,
            mass_spec,
            remove_transrot=True,
            dimensionless=False,
            zero_freq_cutoff=None,
            mass_weighted=None,
            origin=None,
            projector=None,
            zero_gradient_cutoff=None,
            return_status=False,
            **opts
    ):
        """
        Generates normal modes from the specified F and G matrices

        :param basis:
        :param f_matrix: second derivatives of the potential
        :param mass_spec:
        :param mass_units:
        :param remove_transrot:
        :param opts:
        :return:
        """

        mode_data = cls.get_rp_modes(
            gradient, f_matrix, mass_spec,
            remove_transrot=remove_transrot,
            dimensionless=dimensionless,
            zero_freq_cutoff=zero_freq_cutoff,
            mass_weighted=mass_weighted,
            return_gmatrix=True,
            projector=projector,
            zero_gradient_cutoff=zero_gradient_cutoff,
            return_indices=True
        )
        if mode_data is None:
            return None

        (freqs, modes, inv), g_matrix, inds = mode_data

        mass_weighted = mass_weighted or dimensionless

        if mass_weighted and origin is not None:
            origin = np.asanyarray(origin)
            origin = (origin.flatten() @ nput.fractional_power(g_matrix, -1 / 2)).reshape(origin.shape)

        new = cls(basis, modes, inverse=inv, freqs=freqs,
                   mass_weighted=mass_weighted,
                   frequency_scaled=dimensionless,
                   g_matrix=g_matrix,
                   origin=origin,
                   **opts
                   )
        if return_status:
            reg, rp = inds
            return new, len(rp) > 0
        else:
            return new

    @classmethod
    def from_molecule(cls, mol,
                      dimensionless=False,
                      use_internals=None,
                      potential_derivatives=None,
                      project_transrot=True,
                      zero_freq_cutoff=None,
                      masses=None,
                      zero_gradient_cutoff=None,
                      return_status=False,
                      **opts
                      ):
        from ..Molecools import Molecule
        #TODO: cut down on duplication

        mol = mol  # type:Molecule
        if use_internals is None:
            use_internals = mol.internal_coordinates is not None
        if potential_derivatives is None:
            potential_derivatives = mol.potential_derivatives
            if potential_derivatives is None and mol.energy_evaluator is not None:
                potential_derivatives = mol.calculate_energy(order=2)[1:]
        if use_internals:
            hess = nput.tensor_reexpand(
                mol.get_cartesians_by_internals(1, strip_embedding=True, reembed=True),
                [0, potential_derivatives[1]],
                order=2
            )[1]
            if zero_freq_cutoff is None:
                zero_freq_cutoff = cls.default_projected_zero_freq_cutoff
                if zero_freq_cutoff is None:
                    zero_freq_cutoff = 1 * UnitsData.convert("Wavenumbers", "Hartrees")
            return cls.from_grad_fg(
                mol.internal_coordinates.system,
                potential_derivatives[0],
                hess,
                mol.get_gmatrix(masses=masses),
                dimensionless=dimensionless,
                origin=mol.get_internals(strip_embedding=True),
                zero_freq_cutoff=zero_freq_cutoff,
                zero_gradient_cutoff=zero_gradient_cutoff,
                return_status=return_status,
                **opts
            )
        else:
            if not project_transrot and mol.normal_modes.modes is not None:
                return mol.normal_modes.modes.basis.to_new_modes()

            if masses is None:
                masses = mol.atomic_masses
            hess = potential_derivatives[1]
            if project_transrot:
                proj = mol.get_translation_rotation_projector(mass_weighted=True)
                mw = np.diag(np.repeat(1 / np.sqrt(masses), 3))
                mw_inv = np.diag(np.repeat(np.sqrt(masses), 3))
                proj = mw_inv @ proj @ mw
                hess = proj @ hess @ proj.T
                if zero_freq_cutoff is None:
                    zero_freq_cutoff = cls.default_projected_zero_freq_cutoff
                    if zero_freq_cutoff is None:
                        zero_freq_cutoff = 1 * UnitsData.convert("Wavenumbers", "Hartrees")

            return cls.from_grad_fg(
                mol.coords.system,
                potential_derivatives[0],
                hess,
                masses,
                dimensionless=dimensionless,
                masses=masses,
                origin=mol.coords,
                zero_freq_cutoff=zero_freq_cutoff,
                zero_gradient_cutoff=zero_gradient_cutoff,
                return_status=return_status,
                **opts
            )
