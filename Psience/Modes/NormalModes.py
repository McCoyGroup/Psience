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
                         return_gmatrix=False
                         ):

        f_matrix = np.asanyarray(f_matrix)
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
        freq2, modes = slag.eigh(f_matrix, mass_spec, type=3)

        gi12 = nput.fractional_power(mass_spec, -1 / 2)
        g12 = nput.fractional_power(mass_spec, 1 / 2)
        V = gi12 @ modes
        # modes = g12 @ V
        # therefore, inv = V.T @ gi12
        # and modes[:, nz] = g12 @ V[:, nz]; inv[nz, :] = (V.T)[nz, :] @ gi12 (this is only the _left_ inverse)

        if remove_transrot:
            if zero_freq_cutoff is None:
                zero_freq_cutoff = cls.default_zero_freq_cutoff
            nonzero = np.abs(freq2) > zero_freq_cutoff**2
        # if len(nonzero) < len(modes):
        #     if np.linalg.det(modes) != 1 or not np.allclose(modes.T @ modes, np.eyelen(modes)):
        #         # we're working in a non-invertible subspace...
        #         raise ValueError("non-invertible subspace of normal modes found, inverse can't be evaluated")
        #     else:
        #         inv = modes.T
        # else:
        #     inv = np.linalg.inv(modes)
        # nonzero = np.array([0, 3, 4, 5])
            freq2 = freq2[nonzero]
            modes = modes[:, nonzero]
            inv = V.T[nonzero, :] @ gi12
        else:
            inv = V.T @ gi12
        # else:
        #     if np.linalg.det(modes) == 1 and np.allclose(modes.T @ modes, np.eyelen(modes)):
        #         inv = modes.T
        #     else:
        #         inv = np.linalg.inv(modes)

        freqs = np.sign(freq2) * np.sqrt(np.abs(freq2))
        if mass_weighted:
            # gi12 = slag.fractional_matrix_power(mass_spec, -1 / 2)
            # g12 = slag.fractional_matrix_power(mass_spec, 1 / 2)
            inv = inv @ g12
            modes = gi12 @ modes
        if dimensionless:
            weighting = np.sqrt(np.abs(freqs))
            modes = modes / weighting[np.newaxis, :]
            inv = inv * weighting[:, np.newaxis]
        # if mode == 'reasonable':
        # modes, inv = inv.T, modes.T

        sorting = np.argsort(freqs)
        modes = modes[:, sorting]
        inv = inv[sorting, :]

        mode_data = cls.ModeData(freqs, inv.T, modes.T)
        if return_gmatrix:
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
                                                 return_gmatrix=True
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

    @classmethod
    def _atom_projector(cls, n, i):
        if nput.is_numeric(i):
            i = [i]
        z = np.zeros((3 * n, 3 * n))
        for i in i:
            x = np.arange(3 * i, 3 * (i + 1))
            z[x, x] = 1
        return z
    def get_nearest_mode_transform(self,
                                   alternate_modes:np.ndarray,
                                   mass_weighted=False,
                                   atoms=None,
                                   maximum_similarity=True,
                                   unitarize=True
                                   ):
        if not mass_weighted:
            modes = self.remove_mass_weighting()
        else:
            modes = self.make_mass_weighted()

        inv = modes.inverse
        modes = modes.matrix
        if atoms is not None:
            proj = self._atom_projector(modes.shape[0], atoms)
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

    def get_atom_localized_mode_transformation(self,
                                               atoms,
                                               masses=None, origin=None,
                                               localization_type='ned',
                                               allow_mode_mixing=False,
                                               maximum_similarity=False,
                                               unitarize=True
                                               ):
        if masses is None:
            masses = self.masses
            mw = self.make_mass_weighted()
        else:
            mw = self.make_mass_weighted(masses=masses)

        if origin is None:
            origin = mw.remove_mass_weighting().origin

        if nput.is_numeric(atoms):
            atoms = [atoms]

        f = mw.matrix @ np.diag(self.freqs**2) @ mw.matrix.T

        if localization_type == 'direct':
            tr_tf, tr_inv = nput.translation_rotation_invariant_transformation(
                origin.reshape(-1, 3),
                masses,
                mass_weighted=True
            )
        else:
            tr_tf = tr_inv = np.eye(mw.matrix.shape[0])

        g_matrix = tr_tf.T @ tr_tf
        if allow_mode_mixing:
            proj = self._atom_projector(len(masses), atoms)
            freqs, modes, inv = self.get_normal_modes(
                tr_inv @ proj @ f @ proj @ tr_inv.T,
                g_matrix,
                remove_transrot=True
            )
        else:
            # freqs = []
            modes = []
            # inv = []
            for atom in atoms:
                proj = self._atom_projector(len(masses), atom)
                sub_freqs, sub_modes, sub_inv = self.get_normal_modes(
                    tr_inv @ proj @ f @ proj @ tr_inv.T,
                    g_matrix,
                    remove_transrot=True
                )
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

    def get_internal_localized_mode_transformation(
            self,
            expansion_coordinates: "Iterable[Iterable[int]|dict]",
            fixed_atoms=None,
            mass_weighted=False,
            project_transrot=True,
            atoms=None,
            maximum_similarity=False,
            unitarize=True
    ):
        # from .ObliqueModes import ObliqueModeGenerator

        nmw = self.remove_mass_weighting()
        base_derivs = nput.internal_coordinate_tensors(
            nmw.origin,
            expansion_coordinates,
            fixed_atoms=fixed_atoms
        )

        if mass_weighted:
            base_derivs = np.diag(np.repeat(1/np.sqrt(self.masses), 3)) @ base_derivs

        if project_transrot:
            origin = self.origin
            masses = self.masses
            projector = nput.translation_rotation_projector(
                    origin.reshape(-1, 3),
                    masses,
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

    @property
    def localizer_dispatch(self):
        return {
            self.LocalizationMethods.MaximumSimilarity.value:(self.get_nearest_mode_transform, 'target_modes'),
            self.LocalizationMethods.AtomLocalized.value:(self.get_atom_localized_mode_transformation, 'atoms'),
            self.LocalizationMethods.DisplacmentMinimized.value:(self.get_displacement_localized_mode_transformation, 'mode_blocks'),
            self.LocalizationMethods.Internals.value:(self.get_internal_localized_mode_transformation, 'internals')
        }

    def localize(self,
                 method=None,
                 *,
                 atoms=None,
                 target_modes=None,
                 internals=None,
                 mode_blocks=None,
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
    def get_rp_modes(cls,
                  gradient,
                  f_matrix,
                  mass_spec,
                  mass_units="AtomicMassUnits",
                  remove_transrot=True,
                  dimensionless=False
                  ):
        ...