
import numpy as np, scipy.linalg as slag, itertools

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

    default_zero_freq_cutoff = 2.5e-4 # 50 wavenumbers...
    @classmethod
    def get_normal_modes(cls,
                         f_matrix,
                         mass_spec,
                         # mass_units="AtomicMassUnits",
                         remove_transrot=True,
                         dimensionless=False,
                         mass_weighted=None,
                         zero_freq_cutoff=None,
                         mode='default'
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
        if remove_transrot:
            gi12 = slag.fractional_matrix_power(mass_spec, -1 / 2)
            g12 = slag.fractional_matrix_power(mass_spec, 1 / 2)
            V = gi12 @ modes
            # modes = g12 @ V
            # therefore, inv = V.T @ gi12
            # and modes[:, nz] = g12 @ V[:, nz]; inv[nz, :] = (V.T)[nz, :] @ gi12 (this is only the _left_ inverse)
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
            if np.linalg.det(modes) == 1 and np.allclose(modes.T @ modes, np.eyelen(modes)):
                inv = modes.T
            else:
                inv = np.linalg.inv(modes)

        freqs = np.sign(freq2) * np.sqrt(np.abs(freq2))
        if mode == 'reasonable':
            modes, inv = inv.T, modes.T
            if mass_weighted:
                gi12 = slag.fractional_matrix_power(mass_spec, -1 / 2)
                g12 = slag.fractional_matrix_power(mass_spec, 1 / 2)
                inv = inv @ gi12
                modes = g12 @ modes
            if dimensionless:
                weighting = np.sqrt(np.abs(freqs))
                modes = modes * weighting[np.newaxis, :]
                inv = inv / weighting[:, np.newaxis]
        else:
            if mass_weighted:
                gi12 = slag.fractional_matrix_power(mass_spec, -1 / 2)
                g12 = slag.fractional_matrix_power(mass_spec, 1 / 2)
                inv = inv @ g12
                modes = gi12 @ modes
            if dimensionless:
                weighting = np.sqrt(np.abs(freqs))
                modes = modes / weighting[np.newaxis, :]
                inv = inv * weighting[:, np.newaxis]

        sorting = np.argsort(freqs)
        modes = modes[:, sorting]
        inv = inv[sorting, :]

        return freqs, modes, inv
    @classmethod
    def from_fg(cls,
                basis,
                f_matrix,
                mass_spec,
                remove_transrot=True,
                dimensionless=False,
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

        freqs, modes, inv = cls.get_normal_modes(f_matrix, mass_spec,
                                          remove_transrot=remove_transrot,
                                          dimensionless=dimensionless)

        return cls(basis, modes, inverse=inv, freqs=freqs, **opts)

    def get_nearest_mode_transform(self, alternate_modes):
        modes = self.matrix

        ls_coords = np.linalg.inv(modes.T @ modes) @ modes.T @ alternate_modes

        U, s, V = np.linalg.svd(ls_coords)
        good_s = np.where(s > 1e-8)[0]
        R = (U[:, good_s] @ V[good_s, :])

        return R

    def get_localized_modes(self, target_coords,
                            symmetrizer=None,
                            mode_selection=None,
                            make_oblique=False,
                            rediagonalize=False
                            ):
        from .ObliqueModes import ObliqueModeGenerator

        carts = self.origin
        if carts.ndim != 2:
            raise NotImplementedError("can't get local modes for non-Cartesian basis")

        targets = []
        for idx in target_coords:
            nidx = len(idx)
            if nidx == 2:
                targets.append(nput.dist_vec(carts, *idx))
            elif nidx == 3:
                targets.append(nput.angle_vec(carts, *idx))
            elif nidx == 4:
                targets.append(nput.dihed_vec(carts, *idx))
            else:
                raise ValueError("can't parse coordinate spec {}".format(idx))
        target_modes = np.array(targets).T

        if symmetrizer is not None:
            if isinstance(symmetrizer, (np.ndarray, list)):
                symmetrizer = np.asanyarray(symmetrizer)
                target_modes = target_modes @ symmetrizer.T
            else:
                target_modes = symmetrizer(target_modes)

        ltf = self.get_nearest_mode_transform(target_modes)
        ltf_inv = ltf.T

        if make_oblique:
            new_f = ltf.T @ np.diag(self.freqs ** 2) @ ltf
            new_g = np.eye(len(new_f))
            ob_f, ob_g, ob_u, ob_ui = ObliqueModeGenerator(new_f, new_g).run()
            if mode_selection is not None:
                ltf = ltf @ ob_u[:, mode_selection]
                ltf_inv = ob_ui[mode_selection, :] @ ltf_inv
            else:
                ltf = ltf @ ob_u
                ltf_inv = ob_ui @ ltf_inv
        elif mode_selection is not None:
            ltf = ltf[:, mode_selection]
            ltf_inv = ltf.T

        if rediagonalize:
            new_f = ltf_inv @ np.diag(self.freqs ** 2) @ ltf_inv.T
            new_freq2, mode_tf = np.linalg.eigh(new_f)
            mode_inv = mode_tf.T
            if not make_oblique:
                new_freqs = np.sqrt(new_freq2)
            else:
                new_freqs = new_freq2
                mode_tf = mode_tf @ np.diag(1/np.sqrt(new_freqs))
                mode_inv = np.diag(np.sqrt(new_freqs)) @ mode_inv

            full_tf = ltf @ mode_tf
            full_inv = mode_inv @ ltf_inv
            new_modes = self.matrix @ full_tf
            new_inv = full_inv @ self.inverse

            return (full_tf, full_inv), type(self)(
                self.basis,
                new_modes,
                inverse=new_inv,
                origin=self.origin,
                freqs=new_freqs
            )
        else:
            return ltf

    def make_dimensionless(self, masses):
        L = self.matrix.T
        Linv = self.inverse
        freqs = self.freqs
        conv = np.sqrt(np.broadcast_to(freqs[:, np.newaxis], L.shape))
        if masses is not None:
            trip_mass = np.broadcast_to(
                np.sqrt(masses)[:, np.newaxis],
                (len(masses), L.shape[1]//len(masses))
            ).flatten()
            mass_conv = np.broadcast_to(trip_mass[np.newaxis, :], L.shape)
            conv = conv * mass_conv
        L = L * conv
        Linv = Linv / conv
        return type(self)(self.basis, L.T, inverse=Linv, origin=self.origin, freqs=self.freqs)

    @classmethod
    def from_molecule(cls, mol, dimensionless=False, use_internals=None):
        from ..Molecools import Molecule

        mol = mol  # type:Molecule
        if use_internals is None:
            use_internals = mol.internal_coordinates is not None
        if use_internals:
            return cls.from_fg(
                mol.internal_coordinates.system,
                mol.get_internal_potential_derivatives(2)[1],
                mol.g_matrix,
                dimensionless=dimensionless,
                origin=mol.coords
            )
        else:
            return cls.from_fg(
                mol.coords.system,
                mol.potential_derivatives[1],
                mol.atomic_masses,
                dimensionless=dimensionless,
                origin=mol.coords
            )

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