
import numpy as np, scipy.linalg as slag
import scipy.linalg

from McUtils.Data import AtomData, UnitsData

from .MixtureModes import MixtureModes

__all__ = [
    "NormalModes",
    "ReactionPathModes"
]

__reload_hook__ = [".MixtureModes"]

class NormalModes(MixtureModes):
    name="NormalModes"

    default_zero_freq_cutoff = 4.5e-5 # 10 wavenumbers...
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
            nonzero = np.abs(freq2) > zero_freq_cutoff
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