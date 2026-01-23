
import numpy as np
from .BaseSpectrum import DiscreteSpectrum

__all__ = [
    "HarmonicSpectrum"
]

class HarmonicSpectrum(DiscreteSpectrum):

    @classmethod
    def from_normal_modes(cls, nms, dipole_derivatives, **opts):
        nms = nms.remove_mass_weighting().make_frequency_scaled()
        freqs = nms.freqs
        tms = 1/np.sqrt(2) * nms.coords_by_modes @ dipole_derivatives
        return cls.from_transition_moments(
            freqs,
            tms,
            **opts
        )

    @classmethod
    def raman_from_modes(cls, nms, polarizability_derivatives, **opts):
        nms = nms.remove_mass_weighting().make_frequency_scaled()
        freqs = nms.freqs
        pol_moms = 1/np.sqrt(2) * np.tensordot(nms.coords_by_modes, polarizability_derivatives, axes=[-1, 0])
        return cls.from_raman_moments(
            freqs,
            pol_moms,
            **opts
        )

    @classmethod
    def from_mol(cls, mol, modes=None, dipole_derivatives=None, **opts):
        return cls.from_normal_modes(
            mol.get_normal_modes() if modes is None else modes,
            mol.get_cartesian_dipole_derivatives(1)[0] if dipole_derivatives is None else dipole_derivatives[1],
            **opts
        )

    @classmethod
    def raman_from_mol(cls, mol, modes=None, polarizability_derivatives=None, **opts):
        if polarizability_derivatives is None:
            polarizability_derivatives = mol.get_cartesian_polarizability_derivatives(1)
        return cls.raman_from_modes(
            mol.get_normal_modes() if modes is None else modes,
            polarizability_derivatives[1],
            **opts
        )