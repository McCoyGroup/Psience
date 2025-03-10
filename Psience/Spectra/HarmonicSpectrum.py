
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
    def from_mol(cls, mol, **opts):
        return cls.from_normal_modes(
            mol.get_normal_modes(),
            mol.get_cartesian_dipole_derivatives(1)[0],
            **opts
        )