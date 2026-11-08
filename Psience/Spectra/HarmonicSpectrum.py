
import numpy as np
from McUtils.Data import UnitsData
from .BaseSpectrum import DiscreteSpectrum

__all__ = [
    "HarmonicSpectrum"
]

class HarmonicSpectrum(DiscreteSpectrum):

    @classmethod
    def from_normal_modes(cls, nms, dipole_derivatives, **opts):
        """
        **LLM Docstring**

        Build a harmonic (double-harmonic-approximation) IR spectrum from a set of normal modes and their Cartesian dipole derivatives, converting the modes to a non-mass-weighted, frequency-scaled basis and computing the transition moments from the dipole derivatives projected onto the modes.

        :param nms: the normal modes to build the spectrum from
        :type nms: MixtureModes
        :param dipole_derivatives: the Cartesian dipole-derivative tensor (first derivatives of the dipole with respect to Cartesian displacements)
        :type dipole_derivatives: np.ndarray
        :param opts: extra options forwarded to `from_transition_moments`
        :type opts: dict
        :return: the constructed harmonic spectrum
        :rtype: HarmonicSpectrum
        """
        nms = nms.remove_mass_weighting().make_frequency_scaled()
        freqs = nms.freqs
        tms = 1/np.sqrt(2) * nms.coords_by_modes @ dipole_derivatives
        return cls.from_transition_moments(
            freqs,
            tms,
            **opts
        )

    @classmethod
    def raman_from_modes(cls, nms, polarizability_derivatives, pump_frequency=0, **opts):
        """Build harmonic Raman activities, or laser-dependent cross sections.

        Cartesian polarizability derivatives are in atomic units (a0^2).
        With no pump, intensities are the conventional activities
        S = 45 (d alpha_iso/dQ)^2 + 7 (d gamma/dQ)^2 in Angstroms**4/amu,
        where Q is a mass-weighted, frequency-unscaled normal coordinate.

        A positive pump photon energy in Hartrees instead returns ground-state
        Stokes differential cross sections in cm^2/sr, using the geometry and
        convention documented by DiscreteSpectrum.from_raman_moments.
        """
        nms = nms.remove_mass_weighting().remove_frequency_scaling()
        freqs = np.asarray(nms.freqs)
        if np.iscomplexobj(freqs) or not np.all(np.isfinite(freqs)):
            raise ValueError("harmonic Raman frequencies must be finite and real")
        pol_derivs = np.tensordot(nms.coords_by_modes, polarizability_derivatives, axes=[-1, 0])
        if pump_frequency is None or (np.isscalar(pump_frequency) and pump_frequency == 0):
            isotropic, anisotropic = cls._raman_invariants(pol_derivs)
            units = (UnitsData.convert("BohrRadius", "Angstroms") ** 4
                     / UnitsData.convert("ElectronMass", "AtomicMassUnits"))
            opts.update(frequency_unit="Wavenumbers", intensity_unit="Angstroms**4/AtomicMassUnits",
                        intensity_label="Raman activity")
            return cls(freqs * UnitsData.convert("Hartrees", "Wavenumbers"),
                       units * (45 * isotropic + 7 * anisotropic), **opts)
        if np.any(freqs <= 0):
            raise ValueError("harmonic Raman cross sections require positive vibrational frequencies")
        pol_moms = pol_derivs / np.sqrt(2 * freqs)[..., np.newaxis, np.newaxis]
        return cls.from_raman_moments(
            freqs,
            pol_moms,
            pump_frequency=pump_frequency,
            **opts
        )

    @classmethod
    def from_mol(cls, mol, modes=None, dipole_derivatives=None, **opts):
        """
        **LLM Docstring**

        Build a harmonic IR spectrum for a molecule, using its normal modes and Cartesian dipole derivatives (computed from the molecule if not given directly), via `from_normal_modes`.

        :param mol: the molecule to build the spectrum for
        :type mol: Molecule
        :param modes: the normal modes to use; computed via `mol.get_normal_modes()` if not given
        :type modes: MixtureModes | None
        :param dipole_derivatives: the dipole derivatives to use; if given as a full derivative expansion, the first-order term (`[1]`) is used, otherwise computed via `mol.get_cartesian_dipole_derivatives(1)[0]`
        :type dipole_derivatives: object | None
        :param opts: extra options forwarded to `from_normal_modes`
        :type opts: dict
        :return: the constructed harmonic spectrum
        :rtype: HarmonicSpectrum
        """
        return cls.from_normal_modes(
            mol.get_normal_modes() if modes is None else modes,
            mol.get_cartesian_dipole_derivatives(1)[0] if dipole_derivatives is None else dipole_derivatives[1],
            **opts
        )

    @classmethod
    def raman_from_mol(cls, mol, modes=None, polarizability_derivatives=None, **opts):
        """
        **LLM Docstring**

        Build harmonic Raman activities (or cross sections when a pump is given)
        from the molecule's Cartesian modes and polarizability derivatives.

        :param mol: the molecule to build the spectrum for
        :type mol: Molecule
        :param modes: the Cartesian normal modes to use; computed with `use_internals=False` if not given
        :type modes: MixtureModes | None
        :param polarizability_derivatives: the full polarizability-derivative expansion to use (the first-order term is used); computed via `mol.get_cartesian_polarizability_derivatives(1)` if not given
        :type polarizability_derivatives: list[np.ndarray] | None
        :param opts: extra options forwarded to `raman_from_modes`
        :type opts: dict
        :return: the constructed harmonic Raman spectrum
        :rtype: HarmonicSpectrum
        """
        if polarizability_derivatives is None:
            polarizability_derivatives = mol.get_cartesian_polarizability_derivatives(1)
        if polarizability_derivatives is None or len(polarizability_derivatives) < 2:
            raise ValueError("harmonic Raman spectra require first Cartesian polarizability derivatives")
        return cls.raman_from_modes(
            mol.get_normal_modes(use_internals=False) if modes is None else modes,
            polarizability_derivatives[1],
            **opts
        )