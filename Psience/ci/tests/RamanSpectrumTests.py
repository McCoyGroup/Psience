from pathlib import Path
import tempfile
import unittest

import numpy as np
import Psience
from Psience.Modes import NormalModes
from Psience.Molecools import Molecule
from Psience.Spectra import DiscreteSpectrum, HarmonicSpectrum
from McUtils.Coordinerds import CartesianCoordinates3D
from McUtils.Data import UnitsData


class RamanSpectrumTests(unittest.TestCase):
    @staticmethod
    def modes():
        return NormalModes(CartesianCoordinates3D, np.eye(3), inverse=np.eye(3),
                           freqs=np.array([.01, .02, .03]), masses=np.ones(1),
                           origin=np.zeros((1, 3)))

    @staticmethod
    def tensors():
        return np.array([np.eye(3), np.diag([1., -1., 0.]),
                         [[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]]])

    def test_harmonic_activities_and_units(self):
        # Isotropic I -> S=45; traceless diag(1,-1,0) and xy -> S=21.
        spec = HarmonicSpectrum.raman_from_modes(self.modes(), self.tensors())
        # Independent SI dimensional conversion: a0^4 / electron mass to A^4/u.
        conversion = (.529177210903 ** 4) * (1.66053906660e-27 / 9.1093837015e-31)
        np.testing.assert_allclose(spec.intensities, np.array([45, 21, 21]) * conversion)
        self.assertEqual(spec.meta["intensity_unit"], "Angstroms**4/AtomicMassUnits")
        self.assertEqual(spec.meta["frequency_unit"], "Wavenumbers")
        for modes in (self.modes().make_frequency_scaled(), self.modes().make_mass_weighted(),
                      self.modes().make_dimensionless()):
            np.testing.assert_allclose(HarmonicSpectrum.raman_from_modes(modes, self.tensors()).intensities,
                                       spec.intensities)

    def test_rotations_and_complex_phase_preserve_strengths(self):
        angle = .49
        rot = np.array([[np.cos(angle), -np.sin(angle), 0],
                        [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        tensors = self.tensors()
        freqs = self.modes().freqs
        expected = DiscreteSpectrum.from_raman_moments(freqs, tensors)
        rotated = DiscreteSpectrum.from_raman_moments(freqs, rot @ tensors @ rot.T)
        phased = DiscreteSpectrum.from_raman_moments(freqs, tensors * np.exp(.7j))
        np.testing.assert_allclose(rotated.intensities, expected.intensities)
        np.testing.assert_allclose(phased.intensities, expected.intensities)
        np.testing.assert_allclose(expected.intensities, np.array([45, 21, 21]) * .529177210903**6)
        self.assertEqual(expected.meta["intensity_unit"], "Angstroms**6")

    def test_absolute_cross_section_and_stokes_sign(self):
        # An isotropic transition tensor I has d sigma/d Omega = k_s^4 a0^6.
        shifts = np.array([.01, -.01, .10, .11])
        spec = DiscreteSpectrum.from_raman_moments(shifts, np.tile(np.eye(3), (4, 1, 1)), pump_frequency=.1)
        expected = (2*np.pi * np.array([.09, .11, 0., 0.]) * 219474.6313632)**4 * (5.29177210903e-9)**6
        np.testing.assert_allclose(spec.intensities, expected, rtol=1e-9)
        self.assertEqual(spec.meta["intensity_unit"], "Centimeters**2/Steradians")
        self.assertEqual(UnitsData.unit_label(spec.meta["intensity_unit"]), "cm²/sr")
        self.assertLess(spec.intensities[0], spec.intensities[1])

    def test_harmonic_cross_section_contains_zero_point_matrix_element(self):
        modes = self.modes()
        spec = HarmonicSpectrum.raman_from_modes(modes, self.tensors(), pump_frequency=.1)
        moments = self.tensors() / np.sqrt(2*modes.freqs)[:, None, None]
        expected = DiscreteSpectrum.from_raman_moments(modes.freqs, moments, pump_frequency=.1)
        np.testing.assert_allclose(spec.intensities, expected.intensities)

    def test_invalid_inputs(self):
        for pump in (-1, float("nan"), float("inf"), 1j, [1, 2]):
            with self.subTest(pump=pump), self.assertRaises(ValueError):
                DiscreteSpectrum.from_raman_moments([.01], [np.eye(3)], pump_frequency=pump)
        for tensors in (np.ones((1, 6)), np.full((1, 3, 3), np.nan),
                        np.array([[[0, 1, 0], [0, 0, 0], [0, 0, 0]]])):
            with self.subTest(shape=tensors.shape), self.assertRaises(ValueError):
                DiscreteSpectrum.from_raman_moments([.01], tensors)
        with self.assertRaises(ValueError):
            DiscreteSpectrum.from_raman_moments([.01, .02], [np.eye(3)])
        with self.assertRaises(ValueError):
            HarmonicSpectrum.raman_from_modes(self.modes().modify(freqs=np.array([0., .02, .03])),
                                             self.tensors(), pump_frequency=.1)

    def test_units_follow_spectral_operations(self):
        spec = HarmonicSpectrum.raman_from_modes(self.modes(), self.tensors())
        for result in (spec[:2], spec.broaden(), spec.intensity_filter(0, np.inf)):
            self.assertEqual(result.meta["intensity_unit"], spec.meta["intensity_unit"])
        self.assertIsNone(spec.normalize().meta["intensity_unit"])
        ir = DiscreteSpectrum.from_transition_moments(np.array([.01]), np.array([[1., 0., 0.]]))
        self.assertAlmostEqual(ir.intensities[0], 35542.06329390961)
        self.assertEqual(ir.meta["intensity_unit"], "KilometersPerMole")
        self.assertEqual(UnitsData.unit_label(spec.meta["intensity_unit"]), "Å⁴/amu")

    def test_water_molecule_from_gaussian_numerical_derivatives(self):
        data = Path(Psience.__file__).resolve().parent.parent / "ci/tests/TestData/water_freq.fchk"
        mol = Molecule.from_file(str(data))
        spec = mol.get_harmonic_raman_spectrum()
        # Independently unpack the first 18 numbers of the checkpoint's numerical
        # polarizability block (six components per mode), in sqrt(amu)*bohr.
        np.testing.assert_allclose(spec.frequencies, [1622.30295, 3803.30058, 3937.52557], atol=2e-5)
        np.testing.assert_allclose(spec.intensities, [1.95602618, 104.09028427, 24.17356499], rtol=2e-6)
        np.testing.assert_allclose(mol.get_harmonic_spectrum().intensities,
                                   [67.4560596, 4.14281444, 67.0202811], rtol=5e-6)
        internal_mol = Molecule.from_file(str(data), internals=[[0, -1, -1, -1],
                                                              [1, 0, -1, -1], [2, 0, 1, -1]])
        np.testing.assert_allclose(internal_mol.get_harmonic_raman_spectrum().intensities, spec.intensities)
        amu_to_me = 1822.888486209
        self.assertAlmostEqual(mol.polarizability_derivatives[2][0, 0, 0, 0], .918608324 / amu_to_me, places=11)

    def test_constant_only_expansion_requests_the_missing_first_derivative(self):
        from unittest.mock import patch
        mol = Molecule(["O", "H", "H"], [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
                       polarizability_derivatives=[np.eye(3)])
        derivatives = [np.eye(3), np.ones((9, 3, 3))]
        with patch.object(mol, "calculate_dipole_polarizability", return_value=(None, derivatives)) as calculate:
            actual = mol.get_cartesian_polarizability_derivatives(1, evaluator=object())
        calculate.assert_called_once()
        self.assertEqual(calculate.call_args.kwargs["order"], 1)
        np.testing.assert_array_equal(actual[1], derivatives[1])

    def test_gaussian_cartesian_derivatives_and_missing_derivatives(self):
        # This scalar polarizability plus Cartesian derivative block must load
        # without requiring numerical normal-coordinate derivatives.
        from Psience.Molecools.Properties import DipoleSurfaceManager
        values = np.arange(54.)
        text = "title\nmethod\n" + f"{'Number of atoms':<43}I                3\n"
        text += f"{'Polarizability':<43}R   N=           6\n"
        text += "".join(f"{v:16.8E}" for v in [1, 0, 1, 0, 0]) + "\n" + f"{1:16.8E}\n"
        text += f"{'Polarizability Derivatives':<43}R   N=          54\n"
        text += "\n".join("".join(f"{v:16.8E}" for v in values[i:i+5]) for i in range(0, 54, 5)) + "\n"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "analytic.fchk"
            path.write_text(text)
            derivs = DipoleSurfaceManager(None).load_polarizability_derivatives(str(path))
        self.assertEqual(derivs[1].shape, (9, 3, 3))
        np.testing.assert_array_equal(derivs[1][0], [[0, 1, 3], [1, 2, 4], [3, 4, 5]])
        mol = Molecule(["O", "H", "H"], [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
                       polarizability_derivatives=[np.eye(3)])
        with self.assertRaisesRegex(ValueError, "first Cartesian polarizability derivatives"):
            mol.get_harmonic_raman_spectrum()
