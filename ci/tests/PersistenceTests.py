"""
Tests for persisting Psience objects through `McUtils.Scaffolding` checkpointers.

Plain `unittest` tests that only write inside temp directories. Gaps that aren't fixed yet are
marked `@unittest.expectedFailure` with the phase of the persistence plan
(`claude_drafts/persistence/PLAN.md`) that fixes them; remove the marker in the patch that lands the fix.
"""

import os, shutil, tempfile, unittest
import numpy as np

import warnings
from McUtils.Scaffolding import Checkpointer, PseudoPickler
from McUtils.Coordinerds import CoordinateSet, CartesianCoordinates3D
from Psience.Molecools import Molecule
from Psience.Modes import NormalModes
from Psience.BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from Psience.Molecools.Serializers import MoleculePropertyCache

__all__ = ["PersistenceTests"]

WATER = (['O', 'H', 'H'], np.array([[0., 0., 0.], [0.96, 0., 0.], [-0.24, 0.93, 0.]]))


class PersistenceTests(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.mkdtemp(prefix="psience_persistence_tests_")

    def tearDown(self):
        shutil.rmtree(self._tmp, ignore_errors=True)

    def path(self, name):
        return os.path.join(self._tmp, name)

    #region Molecule round trips
    def _roundtrip_molecule(self, ext):
        mol = Molecule(*WATER)
        fn = self.path('mol' + ext)
        with Checkpointer.from_file(fn) as chk:
            chk['mol'] = mol
        with Checkpointer.from_file(fn) as chk:
            loaded = chk['mol']
        self.assertIsInstance(loaded, Molecule)
        self.assertEqual(tuple(loaded.atoms), tuple(mol.atoms))
        self.assertTrue(np.allclose(loaded.coords, mol.coords))
        self.assertTrue(np.allclose(loaded.masses, mol.masses))

    def test_MoleculeRoundTripJSON(self):
        self._roundtrip_molecule('.json')

    def test_MoleculeRoundTripNPZ(self):
        self._roundtrip_molecule('.npz')

    def test_MoleculeRoundTripHDF5(self):
        self._roundtrip_molecule('.hdf5')
    #endregion

    #region Object round trips through every checkpoint format
    EXTENSIONS = ['.json', '.npz', '.hdf5']

    def _roundtrip(self, name, value, check):
        for ext in self.EXTENSIONS:
            fn = self.path(name + ext)
            with Checkpointer.from_file(fn) as chk:
                chk['obj'] = value
            with Checkpointer.from_file(fn) as chk:
                loaded = chk['obj']
            check(loaded, ext)

    def test_MoleculeWithInternalsRoundTrip(self):
        mol = Molecule(*WATER, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        def check(loaded, ext):
            self.assertIsInstance(loaded, Molecule, ext)
            self.assertEqual(tuple(loaded.atoms), tuple(mol.atoms), ext)
            self.assertTrue(np.allclose(loaded.coords, mol.coords), ext)
            self.assertIsInstance(loaded.coords, CoordinateSet, ext)
            self.assertEqual(np.asarray(loaded.internals['zmatrix']).tolist(), np.asarray(mol.internals['zmatrix']).tolist(), ext)
        self._roundtrip('mol_int', mol, check)

    def test_NormalModesRoundTrip(self):
        modes = NormalModes(None, np.eye(3), freqs=np.array([1., 2., 3.]))
        def check(loaded, ext):
            self.assertIsInstance(loaded, NormalModes, ext)
            self.assertTrue(np.allclose(loaded.matrix, modes.matrix), ext)
            self.assertTrue(np.allclose(loaded.freqs, modes.freqs), ext)
        self._roundtrip('modes', modes, check)

    def test_HarmonicOscillatorBasisAndStateSpaceRoundTrip(self):
        basis = HarmonicOscillatorProductBasis(3)
        space = BasisStateSpace(basis, [[0, 0, 0], [1, 0, 0], [0, 1, 1]])
        def check_basis(loaded, ext):
            self.assertIsInstance(loaded, HarmonicOscillatorProductBasis, ext)
            self.assertEqual(loaded.ndim, basis.ndim, ext)
        def check_space(loaded, ext):
            self.assertIsInstance(loaded, BasisStateSpace, ext)
            self.assertEqual(np.asarray(loaded.excitations).tolist(), np.asarray(space.excitations).tolist(), ext)
        self._roundtrip('basis', basis, check_basis)
        self._roundtrip('space', space, check_space)

    def test_CoordinateSetRoundTrip(self):
        coords = CoordinateSet(np.random.rand(4, 3), CartesianCoordinates3D)
        def check(loaded, ext):
            self.assertIsInstance(loaded, CoordinateSet, ext)
            self.assertTrue(np.allclose(loaded, coords), ext)
            self.assertEqual(type(loaded.system), type(coords.system), ext)
        self._roundtrip('coords', coords, check)

    def test_MoleculeWithDerivativesRoundTrip(self):
        import McUtils
        fchk = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(McUtils.__file__))),
                            'ci', 'tests', 'TestData', 'HOD_freq.fchk')
        if not os.path.exists(fchk):
            self.skipTest("McUtils test data not found")
        mol = Molecule.from_file(fchk, 'fchk')  # (the D atom is why `to_state` can't ask RDKit for bonds)
        freqs = mol.get_normal_modes().freqs
        def check(loaded, ext):
            self.assertTrue(np.allclose(loaded.coords, mol.coords), ext)
            for a, b in zip(loaded.potential_derivatives, mol.potential_derivatives):
                self.assertTrue(np.allclose(a, b), ext)
            self.assertTrue(np.allclose(loaded.get_normal_modes().freqs, freqs), ext)
        self._roundtrip('mol_fchk', mol, check)

    def test_SharedMoleculeReference(self):
        mol = Molecule(*WATER)
        def check(loaded, ext):
            self.assertIs(loaded['a'], loaded['b'], ext)
        self._roundtrip('shared', {'a': mol, 'b': mol}, check)

    def test_ProtocolOneFileStillLoads(self):
        fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TestData', 'protocol1_molecule.json')
        from McUtils.Scaffolding import JSONSerializer
        with self.assertWarns(DeprecationWarning):
            with Checkpointer.from_file(fn, serializer=JSONSerializer(allow_pickle=True)) as chk:
                loaded = chk['mol']
        self.assertEqual(tuple(loaded.atoms), tuple(WATER[0]))
        self.assertTrue(np.allclose(loaded.coords, WATER[1]))
        with open(fn) as f:
            import json
            spec = json.load(f)['mol']
        from McUtils.Scaffolding import UnsafeLoadError
        with self.assertRaises(UnsafeLoadError):
            PseudoPickler().deserialize(spec)
    #endregion

    #region Molecule property cache
    def test_MoleculeCachedEvalPersistsForSameGeometry(self):
        fn = self.path('props.hdf5')
        calls = []

        def generate():
            calls.append(1)
            return np.arange(4)

        for _ in range(2):
            mol = Molecule(*WATER, checkpoint_file=fn)
            val = mol.cached_eval('values', generate)
        self.assertEqual(len(calls), 1)
        self.assertEqual(np.asarray(val).tolist(), [0, 1, 2, 3])

    @unittest.expectedFailure  # phase 7a: the cache is keyed by name only, so a new geometry gets the old value
    def test_MoleculeCachedEvalChangedGeometry(self):
        fn = self.path('props.hdf5')
        Molecule(*WATER, checkpoint_file=fn).cached_eval('values', lambda: np.arange(4))
        moved = Molecule(WATER[0], WATER[1] * 1.1, checkpoint_file=fn)
        val = moved.cached_eval('values', lambda: np.zeros(4))
        self.assertEqual(np.asarray(val).tolist(), [0., 0., 0., 0.])

    def test_MoleculePropertyCacheDelete(self):
        cache = MoleculePropertyCache(Molecule(*WATER), {})
        cache['a'] = 1
        del cache['a']
        with self.assertRaises(KeyError):
            cache['a']
    #endregion


if __name__ == '__main__':
    unittest.main()
