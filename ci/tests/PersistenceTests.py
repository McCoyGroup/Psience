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

    @unittest.expectedFailure  # by design: without a `checkpoint_id` entries are keyed by name only; the identifier is yours to supply
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

    #region Compiled expression store on the shared entry store (phase 6)
    def test_CompiledExpressionStoreUsesEntryStore(self):
        import types
        from McUtils.Devutils import EntryStore
        from Psience.VPT2.CompiledExpressions import AnalyticExpressionStore
        store = AnalyticExpressionStore(self.path('compiled'))
        self.assertIsInstance(store.entries, EntryStore)
        wrapper = types.SimpleNamespace(expr=0, change=(), num_fixed=0, op=None)
        compiled = store.compile_factory(lambda: wrapper, context='zero', max_terms=10)
        self.assertTrue(compiled.manifest['complete'])
        signature = compiled.folder.name
        self.assertEqual(store.entries.keys(), [signature])  # folder/<signature>/manifest.json, no shards
        self.assertTrue(os.path.isfile(self.path(os.path.join('compiled', signature, 'manifest.json'))))
        # a second compile hits the stored entry without calling the factory
        def fail():
            raise AssertionError('should have been served from the store')
        again = store.compile_factory(fail, context='zero', max_terms=10)
        self.assertEqual(again.folder, compiled.folder)
        # aliases keep their old layout and point at the compiled entry
        aliases = os.listdir(self.path(os.path.join('compiled', 'aliases')))
        self.assertEqual(len(aliases), 1)
        import json
        with open(self.path(os.path.join('compiled', 'aliases', aliases[0]))) as f:
            self.assertEqual(json.load(f), {'compiled': signature})
        self.assertEqual([n for n in os.listdir(self.path('compiled')) if n.startswith('.')], [])

    def test_CompiledExpressionStoreFailedCompileLeavesNothing(self):
        from Psience.VPT2.CompiledExpressions import AnalyticExpressionStore
        store = AnalyticExpressionStore(self.path('compiled2'))
        def broken():
            raise RuntimeError('compile failed')
        with self.assertRaises(RuntimeError):
            store.compile_factory(broken, context='zero', max_terms=10)
        self.assertEqual(os.listdir(self.path('compiled2')), [])
    #endregion

    #region Persistent molecule property cache with a caller-supplied identifier (phase 7a)
    class _CountingEvaluator:
        """An energy evaluator (any object with `evaluate`) that counts its calls"""
        distance_units = 'BohrRadius'
        def __init__(self, k=1.0):
            self.k, self.calls = k, 0
        def evaluate(self, coords, order=0):
            self.calls += 1
            x = np.asarray(coords, dtype=float).reshape(-1)
            out = [np.array(self.k * np.sum(x ** 2) / 2)]
            if order >= 1:
                out.append(self.k * x)
            if order >= 2:
                out.append(self.k * np.eye(len(x)))
            return out[:order + 1]

    def test_MoleculeCachedEvalInObjectStore(self):
        from McUtils.Scaffolding import ObjectStore
        store = ObjectStore(self.path('mol_store'))
        calls = []
        def generate():
            calls.append(1)
            return {'a': np.arange(3.), 'n': 4}
        mol = Molecule(*WATER, checkpoint_file=store, checkpoint_id='water-1')
        first = mol.cached_eval('props', generate)
        again = Molecule(*WATER, checkpoint_file=store, checkpoint_id='water-1')  # same identifier: a hit
        self.assertEqual(again.cached_eval('props', generate)['n'], 4)
        self.assertEqual(len(calls), 1)
        Molecule(*WATER, checkpoint_file=store, checkpoint_id='water-2').cached_eval('props', generate)  # other identifier
        mol.cached_eval('props', generate, options={'order': 2})  # different options: a miss
        self.assertEqual(len(calls), 3)
        self.assertEqual(len(store.keys()), 3)
        meta = store.meta(store.find(provenance__property='props', provenance__identifier='water-1')[0])
        self.assertEqual(meta['provenance']['identifier'], 'water-1')
        self.assertEqual(len(store.find(provenance__identifier='water-2')), 1)
        self.assertTrue(np.array_equal(first['a'], np.arange(3.)))
        # the identifier is the caller's: a moved geometry under the same identifier is the same entry
        moved = Molecule(WATER[0], WATER[1] * 1.1, checkpoint_file=store, checkpoint_id='water-1')
        self.assertEqual(moved.cached_eval('props', lambda: 'recomputed')['n'], 4)
        def bad():
            return {1, 2}  # not serializable: returned, not cached
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            self.assertEqual(mol.cached_eval('odd', bad), {1, 2})
        self.assertTrue(any("couldn't cache" in str(w.message) for w in caught))
        self.assertEqual(mol.cached_eval('odd', lambda: 5, condition=lambda k: False), 5)  # bypassed

    def test_ObjectStoreCacheNeedsAnIdentifier(self):
        from McUtils.Scaffolding import ObjectStore
        store = ObjectStore(self.path('anon_store'))
        mol = Molecule(*WATER, checkpoint_file=store)
        with self.assertRaises(ValueError):
            mol.cached_eval('props', lambda: 1)
        self.assertEqual(store.keys(), [])
        # derived properties are only cached automatically when there is an identifier
        evaluator = self._CountingEvaluator()
        anon = Molecule(*WATER, energy_evaluator=evaluator, checkpoint_file=store)
        anon.calculate_energy(order=1); anon.calculate_energy(order=1)
        self.assertEqual(evaluator.calls, 2)
        self.assertEqual(store.keys(), [])

    def test_CheckpointerCacheScopedByIdentifier(self):
        fn = self.path('ident.json')
        calls = []
        def generate():
            calls.append(1)
            return [1, 2]
        for ident in ('a', 'a', 'b'):
            Molecule(*WATER, checkpoint_file=fn, checkpoint_id=ident).cached_eval('values', generate)
        self.assertEqual(len(calls), 2)
        Molecule(*WATER, checkpoint_file=fn).cached_eval('values', generate)  # no identifier: the plain key, as before
        self.assertEqual(len(calls), 3)

    def test_MoleculeCheckpointTrueIsRejected(self):
        with self.assertRaises(ValueError):
            Molecule(*WATER, checkpoint_file=True)  # no implicit store location
        with self.assertRaises(ValueError):
            Molecule(*WATER, energy_evaluator=self._CountingEvaluator(), checkpoint_file=True)

    def test_MoleculeEnergyCalculationIsCached(self):
        from McUtils.Scaffolding import ObjectStore
        store = ObjectStore(self.path('energy_store'))
        evaluator = self._CountingEvaluator()
        mol = Molecule(*WATER, energy_evaluator=evaluator, checkpoint_file=store, checkpoint_id='w')
        e1 = mol.calculate_energy(order=2)
        self.assertEqual(evaluator.calls, 1)
        e2 = mol.calculate_energy(order=2)
        again = Molecule(*WATER, energy_evaluator=evaluator, checkpoint_file=store, checkpoint_id='w').calculate_energy(order=2)
        self.assertEqual(evaluator.calls, 1)  # served from the store, even by a new molecule object
        for a, b in zip(e1, again):
            self.assertTrue(np.allclose(a, b))
        mol.calculate_energy(order=1)  # a different order is a different entry
        self.assertEqual(evaluator.calls, 2)
        Molecule(*WATER, energy_evaluator=evaluator, checkpoint_file=store, checkpoint_id='w2').calculate_energy(order=2)
        self.assertEqual(evaluator.calls, 3)  # another identifier misses
        # an evaluator *object* is not part of the key (the identifier covers it): the same identifier hits
        other = Molecule(*WATER, energy_evaluator=self._CountingEvaluator(k=2.0), checkpoint_file=store, checkpoint_id='w')
        other.calculate_energy(order=2)
        self.assertEqual(len(store.find(provenance__property='energy')), 3)
        self.assertEqual(other.energy_evaluator.calls, 0)
        # coordinates given explicitly aren't this molecule's geometry: never cached here
        before = evaluator.calls
        mol.calculate_energy(WATER[1], order=2)
        mol.calculate_energy(WATER[1], order=2)
        self.assertEqual(evaluator.calls, before + 2)
        # without a persistent checkpoint nothing is cached at all
        plain = Molecule(*WATER, energy_evaluator=evaluator)
        before = evaluator.calls
        plain.calculate_energy(order=2); plain.calculate_energy(order=2)
        self.assertEqual(evaluator.calls, before + 2)

    def test_MoleculeNormalModesAndPointGroupAreCached(self):
        from McUtils.Scaffolding import ObjectStore
        store = ObjectStore(self.path('modes_store'))
        evaluator = self._CountingEvaluator()
        mol = Molecule(*WATER, energy_evaluator=evaluator, checkpoint_file=store, checkpoint_id='w')
        mol.potential_derivatives = [np.zeros(9), np.eye(9) * 0.5]
        modes = mol.get_normal_modes()
        self.assertEqual(len(store.find(provenance__property='normal_modes')), 1)
        other = Molecule(*WATER, checkpoint_file=store, checkpoint_id='w')
        other.potential_derivatives = [np.zeros(9), np.eye(9) * 0.5]
        cached = other.get_normal_modes()
        self.assertEqual(len(store.find(provenance__property='normal_modes')), 1)
        self.assertTrue(np.allclose(np.asarray(modes.freqs), np.asarray(cached.freqs)))
        stiffer = Molecule(*WATER, checkpoint_file=store, checkpoint_id='w')  # the molecule's own derivatives are part of the key
        stiffer.potential_derivatives = [np.zeros(9), np.eye(9) * 2.0]
        stiffer.get_normal_modes()
        self.assertEqual(len(store.find(provenance__property='normal_modes')), 2)
        pg = mol.point_group
        self.assertEqual(Molecule(*WATER, checkpoint_file=store, checkpoint_id='w').point_group.name, pg.name)
        self.assertEqual(len(store.find(provenance__property='point_group')), 1)

    def _id_mol(self, store, ident='w', evaluator=None, **kw):
        evaluator = self._CountingEvaluator() if evaluator is None else evaluator
        return Molecule(*WATER, energy_evaluator=evaluator, checkpoint_file=store, checkpoint_id=ident, **kw), evaluator

    def test_ModifyInheritsIdentifierAndSuspendsWhatChanged(self):
        from McUtils.Scaffolding import ObjectStore
        store = ObjectStore(self.path('modify_store'))
        mol, evaluator = self._id_mol(store)
        mol.calculate_energy(order=1)
        pg = mol.point_group
        self.assertEqual(evaluator.calls, 1)
        # a copy, or a change that nothing cached depends on, shares the cache and the identifier
        for same in (mol.copy(), mol.modify(display_mode='plain'), mol.modify(meta={'note': 1})):
            self.assertEqual(same.checkpoint.identifier, 'w')
            same.calculate_energy(order=1)
            self.assertEqual(same.point_group.name, pg.name)
        self.assertEqual(evaluator.calls, 1)
        # new coordinates: nothing cached under the old identifier is reused or written, whatever its key
        moved = mol.modify(coords=WATER[1] * 1.1)
        self.assertEqual(moved.checkpoint.identifier, 'w')
        before = len(store.keys())
        moved.calculate_energy(order=1); moved.calculate_energy(order=1)
        self.assertEqual(evaluator.calls, 3)
        self.assertEqual(moved.cached_eval('anything', lambda: 5), 5)
        moved.point_group
        self.assertEqual(len(store.keys()), before)
        # ...until the caller says the changed molecule has its own entries
        named = mol.modify(coords=WATER[1] * 1.1, checkpoint_id='w-stretched')
        named.calculate_energy(order=1); named.calculate_energy(order=1)
        self.assertEqual(evaluator.calls, 4)
        self.assertEqual(len(store.find(provenance__identifier='w-stretched')), 1)
        # only the values that depend on the change are suspended: a new evaluator leaves the point group alone
        other = mol.modify(energy_evaluator=self._CountingEvaluator(k=3.0))
        other.energy_evaluator.calls = 0
        other.calculate_energy(order=1); other.calculate_energy(order=1)
        self.assertEqual(other.energy_evaluator.calls, 2)
        before = len(store.find(provenance__property='point_group'))
        self.assertEqual(other.point_group.name, pg.name)
        self.assertEqual(len(store.find(provenance__property='point_group')), before)
        # a change record is carried forward through further modify calls
        again = moved.modify(display_mode='plain')
        again.calculate_energy(order=1); again.calculate_energy(order=1)
        self.assertEqual(again.checkpoint.changed, {'coords'})
        # `None` drops caching, and a molecule without an identifier passes nothing on
        dropped = mol.modify(checkpoint_id=None)
        self.assertFalse(dropped.checkpoint.active)
        plain = Molecule(*WATER, checkpoint_file=self.path('legacy.json'))
        self.assertFalse(plain.modify().checkpoint.persistent)

    def test_SettersSuspendDependentValues(self):
        from McUtils.Scaffolding import ObjectStore
        store = ObjectStore(self.path('setter_store'))
        mol, evaluator = self._id_mol(store)
        mol.potential_derivatives = [np.zeros(9), np.eye(9) * 0.5]
        mol.get_normal_modes(); mol.calculate_energy(order=1); mol.calculate_energy(order=1)
        self.assertEqual(evaluator.calls, 1)  # computing the properties doesn't suspend them
        self.assertEqual(len(store.find(provenance__property='normal_modes')), 1)
        # the derivatives are part of the normal modes' key, so changing them just means a new entry
        mol.potential_derivatives = [np.zeros(9), np.eye(9) * 2.0]
        mol.get_normal_modes()
        self.assertEqual(len(store.find(provenance__property='normal_modes')), 2)
        mol.calculate_energy(order=1)
        self.assertEqual(evaluator.calls, 1)  # the energy doesn't depend on the derivatives
        # new coordinates are not part of any key: everything depending on them is suspended
        mol.coords = WATER[1] * 1.1
        mol.get_normal_modes()
        self.assertEqual(len(store.find(provenance__property='normal_modes')), 2)
        mol.calculate_energy(order=1); mol.calculate_energy(order=1)
        self.assertEqual(evaluator.calls, 3)
    #endregion

    #region 7c: AIMD trajectory files
    @staticmethod
    def _aimd(path=None, steps_kw=None, **kw):
        from Psience.AIMD import AIMDSimulator
        masses = [1836., 1836., 1836.]
        x0 = np.array([[0., 0., 0.], [1.8, 0., 0.], [0., 2.0, 0.2]])
        v0 = np.array([[0., 0., 0.], [0., .002, 0.], [.001, 0., 0.]])
        force = lambda x: -0.3 * (x - x0)  # harmonic about the start
        opts = dict(timestep=5., velocities=v0, track_kinetic_energy=True, track_velocities=True)
        opts.update(kw)
        return AIMDSimulator(masses, x0 + .1, force, trajectory_file=path, **opts), force

    def test_AIMDStreamsFramesToFile(self):
        ref, _ = self._aimd()
        ref.propagate(25)
        path = self.path('aimd.nptar')
        sim, _ = self._aimd(path, flush_every=10)
        sim.propagate(25)
        self.assertEqual(len(sim.trajectory), 6)  # 26 frames, 20 flushed
        self.assertEqual(sim.num_frames, 26)
        sim.flush()
        self.assertEqual(len(sim.trajectory), 0)
        for a, b in zip(ref.extract_trajectory(embed=False, extract_velocities=True),
                        sim.extract_trajectory(embed=False, extract_velocities=True)):
            self.assertTrue(np.allclose(a, b))
        with self.assertRaises(ValueError):
            self._aimd(path)  # won't clobber a run

    def test_AIMDRestartReproducesRemainingSteps(self):
        from Psience.AIMD import AIMDSimulator
        for kw in (dict(), dict(step_predictor='beeman'),
                   dict(thermostat=dict(type='langevin', seed=7, friction=.05), target_temperature=300.)):
            with self.subTest(**{k: str(v) for k, v in kw.items()}):
                tag = str(abs(hash(str(kw))))
                ref, _ = self._aimd(**kw)
                ref.propagate(40)
                path = self.path('restart_{}.nptar'.format(tag))
                sim, force = self._aimd(path, flush_every=7, **kw)
                sim.propagate(23)
                sim.flush()
                del sim  # "crash"
                restarted = AIMDSimulator.restart(path, force, **{k: v for k, v in kw.items()
                                                                  if k in ('thermostat', 'step_predictor')})
                self.assertEqual(restarted.steps, 23)
                restarted.propagate(17)
                restarted.flush()
                a = ref.extract_trajectory(embed=False, extract_velocities=True)
                b = restarted.extract_trajectory(embed=False, extract_velocities=True)
                self.assertEqual(a[0].shape, b[0].shape)
                for x, y in zip(a, b):
                    self.assertTrue(np.allclose(x, y, atol=1e-10))

    def test_AIMDRestartChecksWhatItCannotStore(self):
        from Psience.AIMD import AIMDSimulator
        path = self.path('check.nptar')
        sim, force = self._aimd(path, thermostat=dict(type='langevin', seed=1), target_temperature=300.)
        sim.propagate(3)
        sim.flush()
        with self.assertRaises(ValueError):
            AIMDSimulator.restart(path, force)  # thermostat missing
        with self.assertRaises(ValueError):
            AIMDSimulator.restart(self.path('nothing.nptar'), force)
    #endregion

    #region 7d: persistent VPT2 polynomial cache
    @staticmethod
    def _poly_values(n, seed=0):
        from Psience.VPT2.Analytic import _HashableNumPyIndex
        rng = np.random.default_rng(seed)
        return {_HashableNumPyIndex(rng.integers(0, 5, size=3)): float(rng.normal()) for _ in range(n)}

    def test_PersistentPolyCacheRoundTrip(self):
        from Psience.VPT2.PersistentCaches import PersistentEvaluationCache
        root = self.path('vpt2_store')
        poly_key = (('poly', 1, 2), (1, -1), None)
        cache = PersistentEvaluationCache(root, signature='impl-1')
        self.assertNotIn(poly_key, cache)
        cache[poly_key] = values = self._poly_values(20)
        self.assertEqual(cache.flush(), 1)
        self.assertEqual(cache.flush(), 0)  # nothing new
        fresh = PersistentEvaluationCache(root, signature='impl-1')  # a new process
        self.assertIn(poly_key, fresh)
        self.assertEqual(fresh[poly_key], values)
        self.assertEqual(fresh.disk_hits, 1)
        values.update(self._poly_values(5, seed=1))  # the dict grows as new states are evaluated
        self.assertEqual(cache.flush(), 1)
        self.assertEqual(PersistentEvaluationCache(root, signature='impl-1')[poly_key], values)

    def test_PersistentPolyCacheKeysAndFallbacks(self):
        from Psience.VPT2.PersistentCaches import PersistentEvaluationCache
        root = self.path('vpt2_store2')
        a = PersistentEvaluationCache(root, signature='impl-1')
        a[('p', 1)] = self._poly_values(4)
        a.flush()
        self.assertNotIn(('p', 2), PersistentEvaluationCache(root, signature='impl-1'))  # other polynomial
        self.assertNotIn(('p', 1), PersistentEvaluationCache(root, signature='impl-2'))  # other implementation
        self.assertIn(('p', 1), PersistentEvaluationCache(root, signature='impl-1'))

        def unfingerprintable(key):
            raise TypeError('no identity')
        b = PersistentEvaluationCache(root, fingerprint=unfingerprintable, signature='x')
        b[('q',)] = self._poly_values(3)
        b['plain'] = {'not': 'states'}
        self.assertEqual(b.flush(), 0)  # memory only, no error
        self.assertEqual(len(b), 2)

    def test_EvaluatorPersistentCacheHooks(self):
        from Psience.VPT2.Analytic import PerturbationTheoryExpressionEvaluator as Ev
        from Psience.VPT2.PersistentCaches import PersistentEvaluationCache
        original = Ev._poly_cache
        try:
            with self.assertRaises(ValueError):
                Ev.use_persistent_cache(None)
            with self.assertRaises(ValueError):
                Ev.use_persistent_cache(True)
            self.assertIs(Ev._poly_cache, original)
            layer = Ev.use_persistent_cache(self.path('hook_store'))
            self.assertIsInstance(Ev._poly_cache, PersistentEvaluationCache)
            self.assertIs(layer.memory, original)
            Ev._poly_cache[('k', 0)] = self._poly_values(3)
            self.assertEqual(Ev.flush_persistent_cache(), 1)
            Ev.use_persistent_cache(False)
            self.assertIs(Ev._poly_cache, original)
            self.assertEqual(Ev.flush_persistent_cache(), 0)
        finally:
            Ev._poly_cache = original
    #endregion

    #region 7e: scan results in an object store
    class _ScanLoader:
        def __call__(self, path):
            return np.loadtxt(path)

    class _ScanExtractor:
        __fingerprint_exclude__ = ('calls',)

        def __init__(self, scale=1.0):
            self.scale = scale
            self.calls = 0

        def __call__(self, data):
            self.calls += 1
            return {'energy': np.array(self.scale * data.sum()), 'vec': self.scale * data}

    def _make_scan(self):
        import json
        from Psience.Data.ScanManager import ScanManager
        mgr = ScanManager(self.path('scan_out'), scan_id='t')
        os.makedirs(mgr.scan_dir)
        steps = []
        for i in range(3):
            for j in range(2):
                name = 'scan_{:03d}_{:03d}'.format(i, j)
                with open(os.path.join(mgr.scan_dir, name + '.out'), 'w') as f:
                    f.write('{} {}\n'.format(i + 1, j + 1))
                steps.append({'index': [i, j], 'values': [float(i), float(j)], 'file': name + '.in'})
        json.dump({'scan_id': 't', 'shape': [3, 2], 'steps': steps}, open(mgr.scan_info_file, 'w'))
        return mgr

    def test_ScanManagerStoresPointsAndSkipsReruns(self):
        from McUtils.Scaffolding import ObjectStore
        mgr = self._make_scan()
        gen = lambda f: os.path.splitext(f)[0] + '.out'
        store = ObjectStore(self.path('scan_store'))
        extractor = self._ScanExtractor()
        first = mgr.parse(extractor, output_file_generator=gen, molecule_loader=self._ScanLoader(), store=store)
        self.assertEqual(extractor.calls, 6)
        self.assertEqual(len(store.keys()), 6)
        self.assertEqual(sorted(tuple(store.meta(k)['provenance']['grid_index']) for k in store.keys()),
                         [(i, j) for i in range(3) for j in range(2)])
        self.assertEqual(first['energy'].shape, (3, 2))
        self.assertEqual(first['energy'][2, 1], 5.0)
        again = self._ScanExtractor()
        second = mgr.parse(again, output_file_generator=gen, molecule_loader=self._ScanLoader(), store=store)
        self.assertEqual(again.calls, 0)  # every point from the store
        for k in first:
            self.assertTrue(np.array_equal(first[k], second[k]))
        # a changed output re-reads that point only; a changed extractor re-reads all
        with open(os.path.join(mgr.scan_dir, 'scan_001_001.out'), 'w') as f:
            f.write('10 10\n')
        changed = self._ScanExtractor()
        third = mgr.parse(changed, output_file_generator=gen, molecule_loader=self._ScanLoader(), store=store)
        self.assertEqual(changed.calls, 1)
        self.assertEqual(third['energy'][1, 1], 20.0)
        scaled = self._ScanExtractor(scale=2.0)
        mgr.parse(scaled, output_file_generator=gen, molecule_loader=self._ScanLoader(), store=store)
        self.assertEqual(scaled.calls, 6)

    def test_ScanManagerStoreMissingOutputsAndUnfingerprintable(self):
        mgr = self._make_scan()
        gen = lambda f: os.path.splitext(f)[0] + '.out'
        os.remove(os.path.join(mgr.scan_dir, 'scan_000_000.out'))
        extractor = self._ScanExtractor()
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always')
            found = mgr.parse(extractor, output_file_generator=gen, molecule_loader=self._ScanLoader(),
                              store=self.path('scan_store2'))
        self.assertTrue(np.isnan(found['energy'][0, 0]))
        self.assertEqual(extractor.calls, 5)
        extractor.opaque = object()  # can't be identified: parse still works, just without the store
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            mgr.parse(extractor, output_file_generator=gen, molecule_loader=self._ScanLoader(),
                      store=self.path('scan_store3'))
        self.assertTrue(any("aren't cached" in str(w.message) for w in caught))
        self.assertEqual(extractor.calls, 10)

    def test_ScanManagerManifestIsWrittenAtomically(self):
        import json
        mgr = self._make_scan()
        before = os.listdir(mgr.scan_dir)
        info = mgr.load_scan_info()
        info['note'] = 'x'
        from McUtils import Devutils as dev
        dev.atomic_write_json(mgr.scan_info_file, info)
        self.assertEqual(mgr.load_scan_info()['note'], 'x')
        self.assertEqual(sorted(os.listdir(mgr.scan_dir)), sorted(before))  # no temp files left behind
    #endregion


if __name__ == '__main__':
    unittest.main()
