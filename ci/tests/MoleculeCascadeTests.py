"""
Tests for the property cascade on `Molecule`, `MolecularEmbedding` and `MolecularTopology`
(`McUtils.Devutils.CascadingPropertyManager`): assigning an input (coordinates, evaluators, derivatives, normal modes, ...)
resets the values derived from it, and tells the molecule's persistent cache to stop vouching for the
entries that depend on it.
"""

import os, shutil, tempfile, unittest
import numpy as np

from McUtils.Scaffolding import ObjectStore
from McUtils.Devutils import CascadingPropertyManager
from Psience.Molecools import Molecule

__all__ = ["MoleculeCascadeTests"]

WATER = (['O', 'H', 'H'], np.array([[0., 0., 0.], [0.96, 0., 0.], [-0.24, 0.93, 0.]]))


class _Counting:
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


class MoleculeCascadeTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.mkdtemp(prefix='cascade_')

    def tearDown(self):
        shutil.rmtree(self._tmp, ignore_errors=True)

    def path(self, *p):
        return os.path.join(self._tmp, *p)

    def test_EvaluatorsAreCascadingProperties(self):
        for name in ('energy_evaluator', 'dipole_evaluator', 'charge_evaluator', 'polarizability_evaluator'):
            self.assertIn(name, CascadingPropertyManager.graph_for(Molecule))
        a, b = _Counting(1.0), _Counting(3.0)
        mol = Molecule(*WATER, energy_evaluator=a)
        e0 = mol.energy
        self.assertEqual((mol.energy, a.calls), (e0, 1))
        mol.energy_evaluator = a  # the same object: nothing to reset
        self.assertEqual((mol.energy, a.calls), (e0, 1))
        mol.energy_evaluator = b  # reassigned in place: the energy belonged to the old evaluator
        self.assertAlmostEqual(mol.energy, 3 * e0)
        self.assertEqual(b.calls, 1)

    def test_EvaluatorAssignmentSuspendsPersistentEntries(self):
        store = ObjectStore(self.path('store'))
        a, b = _Counting(1.0), _Counting(3.0)
        mol = Molecule(*WATER, energy_evaluator=a, checkpoint_file=store, checkpoint_id='w')
        mol.calculate_energy(order=1); mol.calculate_energy(order=1)
        pg = mol.point_group
        self.assertEqual(a.calls, 1)
        n = len(store.keys())
        mol.energy_evaluator = b
        mol.calculate_energy(order=1); mol.calculate_energy(order=1)
        self.assertEqual(b.calls, 2)  # not cached under an identifier that was given for the old evaluator
        self.assertEqual(len(store.keys()), n)
        self.assertEqual(mol.point_group.name, pg.name)  # the point group doesn't depend on the evaluator
        self.assertEqual(len(store.keys()), n)

    def test_GeometryChangeRefreshesEnergyAndPointGroup(self):
        mol = Molecule(*WATER, energy_evaluator=_Counting())
        e0, pg = mol.energy, mol.point_group
        self.assertIs(mol.point_group, pg)
        mol.coords = WATER[1] * 1.1
        self.assertNotAlmostEqual(mol.energy, e0)  # was returned stale before
        self.assertIsNot(mol.point_group, pg)

    def test_DerivedObjectsFollowTheirInputs(self):
        mol = Molecule(*WATER)
        everything = ('_hamiltonian', '_evaluator', '_mode_embedding')
        for change, reset in (
                (lambda: setattr(mol, 'potential_derivatives', [np.zeros(9), np.eye(9)]), everything),
                (lambda: setattr(mol.potential_surface, 'derivatives', [np.zeros(9), 2 * np.eye(9)]), everything),
                (lambda: setattr(mol.dipole_surface, 'derivatives', [np.zeros(3), np.zeros((9, 3))]), ('_hamiltonian',)),
                (lambda: setattr(mol, 'coords', WATER[1] * 1.05), everything),
                (lambda: setattr(mol, 'embedding', mol.embedding), everything),
        ):
            mol.hamiltonian, mol.evaluator, mol.mode_embedding = 'h', 'e', 'm'
            change()
            for slot in everything:
                self.assertEqual(getattr(mol, slot) is None, slot in reset, slot)
        # what doesn't depend on a change is left alone
        mol.hamiltonian, mol.evaluator, mol.mode_embedding = 'h', 'e', 'm'
        mol.dipole_evaluator = _Counting()
        self.assertEqual((mol._hamiltonian, mol._evaluator, mol._mode_embedding), ('h', 'e', 'm'))

    def test_ComputedFillsAreNotChanges(self):
        store = ObjectStore(self.path('store'))
        mol = Molecule(*WATER, energy_evaluator=_Counting(), checkpoint_file=store, checkpoint_id='w')
        calls = []
        mol.cached_eval('x', lambda: calls.append(1) or 5)
        mol.hamiltonian = 'h'
        mol.get_cartesian_potential_derivatives(order=2)  # stores what it computed on the molecule
        self.assertIsNotNone(mol.potential_derivatives)
        self.assertEqual(mol.checkpoint.changed, set())  # computed from the current state, so nothing was invalidated
        self.assertIsNone(mol._hamiltonian)  # ...but what was built without them is rebuilt
        mol.cached_eval('x', lambda: calls.append(1) or 5)
        self.assertEqual(len(calls), 1)
        mol.potential_derivatives = [np.zeros(9), np.eye(9)]  # while this one is a change
        self.assertEqual(mol.checkpoint.changed, {'potential_derivatives'})

    def test_EmbeddingCachesFollowTheCoordinates(self):
        mol = Molecule(*WATER, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        emb = mol.embedding
        tr, frame = emb.translation_rotation_modes, emb.inertial_frame
        self.assertIs(emb.translation_rotation_modes, tr)
        self.assertIs(emb.inertial_frame, frame)  # the frame was recomputed on every access before
        ints = emb.internal_coordinates
        jacobians = emb._jacobians
        mol.coords = WATER[1] * 1.1
        self.assertIsNot(emb.translation_rotation_modes, tr)  # was stale before
        self.assertFalse(np.allclose(emb.inertial_frame[0], frame[0]))
        self.assertFalse(np.allclose(np.asarray(emb.internal_coordinates), np.asarray(ints)))
        self.assertIsNot(emb._jacobians, jacobians)
        # a new specification (or Z-matrix) resets what was derived from the old one
        ints, jacobians = emb.internal_coordinates, emb._jacobians
        emb.zmatrix = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 1, 0, -1]]
        self.assertIsNot(emb._jacobians, jacobians)
        self.assertIsNot(emb.internal_coordinates, ints)

    def test_TopologyGraphsFollowTheBonds(self):
        mol = Molecule(*WATER, bonds=[[0, 1], [0, 2]])
        topo = mol.topology
        graph, frags = topo.edge_graph, topo.fragment_indices
        self.assertIs(topo.edge_graph, graph)
        self.assertEqual(len(frags), 1)
        topo.bonds = [[0, 1]]
        self.assertIsNot(topo.edge_graph, graph)
        self.assertEqual(len(topo.fragment_indices), 2)


if __name__ == '__main__':
    unittest.main()
