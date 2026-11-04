"""Frontier lifetime and exact stable-union regressions."""
import unittest,gc,weakref,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace
from Psience.BasisReps.ExternalStateSpaces import StateSpaceDiskStore

class PropertyFrontierStorageTests(unittest.TestCase):
    def test_unreferenced_frontier_mappings_are_released(self):
        store=StateSpaceDiskStore(block_size=7);self.addCleanup(store.close)
        refs=[]
        for i in range(300):
            array=store.array('tiny',np.int64,(1,));array[0]=i
            refs.append(weakref.ref(array._mmap));del array
        gc.collect()
        self.assertLessEqual(sum(r() is not None for r in refs),1)
        self.assertLessEqual(len(store._mappings),1)

    def test_live_numpy_views_keep_their_mapping(self):
        store=StateSpaceDiskStore();self.addCleanup(store.close)
        array=store.array('values',np.int64,(12,));array[:]=np.arange(12)
        mapping=weakref.ref(array._mmap);view=np.asarray(array)[2:8]
        del array;gc.collect()
        self.assertIsNotNone(mapping())
        np.testing.assert_array_equal(view,np.arange(2,8))
        del view;gc.collect()
        self.assertIsNone(mapping())

    def test_explicit_close_closes_live_mappings_before_unlink(self):
        store=StateSpaceDiskStore();array=store.array('values',np.int64,(12,));mapping=array._mmap
        path=store.directory;store.close()
        self.assertTrue(mapping.closed);self.assertFalse(path.exists())

    def test_bulk_union_preserves_sequential_values_and_order(self):
        rng=np.random.default_rng(813);basis=HarmonicOscillatorProductBasis(3)
        for sorted_output in [False,True]:
            with tempfile.TemporaryDirectory() as directory:
                store=StateSpaceDiskStore(directory,block_size=7)
                try:
                    resident=[BasisStateSpace(basis,rng.integers(0,80,size=19),mode='indices') for _ in range(9)]
                    for first_disk in [False,True]:
                        spaces=resident.copy()
                        if first_disk:spaces[0]=store.frontier(basis,[resident[0].indices],sort=False)
                        expected=store.coerce(spaces[0])
                        for space in spaces[1:]:expected=expected.union(space,sort=sorted_output)
                        actual=store.union_spaces(iter(spaces),sort=sorted_output)
                        np.testing.assert_array_equal(actual.indices,expected.indices)
                finally:store.close()

    def test_bulk_union_single_pass_and_empty_inputs(self):
        basis=HarmonicOscillatorProductBasis(3);store=StateSpaceDiskStore(block_size=3)
        self.addCleanup(store.close)
        spaces=[store.frontier(basis,[np.array([i,i+1])],sort=False) for i in range(8)]
        self.assertIsNone(store.union_spaces([]))
        self.assertIs(store.union_spaces([spaces[0]]),spaces[0])
        original=store.frontier
        with patch.object(store,'frontier',wraps=original) as build:
            actual=store.union_spaces(spaces)
            self.assertEqual(build.call_count,1)
        np.testing.assert_array_equal(actual.indices,np.arange(9))

if __name__=='__main__':unittest.main()
