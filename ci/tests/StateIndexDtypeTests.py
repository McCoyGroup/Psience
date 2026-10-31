"""Regressions for exact integer IDs at state-space construction boundaries."""
import unittest
import numpy as np
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace
from McUtils.Numputils import external_unique

class IndexDtypeTests(unittest.TestCase):
    def setUp(self):self.basis=HarmonicOscillatorProductBasis(72)

    def test_constructor_and_indexer_use_declared_dtype(self):
        for states,mode in [([0,1],'indices'),([], 'indices'),
                            (np.zeros((1,72),dtype=int),'excitations'),
                            (np.empty((0,72),dtype=int),'excitations')]:
            space=BasisStateSpace(self.basis,states,mode=mode)
            self.assertEqual(space.indices.dtype,space.indices_dtype)

    def test_concatenate_preserves_keys_above_float_precision(self):
        first=2**53+1
        a=BasisStateSpace(self.basis,[first],mode='indices')
        b=BasisStateSpace(self.basis,[first+1],mode='indices')
        b.indices=[first+1]
        c=a.concatenate(b,track_excitations=False)
        self.assertEqual([int(x) for x in c.indices],[first,first+1])
        self.assertEqual(c.indices.dtype,c.indices_dtype)
        np.testing.assert_array_equal(external_unique(c.indices,sort=False),[first,first+1])

    def test_union_preserves_keys_above_float_precision(self):
        first=2**53+1
        a=BasisStateSpace(self.basis,[first],mode='indices')
        b=BasisStateSpace(self.basis,[first+1],mode='indices')
        b.indices=[first+1]
        c=a.union(b,track_excitations=False)
        self.assertEqual([int(x) for x in c.indices],[first,first+1])
        self.assertEqual(c.indices.dtype,c.indices_dtype)

if __name__=='__main__':unittest.main()
