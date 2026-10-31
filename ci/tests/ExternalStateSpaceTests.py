from pathlib import Path
import tempfile,unittest
import numpy as np
import scipy.sparse as sp
from unittest.mock import patch
from McUtils.Numputils import external_unique,save_sparse_matrix,load_sparse_matrix
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace,BraKetSpace,SelectionRuleStateSpace

class ExternalTests(unittest.TestCase):
    def test_stable_and_sorted_unique_rows(self):
        rng=np.random.default_rng(721)
        for width in [1,2,4]:
            x=rng.integers(-10,11,size=(117,width))
            if width==1:x=x[:,0]
            unique,first=np.unique(x,axis=0,return_index=True)
            for sorting in [True,False]:
                for chunk in [1,7,1000]:
                    out,where=external_unique([x[:13],x[13:60],x[60:]],sort=sorting,
                        return_index=True,chunk_size=chunk)
                    order=np.arange(len(unique)) if sorting else np.argsort(first)
                    np.testing.assert_array_equal(out,unique[order])
                    np.testing.assert_array_equal(where,first[order])
                    if chunk==7:
                        values,where=external_unique(x,sort=sorting,return_index=True,
                            return_values=False,chunk_size=chunk)
                        self.assertIsNone(values)
                        np.testing.assert_array_equal(where,first[order])

    def test_empty_and_file_output(self):
        np.testing.assert_array_equal(external_unique([np.empty(0,dtype=int)]),[])
        with tempfile.TemporaryDirectory() as temp:
            x=external_unique(np.array([7,3,7,2,1]),sort=False,output_file=Path(temp)/'ids.npy',chunk_size=2)
            self.assertIsInstance(x,np.memmap)
            np.testing.assert_array_equal(x,[7,3,2,1])

    def test_pair_order_and_orientation(self):
        rng=np.random.default_rng(831);basis=HarmonicOscillatorProductBasis(3)
        for size in [0,1,50,250]:
            x=rng.integers(0,15,size=(2,size))
            pairs=BraKetSpace(BasisStateSpace(basis,x[0],mode='indices'),
                              BasisStateSpace(basis,x[1],mode='indices'))
            memory=pairs.remove_duplicates()
            external=pairs.remove_duplicates(method='external',chunk_size=13)
            np.testing.assert_array_equal(memory.bras.indices,external.bras.indices)
            np.testing.assert_array_equal(memory.kets.indices,external.kets.indices)

    def test_streamed_nested_indices(self):
        basis=HarmonicOscillatorProductBasis(2)
        base=BasisStateSpace(basis,[0,1],mode='indices')
        blocks=[BasisStateSpace(basis,v,mode='indices') for v in [[2,3],[1,4]]]
        space=SelectionRuleStateSpace(base,blocks)
        np.testing.assert_array_equal(np.concatenate(list(space.iter_index_blocks())),space.indices)

    def test_mapped_sparse_product(self):
        a=sp.diags([np.arange(1,100),np.full(98,.7)],[0,1],format='csr')
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'csr.npz';save_sparse_matrix(path,a);mapped=load_sparse_matrix(path)
            for name in ['data','indices','indptr']:
                self.assertTrue(np.shares_memory(getattr(mapped.data,name),mapped._mapped_archive[name]))
            x=np.arange(a.shape[1],dtype=float)
            np.testing.assert_array_equal(mapped.dot(x),a.dot(x))
            # Replacing a cache file must not mutate an already mapped matrix.
            save_sparse_matrix(path,2*a)
            np.testing.assert_array_equal(mapped.dot(x),a.dot(x))
            np.testing.assert_array_equal(load_sparse_matrix(path).dot(x),2*a.dot(x))
            mapped._mapped_archive.close()
            del mapped

    def test_temporary_maps_closed_before_cleanup(self):
        import McUtils.Numputils.ExternalMemory as module
        original=module._close_mapping;seen=[]
        def close(array):seen.append(array);original(array)
        with patch.object(module,'_close_mapping',new=close):
            out=external_unique(np.array([7,3,2,7,4,1]),sort=False,chunk_size=2)
        np.testing.assert_array_equal(out,[7,3,2,4,1])
        self.assertTrue(seen)
        self.assertTrue(all(a._mmap.closed for a in seen))

    def test_numpy_only_merge_fallback(self):
        import McUtils.Numputils.ExternalMemory as module
        with patch.object(module,'_merge_function',return_value=module._merge):
            actual=external_unique(np.array([9,1,3,9,2,1]),sort=False,chunk_size=2)
        np.testing.assert_array_equal(actual,[9,1,3,2])

if __name__=='__main__':unittest.main()
