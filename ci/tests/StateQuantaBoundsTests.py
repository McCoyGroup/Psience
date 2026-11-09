import itertools
import unittest
import numpy as np
from Psience.BasisReps import BasisStateSpace


class StateQuantaBoundsTests(unittest.TestCase):
    def test_inclusive_quanta_and_energy_bounds(self):
        frequencies = np.array([.29, .1, .17])
        for minimum, maximum in ((None, 0), (None, 6), (6, 6), (3, 4)):
            actual = BasisStateSpace.states_under_freq_threshold(frequencies, 1., min_quanta=minimum, max_quanta=maximum)
            expected = [state for state in itertools.product(range(maximum + 1), repeat=3)
                        if (minimum is None or sum(state) >= minimum) and sum(state) <= maximum
                        and np.dot(state, frequencies) <= 1. + 1e-14]
            self.assertEqual(set(map(tuple, actual)), set(expected))

    def test_initial_state_respects_both_bounds(self):
        cases = [dict(initial_state=[1, 1], min_quanta=3, max_quanta=3),
                 dict(initial_state=[2, 1], min_quanta=2, max_quanta=2),
                 dict(initial_state=[2, 1], min_quanta=3, max_quanta=3)]
        for options in cases:
            actual = BasisStateSpace.states_under_freq_threshold(np.array([.1, .2]), 1., **options)
            for state in actual:
                self.assertGreaterEqual(sum(state), options['min_quanta'])
                self.assertLessEqual(sum(state), options['max_quanta'])

    def test_initial_state_above_energy_window_is_excluded(self):
        actual = BasisStateSpace.states_in_windows(np.array([1., 2.]), [[0., 1.]], initial_state=[0, 1], max_quantas=[2])[0]
        self.assertEqual(len(actual), 0)

    def test_window_iterator_matches_independent_occupation_enumeration(self):
        frequencies = np.array([.25, .5, .75])
        windows = [[0.,1.5],[.5,2.],[0.,1.5],[3.,2.]]
        initial = np.array([1,0,0])
        maximum = np.array([3,2,2])
        expected = [[] for _ in windows]
        rows = [np.array(s) for s in itertools.product(range(1,4),range(3),range(3))]
        rows.sort(key=lambda s:(sum(s),tuple(i for i,n in enumerate(s-initial) for _ in range(n))))
        for row in rows:
            if row[1] != initial[1]:
                continue
            energy = float(np.dot(initial,frequencies))
            for mode,n in enumerate(row-initial):
                for _ in range(n):
                    energy += frequencies[mode]
            for i,(low,high) in enumerate(windows):
                if 1 <= sum(row) <= 4 and low <= energy <= high:
                    expected[i].append(row)
        expected = [np.asarray(a,dtype=int).reshape((-1,3)) for a in expected]
        options = dict(initial_state=initial,max_state=maximum,fixed_modes=[1],
                       min_quantas=[1]*4,max_quantas=[4]*4,block_size=2,max_bytes=512)
        for sparse in [False,True]:
            pieces = [[] for _ in windows]
            for tags,block in BasisStateSpace.states_in_windows(
                    frequencies,windows,sparse=sparse,**options, return_iterator=True):
                values = block.asarray() if sparse else block
                for tag,row in zip(tags,values):
                    pieces[int(tag)].append(row.copy())
            for actual,reference in zip(pieces,expected):
                np.testing.assert_array_equal(np.asarray(actual,dtype=int).reshape((-1,3)),reference)
            actual = BasisStateSpace.states_in_windows(frequencies,windows,sparse=sparse,**options)
            for value,reference in zip(actual,expected):
                np.testing.assert_array_equal(value.asarray() if sparse else value,reference)
        np.testing.assert_array_equal(initial,[1,0,0])

    def test_window_generation_filter_matches_postfilter_and_preserves_readonly_input(self):
        frequencies = np.array([.25,.5,.75])
        windows = [[0.,2.],[0.,2.]]
        reference = BasisStateSpace.states_in_windows(frequencies,windows,max_quantas=3)
        calls = []
        def predicate(state,window):
            self.assertFalse(state.flags.writeable)
            calls.append(window)
            return state[1] % 2 == window
        actual = BasisStateSpace.states_in_windows(frequencies,windows,
            max_quantas=3,block_size=2,max_bytes=256,filter=predicate)
        for i,(old,new) in enumerate(zip(reference,actual)):
            np.testing.assert_array_equal(new,old[old[:,1]%2==i])
        self.assertTrue(calls)
        from unittest import mock
        with mock.patch('Psience.BasisReps.StateSpaces.SparseArray.from_data',
                        side_effect=AssertionError('rejected state allocated an output block')):
            self.assertEqual(list(BasisStateSpace.states_in_windows(
                frequencies,windows,max_quantas=3,sparse=True,filter=lambda state,window:False, return_iterator=True)),[])

    def test_iterator_is_lazy_and_blocks_own_their_storage(self):
        calls = []
        def predicate(state,window):
            calls.append(int(sum(state)))
            return True
        iterator = BasisStateSpace.states_in_windows(
            np.ones(72),[[0,6]],max_quantas=6,block_size=3,filter=predicate, return_iterator=True)
        tags,first = next(iterator)
        np.testing.assert_array_equal(first,np.zeros((1,72),dtype=int))
        self.assertEqual(calls,[0])
        saved = first.copy()
        tags,second = next(iterator)
        self.assertEqual(second.shape,(3,72))
        self.assertEqual(calls,[0,1,1,1])
        np.testing.assert_array_equal(first,saved)
        iterator.close()

    def test_dense_and_sparse_block_payloads_obey_byte_target(self):
        from McUtils.Numputils import SparseArray
        for sparse in [False,True]:
            blocks = list(BasisStateSpace.states_in_windows(
                [.25,.5,.75],[[0,2],[0,2]],max_quantas=4,
                block_size=100,max_bytes=512,sparse=sparse, return_iterator=True))
            self.assertTrue(len(blocks)>1)
            for tags,block in blocks:
                self.assertLessEqual(len(tags),100)
                self.assertEqual(len(tags),block.shape[0])
                if sparse:
                    self.assertIsInstance(block,SparseArray)
                    csr=block.data
                    payload=csr.data.nbytes+csr.indices.nbytes+csr.indptr.nbytes
                else:
                    payload=block.nbytes
                self.assertLessEqual(tags.nbytes+payload,512)

    def test_mixed_bounded_and_energy_bounded_windows(self):
        values=BasisStateSpace.states_in_windows([1.,2.],[[0.,4.],[0.,4.]],
                                                max_quantas=[0,None],block_size=2)
        np.testing.assert_array_equal(values[0],[[0,0]])
        expected={(a,b) for a in range(5) for b in range(3) if a+2*b<=4}
        self.assertEqual(set(map(tuple,values[1])),expected)

    def test_zero_frequency_with_finite_quanta_and_large_initial_occupations(self):
        for sparse in [False,True]:
            actual=BasisStateSpace.states_in_windows([0.,1.],[[0.,2.]],
                max_quantas=2,sparse=sparse)[0]
            values=actual.asarray() if sparse else actual
            expected={(a,b) for a in range(3) for b in range(3-a)}
            self.assertEqual(set(map(tuple,values)),expected)
            actual=BasisStateSpace.states_in_windows([1.,2.],[[130.,134.]],
                initial_state=[130,0],max_state=[132,1],max_quantas=134,sparse=sparse)[0]
            values=actual.asarray() if sparse else actual
            self.assertGreaterEqual(values[:,0].min(),130)
            self.assertIn((132,1),set(map(tuple,values)))

    def test_window_iterator_validates_memory_and_finiteness(self):
        invalid=[dict(block_size=0),dict(max_bytes=0),dict(max_bytes=1),
                 dict(initial_state=[-1,0]),dict(initial_state=[.5,0])]
        for options in invalid:
            with self.assertRaises(ValueError):
                list(BasisStateSpace.states_in_windows([1.,2.],[[0.,2.]],max_quantas=2,**options, return_iterator=True))
        with self.assertRaises(ValueError):
            list(BasisStateSpace.states_in_windows([0.,1.],[[0.,2.]], return_iterator=True))
        with self.assertRaises(ValueError):
            list(BasisStateSpace.states_in_windows([float('nan')],[[0.,2.]], return_iterator=True))
        self.assertEqual(BasisStateSpace.states_in_windows([1.],[]),[])

    def test_sparse_window_dtype_promotes_without_losing_large_occupations(self):
        for occupation,itemsize in [(0,1),(130,2),(40000,4)]:
            value=BasisStateSpace.states_in_windows([1.],[[occupation,occupation]],
                initial_state=[occupation],max_quantas=occupation,sparse=True)[0]
            np.testing.assert_array_equal(value.asarray(),[[occupation]])
            self.assertEqual(value.data.data.dtype.itemsize,itemsize)

    def test_storage_traversal_and_iterator_are_independent(self):
        from McUtils.Numputils import SparseArray
        from unittest import mock
        frequencies=[.25,.5,.75]
        windows=[[0.,2.],[.5,1.5],[0.,2.]]
        options=dict(initial_state=[1,0,0],max_state=[3,2,2],
                     min_quantas=[1,2,0],max_quantas=[4,3,2],
                     filter=lambda state,window:(state[1]+window)%2==0)
        reference=BasisStateSpace.states_in_windows(frequencies,windows,**options)
        for sparse,sparse_traversal,return_iterator in itertools.product([False,True],repeat=3):
            with self.subTest(sparse=sparse,traversal=sparse_traversal,iterator=return_iterator):
                output=BasisStateSpace.states_in_windows(frequencies,windows,**options,
                    sparse=sparse,sparse_traversal=sparse_traversal,
                    return_iterator=return_iterator,block_size=2,max_bytes=512)
                if return_iterator:
                    rows=[[] for _ in windows]
                    for tags,block in output:
                        self.assertEqual(isinstance(block,SparseArray),sparse)
                        dense=block.asarray() if sparse else block
                        for tag,row in zip(tags,dense):
                            rows[tag].append(row.copy())
                    output=[np.asarray(row,dtype=int).reshape((-1,3)) for row in rows]
                else:
                    self.assertTrue(all(isinstance(block,SparseArray)==sparse for block in output))
                    if sparse: output=[block.asarray() for block in output]
                for actual,expected in zip(output,reference):
                    np.testing.assert_array_equal(actual,expected)
        # Dense output must never pass through sparse packing/densification.
        with mock.patch('Psience.BasisReps.StateSpaces.SparseArray.from_data',
                        side_effect=AssertionError('dense result used sparse storage')):
            eager=BasisStateSpace.states_in_windows(frequencies,windows,**options)
            streamed=list(BasisStateSpace.states_in_windows(frequencies,windows,**options,
                                                           return_iterator=True,block_size=2))
        for actual,expected in zip(eager,reference):
            np.testing.assert_array_equal(actual,expected)
        self.assertTrue(streamed)
        self.assertFalse(hasattr(BasisStateSpace,'iter_states_in_windows'))

    def test_eager_sparse_can_construct_sparse_basis_without_dense_cache(self):
        from Psience.BasisReps import HarmonicOscillatorProductBasis
        sparse=BasisStateSpace.states_in_windows([1.,2.,3.],[[0.,4.]],
                                                max_quantas=3,sparse=True)[0]
        basis=HarmonicOscillatorProductBasis(3)
        space=BasisStateSpace(basis,sparse)
        self.assertTrue(space.has_sparse_excitations)
        self.assertIs(space._sparse_excitations,sparse)
        self.assertIsNone(space._excitations)
        expected=BasisStateSpace.states_in_windows([1.,2.,3.],[[0.,4.]],max_quantas=3)[0]
        np.testing.assert_array_equal(sparse.asarray(),expected)

    def test_sparse_traversal_handles_negative_frequencies_with_finite_bounds(self):
        for sparse in [False,True]:
            actual=BasisStateSpace.states_in_windows([-1.,2.],[[0.,1.]],max_quantas=4,
                                                     sparse=sparse,sparse_traversal=True)[0]
            values=actual.asarray() if sparse else actual
            expected={(a,b) for a in range(5) for b in range(5-a) if 0<=-a+2*b<=1}
            self.assertEqual(set(map(tuple,values)),expected)
