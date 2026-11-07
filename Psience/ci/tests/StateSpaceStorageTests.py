"""Sparse occupations in BasisStateSpace, using recorded Pyrene states."""
import json, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import scipy.sparse as sp
from McUtils.Numputils import SparseArray
from McUtils.Scaffolding import PseudoPickler
from Psience.BasisReps import BasisStateSpace, BasisStateSpaceView, HarmonicOscillatorProductBasis
from Psience.BasisReps.ExternalStateSpaces import StateSpaceDiskStore, DiskSelectionRuleStateSpace
from McUtils.Combinatorics import SymmetricGroupGenerator

class StateSpaceStorageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data=json.loads((Path(__file__).parent/'TestData/pyrene_vpt2_block_handoff.json').read_text())
        cls.states=np.asarray(data['state_space_probe_states'],dtype=np.int64)
        cls.basis=HarmonicOscillatorProductBasis(cls.states.shape[1])
        cls.ranks=cls.basis.ravel_state_inds(cls.states)

    def space(self):return BasisStateSpace(self.basis,sp.csr_matrix(self.states))

    def test_sparse_input_preserves_real_states_without_dense_cache(self):
        space=self.space()
        np.testing.assert_array_equal(space.indices,self.ranks)
        np.testing.assert_array_equal(space.sparse_excitations.asarray(),self.states)
        self.assertIsNone(space._excitations)
        self.assertIsInstance(space._init_states,SparseArray)

    def test_existing_sparse_wrapper_is_accepted_without_rewrapping(self):
        values=SparseArray.from_data(sp.csr_matrix(self.states))
        space=BasisStateSpace(self.basis,values)
        np.testing.assert_array_equal(space.indices,self.ranks)
        self.assertIs(space.sparse_excitations,values)

    def test_subsets_preserve_duplicates_empty_and_singleton_axes(self):
        space=self.space()
        for selected in [[2,0,2],[0],[],[-1,0]]:
            subset=space.take_subspace(selected)
            np.testing.assert_array_equal(subset.sparse_excitations.asarray(),self.states[selected])
            self.assertEqual(subset.sparse_excitations.shape,(len(selected),self.states.shape[1]))
            self.assertIsNone(subset._excitations)

    def test_index_to_sparse_conversion_is_blocked_and_not_retained(self):
        space=BasisStateSpace(self.basis,self.ranks,mode='indices',keep_excitations=False)
        with patch.object(self.basis,'unravel_state_inds',wraps=self.basis.unravel_state_inds) as decode:
            result=space.as_sparse_excitations(block_size=3)
        np.testing.assert_array_equal(result.asarray(),self.states)
        self.assertTrue(decode.call_count>1)
        self.assertTrue(all(len(call.args[0])<=3 for call in decode.call_args_list))
        self.assertIsNone(space._excitations)
        self.assertIsNone(space._sparse_excitations)

    def test_dense_input_can_be_discarded_after_indexing(self):
        space=BasisStateSpace(self.basis,self.states,keep_excitations=False)
        self.assertEqual(space._init_states.ndim,1)
        self.assertIsNone(space._excitations)
        np.testing.assert_array_equal(space.indices,self.ranks)
        np.testing.assert_array_equal(space.excitations,self.states)
        self.assertIsNone(space._excitations)

    def test_sparse_and_dense_queries_do_not_decode_the_stored_space(self):
        space=self.space()
        with patch.object(space,'as_excitations',side_effect=AssertionError('whole dense table requested')):
            np.testing.assert_array_equal(space.find(self.states[[4,1]]),[4,1])
            np.testing.assert_array_equal(space.find(sp.csr_matrix(self.states[[4,1]])),[4,1])

    def test_rank_set_operations_preserve_sparse_cache_policy(self):
        first=self.space().take_subspace([0,1,2,3])
        second=self.space().take_subspace([2,3,4,5])
        with patch.object(BasisStateSpace,'as_excitations',side_effect=AssertionError('whole dense table requested')):
            union=first.union(second,sort=True)
            intersection=first.intersection(second,sort=True)
            difference=first.difference(second)
            np.testing.assert_array_equal(union.indices,np.unique(self.ranks[:6]))
            np.testing.assert_array_equal(intersection.indices,np.intersect1d(self.ranks[:4],self.ranks[2:6]))
            self.assertEqual(set(map(int,difference.indices)),set(map(int,self.ranks[:2])))
            self.assertFalse(union.keep_excitations)

    def test_uncached_sparse_sets_never_select_a_dense_algorithm(self):
        # No cached indices: the old heuristics otherwise prefer occupations.
        with patch.object(BasisStateSpace,'as_excitations',side_effect=AssertionError('whole dense table requested')):
            for operation in ['union','intersection','difference']:
                first=BasisStateSpace(self.basis,sp.csr_matrix(self.states[:4]))
                second=BasisStateSpace(self.basis,sp.csr_matrix(self.states[2:6]))
                self.assertFalse(first.has_indices)
                result=getattr(first,operation)(second)
                self.assertFalse(result.keep_excitations)
                self.assertEqual(first.unique_len,4)

    def test_column_selection_uses_sparse_occupations(self):
        space=self.space()
        modes=np.flatnonzero(np.any(self.states,axis=0))
        with patch.object(space,'as_excitations',side_effect=AssertionError('whole dense table requested')):
            subset=space.take_subdimensions(modes)
        np.testing.assert_array_equal(subset.sparse_excitations.asarray(),self.states[:,modes])
        self.assertIsNone(subset._excitations)

    def test_native_serializer_preserves_sparse_storage(self):
        space=self.space();space.indices
        pickler=PseudoPickler();restored=pickler.deserialize(pickler.serialize(space))
        np.testing.assert_array_equal(restored.indices,self.ranks)
        np.testing.assert_array_equal(restored.sparse_excitations.asarray(),self.states)
        self.assertIsNone(restored._excitations)
        self.assertFalse(restored.keep_excitations)

    def test_empty_sparse_state_space(self):
        space=BasisStateSpace(self.basis,sp.csr_matrix((0,self.basis.ndim),dtype=np.int8))
        self.assertEqual(len(space),0)
        self.assertEqual(space.sparse_excitations.shape,(0,self.basis.ndim))
        self.assertEqual(space.indices.shape,(0,))

    def test_shared_view_and_mapped_back_conversion(self):
        parent=BasisStateSpace(self.basis,self.ranks,mode='indices',keep_excitations=False)
        positions=np.array([3,1,3])
        view=BasisStateSpaceView(self.basis,positions,mode='indices',index_parent=parent)
        np.testing.assert_array_equal(view.as_sparse_excitations(block_size=2).asarray(),self.states[positions])
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=3)
            try:
                mapped=store.frontier(self.basis,[self.ranks],sort=True)
                np.testing.assert_array_equal(mapped.as_sparse_excitations(block_size=2).asarray(),
                    self.basis.unravel_state_inds(np.sort(self.ranks)))
            finally:store.close()

    def test_shared_view_concatenation_returns_an_independent_space(self):
        parent=BasisStateSpace(self.basis,self.ranks,mode='indices',keep_excitations=False)
        first=BasisStateSpaceView(self.basis,np.array([3,1]),mode='indices',index_parent=parent)
        second=BasisStateSpaceView(self.basis,np.array([1,4]),mode='indices',index_parent=parent)
        with patch.object(BasisStateSpace,'as_excitations',side_effect=AssertionError('whole dense table requested')):
            result=first.concatenate(second)
        self.assertIs(type(result),BasisStateSpace)
        self.assertFalse(result.keep_excitations)
        np.testing.assert_array_equal(result.indices,self.ranks[[3,1,1,4]])

    def test_ordinary_dense_api_is_unchanged(self):
        space=BasisStateSpace(self.basis,self.states)
        self.assertIsInstance(space.excitations,np.ndarray)
        self.assertTrue(space.keep_excitations)
        np.testing.assert_array_equal(space.excitations,self.states)
        np.testing.assert_array_equal(space.indices,self.ranks)

    def test_sparse_shape_and_noninteger_input_rejected(self):
        with self.assertRaises(ValueError):
            BasisStateSpace(self.basis,sp.csr_matrix((2,3),dtype=np.int8))
        with self.assertRaises(TypeError):
            BasisStateSpace(self.basis,sp.csr_matrix(self.states.astype(float)+.5))

    def test_integral_float_sparse_and_sparse_setter(self):
        space=BasisStateSpace(self.basis,sp.csr_matrix(self.states.astype(float)))
        np.testing.assert_array_equal(space.indices,self.ranks)
        space.excitations=SparseArray.from_data(sp.csr_matrix(self.states))
        np.testing.assert_array_equal(space.sparse_excitations.asarray(),self.states)
        space.excitations=self.states
        self.assertFalse(space.has_sparse_excitations)
        np.testing.assert_array_equal(space.excitations,self.states)

    def test_split_and_concatenate_do_not_require_dense_occupations(self):
        for space in [self.space(),BasisStateSpace(self.basis,self.ranks,mode='indices',keep_excitations=False)]:
            with patch.object(BasisStateSpace,'as_excitations',side_effect=AssertionError('whole dense table requested')):
                parts=space.split(4)
                result=parts[0]
                for other in parts[1:]:result=result.concatenate(other)
                np.testing.assert_array_equal(result.indices,self.ranks)
                self.assertFalse(result.keep_excitations)

    def test_disk_rule_pairs_match_vector_reference_with_filters(self):
        # Actual Pyrene states, including overtones and combinations. Expected
        # pairs use the established vector iterator plus directed filtering.
        states=self.states
        group=SymmetricGroupGenerator(self.basis.ndim)
        expanded=[]
        for rows,values in group.iter_permutation_rule_direct_sum(states,[(1,)],block_size=101):
            expanded.extend(self.basis.ravel_state_inds(values).tolist())
        expanded=np.unique(expanded)
        excluded=np.zeros(self.basis.ndim,dtype=int);excluded[-1]=-1
        variants=[(self.ranks,None,[(1,),(-1,),()],None,[]),
                  (expanded,self.ranks,[(-1,)],None,[excluded]),
                  (self.ranks,None,[(1,),(-1,)],np.array([0,1,70,71]),[]),
                  (self.ranks,None,[(1,),(-1,)],None,[-excluded])]
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                for source_ids,target_ids,rules,dimensions,exclusions in variants:
                    source=store.frontier(self.basis,[source_ids],sort=True)
                    target=None if target_ids is None else store.frontier(self.basis,[target_ids],sort=True)
                    graph=DiskSelectionRuleStateSpace(source,rules,filter_space=target,target_dimensions=dimensions)
                    if exclusions:graph=graph.filter_transitions(exclusions)
                    blocks=list(graph.iter_pair_blocks());actual=np.concatenate(blocks) if blocks else np.empty((0,2),dtype=int)
                    vectors=self.basis.unravel_state_inds(source.indices).astype(int)
                    selected=vectors if dimensions is None else vectors[:,dimensions]
                    expected=[]
                    for rows,outputs in group.__class__(selected.shape[1]).iter_permutation_rule_direct_sum(selected,rules,block_size=101):
                        if dimensions is not None:
                            restored=vectors[rows].copy();restored[:,dimensions]=outputs;outputs=restored
                        valid=np.ones(len(rows),dtype=bool)
                        for rule in exclusions:valid &= ~np.all(outputs-vectors[rows]==rule,axis=1)
                        ids=self.basis.ravel_state_inds(outputs)
                        if target_ids is not None:valid &= np.isin(ids,target.indices)
                        expected.extend(zip(source.indices[rows[valid]],ids[valid]))
                    self.assertEqual(sorted(map(tuple,actual)),sorted(expected))
            finally:store.close()


    def test_degeneracy_coverage_preserves_real_order_and_duplicates(self):
        from Psience.VPT2.DegeneracySpecs import DegenerateMultiStateSpace
        chosen=[2,0,1,2,3]
        states=BasisStateSpace(self.basis,self.ranks[chosen],mode='indices',keep_excitations=False)
        result=DegenerateMultiStateSpace.from_spec(None,states=states)
        self.assertEqual([list(g.indices) for g in result.spaces.flat],[[self.ranks[i]] for i in chosen])

    def test_degeneracy_coverage_restores_filtered_states_once_in_input_order(self):
        from Psience.VPT2.DegeneracySpecs import DegenerateMultiStateSpace,GroupsDegeneracySpec
        chosen=[2,0,1,2,3]
        states=BasisStateSpace(self.basis,self.ranks[chosen],mode='indices',keep_excitations=False)
        # Overlapping groups and a member outside the requested input.
        spec=GroupsDegeneracySpec([self.states[[0,1,5]],self.states[[1,3]],self.states[[4]]])
        def filtered(group):return group.take_subspace(np.flatnonzero(group.indices!=self.ranks[1]))
        result=DegenerateMultiStateSpace.from_spec(spec,states=states,group_filter=filtered)
        expected=[self.ranks[[0,5]],self.ranks[[3]],self.ranks[[2]],self.ranks[[2]],self.ranks[[1]]]
        self.assertEqual([list(g.indices) for g in result.spaces.flat],[list(g) for g in expected])

    def test_degeneracy_first_matching_explicit_group_is_preserved(self):
        from Psience.VPT2.DegeneracySpecs import DegenerateMultiStateSpace,GroupsDegeneracySpec
        states=BasisStateSpace(self.basis,self.ranks[[1]],mode='indices',keep_excitations=False)
        spec=GroupsDegeneracySpec([self.states[[0,1]],self.states[[1,3]]])
        result=DegenerateMultiStateSpace.from_spec(spec,states=states)
        self.assertEqual([list(g.indices) for g in result.spaces.flat],[list(self.ranks[[0,1]])])

    def test_degeneracy_split_filter_retains_group_order(self):
        from Psience.VPT2.DegeneracySpecs import DegenerateMultiStateSpace,GroupsDegeneracySpec
        states=BasisStateSpace(self.basis,self.ranks[[2,0,1]],mode='indices',keep_excitations=False)
        spec=GroupsDegeneracySpec([self.states[[1,0]]])
        result=DegenerateMultiStateSpace.from_spec(spec,states=states,group_filter=lambda group:group.split(1))
        self.assertEqual([list(g.indices) for g in result.spaces.flat],[[self.ranks[1]],[self.ranks[0]],[self.ranks[2]]])

    def test_degeneracy_coverage_reads_shared_views_in_blocks(self):
        from Psience.VPT2.DegeneracySpecs import DegenerateMultiStateSpace
        parent=BasisStateSpace(self.basis,self.ranks,mode='indices',keep_excitations=False)
        def view(group):
            return BasisStateSpaceView(self.basis,parent.find(group.indices),mode='indices',index_parent=parent)
        with patch.object(BasisStateSpaceView,'as_indices',side_effect=AssertionError('whole shared-view ranks requested')):
            result=DegenerateMultiStateSpace.from_spec(None,states=parent,group_filter=view)
        self.assertEqual([list(g.indices) for g in result.spaces.flat],[[x] for x in self.ranks])

    def test_degeneracy_coverage_accesses_grow_linearly(self):
        from Psience.VPT2.DegeneracySpecs import DegenerateMultiStateSpace
        states=BasisStateSpace(self.basis,self.ranks,mode='indices',keep_excitations=False)
        original=BasisStateSpace._pull_inds;calls=[]
        def observed(space):calls.append(type(space));return original(space)
        with patch.object(BasisStateSpace,'_pull_inds',observed):
            result=DegenerateMultiStateSpace.from_spec(None,states=states)
        self.assertEqual(len(result.spaces),len(states))
        self.assertLess(len(calls),10*len(states)+20)

    def test_degeneracy_empty_filtered_groups_restore_missing_duplicates_once(self):
        from Psience.VPT2.DegeneracySpecs import DegenerateMultiStateSpace
        chosen=[2,0,2,1]
        states=BasisStateSpace(self.basis,self.ranks[chosen],mode='indices',keep_excitations=False)
        result=DegenerateMultiStateSpace.from_spec(None,states=states,group_filter=lambda group:group.take_subspace([]))
        self.assertEqual([list(g.indices) for g in result.spaces.flat],[[],[],[],[],[self.ranks[2]],[self.ranks[0]],[self.ranks[1]]])

if __name__=='__main__':unittest.main()
