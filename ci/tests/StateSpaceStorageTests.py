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

    def _frontier_cache_reference(self, graph):
        def ranks():
            yield from graph.representative_space.iter_index_blocks()
            for pairs in graph.iter_pair_blocks():yield pairs.reshape(-1)
        return np.array(graph.store.frontier(graph.basis,ranks()).indices,copy=True)

    def test_frontier_cache_union_replays_only_uncached_recipes(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=store.frontier(self.basis,[self.ranks])
                first=base.apply_selection_rules([(1,)],target_dimensions=np.array([0,1,70,71]))
                second=base.apply_selection_rules([(-1,),()],target_dimensions=np.array([0,1,70,71]))
                first.to_single();combined=first.union(second)
                expected=self._frontier_cache_reference(combined)
                cached={id(r) for r in first.recipes};seen=[]
                original=DiskSelectionRuleStateSpace.iter_pair_blocks
                def observed(graph,*args,**kwargs):
                    seen.extend(id(r) for r in graph.recipes)
                    yield from original(graph,*args,**kwargs)
                with patch.object(DiskSelectionRuleStateSpace,'iter_pair_blocks',observed):
                    result=combined.to_single()
                np.testing.assert_array_equal(result.indices,expected)
                self.assertTrue(seen);self.assertFalse(cached.intersection(seen))
                self.assertEqual(seen,[id(r) for r in second.recipes])
                self.assertIsNone(combined._frontier_components)
                # The complete directed edge recipes remain available.
                actual=list(map(tuple,np.concatenate(list(combined.iter_pair_blocks()))))
                ordinary=DiskSelectionRuleStateSpace(combined.representative_space,recipes=combined.recipes)
                reference=list(map(tuple,np.concatenate(list(ordinary.iter_pair_blocks()))))
                self.assertEqual(actual,reference)
            finally:store.close()

    def test_frontier_cache_both_cached_are_lazy_and_need_no_edges(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=store.frontier(self.basis,[self.ranks])
                first=base.apply_selection_rules([(1,)],target_dimensions=np.array([70,71]))
                second=base.apply_selection_rules([(-1,),()],target_dimensions=np.array([70,71]))
                left=first.to_single();right=second.to_single()
                expected=np.union1d(left.indices,right.indices)
                with patch.object(DiskSelectionRuleStateSpace,'iter_pair_blocks',side_effect=AssertionError('cached edges replayed')):
                    combined=first.union(second)
                    self.assertIsNone(combined._frontier)
                    self.assertIsNotNone(combined._frontier_components)
                    np.testing.assert_array_equal(combined.to_single().indices,expected)
            finally:store.close()

    def test_frontier_cache_uncached_union_keeps_single_replay_path(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=store.frontier(self.basis,[self.ranks])
                first=base.apply_selection_rules([(1,)],target_dimensions=np.array([70,71]))
                second=base.apply_selection_rules([(-1,),()],target_dimensions=np.array([70,71]))
                combined=first.union(second)
                self.assertIsNone(combined._frontier_components)
                expected=self._frontier_cache_reference(combined)
                np.testing.assert_array_equal(combined.to_single().indices,expected)
            finally:store.close()

    def test_frontier_cache_nested_unions_preserve_recipe_contributions(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=store.frontier(self.basis,[self.ranks])
                first=base.apply_selection_rules([(1,)],target_dimensions=np.array([70,71]))
                middle=base.apply_selection_rules([(-1,)],target_dimensions=np.array([70,71]))
                last=base.apply_selection_rules([()],target_dimensions=np.array([70,71]))
                first.to_single();last.to_single()
                combined=first.union(middle).union(last)
                self.assertEqual(len(combined._frontier_components),3)
                expected=self._frontier_cache_reference(combined)
                seen=[];original=DiskSelectionRuleStateSpace.iter_pair_blocks
                def observed(graph,*args,**kwargs):
                    seen.extend(id(r) for r in graph.recipes)
                    yield from original(graph,*args,**kwargs)
                with patch.object(DiskSelectionRuleStateSpace,'iter_pair_blocks',observed):
                    np.testing.assert_array_equal(combined.to_single().indices,expected)
                self.assertEqual(seen,[id(r) for r in middle.recipes])
            finally:store.close()

    def test_frontier_cache_pending_snapshot_survives_original_recipe_edits(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=store.frontier(self.basis,[self.ranks])
                pending=base.apply_selection_rules([(1,)],target_dimensions=np.array([70,71]))
                cached=base.apply_selection_rules([()],target_dimensions=np.array([70,71]));cached.to_single()
                combined=cached.union(pending);expected=self._frontier_cache_reference(combined)
                excluded=np.zeros(self.basis.ndim,dtype=int);excluded[-1]=1
                pending.filter_transitions([excluded],in_place=True)
                np.testing.assert_array_equal(combined.to_single().indices,expected)
                self.assertIsNot(combined.recipes[-1],pending.recipes[-1])
            finally:store.close()

    def test_frontier_cache_in_place_exclusions_invalidate_all_node_caches(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=store.frontier(self.basis,[np.array([self.ranks.min()],dtype=np.int64)])
                first=base.apply_selection_rules([(1,)],target_dimensions=np.array([70,71]));first.to_single()
                second=base.apply_selection_rules([()],target_dimensions=np.array([70,71]));second.to_single()
                for materialize in [False,True]:
                    combined=first.union(second)
                    if materialize:combined.to_single()
                    excluded=np.zeros(self.basis.ndim,dtype=int);excluded[-1]=1
                    combined.filter_transitions([excluded],in_place=True)
                    self.assertIsNone(combined._frontier);self.assertIsNone(combined._frontier_components)
                    self.assertIsNone(combined._frontier_block_size)
                    expected=self._frontier_cache_reference(combined)
                    np.testing.assert_array_equal(combined.to_single().indices,expected)
            finally:store.close()

    def test_frontier_cache_map_and_source_intersection_recompute_nodes(self):
        from Psience.BasisReps import SelectionRuleStateSpace
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=store.frontier(self.basis,[self.ranks])
                first=base.apply_selection_rules([(1,)],target_dimensions=np.array([70,71]));first.to_single()
                second=base.apply_selection_rules([()],target_dimensions=np.array([70,71]))
                combined=first.union(second)
                selected=base.take_subspace([0])
                source_filter=SelectionRuleStateSpace(selected,[],ignore_shapes=True)
                variants=[combined.map(lambda block:block.take_subspace([0])),
                          combined.intersection(source_filter,handle_subspaces=False)]
                for graph in variants:
                    self.assertIsNone(graph._frontier_components)
                    expected=self._frontier_cache_reference(graph)
                    np.testing.assert_array_equal(graph.to_single().indices,expected)
            finally:store.close()

    def test_frontier_cache_other_store_preserves_callback_block_context(self):
        with tempfile.TemporaryDirectory() as directory:
            first_store=StateSpaceDiskStore(directory,block_size=101)
            other_store=StateSpaceDiskStore(directory,block_size=5)
            try:
                left=first_store.frontier(self.basis,[self.ranks])
                right=other_store.frontier(self.basis,[self.ranks])
                first=left.apply_selection_rules([()]);first.to_single()
                second=right.apply_selection_rules([(1,)],target_dimensions=np.array([0,1,70,71])).map(
                    lambda block:block.take_subspace([0]))
                second.to_single()
                combined=first.union(second)
                self.assertIsNone(combined._frontier_components[1][0])
                expected=self._frontier_cache_reference(combined)
                np.testing.assert_array_equal(combined.to_single().indices,expected)
            finally:first_store.close();other_store.close()

    def test_frontier_cache_changed_block_size_does_not_reuse_old_callbacks(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=5)
            try:
                base=store.frontier(self.basis,[self.ranks])
                first=base.apply_selection_rules([(1,)],target_dimensions=np.array([0,1,70,71])).map(
                    lambda block:block.take_subspace([0]))
                second=base.apply_selection_rules([()]);first.to_single();second.to_single()
                store.block_size=101
                combined=first.union(second)
                self.assertIsNone(combined._frontier_components)
                expected=self._frontier_cache_reference(combined)
                np.testing.assert_array_equal(combined.to_single().indices,expected)
            finally:store.close()

    def test_frontier_cache_releases_intermediate_frontier_owners(self):
        import gc,weakref
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=store.frontier(self.basis,[self.ranks])
                first=base.apply_selection_rules([(1,)],target_dimensions=np.array([70,71]));first.to_single()
                second=base.apply_selection_rules([()],target_dimensions=np.array([70,71]))
                previous=weakref.ref(first._frontier)
                combined=first.union(second);del first
                gc.collect();self.assertIsNotNone(previous())
                combined.to_single();gc.collect()
                self.assertIsNone(combined._frontier_components)
                self.assertIsNone(previous())
            finally:store.close()

    def test_frontier_cache_explicit_nested_polyad_nodes_are_preserved(self):
        from Psience.BasisReps import SelectionRuleStateSpace
        source=BasisStateSpace(self.basis,self.ranks[:1],mode='indices')
        internal=BasisStateSpace(self.basis,self.ranks[1:2],mode='indices')
        target=BasisStateSpace(self.basis,self.ranks[2:3],mode='indices')
        resident=SelectionRuleStateSpace(source,[SelectionRuleStateSpace(internal,[target])])
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                first=DiskSelectionRuleStateSpace.from_resident(resident,store);first.to_single()
                base=store.frontier(self.basis,[self.ranks[:2]])
                second=base.apply_selection_rules([(1,)],target_dimensions=np.array([70,71]))
                combined=first.union(second);expected=self._frontier_cache_reference(combined)
                np.testing.assert_array_equal(combined.to_single().indices,expected)
                self.assertIn(int(self.ranks[1]),combined.indices)
            finally:store.close()

    def _class_frontier_reference(self,graph):
        def blocks():
            yield from graph.representative_space.iter_index_blocks()
            for pairs in graph.iter_pair_blocks():yield pairs.reshape(-1)
        return np.array(graph.store.frontier(graph.basis,blocks()).indices,copy=True)

    def _class_fundamentals(self,store):
        # The full fundamental manifold in the actual 72-mode Pyrene basis.
        group=self.basis.indexer.symmetric_group
        data,offsets=group._get_partition_perms([1])
        return store.frontier(self.basis,[np.arange(int(offsets[0]),int(offsets[0])+int(data[0].num_elements))])

    def test_class_frontier_coalesces_split_support_per_rule(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=self._class_fundamentals(store);values=base.indices;middle=len(values)//2
                first=store.frontier(self.basis,[values[:middle]])
                second=store.frontier(self.basis,[values[middle:]])
                graph=DiskSelectionRuleStateSpace(base,recipes=[(s,[(1,),(-1,)],None,None,(),()) for s in [first,second]])
                expected=self._class_frontier_reference(graph)
                with patch.object(DiskSelectionRuleStateSpace,'iter_pair_blocks',side_effect=AssertionError('full class generated edges')):
                    result=graph.to_single()
                np.testing.assert_array_equal(result.indices,expected)
                self.assertGreater(store.metrics['class_complete_sources'],0)
                self.assertEqual(store.metrics['class_pending_sources'],0)
                self.assertEqual(len(graph.recipes),2)
            finally:store.close()

    def test_class_frontier_does_not_cross_join_sources_and_rules(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=self._class_fundamentals(store);values=base.indices;middle=len(values)//2
                first=store.frontier(self.basis,[values[:middle]])
                second=store.frontier(self.basis,[values[middle:]])
                graph=DiskSelectionRuleStateSpace(base,recipes=[(first,[(1,)],None,None,(),()),(second,[(-1,)],None,None,(),())])
                expected=self._class_frontier_reference(graph)
                np.testing.assert_array_equal(graph.to_single().indices,expected)
            finally:store.close()

    def test_class_frontier_keeps_partial_source_support_exact(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                actual=self.ranks[np.sum(self.states,axis=1)<=1]
                base=store.frontier(self.basis,[np.r_[0,actual]])
                graph=DiskSelectionRuleStateSpace(base,[(1,),(-1,),()])
                expected=self._class_frontier_reference(graph)
                np.testing.assert_array_equal(graph.to_single().indices,expected)
                self.assertGreater(store.metrics.get('class_pending_sources',0),0)
            finally:store.close()

    def test_class_frontier_intersects_explicit_target_ID_filter(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=self._class_fundamentals(store)
                target=store.frontier(self.basis,[self.ranks])
                graph=DiskSelectionRuleStateSpace(base,[(1,),(-1,),()],filter_space=target)
                expected=self._class_frontier_reference(graph)
                with patch.object(DiskSelectionRuleStateSpace,'iter_pair_blocks',side_effect=AssertionError('complete filtered class generated edges')):
                    np.testing.assert_array_equal(graph.to_single().indices,expected)
            finally:store.close()

    def test_class_frontier_quanta_filter_matches_explicit_manifold(self):
        from Psience.BasisReps.ExternalStateSpaces import QuantaStateSpace
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=self._class_fundamentals(store);target=QuantaStateSpace(self.basis,2)
                explicit=store.frontier(self.basis,[np.arange(target.start,target.stop)])
                rules=[(1,),(-1,),()]
                expected=self._class_frontier_reference(DiskSelectionRuleStateSpace(base,rules,filter_space=explicit))
                graph=DiskSelectionRuleStateSpace(base,rules,filter_space=target)
                with patch.object(DiskSelectionRuleStateSpace,'iter_pair_blocks',side_effect=AssertionError('quanta class generated edges')):
                    np.testing.assert_array_equal(graph.to_single().indices,expected)
            finally:store.close()

    def test_class_frontier_preserves_block_dependent_transforms(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=self._class_fundamentals(store)
                graph=DiskSelectionRuleStateSpace(base,[(1,)]).map(lambda space:space.take_subspace(np.arange(0,len(space),2)))
                expected=self._class_frontier_reference(graph)
                np.testing.assert_array_equal(graph.to_single().indices,expected)
                self.assertNotIn('class_frontier_queries',store.metrics)
            finally:store.close()

    def test_class_frontier_preserves_mode_specific_edge_exclusions(self):
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=self._class_fundamentals(store)
                excluded=np.zeros(self.basis.ndim,dtype=int);excluded[0]=1
                graph=DiskSelectionRuleStateSpace(base,[(1,)]).filter_transitions([excluded])
                expected=self._class_frontier_reference(graph)
                np.testing.assert_array_equal(graph.to_single().indices,expected)
                self.assertNotIn('class_frontier_queries',store.metrics)
            finally:store.close()

    def test_class_quanta_filter_keeps_directed_edges_exact_after_frontier(self):
        from Psience.BasisReps.ExternalStateSpaces import QuantaStateSpace
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=self._class_fundamentals(store)
                for q in [0,2]:
                    implicit=QuantaStateSpace(self.basis,q)
                    explicit=store.frontier(self.basis,[np.arange(implicit.start,implicit.stop)])
                    rules=[(1,),(-1,),()]
                    graph=DiskSelectionRuleStateSpace(base,rules,filter_space=implicit)
                    reference=DiskSelectionRuleStateSpace(base,rules,filter_space=explicit)
                    expected=np.concatenate(list(reference.iter_pair_blocks()))
                    graph.to_single()
                    actual=np.concatenate(list(graph.iter_pair_blocks()))
                    order=lambda pairs:pairs[np.lexsort(pairs.T[::-1])]
                    np.testing.assert_array_equal(order(actual),order(expected))
            finally:store.close()

    def test_class_implicit_quanta_blocks_obey_requested_cap(self):
        from Psience.BasisReps.ExternalStateSpaces import QuantaStateSpace,index_blocks
        target=QuantaStateSpace(self.basis,1)
        blocks=list(index_blocks(target,7))
        self.assertTrue(all(len(block)<=7 for block in blocks))
        np.testing.assert_array_equal(np.concatenate(blocks),np.arange(target.start,target.stop))
        with self.assertRaises(RuntimeError):target.indices


    def test_structured_Pyrene_class_filters_accept_subtrees_without_leaves(self):
        from McUtils.Combinatorics import UniquePermutations
        # This five-quanta permutation class is part of the recorded complete
        # 19,757,815-state Pyrene frontier, not a substitute molecular benchmark.
        partition=np.zeros(self.basis.ndim,dtype=int);partition[:5]=1
        generator=UniquePermutations(partition)
        caps=np.ones(self.basis.ndim,dtype=int);caps[0]=0
        statistics={};ranges=list(generator.iter_index_ranges({'maximum':caps},block_size=3,statistics=statistics))
        self.assertEqual(statistics['leaf_checks'],0)
        self.assertLess(statistics['visited_nodes'],10)
        self.assertEqual(sum(int(b-a) for block in ranges for a,b in block),
                         int(generator.num_permutations)*(self.basis.ndim-5)//self.basis.ndim)

    def test_structured_graph_nodes_and_edges_match_pointwise_filters(self):
        from Psience.BasisReps.StateFilters import Postfilter
        with tempfile.TemporaryDirectory() as directory:
            store=StateSpaceDiskStore(directory,block_size=101)
            try:
                base=self._class_fundamentals(store)
                caps=np.full(self.basis.ndim,2);caps[0]=0;caps[-1]=1
                specs=[{'max_quanta':caps},{'max_quanta':2,'parity':'even','partitions':[(),(2,),(1,1)]}]
                for spec in specs:
                    post=Postfilter.construct(base,spec)
                    graph=DiskSelectionRuleStateSpace(base,[(1,),(-1,),()]).map(post.apply)
                    expected=self._class_frontier_reference(graph)
                    actual=graph.to_single()
                    np.testing.assert_array_equal(actual.indices,expected)
                    self.assertGreater(store.metrics.get('constraint_visited_nodes',0),0)
                    # Frontier acceleration must not change the later edge recipes.
                    plain=DiskSelectionRuleStateSpace(base,[(1,),(-1,),()]).map(post.apply)
                    key=lambda pairs:pairs[np.lexsort(pairs.T[::-1])]
                    np.testing.assert_array_equal(key(np.concatenate(list(graph.iter_pair_blocks()))),
                                                  key(np.concatenate(list(plain.iter_pair_blocks()))))
            finally:store.close()

    def test_structured_run_simple_binds_selected_frequencies_without_mutation(self):
        from types import SimpleNamespace
        from inspect import signature
        from Psience.VPT2 import VPTRunner,VPTStateSpace
        from Psience.BasisReps import BasisStateSpaceFilter
        from McUtils.Data import UnitsData
        frequencies=np.array([0.0022091806829265537,0.0026866989790004497,
            0.0030336185319669155,0.0035333560019367226,0.0040757163706985075,
            0.005995761249610487,0.006121095689537073,0.008412070536502078,0.016940751588948696])
        spec={'max_quanta':2,'energy_windows':[[0,4000]],'energy_units':'Wavenumbers'}
        bound=BasisStateSpaceFilter.bind_postfilters(spec,modes=SimpleNamespace(freqs=frequencies),mode_selection=(8,0,4))
        np.testing.assert_array_equal(bound['frequencies'],frequencies[[8,0,4]])
        self.assertNotIn('frequencies',spec)
        self.assertIn('basis_filters',signature(VPTRunner.run_simple).parameters)
        # Explicit postfilters must also work for the unfiltered WFC target.
        filters=VPTStateSpace.get_state_space_filter(self.states[:2],target='wavefunctions',order=2,postfilters={'parity':'even'})
        self.assertTrue(filters)
        self.assertTrue(all(f.postfilters for f in filters.values()))

    def test_sparse_window_stream_preserves_recorded_pyrene_states(self):
        expected={tuple(state) for state in self.states if sum(state)<=2}
        predicate=lambda state,window:tuple(state) in expected
        blocks=list(BasisStateSpace.states_in_windows(
            np.ones(self.basis.ndim),[[0.,2.]],max_quantas=2,
            sparse=True,filter=predicate,block_size=3,max_bytes=1024, return_iterator=True))
        actual=[]
        for tags,block in blocks:
            np.testing.assert_array_equal(tags,np.zeros(len(tags),dtype=int))
            space=BasisStateSpace(self.basis,block)
            self.assertTrue(space.has_sparse_excitations)
            self.assertIsNone(space._excitations)
            self.assertFalse(space.keep_excitations)
            np.testing.assert_array_equal(space.indices,self.basis.ravel_state_inds(block.asarray()))
            actual.extend(map(tuple,space.sparse_excitations.asarray()))
        self.assertEqual(set(actual),expected)
        self.assertEqual(len(actual),len(expected))
        key=lambda row:(sum(row),tuple(i for i,n in enumerate(row) for _ in range(n)))
        self.assertEqual(actual,sorted(expected,key=key))

    def test_sparse_window_collection_never_decodes_dense_occupations(self):
        with patch.object(type(self.space().sparse_excitations),'asarray',
                          side_effect=AssertionError('dense occupation table requested')):
            result=BasisStateSpace.states_in_windows(
                np.ones(self.basis.ndim),[[0.,1.],[2.,1.]],
                max_quantas=1,sparse=True,block_size=3,max_bytes=1024)
            first=BasisStateSpace(self.basis,result[0])
            empty=BasisStateSpace(self.basis,result[1])
            self.assertTrue(first.has_sparse_excitations)
            self.assertEqual(len(first),self.basis.ndim+1)
            self.assertEqual(empty.sparse_excitations.shape,(0,self.basis.ndim))


if __name__=='__main__':unittest.main()
