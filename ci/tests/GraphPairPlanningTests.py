"""Exact directed-graph checks, independent of VPT matrix coefficients."""
import itertools
import unittest
import numpy as np
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace
from Psience.BasisReps.ExternalStateSpaces import (
    StateSpaceDiskStore, DiskBasisStateSpace, DiskSelectionRuleStateSpace,
    _StoredCouplings)


def collect(graph, other=None):
    blocks=list(graph.iter_pair_blocks(filter_space=other))
    assert all(len(b)<=graph.store.block_size for b in blocks)
    pairs=np.concatenate(blocks) if blocks else np.empty((0,2),dtype=np.int64)
    return sorted(map(tuple,pairs.tolist()))


def brute(basis, source, rules, targets=None, dimensions=None, excluded=()):
    result=[]
    for bra in source:
        original=basis.unravel_state_inds(np.array([bra]))[0]
        dims=list(range(basis.ndim)) if dimensions is None else dimensions
        changes=set()
        for rule in rules:
            if len(rule)<=len(dims):
                changes.update(itertools.permutations(tuple(rule)+(0,)*(len(dims)-len(rule))))
        for change in changes:
            output=original.copy();output[dims]+=change
            if np.any(output<0):continue
            if tuple(output-original) in excluded:continue
            ket=int(basis.ravel_state_inds(output[np.newaxis,:])[0])
            if targets is None or ket in targets:result.append((int(bra),ket))
    return sorted(result)


class GraphPairPlanningTests(unittest.TestCase):
    def setUp(self):
        self.basis=HarmonicOscillatorProductBasis([8]*3)
        self.store=StateSpaceDiskStore(block_size=7)
        self.addCleanup(self.store.close)
        occupations=np.array([s for s in itertools.product(range(5),repeat=3) if sum(s)<=4])
        self.source=self.store.coerce(BasisStateSpace(self.basis,occupations))
        self.small=self.store.coerce(BasisStateSpace(self.basis,[[0,0,0],[1,0,0]]))

    def test_reverse_preserves_one_sided_direction(self):
        rules=[(1,), (1,1), (2,)]
        graph=DiskSelectionRuleStateSpace(self.source,rules,filter_space=self.small)
        self.assertEqual(collect(graph),brute(self.basis,self.source.indices,rules,set(self.small.indices)))
        self.assertLessEqual(self.store.metrics['peak_source_states'],len(self.small))
        self.assertIsNone(self.source._lookup)

    def test_two_target_filters_and_canonical_pairs(self):
        rules=[(),(-1,),(1,),(-1,1),(-1,-1,1),(1,1,1)]
        subset=self.store.coerce(BasisStateSpace(self.basis,[[0,0,0],[0,1,0],[1,0,0]]))
        graph=DiskSelectionRuleStateSpace(self.source,rules,filter_space=subset)
        expected=brute(self.basis,self.source.indices,rules,set(subset.indices)&set(self.small.indices))
        self.assertEqual(collect(graph,self.small),expected)
        self.assertEqual(sorted(map(tuple,graph.canonical_pairs(self.small).tolist())),
                         sorted(set(tuple(sorted(pair)) for pair in expected)))

    def test_subset_modes_and_signed_exclusion(self):
        rules=[(-1,), (1,), (-1,1)]
        excluded=((0,1,0),)
        graph=DiskSelectionRuleStateSpace(self.source,rules,filter_space=self.small,
                                         target_dimensions=[0,1]).filter_transitions(excluded)
        self.assertEqual(collect(graph),brute(self.basis,self.source.indices,rules,
                         set(self.small.indices),[0,1],excluded))

    def test_none_keeps_original_cartesian_target_domain(self):
        graph=DiskSelectionRuleStateSpace(self.small,None)
        outside=self.store.coerce(BasisStateSpace(self.basis,[[0,0,0],[2,0,0]]))
        wanted=[(int(s),int(t)) for s in self.small.indices for t in self.small.indices if t in outside.indices]
        self.assertEqual(collect(graph,outside),sorted(wanted))

    def test_noncanonical_bounds_and_duplicate_sources(self):
        ranks=np.array([*self.source.indices[::-1],self.source.indices[0]])
        source=DiskBasisStateSpace(self.basis,ranks,self.store,canonical=False,unique=False)
        self.assertEqual(source.quanta_bounds,(0,4))
        graph=DiskSelectionRuleStateSpace(source,[(-1,),(1,)],filter_space=self.small)
        self.assertEqual(collect(graph),brute(self.basis,ranks,[(-1,),(1,)],set(self.small.indices)))

    def test_owned_records_keep_outer_ownership(self):
        ranks=self.source.indices
        records=np.array([[ranks[0],ranks[3],ranks[1]],
                          [ranks[1],ranks[4],ranks[0]],
                          [ranks[2],ranks[0],ranks[0]]])
        graph=DiskSelectionRuleStateSpace(self.small,recipes=[
            (self.small,_StoredCouplings(records),None,None,(),())])
        expected=[tuple(p[1:]) for p in records if p[0] in self.small.indices and p[2] in self.small.indices]
        self.assertEqual(collect(graph,self.small),sorted(expected))

    def test_empty_filter_and_single_row_blocks(self):
        empty=self.store.coerce(BasisStateSpace(self.basis,[],mode='indices'))
        graph=DiskSelectionRuleStateSpace(self.source,[(1,),(-1,)])
        self.assertEqual(collect(graph,empty),[])
        self.store.block_size=1
        self.assertEqual(collect(graph,self.small),brute(self.basis,self.source.indices,
                         [(1,),(-1,)],set(self.small.indices)))


if __name__=='__main__':unittest.main()
