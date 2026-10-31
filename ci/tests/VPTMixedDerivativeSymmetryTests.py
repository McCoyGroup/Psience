import itertools
import types
import unittest
from unittest.mock import patch
import numpy as np
from McUtils.Combinatorics import UniquePermutations
from McUtils.Scaffolding import NullLogger
from Psience.VPT2.Terms import PotentialTerms, ExpansionTerms
from Psience.VPT2 import AnalyticVPTRunner
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace


def sampled_tensor(size, rank, repeated=2):
    truth = np.zeros((size,) * rank)
    raw = np.zeros_like(truth)
    for number, key in enumerate(itertools.combinations_with_replacement(range(size), rank), 1):
        for pos in UniquePermutations(key).permutations():
            truth[tuple(pos)] = number
            if len(set(pos[:repeated])) == 1:
                raw[tuple(pos)] = number
    return raw, truth


class VPTMixedDerivativeSymmetryTests(unittest.TestCase):
    def sym(self, raw, mode='analytical', mode_axes=2, **opts):
        return PotentialTerms._symmetrize_mixed_derivatives(raw, mode, mode_axes,
            logger=NullLogger(), **opts)

    def test_all_available_quartic_orbits_are_completed(self):
        raw, truth = sampled_tensor(4, 4)
        for mode in ['analytical', 'numerical', 'averaged', 'old']:
            result = self.sym(raw, mode, term_id='v4_cart')
            for pos in np.ndindex(raw.shape):
                expected = truth[pos] if len(set(pos)) < 4 else 0
                self.assertEqual(result[pos], expected, (mode, pos))
        np.testing.assert_array_equal(raw, sampled_tensor(4, 4)[0])

    def test_repeated_index_can_be_smallest_middle_or_largest(self):
        for sample in [(0,0,1,2), (1,1,0,2), (2,2,0,1), (1,1,0,1)]:
            raw = np.zeros((3,) * 4)
            raw[sample] = raw[sample[:2] + sample[2:][::-1]] = 7
            for diagonal in [True, False]:
                result = self.sym(raw, diagonal=diagonal)
                for pos in UniquePermutations(sample).permutations():
                    self.assertEqual(result[tuple(pos)], 7)

    def test_two_mode_only_retains_the_legacy_iijj_subset(self):
        raw, truth = sampled_tensor(4, 4)
        result = self.sym(raw, two_mode_only=True)
        for pos in np.ndindex(raw.shape):
            counts = [pos.count(value) for value in set(pos)]
            self.assertEqual(result[pos], truth[pos] if all(count % 2 == 0 for count in counts) else 0)

    def test_selection_policy_and_zero_fallback(self):
        raw = np.zeros((2,) * 4)
        raw[0,0,1,1], raw[1,1,0,0] = 2, 8
        for mode, expected in [('analytical',2),('numerical',8),('averaged',5)]:
            self.assertEqual(self.sym(raw, mode)[0,1,0,1], expected)
        raw[0,0,1,1] = 0
        for mode in ['analytical','numerical','averaged']:
            self.assertEqual(self.sym(raw, mode)[0,1,0,1], 8)
        self.assertEqual(self.sym(raw, handle_zeros=False)[0,1,0,1], 0)
        self.assertEqual(self.sym(raw, 'averaged', handle_zeros=False)[0,1,0,1], 4)

    def test_value_axes_select_available_samples_independently(self):
        raw = np.zeros((2,3,2,2,2,2))
        raw[:,:,0,0,1,1] = [[0,2,0],[3,4,5]]
        raw[:,:,1,1,0,0] = [[8,9,7],[9,9,9]]
        result = self.sym(raw, val_axes=2)
        np.testing.assert_array_equal(result[:,:,0,1,0,1], [[8,2,7],[3,4,5]])

    def test_pre_transform_symmetry_stays_within_each_basis(self):
        raw = np.arange(2*2*3*3,dtype=float).reshape(2,2,3,3)
        result = self.sym(raw, term_id='v4_int_pre', handle_zeros=False)
        self.assertEqual(result[1,0,2,1], raw[0,1,1,2])
        np.testing.assert_array_equal(result, result.swapaxes(0,1))
        np.testing.assert_array_equal(result, result.swapaxes(2,3))

    def test_zero_rest_false_preserves_unsampled_orbits(self):
        raw = np.zeros((4,) * 4)
        raw[0,1,2,3] = 19
        self.assertEqual(self.sym(raw)[0,1,2,3], 0)
        self.assertEqual(self.sym(raw, zero_rest=False)[0,1,2,3], 19)

    def test_unhandled_keeps_sample_values_without_inventing_symmetry(self):
        raw = np.zeros((3,) * 4);raw[1,1,0,2] = 7
        result = self.sym(raw, 'unhandled')
        np.testing.assert_array_equal(result, raw)

    def test_higher_rank_sampled_symmetry(self):
        raw, truth = sampled_tensor(3, 5)
        np.testing.assert_array_equal(self.sym(raw), truth)

    def test_dipole_semidiagonal_cubic_keeps_repeated_high_index(self):
        raw = np.zeros((3,3,3));raw[2,2,0] = 9
        result = self.sym(raw)
        for pos in UniquePermutations((2,2,0)).permutations():
            self.assertEqual(result[tuple(pos)], 9)
        self.assertEqual(result[0,1,2], 0)

    def test_quartic_filter_preserves_kinetic_index_roles(self):
        tensor = np.arange(4**4).reshape((4,) * 4) + 1
        for two_mode in [False, True]:
            v, g, extra = ExpansionTerms.restrict_quartic_expansion([tensor,2*tensor,17],two_mode_only=two_mode)
            for pos in np.ndindex(tensor.shape):
                counts = [pos.count(value) for value in set(pos)]
                keep = all(count % 2 == 0 for count in counts) if two_mode else max(counts) >= 2
                self.assertEqual(v[pos], tensor[pos] if keep else 0)
                self.assertEqual(g[pos], 2*tensor[pos] if keep else 0)
            self.assertEqual(extra,17)
        result = ExpansionTerms.restrict_quartic_expansion([np.full((1,)*4,7.),0])
        self.assertEqual(result[0].item(),7.)
        self.assertEqual(result[1],0)

    def test_runner_options_select_semidiagonal_legacy_or_full_input(self):
        tensor = np.ones((4,)*4)
        ham = types.SimpleNamespace(molecule=types.SimpleNamespace(internals=None),
            logger=NullLogger(),checkpointer=None,local_mode_couplings=False,local_mode_coupling_order=None)
        ham._get_expansion_orders = lambda *a: {}
        ham.get_perturbations = lambda *a, **kw: [[np.eye(4),np.eye(4)],[np.zeros((4,)*3),0],[tensor,tensor,0,0]]
        class Capture:
            def __init__(self,expansions,*a,**kw):self.expansions=expansions
        construct = AnalyticVPTRunner.from_hamiltonian.__func__
        default = construct(Capture,ham,2).expansions[2][0]
        self.assertEqual(default[0,1,1,2],1)
        self.assertEqual(default[0,1,2,3],0)
        legacy = construct(Capture,ham,2,two_mode_only=True).expansions[2][0]
        self.assertEqual(legacy[0,1,1,2],0)
        self.assertEqual(legacy[0,1,0,1],1)
        full = construct(Capture,ham,2,take_diagonal_v4_terms=False).expansions[2][0]
        self.assertEqual(full[0,1,2,3],1)

    def test_public_entrypoints_forward_analytic_filter_options(self):
        class Capture(AnalyticVPTRunner):
            @classmethod
            def from_hamiltonian(cls, *args, **opts):
                return opts
        classic=types.SimpleNamespace(hamiltonian=object(),pt_opts=types.SimpleNamespace(opts={'order':2}))
        for method in [Capture.construct,Capture.from_file]:
            with patch('Psience.VPT2.Runner.VPTRunner.construct',return_value=(classic,None)) as native:
                result=method('synthetic',take_diagonal_v4_terms=False,two_mode_only=True)
            self.assertFalse(result['take_diagonal_v4_terms'])
            self.assertTrue(result['two_mode_only'])
            self.assertNotIn('two_mode_only',native.call_args.kwargs)
            self.assertNotIn('take_diagonal_v4_terms',native.call_args.kwargs)

    def test_preserved_mixed_quartic_supplies_fundamental_coupling(self):
        raw = np.zeros((3,)*4);raw[1,1,0,2] = raw[1,1,2,0] = .008
        completed = self.sym(raw)
        basis = HarmonicOscillatorProductBasis(3)
        states = BasisStateSpace(basis,[[1,0,0],[0,0,1]])
        def matrix(tensor):
            rep = basis.representation('x','x','x','x',coeffs=tensor/24,axes=[[0,1,2,3],[0,1,2,3]])
            value = rep.get_representation_matrix(states,states,zero_element_warning=False)
            return value.asarray() if hasattr(value,'asarray') else np.asarray(value)
        self.assertAlmostEqual(matrix(ExpansionTerms.restrict_quartic_expansion([completed,0])[0])[0,1],.008/8)
        self.assertEqual(matrix(ExpansionTerms.restrict_quartic_expansion([completed,0],two_mode_only=True)[0])[0,1],0)


if __name__ == '__main__':
    unittest.main()
