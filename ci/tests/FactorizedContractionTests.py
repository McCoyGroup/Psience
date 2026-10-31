"""Independent native-tensor checks for exact coefficient product contractions."""
import itertools,unittest
from unittest.mock import patch
import numpy as np
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace,BraKetSpace
from Psience.BasisReps.DirectContractions import FactorizedCoefficients,QuadraticProductCoefficients

class GeneralFactorTests(unittest.TestCase):
    def pairs(self,nmodes=2,count=24):
        rng=np.random.default_rng(712)
        basis=HarmonicOscillatorProductBasis(nmodes)
        a=rng.integers(0,3,size=(count,nmodes));b=rng.integers(0,3,size=(count,nmodes))
        a[:count//3]=b[:count//3]
        return basis,BraKetSpace(BasisStateSpace(basis,a),BasisStateSpace(basis,b))

    def compare(self,word,coefficients,**opts):
        basis,pairs=self.pairs()
        dense=coefficients.asarray()
        reference=basis.representation(*word,coeffs=dense,contraction_method='tensor',**opts).get_brakets(pairs)
        reference=reference.asarray() if hasattr(reference,'asarray') else reference
        with patch.object(coefficients,'asarray',side_effect=AssertionError('unexpected coefficient densification')):
            actual=basis.representation(*word,coeffs=coefficients,contraction_method='direct',
                direct_block_size=7,**opts).get_brakets(pairs)
        np.testing.assert_allclose(actual,reference,atol=2e-10,rtol=2e-12)
        if len(word)>6:
            dense_direct=basis.representation(*word,coeffs=dense,contraction_method='direct',**opts).get_brakets(pairs)
            np.testing.assert_allclose(dense_direct,reference,atol=2e-10,rtol=2e-12)

    def test_ordered_words_and_degrees_above_six(self):
        rng=np.random.default_rng(122)
        for word in [('x','x','p','p'),('p','x','x','p'),
                     tuple('x'*6),tuple('xp'+'x'*4+'xp'),tuple('x'*8),
                     tuple('xp'+'x'*6+'xp')]:
            degree=len(word)
            factors=[rng.normal(size=(3,2)) for _ in word]
            coeff=FactorizedCoefficients(factors,[(i,) for i in range(degree)],rng.normal(size=3))
            with self.subTest(word=word):self.compare(word,coeff)

    def test_shared_axes_rank_broadcast_components_and_axis_permutation(self):
        rng=np.random.default_rng(21)
        coeff=FactorizedCoefficients([rng.normal(size=(3,2,2)),rng.normal(size=(1,2,2)),
            rng.normal(size=(3,2))],[(0,2),(1,2),(3,)],rng.normal(size=(3,3)))
        self.compare(tuple('pxxp'),coeff,axes=[[0,1,2,3],[2,0,3,1]])

    def test_explicit_exclusions_remain_guarded(self):
        rng=np.random.default_rng(18)
        coeff=FactorizedCoefficients([rng.normal(size=(2,2,2)),rng.normal(size=(2,2,2))],
            [(0,2),(1,3)],np.array([.2,-.4]))
        with self.assertRaisesRegex(ValueError,'unfiltered'):
            self.compare(tuple('pxxp'),coeff,axes=[[0,1,2,3],[2,0,3,1]],
                skipped_indices=np.array([[0,0,1,1],[1,0,0,0]]))

    def test_coefficient_threshold_remains_guarded(self):
        rng=np.random.default_rng(55)
        coeff=FactorizedCoefficients([rng.normal(size=(2,2,2)),rng.normal(size=(2,2,2))],
            [(0,1),(2,3)],np.array([.2,-.4]))
        with self.assertRaisesRegex(ValueError,'unfiltered'):
            self.compare(tuple('xpxp'),coeff,skipped_coefficient_threshold=.2)

    def test_fourteenth_degree_single_mode_has_no_degree_ceiling(self):
        basis=HarmonicOscillatorProductBasis(1)
        states=BasisStateSpace(basis,[[0]])
        coeff=FactorizedCoefficients([np.ones((1,1)) for _ in range(14)],[(i,) for i in range(14)])
        actual=basis.representation(*tuple('x'*14),coeffs=coeff,contraction_method='direct').get_brakets(BraKetSpace(states,states))
        # Gaussian ground-state moment: <q^14> = 13!! / 2^7.
        expected=np.prod(np.arange(1,14,2))/2**7
        np.testing.assert_allclose(actual,[expected],atol=2e-10)
        dense=basis.representation(*tuple('x'*14),coeffs=coeff.asarray(),contraction_method='direct').get_brakets(BraKetSpace(states,states))
        np.testing.assert_allclose(dense,[expected],atol=2e-10)

    def test_high_degree_dense_component_axes(self):
        rng=np.random.default_rng(29)
        coeff=FactorizedCoefficients([rng.normal(size=(2,2)) for _ in range(8)],
            [(i,) for i in range(8)],rng.normal(size=(2,3)))
        self.compare(tuple('xxxxxxpp'),coeff,axes=[list(range(8)),[2,0,7,3,4,1,6,5]])

    def test_zero_degree_component_and_empty_pairs(self):
        basis,pairs=self.pairs()
        coeff=FactorizedCoefficients([np.array([2.,3.])],[()],np.array([[.1,.2],[-.4,.3]]))
        actual=basis.representation(coeffs=coeff,axes=[[],[]],contraction_method='direct').get_brakets(pairs)
        bra,ket=pairs.state_pairs
        expected=np.all(bra==ket,axis=0)[:,None]*coeff.asarray()
        np.testing.assert_allclose(actual,expected,atol=1e-13)
        empty=pairs.take_subspace(slice(0,0))
        c=FactorizedCoefficients([np.ones((1,2,2))],[(0,1)])
        self.assertEqual(basis.representation('x','x',coeffs=c,contraction_method='direct').get_brakets(empty).shape,(0,))

    def test_higher_coriolis_factors_match_existing_native_terms(self):
        from Psience.VPT2.Terms import CoriolisTerm
        from Psience.VPT2.Hamiltonian import PerturbationTheoryHamiltonian
        from McUtils.Scaffolding import NullCheckpointer
        class Terms(CoriolisTerm):
            @property
            def freqs(self):return np.array([.03,.04])
            def _freq_sqrt(self,freqs):return np.sqrt(freqs)
            def get_zetas_and_momi(self):return np.arange(12).reshape(3,2,2)*.003,np.array([.1,.2,.3])
            def moment_of_inertia_derivs(self,order):
                return [np.arange(9*2**n).reshape((3,3)+(2,)*n)*.002 for n in range(order+1)]
        terms=object.__new__(Terms);terms.checkpointer=NullCheckpointer()
        getter=PerturbationTheoryHamiltonian.CoriolisTermGetter(terms,None)
        for order in range(5):
            raw=terms.get_terms(order)[order]
            expected=raw[0,0]+raw[1,1]+raw[2,2]
            actual=getter.get_factored_term(order)
            np.testing.assert_allclose(actual.asarray(),expected,atol=2e-17)
            self.compare(tuple('xp'+'x'*order+'xp'),actual)

    def test_validation(self):
        for factors,axes in [([np.ones((1,2))],[(1,)]),
                             ([np.ones((1,2,2))],[(0,0)]),
                             ([np.ones((2,2)),np.ones((3,2))],[(0,),(1,)])]:
            with self.assertRaises(ValueError):FactorizedCoefficients(factors,axes)
        with self.assertRaises(TypeError):FactorizedCoefficients([np.ones((1,2),complex)],[(0,)])

    def test_coefficient_updates_are_seen_on_repeated_calls(self):
        basis,pairs=self.pairs()
        c=FactorizedCoefficients([np.ones((1,2,2))],[(0,1)])
        rep=basis.representation('x','x',coeffs=c,contraction_method='direct')
        before=rep.get_brakets(pairs)
        c.weights*=2
        after=rep.get_brakets(pairs)
        np.testing.assert_allclose(after,2*before,atol=2e-13)

if __name__=='__main__':unittest.main()
