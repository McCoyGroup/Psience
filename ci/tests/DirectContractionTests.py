import itertools,unittest
import numpy as np
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace,BraKetSpace
from Psience.BasisReps.DirectContractions import QuadraticProductCoefficients

class DirectTests(unittest.TestCase):
    def test_words_and_asymmetric_coefficients(self):
        rng=np.random.default_rng(427)
        basis=HarmonicOscillatorProductBasis(3)
        b=rng.integers(0,4,size=(140,3));k=rng.integers(0,4,size=(140,3))
        b[:30]=k[:30]
        pairs=BraKetSpace(BasisStateSpace(basis,b),BasisStateSpace(basis,k))
        for word in [('x',),('x','x'),('p','p'),('x','x','x'),('p','x','p'),
                     ('x','x','x','x'),('x','p','x','p'),('p','x','x','p')]:
            with self.subTest(word=word):
                c=rng.normal(size=(3,)*len(word));c.flat[::3]=0
                out=[]
                for method in ['tensor','direct']:
                    v=basis.representation(*word,coeffs=c,contraction_method=method,
                                            direct_block_size=37).get_brakets(pairs)
                    out.append(v.asarray() if hasattr(v,'asarray') else v)
                np.testing.assert_allclose(*out,atol=3e-12,rtol=3e-12)

    def test_component_axes_and_axis_permutation(self):
        rng=np.random.default_rng(6);basis=HarmonicOscillatorProductBasis(2)
        states=BasisStateSpace(basis,np.array(list(itertools.product(range(3),repeat=2))))
        i,j=np.indices((len(states),len(states)))
        pairs=BraKetSpace(states.take_subspace(i.ravel()),states.take_subspace(j.ravel()))
        for word,c,axes in [(('x','x','x'),rng.normal(size=(2,2,2,3)),[[0,1,2],[0,1,2]]),
                            (('p','x','p'),rng.normal(size=(2,2,2)),[[0,1,2],[1,0,2]])]:
            out=[]
            for method in ['tensor','direct']:
                v=basis.representation(*word,coeffs=c,axes=axes,contraction_method=method).get_brakets(pairs)
                out.append(v.asarray() if hasattr(v,'asarray') else v)
            np.testing.assert_allclose(*out,atol=3e-12)

    def test_factored_quartic(self):
        rng=np.random.default_rng(8);basis=HarmonicOscillatorProductBasis(3)
        states=BasisStateSpace(basis,rng.integers(0,4,size=(100,3)))
        pairs=BraKetSpace(states,states.take_subspace(np.arange(100)[::-1]))
        factors=QuadraticProductCoefficients(rng.normal(size=(3,3,3)),np.array([.2,.7,-.1]))
        out=[]
        for c,method in [(factors.asarray(),'tensor'),(factors,'direct')]:
            v=basis.representation('x','p','x','p',coeffs=c,contraction_method=method).get_brakets(pairs)
            out.append(v.asarray() if hasattr(v,'asarray') else v)
        np.testing.assert_allclose(*out,atol=3e-12)

    def test_coriolis_factors_reproduce_terms(self):
        from McUtils.Scaffolding import NullCheckpointer
        from Psience.VPT2.Terms import CoriolisTerm
        from Psience.VPT2.Hamiltonian import PerturbationTheoryHamiltonian
        class Terms(CoriolisTerm):
            @property
            def freqs(self):return np.array([.03,.04,.05])
            def _freq_sqrt(self,freqs):return np.sqrt(freqs)
            def get_zetas_and_momi(self):return np.arange(27).reshape(3,3,3)*.003,np.array([.1,.2,.3])
            def moment_of_inertia_derivs(self,order):return [2*np.diag([.1,.2,.3])]
        terms=object.__new__(Terms);terms.checkpointer=NullCheckpointer()
        getter=PerturbationTheoryHamiltonian.CoriolisTermGetter(terms,None)
        actual=getter.get_factored_term(0).asarray()
        raw=terms.get_terms(0)[0];expected=raw[0,0]+raw[1,1]+raw[2,2]
        np.testing.assert_allclose(actual,expected,atol=1e-15)
        self.assertIsNone(PerturbationTheoryHamiltonian.CoriolisTermGetter(terms, []).get_factored_term(1))

    def test_empty_and_lazy_braket_length(self):
        basis=HarmonicOscillatorProductBasis(2)
        states=BasisStateSpace(basis,np.arange(10),mode='indices')
        pairs=BraKetSpace(states,states)
        self.assertEqual(len(pairs),10)
        self.assertIsNone(pairs._state_pairs)
        self.assertIsNone(states._excitations)

    def test_mismatched_mode_count_rejected(self):
        basis=HarmonicOscillatorProductBasis(3)
        wrong=HarmonicOscillatorProductBasis(2)
        states=BasisStateSpace(wrong,[[0,0],[1,0]])
        rep=basis.representation('x','x',coeffs=np.eye(3),contraction_method='direct')
        with self.assertRaisesRegex(ValueError,'mode count'):
            rep.get_brakets(BraKetSpace(states,states))

    def test_constant_scalar_and_component(self):
        basis=HarmonicOscillatorProductBasis(2)
        states=BasisStateSpace(basis,[[0,0],[1,0],[0,0]])
        pairs=BraKetSpace(states,states.take_subspace([0,2,1]))
        for c in [.3,np.array([.2,.1,-.4])]:
            out=[]
            for method in ['tensor','direct']:
                value=basis.representation(coeffs=c,axes=[[],[]],contraction_method=method).get_brakets(pairs)
                out.append(value.asarray() if hasattr(value,'asarray') else value)
            np.testing.assert_allclose(*out,atol=2e-14)

    @staticmethod
    def _molecular_quadratic_fixture():
        from pathlib import Path
        from Psience.VPT2.Hamiltonian import PerturbationTheoryHamiltonian
        path = Path(__file__).parent / 'TestData/HOH_freq.fchk'
        ham = PerturbationTheoryHamiltonian.from_fchk(
            str(path), logger=False, operator_contraction_method='direct')
        states = BasisStateSpace.from_quanta(ham.basis, range(3))
        return ham, states

    def test_molecular_quadratic_diagonal_shortcut_is_exact(self):
        ham, states = self._molecular_quadratic_fixture()
        representation = ham.get_perturbations(0)[0]
        values = []
        for enabled in [False, True]:
            for term in representation.computers:
                term.operator.use_analytic_harmonic_diagonal = enabled
            matrix = representation.get_representation_matrix(
                states, states, diagonal=True, zero_element_warning=False)
            values.append(matrix.asarray())
        self.assertEqual(values[0].tobytes(), values[1].tobytes())
        self.assertTrue(np.all(np.isfinite(values[1])))
        # Operators restored from older checkpoints lack the new option.
        for term in representation.computers:
            del term.operator.use_analytic_harmonic_diagonal
        restored = representation.get_representation_matrix(
            states, states, diagonal=True, zero_element_warning=False)
        self.assertEqual(values[1].tobytes(), restored.asarray().tobytes())

    def test_molecular_quadratic_diagonal_components_masks_and_thresholds(self):
        ham, states = self._molecular_quadratic_fixture()
        pairs = BraKetSpace(states, states)
        for word, raw in [(('x', 'x'), ham.V_terms[0]), (('p', 'p'), ham.G_terms[0])]:
            coefficients = np.stack([raw, raw.T], axis=1)
            for threshold, skipped in [(1e-14, None), (1., None), (1e-14, [(0, 0)])]:
                with self.subTest(word=word, threshold=threshold, skipped=skipped):
                    values = []
                    for enabled in [False, True]:
                        rep = ham.basis.representation(
                            *word, coeffs=coefficients, axes=[[1, 0], [2, 0]],
                            contraction_method='direct', zero_threshold=threshold,
                            skipped_indices=skipped, use_analytic_harmonic_diagonal=enabled)
                        values.append(rep.get_brakets(pairs))
                    self.assertEqual(values[0].tobytes(), values[1].tobytes())

    def test_molecular_mixed_diagonal_and_off_diagonal_pairs_are_exact(self):
        ham, states = self._molecular_quadratic_fixture()
        order = np.arange(len(states))[::-1].copy()
        order[::3] = np.arange(len(states))[::3]
        pairs = BraKetSpace(states, states.take_subspace(order))
        for word, raw in [(('x', 'x'), ham.V_terms[0]), (('p', 'p'), ham.G_terms[0])]:
            values = [ham.basis.representation(
                *word, coeffs=raw, contraction_method='direct',
                use_analytic_harmonic_diagonal=enabled).get_brakets(pairs)
                for enabled in [False, True]]
            self.assertEqual(values[0].tobytes(), values[1].tobytes())

    def test_run_simple_quadratic_diagonal_option_defaults_on_and_can_disable(self):
        from Psience.VPT2 import VPTRunner
        ham, states = self._molecular_quadratic_fixture()
        for options, expected in [({}, True), ({'use_analytic_harmonic_diagonal': True}, True),
                                  ({'use_analytic_harmonic_diagonal': False}, False)]:
            runner, _ = VPTRunner.construct(ham.molecule, states.excitations,
                                           logger=False, degeneracy_specs='auto', **options)
            self.assertIs(runner.hamiltonian.operator_settings['use_analytic_harmonic_diagonal'], expected)
            h0 = runner.hamiltonian.get_perturbations(0)[0]
            self.assertTrue(all(term.operator.use_analytic_harmonic_diagonal is expected
                                for term in h0.computers))

    def test_nonharmonic_basis_keeps_its_native_diagonal(self):
        from Psience.BasisReps import SimpleProductBasis, LegendreBasis
        ham, actual_states = self._molecular_quadratic_fixture()
        basis = SimpleProductBasis(LegendreBasis, [8, 8])
        states = BasisStateSpace(basis, actual_states.excitations[:, :2])
        pairs = BraKetSpace(states, states)
        for method in ['tensor', 'auto']:
            values = []
            for enabled in [False, True]:
                rep = basis.representation('x', 'x', coeffs=ham.V_terms[0][:2, :2],
                    contraction_method=method, use_analytic_harmonic_diagonal=enabled)
                value = rep.get_brakets(pairs)
                values.append(value.asarray() if hasattr(value, 'asarray') else value)
            self.assertEqual(values[0].tobytes(), values[1].tobytes())

if __name__=='__main__':unittest.main()
