import unittest
from types import SimpleNamespace
from unittest import mock
import numpy as np
from Psience.BasisReps import BasisStateSpace,HarmonicOscillatorProductBasis,PropertyCorrectionFilter
from Psience.VPT2 import VPTRunner,VPTSystem,VPTStateSpace
from Psience.VPT2.Runner import VPTHamiltonianOptions,VPTRuntimeOptions,VPTSolverOptions
from Psience.VPT2.Runner import resolve_evaluator_defaults
from Psience.VPT2.DegeneracySpecs import StronglyCoupledDegeneracySpec,DegeneracySpec
from Psience.VPT2.Corrections import PerturbationTheoryCorrections
from Psience.VPT2.Terms import ExpansionTerms
from McUtils.Numputils import SparseArray

class EvaluatorDefaultsTests(unittest.TestCase):
    def setUp(self):
        self.basis=HarmonicOscillatorProductBasis(2)
        self.targets=BasisStateSpace.from_quanta(self.basis,range(3))
        self.system=object.__new__(VPTSystem)
        self.system.mode_selection=[0,1];self.system.use_local_modes=False
        self.system.mol=SimpleNamespace(normal_modes=SimpleNamespace(modes=SimpleNamespace(freqs=np.array([.007,.011]))))
        self.space=VPTStateSpace(self.targets,system=self.system)
        self.initial=VPTStateSpace(BasisStateSpace(self.basis,[[0,0]]),system=self.system)

    def test_profiles_and_default_constructor(self):
        for profile in ['new','old']:
            runner,_=VPTRunner.construct(self.system,self.targets.excitations,evaluator_defaults=profile)
            opts=dict(runner.ham_opts.opts,**runner.runtime_opts.ham_opts,
                **runner.runtime_opts.solver_opts,**runner.pt_opts.opts)
            for k,v in resolve_evaluator_defaults(profile).items():self.assertEqual(opts[k],v)
        runner=VPTRunner(self.system,self.space,initial_states=self.initial)
        self.assertEqual(runner.evaluator_defaults,'new')
        with self.assertRaises(ValueError):VPTRunner.construct(self.system,[[0,0]],evaluator_defaults='typo')

    def test_overrides_do_not_mutate_supplied_options(self):
        h=VPTHamiltonianOptions(cartesian_analytic_deriv_order=0)
        r=VPTRuntimeOptions(operator_contraction_method='tensor')
        s=VPTSolverOptions(use_block_recurrence=False)
        before=(dict(h.opts),dict(r.ham_opts),dict(r.solver_opts),dict(s.opts))
        runner=VPTRunner(self.system,self.space,initial_states=self.initial,
            hamiltonian_options=h,runtime_options=r,solver_options=s)
        self.assertEqual(before,(h.opts,r.ham_opts,r.solver_opts,s.opts))
        self.assertFalse(runner.pt_opts.opts['use_block_recurrence'])
        self.assertEqual(runner.ham_opts.opts['cartesian_analytic_deriv_order'],0)

    def test_native_derivative_profiles(self):
        class Geometry:
            internal_coordinates=None
            coords=np.zeros((2,3))
            def _atomic_masses(self):return np.ones(2)
        for mode,expected in [('new',('fast',-1,'analytical')),('old',('old',0,'old'))]:
            exp=ExpansionTerms(Geometry(),modes=SimpleNamespace(freqs=np.array([.007,.011])),
                undimensionalize=False,expansion_handling_mode=mode)
            self.assertEqual((exp.cartesian_by_internal_derivative_method,
                exp.cartesian_analytic_deriv_order,exp.mixed_derivative_handling_mode.value),expected)

    def test_property_use_retains_closed_paths_through_sixth_order(self):
        for order in range(1,7):
            plan=PropertyCorrectionFilter.from_property(self.targets,order=order)
            self.assertTrue(plan.uses)
        plan=PropertyCorrectionFilter.from_property(self.targets,order=4,target_property='frequencies')
        n=BasisStateSpace(self.basis,[[2,0]])
        self.assertEqual(len(plan(n,BasisStateSpace(self.basis,[[2,2]]),2)),1)
        plan=PropertyCorrectionFilter.from_property(self.targets,order=2)
        self.assertEqual(len(plan(n,BasisStateSpace(self.basis,[[2,2]]),2)),0)

    def test_classic_coefficient_filter_is_order_aware_and_nonmutating(self):
        n=BasisStateSpace(self.basis,[[2,0]])
        ms=BasisStateSpace(self.basis,[[1,1],[2,2],[1,2],[1,0]])
        arrays=[SparseArray.from_data(np.zeros((1,4)))]+[SparseArray.from_data(np.full((1,4),.5)) for _ in range(2)]
        corr=PerturbationTheoryCorrections(n,None,ms,np.zeros((1,3)),arrays)
        keep=lambda n,ms:ms
        old=corr.find_strong_couplings(.3,state_filter=keep)
        new=corr.find_strong_couplings(.3,state_filter=keep,
            correction_filter=PropertyCorrectionFilter.from_property(self.targets,order=2))
        key=int(n.indices[0])
        self.assertEqual(len(old[key][2]),4)
        self.assertEqual(new[key][2].excitations.tolist(),[[1,1]])
        self.assertEqual(len(new[key][1]),2)
        for a in arrays[1:]:np.testing.assert_array_equal(a.asarray(),np.full((1,4),.5))

    def test_per_spec_filter_setting_is_parsed(self):
        for enabled in [True,False,None]:
            spec=DegeneracySpec.from_spec(dict(wfc_threshold=.3,filter_by_property=enabled))
            self.assertIs(spec.filter_by_property,enabled)

    def test_unrelated_evaluators_are_not_adapted(self):
        evaluator=mock.Mock()
        spec=StronglyCoupledDegeneracySpec(evaluator=evaluator)
        self.assertIs(spec.evaluator,evaluator)
        evaluator.assert_not_called()

    def test_custom_operator_rules_and_unknown_support(self):
        rules=PropertyCorrectionFilter.operator_rules(self.basis,{'O':[0,np.ones(2)]},2,False)
        self.assertEqual(rules[0],((),))
        self.assertEqual(set(rules[1]),set(self.basis.selection_rules('x')))
        self.assertIsNone(PropertyCorrectionFilter.operator_rules(self.basis,np.eye(4),2))

    def test_solver_factory_keeps_handlers_and_callback(self):
        runner=VPTRunner(self.system,self.space,initial_states=self.initial)
        captured={}
        runner._ham=SimpleNamespace(get_solver=lambda states,**opts:captured.update(opts))
        runner.get_solver()
        self.assertIs(captured['degeneracy_handlers'],self.space.degeneracy_specs)
        self.assertIs(captured['strong_coupling_correction_filter'],runner.pt_opts.opts['strong_coupling_correction_filter'])
