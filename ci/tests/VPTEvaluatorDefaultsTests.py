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


    def test_vectorized_endpoint_screening_is_chunk_invariant(self):
        import itertools
        from Psience.VPT2.Discovery import AnalyticPropertyDiscovery
        discovery=object.__new__(AnalyticPropertyDiscovery)
        rules={1:[(),(0,),(1,),(-1,)],2:[(),(0,),(2,),(-2,),(1,1),(1,-1),(-1,-1)]}
        discovery._wavefunction_rules=lambda order:rules[order]
        discovery.endpoint_quanta=lambda source,order:set(range(7))
        states=np.asarray(list(itertools.product(range(4),repeat=3)),dtype=int)
        for source in [[0,0,0],[1,1,0],[0,0,2]]:
            source=np.asarray(source)
            for order in [1,2]:
                expected=discovery.allowed_endpoints(source,states,order)
                for size in [1,2,9,64]:
                    actual=np.concatenate([discovery.allowed_endpoints(source,states[i:i+size],order)
                                           for i in range(0,len(states),size)])
                    np.testing.assert_array_equal(actual,expected)

    def test_sparse_discovery_stream_matches_eager_vectorized_filter(self):
        from Psience.VPT2.Discovery import AnalyticPropertyDiscovery
        discovery=object.__new__(AnalyticPropertyDiscovery)
        discovery.eval=SimpleNamespace(freqs=np.array([.75,.25,.5]))
        rules={1:[(1,),(-1,)],2:[(),(2,),(-2,),(1,1),(1,-1),(-1,-1)]}
        discovery._wavefunction_rules=lambda order:rules[order]
        discovery.endpoint_quanta=lambda source,order:{0,1,2,3}
        basis=HarmonicOscillatorProductBasis(3)
        inputs=BasisStateSpace(basis,[[1,0,0],[0,1,1],[0,0,0]])
        sources=inputs.excitations
        sorting=np.argsort(discovery.eval.freqs);inverse=np.argsort(sorting)
        energies=sources@discovery.eval.freqs
        windows=[[max(float(e)-1.,0),float(e)+1.] for e in energies]
        expected={}
        for order in [1,2]:
            maximum=np.max(sources,axis=0)+max(abs(c) for rule in rules[order] for c in rule)
            raw=BasisStateSpace.states_in_windows(discovery.eval.freqs[sorting],windows,
                max_state=maximum[sorting],min_quantas=[0]*3,max_quantas=[3]*3)
            for i,rows in enumerate(raw):
                rows=rows[:,inverse]
                expected[i,order]=rows[discovery.allowed_endpoints(sources[i],rows,order)]
        native_windows=BasisStateSpace.states_in_windows
        def streamed_windows(*args,**options):
            self.assertTrue(options.get('return_iterator'))
            self.assertFalse(options.get('sparse'))
            self.assertNotIn('filter',options)
            iterator=native_windows(*args,**options)
            def observed():
                for tags,block in iterator:
                    self.assertLessEqual(tags.nbytes+block.nbytes,512)
                    previous=screen.call_count
                    yield tags,block
                    self.assertGreater(screen.call_count,previous,
                        'a work block was retained without running the vectorized filter')
            return observed()
        with mock.patch.object(discovery,'allowed_endpoints',wraps=discovery.allowed_endpoints) as screen, \
             mock.patch.object(BasisStateSpace,'states_in_windows',side_effect=streamed_windows):
            result=discovery.candidate_spaces(inputs,1.,(1,2),block_size=8,max_bytes=512)
        self.assertTrue(screen.call_args_list)
        self.assertTrue(all(call.args[1].ndim==2 for call in screen.call_args_list))
        self.assertTrue(any(len(call.args[1])>1 for call in screen.call_args_list))
        for i,orders in enumerate(result):
            for order,space in orders.items():
                self.assertTrue(space.has_sparse_excitations)
                self.assertIsNone(space._excitations)
                np.testing.assert_array_equal(space.sparse_excitations.asarray(),expected[i,order])
                np.testing.assert_array_equal(space.indices,basis.ravel_state_inds(expected[i,order]))

    def test_vectorized_discovery_preserves_large_sparse_occupations(self):
        from Psience.VPT2.Discovery import AnalyticPropertyDiscovery
        basis=HarmonicOscillatorProductBasis(1)
        inputs=BasisStateSpace(basis,SparseArray.from_data(np.array([[130]],dtype=np.int16),layout='csr'))
        discovery=object.__new__(AnalyticPropertyDiscovery)
        discovery.eval=SimpleNamespace(freqs=np.array([1.]))
        discovery._wavefunction_rules=lambda order:[(0,)]
        discovery.endpoint_quanta=lambda source,order:{130}
        result=discovery.candidate_spaces(inputs,1.,(1,),max_bytes=512)[0][1]
        self.assertTrue(result.has_sparse_excitations)
        np.testing.assert_array_equal(result.sparse_excitations.asarray(),[[130]])


    def test_window_rule_prefix_traversal_all_storage_and_chunk_modes(self):
        from Psience.BasisReps import BasisStateSpaceFilter
        import itertools
        basis=HarmonicOscillatorProductBasis(4)
        sources=BasisStateSpace(basis,[[2,0,1,0],[0,1,0,1]])
        rules=[[],[0],[1],[-1],[2,-1],[1,1,-1],[1,-1],[2,-2],[1,1],[-1,-1]]
        windows=np.array([[-2.,3.],[-1.,2.]])
        frequencies=np.array([.3,-.2,.1,.4])
        floor=np.array([1,0,0,0]);maximum=np.array([3,2,2,2])
        filt=BasisStateSpaceFilter(sources,[(None,rules)],None)
        def expected(source,window):
            signatures={tuple(sorted(x for x in rule if x)) for rule in rules}
            rows=[]
            for row in itertools.product(*(range(int(low),int(high)+1) for low,high in zip(floor,maximum))):
                if tuple(sorted(int(v) for v in np.asarray(row)-source if v)) not in signatures:continue
                energy=float(np.dot(frequencies,floor))
                for mode in range(4):
                    for _ in range(row[mode]-floor[mode]):energy+=float(frequencies[mode])
                if window[0]<=energy<=window[1] and sum(row)<=5:rows.append(row)
            return np.asarray(sorted(rows,key=lambda row:(sum(row),tuple(-x for x in row))),dtype=int).reshape(-1,4)
        wanted=[expected(source,window) for source,window in zip(sources.excitations,windows)]
        for sparse,traversal,eager in itertools.product([False,True],repeat=3):
            for chunk in [1,7]:
                with self.subTest(sparse=sparse,traversal=traversal,eager=eager,chunk=chunk):
                    options=dict(initial_state=floor,max_state=maximum,max_quantas=5,
                        sparse=sparse,sparse_traversal=traversal,filter=filt,block_size=chunk,max_bytes=2048)
                    if eager:
                        result=BasisStateSpace.states_in_windows(frequencies,windows,**options)
                        actual=[r.asarray() if sparse else r for r in result]
                    else:
                        rows=[[] for _ in windows]
                        for tags,block in BasisStateSpace.states_in_windows(
                                frequencies,windows,return_iterator=True,**options):
                            value=block.asarray() if sparse else block
                            self.assertLessEqual(len(tags),chunk)
                            for window in np.unique(tags):rows[int(window)].append(value[tags==window])
                        actual=[np.concatenate(r) if r else np.empty((0,4),dtype=int) for r in rows]
                    for left,right in zip(actual,wanted):np.testing.assert_array_equal(left,right)

    def test_window_rule_prefix_exact_boundaries_and_resume(self):
        from Psience.BasisReps import BasisStateSpaceFilter
        basis=HarmonicOscillatorProductBasis(3)
        source=BasisStateSpace(basis,[[1,0,0]])
        frequencies=np.array([.1,.2,.3]);floor=np.array([0,0,0])
        target=np.array([2,1,0])
        energy=0.
        for mode,occupation in enumerate(target):
            for _ in range(int(occupation)):energy+=float(frequencies[mode])
        rules=[[1,1],[1],[2,-1],[1,-1],[1,1,-1]]
        filt=BasisStateSpaceFilter(source,[(None,rules)],None)
        windows=[[energy,energy],[np.nextafter(energy,np.inf),energy+1.]]
        actual=BasisStateSpace.states_in_windows(frequencies,windows,filter=filt,max_quantas=3)
        reference=BasisStateSpace.states_in_windows(frequencies,windows,max_quantas=3)
        signatures={tuple(sorted(rule)) for rule in rules}
        for a,r in zip(actual,reference):
            selected=[row for row in r if tuple(sorted(int(x) for x in row-source.excitations[0] if x)) in signatures]
            np.testing.assert_array_equal(a,np.asarray(selected,dtype=int).reshape(-1,3))
        statistics={}
        stream=BasisStateSpace.states_in_windows(frequencies,[[0.,1.]],filter=filt,max_quantas=3,
            return_iterator=True,block_size=1,statistics=statistics)
        tag,first=next(stream);held=first.copy()
        self.assertEqual(statistics['emitted_rows'],1)
        self.assertEqual(statistics['maximum_output_rows'],1)
        next(stream)
        self.assertEqual(statistics['emitted_rows'],2)
        self.assertEqual(statistics['scalar_filter_calls'],0)
        np.testing.assert_array_equal(first,held)
        stream.close()

    def test_window_rule_prefix_prefilter_priority_and_fixed_floor(self):
        from Psience.BasisReps import BasisStateSpaceFilter
        basis=HarmonicOscillatorProductBasis(2)
        sources=BasisStateSpace(basis,[[1,0],[0,1]])
        first=BasisStateSpace(basis,[[1,0]])
        filt=BasisStateSpaceFilter(sources,[(first,[[1]]),(None,[[-1]])],None)
        result=BasisStateSpace.states_in_windows([.1,.2],[[0.,2.]]*2,
            filter=filt,max_quantas=3)
        np.testing.assert_array_equal(result[0],[[2,0],[1,1]])
        np.testing.assert_array_equal(result[1],[[0,0]])
        filt=BasisStateSpaceFilter(first,[(None,[[],[1]])],None)
        result=BasisStateSpace.states_in_windows([.1,.2],[[0.,2.]],filter=filt,
            initial_state=[1,0],max_state=[0,2],fixed_modes=[0],max_quantas=3)[0]
        np.testing.assert_array_equal(result,[[1,0],[1,1]])
        with self.assertRaises(ValueError):
            BasisStateSpace.states_in_windows([.1,.2],[[0.,1.]],
                filter=BasisStateSpaceFilter(first,None,None))

    def test_property_generation_filter_matches_vectorized_candidates(self):
        from Psience.VPT2.Discovery import AnalyticPropertyDiscovery
        basis=HarmonicOscillatorProductBasis(4)
        inputs=BasisStateSpace(basis,[[0,0,0,0],[1,0,0,0],[0,0,1,1],[0,2,0,0]])
        discovery=object.__new__(AnalyticPropertyDiscovery)
        discovery.eval=SimpleNamespace(freqs=np.array([.2,.1,.3,.25]))
        rules={1:[[1],[-1],[3],[2,-1],[1,1,-1]],
               2:[[],[0],[2],[-2],[1,1],[1,-1],[2,-2],[1,1,-1,-1]]}
        discovery._wavefunction_rules=lambda order:rules[order]
        discovery.endpoint_quanta=lambda row,order:{int(np.sum(row))+sum(r) for r in rules[order]
            if int(np.sum(row))+sum(r)>=0}
        wanted=discovery.candidate_spaces(inputs,.5,(1,2),block_size=7,max_bytes=2048)
        with mock.patch.object(discovery,'allowed_endpoints',
                               side_effect=AssertionError('late endpoint filter called')):
            actual=discovery.candidate_spaces(inputs,.5,(1,2),block_size=1,max_bytes=512,
                                               filter_at_generation=True)
        for a,b in zip(actual,wanted):
            for order in a:
                np.testing.assert_array_equal(a[order].excitations,b[order].excitations)
                self.assertTrue(a[order].has_sparse_excitations)
