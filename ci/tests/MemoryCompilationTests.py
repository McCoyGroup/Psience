"""Native regressions for RAM compilation, explicit opt-out and timing accounting."""
import contextlib,itertools,io,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock
import numpy as np
from McUtils.Parallelizers import MultiprocessingParallelizer
from McUtils.Scaffolding import Logger
from Psience.VPT2 import AnalyticVPTRunner
from Psience.VPT2 import CompiledExpressions as ce
from Psience.VPT2.Projected import ProjectedCorrectionEvaluator
from Psience.VPT2.Runner import _timed_vpt_step

class MemoryCompilationTests(unittest.TestCase):
    def runner(self,**opts):
        w=np.array([.0201,.0217,.0313]);v3=np.zeros((3,)*3);v4=np.zeros((3,)*4)
        for indices,value in [((0,0,0),.0003),((0,1,1),.00008),((0,0,1),.0001)]:
            for p in set(itertools.permutations(indices)):v3[p]=value
        v4[0,0,0,0]=.00001
        dipole=[[0.,np.array([.2,.1,.3]),np.diag([.02,-.01,.01]),np.zeros((3,)*3)]]
        dipole+=[[0.,np.zeros(3),np.zeros((3,3)),np.zeros((3,)*3)] for _ in range(2)]
        return AnalyticVPTRunner([[np.diag(w)/2,np.diag(w)/2],[v3,0],[v4,0,0,0,0]],
            order=2,freqs=w,internals=False,dipole_expansion=dipole,logger=False,**opts)
    def energy(self,runner):return np.asarray(runner.eval.get_energy_corrections([[0,0,0],[1,0,0]],order=2))
    def moments(self,runner,block=None):
        states=[[0,0,0],[1,0,0],[0,1,0]]
        space=runner.prep_states(states,degeneracy_specs={'groups':[block or [[1,0,0],[0,1,0]]]})
        return np.asarray(runner.get_transition_moment_corrections(space,order=2,state_specific_degeneracies=True))
    @contextlib.contextmanager
    def no_expression_writes(self):
        with mock.patch.object(ce,'_atomic_json',side_effect=AssertionError('expression JSON write')),\
             mock.patch.object(ce,'_atomic_pickle',side_effect=AssertionError('expression pickle write')),\
             mock.patch.object(ce,'_CompiledEntries',side_effect=AssertionError('disk expression store')):
            yield
    def test_default_compiles_in_memory_without_expression_writes(self):
        r=self.runner()
        with self.no_expression_writes():result=self.energy(r)
        self.assertTrue(np.isfinite(result).all());self.assertIsNone(r.eval.expression_store.folder)
        self.assertGreater(r.eval.compilation_statistics()['compilations'],0)
    def test_memory_disk_and_disabled_energy_match(self):
        reference=self.energy(self.runner(compile=False))
        with tempfile.TemporaryDirectory() as folder:
            for r in [self.runner(),self.runner(expressions_folder=folder)]:
                np.testing.assert_allclose(self.energy(r),reference,rtol=2e-12,atol=2e-14)
    def test_memory_and_disabled_projected_transition_moments_match(self):
        reference=self.moments(self.runner(compile=False));r=self.runner();contexts=[]
        original=ProjectedCorrectionEvaluator.operators
        def observed(*args,**kwargs):contexts.append(ce.current_compilation());return original(*args,**kwargs)
        with self.no_expression_writes(),mock.patch.object(ProjectedCorrectionEvaluator,'operators',observed):
            result=self.moments(r)
        np.testing.assert_allclose(result,reference,rtol=2e-11,atol=2e-13)
        self.assertTrue(contexts);self.assertTrue(all(c is not None for c in contexts))
    def test_compile_false_suppresses_even_an_explicit_disk_folder(self):
        with tempfile.TemporaryDirectory() as parent:
            folder=Path(parent)/'must-not-exist';r=self.runner(compile=False,expressions_folder=folder)
            with self.no_expression_writes():self.energy(r)
            self.assertFalse(folder.exists());self.assertEqual(r.eval.compilation_statistics()['requests'],0)
    def test_compile_false_suppresses_an_outer_context_and_restores_it(self):
        store=ce.AnalyticExpressionStore()
        with ce.compiled_evaluation_context(store):
            self.energy(self.runner(compile=False))
            self.assertIs(ce.current_compilation()[0],store)
        self.assertEqual(store.statistics()['requests'],0);self.assertIsNone(ce.current_compilation())
    def test_warm_memory_cache_reuses_compiled_entries(self):
        r=self.runner();a=self.energy(r);before=r.eval.compilation_statistics();b=self.energy(r);after=r.eval.compilation_statistics()
        np.testing.assert_array_equal(a,b)
        self.assertEqual(after['compilations'],before['compilations'])
        self.assertGreater(after['cache_hits'],before['cache_hits'])
    def test_disk_cache_remains_loadable_and_reuses_parts(self):
        with tempfile.TemporaryDirectory() as folder:
            a=self.energy(self.runner(expressions_folder=folder));parts={str(p):p.stat().st_mtime_ns for p in Path(folder).rglob('part-*.pkl')}
            r=self.runner(expressions_folder=folder);b=self.energy(r)
            np.testing.assert_array_equal(a,b);self.assertTrue(parts)
            self.assertEqual(parts,{str(p):p.stat().st_mtime_ns for p in Path(folder).rglob('part-*.pkl')})
            self.assertEqual(r.eval.compilation_statistics()['compilations'],0)
    def test_memory_cache_entry_and_payload_limits_evict_only_reuse(self):
        store=ce.AnalyticExpressionStore(memory_cache_size=1,memory_cache_bytes=100)
        def entry(n,size):return ce.CompiledAnalyticExpression(None,dict(signature=str(n)*64,parts=[dict(bytes=size)]),parts={})
        a=entry(1,40);b=entry(2,60);store._remember(a);store._remember(b)
        self.assertIsNone(store._get_memory('1'*64));self.assertIs(store._get_memory('2'*64),b)
        self.assertEqual(a.manifest['parts'][0]['bytes'],40)
        large=entry(3,200);self.assertIs(store._remember(large),large);self.assertIsNone(store._get_memory('3'*64))
        store.clear();self.assertEqual(store.statistics()['memory_payload_bytes'],0)
    def test_run_level_override_restores_policy_on_return_and_exception(self):
        r=self.runner()
        def inner(*args,**kwargs):self.assertFalse(r.eval.compile_expressions);return 17
        with mock.patch.object(r,'_run_VPT',side_effect=inner):self.assertEqual(r.run_VPT([],compile=False),17)
        self.assertTrue(r.eval.compile_expressions)
        with mock.patch.object(r,'_run_VPT',side_effect=RuntimeError('fixture')):
            with self.assertRaisesRegex(RuntimeError,'fixture'):r.run_VPT([],compile=False)
        self.assertTrue(r.eval.compile_expressions)
    def test_stage_logging_reports_new_and_reused_compilation_separately(self):
        stream=io.StringIO();logger=Logger(stream);r=self.runner()
        with _timed_vpt_step(logger,'Cold',r.eval):self.energy(r)
        before=r.eval.compilation_statistics()
        with _timed_vpt_step(logger,'Warm',r.eval):self.energy(r)
        after=r.eval.compilation_statistics();log=stream.getvalue()
        self.assertIn('Cold: expression compilation',log);self.assertIn('Warm: expression compilation 0.000s (0 new,',log)
        self.assertIn('cache lookup',log);self.assertIn('numerical evaluation',log)
        for key in ['compilation_seconds','cache_lookup_seconds','preparation_seconds','evaluation_seconds']:
            self.assertGreaterEqual(after[key],before[key]);self.assertTrue(np.isfinite(after[key]))
    def test_failed_compilation_is_timed_without_creating_an_entry(self):
        store=ce.AnalyticExpressionStore()
        with self.assertRaisesRegex(RuntimeError,'build fixture'):
            store.compile_factory(lambda:(_ for _ in ()).throw(RuntimeError('build fixture')),context={'fixture':1})
        stats=store.statistics();self.assertEqual(stats['failures'],1);self.assertEqual(stats['compilations'],0)
        self.assertGreater(stats['compilation_seconds'],0.);self.assertEqual(stats['memory_entries'],0)
    def test_memory_parallel_execution_matches_serial(self):
        expected=self.energy(self.runner())
        with MultiprocessingParallelizer(processes=2,stall_timeout=60) as backend:
            r=self.runner(parallelizer=backend)
            with self.no_expression_writes():actual=self.energy(r)
        np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=2e-14)
    def test_explicit_frequency_compilation_uses_memory_without_a_folder(self):
        r=self.runner();reference=self.energy(self.runner(compile=False))
        with self.no_expression_writes():
            compiled=r.eval.precompile_frequency_corrections([0,0,0],[1,0,0],order=2)
            values=r.eval.get_frequency_corrections([0,0,0],[1,0,0],order=2)
        self.assertTrue(compiled)
        self.assertAlmostEqual(float(np.sum(values)),float(np.sum(reference[:,1]-reference[:,0])),places=12)
        self.assertIsNone(r.eval.expression_store.folder)
    def test_cache_reconfiguration_resets_timing_epoch_without_negative_counts(self):
        r=self.runner();reference=self.energy(r)
        with tempfile.TemporaryDirectory() as parent:
            folder=Path(parent)/'new-cache';r.eval.expressions_folder=folder
            self.assertEqual(r.eval.compilation_statistics()['compilations'],0)
            self.assertFalse(folder.exists())
            np.testing.assert_array_equal(self.energy(r),reference)
        r.eval.expressions_folder=None;r.eval.compilation_options['cache_bytes']=64*1024**2
        self.assertEqual(r.eval.compilation_statistics()['compilations'],0)
        with self.no_expression_writes():np.testing.assert_array_equal(self.energy(r),reference)
        self.assertEqual(r.eval.expression_store.memory_cache_bytes,64*1024**2)
    def test_new_resonance_block_rebinds_equations_without_reusing_values(self):
        r=self.runner();self.moments(r);before=r.eval.compilation_statistics()['cache_hits']
        block=[[1,0,0],[0,0,1]]
        with self.no_expression_writes():actual=self.moments(r,block)
        expected=self.moments(self.runner(compile=False),block)
        np.testing.assert_allclose(actual,expected,rtol=2e-11,atol=2e-13)
        self.assertGreater(r.eval.compilation_statistics()['cache_hits'],before)
    def test_hamiltonian_clone_preserves_policy_and_derives_new_harmonic_frequencies(self):
        r=self.runner(compile=False);new_frequencies=np.array([.04,.06,.08])
        expansions=[[np.diag(new_frequencies)/2,np.diag(new_frequencies)/2],[0,0],[0,0,0,0,0]]
        with mock.patch.object(r.eval,'_prep_ham_corrs',return_value=([],[])),\
             mock.patch.object(r.eval,'_prep_operator_expansion',return_value=expansions),\
             mock.patch.object(r.eval.solver,'modify_hamiltonian',return_value=r.eval.solver):
            changed=r.eval.modify_hamiltonian(object())
        np.testing.assert_array_equal(changed.freqs,new_frequencies)
        self.assertFalse(changed.compile_expressions)



"""Proposed additions to Psience/ci/tests/MemoryCompilationTests.py."""
import contextlib, gc, itertools, unittest, weakref, time
from unittest import mock
import numpy as np
from Psience.VPT2 import AnalyticVPTRunner
from Psience.VPT2.Analytic import PolyPath, PerturbationTheoryExpressionEvaluator
from Psience.VPT2.BoundedStream import StreamBudget, retained_size
from Psience.VPT2.Projected import AnalyticBlockExpansion, ProjectedCorrectionEvaluator

class CompilationOverheadTests(unittest.TestCase):
    def test_weak_comparison_cache_preserves_equal_unequal_and_lifetimes(self):
        a=PolyPath.from_coeffs([np.array([1.,2.])],steps=0)
        b=PolyPath.from_coeffs([np.array([1.,2.])],steps=0)
        c=PolyPath.from_coeffs([np.array([1.,3.])],steps=0)
        with PolyPath.cached_comparisons(2):
            self.assertEqual(a,b);self.assertEqual(a,b);self.assertNotEqual(a,c)
            pairs,limit=PolyPath._comparison_context.get()
            self.assertLessEqual(len(pairs),limit)
            ref=weakref.ref(b);del b;gc.collect();self.assertIsNone(ref())
            with PolyPath.cached_comparisons(0):self.assertEqual(a,a)
        self.assertIsNone(PolyPath._comparison_context.get())
        with self.assertRaisesRegex(ValueError,'nonnegative'):
            with PolyPath.cached_comparisons(-1):pass

    def test_bounded_accounting_preserves_flush_boundaries(self):
        decisions=[]
        for bounded in (False,True):
            budget=StreamBudget(dict(max_terms=5,max_bytes=800,check_interval=2,
                                     bounded_byte_accounting=bounded))
            pending=[];events=[]
            for i in range(20):
                pending.append(np.arange(i+7,dtype=float))
                reason=budget.incoming(pending)
                if reason:
                    events.append((i,reason));pending.clear();budget.flushed(reason)
            decisions.append(events)
        self.assertEqual(*decisions)

    def test_bounded_walk_is_exact_below_limit_and_a_lower_bound_above(self):
        base=np.arange(200);shared={'base':base,'view':base[::2]}
        graph=[shared,shared];graph.append(graph)
        size=retained_size(graph)
        self.assertEqual(retained_size(graph,stop_after=size+1),size)
        cut=retained_size(graph,stop_after=100)
        self.assertGreater(cut,100);self.assertLessEqual(cut,size)

    def test_degree_hint_matches_native_evaluation_and_bypass_is_optional(self):
        from Psience.VPT2 import CompiledExpressions as ce
        from Psience.VPT2.Analytic import PTTensorCoeffProductSum
        poly=PTTensorCoeffProductSum({((1,0,0,0,0),):
            PolyPath.from_coeffs([np.array([1.,.2])],steps=0)})
        manifest=dict(signature='degree-regression',projection=None,options={},num_fixed=1,
                      change=[1],baseline_shift=None,
                      parts=[dict(file='fixture',rank=1,energy_rows=1)])
        expression=ce.CompiledAnalyticExpression(None,manifest,parts={'fixture':poly})
        args=([(np.zeros(2,dtype=int),np.array([[0,1]]))],[[[0,0],[np.ones((2,)*3),0]]],np.array([.02,.03]))
        results=[]
        for cache_degree in (False,True):
            evaluation=expression.prepare_evaluation(*args,resume=False,
                cache_polynomial_order=cache_degree,path_equality_cache_size=0)
            results.append(evaluation.run())
        for a,b in zip(*results):np.testing.assert_array_equal(a,b)

def intermediate_fixture():
    """Two modes; fully symmetric V, valid p/Q-separated G, nonzero G3."""
    rng=np.random.default_rng(20261107);n=2;freqs=np.array([.0201,.0108]);exp=[[np.diag(freqs)/2]*2]
    for order in (1,2,3):
        rank=order+2
        raw=rng.normal(size=(n,)*rank)
        v=sum(raw.transpose(p) for p in itertools.permutations(range(rank)))
        import math
        v*=2e-5/(10**(order-1)*math.factorial(rank))
        raw=rng.normal(size=(n,)*rank)
        g=sum(raw.transpose(p+(order,order+1))+raw.transpose(p+(order+1,order))
              for p in itertools.permutations(range(order)))
        g*=1e-5/(10**(order-1)*2*math.factorial(order))
        exp.append([v,g]+([np.zeros((n,)*(order-2)),0,0] if order>=2 else []))
    return exp,freqs,[[1,0],[0,2]]

class IntermediateBlockTests(unittest.TestCase):
    observations={}
    def test_direct_recurrence_defers_metric_and_recovers_canonical_series(self):
        exp,freqs,states=intermediate_fixture()
        runner=AnalyticVPTRunner(exp,order=3,freqs=freqs,internals=True,logger=False,
                                 compilation_options={'mode_dimension':2})
        normal=AnalyticBlockExpansion(runner.eval,states,zero_cutoff=1e-24)
        direct=AnalyticBlockExpansion(runner.eval,states,zero_cutoff=1e-24,
                                     intermediate_normalization=True)
        from Psience.VPT2 import CompiledExpressions as ce
        with ce.compiled_evaluation_context(runner.eval.expression_store):
            tick=time.perf_counter();normal.extend(3);normal_seconds=time.perf_counter()-tick
            tick=time.perf_counter();direct.extend(3);direct_seconds=time.perf_counter()-tick
        self.assertEqual(len(direct._metric),1)
        for n in range(1,4):np.testing.assert_array_equal(direct.N[n],np.zeros((2,2)))
        with ce.compiled_evaluation_context(runner.eval.expression_store):
            tick=time.perf_counter();transformed=direct.orthonormal_hamiltonian(3);metric_seconds=time.perf_counter()-tick
        comparisons=[]
        for n,(a,b) in enumerate(zip(transformed,normal.F)):
            den=np.abs(b);mask=den>0;delta=np.abs(a-b)
            comparisons.append(dict(order=n,max_absolute_hartree=float(delta.max()),
                max_percent=float((100*delta[mask]/den[mask]).max()) if mask.any() else None,
                denominator='Absolute corresponding nonzero canonical reference entry'))
        type(self).observations=dict(normal_seconds=normal_seconds,direct_seconds=direct_seconds,
            deferred_metric_seconds=metric_seconds,comparisons=comparisons,
            bloch=[a.tolist() for a in direct.F],canonical=[a.tolist() for a in normal.F],
            transformed=[a.tolist() for a in transformed],metric=[a.tolist() for a in direct._metric],
            timing_scope='Small two-mode fixture; shared symbolic/runtime caches, not a full-molecule speed benchmark')
        for a,b in zip(transformed,normal.F):np.testing.assert_allclose(a,b,rtol=3e-10,atol=3e-13)
        before=len(direct._elements)
        with mock.patch.object(direct,'inner',side_effect=AssertionError('repeated overlap contraction')):
            again=direct.orthonormal_hamiltonian(3)
            converted_block=direct.orthonormal_block(3)
        self.assertEqual(before,len(direct._elements))
        for a,b in zip(again,transformed):np.testing.assert_array_equal(a,b)
        for a,b in zip(converted_block.F,normal.F):np.testing.assert_allclose(a,b,rtol=3e-10,atol=3e-13)
        for a,b in zip(converted_block.N,normal.N):np.testing.assert_allclose(a,b,rtol=3e-10,atol=3e-13)
        # Independently recover the intermediate third-order Q columns from
        # the existing normalized series: Omega_3=Y_3-Y_1 N_2-P N_3.
        # All sampled final states lie in Q, so the last term is exactly zero.
        qstates=[[0,0],[0,1],[2,0],[1,1]];wave_comparisons=[]
        with ce.compiled_evaluation_context(runner.eval.expression_store):
            for a in range(len(states)):
                for final in qstates:
                    reference=sum(normal.element(term,origin,final)
                                  for origin,term in normal.Y[3][a].items())
                    canonical_reference=reference
                    reference-=sum(normal.N[2][b,a]*sum(normal.element(term,origin,final)
                        for origin,term in normal.Y[1][b].items()) for b in range(len(states)))
                    actual=sum(direct.element(term,origin,final)
                               for origin,term in direct.Y[3][a].items())
                    converted=sum(converted_block.element(term,origin,final)
                                  for origin,term in converted_block.Y[3][a].items())
                    self.assertAlmostEqual(converted,canonical_reference,
                        delta=3e-12+3e-10*abs(canonical_reference))
                    self.assertAlmostEqual(actual,reference,delta=3e-12+3e-10*abs(reference))
                    wave_comparisons.append(dict(column=a,final=final,actual=float(actual),reference=float(reference),
                        absolute_difference=float(abs(actual-reference)),
                        percent_difference=float(100*abs(actual-reference)/abs(reference)) if reference else None,
                        converted_canonical=float(converted),canonical_reference=float(canonical_reference),
                        canonical_absolute_difference=float(abs(converted-canonical_reference)),
                        canonical_percent_difference=float(100*abs(converted-canonical_reference)/abs(canonical_reference)) if canonical_reference else None,
                        denominator='Absolute normalized-series conversion reference amplitude'))
        type(self).observations['wavefunction_comparisons']=wave_comparisons
        with self.assertRaisesRegex(ValueError,'normalized observables'):
            direct.operator_correction(direct,0,0,runner.eval.solver.hamiltonian_expansion,2)

    def test_intermediate_wavefunctions_do_not_build_normalized_block(self):
        exp,freqs,states=intermediate_fixture()
        runner=AnalyticVPTRunner(exp,order=2,freqs=freqs,internals=True,
                                 intermediate_normalization=True,logger=False,
                                 compilation_options={'mode_dimension':2})
        projected=ProjectedCorrectionEvaluator(runner.eval,blocks=[states],zero_cutoff=1e-24)
        with mock.patch.object(AnalyticBlockExpansion,'overlap_metric',side_effect=AssertionError('metric')):
            projected.wavefunctions([[states[0],[states[1]]]],order=2)
        self.assertTrue(projected.cache)
        self.assertTrue(all(block.intermediate_normalization for block in projected.cache.values()))

    def test_explicit_bloch_api_keeps_default_canonical_path_separate(self):
        exp,freqs,states=intermediate_fixture()
        runner=AnalyticVPTRunner(exp,order=2,freqs=freqs,internals=True,logger=False,
                                 compilation_options={'mode_dimension':2})
        pairs=[[s,states] for s in states]
        raw=runner.eval.get_reexpressed_hamiltonian(pairs,order=2,degenerate_blocks=[states],
             include_diagonal=True,normalization='intermediate',zero_cutoff=1e-24)
        self.assertTrue(raw.corrections)
        reused=ProjectedCorrectionEvaluator(runner.eval,blocks=[states],zero_cutoff=1e-24)
        direct=reused.cache[(tuple(sorted(tuple(s) for s in states)),'intermediate')]
        canonical,_=reused.block(states[0],2)
        self.assertIs(canonical,direct._canonical_expansion)
        self.assertFalse(canonical.intermediate_normalization)
        self.assertIs(reused.block(states[0],2)[0],canonical)
        with self.assertRaisesRegex(ValueError,'contribution triples'):
            projected=ProjectedCorrectionEvaluator(runner.eval,blocks=[states])
            projected.hamiltonian(pairs,2,terms=[(0,1,1)],normalization='intermediate')

class FreeIndexExecutorTests(unittest.TestCase):
    observations = {}

    def fixture(self):
        from Psience.VPT2.Analytic import (PTTensorCoeffProductSum as Coeff,
            ProductPTPolynomial, PTEnergyChangeProductSum,
            PerturbationTheoryExpressionEvaluator as Expression)
        a=Coeff.canonical_key(((1,0,0,0,1),))
        b=Coeff.canonical_key(((1,0,0,1,1),))
        p=ProductPTPolynomial([np.array([1.,.3]),np.array([2.,.5])],steps=0)
        value=PTEnergyChangeProductSum({((0,1,-2),):p,():p.scale(.01)})
        poly=Coeff({a:value,b:value.permute((1,0))},canonicalize=False)
        wrapper=Expression(None,poly,change=(),logger=False)
        v=np.empty((4,)*3)
        for index in np.ndindex(v.shape):
            v[index]=1+sum((j+1)*i for j,i in enumerate(sorted(index)))
        coeffs=[[[0,0],[v,0]]]
        states=[(np.array([0,1,2,3]),np.array([[0,1,2,3],[3,0,2,1]])),
                (np.array([2,0,3,1]),np.array([[1,3,0,2]]))]
        return wrapper,states,coeffs,np.array([.013,.021,.034,.049])

    def compare(self,actual,reference):
        delta=np.abs(np.asarray(actual)-np.asarray(reference))
        den=np.abs(reference);mask=den>0
        np.testing.assert_allclose(actual,reference,rtol=3e-12,atol=3e-12)
        return dict(max_absolute=float(delta.max()),
            max_percent=float((100*delta[mask]/den[mask]).max()) if mask.any() else None,
            denominator='Absolute corresponding nonzero unreduced reference entry')

    def test_shared_plans_preserve_orientations_and_nonuniform_states(self):
        from Psience.VPT2.IndexedEvaluator import IndexedBlockEvaluator
        wrapper,states,coeffs,freqs=self.fixture()
        store=ce.AnalyticExpressionStore()
        baseline=store.compile(wrapper,context={'test':'binding'},mode_dimension=4,
            free_index_reduction=False)
        candidate=store.compile(wrapper,context={'test':'binding'},mode_dimension=4,
            free_index_reduction=True)
        self.assertNotEqual(baseline.manifest['signature'],candidate.manifest['signature'])
        part=next(iter(candidate.parts.values()))
        metadata=part._free_index_plan_data
        self.assertEqual(metadata['stats']['canonical_plans'],1)
        self.assertEqual(metadata['stats']['shared_plans'],1)
        plan=IndexedBlockEvaluator.compile_plan(part,reduced=True)
        keys=list(plan.subexpressions)
        self.assertIs(plan.subexpressions[keys[0]],plan.subexpressions[keys[1]])
        expected=baseline.evaluate(states,coeffs,freqs,resume=False,batch_size=2)
        seen=[];original=IndexedBlockEvaluator.finish
        def finish(executor):
            record=original(executor);seen.append(record);return record
        with mock.patch.object(IndexedBlockEvaluator,'finish',finish):
            actual=candidate.evaluate(states,coeffs,freqs,resume=False,batch_size=2)
        checks=[self.compare(a,b) for a,b in zip(actual,expected)]
        self.assertTrue(any(s['free_index_bound_blocks']>0 for s in seen))
        disabled=candidate.evaluate(states,coeffs,freqs,resume=False,batch_size=3,
            free_index_reduction=False)
        for a,b in zip(disabled,expected):self.compare(a,b)
        type(self).observations['binding']=checks

    def test_disk_roundtrip_resume_and_flag_cache_namespaces(self):
        wrapper,states,coeffs,freqs=self.fixture()
        with tempfile.TemporaryDirectory() as folder:
            store=ce.AnalyticExpressionStore(folder)
            a=store.compile(wrapper,context={'test':'disk'},mode_dimension=4,
                free_index_reduction=False)
            b=store.compile(wrapper,context={'test':'disk'},mode_dimension=4,
                free_index_reduction=True)
            expected=a.evaluate(states,coeffs,freqs,batch_size=2,task_size=2)
            first=b.evaluate(states,coeffs,freqs,batch_size=2,task_size=2)
            loaded=ce.CompiledAnalyticExpression.load(b.folder)
            second=loaded.evaluate(states,coeffs,freqs,batch_size=2,task_size=2)
            for x,y,z in zip(expected,first,second):
                self.compare(y,x);np.testing.assert_array_equal(y,z)
        with self.assertRaisesRegex(TypeError,'boolean'):
            ce.AnalyticExpressionStore().compile(wrapper,context={},free_index_reduction='yes')

    def test_external_predicates_keep_native_plan(self):
        from Psience.VPT2.IndexedEvaluator import IndexedBlockEvaluator
        wrapper,states,coeffs,freqs=self.fixture()
        compiled=ce.AnalyticExpressionStore().compile(wrapper,context={},mode_dimension=4,
            free_index_reduction=True)
        part=next(iter(compiled.parts.values()))
        token=ce._current_free_index_reduction.set(True)
        try:
            executor=IndexedBlockEvaluator(PerturbationTheoryExpressionEvaluator,part)
        finally:ce._current_free_index_reduction.reset(token)
        self.assertIsNotNone(executor.plan.index_bindings)
        # Stop after native-plan selection, before coefficient gathering. This
        # tests the fallback boundary without inventing an external predicate.
        with mock.patch.object(executor,'get_coefficient_plan',side_effect=RuntimeError('stop')):
            with self.assertRaisesRegex(RuntimeError,'stop'):
                executor.eval_perm(part,(),None,(0,1),states,None,None,freqs,(),(),0,
                    object(),[False,False,np.logical_and],1e-24,None,None,None,None,None,None,None,
                    False,None,False)
        self.assertIsNone(executor.plan.index_bindings)
        self.assertEqual(executor.stats.free_index_fallback_calls,1)

    def test_native_projected_internal_blocks_through_third_order(self):
        matrices=[]
        for enabled in (False,True):
            exp,freqs,states=intermediate_fixture()
            runner=AnalyticVPTRunner(exp,order=3,freqs=freqs,internals=True,logger=False,
                compilation_options={'mode_dimension':2,'free_index_reduction':enabled})
            block=AnalyticBlockExpansion(runner.eval,states,zero_cutoff=1e-24)
            with ce.compiled_evaluation_context(runner.eval.expression_store,
                    compilation_options={'mode_dimension':2,'free_index_reduction':enabled}):
                block.extend(3)
            matrices.append(np.asarray(block.F))
        type(self).observations['projected_internal']=self.compare(matrices[1],matrices[0])

    def test_coalescing_switch_remains_independent(self):
        arrays=[]
        for coalescing in (False,True):
            for reduction in (False,True):
                runner=MemoryCompilationTests().runner(compilation_options={
                    'coalesce_combine_nodes':coalescing,'free_index_reduction':reduction})
                arrays.append(MemoryCompilationTests().energy(runner))
        for actual in arrays[1:]:self.compare(actual,arrays[0])


class FreeIndexRemapTests(unittest.TestCase):
    def test_inverse_views_collapse_without_mutating_original_graph(self):
        base=PolyPath.from_coeffs([np.array([1.,.3]),np.array([2.,.5]),np.array([3.,.7])],steps=0).scale(2.)
        permutation=(2,0,1);inverse=(1,2,0)
        native=base.permute(permutation).permute(inverse)
        before=hash(native)
        token=PolyPath._index_reduction_context.set(True)
        try:
            optimized=base.permute(permutation).permute(inverse)
            normalized=ce._normalize_index_paths(native,{})
        finally:PolyPath._index_reduction_context.reset(token)
        self.assertIs(optimized,base)
        self.assertIs(normalized,base)
        self.assertEqual(hash(native),before)
        states=np.array([[0,1,2],[3,0,1],[2,3,0]])
        np.testing.assert_array_equal(native.evaluate_polynomial(states),optimized.evaluate_polynomial(states))

    def test_leaf_padding_and_permutation_are_exact_native_views(self):
        base=PolyPath.from_coeffs([np.array([1.,.3]),np.array([2.,.5])],steps=0)
        expected=base.pad([1,1]).permute((3,1,0,2))
        token=PolyPath._index_reduction_context.set(True)
        try:actual=base.pad([1,1]).permute((3,1,0,2))
        finally:PolyPath._index_reduction_context.reset(token)
        self.assertIsNone(actual._node)
        states=np.array([[0,1,2,3],[3,0,1,2],[2,3,0,1]])
        np.testing.assert_array_equal(actual.evaluate_polynomial(states),expected.evaluate_polynomial(states))

    def test_real_lazy_path_orbits_share_the_lowered_plan(self):
        from Psience.VPT2.Analytic import PTEnergyChangeProductSum,PTTensorCoeffProductSum
        wrapper,states,coeffs,freqs=FreeIndexExecutorTests().fixture()
        a,b=list(wrapper.expr.terms)
        value=wrapper.expr.terms[a]
        value=PTEnergyChangeProductSum({energy:PolyPath.from_polynomial(poly) for energy,poly in value.terms.items()})
        wrapper.expr=PTTensorCoeffProductSum({a:value,b:value.permute((1,0))},canonicalize=False)
        # Production wrappers supply a symbolic operator identity. A bare
        # PolyPath cannot currently be used as an alias fingerprint.
        wrapper.op=('free-index-lazy-binding-fixture',)
        baseline=ce.AnalyticExpressionStore().compile(wrapper,context={},mode_dimension=4)
        candidate=ce.AnalyticExpressionStore().compile(wrapper,context={},mode_dimension=4,free_index_reduction=True)
        metadata=next(iter(candidate.parts.values()))._free_index_plan_data
        self.assertEqual(metadata['stats']['canonical_plans'],1)
        expected=baseline.evaluate(states,coeffs,freqs,resume=False)
        actual=candidate.evaluate(states,coeffs,freqs,resume=False)
        for x,y in zip(actual,expected):FreeIndexExecutorTests().compare(x,y)


class FreeIndexBypassTests(unittest.TestCase):
    observations={}
    def test_known_low_rank_parts_do_no_canonical_or_graph_work(self):
        from Psience.VPT2.Analytic import PTTensorCoeffProductSum as Coeff
        for rank,fixed in ((0,0),(1,0),(1,1),(2,1),(2,2),(3,3),(4,3)):
            key=() if rank==0 else ((1,0)+tuple(range(rank)),)
            expression=Coeff({key:PolyPath.from_coeffs([np.ones(1)]*rank)},canonicalize=False)
            expression._free_index_plan_data={'stale':True}
            with mock.patch.object(ce,'_canonical_free_index_key',side_effect=AssertionError('orbit scan')), \
                 mock.patch.object(ce,'_normalize_index_paths',side_effect=AssertionError('graph walk')), \
                 mock.patch.object(ce,'_free_index_value_key',side_effect=AssertionError('graph hash')):
                self.assertIsNone(ce._prepare_free_index_plan(expression,fixed,rank=rank))
            self.assertFalse(hasattr(expression,'_free_index_plan_data'))

    def test_mixed_ranks_keep_ineligible_values_native(self):
        from Psience.VPT2.Analytic import PTTensorCoeffProductSum as Coeff
        keys=((),((1,0,0),),((1,0,0,1),))
        values=(PolyPath.from_coeffs([]),PolyPath.from_coeffs([np.array([1.,.2])]),
                PolyPath.from_coeffs([np.array([1.,.2]),np.array([2.,.3])]))
        expression=Coeff(dict(zip(keys,values)),canonicalize=False)
        with mock.patch.object(ce,'_canonical_free_index_key',wraps=ce._canonical_free_index_key) as orbit, \
             mock.patch.object(ce,'_normalize_index_paths',wraps=ce._normalize_index_paths) as normalize:
            metadata=ce._prepare_free_index_plan(expression,0)
        self.assertEqual(orbit.call_count,1);self.assertEqual(normalize.call_count,1)
        self.assertEqual(metadata['stats']['native_passthrough_keys'],2)
        for key,value in zip(keys[:2],values[:2]):
            self.assertIs(metadata['values'][key],value);self.assertIsNone(metadata['bindings'][key])

    def test_compiled_bypassed_parts_match_native_values(self):
        from Psience.VPT2.Analytic import (PTTensorCoeffProductSum as Coeff,
            PTEnergyChangeProductSum,PerturbationTheoryExpressionEvaluator as Expression)
        _,states,coeffs,freqs=FreeIndexExecutorTests().fixture()
        records=[]
        for change,keys in (((),(((1,0,0,0,0),),)),((1,),(((1,0,0,0,0),),)),
                            ((1,),(((1,0,0,0,1),),((1,0,0,1,1),)))):
            terms={}
            for key in keys:
                rank=len({i for coefficient in key for i in coefficient[2:]})
                poly=PolyPath.from_coeffs([np.array([1.+i,.2]) for i in range(rank)],steps=0)
                terms[Coeff.canonical_key(key)]=PTEnergyChangeProductSum({():poly})
            wrapper=Expression(('low-free-axis-fixture',change,keys),Coeff(terms,canonicalize=False),change=change,logger=False)
            baseline=ce.AnalyticExpressionStore().compile(wrapper,context={'case':keys,'change':change},mode_dimension=4)
            with mock.patch.object(ce,'_canonical_free_index_key',side_effect=AssertionError('orbit scan')), \
                 mock.patch.object(ce,'_normalize_index_paths',side_effect=AssertionError('graph walk')), \
                 mock.patch.object(ce,'_free_index_value_key',side_effect=AssertionError('graph hash')):
                candidate=ce.AnalyticExpressionStore().compile(wrapper,context={'case':keys,'change':change},mode_dimension=4,free_index_reduction=True)
            self.assertTrue(candidate.parts)
            self.assertTrue(all(not hasattr(part,'_free_index_plan_data') for part in candidate.parts.values()))
            expected=baseline.evaluate(states,coeffs,freqs,resume=False)
            actual=candidate.evaluate(states,coeffs,freqs,resume=False)
            records.extend(FreeIndexExecutorTests().compare(a,b) for a,b in zip(actual,expected))
        type(self).observations['native_low_free_axis_values']=records


class SharedCanonicalValueTests(unittest.TestCase):
    observations={}
    def table(self,physical_blocks=None):
        from Psience.VPT2.IndexedEvaluator import ConcreteBlockTable,IndexedBlockEvaluator
        from Psience.VPT2.Analytic import PerturbationTheoryExpressionEvaluator as Expression,_StatePermutationBlockIdentity
        wrapper,*_=FreeIndexExecutorTests().fixture()
        compiled=ce.AnalyticExpressionStore().compile(wrapper,context={'shared':'table'},mode_dimension=4,free_index_reduction=True)
        part=next(iter(compiled.parts.values()));plan=IndexedBlockEvaluator.compile_plan(part,reduced=True)
        subplan=next(iter(plan.subexpressions.values()));table=ConcreteBlockTable()
        state=np.array([[0,2,1,3]]);freqs=np.array([.013,.021,.034,.049])
        for n,physical in enumerate(physical_blocks if physical_blocks is not None else (np.array([[0,1],[0,2]]),np.array([[0,2],[0,3]]))):
            substates=np.moveaxis(state[:,physical],0,1)
            tuples=[[tuple(row) for row in block] for block in substates]
            table.append(0,np.array([n,n+1]),np.array([n,n+1]),np.array([[2.+n,3.+n]]),
                ('source',n),subplan,physical,substates,tuples,freqs[physical],
                _StatePermutationBlockIdentity(tuples),np.zeros((1,2)),state)
        token=ce._current_free_index_reduction.set(True)
        try:executor=IndexedBlockEvaluator(Expression,part)
        finally:ce._current_free_index_reduction.reset(token)
        return executor,table
    def test_physical_pools_deduplicate_and_restore_original_routes(self):
        executor,table=self.table();joined,routes,aliases=executor._share_canonical_value_blocks(table)
        self.assertEqual(len(joined),1);self.assertEqual(len(joined.permutation_subsets[0]),3)
        for index,(group,inverse) in enumerate(routes):
            np.testing.assert_array_equal(joined.permutation_subsets[group][inverse],table.permutation_subsets[index])
            np.testing.assert_array_equal(joined.permutation_frequencies[group][inverse],table.permutation_frequencies[index])
            np.testing.assert_array_equal(joined.permutation_substates[group][inverse],table.permutation_substates[index])
        self.assertEqual(executor.stats.shared_value_duplicate_points,1)
        self.assertEqual(aliases[0],tuple(table.expression_keys))
    def test_join_budget_keeps_native_blocks_separate(self):
        executor,table=self.table();executor.max_shared_value_permutations=2
        joined,routes,aliases=executor._share_canonical_value_blocks(table)
        self.assertIs(joined,table);self.assertIsNone(routes)
    def test_disjoint_physical_bindings_keep_native_execution(self):
        executor,table=self.table((np.array([[0,1],[0,2]]),np.array([[0,3],[1,3]])))
        joined,routes,aliases=executor._share_canonical_value_blocks(table)
        self.assertIs(joined,table);self.assertIsNone(routes);self.assertIsNone(aliases)
        self.assertEqual(executor.stats.shared_value_blocks_eliminated,0)
        self.assertEqual(executor.stats.shared_value_duplicate_points,0)
        self.assertEqual(executor.stats.shared_value_disjoint_blocks,2)
    def test_actual_contraction_vectors_fall_with_identical_values(self):
        self.evaluate_case(False,False)
    def test_nonunit_prefactors_batched_states_and_cutoffs(self):
        self.evaluate_case(True,True)
    def evaluate_case(self,batched,cutoffs):
        from Psience.VPT2.IndexedEvaluator import IndexedBlockEvaluator
        from Psience.VPT2.Analytic import PTEnergyChangeProductSum,PTTensorCoeffProductSum
        wrapper,states,coeffs,freqs=FreeIndexExecutorTests().fixture()
        if batched:states=[(state,np.array([[0,1,2,3]])) for state,_ in states]
        wrapper.expr=PTTensorCoeffProductSum({key:PTEnergyChangeProductSum(dict(value.terms),prefactor=2.3)
            for key,value in wrapper.expr.terms.items()},canonicalize=False)
        if cutoffs:
            for index in set(itertools.permutations((0,0,1))):coeffs[0][1][0][index]=1e-32
            for index in set(itertools.permutations((0,0,2))):coeffs[0][1][0][index]=0.
        compiled=ce.AnalyticExpressionStore().compile(wrapper,context={'shared':'numeric','batched':batched},mode_dimension=4,free_index_reduction=True)
        outputs=[];statistics=[]
        native_finish=IndexedBlockEvaluator.finish
        for enabled in (False,True):
            captured=[]
            def finish(executor):
                result=native_finish(executor);captured.append(result);return result
            with mock.patch.object(IndexedBlockEvaluator,'use_shared_canonical_values',enabled), \
                 mock.patch.object(IndexedBlockEvaluator,'finish',finish):
                outputs.append(compiled.evaluate(states,coeffs,freqs,resume=False,zero_cutoff=1e-24))
            statistics.append(captured)
        residual=[FreeIndexExecutorTests().compare(a,b) for a,b in zip(outputs[1],outputs[0])]
        before=sum(row['contraction_rows'] for row in statistics[0]);after=sum(row['contraction_rows'] for row in statistics[1])
        if batched:
            self.assertEqual(after,before)
            self.assertEqual(sum(row['shared_value_duplicate_points'] for row in statistics[1]),0)
        else:
            self.assertLess(after,before)
            self.assertGreater(sum(row['shared_value_blocks_eliminated'] for row in statistics[1]),0)
        self.assertLessEqual(sum(row['energy_scalar_lanes'] for row in statistics[1]),sum(row['energy_scalar_lanes'] for row in statistics[0]))
        type(self).observations['batched' if batched else 'routed']=dict(residual=residual,rows_before=before,rows_after=after,statistics=statistics)


class SharedEnergyRowTests(unittest.TestCase):
    observations = {}

    def evaluate(self, capacity=32*1024**2):
        from Psience.VPT2.IndexedEvaluator import IndexedBlockEvaluator
        from Psience.VPT2.Analytic import PTTensorCoeffProductSum as Coeff
        wrapper, states, coeffs, freqs = FreeIndexExecutorTests().fixture()
        value = next(iter(wrapper.expr.terms.values()))
        a = Coeff.canonical_key(((1,0,0,0,1),))
        b = Coeff.canonical_key(((1,1,0,0,1),))
        coeffs[0][1][1] = coeffs[0][1][0] * .37
        wrapper.expr = Coeff({a:value, b:value.scale(1.7)}, canonicalize=False)
        compiled = ce.AnalyticExpressionStore().compile(wrapper, context={'energy_row_reuse':capacity},
            mode_dimension=4, free_index_reduction=True)
        outputs, statistics = [], []
        original = IndexedBlockEvaluator.finish
        for enabled in (False,True):
            captured = []
            def finish(executor):
                result=original(executor);captured.append(result);return result
            with mock.patch.object(IndexedBlockEvaluator,'use_shared_energy_rows',enabled), \
                 mock.patch.object(IndexedBlockEvaluator,'max_energy_row_cache_bytes',capacity), \
                 mock.patch.object(IndexedBlockEvaluator,'finish',finish):
                outputs.append(compiled.evaluate(states,coeffs,freqs,resume=False,zero_cutoff=1e-24))
            statistics.append(captured)
        residual = [FreeIndexExecutorTests().compare(a,b) for a,b in zip(outputs[1],outputs[0])]
        before = sum(r['energy_computed_scalar_lanes'] for r in statistics[0])
        after = sum(r['energy_computed_scalar_lanes'] for r in statistics[1])
        self.assertLess(after,before)
        self.assertGreater(sum(r['energy_rows_reused'] for r in statistics[1]),0)
        self.assertLessEqual(max(r['max_energy_row_cache_bytes'] for r in statistics[1]),capacity)
        type(self).observations[str(capacity)] = dict(residual=residual,scalar_lanes_before=before,
            scalar_lanes_after=after,statistics=statistics)

    def test_distinct_numerators_share_exact_denominators(self):
        self.evaluate()

    def test_bounded_cache_keeps_values_and_temporary_row_reuse(self):
        self.evaluate(1)


if __name__=='__main__':unittest.main()
