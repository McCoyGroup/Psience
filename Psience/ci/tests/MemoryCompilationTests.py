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
    def test_new_resonance_block_does_not_reuse_the_previous_block_expression(self):
        r=self.runner();self.moments(r);before=set(r.eval.expression_store._memory)
        block=[[1,0,0],[0,0,1]]
        with self.no_expression_writes():actual=self.moments(r,block)
        expected=self.moments(self.runner(compile=False),block)
        np.testing.assert_allclose(actual,expected,rtol=2e-11,atol=2e-13)
        self.assertTrue(set(r.eval.expression_store._memory)-before)
    def test_hamiltonian_clone_preserves_policy_and_derives_new_harmonic_frequencies(self):
        r=self.runner(compile=False);new_frequencies=np.array([.04,.06,.08])
        expansions=[[np.diag(new_frequencies)/2,np.diag(new_frequencies)/2],[0,0],[0,0,0,0,0]]
        with mock.patch.object(r.eval,'_prep_ham_corrs',return_value=([],[])),\
             mock.patch.object(r.eval,'_prep_operator_expansion',return_value=expansions),\
             mock.patch.object(r.eval.solver,'modify_hamiltonian',return_value=r.eval.solver):
            changed=r.eval.modify_hamiltonian(object())
        np.testing.assert_array_equal(changed.freqs,new_frequencies)
        self.assertFalse(changed.compile_expressions)

if __name__=='__main__':unittest.main()
