"""Regress projected operator compilation and resonance-aware cache reuse."""
import contextlib,itertools,tempfile,unittest,pickle
from pathlib import Path
from unittest import mock
import numpy as np
from Psience.VPT2 import AnalyticVPTRunner
from Psience.VPT2 import CompiledExpressions as ce
from Psience.VPT2.Projected import ProjectedCorrectionEvaluator
from pathlib import Path
import hashlib,json,collections
def inventory(folder):
    root=Path(folder);rows=[];groups=collections.defaultdict(lambda:dict(files=0,logical_bytes=0,allocated_bytes=0))
    for path in sorted(root.rglob('*')):
        if not path.is_file():continue
        rel=path.relative_to(root);stat=path.stat()
        if 'evaluations' in rel.parts:kind='numerical_evaluation_checkpoints'
        elif path.name.startswith('part-') and path.suffix=='.pkl':kind='compiled_expression_parts'
        elif path.name=='manifest.json':kind='expression_manifests'
        elif rel.parts[0]=='aliases':kind='expression_aliases'
        else:kind='other'
        data=path.read_bytes();allocated=stat.st_blocks*512
        rows.append(dict(path=str(rel),kind=kind,logical_bytes=stat.st_size,allocated_bytes=allocated,mtime_ns=stat.st_mtime_ns,sha256=hashlib.sha256(data).hexdigest()))
        item=groups[kind];item['files']+=1;item['logical_bytes']+=stat.st_size;item['allocated_bytes']+=allocated
    return dict(files=len(rows),logical_bytes=sum(r['logical_bytes'] for r in rows),allocated_bytes=sum(r['allocated_bytes'] for r in rows),
        categories=dict(groups),largest_file_bytes=max((r['logical_bytes'] for r in rows),default=0),file_records=rows,
        expression_signatures=sorted(p.parent.name for p in root.glob('*/manifest.json')))


class ProjectedOperatorCompilationTests(unittest.TestCase):
    observations={}
    def test_projection_factor_is_only_a_formal_role(self):
        from Psience.VPT2 import Projected as projected
        factor=projected.ProjectionFactor('R',3)
        self.assertEqual(vars(factor),dict(kind='R',slot=3))
        self.assertEqual(vars(pickle.loads(pickle.dumps(factor))),vars(factor))
        self.assertFalse(hasattr(projected,'StateProjectionFactor'))
        with self.assertRaises(ValueError):projected.ProjectionFactor('R',[(1,0)])

    def test_equations_rebind_states_frequencies_and_block_coefficients(self):
        from Psience.VPT2.Projected import AnalyticBlockExpansion
        r=self.runner(None);h=r.eval.solver.hamiltonian_expansion[1]
        a=AnalyticBlockExpansion(r.eval,[[0,0,0]])
        b=AnalyticBlockExpansion(r.eval,[[0,0,0],[1,0,0]],freqs=[.022,.024,.033])
        ta=a._scalar(.7)*a._resolvent(h,(0,0,0),(0,0,0))
        tb=b._scalar(1.3)*b._resolvent(h,(0,0,0),(0,0,0))
        ga,ba,va=a._generic_term(ta,len(a.expansions))
        gb,bb,vb=b._generic_term(tb,len(b.expansions))
        self.assertEqual(ce._fingerprint(ga),ce._fingerprint(gb))
        self.assertNotEqual(ba,bb);self.assertNotEqual(va,vb)
        before=ce._fingerprint(ga)
        generic_polynomial=pickle.dumps(ga.get_poly_terms((1,)))
        self.assertNotIn(b'StateProjectionFactor',generic_polynomial)
        self.assertNotIn(b'_ProjectionEvaluationContext',generic_polynomial)
        expected=a.element(ta,[0,0,0],[1,0,0])
        self.assertGreater(abs(expected),1e-10)
        self.assertEqual(b.element(tb,[0,0,0],[1,0,0]),0.)
        with tempfile.TemporaryDirectory() as folder:
            store=ce.AnalyticExpressionStore(folder)
            with ce.compiled_evaluation_context(store):
                aa=AnalyticBlockExpansion(r.eval,[[0,0,0]])
                ab=AnalyticBlockExpansion(r.eval,[[0,0,0],[1,0,0]],freqs=b.freqs)
                actual=aa.element(aa._scalar(.7)*aa._resolvent(h,(0,0,0),(0,0,0)),[0,0,0],[1,0,0])
                first=inventory(folder)
                self.assertAlmostEqual(actual,expected,places=13)
                excluded=ab.element(ab._scalar(1.3)*ab._resolvent(h,(0,0,0),(0,0,0)),[0,0,0],[1,0,0])
                self.assertEqual(excluded,0.)
                repeated=aa.element(aa._scalar(1.4)*aa._resolvent(h,(0,0,0),(0,0,0)),[0,0,0],[1,0,0])
                self.assertAlmostEqual(repeated,2*expected,places=13)
                zero=aa.element(aa._scalar(0.)*aa._resolvent(h,(0,0,0),(0,0,0)),[0,0,0],[1,0,0])
                self.assertEqual(zero,0.)
            changed=inventory(folder)
            self.assertEqual(first['expression_signatures'],changed['expression_signatures'])
            self.assertGreater(changed['categories']['numerical_evaluation_checkpoints']['files'],
                               first['categories']['numerical_evaluation_checkpoints']['files'])
        self.assertEqual(before,ce._fingerprint(ga))

    def test_projector_row_batch_matches_scalar_reference(self):
        from Psience.VPT2.Projected import ProjectionFactor,projected_path_factors,projected_energy_rows
        bindings=((((1,0,0),(0,1,0)),(.02,.023,.031),(1,0,0),False),)
        permutations=np.asarray(list(itertools.permutations(range(3))))
        state=np.array([[0,0,0]]);freqs=np.asarray(bindings[0][1])[permutations]
        keys=[((ProjectionFactor(kind,0),1,0,0),) for kind in ['P','Q','R','V']]
        expected=[projected_path_factors(key,state,permutations,freqs,projection_bindings=bindings) for key in keys]
        masks,denominators=projected_energy_rows(keys,[state]*4,[permutations]*4,
            np.stack([freqs]*4),projection_bindings=bindings)
        for i,(mask,denominator) in enumerate(expected):
            np.testing.assert_array_equal(masks[i],mask)
            np.testing.assert_allclose(denominators[i],denominator,atol=1e-15,rtol=1e-14)

    def test_trivial_resolvent_retains_shifted_reference_energy(self):
        from Psience.VPT2.Projected import ProjectionFactor,projected_path_factors,projected_energy_rows
        binding=((((0,1,0),),(.02,.023,.031),(0,1,0),True),)
        key=((ProjectionFactor('R'),0,0),)
        state=np.array([[0,0,0]]);perms=np.array([[1,0]]);freqs=np.array([[.023,.02]])
        mask,denominator=projected_path_factors(key,state,perms,freqs,projection_bindings=binding)
        self.assertTrue(mask[0]);self.assertAlmostEqual(denominator[0],.023,places=15)
        masks,denominators=projected_energy_rows([key],[state],[perms],freqs[None],projection_bindings=binding)
        np.testing.assert_array_equal(masks,[mask]);np.testing.assert_allclose(denominators,[denominator])

    def test_projection_bindings_cross_compiled_worker_processes(self):
        from McUtils.Parallelizers import MultiprocessingParallelizer
        from Psience.VPT2.Projected import AnalyticBlockExpansion
        r=self.runner(None);h=r.eval.solver.hamiltonian_expansion[1]
        a=AnalyticBlockExpansion(r.eval,[[0,0,0]])
        expected=a.element(a._scalar(.7)*a._resolvent(h,(0,0,0),(0,0,0)),[0,0,0],[1,0,0])
        with MultiprocessingParallelizer(processes=2,stall_timeout=60) as backend:
            r.eval.parallelizer=backend
            b=AnalyticBlockExpansion(r.eval,[[0,0,0]])
            with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()):
                actual=b.element(b._scalar(.7)*b._resolvent(h,(0,0,0),(0,0,0)),[0,0,0],[1,0,0])
        self.assertAlmostEqual(actual,expected,places=13)

    def runner(self,folder):
        w=np.array([.0201,.0217,.0313]);v3=np.zeros((3,)*3);v4=np.zeros((3,)*4)
        for indices,value in [((0,0,0),.0003),((0,1,1),.00008),((0,0,1),.0001),((0,1,2),.00007)]:
            for p in set(itertools.permutations(indices)):v3[p]=value
        for p in set(itertools.permutations((0,0,1,1))):v4[p]=.00002
        mu1=np.array([.2,.1,.3]);mu2=np.diag([.02,-.01,.01]);mu3=np.zeros((3,)*3);mu3[0,0,0]=.0004
        dipole=[[0.,mu1,mu2,mu3]]+[[0.,np.zeros(3),np.zeros((3,3)),np.zeros((3,)*3)] for _ in range(2)]
        return AnalyticVPTRunner([[np.diag(w)/2,np.diag(w)/2],[v3,0],[v4,0,0,0,0]],
            order=2,freqs=w,internals=False,dipole_expansion=dipole,expressions_folder=folder,
            compile=folder is not None,logger=False)
    def moments(self,runner,block):
        states=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
        space=runner.prep_states(states,degeneracy_specs={'groups':[block]})
        return np.asarray(runner.get_transition_moment_corrections(space,order=2,state_specific_degeneracies=True))
    def test_projected_branch_uses_compilation_and_matches_uncompiled(self):
        block=[[1,0,0],[0,1,0]]
        reference=self.moments(self.runner(None),block)
        native=ProjectedCorrectionEvaluator.operators;contexts=[]
        def observed(obj,*args,**kwargs):
            contexts.append(ce.current_compilation());return native(obj,*args,**kwargs)
        with tempfile.TemporaryDirectory() as folder, mock.patch.object(ProjectedCorrectionEvaluator,'operators',observed):
            actual=self.moments(self.runner(folder),block)
            self.assertTrue(contexts);self.assertTrue(all(c is not None for c in contexts))
            self.assertTrue(list(Path(folder).glob('*/manifest.json')))
        np.testing.assert_allclose(actual,reference,rtol=2e-11,atol=2e-13)
        self.assertIsNone(ce.current_compilation())
    def test_without_cache_folder_keeps_uncompiled_projected_path(self):
        native=ProjectedCorrectionEvaluator.operators;contexts=[]
        def observed(obj,*args,**kwargs):
            contexts.append(ce.current_compilation());return native(obj,*args,**kwargs)
        with mock.patch.object(ProjectedCorrectionEvaluator,'operators',observed):
            result=self.moments(self.runner(None),[[1,0,0],[0,1,0]])
        self.assertTrue(np.isfinite(result).all());self.assertTrue(contexts)
        self.assertTrue(all(c is None for c in contexts))
    def test_existing_context_is_reused(self):
        native=ProjectedCorrectionEvaluator.operators;contexts=[]
        def observed(obj,*args,**kwargs):
            contexts.append(ce.current_compilation());return native(obj,*args,**kwargs)
        with tempfile.TemporaryDirectory() as outer, tempfile.TemporaryDirectory() as inner:
            store=ce.AnalyticExpressionStore(outer)
            with ce.compiled_evaluation_context(store),mock.patch.object(ProjectedCorrectionEvaluator,'operators',observed):
                self.moments(self.runner(inner),[[1,0,0],[0,1,0]])
            self.assertTrue(contexts);self.assertTrue(all(c[0] is store for c in contexts))
            self.assertFalse(list(Path(inner).iterdir()));self.assertTrue(list(Path(outer).glob('*/manifest.json')))
        self.assertIsNone(ce.current_compilation())
    def test_context_restored_when_projected_evaluation_raises(self):
        with tempfile.TemporaryDirectory() as folder,mock.patch.object(ProjectedCorrectionEvaluator,'operators',side_effect=RuntimeError('fixture failure')):
            with self.assertRaisesRegex(RuntimeError,'fixture failure'):
                self.moments(self.runner(folder),[[1,0,0],[0,1,0]])
        self.assertIsNone(ce.current_compilation())
    def test_resonance_block_changes_keep_common_equations_immutable(self):
        with tempfile.TemporaryDirectory() as folder:
            a=self.moments(self.runner(folder),[[1,0,0],[0,1,0]]);first=inventory(folder)
            again=self.moments(self.runner(folder),[[1,0,0],[0,1,0]]);repeat=inventory(folder)
            b=self.moments(self.runner(folder),[[1,0,0],[0,0,1]]);changed=inventory(folder)
            np.testing.assert_array_equal(again,a)
            def parts(value):return {r['path']:(r['sha256'],r['mtime_ns']) for r in value['file_records'] if r['kind']=='compiled_expression_parts'}
            self.assertEqual(parts(first),parts(repeat))
            added=set(changed['expression_signatures'])-set(first['expression_signatures'])
            self.assertTrue(set(first['expression_signatures']) & set(changed['expression_signatures']))
            self.assertTrue(parts(first).items() <= parts(changed).items())
            self.observations.update(same_resonance_reuses_expression_parts=True,
                different_resonance_added_signatures=len(added),first=first,repeat=repeat,changed=changed,
                max_transition_moment_change=float(np.max(np.abs(a-b))))

if __name__=='__main__':unittest.main()
