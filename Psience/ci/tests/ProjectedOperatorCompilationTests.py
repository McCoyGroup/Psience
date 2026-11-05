"""Regress projected operator compilation and resonance-aware cache reuse."""
import contextlib,itertools,tempfile,unittest
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
    def runner(self,folder):
        w=np.array([.0201,.0217,.0313]);v3=np.zeros((3,)*3);v4=np.zeros((3,)*4)
        for indices,value in [((0,0,0),.0003),((0,1,1),.00008),((0,0,1),.0001),((0,1,2),.00007)]:
            for p in set(itertools.permutations(indices)):v3[p]=value
        for p in set(itertools.permutations((0,0,1,1))):v4[p]=.00002
        mu1=np.array([.2,.1,.3]);mu2=np.diag([.02,-.01,.01]);mu3=np.zeros((3,)*3);mu3[0,0,0]=.0004
        dipole=[[0.,mu1,mu2,mu3]]+[[0.,np.zeros(3),np.zeros((3,3)),np.zeros((3,)*3)] for _ in range(2)]
        return AnalyticVPTRunner([[np.diag(w)/2,np.diag(w)/2],[v3,0],[v4,0,0,0,0]],
            order=2,freqs=w,internals=False,dipole_expansion=dipole,expressions_folder=folder,logger=False)
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
    def test_resonance_block_changes_create_distinct_cache_entries(self):
        with tempfile.TemporaryDirectory() as folder:
            a=self.moments(self.runner(folder),[[1,0,0],[0,1,0]]);first=inventory(folder)
            again=self.moments(self.runner(folder),[[1,0,0],[0,1,0]]);repeat=inventory(folder)
            b=self.moments(self.runner(folder),[[1,0,0],[0,0,1]]);changed=inventory(folder)
            np.testing.assert_array_equal(again,a)
            def parts(value):return {r['path']:(r['sha256'],r['mtime_ns']) for r in value['file_records'] if r['kind']=='compiled_expression_parts'}
            self.assertEqual(parts(first),parts(repeat))
            added=set(changed['expression_signatures'])-set(first['expression_signatures'])
            self.assertTrue(added)
            self.observations.update(same_resonance_reuses_expression_parts=True,
                different_resonance_added_signatures=len(added),first=first,repeat=repeat,changed=changed,
                max_transition_moment_change=float(np.max(np.abs(a-b))))

if __name__=='__main__':unittest.main()
