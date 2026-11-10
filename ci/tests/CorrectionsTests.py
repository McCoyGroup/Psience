"""Real-data tests for correction storage, filtering, and block evaluation."""
import ast, gc, hashlib, importlib, json, os, tempfile, unittest, weakref
from pathlib import Path
from unittest.mock import patch
import numpy as np
from McUtils.Numputils import RowSparseArray, SparseArray, matrix_dot_into
from McUtils.Scaffolding import NullLogger
from Psience.Molecools import Molecule
from Psience.BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis, BasisStateSpaceFilter
from Psience.VPT2 import VPTRunner
from Psience.VPT2.Corrections import PerturbationTheoryCorrections
from Psience.VPT2.Solver import PerturbationTheorySolver
from Psience.VPT2.Runner import VPTStateSpace
from Psience.VPT2.DegeneracySpecs import (
    DegenerateMultiStateSpace, GroupsDegeneracySpec, MartinTestDegeneracySpec)

# Tests intentionally inspect the internal collectors' buffer ownership.
SparseCorrectionArray = PerturbationTheoryCorrections._SparseCorrectionArray
CompactCorrectionArray = PerturbationTheoryCorrections._CompactCorrectionArray

def fixture(name):
    root=os.environ.get("PSIENCE_REAL_VPT_INPUTS")
    roots=([Path(root)] if root else [])+[Path(__file__).parent/"TestData",Path(__file__).parent/"validation_data"]
    return next((p/name for p in roots if (p/name).is_file()), None)


@unittest.skipUnless(fixture("hno3.fchk"), "Needs actual hno3.fchk; set PSIENCE_REAL_VPT_INPUTS to its directory")
class CompactCorrectionStorageTests(unittest.TestCase):
    def test_real_two_quantum_corrections_and_intensities(self):
        results=[]
        for storage in ["sparse", "compact"]:
            molecule=Molecule.from_file(str(fixture("hno3.fchk")))
            wave=VPTRunner.run_simple(molecule,2,order=2,expansion_order=2,
                evaluator_defaults="new",degeneracy_specs="auto",
                correction_storage=storage,correction_block_size=4096,
                state_space_storage="memory",mixed_derivative_handling_mode="analytical",
                logger=NullLogger())
            results.append(wave)
        reference,candidate=results
        np.testing.assert_array_equal(reference.corrs.states.excitations,candidate.corrs.states.excitations)
        np.testing.assert_array_equal(reference.corrs.total_basis.indices,candidate.corrs.total_basis.indices)
        np.testing.assert_array_equal(reference.corrs.energy_corrs,candidate.corrs.energy_corrs)
        np.testing.assert_allclose(reference.energies,candidate.energies,rtol=2e-11,atol=2e-12)
        a,b=reference.intensities,candidate.intensities
        if isinstance(a,dict):
            self.assertEqual(set(a),set(b))
            for key in a:np.testing.assert_allclose(a[key],b[key],rtol=2e-11,atol=2e-12)
        else:np.testing.assert_allclose(a,b,rtol=2e-11,atol=2e-12)
        def entries(matrix):
            blocks=list(matrix.iter_nonzero_blocks(block_size=4096))
            values=np.concatenate([v for v,inds in blocks]) if blocks else np.empty(0)
            indices=[np.concatenate([inds[d] for v,inds in blocks]) if blocks else np.empty(0,dtype=int) for d in range(2)]
            return values,indices
        for expected,actual in zip(reference.corrs.wfn_corrections,candidate.corrs.wfn_corrections):
            self.assertIsInstance(actual,RowSparseArray)
            av,ai=entries(expected);bv,bi=entries(actual)
            np.testing.assert_array_equal(av,bv)
            for a,b in zip(ai,bi):np.testing.assert_array_equal(a,b)
            selected=actual[[0,1,1],:]
            native=selected.asarray()
            np.testing.assert_allclose(selected.dot(selected.T).asarray(),native.dot(native.T),rtol=2e-12,atol=2e-13)
            np.testing.assert_allclose(SparseArray.from_data(native).dot(selected.T).asarray(),native.dot(native.T),rtol=2e-12,atol=2e-13)

    def test_sparse_lower_orders_and_discovery_release_on_actual_HNO3(self):
        import scipy.sparse as sp
        captures=[];owners=[];extensions=[]
        native_groups=PerturbationTheorySolver._strong_coupling_groups
        native_extend=PerturbationTheorySolver._extend_strong_coupling_spaces
        native_block=PerturbationTheorySolver._apply_VPT_block_equations_in_place
        def groups(solver,spec,sc,corrs,threshold):
            matrices=corrs.wfn_corrections
            owners[:]=[weakref.ref(m) for m in matrices]
            for matrix in matrices:
                for row in getattr(matrix,'rows',()):
                    owners.extend([weakref.ref(row),weakref.ref(row.values)])
            return native_groups(solver,spec,sc,corrs,threshold)
        def extend(solver,*args,**kwargs):
            extensions.append((len(owners),sum(owner() is not None for owner in owners)))
            return native_extend(solver,*args,**kwargs)
        def block(solver,group,*args,**kwargs):
            result=native_block(solver,group,*args,**kwargs)
            if solver.sparse_low_order_corrections:
                self.assertIsInstance(result[2],tuple)
                self.assertTrue(sp.issparse(result[2][0]))
                self.assertTrue(sp.issparse(result[2][1]))
                self.assertEqual(result[2][1].dtype,np.dtype(float))
                if not captures:
                    H=list(kwargs.get('perturbations') or solver.representations)
                    # Same actual molecular block/denominators, narrower input
                    # storage. This checks accumulation dtype, not a new model.
                    H[1]=SparseArray.from_data(H[1].data.astype(np.float32))
                    solver.sparse_low_order_corrections=False
                    dense=native_block(solver,group,perturbations=H)
                    solver.sparse_low_order_corrections=True
                    mixed=native_block(solver,group,perturbations=H)
                    actual=np.stack([x.toarray() if sp.issparse(x) else x for x in mixed[2]],axis=1)
                    np.testing.assert_allclose(actual,dense[2],rtol=2e-11,atol=2e-12)
                    np.testing.assert_allclose(mixed[0],dense[0],rtol=2e-11,atol=2e-12)
                captures.append(tuple(x.shape for x in result[2]))
            else:self.assertIsInstance(result[2],np.ndarray)
            return result
        waves=[]
        with patch.object(PerturbationTheorySolver,'_strong_coupling_groups',groups), \
                patch.object(PerturbationTheorySolver,'_extend_strong_coupling_spaces',extend), \
                patch.object(PerturbationTheorySolver,'_apply_VPT_block_equations_in_place',block):
            for sparse in [False,True]:
                molecule=Molecule.from_file(str(fixture('hno3.fchk')))
                molecule.internals=None
                waves.append(VPTRunner.run_simple(molecule,2,order=2,expansion_order=2,
                    evaluator_defaults='new',degeneracy_specs='auto',use_block_recurrence=True,
                    mixed_derivative_handling_mode='analytical',memory_constrained=True,
                    state_space_storage='memory',correction_storage='compact',
                    sparse_low_order_corrections=sparse,logger=False))
        self.assertTrue(captures,'No actual multi-state recurrence returned')
        self.assertTrue(extensions,'No automatic resonance pass was checked')
        self.assertTrue(all(count>0 and alive==0 for count,alive in extensions),extensions)
        a,b=waves
        np.testing.assert_array_equal(a.corrs.states.indices,b.corrs.states.indices)
        np.testing.assert_array_equal(a.corrs.total_basis.indices,b.corrs.total_basis.indices)
        ga=[g.indices.tolist() for g in a.corrs.degenerate_states]
        gb=[g.indices.tolist() for g in b.corrs.degenerate_states]
        self.assertEqual(ga,gb)
        np.testing.assert_allclose(a.energies,b.energies,rtol=2e-11,atol=2e-12)
        np.testing.assert_allclose(a.corrs.energy_corrs,b.corrs.energy_corrs,rtol=2e-11,atol=2e-12)
        for left,right in zip(a.corrs.wfn_corrections,b.corrs.wfn_corrections):
            np.testing.assert_allclose(left.asarray(),right.asarray(),rtol=2e-11,atol=2e-12)
        if isinstance(a.intensities,dict):
            self.assertEqual(set(a.intensities),set(b.intensities))
            for key in a.intensities:
                np.testing.assert_allclose(a.intensities[key],b.intensities[key],rtol=2e-11,atol=2e-12)
        else:np.testing.assert_allclose(a.intensities,b.intensities,rtol=2e-11,atol=2e-12)



@unittest.skipUnless(fixture('hno3.fchk'),'Needs actual hno3.fchk; set PSIENCE_REAL_VPT_INPUTS')
class FileBackedCorrectionTests(unittest.TestCase):
    def test_complete_two_quantum_spectrum_and_coefficients(self):
        from Psience.Molecools import Molecule
        from Psience.VPT2 import VPTRunner
        from Psience.BasisReps import BasisStateSpace,HarmonicOscillatorProductBasis
        from McUtils.Scaffolding import NullLogger
        with tempfile.TemporaryDirectory() as directory:
            waves=[]
            for spill in [False,True]:
                opts=dict(order=2,expansion_order=2,evaluator_defaults='new',degeneracy_specs='auto',
                    mixed_derivative_handling_mode='analytical',state_space_storage='external',
                    state_space_scratch=str(Path(directory)/str(spill)),correction_storage='compact',logger=NullLogger())
                if spill:opts.update(correction_spill_bytes=1024**3,correction_resident_bytes=1,
                    correction_spill_segment_bytes=4*1024**2)
                waves.append(VPTRunner.run_simple(Molecule.from_file(str(fixture('hno3.fchk'))),2,**opts))
            reference,candidate=waves
            for attr in ['energies']:
                np.testing.assert_array_equal(getattr(reference,attr),getattr(candidate,attr))
            np.testing.assert_array_equal(reference.corrs.states.excitations,candidate.corrs.states.excitations)
            np.testing.assert_array_equal(reference.corrs.energy_corrs,candidate.corrs.energy_corrs)
            a,b=reference.intensities,candidate.intensities
            if isinstance(a,dict):
                self.assertEqual(set(a),set(b))
                for key in a:np.testing.assert_array_equal(a[key],b[key])
            else:np.testing.assert_array_equal(a,b)
            expected=BasisStateSpace.from_quanta(HarmonicOscillatorProductBasis(9),[1,2]).excitations
            self.assertTrue({tuple(x) for x in expected}<={tuple(x) for x in candidate.corrs.states.excitations})
            def hashes(matrix):
                digests=[hashlib.sha256() for _ in range(3)]
                for values,(rows,columns) in matrix.iter_nonzero_blocks(block_size=4096):
                    for digest,array in zip(digests,[values,rows,columns]):digest.update(np.asarray(array).tobytes())
                return [d.hexdigest() for d in digests]
            for a,b in zip(reference.corrs.wfn_corrections,candidate.corrs.wfn_corrections):
                self.assertEqual(hashes(a),hashes(b))
                self.assertGreater(b.storage.snapshot()['spilled_rows'],0)

            # Exercise the shared coupling-entry filter on actual calculated
            # coefficients, with the native state-space prefilter interface.
            corrs=reference.corrs
            dense=np.stack([matrix.asarray() for matrix in corrs.wfn_corrections],axis=1)
            nstates,norders,nbasis=dense.shape
            zero_rule=np.zeros(corrs.states.ndim,dtype=int)
            rule=BasisStateSpaceFilter(corrs.states,[(corrs.states,[zero_rule])],None)
            filters={(level,):rule for level in range(norders)}
            expected,expected_spaces=PerturbationTheoryCorrections.create_coupling_matrix(
                dense,corrs.states,corrs.total_basis,nstates,norders-1,filters=filters)
            for mode in ['sparse','compact']:
                packed=PerturbationTheoryCorrections.create_correction_array(
                    mode,nstates,norders-1,nbasis,4096)
                packed[:]=dense
                actual,actual_spaces=PerturbationTheoryCorrections.create_coupling_matrix(
                    packed,corrs.states,corrs.total_basis,nstates,norders-1,filters=filters)
                for a,b in zip(expected,actual):
                    np.testing.assert_array_equal(a.asarray(),b.asarray())
                for a,b in zip(expected_spaces,actual_spaces):
                    for x,y in zip(a.spaces,b.spaces):
                        np.testing.assert_array_equal(x.indices,y.indices)


class VPTBlockHandoffTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(__file__).parent/'TestData'/'pyrene_vpt2_block_handoff.json'
        cls.fixture=json.loads(path.read_text())
        cls.coefficients=np.asarray(cls.fixture['coefficients'],dtype=float)
        assert cls.coefficients.ndim==3 and len(cls.coefficients)>=2

    def data(self):
        # Restore the solver's order/state/basis backing layout from real data.
        return self.coefficients.transpose(1,0,2).copy().transpose(1,0,2)

    def collector(self, cls, data):
        reference=weakref.ref(data)
        class RecordingRows(cls):
            def __setitem__(self, selected, corrections):
                original=reference()
                self.shared.append(bool(np.shares_memory(np.asarray(corrections),original)))
                del original
                return super().__setitem__(selected,corrections)
        store=RecordingRows(data.shape[0],data.shape[1]-1,data.shape[2],64,cutoff=1e-14)
        store.shared=[]
        return store

    def dense(self, store):
        output=np.zeros(store.shape)
        for i,levels in enumerate(store.rows):
            for k,row in enumerate(levels):
                if isinstance(store,CompactCorrectionArray):
                    native=row.native()
                    output[i,k]=np.asarray(native.toarray() if hasattr(native,'toarray') else native).reshape(-1)
                else:
                    values,columns=row
                    output[i,k,columns]=values
        return output

    def test_all_present_borrows_block_then_releases_owner(self):
        for cls in [SparseCorrectionArray,CompactCorrectionArray]:
            with self.subTest(storage=cls.__name__):
                data=self.data()
                returned=weakref.ref(data);owner=weakref.ref(data.base)
                store=self.collector(cls,data)
                selected=np.arange(len(data))[::-1]
                expected=np.zeros_like(data)
                expected[selected]=np.where(np.abs(data)>1e-14,data,0)
                store.set_rows(selected,data)
                self.assertEqual(store.shared,[True])
                np.testing.assert_array_equal(self.dense(store),expected)
                del data
                gc.collect()
                self.assertIsNone(returned())
                self.assertIsNone(owner())

    def test_partial_and_empty_selections_borrow_individual_rows(self):
        for cls in [SparseCorrectionArray,CompactCorrectionArray]:
            for empty in [False,True]:
                with self.subTest(storage=cls.__name__,empty=empty):
                    data=self.data();store=self.collector(cls,data)
                    selected=np.arange(len(data))[::-1]
                    selected[::2]=-1
                    if empty:selected[:]=-1
                    expected=np.zeros_like(data)
                    for source,destination in enumerate(selected):
                        if destination>=0:expected[destination]=np.where(np.abs(data[source])>1e-14,data[source],0)
                    store.set_rows(selected,data)
                    self.assertEqual(store.shared,[True]*int(np.count_nonzero(selected>=0)))
                    np.testing.assert_array_equal(self.dense(store),expected)

    def group(self):
        states = self.fixture['source_group_states']
        return BasisStateSpace(HarmonicOscillatorProductBasis(len(states[0])), states)

    def test_corrections_owns_internal_storage_and_factory(self):
        module=importlib.import_module('Psience.VPT2.Corrections')
        self.assertFalse((Path(module.__file__).parent/'CorrectionStorage.py').exists())
        for kind in [SparseCorrectionArray,CompactCorrectionArray]:
            self.assertEqual(kind.__module__,module.__name__)
            self.assertTrue(kind.__qualname__.startswith('PerturbationTheoryCorrections._'))
        for mode in ['dense','sparse','compact']:
            nstates,norders,nbasis=self.coefficients.shape
            packed=PerturbationTheoryCorrections.create_correction_array(
                mode,nstates,norders-1,nbasis,64)
            self.assertEqual(packed.shape,self.coefficients.shape)
            packed[:]=self.coefficients
            if mode=='dense':np.testing.assert_array_equal(packed,self.coefficients)
            else:np.testing.assert_array_equal(self.dense(packed),self.coefficients)

    def test_degeneracy_specs_has_no_solver_back_import(self):
        module = importlib.import_module('Psience.VPT2.DegeneracySpecs')
        tree = ast.parse(Path(module.__file__).read_text())
        imports = [node for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
        self.assertFalse(any(node.level == 1 and node.module == 'Solver' for node in imports))

    def test_injected_filter_preserves_martin_context(self):
        # Check configuration dispatch using recorded Pyrene states/coefficients;
        # this does not recalculate a molecular wavefunction.
        group = self.group()
        spec = MartinTestDegeneracySpec()
        spec._states = spec._basis = group
        spec._matrix = self.coefficients[:, 1, :]
        calls = []
        def record(selected, **options):
            calls.append(options)
            return selected
        selected = DegenerateMultiStateSpace.construct_filer(
            spec=spec, filter_func=record)(group)
        self.assertIs(selected, group)
        self.assertEqual(len(calls), 1)
        self.assertIs(calls[0]['states'], group)
        self.assertIs(calls[0]['basis'], group)
        self.assertIs(calls[0]['corrections'], spec._matrix)
        self.assertEqual(calls[0]['threshold'], spec.threshold)
        self.assertEqual(calls[0]['threshold_step_size'], 1 / 219475)

    def test_custom_and_unfiltered_specs_bypass_default_implementation(self):
        group = self.group()
        def unexpected(*args, **kwargs):
            self.fail('Explicit custom/unfiltered spec reached default implementation')
        calls = []
        def custom(selected):
            calls.append(selected)
            return selected
        for override in [custom, 'unfiltered']:
            spec = GroupsDegeneracySpec([group.excitations], group_filter=override)
            result = DegenerateMultiStateSpace.from_spec(
                spec, states=group,
                group_filter=spec.get_degenerate_group_filter(),
                default_group_filter=unexpected)
            np.testing.assert_array_equal(result.to_single().indices, group.indices)
        self.assertEqual(len(calls), 1)

    def test_pre_solver_setup_supplies_bounded_filter(self):
        group = self.group()
        spec = dict(groups=[group.excitations], group_filter='default')
        # Explicit groups need no frequency-based test for this dispatch check.
        calls = []
        def record(selected, **options):
            calls.append(options)
            return selected
        with patch.object(PerturbationTheorySolver, 'filter_degenerate_group', side_effect=record):
            _, groups = VPTStateSpace.build_degenerate_state_spaces(
                object.__new__(VPTStateSpace), spec, group, freqs=np.empty(0))
        self.assertEqual(len(calls), 1)
        np.testing.assert_array_equal(groups[0], group.excitations)

    def test_classic_solver_selects_bounded_or_dense_filter(self):
        group = self.group()
        spec = GroupsDegeneracySpec([group.excitations], group_filter='default')
        # Only exercise group construction, without representations or VPT.
        class ContextSolver(PerturbationTheorySolver):
            @property
            def zero_order_energies(self): return None
            @property
            def high_frequency_modes(self): return None
        calls = []
        def record(selected, **options):
            calls.append(options)
            return selected
        for constrained in [True, False]:
            solver = object.__new__(ContextSolver)
            solver._deg_states = None
            solver.degenerate_states = None
            solver.degeneracy_handlers = [spec]
            solver.states = group
            solver.full_basis = None
            solver.logger = NullLogger()
            solver.resonance_memory_constrained = constrained
            solver.resonance_block_size = 17
            calls.clear()
            with patch.object(ContextSolver, 'filter_degenerate_group', side_effect=record):
                solver.degenerate_spaces
            self.assertEqual(len(calls), int(constrained))
            if constrained:
                self.assertEqual(calls[0]['block_size'], 17)


    def test_actual_Pyrene_mixed_rows_pack_without_dense_block(self):
        import scipy.sparse as sp
        from pathlib import Path
        from Psience.VPT2.Corrections import PerturbationTheoryCorrections
        data=json.loads((Path(__file__).parent/'TestData/pyrene_vpt2_block_handoff.json').read_text())
        values=np.asarray(data['coefficients'],dtype=float)
        for mode in ['sparse','compact']:
            expected=PerturbationTheoryCorrections.create_correction_array(mode,len(values),values.shape[1]-1,values.shape[2],127,cutoff=1e-14)
            actual=PerturbationTheoryCorrections.create_correction_array(mode,len(values),values.shape[1]-1,values.shape[2],127,cutoff=1e-14)
            expected.set_rows(np.arange(len(values)),values)
            levels=tuple(sp.csr_matrix(values[:,k]) if k<2 else values[:,k] for k in range(values.shape[1]))
            actual.set_rows(np.arange(len(values)),levels)
            for i in range(len(values)):
                for k in range(values.shape[1]):
                    if mode=='compact':
                        a,b=actual.rows[i][k],expected.rows[i][k]
                        np.testing.assert_array_equal(a.values,b.values)
                        np.testing.assert_array_equal(np.asarray(a.columns),np.asarray(b.columns))
                    else:
                        for a,b in zip(actual.rows[i][k],expected.rows[i][k]):np.testing.assert_array_equal(a,b)



class VPTWorkspaceTests(unittest.TestCase):
    groups={'groups':[[[2,0,0],[0,1,0]],[[0,2,0],[0,0,2]]]}

    @classmethod
    def setUpClass(cls):
        cls.input=Path(__file__).parent/'TestData/HOH_freq.fchk'
        if not cls.input.is_file():
            cls.input=Path(__file__).resolve().parents[2]/'ci/tests/TestData/HOH_freq.fchk'
        if not cls.input.is_file():raise FileNotFoundError('Actual water FChk fixture is required')

    def run_water(self,order,memory,**options):
        molecule=Molecule.from_file(str(self.input));molecule.internals=None
        return VPTRunner.run_simple(molecule,2,order=order,expansion_order=2,
            calculate_intensities=order==2,degeneracy_specs=self.groups,
            evaluator_defaults='new',memory_constrained=memory,
            wavefunction_workspace_bytes=4096,correction_storage='dense',logger=False,**options)

    def assert_same_corrections(self,left,right):
        np.testing.assert_array_equal(left.corrs.states.excitations,right.corrs.states.excitations)
        np.testing.assert_array_equal(left.corrs.total_basis.indices,right.corrs.total_basis.indices)
        np.testing.assert_allclose(left.energies,right.energies,rtol=2e-11,atol=2e-12)
        for a,b in zip(left.corrs.wfn_corrections,right.corrs.wfn_corrections):
            np.testing.assert_allclose(a.asarray(),b.asarray(),rtol=2e-11,atol=2e-12)
        for a,b in zip(left.corrs.degenerate_hamiltonians,right.corrs.degenerate_hamiltonians):
            np.testing.assert_allclose(a,b,rtol=2e-11,atol=2e-12)

    def test_arbitrary_order_and_intermediate_normalization(self):
        # Higher orders use the actual quartic-truncated Hamiltonian; they do
        # not purport to include a full higher-order molecular force field.
        for order,intermediate in [(2,False),(4,False),(6,False),(4,True)]:
            with self.subTest(order=order,intermediate=intermediate):
                a=self.run_water(order,False,intermediate_normalization=intermediate)
                b=self.run_water(order,True,intermediate_normalization=intermediate)
                self.assert_same_corrections(a,b)
                if order==2:
                    ai,bi=a.intensities,b.intensities
                    if isinstance(ai,dict):
                        self.assertEqual(set(ai),set(bi))
                        for key in ai:np.testing.assert_allclose(ai[key],bi[key],rtol=2e-11,atol=2e-12)
                    else:np.testing.assert_allclose(ai,bi,rtol=2e-11,atol=2e-12)

    def test_real_operator_accumulation_and_projection(self):
        captured=[];method=PerturbationTheorySolver.apply_VPT_block_equations
        def capture(solver,group,*args,**kwargs):
            if not captured:captured.append((solver,group,kwargs.get('perturbations') or solver.representations))
            return method(solver,group,*args,**kwargs)
        PerturbationTheorySolver.apply_VPT_block_equations=capture
        try:wave=self.run_water(2,False)
        finally:PerturbationTheorySolver.apply_VPT_block_equations=method
        solver,group,H=captured[0]
        positions=wave.corrs.states.find(group)
        b=np.ascontiguousarray(wave.corrs.wfn_corrections[1][positions,:].asarray()).T
        target=np.ascontiguousarray(wave.corrs.wfn_corrections[2][positions,:].asarray()).T
        expected=target-H[1].dot(b)
        rows=solver.flat_total_space.find(group);projection=np.empty((len(rows),b.shape[1]))
        chunks=[]
        def inspect(rs,cs,product):
            chunks.append(product.nbytes)
            selected=np.flatnonzero((rows>=rs.start)&(rows<rs.stop))
            projection[selected,cs]=product[rows[selected]-rs.start]
        returned=matrix_dot_into(H[1],b,target,alpha=-1,beta=1,working_memory=4096,block_callback=inspect)
        self.assertIs(returned,target)
        self.assertGreater(len(chunks),1)
        self.assertLessEqual(max(chunks),4096)
        np.testing.assert_allclose(target,expected,rtol=2e-11,atol=2e-12)
        np.testing.assert_allclose(projection,H[1].dot(b)[rows],rtol=2e-11,atol=2e-12)
        before=b.copy()
        with self.assertRaises(ValueError):matrix_dot_into(H[1],b,b)
        np.testing.assert_array_equal(before,b)

    def test_unprojected_near_zero_error_unchanged(self):
        captured=[];method=PerturbationTheorySolver.apply_VPT_block_equations
        def capture(solver,group,*args,**kwargs):
            if not captured:captured.append((solver,group,kwargs.get('perturbations') or solver.representations))
            return method(solver,group,*args,**kwargs)
        PerturbationTheorySolver.apply_VPT_block_equations=capture
        try:self.run_water(2,False)
        finally:PerturbationTheorySolver.apply_VPT_block_equations=method
        solver,group,H=captured[0];inds=solver.flat_total_space.find(group)
        gaps=np.abs(solver.zero_order_energies-solver.zero_order_energies[inds[0]])
        gaps[inds]=np.inf;cutoff=float(np.min(gaps)*1.01+np.finfo(float).eps)
        errors=[]
        for flag in [False,True]:
            solver.memory_constrained=flag
            with self.assertRaises(ValueError) as error:
                solver.apply_VPT_block_equations(group,perturbations=H,non_zero_cutoff=cutoff)
            errors.append(str(error.exception))
        self.assertEqual(*errors)


if __name__=='__main__': unittest.main()
