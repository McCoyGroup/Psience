from pathlib import Path
import tempfile,unittest
import numpy as np
import scipy.sparse as sp
from unittest.mock import patch
from McUtils.Numputils import external_unique,save_sparse_matrix,load_sparse_matrix
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace,BraKetSpace,SelectionRuleStateSpace

class ExternalTests(unittest.TestCase):
    def test_stable_and_sorted_unique_rows(self):
        rng=np.random.default_rng(721)
        for width in [1,2,4]:
            x=rng.integers(-10,11,size=(117,width))
            if width==1:x=x[:,0]
            unique,first=np.unique(x,axis=0,return_index=True)
            for sorting in [True,False]:
                for chunk in [1,7,1000]:
                    out,where=external_unique([x[:13],x[13:60],x[60:]],sort=sorting,
                        return_index=True,chunk_size=chunk)
                    order=np.arange(len(unique)) if sorting else np.argsort(first)
                    np.testing.assert_array_equal(out,unique[order])
                    np.testing.assert_array_equal(where,first[order])
                    if chunk==7:
                        values,where=external_unique(x,sort=sorting,return_index=True,
                            return_values=False,chunk_size=chunk)
                        self.assertIsNone(values)
                        np.testing.assert_array_equal(where,first[order])

    def test_empty_and_file_output(self):
        np.testing.assert_array_equal(external_unique([np.empty(0,dtype=int)]),[])
        with tempfile.TemporaryDirectory() as temp:
            x=external_unique(np.array([7,3,7,2,1]),sort=False,output_file=Path(temp)/'ids.npy',chunk_size=2)
            self.assertIsInstance(x,np.memmap)
            np.testing.assert_array_equal(x,[7,3,2,1])

    def test_pair_order_and_orientation(self):
        rng=np.random.default_rng(831);basis=HarmonicOscillatorProductBasis(3)
        for size in [0,1,50,250]:
            x=rng.integers(0,15,size=(2,size))
            pairs=BraKetSpace(BasisStateSpace(basis,x[0],mode='indices'),
                              BasisStateSpace(basis,x[1],mode='indices'))
            memory=pairs.remove_duplicates()
            external=pairs.remove_duplicates(method='external',chunk_size=13)
            np.testing.assert_array_equal(memory.bras.indices,external.bras.indices)
            np.testing.assert_array_equal(memory.kets.indices,external.kets.indices)

    def test_streamed_nested_indices(self):
        basis=HarmonicOscillatorProductBasis(2)
        base=BasisStateSpace(basis,[0,1],mode='indices')
        blocks=[BasisStateSpace(basis,v,mode='indices') for v in [[2,3],[1,4]]]
        space=SelectionRuleStateSpace(base,blocks)
        np.testing.assert_array_equal(np.concatenate(list(space.iter_index_blocks())),space.indices)

    def test_mapped_sparse_product(self):
        a=sp.diags([np.arange(1,100),np.full(98,.7)],[0,1],format='csr')
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'csr.npz';save_sparse_matrix(path,a);mapped=load_sparse_matrix(path)
            for name in ['data','indices','indptr']:
                self.assertTrue(np.shares_memory(getattr(mapped.data,name),mapped._mapped_archive[name]))
            x=np.arange(a.shape[1],dtype=float)
            np.testing.assert_array_equal(mapped.dot(x),a.dot(x))
            # Replacing a cache file must not mutate an already mapped matrix.
            save_sparse_matrix(path,2*a)
            np.testing.assert_array_equal(mapped.dot(x),a.dot(x))
            np.testing.assert_array_equal(load_sparse_matrix(path).dot(x),2*a.dot(x))
            mapped._mapped_archive.close()
            del mapped

    def test_temporary_maps_closed_before_cleanup(self):
        import McUtils.Numputils.ExternalMemory as module
        original=module._close_mapping;seen=[]
        def close(array):seen.append(array);original(array)
        with patch.object(module,'_close_mapping',new=close):
            out=external_unique(np.array([7,3,2,7,4,1]),sort=False,chunk_size=2)
        np.testing.assert_array_equal(out,[7,3,2,4,1])
        self.assertTrue(seen)
        self.assertTrue(all(a._mmap.closed for a in seen))

    def test_numpy_only_merge_fallback(self):
        import McUtils.Numputils.ExternalMemory as module
        with patch.object(module,'_merge_function',return_value=module._merge):
            actual=external_unique(np.array([9,1,3,9,2,1]),sort=False,chunk_size=2)
        np.testing.assert_array_equal(actual,[9,1,3,2])

if __name__=='__main__':unittest.main()


class EvaluatedRepresentationTests(unittest.TestCase):
    """Actual quartic-water operator tests; all molecular evaluation is on HPC."""
    @staticmethod
    def _fixture():
        from Psience.VPT2.Hamiltonian import PerturbationTheoryHamiltonian
        path = Path(__file__).parent/'TestData/HOH_freq.fchk'
        ham = PerturbationTheoryHamiltonian.from_fchk(str(path), logger=False,
                                                    operator_contraction_method='direct')
        return ham, BasisStateSpace.from_quanta(ham.basis, range(3))

    def test_molecular_resident_and_deferred_modes(self):
        from Psience.BasisReps.ExternalStateSpaces import StateSpaceDiskStore, DiskSelectionRuleStateSpace
        from McUtils.Numputils import EvaluatedSparseMatrix, matrix_dot_into
        ham, states = self._fixture(); rep = ham.get_perturbations(2)[2]
        resident = rep.get_transformed_space(states)
        for external in [False, True]:
            with tempfile.TemporaryDirectory() as temp:
                store = StateSpaceDiskStore(temp, block_size=23)
                try:
                    graph = (DiskSelectionRuleStateSpace(store.coerce(states), rep.selection_rules)
                             if external else resident)
                    total = (store.union_spaces([store.coerce(states), graph.to_single()], sort=True)
                             if external else states.union(graph.to_single()).take_unique(sort=True))
                    stored = rep.get_representation_matrix(graph, total, zero_element_warning=False)
                    expected = stored.asarray()
                    for mode in ['positions', 'runtime']:
                        matrix = rep.get_representation_matrix(graph, total, evaluation_mode=mode,
                            zero_element_warning=False, evaluation_block_size=19, evaluation_working_memory=2**20)
                        self.assertIsInstance(matrix, EvaluatedSparseMatrix)
                        self.assertEqual(matrix.evaluation_mode, mode)
                        if mode == 'runtime': self.assertEqual(matrix.nbytes, 0)
                        np.testing.assert_array_equal(matrix.asarray(), expected)
                        rhs = np.zeros((len(total), 2)); rhs[[0, 2, 7], :] = [[1, .1], [.3, 0], [0, -.2]]
                        np.testing.assert_allclose(matrix.dot(rhs), expected.dot(rhs), rtol=2e-13, atol=2e-14)
                        np.testing.assert_array_equal(matrix[[2, 0, 2], [4, 1, 4]], expected[[2, 0, 2], [4, 1, 4]])
                finally: store.close()

    def test_deferred_endpoint_restrictions_and_postfilter_semantics(self):
        from Psience.BasisReps.ExternalStateSpaces import StateSpaceDiskStore, DiskSelectionRuleStateSpace
        ham, states = self._fixture(); rep = ham.get_perturbations(2)[2]
        with tempfile.TemporaryDirectory() as temp:
            store = StateSpaceDiskStore(temp, block_size=11)
            try:
                source = store.coerce(states)
                # A content-dependent postfilter must see the original batches.
                def postfilter(space): return space.take_subspace(slice(0, len(space)//2))
                recipes = [(source, rep.selection_rules, None, None, (postfilter,), ())]
                graph = DiskSelectionRuleStateSpace(source, recipes=recipes)
                full = np.concatenate(list(graph.iter_pair_blocks()))
                selected = source.take_subspace([0, 2])
                actual = np.concatenate(list(graph.iter_pair_blocks(source_filter=selected)))
                expected = full[np.isin(full[:, 0], selected.indices)]
                np.testing.assert_array_equal(actual, expected)
                # Explicit owned nested records must filter actual bras, not
                # confuse the outer representative/owner with the edge source.
                child = SelectionRuleStateSpace(states.take_subspace([1]), [states.take_subspace([2, 3])])
                nested = SelectionRuleStateSpace(states.take_subspace([0]), [child])
                owned = DiskSelectionRuleStateSpace.from_resident(nested, store)
                full = np.concatenate(list(owned.iter_pair_blocks()))
                chosen = states.take_subspace([1])
                actual = np.concatenate(list(owned.iter_pair_blocks(source_filter=chosen)))
                np.testing.assert_array_equal(actual, full[np.isin(full[:, 0], chosen.indices)])
            finally: store.close()

    def test_molecular_h0_vector_and_policy_options(self):
        from Psience.VPT2 import VPTSolverOptions
        from McUtils.Numputils import EvaluatedSparseMatrix
        ham, states = self._fixture(); rep = ham.get_perturbations(0)[0]
        stored = rep.get_representation_matrix(states, states, diagonal=True, zero_element_warning=False)
        vector = rep.get_representation_matrix(states, states, diagonal=True, evaluation_mode='diagonal', zero_element_warning=False)
        implicit = rep.get_representation_matrix(states, states, diagonal=True, evaluation_mode='runtime', zero_element_warning=False)
        np.testing.assert_array_equal(vector.diag, stored.diag)
        np.testing.assert_array_equal(implicit.diag, stored.diag)
        self.assertIsNone(vector.positions)
        self.assertEqual(vector.nbytes, vector.diag.nbytes)
        self.assertEqual(implicit.nbytes, 0)
        opts = VPTSolverOptions(evaluation_mode={1: 'stored', 2: 'runtime'})
        self.assertEqual(opts.opts['evaluation_mode'], {1: 'stored', 2: 'runtime'})

    def test_deferred_directed_restricted_modes_and_noncanonical_basis(self):
        from Psience.BasisReps.ExternalStateSpaces import StateSpaceDiskStore, DiskSelectionRuleStateSpace, DiskBasisStateSpace
        ham, states = self._fixture(); rep = ham.get_perturbations(2)[2]
        with tempfile.TemporaryDirectory() as temp:
            store = StateSpaceDiskStore(temp, block_size=31)
            try:
                source = store.coerce(states)
                # Deliberately directed rules, a mode subset, and a signed
                # exclusion must keep the stored symmetric graph semantics.
                rules = ((2,), (-2,), ())
                recipes = [(source, rules, None, (1,), (), ((0, 2, 0),))]
                graph = DiskSelectionRuleStateSpace(source, recipes=recipes)
                canonical = store.union_spaces([source, graph.to_single()], sort=True)
                total = DiskBasisStateSpace(ham.basis, canonical.indices[::-1], store, canonical=False, unique=True)
                target_filter = source.take_subspace([0, 2, 4])
                expected = rep.get_representation_matrix(graph, total, filter_space=target_filter,
                                                        zero_element_warning=False).asarray()
                for mode in ['positions', 'runtime']:
                    actual = rep.get_representation_matrix(graph, total, filter_space=target_filter,
                        evaluation_mode=mode, evaluation_block_size=13, zero_element_warning=False)
                    np.testing.assert_array_equal(actual.asarray(), expected)
            finally: store.close()

    def test_bounded_runtime_avoids_forced_global_collection(self):
        import Psience.BasisReps.Representations as module
        ham, states = self._fixture(); rep = ham.get_perturbations(2)[2]
        graph = rep.get_transformed_space(states)
        total = states.union(graph.to_single()).take_unique(sort=True)
        runtime = rep.get_representation_matrix(graph, total, evaluation_mode='runtime', zero_element_warning=False)
        pairs = BraKetSpace(total.take_subspace([0]), total.take_subspace([0]))
        with patch.object(module.gc, 'collect', return_value=0) as collect:
            actual = runtime[0, 0]
            self.assertEqual(collect.call_count, 0)
            expected = rep.get_brakets(pairs, check_orthogonality=True, memory_constrained=True)
            self.assertGreater(collect.call_count, 0)
        np.testing.assert_array_equal(actual, np.asarray(expected)[0])

    def test_runtime_coordinates_outlive_temporary_query_maps(self):
        from Psience.BasisReps.ExternalStateSpaces import StateSpaceDiskStore, DiskSelectionRuleStateSpace
        ham, states = self._fixture(); rep = ham.get_perturbations(2)[2]
        with tempfile.TemporaryDirectory() as temp:
            store = StateSpaceDiskStore(temp, block_size=17)
            try:
                source = store.coerce(states)
                graph = DiskSelectionRuleStateSpace(source, rep.selection_rules)
                total = store.union_spaces([source, graph.to_single()], sort=True)
                runtime = rep.get_representation_matrix(graph, total, evaluation_mode='runtime',
                    evaluation_block_size=7, zero_element_warning=False)
                expected = rep.get_representation_matrix(graph, total, zero_element_warning=False).asarray()
                blocks = list(runtime.iter_evaluated_blocks(block_size=7))
                self.assertFalse(list(store.directory.glob('operator_pairs_*')))
                # Every query map is now closed/unlinked. Retained coordinates
                # must own ordinary bounded arrays, not dangling mapped views.
                actual = np.zeros_like(expected)
                for values, indices in blocks: actual[indices] = values
                np.testing.assert_array_equal(actual, expected)
            finally: store.close()

    def test_zero_perturbation_modes(self):
        ham, states = self._fixture(); rep = 0*ham.get_perturbations(2)[2]
        metadata_before = rep.is_zero
        graph = ham.get_perturbations(2)[2].get_transformed_space(states)
        total = states.union(graph.to_single()).take_unique(sort=True)
        for mode in ['positions', 'runtime']:
            matrix = rep.get_representation_matrix(graph, total, evaluation_mode=mode)
            np.testing.assert_array_equal(matrix.dot(np.ones(len(total))), np.zeros(len(total)))
            self.assertEqual(matrix[0, 0], 0.)
            self.assertEqual(rep.is_zero, metadata_before)
            if mode == 'runtime': self.assertEqual(matrix.nbytes, 0)

    def test_default_options_leave_the_evaluation_policy_absent(self):
        from Psience.VPT2 import VPTSolverOptions
        self.assertNotIn('evaluation_mode', VPTSolverOptions().opts)
        self.assertNotIn('evaluation_mode', VPTSolverOptions(evaluation_mode=None).opts)
        ham, states = self._fixture()
        solver = ham.get_solver(states, order=2, expansion_order=2)
        self.assertIsNone(solver.evaluation_mode)
        self.assertNotIn('evaluation_mode', vars(solver))
        explicit = ham.get_solver(states, order=2, expansion_order=2, evaluation_mode={2:'runtime'})
        self.assertEqual(explicit.evaluation_mode, {2:'runtime'})
        self.assertNotIn(0, explicit.evaluation_mode)
        with self.assertRaises(ValueError):
            ham.get_solver(states, order=2, expansion_order=2, evaluation_mode='runtime',
                           use_cached_representations=True)

    def test_default_representation_cannot_enter_the_evaluated_backend(self):
        from McUtils.Numputils import EvaluatedSparseMatrix
        ham, states = self._fixture(); rep = ham.get_perturbations(2)[2]
        graph = rep.get_transformed_space(states)
        total = states.union(graph.to_single()).take_unique(sort=True)
        expected = rep.get_representation_matrix(graph, total, zero_element_warning=False).asarray()
        with patch.object(rep, 'get_evaluated_representation', side_effect=AssertionError('Implicit backend is opt-in')), \
             patch.object(EvaluatedSparseMatrix, '__init__', side_effect=AssertionError('Implicit matrix is opt-in')):
            for options in [{}, {'evaluation_mode':None}, {'evaluation_mode':'stored'}]:
                actual = rep.get_representation_matrix(graph, total, zero_element_warning=False, **options)
                self.assertNotIsInstance(actual, EvaluatedSparseMatrix)
                np.testing.assert_array_equal(actual.asarray(), expected)

    def test_default_cached_h0_is_returned_without_conversion(self):
        from McUtils.Numputils import EvaluatedSparseMatrix
        ham, states = self._fixture(); h0 = ham.get_perturbations(0)[0]
        legacy = h0.get_representation_matrix(states, states, diagonal=True, zero_element_warning=False)
        solver = ham.get_solver(states, order=2, expansion_order=2)
        class Cache(dict):
            def __enter__(self): return self
            def __exit__(self, *args): return False
        cached = [legacy]
        solver.checkpointer = Cache(representations=cached)
        solver.use_cached_representations = True
        with patch.object(EvaluatedSparseMatrix, '__init__', side_effect=AssertionError('Cached H0 must stay legacy')):
            self.assertIs(solver.get_VPT_representations(), cached)
            self.assertIs(solver.get_VPT_representations()[0], legacy)
