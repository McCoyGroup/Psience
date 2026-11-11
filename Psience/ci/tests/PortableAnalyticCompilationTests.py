"""Compiled artifacts must execute without expression planning or construction."""
import json
import pickle
from pathlib import Path
import shutil
import subprocess
import sys
from unittest import mock

import numpy as np
import pytest
from McUtils.Zachary import Symbols

from Psience.VPT2 import Analytic as a
from Psience.VPT2 import Evaluators as ce
from Psience.VPT2 import AnalyticVPTRunner


def simple_run_inputs(nmodes=1, scale=1.):
    """Explicit synthetic surfaces exercise the real public system constructor."""
    frequencies = .024*np.power(np.arange(1, nmodes+1), .71)
    system = [['H', 'H', 'H'], np.array([[0., 0., 0.], [1.4, 0., 0.], [0., 1.4, 0.]])]
    states = [np.zeros(nmodes, dtype=int), np.eye(nmodes, dtype=int)[0]]
    options = dict(order=2, logger=False, clear_caches=True, force_return_on_crash=False,
        modes=dict(matrix=np.eye(9, nmodes), inverse=np.eye(nmodes, 9), freqs=frequencies),
        potential_terms=[np.diag(frequencies), np.ones((nmodes,)*3)*.0003*scale,
                         np.ones((nmodes,)*4)*.00001*scale],
        kinetic_terms=[np.diag(frequencies), 0, 0], coriolis_terms=[0, 0, 0],
        pseudopotential_terms=[0, 0, 0],
        dipole_terms=[[0., np.ones(nmodes)*factor, np.ones((nmodes,)*2)*.01,
                       np.zeros((nmodes,)*3)] for factor in (1., .3, .1)])
    return system, states, options


@pytest.mark.parametrize('backend', ['portable', 'pickle'])
def test_run_simple_builds_and_reuses_single_libraries_across_systems(tmp_path, backend):
    system, states, options = simple_run_inputs()
    reference = AnalyticVPTRunner.run_simple(system, states, compile=False, **options)
    ir_file, compiled_file = tmp_path/'expressions/library.ir.tar', tmp_path/'compiled/library.plan.tar'
    runner, actual = AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file,
        compiled_file=compiled_file, storage_backend=backend, return_runner=True, **options)
    for name in ('energy_corrections', 'transition_moment_corrections'):
        np.testing.assert_allclose(getattr(actual, name), getattr(reference, name), atol=2e-13, rtol=2e-11)
    assert ce.AnalyticExpressionIR.expression_names(ir_file)
    names = ce.CompiledAnalyticExpression.expression_names(compiled_file)
    assert names
    for name in names:
        with ce._CompiledAnalyticArchive.load(compiled_file, expression=name) as entry:
            assert entry.storage_backend == backend
            assert entry.manifest['specializations'] == dict(active_operators=None, mode_dimension=None)
    assert sorted(path.name for path in compiled_file.parent.iterdir()) == [compiled_file.name]
    assert not list(tmp_path.rglob('evaluations'))
    before = compiled_file.read_bytes()
    # Transfer only the compiled file. A different number of modes and changed
    # physical tensors must not construct expressions or re-plan saved programs.
    recipient = tmp_path/'recipient/library.plan.tar'
    recipient.parent.mkdir()
    shutil.copyfile(compiled_file, recipient)
    ir_file.unlink()
    system, states, options = simple_run_inputs(2, scale=.7)
    reference = AnalyticVPTRunner.run_simple(system, states, compile=False, **options)
    with mock.patch.object(a.PerturbationTheoryTerm, 'get_poly_terms', side_effect=AssertionError('symbolic rebuild')), \
         mock.patch.object(a._CompiledPolynomial, 'compile', side_effect=AssertionError('polynomial compile')), \
         mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan', side_effect=AssertionError('plan compile')), \
         mock.patch.object(ce, '_implementation_signature', side_effect=AssertionError('producer source gate')):
        warm_runner, actual = AnalyticVPTRunner.run_simple(system, states, compiled_file=recipient,
            compilation_options={'reuse_only':True}, return_runner=True, **options)
    for name in ('energy_corrections', 'transition_moment_corrections'):
        np.testing.assert_allclose(getattr(actual, name), getattr(reference, name), atol=2e-13, rtol=2e-11)
    assert warm_runner.eval.compilation_statistics()['compilations'] == 0
    assert recipient.read_bytes() == before
    assert sorted(path.name for path in recipient.parent.iterdir()) == [recipient.name]
    # Repeat in an independent interpreter, without inherited expression caches.
    script = '''
import json, sys
from pathlib import Path
from unittest import mock
import numpy as np
sys.path.insert(0, sys.argv[1])
from PortableAnalyticCompilationTests import simple_run_inputs
from Psience.VPT2 import AnalyticVPTRunner, Analytic as a, Evaluators as ce
system, states, options = simple_run_inputs(2, scale=.7)
with mock.patch.object(a.PerturbationTheoryTerm, 'get_poly_terms', side_effect=AssertionError('rebuild')), \\
     mock.patch.object(a._CompiledPolynomial, 'compile', side_effect=AssertionError('compile')), \\
     mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan', side_effect=AssertionError('planning')), \\
     mock.patch.object(ce, '_implementation_signature', side_effect=AssertionError('source gate')):
    result = AnalyticVPTRunner.run_simple(system, states, compiled_file=sys.argv[2],
        compilation_options={'reuse_only': True}, **options)
print(json.dumps({name:np.asarray(getattr(result,name)).tolist() for name in
                 ('energy_corrections','transition_moment_corrections')}))
'''
    child = subprocess.run([sys.executable, '-c', script, str(Path(__file__).parent), str(recipient)],
                           capture_output=True, text=True)
    assert child.returncode == 0, child.stderr
    for name, values in json.loads(child.stdout).items():
        np.testing.assert_allclose(values, getattr(reference,name), atol=2e-13, rtol=2e-11)


@pytest.mark.parametrize('input_backend,output_backend', [('portable','pickle'), ('pickle','portable')])
def test_run_simple_recompiles_retained_ir_without_symbolic_rebuild(tmp_path, input_backend, output_backend):
    system, states, options = simple_run_inputs()
    options['calculate_intensities'] = False
    ir_file, compiled_file = tmp_path/'ir/library.tar', tmp_path/'compiled/library.tar'
    reference = AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file,
        compiled_file=compiled_file, storage_backend=input_backend, **options)
    compiled_file.unlink()
    with mock.patch.object(a.PerturbationTheoryTerm, 'get_poly_terms', side_effect=AssertionError('symbolic rebuild')):
        actual = AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file,
            compiled_file=compiled_file, storage_backend=output_backend, **options)
    np.testing.assert_allclose(actual.energy_corrections, reference.energy_corrections, atol=2e-14)
    for name in ce.CompiledAnalyticExpression.expression_names(compiled_file):
        with ce._CompiledAnalyticArchive.load(compiled_file, expression=name) as entry:
            assert entry.storage_backend == output_backend


def test_run_simple_rejects_disabled_compilation_with_artifact_paths():
    with mock.patch.object(AnalyticVPTRunner, 'construct', side_effect=AssertionError('system construction')):
        with pytest.raises(ValueError, match='compile=True'):
            AnalyticVPTRunner.run_simple('unused', [], compile=False, compiled_file='unused.tar')


def test_run_simple_accepts_offline_compilation_of_its_retained_ir(tmp_path):
    system, states, options = simple_run_inputs()
    options['calculate_intensities'] = False
    ir_file, compiled_file = tmp_path/'ir/library.tar', tmp_path/'compiled/library.tar'
    reference = AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file,
        compiled_file=compiled_file, **options)
    compiled_file.unlink()
    with mock.patch.object(a.PerturbationTheoryTerm, 'get_poly_terms', side_effect=AssertionError('rebuild')):
        ce.AnalyticExpressionStore().compile_ir(ir_file, compiled_file, storage_backend='pickle')
        ir_file.unlink()
        actual = AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file,
            compilation_options={'reuse_only':True}, **options)
    np.testing.assert_array_equal(actual.energy_corrections, reference.energy_corrections)


def test_run_simple_extends_and_converts_libraries_preserving_existing_entries(tmp_path):
    system, states, options = simple_run_inputs()
    ir_file, compiled_file = tmp_path/'ir/library.tar', tmp_path/'compiled/library.tar'
    AnalyticVPTRunner.run_simple(system, states, calculate_intensities=False,
        ir_file=ir_file, compiled_file=compiled_file, **options)
    before = ce.CompiledAnalyticExpression.expression_names(compiled_file)
    reference = AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file,
        compiled_file=compiled_file, **options)
    after = ce.CompiledAnalyticExpression.expression_names(compiled_file)
    assert set(before) < set(after)
    with mock.patch.object(a.PerturbationTheoryTerm, 'get_poly_terms', side_effect=AssertionError('rebuild')), \
         mock.patch.object(a._CompiledPolynomial, 'compile', side_effect=AssertionError('compile')):
        actual = AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file, compiled_file=compiled_file,
            storage_backend='pickle', compilation_options={'reuse_only':True}, **options)
    assert ce.CompiledAnalyticExpression.expression_names(compiled_file) == after
    for name in ('energy_corrections', 'transition_moment_corrections'):
        np.testing.assert_array_equal(getattr(actual,name),getattr(reference,name))
    for reader, path in ((ce.AnalyticExpressionIR,ir_file),(ce._CompiledAnalyticArchive,compiled_file)):
        for name in reader.expression_names(path):
            with reader.load(path,expression=name) as entry: assert entry.storage_backend == 'pickle'


def test_run_simple_checkpoints_only_in_independent_evaluation_folder(tmp_path):
    system, states, options = simple_run_inputs()
    options['calculate_intensities'] = False
    compiled_file = tmp_path/'compiled/library.tar'
    evaluation_folder = tmp_path/'evaluations/run-001'
    reference = AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file,
        evaluation_folder=evaluation_folder, **options)
    assert list(evaluation_folder.rglob('*.pkl'))
    before = compiled_file.read_bytes()
    with mock.patch.object(a._CompiledPolynomial, 'evaluate_polynomial', side_effect=AssertionError('checkpoint missed')):
        actual = AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file,
            evaluation_folder=evaluation_folder, **options)
    np.testing.assert_array_equal(actual.energy_corrections,reference.energy_corrections)
    assert compiled_file.read_bytes() == before
    assert sorted(p.name for p in compiled_file.parent.iterdir()) == [compiled_file.name]
    with pytest.raises(ValueError, match='outside'):
        AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file,
            evaluation_folder=compiled_file.parent/'evaluations', **options)


def test_run_simple_failed_extension_leaves_published_files_unchanged(tmp_path):
    system, states, options = simple_run_inputs()
    ir_file, compiled_file = tmp_path/'ir/library.tar', tmp_path/'compiled/library.tar'
    AnalyticVPTRunner.run_simple(system, states, calculate_intensities=False,
        ir_file=ir_file, compiled_file=compiled_file, **options)
    before = ir_file.read_bytes(), compiled_file.read_bytes()
    with mock.patch.object(ce.AnalyticExpressionStore, '_lower_ir_part', side_effect=RuntimeError('interrupted compile')):
        with pytest.raises(RuntimeError, match='interrupted compile'):
            AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file,
                compiled_file=compiled_file, **options)
    assert (ir_file.read_bytes(),compiled_file.read_bytes()) == before


def test_run_simple_automatic_action_policy_uses_saved_full_programs(tmp_path):
    system, states, options = simple_run_inputs(7)
    options.update(calculate_intensities=False)
    reference = AnalyticVPTRunner.run_simple(system, states, compile=False, evaluation_mode='direct', **options)
    with mock.patch.object(ce.AnalyticCorrectionCompiler, 'evaluate_action_expansion',
                           side_effect=AssertionError('controller symbolic projection')):
        actual = AnalyticVPTRunner.run_simple(system, states,
            compiled_file=tmp_path/'compiled/library.tar', **options)
    np.testing.assert_allclose(actual.energy_corrections,reference.energy_corrections,atol=2e-14,rtol=2e-12)


@pytest.mark.parametrize('handling', ['projected', 'conditional', 'filtered'])
def test_run_simple_resonance_libraries_bind_membership_at_runtime(tmp_path, handling):
    system, _, options = simple_run_inputs(2)
    states = [[0,0], [1,0], [0,1]]
    options.update(degeneracy_specs={'groups':[[[1,0],[0,1]]]}, degeneracy_handling=handling,
                   evaluation_mode='direct')
    reference = AnalyticVPTRunner.run_simple(system, states, compile=False, **options)
    compiled_file = tmp_path/'compiled/library.tar'
    actual = AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file, **options)
    for name in ('energy_corrections', 'transition_moment_corrections'):
        np.testing.assert_allclose(getattr(actual,name),getattr(reference,name),atol=2e-12,rtol=2e-10)
    with mock.patch.object(a.PerturbationTheoryTerm, 'get_poly_terms', side_effect=AssertionError('rebuild')), \
         mock.patch.object(a._CompiledPolynomial, 'compile', side_effect=AssertionError('compile')):
        actual = AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file,
            compilation_options={'reuse_only':True}, **options)
    for name in ('energy_corrections', 'transition_moment_corrections'):
        np.testing.assert_allclose(getattr(actual,name),getattr(reference,name),atol=2e-12,rtol=2e-10)


def test_run_simple_saves_additional_operator_corrections(tmp_path):
    system, states, options = simple_run_inputs()
    options.update(calculate_intensities=False,
        operator_terms={'position':[np.ones(1), np.ones((1,)*2)*.01, np.zeros((1,)*3)]})
    reference = AnalyticVPTRunner.run_simple(system, states, compile=False, **options)
    compiled_file = tmp_path/'compiled/library.tar'
    actual = AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file, **options)
    np.testing.assert_allclose(actual.operator_corrections,reference.operator_corrections,atol=2e-13,rtol=2e-11)
    with mock.patch.object(a.PerturbationTheoryTerm, 'get_poly_terms', side_effect=AssertionError('rebuild')), \
         mock.patch.object(a._CompiledPolynomial, 'compile', side_effect=AssertionError('compile')):
        actual = AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file,
            compilation_options={'reuse_only':True}, **options)
    np.testing.assert_allclose(actual.operator_corrections,reference.operator_corrections,atol=2e-13,rtol=2e-11)


def test_artifact_configuration_rejects_overlapping_paths_and_missing_requests(tmp_path):
    store = ce.AnalyticExpressionStore()
    with pytest.raises(ValueError, match='separate'):
        with store.artifact_session(ir_file=tmp_path/'same.tar', compiled_file=tmp_path/'same.tar'): pass
    with pytest.raises(ValueError, match='requires'):
        with store.artifact_session(ir_file=tmp_path/'source.tar'): pass
    with pytest.raises(ValueError, match='compatible expression'):
        with store.artifact_session(compiled_file=tmp_path/'compiled/library.tar',reuse_only=True):
            store.compile_factory(lambda: (_ for _ in ()).throw(AssertionError('unexpected build')),context='missing')
    assert not (tmp_path/'compiled/library.tar').exists()


def test_run_simple_can_add_retained_ir_to_an_existing_compiled_library(tmp_path):
    system, states, options = simple_run_inputs()
    options['calculate_intensities'] = False
    ir_file, compiled_file = tmp_path/'ir/library.tar', tmp_path/'compiled/library.tar'
    reference = AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file, **options)
    before = compiled_file.read_bytes()
    with mock.patch.object(a._CompiledPolynomial, 'compile', side_effect=AssertionError('unnecessary recompile')):
        actual = AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file,
            compiled_file=compiled_file, **options)
    assert ce.AnalyticExpressionIR.expression_names(ir_file)
    assert compiled_file.read_bytes() == before
    np.testing.assert_array_equal(actual.energy_corrections,reference.energy_corrections)


def test_run_simple_reuse_only_reports_missing_requested_ir(tmp_path):
    system, states, options = simple_run_inputs()
    options['calculate_intensities'] = False
    ir_file, compiled_file = tmp_path/'ir/library.tar', tmp_path/'compiled/library.tar'
    AnalyticVPTRunner.run_simple(system, states, compiled_file=compiled_file, **options)
    before = compiled_file.read_bytes()
    with pytest.raises(ValueError, match='mathematical source'):
        AnalyticVPTRunner.run_simple(system, states, ir_file=ir_file, compiled_file=compiled_file,
            compilation_options={'reuse_only':True}, **options)
    assert not ir_file.exists()
    assert compiled_file.read_bytes() == before


def expression():
    return a.AnalyticPerturbationTheorySolver.from_order(2, internals=False, logger=False).energy_correction(2)([])


def inputs(nmodes, seed=97):
    rng = np.random.default_rng(seed)
    frequencies = .017*np.sqrt(np.arange(2, nmodes+2))
    coefficients = [[np.diag(frequencies)/2, np.diag(frequencies)/2],
                    [rng.normal(size=(nmodes,)*3)*1.e-4, 0],
                    [rng.normal(size=(nmodes,)*4)*1.e-5, 0, 0, 0, 0]]
    states = [(np.zeros(nmodes, dtype=int), np.arange(nmodes)[None]),
              (np.arange(nmodes) % 3, np.array([np.arange(nmodes), np.arange(nmodes)[::-1]]))]
    return states, [coefficients], frequencies


def resumable_expression():
    key = a.PTTensorCoeffProductSum.canonical_key(((1, 0, 0, 0, 1),))
    polynomial = a.ProductPTPolynomial([np.array([1., .3]), np.array([2., .5])], steps=0)
    value = a.PTEnergyChangeProductSum({((0, 1, -2),): polynomial, (): polynomial.scale(.01)})
    poly = a.PTTensorCoeffProductSum({key: value}, canonicalize=False)
    return a.PerturbationTheoryExpressionEvaluator(None, poly, change=(), logger=False)


def test_polynomial_program_matches_native_transformations():
    leaf = a.PolyPath.from_coeffs([np.array([1., 2., -.3]), np.array([2., -1.])], steps=3)
    polynomials = [leaf, leaf.shift((1, -2)), leaf.permute((1, 0)),
                   leaf*leaf.shift((1, 1)), (leaf+leaf.shift((0, 1)))*leaf,
                   leaf.permutation_sum([(0, 1), (1, 0)])]
    states = np.array([[0, 0], [1, 3], [2, 4], [7, 1]])
    for poly in polynomials:
        expected = poly.evaluate_polynomial(states)
        compiled = a._CompiledPolynomial.compile(poly)
        with mock.patch.object(a.PolyAxis, 'coeffs', new_callable=mock.PropertyMock,
                               side_effect=AssertionError('factor convolution at runtime')):
            np.testing.assert_array_equal(compiled.evaluate_polynomial(states), expected)


def test_deep_shared_program_roundtrip_needs_no_recursive_symbolic_execution(tmp_path):
    common = a.PolyPath.from_coeffs([np.array([1., 2.])])
    deep = common
    for _ in range(1500):
        deep = a.PolyPath((), node=('scale', (deep, 1.001)), ndim=1, order=(1,))
    cache = {}
    first = a._CompiledPolynomial.compile(common, cache=cache)
    second = a._CompiledPolynomial.compile(deep, cache=cache)
    names = ce._CompiledAnalyticCodec(tmp_path, 'deep').write((first, second))
    record = json.loads((tmp_path/names[0]).read_text())
    restored = ce._CompiledAnalyticCodec.read(record, lambda node:
        np.load(tmp_path/node['file'], mmap_mode='r', allow_pickle=False))
    node = restored[1]
    for _ in range(1500): node = node.instructions[0][1][0]
    assert node is restored[0]
    states = np.array([[0], [2], [5]])
    expected = common.evaluate_polynomial(states)
    for _ in range(1500): expected = 1.001*expected
    with mock.patch.object(a.PolyPath, 'evaluate_polynomial', side_effect=AssertionError('symbolic execution')):
        np.testing.assert_array_equal(restored[1].evaluate_polynomial(states), expected)


def test_previous_flat_numeric_program_record_remains_readable():
    def previous(value):
        if type(value) is tuple: return {'tuple': [previous(item) for item in value]}
        if type(value) is int: return {'integer': str(value)}
        if type(value) is float: return {'float': value.hex()}
        return value
    fields = (((('constant', 3.25, None, ()),)), 1, (0,), 'legacy-constant', False, None)
    record = dict(format='psience-compiled-expression', version=1, root={'ref':0},
                  nodes=[dict(kind='polynomial-program', fields=previous(fields))])
    restored = ce._CompiledAnalyticCodec.read(record, lambda node: None)
    np.testing.assert_array_equal(restored.evaluate_polynomial(np.array([[0], [3]])), [3.25, 3.25])


@pytest.mark.parametrize('reduce', [False, True])
def test_single_compiled_file_evaluates_multiple_systems_without_planning(tmp_path, reduce):
    wrapper = expression()
    store = ce.AnalyticExpressionStore()
    reference = store.compile(wrapper, free_index_reduction=reduce)
    ir_path, compiled_path = tmp_path/'expression.ir.tar', tmp_path/'expression.plan.tar'
    with store.export_ir(wrapper, ir_path, max_terms=3): pass
    compiled = store.compile_ir(ir_path, compiled_path, free_index_reduction=reduce)
    assert compiled.artifact == compiled_path
    assert compiled.manifest['kind'] == 'compiled-analytic-expression'
    recipient = tmp_path/'recipient'
    recipient.mkdir()
    shutil.copyfile(compiled_path, recipient/compiled_path.name)
    assert sorted(p.name for p in recipient.iterdir()) == [compiled_path.name]
    for nmodes in [2, 3, 4]:
        arguments = inputs(nmodes)
        expected = reference.evaluate(*arguments)
        with mock.patch.object(ce, '_implementation_signature', side_effect=AssertionError('source gate')), \
             mock.patch.object(ce, '_prepare_free_index_plan', side_effect=AssertionError('free-index planning')), \
             mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan', side_effect=AssertionError('structural planning')), \
             mock.patch.object(ce.IndexedCoefficientGatherPlan, '_get_request_pattern_map',
                               side_effect=AssertionError('pattern planning')), \
             mock.patch.object(a.PerturbationTheoryExpressionEvaluator, '_get_max_order',
                               side_effect=AssertionError('degree planning')), \
             mock.patch.object(a._CompiledPolynomial, 'compile', side_effect=AssertionError('polynomial compilation')), \
             mock.patch.object(a.PolyAxis, 'coeffs', new_callable=mock.PropertyMock,
                               side_effect=AssertionError('factor convolution')):
            restored = ce.CompiledAnalyticExpression.load(recipient/compiled_path.name)
            actual = restored.evaluate(*arguments)
        for expected_part, actual_part in zip(expected, actual):
            np.testing.assert_allclose(actual_part, expected_part, atol=2.e-14, rtol=2.e-12)
    assert sorted(p.name for p in recipient.iterdir()) == [compiled_path.name]


def test_compiler_accepts_only_ir_and_rejects_output_replacing_input(tmp_path):
    source = tmp_path/'expression.ir.tar'
    store = ce.AnalyticExpressionStore()
    with store.export_ir(expression(), source): pass
    with pytest.raises(ValueError, match='separate'):
        store.compile_ir(source, source)
    with store.export_ir(expression(), staging_folder=tmp_path/'ir') as ir:
        compiled = store.compile_ir(ir, staging_folder=tmp_path/'lowered')
        assert compiled.folder is not None
        assert not (compiled.folder/'evaluations').exists()
        with pytest.raises(ValueError, match='unsupported'):
            ce.AnalyticExpressionIR.load(compiled.folder)


def test_existing_compile_interface_accepts_ir_and_native_file_output(tmp_path):
    store = ce.AnalyticExpressionStore()
    with store.export_ir(expression(), tmp_path/'source.ir.tar'): pass
    from_ir = store.compile(tmp_path/'source.ir.tar', file=tmp_path/'from-ir.plan.tar')
    from_native = store.compile(expression(), file=tmp_path/'from-native.plan.tar')
    for expected, actual in zip(from_ir.evaluate(*inputs(3)), from_native.evaluate(*inputs(3))):
        np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=2e-12)
    assert store.statistics()['compilations'] == 2


def test_all_ir_and_compiled_backend_combinations_share_expression_identity(tmp_path):
    store = ce.AnalyticExpressionStore()
    wrapper = expression()
    expected = store.compile(wrapper).evaluate(*inputs(4))
    signatures, contents = set(), set()
    for ir_backend in ('portable', 'pickle'):
        source = tmp_path/(ir_backend+'.ir.tar')
        with store.export_ir(wrapper, source, storage_backend=ir_backend) as ir:
            assert ir.storage_backend == ir_backend
            contents.add(tuple(part['content_sha256'] for part in ir.manifest['parts']))
        for compiled_backend in ('portable', 'pickle'):
            target = tmp_path/(ir_backend+'-'+compiled_backend+'.plan.tar')
            compiled = store.compile_ir(source, target, storage_backend=compiled_backend)
            signatures.add(compiled.manifest['signature'])
            assert compiled.manifest['storage_backend'] == compiled_backend
            for left, right in zip(expected, compiled.evaluate(*inputs(4))):
                np.testing.assert_allclose(left, right, atol=2e-14, rtol=2e-12)
    assert len(signatures) == len(contents) == 1
    assert not list(tmp_path.rglob('evaluations'))


def test_backend_conversion_and_mixed_collection_need_no_compilation(tmp_path):
    store = ce.AnalyticExpressionStore()
    with store.export_ir(expression(), tmp_path/'source.ir.tar') as source:
        with source.save(tmp_path/'python.ir.tar', storage_backend='pickle') as python_ir:
            assert python_ir.storage_backend == 'pickle'
    compiled = store.compile_ir(tmp_path/'source.ir.tar', tmp_path/'source.plan.tar')
    with mock.patch.object(a._CompiledPolynomial, 'compile', side_effect=AssertionError('recompile')), \
         mock.patch.object(ce, '_prepare_free_index_plan', side_effect=AssertionError('replan')):
        python = compiled.save(tmp_path/'python.plan.tar', storage_backend='pickle')
        restored = python.save(tmp_path/'roundtrip.plan.tar', storage_backend='portable')
        collection = store.pack_ir({'portable':tmp_path/'source.ir.tar', 'python':tmp_path/'python.ir.tar'},
                                   tmp_path/'mixed.ir.tar')
    assert python.manifest['signature'] == restored.manifest['signature'] == compiled.manifest['signature']
    converted = store.compile_ir(tmp_path/'python.ir.tar', tmp_path/'from-saved-ir.plan.tar')
    assert isinstance(converted, ce.CompiledAnalyticExpression)
    assert converted.manifest['signature'] == compiled.manifest['signature']
    for actual in (python, restored):
        for left, right in zip(compiled.evaluate(*inputs(3)), actual.evaluate(*inputs(3))):
            np.testing.assert_array_equal(left, right)
    for name, backend in [('portable', 'portable'), ('python', 'pickle')]:
        with ce.AnalyticExpressionIR.load(collection, expression=name) as entry:
            assert entry.storage_backend == backend
    converted = store.pack_ir({'first label':tmp_path/'source.ir.tar', 'second label':tmp_path/'python.ir.tar'},
                              tmp_path/'uniform.ir.tar', storage_backend='pickle')
    assert ce.AnalyticExpressionIR.expression_names(converted) == ('first label', 'second label')
    for name in ('first label', 'second label'):
        with ce.AnalyticExpressionIR.load(converted, expression=name) as entry:
            assert entry.manifest['storage_version'] == 2
            assert list(entry.iter_parts())
    outputs = store.compile_ir(collection, tmp_path/'library.plan.tar', storage_backend='pickle')
    for entry in outputs.values():
        assert entry.manifest['storage_backend'] == 'pickle'
        assert entry.evaluate(*inputs(3))


def test_storage_backend_validation_and_manifest_dispatch(tmp_path):
    store = ce.AnalyticExpressionStore()
    with pytest.raises(ValueError, match='storage_backend'):
        store.export_ir(expression(), tmp_path/'invalid.ir.tar', storage_backend='unknown')
    with pytest.raises(ValueError, match='storage_backend'):
        store.compile_ir(tmp_path/'missing', tmp_path/'invalid.plan.tar', storage_backend='unknown')
    with pytest.raises(ValueError, match='output'):
        store.compile(expression(), storage_backend='pickle')
    with pytest.raises(ValueError, match='storage_backend'):
        store.compile_factory(expression, context='fixture', storage_backend='pickle')
    with store.export_ir(expression(), staging_folder=tmp_path/'ir', storage_backend='pickle') as ir:
        manifest = dict(ir.manifest, storage_backend='portable', storage_version=1)
        (ir.source/'manifest.json').write_text(json.dumps(manifest))
        with ce.AnalyticExpressionIR.load(ir.source) as wrong:
            with pytest.raises(ValueError, match='storage backend'): next(wrong.iter_parts())


def test_older_record_pickle_storage_upgrades_without_changing_expression(tmp_path):
    store = ce.AnalyticExpressionStore()
    with store.export_ir(expression(), staging_folder=tmp_path/'old') as source:
        parts = []
        for part in source.manifest['parts']:
            record = source._read_record(part)
            filename = Path(part['file']).stem+'.pkl'
            ce._AnalyticIRCodec.write_record(source.source/filename, record, 'pickle')
            parts.append(dict(part, file=filename, sha256=ce._file_hash(source.source/filename)))
        (source.source/'manifest.json').write_text(json.dumps(dict(source.manifest,
            storage_backend='pickle', storage_version=1, parts=parts)))
        folder = source.source
    with ce.AnalyticExpressionIR.load(folder) as old:
        with old.save(tmp_path/'upgraded.ir.tar', storage_backend='pickle') as upgraded:
            assert upgraded.manifest['storage_version'] == 2
            assert [p['content_sha256'] for p in upgraded.manifest['parts']] == [p['content_sha256'] for p in parts]
            compiled = store.compile_ir(upgraded, tmp_path/'upgraded.plan.tar', storage_backend='pickle')
    reference = store.compile(expression()).evaluate(*inputs(3))
    for left, right in zip(reference, compiled.evaluate(*inputs(3))):
        np.testing.assert_allclose(left, right, atol=2e-14, rtol=2e-12)


@pytest.mark.parametrize('backend', ['portable', 'pickle'])
def test_repacking_checks_part_integrity_without_deserializing(tmp_path, backend):
    store = ce.AnalyticExpressionStore()
    with store.export_ir(expression(), staging_folder=tmp_path/'source', storage_backend=backend) as reader:
        part = reader.source/reader.manifest['parts'][0]['file']
        part.write_bytes(part.read_bytes()+b'damaged')
        with mock.patch.object(reader, '_read_record', side_effect=AssertionError('unnecessary decoding')):
            with pytest.raises(ValueError, match='checksum'): reader.save(tmp_path/'invalid.tar')
    assert not (tmp_path/'invalid.tar').exists()


@pytest.mark.parametrize('backend', ['portable', 'pickle'])
def test_empty_expression_and_readonly_numeric_data_in_each_backend(tmp_path, backend):
    store = ce.AnalyticExpressionStore()
    zero = a.PerturbationTheoryExpressionEvaluator(None, a.PTTensorCoeffProductSum({}), change=(), logger=False)
    empty = store.compile(zero, file=tmp_path/'zero.tar', storage_backend=backend)
    for value in empty.evaluate(*inputs(3)): assert not np.any(value)
    compiled = store.compile(expression(), file=tmp_path/'numeric.tar', storage_backend=backend)
    part = next(part for part in compiled.iter_parts() if part.terms)
    value = next(iter(part.terms.values()))
    program = value if isinstance(value, a._CompiledPolynomial) else next(iter(value.terms.values()))
    assert not program.table.parameters.flags.writeable
    with pytest.raises(ValueError): program.table.parameters.flags.writeable = True


def test_existing_pickle_cache_migrates_without_regenerating_symbolic_expression(tmp_path):
    legacy = ce.AnalyticExpressionStore(tmp_path/'legacy').compile(expression())
    expected = legacy.evaluate(*inputs(3))
    with mock.patch.object(a.AnalyticPerturbationTheorySolver, '__init__',
                           side_effect=AssertionError('symbolic regeneration')):
        migrated = legacy.save(tmp_path/'migrated.plan.tar')
    with pytest.raises(ValueError, match='retain.*IR|mathematical IR'):
        ce.AnalyticExpressionStore().export_ir(migrated, tmp_path/'invalid.ir.tar')
    legacy_folder = legacy.folder
    shutil.rmtree(tmp_path/'legacy')
    assert not legacy_folder.exists()
    with mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan', side_effect=AssertionError('recompilation')):
        repacked = migrated.save(tmp_path/'repacked.plan.tar')
    (tmp_path/'migrated.plan.tar').unlink()
    for left, right in zip(expected, repacked.evaluate(*inputs(3))):
        np.testing.assert_allclose(left, right, atol=2e-14, rtol=2e-12)


def test_data_migration_accepts_earlier_producer_and_still_checks_parts(tmp_path):
    store = ce.AnalyticExpressionStore()
    legacy = ce.AnalyticExpressionStore(tmp_path/'legacy').compile(expression())
    expected = legacy.evaluate(*inputs(3))
    manifest = dict(legacy.manifest, implementation='earlier-producer-source')
    (legacy.folder/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='source version changed'):
        ce.CompiledAnalyticExpression.load(legacy.folder)
    with mock.patch.object(a.AnalyticPerturbationTheorySolver, '__init__',
                           side_effect=AssertionError('symbolic regeneration')):
        with store.export_ir(legacy.folder, tmp_path/'migrated.ir.tar'): pass
        migrated = store.compile_ir(tmp_path/'migrated.ir.tar', tmp_path/'migrated.plan.tar')
    for left, right in zip(expected, migrated.evaluate(*inputs(3))):
        np.testing.assert_allclose(left, right, atol=2e-14, rtol=2e-12)
    part = legacy.folder/manifest['parts'][0]['file']
    part.write_bytes(part.read_bytes()+b'damaged')
    with pytest.raises(ValueError, match='checksum'):
        store.export_ir(legacy.folder, tmp_path/'damaged.ir.tar')


@pytest.mark.parametrize('legacy_journal', [False, True])
def test_portable_evaluation_resumes_external_partial_checkpoint_exactly(tmp_path, legacy_journal):
    compiled = ce.AnalyticExpressionStore().compile(resumable_expression(), file=tmp_path/'expression.plan.tar')
    states, coefficients, frequencies = inputs(4)
    expected = compiled.evaluate(states, coefficients, frequencies, batch_size=1, task_size=32)
    evaluation = compiled.prepare_evaluation(states, coefficients, frequencies, batch_size=1, task_size=32,
                                             evaluation_folder=tmp_path/'results')
    native = a.PerturbationTheoryExpressionEvaluator.evaluate_polynomial_expression
    calls = 0
    def interrupted(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 3: raise RuntimeError('interrupted evaluation')
        return native(*args, **kwargs)
    with mock.patch.object(a.PerturbationTheoryExpressionEvaluator, 'evaluate_polynomial_expression', interrupted):
        with pytest.raises(RuntimeError, match='interrupted'): evaluation.run()
    target = evaluation.folder/evaluation.tasks[0]['file']
    record = pickle.loads(target.read_bytes())
    assert not record['complete'] and record['batch_count'] == 2
    journal = target.with_suffix('.batches')
    if legacy_journal:
        batches = []
        for ordinal in range(record.pop('batch_count')):
            path = journal/('batch-{:08d}.pkl'.format(ordinal))
            batches.append(pickle.loads(path.read_bytes())['batch'])
            path.unlink()
        journal.rmdir()
        record['batches'] = batches
        ce._atomic_pickle(target, record)
    for actual in [evaluation.run(), evaluation.run()]:
        for left, right in zip(expected, actual): np.testing.assert_array_equal(left, right)
    assert not journal.exists()
    assert not (tmp_path/'evaluations').exists()


@pytest.mark.parametrize('damage', ['range', 'signature'])
def test_portable_external_partial_checkpoint_rejects_corruption(tmp_path, damage):
    compiled = ce.AnalyticExpressionStore().compile(resumable_expression(), file=tmp_path/'expression.plan.tar')
    states, coefficients, frequencies = inputs(4)
    evaluation = compiled.prepare_evaluation(states, coefficients, frequencies, batch_size=1, task_size=32,
                                             evaluation_folder=tmp_path/'results')
    native = a.PerturbationTheoryExpressionEvaluator.evaluate_polynomial_expression
    calls = 0
    def interrupted(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2: raise RuntimeError('stop')
        return native(*args, **kwargs)
    with mock.patch.object(a.PerturbationTheoryExpressionEvaluator, 'evaluate_polynomial_expression', interrupted):
        with pytest.raises(RuntimeError, match='stop'): evaluation.run()
    target = evaluation.folder/evaluation.tasks[0]['file']
    path = target.with_suffix('.batches')/'batch-00000000.pkl'
    record = pickle.loads(path.read_bytes())
    if damage == 'range': record['batch']['start'] += 1
    else: record['signature'] = 'incorrect'
    ce._atomic_pickle(path, record)
    with pytest.raises(ValueError, match='coverage|context mismatch'): evaluation.run()


def test_compiled_specializations_validate_active_operators_and_mode_capacity(tmp_path):
    store = ce.AnalyticExpressionStore()
    with store.export_ir(expression(), tmp_path/'pruned.ir.tar', active_operators=[(1, 0)]): pass
    compiled = store.compile_ir(tmp_path/'pruned.ir.tar', tmp_path/'pruned.plan.tar')
    with pytest.raises(ValueError, match='excludes active'):
        compiled.evaluate(*inputs(3))
    arguments = inputs(3)
    arguments[1][0][2] = [0, 0, 0, 0, 0]
    actual = compiled.evaluate(*arguments)
    expected = store.compile(expression()).evaluate(*arguments)
    for expected_part, actual_part in zip(expected, actual):
        np.testing.assert_allclose(actual_part, expected_part, rtol=2e-12, atol=2e-14)
    with store.export_ir(expression(), tmp_path/'small.ir.tar', mode_dimension=2): pass
    small = store.compile_ir(tmp_path/'small.ir.tar', tmp_path/'small.plan.tar')
    with pytest.raises(ValueError, match='fewer modes'):
        small.evaluate(*inputs(3))


def test_compiled_action_coefficient_uses_prepared_plan(tmp_path):
    store = ce.AnalyticExpressionStore()
    projection = ce._OccupationCoefficientProjection((1, 1))
    wrapper = expression()
    reference = store.compile(wrapper, projection=projection)
    with store.export_ir(wrapper, tmp_path/'action.ir.tar', projection=projection): pass
    compiled = store.compile_ir(tmp_path/'action.ir.tar', tmp_path/'action.plan.tar')
    _, coefficients, frequencies = inputs(4)
    expected = reference.evaluate_coefficient(coefficients, frequencies, indices=(0, 2))
    with mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan', side_effect=AssertionError('runtime planning')):
        actual = compiled.evaluate_coefficient(coefficients, frequencies, indices=(0, 2))
    np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=2e-12)


def test_expression_collections_compile_and_transfer_as_one_file(tmp_path):
    store = ce.AnalyticExpressionStore()
    solver = a.AnalyticPerturbationTheorySolver.from_order(2, internals=False, logger=False)
    wrappers = {'energy/order 2': expression(), 'transition Δn=1': solver.wavefunction_correction(1)((1,))}
    references, sources = {}, {}
    for index, (name, wrapper) in enumerate(wrappers.items()):
        references[name] = store.compile(wrapper)
        sources[name] = tmp_path/('source-{}.ir.tar'.format(index))
        with store.export_ir(wrapper, sources[name]): pass
    library = store.pack_ir(sources, tmp_path/'library.ir.tar')
    assert ce.AnalyticExpressionIR.expression_names(library) == tuple(wrappers)
    with pytest.raises(ValueError, match='select an expression'):
        ce.AnalyticExpressionIR.load(library)
    with ce.AnalyticExpressionIR.load(library, expression='energy/order 2') as selected:
        assert list(selected.iter_parts())
    compiled = store.compile_ir(library, tmp_path/'library.plan.tar')
    assert set(compiled) == set(wrappers)
    assert {value.artifact for value in compiled.values()} == {tmp_path/'library.plan.tar'}
    recipient = tmp_path/'recipient'
    recipient.mkdir()
    shutil.copyfile(tmp_path/'library.plan.tar', recipient/'library.plan.tar')
    assert ce.CompiledAnalyticExpression.expression_names(recipient/'library.plan.tar') == tuple(wrappers)
    with pytest.raises(ValueError, match='select an expression'):
        ce.CompiledAnalyticExpression.load(recipient/'library.plan.tar')
    for name in wrappers:
        expected = references[name].evaluate(*inputs(3))
        restored = ce.CompiledAnalyticExpression.load(recipient/'library.plan.tar', expression=name)
        with mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan', side_effect=AssertionError('runtime planning')):
            actual = restored.evaluate(*inputs(3))
        for expected_part, actual_part in zip(expected, actual):
            np.testing.assert_allclose(actual_part, expected_part, atol=2e-14, rtol=2e-12)
    assert sorted(path.name for path in recipient.iterdir()) == ['library.plan.tar']


@pytest.mark.parametrize('condition_kind', ['invariant', 'relationship', 'lattice', 'window'])
@pytest.mark.parametrize('backend', ['portable', 'pickle'])
def test_conditional_compiled_expression_uses_prepared_functions_and_denominators(tmp_path, condition_kind, backend):
    q = Symbols('q').q
    condition = {'invariant':lambda: a.PolynomialInvariantCondition(q*(q-1), [0]),
                 'relationship':lambda: a.PolynomialRelationshipCondition(q, [0], comparison='ge'),
                 'lattice':lambda: a.IntegerLatticeCondition([[1, -1]]),
                 'window':lambda: a.EnergyWindowCondition(.001)}[condition_kind]()
    key = ((a.ConditionEndpoint(condition, 'resolvent', False, 0), 0, 0),
           (a.ConditionEndpoint(condition, 'resolvent', False, 1), 1, -1))
    polynomial = a.PolyPath.from_coeffs([np.array([1., 2.]), np.array([1.])])
    energy = a.PTEnergyChangeProductSum({key: polynomial}, canonicalize=False)
    part = a.PTTensorCoeffProductSum({((0, 0, 0, 1),): energy}, canonicalize=False)
    wrapper = a.PerturbationTheoryExpressionEvaluator(None, a.SqrtChangePoly(part, [], []), [], logger=False)
    store = ce.AnalyticExpressionStore()
    reference = store.compile(wrapper)
    with store.export_ir(wrapper, tmp_path/'condition.ir.tar', storage_backend=backend): pass
    compiled = store.compile_ir(tmp_path/'condition.ir.tar', tmp_path/'condition.plan.tar', storage_backend=backend)
    arguments = inputs(3)
    arguments[1][0][0][0] = np.ones((3, 3))*.1
    expected = reference.evaluate(*arguments)
    from McUtils.Zachary.Symbolic.ElementaryFunctions import Functionlike
    with mock.patch.object(a, '_compile_energy_program', side_effect=AssertionError('denominator compilation')), \
         mock.patch.object(Functionlike, 'compile', side_effect=AssertionError('function compilation')), \
         mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan', side_effect=AssertionError('expression compilation')):
        restored = ce.CompiledAnalyticExpression.load(compiled.artifact)
        actual = restored.evaluate(*arguments)
    for expected_part, actual_part in zip(expected, actual):
        np.testing.assert_allclose(actual_part, expected_part, rtol=2e-12, atol=2e-14)


def test_compiled_projection_program_retains_batched_projector_semantics():
    key = ((a.ProjectionFactor('Q'), 1, -1), (a.ProjectionFactor('R'), 1, -1), (0, 1, 0))
    states = [np.array([[1, 1]])]
    permutations = [np.array([[0, 1], [1, 0]])]
    frequencies = np.array([[[.021, .032], [.032, .021]]])
    bindings = ((tuple([(1, 1), (2, 0)]), (.021, .032), (1, 1), False),)
    expected = a.projected_energy_rows([key], states, permutations, frequencies, projection_bindings=bindings)
    program = a._compile_energy_program(key)
    with mock.patch.object(a, '_bound_projection_row', side_effect=AssertionError('structural row binding')):
        actual = a.projected_energy_rows([key], states, permutations, frequencies,
                                         projection_bindings=bindings, programs=[program])
    for expected_part, actual_part in zip(expected, actual):
        np.testing.assert_array_equal(actual_part, expected_part)


@pytest.mark.parametrize('method', ['linear', 'perfect'])
def test_runtime_resonance_binding_uses_compiled_energy_requirements(tmp_path, method):
    wrapper = expression()
    store = ce.AnalyticExpressionStore()
    reference = store.compile(wrapper)
    with store.export_ir(wrapper, tmp_path/'resonant.ir.tar'): pass
    compiled = store.compile_ir(tmp_path/'resonant.ir.tar', tmp_path/'resonant.plan.tar')
    arguments = inputs(3)
    changes = [((0,), (1,)), ((1,), (-1,)), ((0, 1), (1, -1))]
    with mock.patch.object(a.PerturbationTheoryExpressionEvaluator, 'default_deg_id_method', method):
        expected = reference.evaluate(*arguments, degenerate_changes=changes)
        with mock.patch.object(a.DegeneracyIdentificationContext, '_requirement',
                               side_effect=AssertionError('energy requirement planning')), \
             mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan',
                               side_effect=AssertionError('structural fallback planning')):
            actual = compiled.evaluate(*arguments, degenerate_changes=changes)
    for expected_part, actual_part in zip(expected, actual):
        np.testing.assert_allclose(actual_part, expected_part, atol=2e-14, rtol=2e-12)


@pytest.mark.parametrize('change', [(1,), (-1,)])
def test_transition_program_uses_compiled_occupation_offsets(tmp_path, change):
    solver = a.AnalyticPerturbationTheorySolver.from_order(2, internals=False, logger=False)
    wrapper = solver.wavefunction_correction(1)(change)
    store = ce.AnalyticExpressionStore()
    reference = store.compile(wrapper)
    with store.export_ir(wrapper, tmp_path/'transition.ir.tar'): pass
    compiled = store.compile_ir(tmp_path/'transition.ir.tar', tmp_path/'transition.plan.tar')
    arguments = inputs(3)
    expected = reference.evaluate(*arguments)
    actual = compiled.evaluate(*arguments)
    for expected_part, actual_part in zip(expected, actual):
        np.testing.assert_allclose(actual_part, expected_part, atol=2e-14, rtol=2e-12)


@pytest.mark.parametrize('ir_backend,compiled_backend', [
    ('portable', 'portable'), ('pickle', 'portable'), ('portable', 'pickle'), ('pickle', 'pickle')])
def test_fresh_process_compiles_transferred_ir_then_evaluates_without_source_ir(tmp_path, ir_backend, compiled_backend):
    source, target = tmp_path/'source.ir.tar', tmp_path/'transferred.plan.tar'
    with ce.AnalyticExpressionStore().export_ir(expression(), source, storage_backend=ir_backend): pass
    reference = ce.AnalyticExpressionStore().compile(expression())
    expectations = [[array.tolist() for array in reference.evaluate(*inputs(n))] for n in [2, 3, 4]]
    expected_file = tmp_path/'reference.json'
    expected_file.write_text(json.dumps(expectations))
    compile_code = '''
import sys
from unittest import mock
from Psience.VPT2 import Analytic as a, Evaluators as ce
with mock.patch.object(a.AnalyticPerturbationTheorySolver, '__init__', side_effect=AssertionError('solver construction')):
    ce.AnalyticExpressionStore().compile_ir(sys.argv[1], sys.argv[2], storage_backend=sys.argv[3])
'''
    completed = subprocess.run([sys.executable, '-c', compile_code, str(source), str(target), compiled_backend],
                               cwd=tmp_path, text=True, capture_output=True, timeout=60)
    assert completed.returncode == 0, completed.stderr
    source.unlink()
    recipient = tmp_path/'recipient'
    recipient.mkdir()
    shutil.copyfile(target, recipient/target.name)
    evaluation_code = '''
import json, sys
from unittest import mock
import numpy as np
from Psience.VPT2 import Analytic as a, Evaluators as ce
from McUtils.Zachary.Symbolic.ElementaryFunctions import Functionlike
expected = json.loads(open(sys.argv[2]).read())
def fail(*args, **kwargs): raise AssertionError('structural compiler called during evaluation')
with mock.patch.object(ce, '_implementation_signature', fail), \
     mock.patch.object(ce, '_prepare_free_index_plan', fail), \
     mock.patch.object(ce.IndexedBlockEvaluator, 'compile_plan', fail), \
     mock.patch.object(ce.IndexedCoefficientGatherPlan, '_get_request_pattern_map', fail), \
     mock.patch.object(a._CompiledPolynomial, 'compile', fail), \
     mock.patch.object(a._PolynomialProgramTable, 'pack', fail), \
     mock.patch.object(a.PolyPath, 'evaluate_polynomial', fail), \
     mock.patch.object(a, '_compile_energy_program', fail), \
     mock.patch.object(a.PerturbationTheoryExpressionEvaluator, '_get_max_order', fail), \
     mock.patch.object(a.PolyAxis, 'coeffs', new_callable=mock.PropertyMock, side_effect=fail), \
     mock.patch.object(Functionlike, 'compile', fail), \
     mock.patch.object(a.AnalyticPerturbationTheorySolver, '__init__', fail):
    expression = ce.CompiledAnalyticExpression.load(sys.argv[1])
    for index, n in enumerate([2, 3, 4]):
        rng = np.random.default_rng(97)
        frequencies = .017*np.sqrt(np.arange(2, n+2))
        coefficients = [[np.diag(frequencies)/2, np.diag(frequencies)/2],
                        [rng.normal(size=(n,)*3)*1.e-4, 0],
                        [rng.normal(size=(n,)*4)*1.e-5, 0, 0, 0, 0]]
        states = [(np.zeros(n, dtype=int), np.arange(n)[None]),
                  (np.arange(n)%3, np.array([np.arange(n), np.arange(n)[::-1]]))]
        actual = expression.evaluate(states, [coefficients], frequencies)
        for result, reference in zip(actual, expected[index]):
            np.testing.assert_allclose(result, reference, atol=2e-14, rtol=2e-12)
print(json.dumps({'systems':3, 'planners_disabled':True}))
'''
    completed = subprocess.run([sys.executable, '-c', evaluation_code, str(recipient/target.name), str(expected_file)],
                               cwd=recipient, text=True, capture_output=True, timeout=60)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == dict(systems=3, planners_disabled=True)
    assert sorted(p.name for p in recipient.iterdir()) == [target.name]
