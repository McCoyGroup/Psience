"""Portable semantic expression records and independent numerical storage."""
import json
import pickle
from pathlib import Path
import shutil
from unittest import mock

import numpy as np
import pytest

from Psience.VPT2 import Analytic as a
from Psience.VPT2 import Evaluators as ce


def codec_roundtrip(tmp_path, value):
    names = ce._AnalyticIRCodec(tmp_path, 'part').write(value)
    record = json.loads((tmp_path / names[0]).read_text())
    def load_array(node):
        return np.load(tmp_path / node['file'], mmap_mode='r', allow_pickle=False)
    return ce._AnalyticIRCodec.read(record, load_array), record


def test_native_polynomial_graph_sharing_tags_and_numerics(tmp_path):
    leaf = a.PolyPath.from_coeffs([np.array([1., 2.], dtype='>f8'), np.array([3., -.5])])
    graph = (leaf.shift((1, 0)) * leaf).permute((1, 0))
    condition = a.IntegerLatticeCondition([[1, -1]]) & a.EnergyWindowCondition(.03)
    endpoint = a.ConditionEndpoint(condition, 'resolvent', False, 0)
    roles = (a.ProjectionFactor('Q', 2), a._ChangeDomainIndex(True, 1), endpoint)
    restored, record = codec_roundtrip(tmp_path, [graph, graph, roles, complex(-0., .5), np.float32(-0.)])
    assert restored[0] is restored[1]
    assert type(restored[2][0]) is a.ProjectionFactor
    assert restored[2][0].kind == 'Q' and restored[2][0].slot == 2
    assert type(restored[2][1]) is a._ChangeDomainIndex
    assert restored[2][1].retain and restored[2][1].endpoint == 1
    assert type(restored[2][2]) is a.ConditionEndpoint
    assert restored[2][2].condition == condition
    initial = np.array([[1, 0], [0, 2], [4, 3]])
    final = np.array([[0, 1], [2, 0], [4, 4]])
    np.testing.assert_array_equal(restored[2][2].condition.matches(initial, final, [.01, .02]),
                                  condition.matches(initial, final, [.01, .02]))
    assert type(restored[4]) is np.float32 and np.signbit(restored[4])
    states = np.array([[0, 0], [1, 3], [2, 4], [7, 1]])
    np.testing.assert_array_equal(restored[0].evaluate_polynomial(states), graph.evaluate_polynomial(states))
    assert 'pickle' not in json.dumps(record)
    assert len(list(tmp_path.glob('*.npy'))) == 1


def test_deep_shared_graph_roundtrip_is_iterative(tmp_path):
    leaf = a.PolyPath.from_coeffs([np.array([1., 1.])])
    graph = leaf
    for _ in range(1500):
        graph = a.PolyPath((), node=('add', (graph, leaf)), ndim=1, order=(1,), reduced=True)
    restored, _ = codec_roundtrip(tmp_path, graph)
    count = 0
    while restored._node is not None:
        restored = restored._node[1][0]
        count += 1
    assert count == 1500


def test_tensor_sum_roundtrip_excludes_execution_and_runtime_caches(tmp_path):
    path = a.PolyPath.from_coeffs([np.array([1., 1.]), np.array([1.])])
    value = a.PTEnergyChangeProductSum({((0, 1, -1),): path}, canonicalize=False)
    part = a.PTTensorCoeffProductSum({((1, 0, 0, 1),): value}, canonicalize=False, reduced=True)
    part._free_index_plan_data = {'must': 'not enter simplified IR'}
    part._cache = object()
    restored, record = codec_roundtrip(tmp_path, part)
    assert not hasattr(restored, '_free_index_plan_data')
    assert not hasattr(restored, '_cache')
    assert tuple(restored.terms) == tuple(part.terms)
    assert tuple(next(iter(restored.terms.values())).terms) == tuple(value.terms)


def test_unsupported_custom_objects_and_cyclic_or_unknown_records_fail(tmp_path):
    with pytest.raises(TypeError, match='unsupported analytic IR value'):
        ce._AnalyticIRCodec(tmp_path, 'bad').write(lambda x: x)
    _, record = codec_roundtrip(tmp_path, a.PolyPath.from_coeffs([np.array([1.])]))
    record['nodes'][0] = {'kind': 'polynomial-path', 'fields': {'ref': 0}}
    with pytest.raises(ValueError, match='cyclic'):
        ce._AnalyticIRCodec.read(record, lambda n: None)
    record['nodes'][0] = {'kind': 'unknown', 'fields': {'tuple': []}}
    with pytest.raises(ValueError, match='unknown analytic IR node'):
        ce._AnalyticIRCodec.read(record, lambda n: None)


def test_primitive_buffer_shapes_precision_and_repeated_arrays(tmp_path):
    arrays = [np.empty((0, 3)), np.array(2.5), np.asfortranarray(np.arange(24).reshape(4, 6)),
              np.array([-0., np.inf, -np.inf, np.nan]), np.array([1+2j]), np.array(['abc', 'δ'])]
    restored, record = codec_roundtrip(tmp_path, arrays + [arrays[2], 2**100])
    for expected, actual in zip(arrays, restored):
        np.testing.assert_array_equal(actual, expected)
    assert restored[2] is restored[-2]
    assert restored[-1] == 2**100
    assert np.signbit(restored[3][0])
    with pytest.raises(TypeError, match='non-object'):
        ce._AnalyticIRCodec(tmp_path, 'object').write(np.array([object()], dtype=object))


def test_compact_values_preserve_large_integers_and_read_previous_records(tmp_path):
    expected = (7, -3, 2**100, [1, 2], ('x', 4))
    actual, record = codec_roundtrip(tmp_path, expected)
    assert actual == expected
    assert record['version'] == 2
    assert record['root'] == [7, -3, {'integer': str(2**100)}, {'list': [1, 2]}, ['x', 4]]
    previous = dict(format='psience-analytic-expression', version=1,
                    root={'tuple': [{'integer': '7'}, {'integer': str(2**100)}]}, nodes=[])
    assert ce._AnalyticIRCodec.read(previous, lambda node: None) == (7, 2**100)


def test_equal_immutable_buffers_share_storage_without_merging_signed_zero(tmp_path):
    values = [np.array([1., 2.], dtype='>f8'), np.array([1., 2.], dtype='<f8'),
              np.array([-0.]), np.array([0.])]
    restored, record = codec_roundtrip(tmp_path, values)
    assert restored[0] is restored[1]
    assert len([node for node in record['nodes'] if node['kind'] == 'array']) == 3
    assert not restored[0].flags.writeable
    assert np.signbit(restored[2][0]) and not np.signbit(restored[3][0])


def test_pickle_backend_preserves_deep_graphs_without_changing_recursion_limit(tmp_path):
    root = a.PolyPath.from_coeffs([np.array([1., 2.])])
    for _ in range(1500): root = a.PolyPath((), node=('scale', (root, 1.001)), ndim=1, order=(1,))
    codec = ce._AnalyticIRCodec(tmp_path, 'deep')
    files = codec.write(root, storage_backend='pickle')
    record = pickle.loads((tmp_path/files[0]).read_bytes())
    assert record['representation'] == 'graph'
    metadata = dict(file=files[0], sha256=ce._file_hash(tmp_path/files[0]), buffers=files[1:])
    (tmp_path/'manifest.json').write_text(json.dumps(dict(kind='analytic-expression-ir', ir_version=1,
        complete=True, storage_backend='pickle', storage_version=2, parts=[metadata])))
    with ce.AnalyticExpressionIR.load(tmp_path) as reader: restored = next(reader.iter_parts())
    program = a._CompiledPolynomial.compile(restored)
    expected = np.array([5.])
    for _ in range(1500): expected = 1.001*expected
    np.testing.assert_array_equal(program.evaluate_polynomial(np.array([[2]])), expected)


def native_expression():
    solver = a.AnalyticPerturbationTheorySolver.from_order(2, internals=False, logger=False)
    return solver.energy_correction(2)([])


def test_export_simplifies_in_bounded_parts_without_execution_planning(tmp_path):
    output = tmp_path / 'expression.ir.tar'
    with mock.patch.object(ce, '_prepare_free_index_plan', side_effect=AssertionError('execution planning')):
        with ce.AnalyticExpressionStore().export_ir(native_expression(), output, max_terms=3) as ir:
            assert ir.manifest['kind'] == 'analytic-expression-ir'
            assert ir.manifest['specializations'] == {'active_operators': None, 'mode_dimension': None}
            parts = list(ir.iter_parts())
            assert parts and all(type(part) is a.PTTensorCoeffProductSum for part in parts)
            assert not any(hasattr(part, '_free_index_plan_data') for part in parts)
    assert sorted(p.name for p in tmp_path.iterdir()) == ['expression.ir.tar']
    recipient = tmp_path / 'recipient'
    recipient.mkdir()
    shutil.copyfile(output, recipient / output.name)
    with mock.patch.object(ce, '_implementation_signature', side_effect=AssertionError('source fingerprint gate')):
        with ce.AnalyticExpressionIR.load(recipient / output.name) as copied:
            assert len(list(copied.iter_parts())) == len(parts)
    assert sorted(p.name for p in recipient.iterdir()) == [output.name]


def test_directory_and_archive_ir_preserve_native_evaluation(tmp_path):
    expression = native_expression()
    store = ce.AnalyticExpressionStore()
    expected = store.compile(expression, max_terms=3, free_index_reduction=False)
    states = [(np.array([0, 1]), np.array([[0, 1]])), (np.array([2, 0]), np.array([[0, 1]]))]
    freqs = np.array([.017, .026])
    coefficients = [[np.diag(freqs) / 2, np.diag(freqs) / 2],
                    [np.ones((2, 2, 2)) * .0001, 0], [np.ones((2,) * 4) * .00001, 0, 0, 0, 0]]
    reference = expected.evaluate(states, [coefficients], freqs)
    with store.export_ir(expression, staging_folder=tmp_path / 'build', max_terms=3) as directory:
        with pytest.raises(ValueError, match='must be compiled'):
            ce.CompiledAnalyticExpression.load(directory.source)
        parts = {meta['file']: part for meta, part in zip(directory.manifest['parts'], directory.iter_parts())}
        # Execute the reconstructed mathematical expression through the unchanged
        # native evaluator. This is not yet the standalone IR-to-plan compiler.
        restored = ce.CompiledAnalyticExpression(None, directory.manifest, parts=parts)
        actual = restored.evaluate(states, [coefficients], freqs)
    for left, right in zip(reference, actual):
        np.testing.assert_allclose(left, right, atol=2e-14, rtol=2e-12)


def test_ir_rejects_corrupt_parts_buffers_and_unsupported_version(tmp_path):
    with ce.AnalyticExpressionStore().export_ir(native_expression(), staging_folder=tmp_path) as ir:
        directory = ir.source
        first = directory / ir.manifest['parts'][0]['file']
        original = first.read_bytes()
        first.write_bytes(original + b' ')
        with pytest.raises(ValueError, match='part checksum'):
            list(ir.iter_parts())
        first.write_bytes(original)
        node = next(node for part in ir.manifest['parts']
                    for node in json.loads((directory / part['file']).read_text())['nodes']
                    if node['kind'] == 'array')
        buffer = directory / node['file']
        data = buffer.read_bytes()
        buffer.write_bytes(data[:-1] + bytes([data[-1] ^ 1]))
        with pytest.raises(ValueError, match='array checksum'):
            list(ir.iter_parts())
        manifest = dict(ir.manifest, ir_version=1000)
        (directory / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='unsupported'):
        ce.AnalyticExpressionIR.load(directory)


def test_evaluation_writes_only_to_explicit_independent_destination(tmp_path):
    compiled = ce.AnalyticExpressionStore(tmp_path / 'compiled').compile(native_expression())
    before = {str(p.relative_to(tmp_path / 'compiled')): p.read_bytes()
              for p in (tmp_path / 'compiled').rglob('*') if p.is_file()}
    # Empty fixtures exercise preparation/publication without introducing physical
    # tensor assumptions. The full native energy check above covers contraction.
    empty = dict(compiled.manifest, parts=[])
    expression = ce.CompiledAnalyticExpression(compiled.folder, empty)
    states = [(np.zeros(2, dtype=int), np.arange(2)[None])]
    for checkpoint in [None, tmp_path / 'evaluations']:
        plan = expression.prepare_evaluation(states, [[[0.]]], [.1, .2], evaluation_folder=checkpoint)
        assert (plan.folder is None) == (checkpoint is None)
        plan.run()
    after = {str(p.relative_to(tmp_path / 'compiled')): p.read_bytes()
             for p in (tmp_path / 'compiled').rglob('*') if p.is_file()}
    assert after == before
    assert not list((tmp_path / 'compiled').rglob('evaluations'))
    assert list((tmp_path / 'evaluations').rglob('result.json'))
    with pytest.raises(ValueError, match='outside'):
        expression.prepare_evaluation(states, [[[0.]]], [.1, .2], evaluation_folder=compiled.folder / 'evaluations')


def test_memory_expression_can_checkpoint_and_context_controls_destination(tmp_path):
    expression = ce.CompiledAnalyticExpression(None, dict(signature='fixture', options={}, projection=None,
                                                       num_fixed=0, parts=[]))
    with ce.compiled_evaluation_context(ce.AnalyticExpressionStore(),
                                       evaluation_options={'evaluation_folder': tmp_path / 'runs'}):
        plan = expression.prepare_evaluation([(np.zeros(1, dtype=int), np.zeros((1, 1), dtype=int))],
                                             [[[0.]]], [.1])
        assert plan.folder.parent == tmp_path / 'runs'
        plan.run()
    assert list((tmp_path / 'runs').rglob('result.json'))


def test_explicit_evaluation_destination_resumes_native_results(tmp_path):
    expression = ce.AnalyticExpressionStore(tmp_path / 'compiled').compile(native_expression(), max_terms=3)
    frequencies = np.array([.017, .026])
    coefficients = [[np.diag(frequencies) / 2, np.diag(frequencies) / 2],
                    [np.ones((2,) * 3) * .0001, 0], [np.ones((2,) * 4) * .00001, 0, 0, 0, 0]]
    states = [(np.array([0, 1]), np.array([[0, 1]]))]
    first = expression.prepare_evaluation(states, [coefficients], frequencies,
                                          evaluation_folder=tmp_path / 'runs')
    reference = first.run()
    with mock.patch.object(a.PerturbationTheoryExpressionEvaluator, 'evaluate_polynomial_expression',
                           side_effect=AssertionError('completed task was reevaluated')):
        actual = expression.evaluate(states, [coefficients], frequencies, evaluation_folder=tmp_path / 'runs')
    for expected, value in zip(reference, actual):
        np.testing.assert_array_equal(expected, value)
    result = json.loads((first.folder / 'result.json').read_text())
    assert result['resumed_tasks'] == len(first.tasks) > 0
