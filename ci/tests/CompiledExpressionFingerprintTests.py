"""
Regression: `Psience.VPT2.CompiledExpressions` fingerprints now come from `McUtils.Devutils`.
The key of every compiled-expression store entry depends on these digests, so they must stay
byte-identical to the original in-module implementation, which is embedded below.
"""

import functools, hashlib, inspect, json, marshal, os, pickle, sys, tempfile, unittest, uuid
from pathlib import Path
import numpy as np

from McUtils.Parallelizers import Parallelizer, SharedMemoryArrayTree
from Psience.VPT2 import CompiledExpressions as ce

__all__ = ["CompiledExpressionFingerprintTests"]

_ORIGINAL = {}
_SOURCE = r'''
def _atomic_json(path, value):
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        temporary.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_pickle(path, value):
    # Preserve shared native PolyPath nodes. JSON/NPZ object conversion expands
    # those graphs; ordinary pickle is the storage format, not a new algebra.
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temporary.open('wb') as stream:
            pickle.dump(value, stream, protocol=5)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _file_hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _freeze(value, seen=None):
    """Fingerprint definitions without including mutable evaluation caches."""
    if seen is None:
        seen = {}
    from McUtils.Parallelizers import SharedMemoryArrayTree
    if isinstance(value, SharedMemoryArrayTree):
        return _freeze(value.tree, seen)
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, np.generic):
        return _freeze(value.item(), seen)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        array = np.ascontiguousarray(value)
        if array.dtype.hasobject:
            return ['object-array', array.shape, _freeze(array.tolist(), seen)]
        return ['array', array.shape, array.dtype.str,
                hashlib.sha256(memoryview(array).cast('B')).hexdigest()]
    if isinstance(value, (list, tuple)):
        return [type(value).__name__, [_freeze(v, seen) for v in value]]
    if isinstance(value, (set, frozenset)):
        return ['set', sorted((_freeze(v, seen) for v in value), key=repr)]
    if isinstance(value, dict):
        return ['dict', sorted([[_freeze(k, seen), _freeze(v, seen)]
                               for k, v in value.items()], key=lambda p:repr(p[0]))]
    if isinstance(value, functools.partial):
        return ['partial', _freeze(value.func,seen), _freeze(value.args,seen), _freeze(value.keywords,seen)]
    if inspect.ismethod(value):
        return ['method', _freeze(value.__self__, seen), _freeze(value.__func__, seen)]
    if inspect.isfunction(value):
        return ['function', value.__module__, value.__qualname__,
                hashlib.sha256(marshal.dumps(value.__code__)).hexdigest(),
                _freeze(value.__defaults__, seen),
                _freeze(tuple(c.cell_contents for c in (value.__closure__ or ())), seen)]
    identity = id(value)
    if identity in seen:
        return ['reference', seen[identity]]
    seen[identity] = len(seen)
    kind = type(value).__module__ + '.' + type(value).__qualname__
    from .Analytic import AnalyticPerturbationTheorySolver
    if isinstance(value, AnalyticPerturbationTheorySolver):
        names = ['hamiltonian_expansion', 'allowed_terms', 'allowed_coefficients',
                 'disallowed_coefficients', 'allowed_energy_changes',
                 'intermediate_normalization', 'polynomial_representation', 'resonance_condition']
        return [kind, {name:_freeze(getattr(value, name), seen) for name in names}]
    if isinstance(value, Parallelizer):
        raise TypeError('a parallelizer is runtime state, not an expression definition')
    attrs = getattr(value, '__dict__', None)
    if attrs is not None:
        runtime = {'logger', 'checkpoint', '_changes', '_exprs', '_combined', '_eager',
                   '_op_maps', '_conditional_stacks', '_cache', '_grouped_terms'}
        return [kind, {name:_freeze(v, seen) for name, v in sorted(attrs.items()) if name not in runtime}]
    slots = {}
    for base in type(value).__mro__:
        names = getattr(base, '__slots__', ())
        if isinstance(names, str):names = (names,)
        for name in names:
            if name not in {'__weakref__', '__dict__'} and hasattr(value, name):
                slots[name] = _freeze(getattr(value, name), seen)
    if slots:
        return [kind, slots]
    raise TypeError('no compilation signature for {}'.format(kind))


def _fingerprint(value):
    return hashlib.sha256(json.dumps(_freeze(value), sort_keys=True).encode()).hexdigest()


@functools.lru_cache(None)
def _implementation_signature():
    import McUtils, Psience
    result = dict(numpy=np.__version__, python=list(sys.version_info[:2]))
    for module in [McUtils, Psience]:
        root=Path(module.__file__).parent;digest=hashlib.sha256()
        for path in sorted(root.rglob('*.py')):
            digest.update(str(path.relative_to(root)).encode());digest.update(bytes.fromhex(_file_hash(path)))
        result[module.__name__]=digest.hexdigest()
    return result



'''
_ORIGINAL.update(dict(
    __name__="Psience.VPT2.CompiledExpressions", __package__="Psience.VPT2",
    functools=functools, hashlib=hashlib, inspect=inspect, json=json, marshal=marshal, os=os,
    pickle=pickle, sys=sys, uuid=uuid, Path=Path, np=np, Parallelizer=Parallelizer
))
exec(compile(_SOURCE, "<original CompiledExpressions helpers>", "exec"), _ORIGINAL)


class _Plain:
    def __init__(self, **kw):
        self.__dict__.update(kw)


class _Slotted:
    __slots__ = ('a', 'b')
    def __init__(self, a, b):
        self.a, self.b = a, b


def _fn(x, y=2):
    return x + y


class CompiledExpressionFingerprintTests(unittest.TestCase):

    def samples(self):
        shared = _Plain(v=np.arange(4.))
        yield {'a': [1, 2.5, 'x', None, True], 'b': (np.arange(6).reshape(2, 3), {'k': {1, 2, 3}})}
        yield _Plain(a=np.eye(3), logger=object(), _cache={'x': 1}, name='n')
        yield [shared, shared, {'again': shared}]
        yield _Slotted(np.ones(2), 'z')
        yield functools.partial(_fn, 1, y=3)
        yield np.array([{'a': 1}, [2]], dtype=object)
        yield Path('/tmp/some/where')
        yield np.float64(3.5)

    def test_ImmutableConditionFingerprintIgnoresPreparedFunctionState(self):
        from McUtils.Zachary import Symbols
        from Psience.VPT2.Conditional import PolynomialInvariantCondition
        q = Symbols('q').q
        condition = PolynomialInvariantCondition(q*(q-1), [0])
        import pickle
        same = pickle.loads(pickle.dumps(condition))
        before = ce._fingerprint(condition)
        self.assertEqual(before, ce._fingerprint(same))
        condition.matches([0], [1])
        condition.functions[0].__dict__['prepared_runtime_note'] = object()
        self.assertEqual(before, ce._fingerprint(condition))
        different = PolynomialInvariantCondition(q*q, [0])
        self.assertNotEqual(before, ce._fingerprint(different))

    def test_FreezeMatchesOriginal(self):
        for sample in self.samples():
            self.assertEqual(json.dumps(ce._freeze(sample), sort_keys=True),
                             json.dumps(_ORIGINAL['_freeze'](sample), sort_keys=True))

    def test_FingerprintMatchesOriginal(self):
        for sample in self.samples():
            self.assertEqual(ce._fingerprint(sample), _ORIGINAL['_fingerprint'](sample))

    def test_SharedMemoryArrayTreeMatchesOriginal(self):
        tree = SharedMemoryArrayTree.create({'a': np.arange(5.), 'b': [np.ones(2)]})
        try:
            self.assertEqual(ce._fingerprint(tree), _ORIGINAL['_fingerprint'](tree))
        finally:
            tree.close()
            tree.unlink()

    def test_ParallelizerRejected(self):
        with self.assertRaises(TypeError):
            ce._fingerprint(_Plain(p=Parallelizer.lookup(None)))

    def test_ImplementationSignatureMatchesOriginal(self):
        self.assertEqual(ce._implementation_signature(), _ORIGINAL['_implementation_signature']())

    def test_FileHashAndAtomicWritersMatchOriginal(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            ce._atomic_json(d / 'new.json', {'b': 1, 'a': [1, 2]})
            _ORIGINAL['_atomic_json'](d / 'old.json', {'b': 1, 'a': [1, 2]})
            self.assertEqual((d / 'new.json').read_bytes(), (d / 'old.json').read_bytes())
            ce._atomic_pickle(d / 'new.pkl', {'x': np.arange(3)})
            _ORIGINAL['_atomic_pickle'](d / 'old.pkl', {'x': np.arange(3)})
            self.assertEqual((d / 'new.pkl').read_bytes(), (d / 'old.pkl').read_bytes())
            self.assertEqual(ce._file_hash(d / 'new.pkl'), _ORIGINAL['_file_hash'](d / 'old.pkl'))
            self.assertEqual(sorted(p.name for p in d.iterdir()), ['new.json', 'new.pkl', 'old.json', 'old.pkl'])

    def _numerical_plan_fixture(self, folder=None):
        manifest=dict(signature='coefficient-binding-fixture',options={},
                      projection=None,num_fixed=0,parts=[])
        return ce.CompiledAnalyticExpression(folder,manifest)

    def test_NumericalContextBindsTensorOnceAndPreservesSignature(self):
        from unittest import mock
        array=np.arange(24.,dtype=float).reshape(4,6).T
        coefficients=[[[array,2.]]]
        states=[(np.array([0,0]),np.array([[0,1]]))]
        expression=self._numerical_plan_fixture()
        expected=expression.prepare_evaluation(states,coefficients,[.1,.2]).signature
        native=np.ascontiguousarray;calls=[]
        def observe(value,*args,**kwargs):
            if value is array:calls.append(value.shape)
            return native(value,*args,**kwargs)
        with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()), \
                mock.patch.object(np,'ascontiguousarray',side_effect=observe):
            for _ in range(4):
                self.assertEqual(expression.prepare_evaluation(states,coefficients,[.1,.2]).signature,expected)
            self.assertEqual(len(calls),1)
        self.assertIsNone(ce.current_compilation())

    def test_NumericalContextKeepsChangingBindingsInSignature(self):
        expression=self._numerical_plan_fixture()
        tensor=np.arange(6.).reshape(2,3)
        scalar=np.array(.7)
        coefficients=[[[tensor,scalar]]]
        state=np.array([0,0]);permutations=np.array([[0,1]])
        frequencies=np.array([.1,.2])
        projector=np.array([0,1])
        def plan():
            return expression.prepare_evaluation([(state,permutations)],coefficients,frequencies,
                projection_bindings={'projector':projector}).signature
        expected=plan()
        with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()):
            self.assertEqual(plan(),expected)
            keys=[plan()]
            state[0]=1;keys.append(plan())
            permutations[:]=[[1,0]];keys.append(plan())
            frequencies[0]=.15;keys.append(plan())
            projector[0]=1;keys.append(plan())
            scalar[...]=1.4;keys.append(plan())
            self.assertEqual(len(set(keys)),len(keys))

    def test_NewContextRebindsChangedCoefficients(self):
        expression=self._numerical_plan_fixture()
        array=np.arange(4.)
        states=[(np.array([0,0]),np.array([[0,1]]))]
        def plan(coefficients):return expression.prepare_evaluation(states,coefficients,[.1,.2]).signature
        coefficients=[[[array]]]
        with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()):
            before=plan(coefficients)
            self.assertEqual(plan([[[array.copy()]]]),before)
            self.assertNotEqual(plan([[[array+1]]]),before)
        array+=10
        outside=plan(coefficients)
        self.assertNotEqual(outside,before)
        with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()):
            self.assertEqual(plan(coefficients),outside)
        array+=1
        self.assertNotEqual(plan(coefficients),outside)

    def test_CoefficientBindingsDoNotRetainDiscardedArrays(self):
        import gc,weakref
        expression=self._numerical_plan_fixture()
        states=[(np.array([0,0]),np.array([[0,1]]))]
        with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()):
            array=np.arange(12.).reshape(3,4).copy()
            reference=weakref.ref(array)
            coefficients=[[[array]]]
            plan=expression.prepare_evaluation(states,coefficients,[.1,.2])
            del plan,coefficients,array
            gc.collect()
            self.assertIsNone(reference())

    def test_CoefficientBindingScopeRestoresAfterNestedFailure(self):
        from unittest import mock
        expression=self._numerical_plan_fixture()
        array=np.arange(6.)
        states=[(np.array([0,0]),np.array([[0,1]]))]
        coefficients=[[[array]]]
        original=np.ascontiguousarray;calls=[]
        def observed(value,*args,**kwargs):
            if value is array:calls.append(1)
            return original(value,*args,**kwargs)
        with mock.patch.object(np,'ascontiguousarray',side_effect=observed):
            with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()):
                outer=ce.current_compilation()
                expression.prepare_evaluation(states,coefficients,[.1,.2])
                with self.assertRaisesRegex(RuntimeError,'inner'):
                    with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()):
                        expression.prepare_evaluation(states,coefficients,[.1,.2])
                        raise RuntimeError('inner')
                self.assertIs(ce.current_compilation(),outer)
                expression.prepare_evaluation(states,coefficients,[.1,.2])
                self.assertEqual(len(calls),2)
            expression.prepare_evaluation(states,coefficients,[.1,.2])
            self.assertEqual(len(calls),3)

    def test_SharedCoefficientBindingsKeepNativeSignature(self):
        tree=SharedMemoryArrayTree.create({'coefficients':[np.arange(6.).reshape(2,3)]})
        expression=self._numerical_plan_fixture()
        states=[(np.array([0,0]),np.array([[0,1]]))]
        try:
            expected=expression.prepare_evaluation(states,tree,[.1,.2]).signature
            with ce.compiled_evaluation_context(ce.AnalyticExpressionStore()):
                self.assertEqual(expression.prepare_evaluation(states,tree,[.1,.2]).signature,expected)
                self.assertEqual(expression.prepare_evaluation(states,tree,[.1,.2]).signature,expected)
        finally:
            tree.close();tree.unlink()


if __name__ == '__main__':
    unittest.main()
