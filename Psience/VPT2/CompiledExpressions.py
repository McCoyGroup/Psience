"""Persistent, bounded analytic compilation and native parallel contractions.

Only final expressions are compiled. Coefficients, frequencies, conditions,
state requests and numerical cutoffs remain inputs to the native evaluator.
The optional occupation projection is applied after every symbolic product.
"""
import collections
import contextlib
import contextvars
import dataclasses
import functools
import hashlib
import inspect
import itertools
import json
import marshal
import math
import os
from pathlib import Path
import pickle
import sys
import time
import uuid
import numpy as np
from McUtils.Parallelizers import Parallelizer, MultiprocessingParallelizer
from McUtils.Devutils import (atomic_write, atomic_write_json, file_hash, freeze, fingerprint,
                              implementation_signature, EntryStore)

from .Conditional import (ConditionalTransformationStack, SameStateCondition,
                          operator_sum, operator_product, ZeroTerm)
from .Analytic import ShiftedEnergyBaseline

__all__ = ['AnalyticExpressionStore', 'CompiledAnalyticExpression',
           'CompiledAnalyticEvaluation', 'AnalyticCorrectionCompiler']

_current_compilation = contextvars.ContextVar('analytic_compilation', default=None)


@contextlib.contextmanager
def compiled_evaluation_context(store, compilation_options=None, evaluation_options=None,
                                parallelizer=None):
    token = _current_compilation.set(None if store is None else
        (store, compilation_options or {}, evaluation_options or {}, parallelizer))
    try:
        yield
    finally:
        _current_compilation.reset(token)


def current_compilation():
    return _current_compilation.get()


def _atomic_json(path, value):
    atomic_write_json(path, value, sort_keys=True)


def _atomic_pickle(path, value):
    # Preserve shared native PolyPath nodes. JSON/NPZ object conversion expands
    # those graphs; ordinary pickle is the storage format, not a new algebra.
    atomic_write(path, lambda stream: pickle.dump(value, stream, protocol=5), fsync=False, keep_mode=False)


_file_hash = file_hash

# Evaluation caches and handles that are not part of an expression's definition.
_RUNTIME_ATTRS = ('logger', 'checkpoint', '_cascade_manager', '_changes', '_exprs', '_combined', '_eager',
                  '_op_maps', '_conditional_stacks', '_cache', '_grouped_terms')


def _freeze_solver(value, kind, rec):
    names = ['hamiltonian_expansion', 'allowed_terms', 'allowed_coefficients',
             'disallowed_coefficients', 'allowed_energy_changes',
             'intermediate_normalization', 'polynomial_representation', 'resonance_condition']
    return [kind, {name: rec(getattr(value, name)) for name in names}]


def _reject_parallelizer(value, kind, rec):
    raise TypeError('a parallelizer is runtime state, not an expression definition')


def _freeze_options():
    from McUtils.Parallelizers import SharedMemoryArrayTree
    from .Analytic import AnalyticPerturbationTheorySolver
    return dict(
        pre_handlers=[(lambda v: isinstance(v, SharedMemoryArrayTree), lambda v, rec: rec(v.tree))],
        handlers=[(lambda v: isinstance(v, AnalyticPerturbationTheorySolver), _freeze_solver),
                  (lambda v: isinstance(v, Parallelizer), _reject_parallelizer)],
        exclude_attrs=_RUNTIME_ATTRS
    )


def _freeze(value, seen=None):
    """Fingerprint definitions without including mutable evaluation caches."""
    return freeze(value, seen=seen, **_freeze_options())


def _fingerprint(value):
    return fingerprint(value, **_freeze_options())


def _implementation_signature():
    return implementation_signature('McUtils', 'Psience')


def _polynomial_degree(expression):
    """Bound operator letters from the actual DAG, independent of VPT order."""
    from .Analytic import PTTensorCoeffProductDAG, PTTensorCoeffProductSum
    def leaf(value):
        if isinstance(value, PTTensorCoeffProductSum):
            return max((sum(len(c)-2 for c in key) for key in value.terms), default=0)
        return 0
    if not isinstance(expression, PTTensorCoeffProductDAG):
        return leaf(expression)
    cache = {};pending = [(expression, False)]
    while pending:
        node, ready = pending.pop()
        if node in cache:continue
        kind, *args = node._node
        children = [arg for arg in args if isinstance(arg, PTTensorCoeffProductDAG)]
        if not ready:
            pending.append((node, True));pending.extend((child, False) for child in children)
        elif kind == 'leaf':cache[node] = leaf(args[0])
        elif kind == 'add':cache[node] = max(cache[child] for child in children)
        elif kind in {'mul_along','mul_simple','rmul_along','rmul_simple'}:
            cache[node] = sum(cache[arg] if isinstance(arg, PTTensorCoeffProductDAG) else leaf(arg)
                              for arg in args[:2])
        else:cache[node] = max((cache[child] for child in children), default=0)
    return cache[expression]


@dataclasses.dataclass(frozen=True)
class _OccupationProjection:
    before: int
    after: int
    background: int = 0

    def compiler(self):
        from .Analytic import PolyPath, PTEnergyChangeProductSum, _BoundedEvaluationCache
        from .Conditional import canonical_condition_key
        projection = self
        class Transform:
            def __init__(self):
                self.cache = _BoundedEvaluationCache(100000, 128*1024**2)
                self.nodes = self.cache.view('projection_nodes');self.axes = self.cache.view('projection_axes')
                self.memo = {};self.constants = {}
            def difference(self, value, rank):
                if np.isscalar(value):return np.longdouble(0)
                key = (value, rank)
                if key not in self.memo:
                    path = PolyPath.from_polynomial(value)
                    states = np.full((2, max(rank, path.ndim)), projection.background, dtype=int)
                    states[:,0] = [projection.before, projection.after]
                    values = np.asarray(path.evaluate_polynomial(states,node_cache=self.nodes,axis_cache=self.axes),dtype=np.longdouble)
                    if len(self.memo)>=10000:self.memo.pop(next(iter(self.memo)))
                    self.memo[key] = values[1]-values[0]
                return self.memo[key]
            def constant(self,value,rank):
                key = (float(value),rank)
                if key not in self.constants:
                    self.constants[key] = PolyPath.from_coeffs([np.ones(1)]*rank,prefactor=key[0],steps=0)
                return self.constants[key]
            def __call__(self,key,value):
                rank = len({i for c in key for i in c[2:]})
                if rank == 0:return None
                if isinstance(value,PTEnergyChangeProductSum):
                    rows = {}
                    for energy,numerator in value.terms.items():
                        if any(getattr(row[0],'_state_projection_factor',False) or row[0] not in {0,1} for row in energy):
                            raise ValueError('occupation projection requires occupation-independent harmonic denominators')
                        canonical,scale = canonical_condition_key(energy)
                        if canonical is None:raise ValueError('invalid harmonic denominator')
                        rows[canonical] = rows.get(canonical,np.longdouble(0))+scale*self.difference(numerator,rank)
                    rows = {energy:self.constant(weight,rank) for energy,weight in rows.items() if weight!=0}
                    return value.mutate(rows,reduced=True) if rows else None
                weight = self.difference(value,rank)
                return self.constant(weight,rank) if weight else None
        return Transform()


class _CompiledEntries(EntryStore):
    # compiled expressions live at folder/<signature>/manifest.json, aliases at folder/aliases/<signature>.json
    meta_name = 'manifest.json'
    shard = 0
    alias_field = 'compiled'
    stage_prefix = '.compile-'


class AnalyticExpressionStore:
    format_version = 1
    def __init__(self, folder=None, *, memory_cache_size=128, memory_cache_bytes=512*1024**2):
        self.folder = None if folder is None else Path(folder).expanduser().resolve()
        self.entries = None
        if self.folder is not None:
            self.folder.mkdir(parents=True,exist_ok=True)
            self.entries = _CompiledEntries(self.folder)
        self.memory_cache_size = max(0, int(memory_cache_size))
        self.memory_cache_bytes = None if memory_cache_bytes is None else max(0, int(memory_cache_bytes))
        self._memory = collections.OrderedDict()
        self._aliases = collections.OrderedDict()
        self.memory_bytes = 0
        self.timings = dict(requests=0, compilations=0, cache_hits=0, failures=0,
            compilation_seconds=0., cache_lookup_seconds=0., preparation_seconds=0., evaluation_seconds=0.)

    def statistics(self):
        return dict(self.timings, storage='memory' if self.folder is None else 'disk',
                    memory_entries=len(self._memory), memory_payload_bytes=self.memory_bytes)

    def clear(self):
        """Release reusable in-memory entries; never remove persistent cache files."""
        self._memory.clear(); self._aliases.clear(); self.memory_bytes=0

    def _remember(self, compiled):
        size = sum(part['bytes'] for part in compiled.manifest['parts'])
        key = compiled.manifest['signature']
        if self.memory_cache_size and (self.memory_cache_bytes is None or size <= self.memory_cache_bytes):
            previous = self._memory.pop(key, None)
            if previous is not None:self.memory_bytes -= sum(p['bytes'] for p in previous.manifest['parts'])
            self._memory[key]=compiled; self.memory_bytes+=size
            while len(self._memory)>self.memory_cache_size or (self.memory_cache_bytes is not None and self.memory_bytes>self.memory_cache_bytes):
                removed, value=self._memory.popitem(last=False)
                self.memory_bytes-=sum(p['bytes'] for p in value.manifest['parts'])
                for alias,target in list(self._aliases.items()):
                    if target==removed:del self._aliases[alias]
        return compiled

    def _get_memory(self, signature):
        value=self._memory.get(signature)
        if value is not None:self._memory.move_to_end(signature)
        return value

    def _signature(self,context,active_operators,projection,options):
        return _fingerprint(dict(format=self.format_version,implementation=_implementation_signature(),
            definition=context,active_operators=active_operators,projection=projection,options=options))

    def _find_alias(self,signature):
        if self.folder is None:
            target=self._aliases.get(signature)
            return None if target is None else self._get_memory(target)
        target=self.entries.resolve_alias(signature)
        return None if target is None else CompiledAnalyticExpression.load(self.folder/target)

    def compile(self, expression, *, context=None, active_operators=None, projection=None,
                max_terms=32768, max_bytes=64*1024**2, cache_size=100000,
                cache_bytes=512*1024**2, node_cache_size=8192, mode_dimension=None, progress=None,
                coalesce_combine_nodes=True):
        """Compile native expressions; early DAG combine-node coalescing is on by default.

        Set coalesce_combine_nodes=False (also accepted through evaluator
        compilation_options) to disable that optimization for diagnostics or
        benchmarking. Bounded final-part accumulation remains active.
        """
        return self.compile_factory(lambda:expression, context=expression if context is None else context,
            active_operators=active_operators,projection=projection,max_terms=max_terms,max_bytes=max_bytes,
            cache_size=cache_size,cache_bytes=cache_bytes,node_cache_size=node_cache_size,
            mode_dimension=mode_dimension,progress=progress,coalesce_combine_nodes=coalesce_combine_nodes)

    def compile_factory(self, factory, **options):
        """Compile or reuse an entry, accounting separately for lookup and actual construction."""
        started=time.perf_counter(); building=[]; completed=False
        self.timings['requests']+=1
        def build():
            building.append(time.perf_counter())
            return factory()
        try:
            result=self._compile_factory(build, **options)
            result.store=self; completed=True
            return result
        finally:
            elapsed=time.perf_counter()-started
            construction=0. if not building else time.perf_counter()-building[0]
            self.timings['compilation_seconds']+=construction
            self.timings['cache_lookup_seconds']+=max(0.,elapsed-construction)
            self.timings['compilations' if building else 'cache_hits']+=int(completed)
            self.timings['failures']+=int(not completed)

    def _compile_factory(self, factory, *, context, active_operators=None, projection=None,
                         max_terms=32768,max_bytes=64*1024**2,cache_size=100000,
                          cache_bytes=512*1024**2,node_cache_size=8192,mode_dimension=None,progress=None,
                          coalesce_combine_nodes=True):
        from .Analytic import (SqrtChangePoly,PTTensorCoeffProductDAG,PTTensorCoeffProductSum,
            PTTensorCoeffProductDAGEvaluationPlan,_BoundedEvaluationCache,_BalancedPolynomialAccumulator)
        from .BoundedStream import StreamBudget
        if not isinstance(coalesce_combine_nodes, (bool, np.bool_)):
            raise TypeError('coalesce_combine_nodes must be boolean')
        coalescing_enabled = bool(coalesce_combine_nodes)
        options = dict(max_terms=max_terms,max_bytes=max_bytes,cache_size=cache_size,
                       cache_bytes=cache_bytes,node_cache_size=node_cache_size,mode_dimension=mode_dimension,
                       coalesce_combine_nodes=coalescing_enabled)
        if max_terms<1 or (max_bytes is not None and max_bytes<1):raise ValueError('invalid compilation budget')
        signature = self._signature(context,active_operators,projection,options)
        directory = None if self.folder is None else self.folder/signature
        if self.folder is None:
            cached=self._get_memory(signature)
            if cached is not None:return cached
        elif self.entries.has(signature):return CompiledAnalyticExpression.load(directory)
        alias=self._find_alias(signature)
        if alias is None and active_operators is not None:
            # A cache compiled without pruning is valid for a narrower active
            # coefficient set. The reverse substitution would lose terms.
            alias=self._find_alias(self._signature(context,None,projection,options))
        if alias is not None:return alias
        started = time.perf_counter()
        wrapper = factory()
        expr = wrapper.expr
        if np.isscalar(expr):
            if expr!=0:raise ValueError('nonzero scalar expression needs a native identity polynomial')
            change = tuple(wrapper.change) if wrapper.change is not None else ()
            baseline = None;poly = PTTensorCoeffProductSum({})
        elif isinstance(expr,SqrtChangePoly):
            change = tuple(int(i) for i in expr.poly_change)
            baseline = tuple(int(i) for i in expr.shift_start);poly = expr.poly_obj
        else:
            change = tuple(wrapper.change) if wrapper.change is not None else ()
            baseline = None;poly = expr
        if not isinstance(poly,PTTensorCoeffProductSum):poly = PTTensorCoeffProductSum({():poly})
        if projection is not None and (change or baseline and any(baseline)):
            raise ValueError('final occupation projection requires a diagonal expression')
        if isinstance(poly,PTTensorCoeffProductDAG) and active_operators is not None:
            poly = poly.prune_operators(poly.operator_keys-set(active_operators))
        degree = _polynomial_degree(poly)
        maximum = len(change)+max(0,(degree-sum(abs(i) for i in change))//2)
        if mode_dimension is not None:
            if mode_dimension<1:raise ValueError('mode_dimension must be positive')
            maximum=min(maximum,mode_dimension)
        if isinstance(poly,PTTensorCoeffProductDAG):
            class CompilationPlan(PTTensorCoeffProductDAGEvaluationPlan):
                coalesce_combine_nodes = coalescing_enabled
            plan = CompilationPlan(poly,
                cache=_BoundedEvaluationCache(cache_size,cache_bytes),max_cached_terms=node_cache_size,
                max_mode_indices=maximum)
            iterator = plan.iter_terms()
            if progress is not None:
                plan.progress_callback=lambda p:progress(dict(stage='generating',**p.stats()))
        else:
            plan = None;iterator = poly.terms.items()
        transform = projection.compiler() if projection is not None else None
        budget = StreamBudget(dict(max_terms=max_terms,max_bytes=max_bytes))
        staging = contextlib.nullcontext(None) if self.folder is None else self.entries.stage(signature)
        memory_parts={}
        with staging as staged:
            stage = None if staged is None else Path(staged.path)
            pending={};parts=[]
            def flush(reason):
                if not pending:return
                compact = PTTensorCoeffProductSum({key:acc.total() for key,acc in pending.items()},
                    canonicalize=False,reduced=False).combine()
                groups = {}
                for key,value in compact.terms.items():
                    rank = len({i for coefficient in key for i in coefficient[2:]})
                    groups.setdefault(rank,{})[key] = value
                for rank,terms in sorted(groups.items()):
                    name = 'part-{:06d}-r{}.pkl'.format(len(parts),rank)
                    part = PTTensorCoeffProductSum(terms,canonicalize=False,reduced=True)
                    if stage is None:
                        payload=pickle.dumps(part,protocol=5)
                        part_hash=hashlib.sha256(payload).hexdigest();part_bytes=len(payload)
                        memory_parts[name]=part
                    else:
                        _atomic_pickle(stage/name,part)
                        part_hash=_file_hash(stage/name);part_bytes=(stage/name).stat().st_size
                    parts.append(dict(file=name,rank=rank,sha256=part_hash,bytes=part_bytes,
                        coefficient_keys=len(terms),energy_rows=sum(len(p.terms) if hasattr(p,'terms') else 1 for p in terms.values())))
                pending.clear();budget.flushed(reason)
            for key,value in iterator:
                if transform is not None:value = transform(key,value)
                if value is None:continue
                pending.setdefault(key,_BalancedPolynomialAccumulator()).append(value)
                reason = budget.incoming(pending)
                if reason:flush(reason)
                if progress is not None:progress(dict(stage='compiling',parts=len(parts),**budget.stats()))
            flush('end')
            manifest = dict(format=self.format_version,complete=True,signature=signature,parts=parts,
                change=change,baseline_shift=baseline,num_fixed=wrapper.num_fixed,
                projection=None if projection is None else dataclasses.asdict(projection),
                active_operators=None if active_operators is None else sorted(active_operators),
                operator_letter_bound=degree,maximum_mode_indices=maximum,options=options,
                seconds=time.perf_counter()-started,stream_stats=budget.stats(),
                plan_stats=None if plan is None else plan.stats(),implementation=_implementation_signature())
            if stage is not None:_atomic_json(stage/'manifest.json',manifest)
        compiled=(self._remember(CompiledAnalyticExpression(None,manifest,parts=memory_parts))
                  if directory is None else CompiledAnalyticExpression.load(directory))
        if projection is None:
            native=[wrapper.op if wrapper.op is not None else poly]
            from .ConditionalEvaluator import CorrectionView
            if isinstance(wrapper.op,CorrectionView) and not any(getattr(wrapper.op,name) is not None for name in
                    ['allowed_terms','allowed_coefficients','disallowed_coefficients','allowed_energy_changes']):
                native.append(wrapper.op.term)
            for op in native:
                key=self._signature(dict(operator=op,change=change,baseline=baseline),active_operators,None,options)
                if self.folder is None:
                    if signature in self._memory:
                        self._aliases[key]=signature
                        while len(self._aliases)>4*self.memory_cache_size:self._aliases.popitem(last=False)
                else:self.entries.set_alias(key,signature)
        return compiled


@dataclasses.dataclass(frozen=True)
class _TemplateSubset:
    rank: int
    def count(self,ndim,num_fixed,free_count):return int(free_count==self.rank-num_fixed)
    def iter_subsets(self,ndim,num_fixed,free_count):
        if free_count==self.rank-num_fixed:yield tuple(range(num_fixed,self.rank))


class CompiledAnalyticExpression:
    def __init__(self,folder,manifest,parts=None):
        self.folder=None if folder is None else Path(folder)
        self.manifest=manifest;self.parts=parts;self.store=None
    @classmethod
    def load(cls,folder):
        folder=Path(folder).resolve();manifest=json.loads((folder/'manifest.json').read_text())
        if not manifest.get('complete') or manifest.get('format')!=AnalyticExpressionStore.format_version:
            raise ValueError('incomplete or unsupported compiled expression')
        if manifest['implementation']!=_implementation_signature():raise ValueError('compiled expression source version changed')
        for part in manifest['parts']:
            if Path(part['file']).name!=part['file']:raise ValueError('invalid expression part path')
            if _file_hash(folder/part['file'])!=part['sha256']:raise ValueError('compiled expression checksum mismatch')
        return cls(folder,manifest)
    def prepare_evaluation(self,state_perms,coeffs,freqs,*,anchored_modes=(),parallelizer=None,
                           batch_size=512,task_size=32768,resume=True,use_shared_memory=None,**numerical_options):
        started=time.perf_counter()
        try:
            return CompiledAnalyticEvaluation(self,state_perms,coeffs,freqs,anchored_modes=anchored_modes,
                parallelizer=parallelizer,batch_size=batch_size,task_size=task_size,resume=resume,
                numerical_options=numerical_options,use_shared_memory=use_shared_memory)
        finally:
            if self.store is not None:self.store.timings['preparation_seconds']+=time.perf_counter()-started
    def evaluate(self,state_perms,coeffs,freqs,**options):
        return self.prepare_evaluation(state_perms,coeffs,freqs,**options).run()


class CompiledAnalyticEvaluation:
    """A reusable numerical task plan; the parallelizer never sees a symbolic DAG."""
    def __init__(self,compiled,state_perms,coeffs,freqs,*,anchored_modes=(),parallelizer=None,
                 batch_size=512,task_size=32768,resume=True,numerical_options=None,use_shared_memory=None):
        self.compiled=compiled
        self.state_perms=[(np.asarray(state,dtype=int),np.asarray(perms,dtype=int)) for state,perms in state_perms]
        self.coeffs=coeffs;self.freqs=np.asarray(freqs);self.anchored_modes=tuple(anchored_modes)
        cap=compiled.manifest['options'].get('mode_dimension')
        if cap is not None and len(self.freqs)>cap:
            raise ValueError('compiled expression was pruned for fewer modes; recompile for the larger basis')
        self.parallelizer=Parallelizer.lookup(parallelizer) if parallelizer is not None else None
        self.batch_size=int(batch_size);self.task_size=int(task_size);self.resume=bool(resume)
        if min(self.batch_size,self.task_size)<1:raise ValueError('batch/task sizes must be positive')
        self.numerical_options=dict(numerical_options or {})
        forbidden={'parallelizer','mode_subset_selector','evaluation_mode'} & self.numerical_options.keys()
        if forbidden:raise TypeError('compiled numerical options override scheduling: {}'.format(sorted(forbidden)))
        self.use_shared_memory=use_shared_memory is not False
        projection=compiled.manifest['projection']
        if projection is not None:
            if len(self.anchored_modes)!=1:raise ValueError('projected transition requires its one anchored physical mode')
            anchor=self.anchored_modes[0]
            for state,perms in self.state_perms:
                expected=np.full(len(self.freqs),projection['background'],dtype=int);expected[anchor]=projection['after']
                if not np.array_equal(state,expected) or not np.array_equal(perms,np.arange(len(self.freqs))[None]):
                    raise ValueError('occupation-specialized expression used outside its target-state domain')
        elif self.anchored_modes:raise ValueError('anchoring requires a matching final occupation projection')
        signature=_fingerprint(dict(expression=compiled.manifest['signature'],states=self.state_perms,
            coefficients=coeffs,frequencies=self.freqs,anchors=self.anchored_modes,numerical=self.numerical_options,
            batch_size=self.batch_size,task_size=self.task_size))
        self.folder=None if compiled.folder is None else compiled.folder/'evaluations'/signature
        if self.folder is not None:self.folder.mkdir(parents=True,exist_ok=True)
        tasks=[];ndim=len(self.freqs);num_fixed=compiled.manifest['num_fixed']
        for part in compiled.manifest['parts']:
            rank=part['rank'];free=rank-num_fixed
            if free<0:raise ValueError('coefficient aliases do not cover the fixed change')
            count=math.comb(ndim-num_fixed-len(self.anchored_modes),free-len(self.anchored_modes)) if free>=len(self.anchored_modes) and free<=ndim-num_fixed else 0
            for start in range(0,count,self.task_size):
                stop=min(start+self.task_size,count);index=len(tasks)
                tasks.append(dict(index=index,part=part,start=start,stop=stop,
                                  file='task-{:08d}.pkl'.format(index)))
        self.tasks=tasks;self.signature=signature
        if self.folder is not None:_atomic_json(self.folder/'plan.json',dict(signature=signature,expression_signature=compiled.manifest['signature'],
            kind='numerical-task-plan',task_count=len(tasks),subset_part_count=sum(t['stop']-t['start'] for t in tasks),
            batch_size=self.batch_size,task_size=self.task_size,tasks=tasks))

    def _context(self,coeffs=None):
        return dict(expression_folder=None if self.compiled.folder is None else str(self.compiled.folder),
            expression_parts=self.compiled.parts,manifest=self.compiled.manifest,
            evaluation_folder=None if self.folder is None else str(self.folder),signature=self.signature,state_perms=self.state_perms,
            coeffs=self.coeffs if coeffs is None else coeffs,freqs=self.freqs,anchors=self.anchored_modes,
            batch_size=self.batch_size,resume=self.resume,numerical_options=self.numerical_options)

    def run_task(self,task):
        if isinstance(task,int):task=self.tasks[task]
        return _run_compiled_tasks([task],self._context())[0]

    def run(self):
        started=time.perf_counter()
        try:
            manager = contextlib.nullcontext() if self.parallelizer is None else self.parallelizer
            with manager:return self._run()
        finally:
            if self.compiled.store is not None:self.compiled.store.timings['evaluation_seconds']+=time.perf_counter()-started

    def _run(self):
        backend=self.parallelizer
        nprocs=1 if backend is None else backend.nprocs
        bins=[[] for _ in range(nprocs)];loads=[0]*nprocs
        for task in sorted(self.tasks,key=lambda t:(t['stop']-t['start'])*max(1,t['part']['energy_rows']),reverse=True):
            index=min(range(nprocs),key=lambda i:loads[i]);bins[index].append(task)
            loads[index]+=(task['stop']-task['start'])*max(1,task['part']['energy_rows'])
        owner=None
        try:
            from McUtils.Parallelizers import SharedMemoryArrayTree
            if isinstance(backend,MultiprocessingParallelizer) and self.use_shared_memory and not isinstance(self.coeffs,SharedMemoryArrayTree):
                from .Analytic import PerturbationTheoryExpressionEvaluator
                owner=SharedMemoryArrayTree.create(
                    PerturbationTheoryExpressionEvaluator._normalize_coefficient_expansions(self.coeffs))
            context=self._context(owner)
            if backend is None:
                records=_run_compiled_tasks(self.tasks,context)
            else:
                records=backend.apply(_map_compiled_bins,bins,context,cleanup=False)
            records=sorted(records,key=lambda r:r['index'])
            if [r['index'] for r in records]!=list(range(len(self.tasks))):raise ValueError('incomplete or duplicated compiled task coverage')
            outputs=[np.zeros((len(self.coeffs),len(perms)),dtype=np.result_type(
                *[r['values'][position].dtype for r in records],float))
                for position,(_,perms) in enumerate(self.state_perms)]
            for position,out in enumerate(outputs):
                for index in np.ndindex(out.shape):
                    values=[r['values'][position][index] for r in records]
                    if any(np.iscomplexobj(v) for v in values):out[index]=sum(values)
                    else:out[index]=math.fsum(values)
            if self.folder is not None:_atomic_json(self.folder/'result.json',dict(complete=True,signature=self.signature,task_count=len(records),
                subset_part_count=sum(t['stop']-t['start'] for t in self.tasks),
                resumed_tasks=sum(int(r.get('resumed',False)) for r in records)))
            return outputs
        finally:
            if owner is not None:
                from .Analytic import PerturbationTheoryExpressionEvaluator
                from .IndexedEvaluator import IndexedBlockEvaluator
                IndexedBlockEvaluator.clear_plan_cache()
                PerturbationTheoryExpressionEvaluator._cached_expansion=None
                owner.dispose()


def _map_compiled_bins(bins,context,parallelizer=None):
    return parallelizer.map(_run_assigned_compiled_bins,bins,extra_args=(context,),vectorized=True,aggregate=True)


def _run_assigned_compiled_bins(bins,context):
    return _run_compiled_tasks([task for group in bins for task in group],context)


def _run_compiled_tasks(tasks,context):
    from .Analytic import PerturbationTheoryExpressionEvaluator as Expression
    from .IndexedEvaluator import IndexedBlockEvaluator
    from McUtils.Parallelizers import SharedMemoryArrayTree
    coefficients=context['coeffs']
    if not isinstance(coefficients,SharedMemoryArrayTree):coefficients=Expression._normalize_coefficient_expansions(coefficients)
    manifest=context['manifest'];num_fixed=manifest['num_fixed'];ndim=len(context['freqs'])
    anchors=context['anchors'];available=tuple(i for i in range(num_fixed,ndim) if i not in anchors)
    cache=collections.OrderedDict();records=[]
    for task in tasks:
        target=None if context['evaluation_folder'] is None else Path(context['evaluation_folder'])/task['file']
        if target is not None and context['resume'] and target.exists():
            with target.open('rb') as stream:record=pickle.load(stream)
            if record['signature']!=context['signature'] or record['index']!=task['index']:raise ValueError('compiled result context mismatch')
            if record['start']!=task['start'] or not task['start']<=record['stop']<=task['stop']:
                raise ValueError('compiled result range mismatch')
            if record['complete']:
                if record['stop']!=task['stop']:raise ValueError('completed compiled task has a missing range')
                record=dict(record,resumed=True);records.append(record);continue
        else:record=None
        part=task['part'];key=part['file']
        if key not in cache:
            if context.get('expression_parts') is not None:
                cache[key]=context['expression_parts'][key]
            else:
                path=Path(context['expression_folder'])/key
                if _file_hash(path)!=part['sha256']:raise ValueError('compiled part changed in worker')
                with path.open('rb') as stream:cache[key]=pickle.load(stream)
            while len(cache)>2:cache.popitem(last=False)
        poly=cache[key];rank=part['rank'];free=rank-num_fixed
        batches=[];position=task['start']
        if record is not None:
            position=record['stop'];batches=record['batches']
            cursor=task['start']
            for batch in batches:
                if batch['start']!=cursor or not cursor<batch['stop']<=task['stop']:
                    raise ValueError('invalid partial compiled task coverage')
                cursor=batch['stop']
            if cursor!=position:raise ValueError('partial task stop does not match its batches')
        iterator=itertools.islice(itertools.combinations(available,free-len(anchors)),position,task['stop'])
        def total_values():
            arrays=[np.zeros((len(coefficients),len(perms))) for _,perms in context['state_perms']]
            for route,array in enumerate(arrays):
                for index in np.ndindex(array.shape):
                    array[index]=math.fsum(batch['values'][route][index] for batch in batches)
            return arrays
        values=total_values()
        for extra_batch in iter(lambda:list(itertools.islice(iterator,context['batch_size'])),[]):
            subsets=[anchors+tuple(extra) for extra in extra_batch]
            requests=[]
            for state,perms in context['state_perms']:
                expanded=[]
                for subset in subsets:
                    order=tuple(range(num_fixed))+subset+tuple(i for i in range(num_fixed,ndim) if i not in subset)
                    expanded.extend(perms[:,order])
                requests.append([state,np.asarray(expanded,dtype=int)])
            options=dict(evaluation_mode='materialized',evaluation_backend='indexed',logger=False,log_scaled=False,
                path_cache_size=10000,path_cache_bytes=64*1024**2,
                state_permutation_cache_size=10000,state_permutation_cache_bytes=16*1024**2)
            options.update(context['numerical_options'])
            result=Expression.evaluate_polynomial_expression(requests,coefficients,context['freqs'],poly,
                manifest['change'],manifest['baseline_shift'],num_fixed,
                mode_subset_selector=_TemplateSubset(rank),**options)
            batch_values=[]
            for array,(_,perms) in zip(result,context['state_perms']):
                array=np.asarray(array).reshape(len(coefficients),len(subsets),len(perms))
                value=np.zeros((len(coefficients),len(perms)))
                for index in np.ndindex(value.shape):value[index]=math.fsum(array[index[0],:,index[1]])
                if not np.isfinite(value).all():raise FloatingPointError('nonfinite compiled contraction')
                batch_values.append(value)
            position+=len(subsets);batches.append(dict(start=position-len(subsets),stop=position,values=batch_values))
            values=total_values()
            if target is not None:_atomic_pickle(target,dict(complete=False,signature=context['signature'],index=task['index'],
                start=task['start'],stop=position,values=values,batches=batches))
        if position!=task['stop']:raise ValueError('compiled subset range incomplete')
        record=dict(complete=True,signature=context['signature'],index=task['index'],start=task['start'],stop=position,
                    values=values,batches=batches,resumed=record is not None)
        if target is not None:_atomic_pickle(target,record)
        records.append(record)
    IndexedBlockEvaluator.clear_plan_cache()
    return records


class AnalyticCorrectionCompiler:
    """Correction selection and final projection belong to the evaluator layer."""
    def __init__(self,evaluator,folder=None):
        self.evaluator=evaluator
        self.store=evaluator.expression_store if folder is None else AnalyticExpressionStore(folder)

    def _generator(self,correction,order,energy_method):
        from .Conditional import SameStateCondition
        if correction=='energy' and energy_method=='auto':
            energy_method='bloch' if self.evaluator.solver.resonance_condition is None or type(self.evaluator.solver.resonance_condition) is SameStateCondition else 'native'
        if correction=='energy' and energy_method=='bloch':
            if self.evaluator.solver.resonance_condition is not None and type(self.evaluator.solver.resonance_condition) is not SameStateCondition:
                raise ValueError('scalar Bloch energies must not replace a nontrivial resonance condition')
            from .Analytic import AnalyticPerturbationTheorySolver
            source=self.evaluator.solver
            prepared=self.evaluator._prep_expansions(self.evaluator.expansions)
            # Remove identically zero Hamiltonian orders before enumerating
            # operator paths, rather than only pruning the resulting tensor DAG.
            expansion=[term if n==0 or n<len(prepared) and any(np.any(c) for c in prepared[n])
                       else ZeroTerm(logger=False) for n,term in enumerate(source.hamiltonian_expansion)]
            solver=AnalyticPerturbationTheorySolver(expansion,logger=source.logger,
                polynomial_representation=source.polynomial_representation,
                resonance_condition=source.resonance_condition,
                allowed_terms=source.allowed_terms,allowed_coefficients=source.allowed_coefficients,
                disallowed_coefficients=source.disallowed_coefficients,
                allowed_energy_changes=source.allowed_energy_changes,
                intermediate_normalization=source.intermediate_normalization)
            terms=AnalyticBlochSeries(solver).hermitian_hamiltonian(order)
            return lambda n:terms[n]
        names={'energy':'energy_correction','overlap':'overlap_correction',
               'wavefunction':'wavefunction_correction','full_wavefunction':'full_wavefunction_correction',
               'operator':'operator_correction'}
        if callable(correction):return correction
        if correction not in names:raise ValueError('unknown analytic correction '+str(correction))
        return getattr(self.evaluator.solver,names[correction])

    def precompile(self,correction='energy',order=None,orders=None,changes=((),),energy_method='native',projection=None,
                   compilation_options=None):
        from .Analytic import PerturbationTheoryExpressionEvaluator
        evaluator=self.evaluator
        if order is None:order=len(evaluator.expansions)-1
        if orders is None:orders=range(0,order+1,2 if correction=='energy' else 1)
        prepared=evaluator._prep_expansions(evaluator.expansions)
        active=tuple((n,k) for n,row in enumerate(prepared) for k,c in enumerate(row) if np.any(c))
        # An operator correction can introduce coefficients absent from the
        # Hamiltonian payload. Their absence here is not evidence of a zero.
        if correction not in {'energy','overlap','wavefunction','full_wavefunction'}:
            active=None
        options=dict(evaluator.compilation_options);options.update(compilation_options or {})
        options.setdefault('mode_dimension',len(evaluator.freqs))
        source=dict(solver=evaluator.solver,correction=correction,energy_method=energy_method)
        generator=[None]
        def expression(n,change):
            if generator[0] is None:generator[0]=self._generator(correction,order,energy_method)
            term=generator[0](n)
            return PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(change),change)
        return {(n,tuple(change)):self.store.compile_factory(lambda n=n,c=tuple(change):expression(n,c),
                    context=dict(source,order=n,change=tuple(change)),active_operators=active,projection=projection,**options)
                for n in orders for change in changes}

    def _frequency_request(self,initial,final):
        from .Conditional import SameStateCondition
        initial=np.asarray(initial,dtype=int);final=np.asarray(final,dtype=int)
        if initial.shape!=final.shape or initial.ndim!=1 or np.any(initial<0) or np.any(final<0):
            raise ValueError('frequency states must be nonnegative vectors of the same dimension')
        if len(initial)!=len(self.evaluator.freqs):raise ValueError('state dimension differs from expansion')
        changed=np.flatnonzero(initial!=final);projection=None;anchors=()
        if len(changed)==1 and (self.evaluator.solver.resonance_condition is None or type(self.evaluator.solver.resonance_condition) is SameStateCondition):
            if len(changed)==1:
                anchor=int(changed[0]);rest=np.delete(initial,anchor)
                if len(rest)==0 or np.all(rest==rest[0]):
                    background=int(rest[0]) if len(rest) else 0
                    projection=_OccupationProjection(int(initial[anchor]),int(final[anchor]),background)
                    anchors=(anchor,)
        return initial,final,projection,anchors

    def precompile_frequencies(self,initial,final,order=None,energy_method='auto',specialize=True,**options):
        _,_,projection,_=self._frequency_request(initial,final)
        if not specialize:projection=None
        return self.precompile('energy',order=order,energy_method=energy_method,projection=projection,**options)

    def frequencies(self,initial,final,order=None,energy_method='auto',parallelizer=None,specialize=True,**options):
        initial,final,projection,anchors=self._frequency_request(initial,final)
        numerical=dict(self.evaluator.evaluation_options);numerical.update(options)
        if not specialize or numerical.get('degenerate_changes') is not None or numerical.get('only_degenerate_terms'):
            projection=None;anchors=()
        compiled=self.precompile('energy',order=order,energy_method=energy_method,projection=projection)
        prepared=self.evaluator._prep_expansions(self.evaluator.expansions)
        identity=np.arange(len(initial))[None]
        state_perms=[(final,identity)] if projection is not None else [(initial,identity),(final,identity)]
        backend=self.evaluator.parallelizer if parallelizer is None else parallelizer
        results=[]
        for key,expression in compiled.items():
            values=expression.evaluate(state_perms,[prepared],self.evaluator.freqs,anchored_modes=anchors,
                parallelizer=backend,**numerical)
            results.append(float(values[0][0,0]) if projection is not None else float(values[1][0,0]-values[0][0,0]))
        return np.asarray(results)


# --------------------------------------------------------------------------------------------- #
#  Bloch series (formerly AnalyticBloch.py)
# Condition-based analytic Bloch wave operator, using native projected terms.
#
# Omega_0=I, D Omega_n=0 for n>0. With V=H-H0:
#   T_n = sum_{r=1}^n H_r Omega_{n-r}
#   K_n = D T_n
#   Omega_n = R_Q (T_n - sum_{r=1}^{n-1} Omega_{n-r} K_r).
# R_Q is the native conditional homological inverse; it cannot propagate through
# the retained relation. No perturbation order is singled out.
# --------------------------------------------------------------------------------------------- #

class AnalyticBlochSeries:
    def __init__(self,solver,condition=None):
        condition=SameStateCondition() if condition is None else condition
        if not condition.equivalence:
            raise ValueError('Bloch model spaces require an equivalence relation')
        self.projector=ConditionalTransformationStack(solver,condition)
        self.condition=self.projector.condition
        self.h0=solver.hamiltonian_expansion[0]
        self.expansion=solver.hamiltonian_expansion
        self.identity,self.zero=self.projector.identity,self.projector.zero
        self.wave_operator=[self.identity]
        self.effective_hamiltonian=[self.h0]

    def extend(self,order):
        if order>=len(self.expansion):raise ValueError('Hamiltonian expansion too short')
        while len(self.wave_operator)<=order:
            n=len(self.wave_operator)
            driving=operator_sum(*(operator_product(self.expansion[r],self.wave_operator[n-r])
                                  for r in range(1,n+1)))
            correction=self.projector.restrict(driving)
            feedback=operator_sum(*(operator_product(self.wave_operator[n-r],self.effective_hamiltonian[r])
                                   for r in range(1,n)))
            self.effective_hamiltonian.append(correction)
            self.wave_operator.append(self.projector.restrict(operator_sum(driving,-feedback),
                                                            retain=False,resolvent=True))
        return self

    def metric(self,order):
        self.extend(order)
        result=[self.identity]
        for n in range(1,order+1):
            result.append(self.projector.restrict(operator_sum(*(
                operator_product(ShiftedEnergyBaseline.lookup(self.wave_operator[i]),self.wave_operator[n-i])
                for i in range(n+1)))))
        return result

    def _product_series(self,left,right,order):
        return [operator_sum(*(operator_product(left[i],right[n-i])
                    for i in range(n+1) if i<len(left) and n-i<len(right))) for n in range(order+1)]

    def _metric_power(self,order,power):
        delta=self.metric(order);delta[0]=self.zero
        result=[self.identity]+[self.zero]*order
        product=[self.identity]+[self.zero]*order;coefficient=1.
        for exponent in range(1,order+1):
            product=self._product_series(product,delta,order)
            coefficient*= (power-exponent+1)/exponent
            for n in range(1,order+1):result[n]=operator_sum(result[n],coefficient*product[n])
        return result

    def normalized_wavefunctions(self,order):
        self.extend(order)
        return self._product_series(self.wave_operator,self._metric_power(order,-.5),order)

    def hermitian_hamiltonian(self,order):
        self.extend(order)
        if isinstance(self.condition,SameStateCondition):
            # Each retained block is a scalar, so the metric similarity cancels.
            return self.effective_hamiltonian[:order+1]
        return self._product_series(self._metric_power(order,.5),
            self._product_series(self.effective_hamiltonian,self._metric_power(order,-.5),order),order)
