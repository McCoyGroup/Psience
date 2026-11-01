"""Mapped state frontiers and replayable, bounded coupling transformations."""
from pathlib import Path
import itertools, tempfile, weakref
import numpy as np
from McUtils.Numputils import external_unique, SparseArray
from McUtils.Combinatorics import SymmetricGroupGenerator
from .StateSpaces import BasisStateSpace, SelectionRuleStateSpace, BraKetSpace

__all__ = ['StateSpaceDiskStore', 'DiskBasisStateSpace', 'DiskSelectionRuleStateSpace',
           'QuantaStateSpace']


def _cleanup_store(mappings, temporary):
    # NFS must see all mappings closed before unlinking their files. NumPy and
    # SciPy may retain ndarray views after the original memmap object is gone.
    for mapping in mappings:
        if not mapping.closed:
            mapping.close()
    temporary.cleanup()


class StateSpaceDiskStore:
    def __init__(self, directory=None, block_size=65536):
        if block_size < 1:
            raise ValueError('state_space_block_size must be positive')
        if directory is not None:
            Path(directory).mkdir(parents=True, exist_ok=True)
        self._temporary = tempfile.TemporaryDirectory(prefix='coupled_states_', dir=directory)
        self.directory = Path(self._temporary.name)
        self.block_size = int(block_size)
        self._counter = itertools.count()
        self._mappings = []
        self._finalize = weakref.finalize(self,_cleanup_store,self._mappings,self._temporary)
        self.metrics = dict(peak_source_states=0, peak_candidate_pairs=0, generated_pairs=0)

    def path(self, label):
        return self.directory / ('{}_{}.npy'.format(label, next(self._counter)))

    def own(self, array):
        if isinstance(array, np.memmap):
            self._mappings.append(array._mmap)
            array._state_space_store = self
        return array

    def array(self, label, dtype, shape):
        return self.own(np.lib.format.open_memmap(self.path(label), mode='w+', dtype=dtype, shape=shape))

    def unique(self, arrays, width=1, sort=True):
        def chunks():
            yield np.empty((0,) if width == 1 else (0, width), dtype=np.int64)
            yield from arrays
        return self.own(external_unique(chunks(), sort=sort, directory=self.directory,
                                       chunk_size=self.block_size, output_file=self.path('unique')))

    def frontier(self, basis, arrays, sort=True):
        return DiskBasisStateSpace(basis, self.unique(arrays,sort=sort), self,
                                  canonical=sort,unique=True)

    def coerce(self, space):
        if isinstance(space, DiskSelectionRuleStateSpace):
            return space.to_single()
        if isinstance(space, DiskBasisStateSpace) and space.store is self:
            return space
        return self.frontier(space.basis, index_blocks(space, self.block_size))

    def close(self):
        self._finalize()


def index_blocks(space, size):
    if hasattr(space, 'iter_index_blocks'):
        for values in space.iter_index_blocks():
            for start in range(0, len(values), size):
                yield values[start:start + size]
    else:
        values = space.indices
        for start in range(0, len(values), size):
            yield values[start:start + size]


class QuantaStateSpace(BasisStateSpace):
    """An implicit fixed-total-quanta prefilter; never enumerates a manifold."""
    def __init__(self, basis, quanta):
        super().__init__(basis, [], mode='indices')
        self.quanta = int(quanta)
        group = SymmetricGroupGenerator(basis.ndim)
        self.start = sum(group.num_terms(q)[0] for q in range(self.quanta))
        self.stop = self.start + group.num_terms(self.quanta)[0]

    def __len__(self):
        return int(self.stop - self.start)

    def contains_indices(self, values):
        if not hasattr(self.basis.indexer,'symmetric_group'):
            return np.sum(self.basis.unravel_state_inds(values),axis=1) == self.quanta
        return (values >= self.start) & (values < self.stop)

    @property
    def indices(self):
        raise RuntimeError('implicit quanta filters must be applied by block membership')


class DiskBasisStateSpace(BasisStateSpace):
    """A sorted unique rank frontier whose canonical data stays mapped."""
    def __init__(self, basis, indices, store, canonical=True, unique=None):
        super().__init__(basis, [], mode='indices')
        self.store = store
        # Bypass the ordinary constructor's dtype conversion/copy and cache.
        self._init_states = self._indices = indices
        self._init_state_types = self.StateSpaceSpec.Indices
        self.keep_excitations = False
        self._sorted = True
        self._canonical = canonical
        self._unique = canonical if unique is None else unique
        self._lookup = None

    def __len__(self):
        return len(self._indices)

    def __repr__(self):
        return '{}(nstates={}, ndim={}, mapped=True)'.format(type(self).__name__, len(self), self.ndim)

    def is_unique(self, *args, **kwargs):
        return self._unique

    def is_sorted(self, *args, **kwargs):
        return self._canonical

    def as_indices(self, *args, **kwargs):
        return self._indices

    def as_unique_indices(self, *args, **kwargs):
        return self.take_unique(*args,**kwargs).indices

    @property
    def unique_indices(self):
        return self.take_unique().indices

    def iter_index_blocks(self):
        for start in range(0, len(self), self.store.block_size):
            yield self._indices[start:start + self.store.block_size]

    def iter_excitations(self):
        for values in self.iter_index_blocks():
            yield self.basis.unravel_state_inds(values)

    def quanta_for_positions(self, positions):
        result = self.store.array('quanta',np.int64,(len(positions),))
        for start in range(0,len(positions),self.store.block_size):
            selected = positions[start:start+self.store.block_size]
            result[start:start+len(selected)] = np.sum(
                self.basis.unravel_state_inds(self.indices[selected]),axis=1)
        return result

    def as_excitations(self, *args, **kwargs):
        if len(self) > self.store.block_size:
            raise RuntimeError('use iter_excitations() for a disk-backed state space')
        return self.basis.unravel_state_inds(self._indices)

    @property
    def excitations(self):
        return self.as_excitations()

    @property
    def indexer(self):
        # An explicit arange(N) is unnecessary for this backend.
        raise RuntimeError('mapped frontiers use block search instead of a resident indexer')

    @indexer.setter
    def indexer(self, value):
        self._indexer = value

    def find(self, values, check=True, minimal_dtype=False, dtype=None, missing_val='raise'):
        if isinstance(values, BasisStateSpace):
            values = values.indices
        values = np.asarray(values)
        if values.ndim > 1:
            values = self.basis.ravel_state_inds(values)
        scalar = values.ndim == 0
        values = values.reshape(-1)
        searched = self._indices
        if not self._canonical:
            if self._lookup is None:
                def records():
                    offset = 0
                    for block in self.iter_index_blocks():
                        yield np.column_stack((block,np.arange(offset,offset+len(block))))
                        offset += len(block)
                self._lookup = self.store.unique(records(),width=2)
            searched = self._lookup[:,0]
        result = np.empty(len(values), dtype=np.int64)
        for start in range(0, len(values), self.store.block_size):
            query = values[start:start + self.store.block_size]
            locations = np.searchsorted(searched, query)
            valid = locations < len(searched)
            if len(self):
                valid &= searched[np.minimum(locations, len(searched) - 1)] == query
            if check and missing_val == 'raise' and not np.all(valid):
                raise IndexError('state not present in the mapped frontier')
            fill = -1 if missing_val == 'raise' else missing_val
            if self._lookup is not None and len(searched):
                locations = self._lookup[np.minimum(locations,len(searched)-1),1]
            result[start:start + len(query)] = np.where(valid, locations, fill)
        return result[0] if scalar else result

    def _contains(self, values):
        return self.find(values, missing_val=-1) >= 0

    def take_unique(self, sort=False, **kwargs):
        if self._unique and (not sort or self._canonical):
            return self
        return self.store.frontier(self.basis,self.iter_index_blocks(),sort=sort)

    def as_sorted(self, *args, **kwargs):
        return self.take_unique(sort=True)

    def to_single(self, *args, **kwargs):
        return self

    def union(self, other, sort=False, **kwargs):
        if self.basis != other.basis:
            raise ValueError('different state-space bases')
        return self.store.frontier(self.basis, itertools.chain(
            self.iter_index_blocks(), index_blocks(other, self.store.block_size)),sort=sort)

    def _select_set(self, other, invert):
        if isinstance(other, QuantaStateSpace):
            contains = other.contains_indices
        else:
            other = self.store.coerce(other)
            contains = other._contains
        def blocks():
            for block in self.iter_index_blocks():
                mask = contains(block)
                yield block[~mask if invert else mask]
        return self.store.frontier(self.basis, blocks())

    def intersection(self, other, *args, **kwargs):
        return self._select_set(other, False)

    def difference(self, other, *args, **kwargs):
        return self._select_set(other, True)

    def take_subspace(self, selection, *args, **kwargs):
        # Subspaces used as bra/ket blocks preserve the requested ordering.
        if not isinstance(selection,slice):
            selection = np.atleast_1d(selection)
            if selection.dtype.kind == 'b':
                selection = np.flatnonzero(selection)
        if isinstance(selection, slice):
            count = len(range(*selection.indices(len(self))))
        else:
            count = len(np.atleast_1d(selection))
        if count > self.store.block_size:
            output = self.store.array('subspace',np.int64,(count,))
            if isinstance(selection,slice):
                first,stop,step = selection.indices(len(self))
            ordered = self._canonical
            previous = None
            for start in range(0,count,self.store.block_size):
                stop_at = min(start+self.store.block_size,count)
                positions = (np.arange(first+start*step,first+stop_at*step,step)
                             if isinstance(selection,slice) else selection[start:stop_at])
                block = self.indices[positions]
                if len(block):
                    ordered &= (previous is None or previous < block[0]) and np.all(block[1:] > block[:-1])
                    previous = block[-1]
                output[start:stop_at] = block
            return type(self)(self.basis,output,self.store,canonical=bool(ordered))
        return BasisStateSpace(self.basis, self._indices[selection], mode='indices')

    def get_representation_brakets(self, other=None, selection_rules=None, return_filter=False, **opts):
        if opts.get('freq_threshold') is not None:
            raise NotImplementedError('frequency-threshold braket filtering')
        rules = [()] if selection_rules is not None and len(selection_rules)==0 else selection_rules
        graph = DiskSelectionRuleStateSpace(self,rules,filter_space=self if other is None else other)
        return (graph,None) if return_filter else graph

    def apply_selection_rules(self, selection_rules, target_dimensions=None, filter_space=None,
                              parallelizer=None, logger=None, iterations=1, **opts):
        if iterations != 1:
            raise NotImplementedError('selection-rule iteration count must be one')
        return DiskSelectionRuleStateSpace(self, selection_rules, filter_space=filter_space,
                                          target_dimensions=target_dimensions)


class _StoredCouplings:
    """Mapped (outer representative, bra, ket) records for explicit graphs."""
    def __init__(self, records):
        self.records = records


class DiskSelectionRuleStateSpace(SelectionRuleStateSpace):
    """A union of replayable coupling recipes, not resident adjacency lists."""
    def __init__(self, source, rules=None, *, filter_space=None, target_dimensions=None, recipes=None):
        self.store = source.store
        self._base_space = source
        self.sel_rules = rules
        self.changes = None
        self.recipes = ([(source, rules, filter_space, target_dimensions, (), ())]
                        if recipes is None else recipes)
        self._frontier = None

    @property
    def basis(self):
        return self._base_space.basis

    @property
    def representative_space(self):
        return self._base_space

    @property
    def full_basis(self):
        return None

    def __repr__(self):
        return '{}(recipes={}, mapped=True)'.format(type(self).__name__, len(self.recipes))

    def __len__(self):
        return len(self.to_single())

    @property
    def indices(self):
        return self.to_single().indices

    def iter_index_blocks(self):
        yield from self.to_single().iter_index_blocks()

    @classmethod
    def from_resident(cls, graph, store):
        """Preserve explicit polyad edges, including nested graph ownership."""
        def records(node, owner=None):
            for bra, child in zip(node.representative_space.indices, node.spaces.flat):
                root = bra if owner is None else owner
                targets = child.representative_space if isinstance(child, SelectionRuleStateSpace) else child
                for block in index_blocks(targets, store.block_size):
                    if np.any(block > np.iinfo(np.int64).max):
                        raise OverflowError('disk state ranks must fit signed int64')
                    block = block.astype(np.int64, copy=False)
                    yield np.column_stack((np.full(len(block), root, dtype=np.int64),
                                           np.full(len(block), bra, dtype=np.int64), block))
                if isinstance(child, SelectionRuleStateSpace):
                    yield from records(child, root)
        source = store.coerce(graph.representative_space)
        saved = _StoredCouplings(store.unique(records(graph), width=3))
        return cls(source, recipes=[(source, saved, None, None, (), ())])

    def iter_pair_blocks(self):
        size = self.store.block_size
        for source, rules, filtered, dimensions, transforms, excluded in self.recipes:
            source_size = min(64, size)
            if filtered is not None:
                filtered = self.store.coerce(filtered)
            if isinstance(rules, _StoredCouplings):
                for start in range(0, len(rules.records), size):
                    block = rules.records[start:start + size]
                    block = block[source._contains(block[:, 0])]
                    pairs = block[:, 1:]
                    valid = np.ones(len(pairs), dtype=bool)
                    if filtered is not None:
                        valid &= filtered._contains(pairs[:, 1])
                    if excluded:
                        delta = (self.basis.unravel_state_inds(pairs[:, 1])
                                 - self.basis.unravel_state_inds(pairs[:, 0]))
                        for rule in excluded:
                            valid &= ~np.all(delta == rule, axis=1)
                    for transform in transforms:
                        kept = transform(BasisStateSpace(self.basis, pairs[valid, 1], mode='indices')).indices
                        valid &= np.isin(pairs[:, 1], kept)
                    if np.any(valid):
                        yield pairs[valid]
                continue
            if rules is None:
                targets = source if filtered is None else filtered
                for rank_block in source.iter_index_blocks():
                    for rank in rank_block:
                        for target_block in targets.iter_index_blocks():
                            selected = target_block
                            if excluded:
                                delta = self.basis.unravel_state_inds(selected)-self.basis.unravel_state_inds(np.array([rank]))
                                valid = np.ones(len(selected),dtype=bool)
                                for rule in excluded:
                                    valid &= ~np.all(delta==rule,axis=1)
                                selected = selected[valid]
                            for transform in transforms:
                                selected = transform(BasisStateSpace(self.basis,selected,mode='indices')).indices
                            if len(selected):
                                yield np.column_stack((np.full(len(selected),rank,dtype=np.int64),selected))
                continue
            for rank_block in source.iter_index_blocks():
                for start in range(0, len(rank_block), source_size):
                    ranks = rank_block[start:start + source_size]
                    original = self.basis.unravel_state_inds(ranks).astype(np.int64,copy=False)
                    self.store.metrics['peak_source_states'] = max(
                        self.store.metrics['peak_source_states'], len(ranks))
                    exc = original if dimensions is None else original[:, dimensions]
                    group = SymmetricGroupGenerator(exc.shape[1])
                    for rows, outputs in group.iter_permutation_rule_direct_sum(exc, rules, block_size=size):
                        self.store.metrics['peak_candidate_pairs'] = max(
                            self.store.metrics['peak_candidate_pairs'], len(outputs))
                        self.store.metrics['generated_pairs'] += len(outputs)
                        if dimensions is not None:
                            restored = original[rows].copy()
                            restored[:, dimensions] = outputs
                            outputs = restored
                        if excluded:
                            delta = outputs - original[rows]
                            mask = np.ones(len(outputs), dtype=bool)
                            for rule in excluded:
                                mask &= ~np.all(delta == rule, axis=1)
                            rows, outputs = rows[mask], outputs[mask]
                        targets = self.basis.ravel_state_inds(outputs)
                        if np.any(targets > np.iinfo(np.int64).max):
                            raise OverflowError('disk state ranks must fit signed int64')
                        targets = targets.astype(np.int64,copy=False)
                        valid = np.ones(len(targets), dtype=bool)
                        if filtered is not None:
                            valid &= filtered._contains(targets)
                        for transform in transforms:
                            block = BasisStateSpace(self.basis, targets[valid], mode='indices')
                            kept = transform(block).indices
                            valid &= np.isin(targets, kept)
                        if np.any(valid):
                            yield np.column_stack((ranks[rows[valid]], targets[valid]))

    def to_single(self, *args, **kwargs):
        if self._frontier is None:
            def ranks():
                yield from self._base_space.iter_index_blocks()
                for pairs in self.iter_pair_blocks():
                    yield pairs.reshape(-1)
            self._frontier = self.store.frontier(self.basis, ranks())
        return self._frontier

    def take_unique(self, *args, **kwargs):
        return self

    def remove_duplicates(self, *args, **kwargs):
        return self

    def concatenate(self, other, *args, **kwargs):
        return self.union(other)

    def union(self, other, *args, **kwargs):
        if not isinstance(other, DiskSelectionRuleStateSpace):
            if not isinstance(other, SelectionRuleStateSpace):
                raise TypeError('coupling graphs can only merge with coupling graphs')
            other = type(self).from_resident(other, self.store)
        if self.basis != other.basis:
            raise ValueError('different state-space bases')
        source = self._base_space.union(other._base_space)
        return type(self)(source, recipes=self.recipes + other.recipes)

    def intersection(self, other, handle_subspaces=True, **opts):
        if handle_subspaces:
            raise NotImplementedError('graph intersection requires source-only selection')
        source = self.store.coerce(other.representative_space)
        recipes = [(s.intersection(source), r, f, d, t, e) for s,r,f,d,t,e in self.recipes]
        return type(self)(self._base_space.intersection(source), recipes=recipes)

    def map(self, function):
        owner = getattr(function,'__self__',None)
        if type(owner).__name__ == 'ExclusionTestFilter' and not getattr(owner.test,'supports_state_blocks',False):
            raise ValueError('custom postfilter tests on disk must declare supports_state_blocks=True')
        recipes = [(s,r,f,d,t + (function,),e) for s,r,f,d,t,e in self.recipes]
        return type(self)(self._base_space, recipes=recipes)

    def filter_transitions(self, excluded_transitions, in_place=False):
        excluded = tuple(tuple(x) for x in excluded_transitions)
        recipes = [(s,r,f,d,t,e + excluded) for s,r,f,d,t,e in self.recipes]
        if in_place:
            self.recipes,self._frontier = recipes,None
            return self
        return type(self)(self._base_space, recipes=recipes)

    def drop_states(self, space):
        return self.map(lambda block: block.difference(space))

    def take_states(self, space):
        return self.map(lambda block: block.intersection(space))

    def canonical_pairs(self, other=None):
        allowed = None if other is None else self.store.coerce(other)
        def blocks():
            for pairs in self.iter_pair_blocks():
                if allowed is not None:
                    pairs = pairs[allowed._contains(pairs[:, 1])]
                yield np.sort(pairs, axis=1)
        return self.store.unique(blocks(), width=2)

    def get_representation_brakets(self, *args, **kwargs):
        raise RuntimeError('consume a deferred coupled space through iter_pair_blocks()')


def _mapped_csr(store, data, indices, indptr, n, extra_shape=()):
    import scipy.sparse as sp
    components = int(np.prod(extra_shape)) if extra_shape else 1
    matrix = sp.csr_matrix((data, indices, indptr), shape=(n,n*components), copy=False)
    matrix._state_space_store = store
    matrix._mapped_components = (data,indices,indptr)
    result = SparseArray.from_data(matrix, shape=(n,n)+extra_shape,cache_block_data=False)
    result._state_space_store = store
    return result


def disk_representation(representation, coupled, total, *, diagonal=False, filter_space=None,
                        zero_element_warning=True):
    """Evaluate unique pairs in blocks directly into mapped CSR components."""
    store = total.store
    size, n = store.block_size, len(total)
    if diagonal:
        dtype = np.int32 if n < 2**31-1 else np.int64
        data = store.array('diagonal', float, (n,))
        indices = store.array('indices', dtype, (n,))
        indptr = store.array('indptr', dtype, (n+1,))
        for start in range(0, n, size):
            stop = min(start+size,n)
            block = BasisStateSpace(total.basis, total.indices[start:stop], mode='indices')
            values = representation.get_brakets(BraKetSpace(block, block), check_orthogonality=False,
                                               memory_constrained=True)
            if isinstance(values, SparseArray):
                values = values.asarray()
            data[start:stop] = values
            indices[start:stop] = np.arange(start,stop)
            indptr[start:stop] = np.arange(start,stop)
        indptr[n] = n
        result = _mapped_csr(store,data,indices,indptr,n)
        result._mapped_diagonal = data
        return result
    if isinstance(coupled, DiskSelectionRuleStateSpace):
        pairs = coupled.canonical_pairs(other=filter_space)
    elif isinstance(coupled,SelectionRuleStateSpace):
        graph = DiskSelectionRuleStateSpace.from_resident(coupled, store)
        pairs = graph.canonical_pairs(other=filter_space)
    elif isinstance(coupled, BraKetSpace):
        def blocks():
            for start in range(0,len(coupled),size):
                yield np.sort(np.column_stack((coupled.bras.indices[start:start+size],
                                               coupled.kets.indices[start:start+size])),axis=1)
        pairs = store.unique(blocks(),width=2)
    else:
        coupled = store.coerce(coupled)
        graph = coupled.apply_selection_rules(representation.selection_rules)
        pairs = graph.canonical_pairs(other=coupled if filter_space is None else filter_space)
    def evaluate(block, warn=True):
        bras = BasisStateSpace(total.basis,block[:,0],mode='indices')
        kets = BasisStateSpace(total.basis,block[:,1],mode='indices')
        values = representation.get_brakets(BraKetSpace(bras,kets),memory_constrained=True)
        if isinstance(values,SparseArray):
            values = values.asarray()
        values = np.asarray(values)
        if values.ndim == 0:
            values = np.full(len(block),values)
        elif len(block)==1 and values.shape[0]!=1:
            values = values[np.newaxis]
        if values.shape[0] != len(block):
            raise ValueError('unexpected streamed matrix-element shape')
        if warn and zero_element_warning and np.any(values==0):
            raise ValueError('zero matrix elements; set zero_element_warning=False if expected')
        representation.clear_cache()
        SparseArray.clear_cache()
        return values
    sample = pairs[:1] if len(pairs) else np.array([[0,0]],dtype=np.int64)
    example = evaluate(sample,warn=False)
    extra_shape = example.shape[1:]
    components = int(np.prod(extra_shape)) if extra_shape else 1
    counts = store.array('counts', np.int64, (n+1,))
    counts[:] = 0
    for start in range(0,len(pairs),size):
        block = pairs[start:start+size]
        rows, cols = total.find(block[:,0]), total.find(block[:,1])
        np.add.at(counts,rows+1,components)
        np.add.at(counts,cols[rows!=cols]+1,components)
    carry = 0
    for start in range(0,n+1,size):
        block = np.cumsum(counts[start:start+size],dtype=np.int64) + carry
        counts[start:start+len(block)] = block
        if len(block):
            carry = int(block[-1])
    dtype = np.int32 if max(n*components,carry) < 2**31-1 else np.int64
    indptr = store.array('indptr',dtype,(n+1,))
    cursor = store.array('cursor',np.int64,(n,))
    for start in range(0,n+1,size):
        indptr[start:start+size] = counts[start:start+size]
    for start in range(0,n,size):
        cursor[start:start+size] = counts[start:min(start+size,n)]
    data = store.array('data',np.result_type(example.dtype,float),(carry,))
    indices = store.array('indices',dtype,(carry,))
    def insert(rows,cols,values):
        if components > 1:
            values = values.reshape(len(rows),components)
            cols = (cols[:,None]*components+np.arange(components)).reshape(-1)
            rows = np.repeat(rows,components)
            values = values.reshape(-1)
        order = np.argsort(rows,kind='stable')
        rows,cols,values = rows[order],cols[order],values[order]
        unique,starts,ct = np.unique(rows,return_index=True,return_counts=True)
        offsets = np.arange(len(rows)) - np.repeat(starts,ct)
        positions = cursor[rows] + offsets
        indices[positions], data[positions] = cols,values
        cursor[unique] += ct
    for start in range(0,len(pairs),size):
        block = pairs[start:start+size]
        values = evaluate(block)
        rows,cols = total.find(block[:,0]),total.find(block[:,1])
        insert(rows,cols,values)
        off = rows!=cols
        insert(cols[off],rows[off],values[off])
        representation.clear_cache()
    return _mapped_csr(store,data,indices,indptr,n,extra_shape=extra_shape)
