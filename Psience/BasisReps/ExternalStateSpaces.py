"""Mapped state frontiers and replayable, bounded coupling transformations."""
from pathlib import Path
import itertools, tempfile, weakref
import numpy as np
from McUtils.Numputils import external_unique, SparseArray
from McUtils.Combinatorics import SymmetricGroupGenerator, UniquePermutations
from .StateIndexers import PermutationStateIndexer
from .StateSpaces import BasisStateSpace, BasisStateSpaceView, SelectionRuleStateSpace, BraKetSpace

__all__ = ['StateSpaceDiskStore', 'DiskBasisStateSpace', 'DiskSelectionRuleStateSpace',
           'QuantaStateSpace']


def _cleanup_store(mappings, temporary):
    # NFS must see all mappings closed before unlinking their files. NumPy and
    # SciPy may retain ndarray views after the original memmap object is gone.
    for mapping in list(mappings):
        if not mapping.closed:
            mapping.close()
    temporary.cleanup()


class StateSpaceDiskStore:
    def __init__(self, directory=None, block_size=65536, *, memory_limit=None):
        if block_size < 1:
            raise ValueError('state_space_block_size must be positive')
        if memory_limit is not None and (not np.isfinite(memory_limit) or memory_limit<0):
            raise ValueError('state-space memory_limit must be nonnegative bytes or None')
        self.memory_limit=memory_limit
        if directory is not None:
            Path(directory).mkdir(parents=True, exist_ok=True)
        self._temporary = tempfile.TemporaryDirectory(prefix='coupled_states_', dir=directory)
        self.directory = Path(self._temporary.name)
        self.block_size = int(block_size)
        self._counter = itertools.count()
        self._mappings = weakref.WeakSet()
        self._finalize = weakref.finalize(self,_cleanup_store,self._mappings,self._temporary)
        self.metrics = dict(peak_source_states=0, peak_candidate_pairs=0, generated_pairs=0)

    def path(self, label):
        return self.directory / ('{}_{}.npy'.format(label, next(self._counter)))

    def own(self, array):
        if isinstance(array, np.memmap):
            # Retain cleanup reachability without owning dead frontier maps.
            if array._mmap is not None:self._mappings.add(array._mmap)
            array._state_space_store = self
        return array

    def array(self, label, dtype, shape):
        # One bounded work chunk does not need a live file/mapping descriptor.
        if np.prod(shape,dtype=object) <= self.block_size:
            return np.zeros(shape,dtype=dtype)
        return self.own(np.lib.format.open_memmap(self.path(label), mode='w+', dtype=dtype, shape=shape))

    def unique(self, arrays, width=1, sort=True):
        def chunks():
            yield np.empty((0,) if width == 1 else (0, width), dtype=np.int64)
            yield from arrays
        iterator=iter(chunks());head=[];count=0
        for block in iterator:
            head.append(block);count+=len(block)
            if count>self.block_size:break
        output=self.path('unique') if count>self.block_size else None
        result=external_unique(itertools.chain(head,iterator),sort=sort,directory=self.directory,
                               chunk_size=self.block_size,output_file=output)
        if isinstance(result,np.memmap) and len(result)<=self.block_size:
            # This result has not escaped: copy it before closing its sole map.
            resident=np.array(result,copy=True)
            result._mmap.close()
            Path(result.filename).unlink()
            result=resident
        return self.own(result)

    def frontier(self, basis, arrays, sort=True):
        return DiskBasisStateSpace(basis, self.unique(arrays,sort=sort), self,
                                  canonical=sort,unique=True)

    def coerce(self, space):
        if isinstance(space, DiskSelectionRuleStateSpace):
            return space.to_single()
        if isinstance(space, DiskBasisStateSpace) and space.store is self:
            return space
        return self.frontier(space.basis, index_blocks(space, self.block_size))

    def _union_views(self, spaces, sort):
        """Union common-parent CSR selections without copying every rank row."""
        parent=spaces[0].index_parent
        if getattr(parent,'store',None) is not self or not parent.is_unique():return None
        if len(spaces)==1:
            first=self.coerce(spaces[0])
            return first.take_unique(sort=True) if sort else first
        # One byte per already-existing parent position, never an adjacency
        # matrix or excitation array. Use it only within the configured state
        # memory budget and when it is no larger than the input index payload.
        size=len(parent)
        payload=sum(space.index_positions.nbytes for space in spaces)
        if size>payload or (self.memory_limit is not None and size>self.memory_limit):return None
        seen=np.zeros(size,dtype=bool)
        count=sum(len(space) for space in spaces)
        ranks=self.array('union_view_ranks',np.int64,(min(size,count),))
        # coerce() sorts the first resident input even for sort=False. Keep
        # that legacy ordering instead of emitting the first CSR row as-is.
        first=self.coerce(spaces[0]);offset=0
        for block in first.iter_index_blocks():
            ranks[offset:offset+len(block)]=block;offset+=len(block)
        for positions in spaces[0].iter_position_blocks(self.block_size):seen[positions]=True
        for space in spaces[1:]:
            for positions in space.iter_position_blocks(self.block_size):
                if space._view_unique is not True:
                    _,first_positions=np.unique(positions,return_index=True)
                    positions=positions[np.sort(first_positions)]
                new=positions[~seen[positions]]
                seen[new]=True
                values=parent.indices[new]
                if values.dtype.kind=='u' and len(values) and values.max()>np.iinfo(np.int64).max:
                    raise OverflowError('disk state ranks must fit signed int64')
                ranks[offset:offset+len(new)]=values
                offset+=len(new)
        result=DiskBasisStateSpace(parent.basis,ranks[:offset],self,canonical=False,unique=True)
        return result.take_unique(sort=True) if sort else result


    def union_spaces(self, spaces, sort=False):
        """Union a stream of frontiers in one bounded stable/sorted pass.

        Preserve the same first-space coercion and first-occurrence order as
        repeated ``coerce(first).union(next)``. Empty input returns None.
        """
        spaces = iter(spaces)
        sentinel = object()
        first = next(spaces, sentinel)
        if first is sentinel:
            return None
        if isinstance(first,BasisStateSpaceView) and first.index_parent is not None:
            views=[first];parent=first.index_parent
            for space in spaces:
                if not isinstance(space,BasisStateSpaceView) or space.index_parent is not parent:
                    spaces=itertools.chain(views,(space,),spaces)
                    break
                views.append(space)
            else:
                result=self._union_views(views,sort)
                if result is not None:return result
                spaces=iter(views)
            first=next(spaces)
        first = self.coerce(first)
        second = next(spaces, sentinel)
        if second is sentinel:
            return first.take_unique(sort=True) if sort else first
        basis = first.basis
        def blocks():
            yield from first.iter_index_blocks()
            for space in itertools.chain((second,), spaces):
                if space.basis != basis:
                    raise ValueError('different state-space bases')
                yield from index_blocks(space, self.block_size)
        return self.frontier(basis, blocks(), sort=sort)

    def close(self):
        self._finalize()


def index_blocks(space, size):
    if isinstance(space, QuantaStateSpace):
        yield from space.iter_index_blocks(size)
    elif hasattr(space, 'iter_index_blocks'):
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

    def _contains(self, values):
        return self.contains_indices(values)

    @property
    def quanta_bounds(self):
        return self.quanta, self.quanta

    def iter_index_blocks(self, block_size=65536):
        if type(self.basis.indexer) is not PermutationStateIndexer:
            raise ValueError('implicit quanta index blocks require the native graded indexer')
        for start in range(self.start, self.stop, block_size):
            yield np.arange(start, min(start+block_size, self.stop), dtype=np.int64)

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

    def iter_excitations(self, block_size=None, *, sparse=False):
        size = self.store.block_size if block_size is None else int(block_size)
        if size < 1:
            raise ValueError("block_size must be positive")
        for values in self.iter_index_blocks():
            for start in range(0,len(values),size):
                block = self.basis.unravel_state_inds(values[start:start+size])
                yield SparseArray.from_data(block) if sparse else block

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

    @property
    def quanta_bounds(self):
        """Bounds of actual target occupations for the native graded indexer."""
        if not getattr(self, '_quanta_bounds_ready', False):
            self._quanta_bounds = None
            if len(self) and hasattr(self.basis.indexer, 'symmetric_group'):
                if self._canonical:
                    endpoints = self._indices[[0, -1]]
                else:
                    # Preserve noncanonical ordering; do not build a lookup.
                    low, high = None, None
                    for block in self.iter_index_blocks():
                        lo, hi = int(np.min(block)), int(np.max(block))
                        low = lo if low is None else min(low, lo)
                        high = hi if high is None else max(high, hi)
                    endpoints = np.array([low, high], dtype=np.int64)
                occupations = self.basis.unravel_state_inds(endpoints)
                self._quanta_bounds = tuple(int(q) for q in np.sum(occupations, axis=1))
            self._quanta_bounds_ready = True
        return self._quanta_bounds



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
        self._frontier_components = None
        self._frontier_block_size = None

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

    def iter_pair_blocks(self, filter_space=None):
        """Generate directed pairs, pushing target restrictions before expansion.

        For native graded bases, generate from the smaller unique endpoint set
        when its inverse rules are cheaper by cardinality. Negating the rules
        and restoring (bra, ket) order preserves nonsymmetric rule directions.
        Explicit owned polyad records keep their original ownership handling.
        """
        size = self.store.block_size
        outer_filter = (filter_space if isinstance(filter_space, QuantaStateSpace) else
                        None if filter_space is None else self.store.coerce(filter_space))
        for source, rules, filtered, dimensions, transforms, excluded in self.recipes:
            source = self.store.coerce(source)
            filtered = (filtered if isinstance(filtered, QuantaStateSpace) else
                        None if filtered is None else self.store.coerce(filtered))
            filters = []
            # A generic postfilter can depend on the contents of its input
            # block. Preserve its original blocks and apply the outer filter
            # afterwards; only plain selection-rule recipes can be reversed.
            for space in (filtered, outer_filter if not transforms else None):
                if space is not None and all(space is not old for old in filters):
                    filters.append(space)
            if any(len(space) == 0 for space in filters):
                continue
            if isinstance(rules, _StoredCouplings):
                for start in range(0, len(rules.records), size):
                    block = rules.records[start:start+size]
                    block = block[source._contains(block[:, 0])]
                    pairs = block[:, 1:]
                    valid = np.ones(len(pairs), dtype=bool)
                    for space in filters:
                        valid &= space._contains(pairs[:, 1])
                    if excluded:
                        delta = (self.basis.unravel_state_inds(pairs[:, 1])
                                 - self.basis.unravel_state_inds(pairs[:, 0]))
                        for rule in excluded:
                            valid &= ~np.all(delta == rule, axis=1)
                    for transform in transforms:
                        kept = transform(BasisStateSpace(self.basis, pairs[valid, 1], mode='indices')).indices
                        valid &= np.isin(pairs[:, 1], kept)
                    if transforms and outer_filter is not None:
                        valid &= outer_filter._contains(pairs[:, 1])
                    if np.any(valid):
                        yield pairs[valid]
                continue
            if rules is None:
                # None means an unrestricted graph on the original targets,
                # not an identity operator or a known polynomial selection rule.
                targets = source if filtered is None else filtered
                for rank_block in source.iter_index_blocks():
                    for rank in rank_block:
                        for target_block in index_blocks(targets, size):
                            selected = target_block
                            if excluded:
                                delta = self.basis.unravel_state_inds(selected)-self.basis.unravel_state_inds(np.array([rank]))
                                valid = np.ones(len(selected), dtype=bool)
                                for rule in excluded:
                                    valid &= ~np.all(delta == rule, axis=1)
                                selected = selected[valid]
                            for transform in transforms:
                                selected = transform(BasisStateSpace(self.basis, selected, mode='indices')).indices
                            if outer_filter is not None:
                                selected = selected[outer_filter._contains(selected)]
                            if len(selected):
                                yield np.column_stack((np.full(len(selected), rank, dtype=np.int64), selected))
                continue

            targets = min(filters, key=len) if filters else None
            # Sorted, unique sources have an existing bounded membership lookup.
            # Reversal must not create a large extra noncanonical lookup table.
            reverse = (not transforms and targets is not None and len(targets) < len(source)
                       and not isinstance(targets, QuantaStateSpace)
                       and source.is_sorted() and source.is_unique() and targets.is_unique()
                       and source.quanta_bounds is not None)
            seeds, membership = (targets, source) if reverse else (source, targets)
            active_rules = tuple(tuple(-int(c) if reverse else int(c) for c in rule) for rule in rules)
            if transforms:
                bounds = None
            elif reverse:
                bounds = membership.quanta_bounds
            else:
                known = [space.quanta_bounds for space in filters if space.quanta_bounds is not None]
                bounds = ((max(q[0] for q in known), min(q[1] for q in known)) if known else None)
            if bounds is not None and bounds[0] > bounds[1]:
                continue
            rule_sums = np.array([sum(rule) for rule in active_rules], dtype=np.int64)
            for rank_block in index_blocks(seeds, size):
                for start in range(0, len(rank_block), min(64, size)):
                    ranks = rank_block[start:start+min(64, size)]
                    if reverse:
                        for space in filters:
                            if space is not seeds:
                                ranks = ranks[space._contains(ranks)]
                    if not len(ranks):
                        continue
                    original = self.basis.unravel_state_inds(ranks).astype(np.int64, copy=False)
                    self.store.metrics['peak_source_states'] = max(self.store.metrics['peak_source_states'], len(ranks))
                    if bounds is None:
                        batches = [(np.arange(len(ranks)), active_rules)]
                    else:
                        quanta = np.sum(original, axis=1)
                        batches = []
                        for q in np.unique(quanta):
                            allowed = (q+rule_sums >= bounds[0]) & (q+rule_sums <= bounds[1])
                            if np.any(allowed):
                                batches.append((np.flatnonzero(quanta == q),
                                                tuple(rule for rule, keep in zip(active_rules, allowed) if keep)))
                    for selected, selected_rules in batches:
                        exc = original[selected] if dimensions is None else original[selected][:, dimensions]
                        index_group = (self.basis.indexer.symmetric_group
                            if type(self.basis.indexer) is PermutationStateIndexer else None)
                        index_only = (dimensions is None and index_group is not None
                                      and len(excluded)*len(original) <= size)
                        group = index_group if index_only else SymmetricGroupGenerator(exc.shape[1])
                        forbidden = []
                        if index_only:
                            # Same directed exclusions, encoded once per source block.
                            # Cache size is bounded by the existing output block budget.
                            for rule in excluded:
                                vectors = original-rule if reverse else original+rule
                                possible = np.all(vectors >= 0,axis=1)
                                values = np.full(len(original),-1,dtype=np.int64)
                                if np.any(possible):
                                    ids = self.basis.ravel_state_inds(vectors[possible])
                                    within = ids <= np.iinfo(np.int64).max
                                    positions = np.flatnonzero(possible)
                                    values[positions[within]] = ids[within]
                                forbidden.append(values)
                        for rows, outputs in group.iter_permutation_rule_direct_sum(
                                exc,selected_rules,block_size=size,return_indices=index_only):
                            rows = selected[rows]
                            self.store.metrics['peak_candidate_pairs'] = max(self.store.metrics['peak_candidate_pairs'], len(outputs))
                            self.store.metrics['generated_pairs'] += len(outputs)
                            if dimensions is not None:
                                restored = original[rows].copy()
                                restored[:, dimensions] = outputs
                                outputs = restored
                            valid = np.ones(len(outputs), dtype=bool)
                            if excluded and not index_only:
                                delta = original[rows]-outputs if reverse else outputs-original[rows]
                                for rule in excluded:
                                    valid &= ~np.all(delta == rule, axis=1)
                            generated = outputs if index_only else self.basis.ravel_state_inds(outputs)
                            if np.any(generated > np.iinfo(np.int64).max):
                                raise OverflowError('disk state ranks must fit signed int64')
                            generated = generated.astype(np.int64, copy=False)
                            for values in forbidden:
                                valid &= generated != values[rows]
                            bras, kets = (generated, ranks[rows]) if reverse else (ranks[rows], generated)
                            if membership is not None:
                                valid &= membership._contains(generated)
                            if not reverse:
                                for space in filters:
                                    if space is not membership:
                                        valid &= space._contains(kets)
                            for transform in transforms:
                                kept = transform(BasisStateSpace(self.basis, kets[valid], mode='indices')).indices
                                valid &= np.isin(kets, kept)
                            if transforms and outer_filter is not None:
                                valid &= outer_filter._contains(kets)
                            if np.any(valid):
                                yield np.column_stack((bras[valid], kets[valid]))


    def _get_structured_frontier(self):
        """Commute only declared pointwise predicates with node-set union."""
        grouped = {}
        for source,rules,filtered,dimensions,transforms,excluded in self.recipes:
            descriptors = []
            for function in transforms:
                owner = getattr(function,'__self__',None)
                describe = getattr(owner,'get_permutation_constraints',None)
                descriptor = None if describe is None else describe()
                if descriptor is None:
                    return None
                descriptors.append(descriptor)
            constraints = UniquePermutations.prepare_constraints(descriptors,self.basis.ndim)
            key = constraints['signature']
            if key not in grouped:
                grouped[key] = (constraints,[])
            grouped[key][1].append((source,rules,filtered,dimensions,(),excluded))
        store = self.store
        def blocks():
            yield from self._base_space.iter_index_blocks()
            for constraints,recipes in grouped.values():
                raw = type(self)(self._base_space,recipes=recipes)
                # Class predicates reject impossible target shapes before
                # table emission. Mode-labelled predicates remain a safe
                # outer bound here and are checked by the prefix traversal.
                frontier = raw._get_class_frontier(class_constraints=constraints)
                if frontier is None:
                    frontier = raw.to_single()
                group = self.basis.indexer.symmetric_group
                statistics = {}
                yield from group.iter_filtered_indices(frontier.indices,constraints,
                                                       store.block_size,statistics=statistics)
                for name,value in statistics.items():
                    key = 'constraint_'+name
                    store.metrics[key] = store.metrics.get(key,0)+value
        return store.frontier(self.basis,blocks())


    def _get_class_frontier(self, class_constraints=None):
        """Evaluate set-valued rule propagation over complete source orbits.

        Edge recipes are preserved. Explicit target membership filters can be
        intersected with class intervals; arbitrary block callbacks and
        mode-specific edge restrictions retain the original pair evaluator.
        """
        if type(self.basis.indexer) is not PermutationStateIndexer:
            return None
        if any(t for s,r,f,d,t,e in self.recipes):
            return self._get_structured_frontier()
        if any(r is None or isinstance(r, _StoredCouplings) or d is not None or t or e
               for s, r, f, d, t, e in self.recipes):
            return None
        if any(f is not None and not isinstance(f, (BasisStateSpace, QuantaStateSpace))
               for s, r, f, d, t, e in self.recipes):
            return None
        store, group = self.store, self.basis.indexer.symmetric_group
        groups, sources, filters = {}, {}, {}
        # A source may be split across recipes. Aggregate support separately
        # for each rule and the SAME target restriction, never a cross product
        # between independent unions of sources and rules.
        for source, rules, filtered, dimensions, transforms, excluded in self.recipes:
            identity = id(source)
            if identity not in sources:
                canonical = store.coerce(source).take_unique(sort=True)
                if any(not np.all(self._base_space._contains(block))
                       for block in canonical.iter_index_blocks()):
                    return None  # accepted bras outside the base need edge ownership
                sources[identity] = canonical
            key = (('quanta', filtered.quanta) if isinstance(filtered, QuantaStateSpace)
                   else ('space', id(filtered)))
            filters[key] = filtered
            for rule in {tuple(sorted(int(c) for c in r if c)) for r in rules}:
                if len(rule) <= self.basis.ndim:
                    groups.setdefault((rule, key), {})[identity] = sources[identity]

        def merge(intervals):
            result = []
            for low, high in sorted(intervals):
                low, high = int(low), int(high)
                if low >= high:
                    continue
                if result and low <= result[-1][1]:
                    result[-1][1] = max(result[-1][1], high)
                else:
                    result.append([low, high])
            return result

        def covered(intervals, coverage):
            j = 0
            for low, high in intervals:
                while j < len(coverage) and coverage[j][1] <= low:
                    j += 1
                if j == len(coverage) or coverage[j][0] > low or coverage[j][1] < high:
                    return False
            return True

        def class_support(space):
            return [(tuple(map(int,partition)),low,high,support)
                    for partition,low,high,support in group.iter_index_classes(space.indices)]

        supports, grouped = {}, []
        for (rule, key), contribution in groups.items():
            identities = tuple(sorted(contribution))
            if identities not in supports:
                space = store.union_spaces(contribution.values(), sort=True)
                supports[identities] = class_support(space)
            grouped.append((rule, filters[key], supports[identities]))
        base = store.coerce(self._base_space).take_unique(sort=True)
        complete = [[low, high] for p, low, high, support in class_support(base)
                    if len(support) == high-low]
        pending, filtered_complete, target_cache = [], [], {}
        full_classes = partial_classes = 0

        def target_ranges(partition, rule, filtered):
            cache_key = partition, rule
            if cache_key not in target_cache:
                target_cache[cache_key] = group.get_direct_sum_class_ranges(partition, [rule], constraints=class_constraints).tolist()
            ranges = target_cache[cache_key]
            if isinstance(filtered, QuantaStateSpace):
                return [[max(low, filtered.start), min(high, filtered.stop)] for low, high in ranges
                        if max(low, filtered.start) < min(high, filtered.stop)]
            return ranges

        for rule, filtered, classes in grouped:
            for partition, low, high, support in classes:
                ranges = target_ranges(partition, rule, filtered)
                if len(support) == high-low:
                    full_classes += 1
                    if filtered is None or isinstance(filtered, QuantaStateSpace):
                        complete.extend(ranges)
                    else:
                        filtered_complete.append((ranges, filtered))
                else:
                    partial_classes += 1
                    pending.append((ranges, support, rule, filtered))
        if not full_classes:
            return None
        complete = merge(complete)
        pending = [item for item in pending if not covered(item[0], complete)]

        def in_ranges(values, ranges):
            if not ranges:
                return np.zeros(len(values), dtype=bool)
            intervals = np.asarray(ranges, dtype=np.int64)
            positions = np.searchsorted(intervals[:, 0], values, side='right')-1
            return (positions >= 0) & (values < intervals[np.maximum(positions, 0), 1])

        def ranges_contain(values):
            return in_ranges(values, complete)

        base_covered = all(np.all(ranges_contain(block)) for block in base.iter_index_blocks())
        store.metrics['class_frontier_queries'] = store.metrics.get('class_frontier_queries', 0)+1
        store.metrics['class_complete_sources'] = store.metrics.get('class_complete_sources', 0)+full_classes
        store.metrics['class_partial_sources'] = store.metrics.get('class_partial_sources', 0)+partial_classes
        store.metrics['class_pending_sources'] = store.metrics.get('class_pending_sources', 0)+len(pending)
        store.metrics['class_output_intervals'] = store.metrics.get('class_output_intervals', 0)+len(complete)
        if not pending and not filtered_complete and base_covered:
            count = sum(high-low for low, high in complete)
            values = store.array('class_frontier', np.int64, (count,))
            position = 0
            for low, high in complete:
                for start in range(low, high, store.block_size):
                    stop = min(start+store.block_size, high)
                    values[position:position+stop-start] = np.arange(start, stop, dtype=np.int64)
                    position += stop-start
            result = DiskBasisStateSpace(self.basis, values, store, canonical=True, unique=True)
        else:
            def blocks():
                yield from base.iter_index_blocks()
                for low, high in complete:
                    for start in range(low, high, store.block_size):
                        yield np.arange(start, min(start+store.block_size, high), dtype=np.int64)
                for ranges, filtered in filtered_complete:
                    for values in index_blocks(filtered, store.block_size):
                        yield values[in_ranges(values, ranges) & ~ranges_contain(values)]
                for ranges, support, rule, filtered in pending:
                    selected = DiskBasisStateSpace(self.basis, support, store, canonical=True, unique=True)
                    # Native quanta filters are handled by the class bounds;
                    # use the same membership rule on each partial output.
                    graph = type(self)(selected, [rule], filter_space=None if isinstance(filtered, QuantaStateSpace) else filtered)
                    for pairs in graph.iter_pair_blocks():
                        values = pairs[:, 1]
                        mask = in_ranges(values, ranges) & ~ranges_contain(values)
                        yield values[mask]
            result = store.frontier(self.basis, blocks())
        return result


    def _get_frontier_components(self, store):
        # Block-dependent postfilters must see the same store/coercion and
        # block boundaries as the computation that produced the cached nodes.
        if self.store is store and self._frontier_block_size == store.block_size:
            if self._frontier is not None:
                return ((self._frontier, None),)
            return self._frontier_components
        return None

    def to_single(self, *args, **kwargs):
        if self._frontier is None:
            self._frontier = self._get_class_frontier()
            if self._frontier is not None:
                self._frontier_block_size = self.store.block_size
                self._frontier_components = None
                return self._frontier
            def ranks():
                yield from self._base_space.iter_index_blocks()
                components = self._get_frontier_components(self.store)
                if components is None:
                    for pairs in self.iter_pair_blocks():
                        yield pairs.reshape(-1)
                else:
                    for frontier, recipes in components:
                        if frontier is not None:
                            yield from frontier.iter_index_blocks()
                        else:
                            # Evaluate a frozen recipe contribution under this
                            # union's store and basis, preserving callback blocks.
                            fragment = type(self)(self._base_space, recipes=recipes)
                            for pairs in fragment.iter_pair_blocks():
                                yield pairs.reshape(-1)
            self._frontier = self.store.frontier(self.basis, ranks())
            self._frontier_block_size = self.store.block_size
            # Once merged, release intermediate node-map ownership. Coupling
            # recipes remain intact for subsequent Hamiltonian evaluation.
            self._frontier_components = None
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
        result = type(self)(source, recipes=self.recipes + other.recipes)
        left = self._get_frontier_components(self.store)
        right = other._get_frontier_components(self.store)
        if left is not None or right is not None:
            # Snapshots protect a pending contribution from later in-place
            # recipe edits to its original graph. No edges are evaluated here.
            if left is None:
                left = ((None, tuple(self.recipes)),)
            if right is None:
                right = ((None, tuple(other.recipes)),)
            result._frontier_components = left + right
            result._frontier_block_size = self.store.block_size
        return result

    def intersection(self, other, handle_subspaces=True, **opts):
        if handle_subspaces:
            raise NotImplementedError('graph intersection requires source-only selection')
        source = self.store.coerce(other.representative_space)
        recipes = [(s.intersection(source), r, f, d, t, e) for s,r,f,d,t,e in self.recipes]
        return type(self)(self._base_space.intersection(source), recipes=recipes)

    def map(self, function):
        owner = getattr(function,'__self__',None)
        if (type(owner).__name__ == 'ExclusionTestFilter'
                and owner.get_permutation_constraints() is None
                and not getattr(owner.test,'supports_state_blocks',False)):
            raise ValueError('custom postfilter tests on disk must declare supports_state_blocks=True')
        recipes = [(s,r,f,d,t + (function,),e) for s,r,f,d,t,e in self.recipes]
        return type(self)(self._base_space, recipes=recipes)

    def filter_transitions(self, excluded_transitions, in_place=False):
        excluded = tuple(tuple(x) for x in excluded_transitions)
        recipes = [(s,r,f,d,t,e + excluded) for s,r,f,d,t,e in self.recipes]
        if in_place:
            self.recipes,self._frontier = recipes,None
            self._frontier_components = self._frontier_block_size = None
            return self
        return type(self)(self._base_space, recipes=recipes)

    def drop_states(self, space):
        return self.map(lambda block: block.difference(space))

    def take_states(self, space):
        return self.map(lambda block: block.intersection(space))

    def canonical_pairs(self, other=None):
        def blocks():
            for pairs in self.iter_pair_blocks(filter_space=other):
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
