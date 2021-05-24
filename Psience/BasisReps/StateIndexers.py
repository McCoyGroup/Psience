"""
Provides utilities for converting state vectors to indices.
Used implicitly by various conversion utilities
"""

import abc, numpy as np, collections

__all__ = [
    "BaseStateIndexer",
    "ArrayStateIndexer",
    "SpaceStateIndexer",
    "PermutationStateIndexer"
]

class BaseStateIndexer(metaclass=abc.ABCMeta):
    """
    Provides a base class for indexing states
    """
    @abc.abstractmethod
    def to_indices(self, states):
        """
        Converts the set of states to numerical indices

        :param states:
        :type states:
        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")
    @abc.abstractmethod
    def from_indices(self, indices):
        """
        Converts the set of states to numerical indices

        :param states:
        :type states:
        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")

class ArrayStateIndexer(BaseStateIndexer):
    """
    Very simple indexer that takes a set of max state dimensions and
    provides the appropriate indices in that space
    """
    def __init__(self, dims):
        self.dims = np.array(dims)

    def to_state(self, serializer=None):
        return {
            'dims': self.dims
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(data['dims'])

    def to_indices(self, states):
        idx = np.asarray(states, dtype=int)
        if len(states) == 0:
            return idx
        if idx.ndim == 1:
            idx = idx[np.newaxis]
        return np.ravel_multi_index(idx.T, self.dims)
    def from_indices(self, indices):
        if len(indices) == 0:
            return np.array([], dtype='int8')
        unravel = np.unravel_index(indices, self.dims)
        return np.array(unravel, dtype='int8').T

class SpaceStateIndexer(BaseStateIndexer):
    """
    Very simple indexer that takes a set of states and indexes based on that.
    Needs some testing to make sure everything works as cleanly as hoped.
    """
    def __init__(self, states):
        self.og_states = np.asarray(states, dtype=int)
        self.indexer = np.lexsort(states.T)
        if self.og_states.ndim != 2:
            raise ValueError("can't index anything other than vectors of excitations")
    def to_indices(self, states):
        if len(states) == 0:
            return np.array([], dtype=int)
        return np.searchsorted(self.og_states, states, indexer=self.indexer)
    def from_indices(self, indices):
        return self.og_states[indices,]

class PermutationStateIndexer(BaseStateIndexer):
    """
    A sophisticated indexer that takes a state dimension and provides
    indices based on the `shortlex` ordering, where ordering is defined
    first by # of quanta of excitation, then by which partitioning of the #quanta,
     it represents, and finally by which permutation of that paritioning it is.
    Is likely about as stable as an indexer can be expected to be over large
    numbers of states. Unlikely to exhaust the max integers available for most
    systems.
    """
    def __init__(self, ndim):
        self.ndim = ndim
        self.partitioning_cache = {}

    def to_state(self, serializer=None):
        return {
            'ndim': self.ndim
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(data['ndim'])

    def integer_partitions(self, num):
        """
        Gives basically the second sort criterion by calculating integer partitions
        (i.e. number of way to split up num quanta across modes)
        Sorted by default, which is our saving grace
        :param num:
        :type num:
        :return:
        :rtype: Tuple[int, int, Iterable[PermutationStateIndexer.PartitionPermutationIndexer]]
        """

        if num not in self.partitioning_cache:
            partitions = np.array(
                [list(reversed(x))+[0]*(self.ndim - len(x)) for x in self._accel_asc(num) if len(x) <= self.ndim],
                dtype='int8'
            )
            sorting = np.lexsort(partitions.T)
            indexers = [ PartitionPermutationIndexer(x) for x in partitions[sorting,] ]
            total_dim = np.sum([x.total_states for x in indexers])
            if num > 0:
                prev_dim = self.integer_partitions(num-1)[1]
                total_dim += prev_dim # recursively we get it all
            else:
                prev_dim = 0
            self.partitioning_cache[num] = (prev_dim, total_dim, indexers)

        return self.partitioning_cache[num]

    def _get_inds_by_quanta(self, states, quanta, assume_sorted=False):
        """
        Pre-sorts states so that they're grouped by quanta and then searches
        for the appropriate int based on the integer partitions at that number of quanta

        :param states:
        :type states:
        :param quanta:
        :type quanta:
        :return:
        :rtype:
        """

        # partitions we use for indexing
        baseline, _, partitions = self.integer_partitions(quanta)
        cls_map = {} # cache of classes so that we can compare fewer things
        if len(states) == 1:
            # no reason to do all this sorting... (and it also breaks down...)
            inds = np.full(len(states), -1)
            s = states[0]
            idx = baseline
            cls, cnt = np.unique(s, return_counts=True)
            cls = np.flip(cls)
            cnt = np.flip(cnt)
            nz = sum(x for x,c in zip(cnt, cls) if c != 0)
            cls_key = PermutationStateKey(nz, tuple(cls))
            for indexer in partitions:
                ikey = indexer.key
                if ikey == cls_key: # it's a match!
                    idx += indexer.get_perm_indices(states)[0]
                    break
                else:
                    # add on number of skipped tree nodes
                    idx += indexer.total_states
            inds[0] = idx
        else:
            # build an array of state classes to sort and stuff
            # to minimize checks down the line
            make_class = lambda cls, cnt: PermutationStateKey(sum(x for x,c in zip(cnt, cls) if c!=0), tuple(np.flip(cls)))
            state_dat = np.array([make_class(*np.unique(s, return_counts=True)) for s in states], dtype=object)
            if not assume_sorted:
                # we do an initial sorting of the states and the data
                # so that we can be sure we have a good ordering for
                # later
                # the sorting is done so that things are first sorted by number of terms in the
                # integer partition, then next by the integer partition itself, and then
                # finally by the actual permutation
                # where the sort by integer partition is done so that a smaller initial value is prioritized (I think)
                # and the sort by actual state is done so a larger initial value is prioritized
                sort_states = np.sort(states, axis=-1)
                sort_classes = np.array([
                    np.concatenate([-np.flip(s), np.flip(z), [c[0]]])
                    for c, s, z in zip(state_dat, states, sort_states)
                ], dtype='int8')
                lexxy = np.lexsort(sort_classes.T)
                state_dat = state_dat[lexxy]
                states = states[lexxy]
            else:
                lexxy = None

            inds = np.full(len(states), -1)
            n = 0
            idxer = partitions[n] #type: PermutationStateIndexer.PartitionPermutationIndexer
            ikey = idxer.key
            idxer_block = []
            i_start = 0
            idx = baseline

            for key, state in zip(state_dat, states):
                # we create the block corresponding
                # to the current indexer, since we know our
                # states and data are sorted appropriately
                key = tuple(key)
                if key in cls_map:
                    idxer_block.append(state)
                    continue
                while key != ikey:
                    i_end = i_start + len(idxer_block)
                    if len(idxer_block) > 0:
                        # print(idxer, idx, idxer_block)
                        # try:
                        block_inds = idxer.get_perm_indices(idxer_block)
                        # except:
                        #
                        #     raise Exception(
                        #         idxer_block,
                        #         np.array(
                        #             [x for i,x in enumerate(sort_classes[lexxy]) if sum(states[i]) == 10]
                        #         )
                        #         )
                        # print(" >", block_inds)
                        inds[i_start:i_end] = idx + block_inds
                    n += 1
                    if n >= len(partitions):
                        break
                    # add on previous indexer states
                    # print("  ", partitions, key, idxer.total_states)
                    idx += idxer.total_states
                    idxer = partitions[n]
                    ikey = idxer.key
                    idxer_block = []
                    i_start = i_end
                else:
                    idxer_block.append(state)

                if n >= len(partitions):
                    raise ValueError("couldn't find indexer for state {} in {}; usually means sorting of states failed".format(state, partitions))

            i_end = i_start + len(idxer_block)
            if len(idxer_block) > 0:
                # print(idxer, idx, idxer_block)
                block_inds = idxer.get_perm_indices(idxer_block)
                # print(" >", block_inds)
                inds[i_start:i_end] = idx + block_inds

            if lexxy is not None:
                inds = inds[np.argsort(lexxy),]

        return inds

    def to_indices(self, states):
        """
        Finds the appropriate integer partitioning for each state
        :param states: 2D array of states as excitations
        :type states: np.ndarray
        :return:
        :rtype:
        """
        # group states by number of quanta
        states = np.asarray(states, dtype=int)
        if len(states) == 0:
            return states

        if states.ndim == 1:
            smol = True
            states = states[np.newaxis]
        else:
            smol = False
        # we'll actually now pull the unique states since that will potentially
        # allow us to be much faster
        states, inverse = np.unique(states, axis=0, return_inverse=True)
        nquants = np.sum(states, axis=1)
        inds = np.full_like(nquants, -1) # set up return value
        for q in np.unique(nquants):
            where_are_we = np.where(nquants == q)[0]
            these_states = states[where_are_we]
            inds[where_are_we] = self._get_inds_by_quanta(these_states, q)
        if smol:
            inds = inds[0]
        return inds[inverse,]

    def _from_quanta_inds(self, inds, indexers, check=True):
        """
        Gets the excitation vectors for inds
        based on a set of permutation indexers

        :param inds:
        :type inds: np.ndarray
        :param indexers:
        :type indexers: tuple[PermutationStateIndexer.PartitionPermutationIndexer]
        :return:
        :rtype:
        """

        og = inds
        res = np.full((len(inds), self.ndim), -1, dtype='int8')
        block_skips = 0
        block_start = 0
        for indexer in indexers:
            # figure out how many states this indexer supports
            block_size = indexer.total_states
            block_end = block_start + block_size
            block_inds = np.where(inds < block_end)[0]
            # we know things are sorted, so we add a constant skip
            # on before sticking back into the res array
            if len(block_inds) > 0:
                res_inds = block_inds + block_skips
                res[res_inds,] = indexer.from_perm_indices(inds[block_inds,] - block_start, assume_sorted=True)
                # finally slice the remaining inds so that we don't do excess where calls
                block_skips += len(block_inds)
                inds = inds[len(block_inds):]
                if len(inds) == 0:
                    break
            block_start = block_end

        if check:
            bad = np.where(np.all(res == -1, axis=1))[0]
            if len(bad) > 0:
                raise ValueError("failed to get permutation for indices {} from indexers {}".format(
                    np.asarray(og, dtype=int)[bad],
                    indexers
                ))

        return res

    max_quants = 1000  # just putting an arbitrary bound on the max quanta to prevent infinite loops
    def from_indices(self, indices, check=False):
        """
        Inverts the index calculation.
        First determines what number of quanta the index corresponds to,
        then which integer partition, and finally just loops through the unique
        permutations of the partition to get the right one.
        This is not assured to be a fast process in any way.

        :param indices:
        :type indices: Iterable[int]
        :return:
        :rtype:
        """

        if len(indices) == 0:
            return np.array([], dtype='int8')

        # we'll pre-sort this so that we can be
        # more efficient in how we apply the algorithm
        rem_inds, inverse = np.unique(np.asanyarray(indices, dtype=int), return_inverse=True)
        res = np.full((len(rem_inds), self.ndim), -1, dtype='int8')
        n = 0
        block_skips=0
        while len(rem_inds)>0 and n < self.max_quants:
            # use successive where calls to chop the
            # set of states down bit by bit
            # then calculate the indices for that number of quanta
            prev_num, block_num, indexers = self.integer_partitions(n)
            # we can be more efficient than this, but in the aggregate
            # it won't be significant I think
            block_inds = np.where(rem_inds < block_num)[0]
            # we know things are sorted, so we add a constant skip
            # on before sticking back into the res array
            res_inds = block_inds + block_skips
            res[res_inds,] = self._from_quanta_inds(rem_inds[block_inds,] - prev_num, indexers)
            # finally slice the remaining inds so that we don't do excess where calls
            block_skips += len(block_inds)
            rem_inds = rem_inds[len(block_inds):]
            n += 1

        res = res[inverse]
        if check:
            bad = np.where(np.all(res == -1, axis=1))[0]
            if len(bad) > 0:
                raise ValueError("failed to get permutation for indices {}".format(
                    np.asarray(indices, dtype=int)[bad]
                ))

        return res


