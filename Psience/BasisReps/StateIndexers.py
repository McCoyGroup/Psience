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
        return np.searchsorted(self.og_states, states, indexer=self.indexer)
    def from_indices(self, indices):
        return self.og_states[indices,]

PermutationStateKey = collections.namedtuple("PermutationStateKey", ['non_zero', 'classes'])
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

    @staticmethod
    def _accel_asc(n):
        """
        Pulled directly from http://jeromekelleher.net/author/jerome-kelleher.html
        Could easily be translated to C++ if speed is crucial (but this will never be a bottleneck)
        :param n:
        :type n:
        :return:
        :rtype:
        """

        a = [0 for i in range(n + 1)]
        k = 1
        y = n - 1
        while k != 0:
            x = a[k - 1] + 1
            k -= 1
            while 2 * x <= y:
                a[k] = x
                y -= x
                k += 1
            l = k + 1
            while x <= y:
                a[k] = x
                a[l] = y
                yield a[:k + 2]
                x += 1
                y -= 1
            a[k] = x + y
            y = x + y - 1
            yield a[:k + 1]

    class PartitionPermutationIndexer:
        """
        An order statistics tree designed to make it easy to
        get the index of a permutation based on the integer partition
        it comes from, which gives the number of character classes,
        and overall permutation length
        """
        def __init__(self, partition):
            self.ndim = len(partition)
            self.partition = partition
            self.classes, self.counts = np.unique(partition, return_counts=True)
            self.classes = np.flip(self.classes)
            self.counts = np.flip(self.counts)
            self.non_zero = sum(x for x,c in zip(self.counts, self.classes) if c != 0)
            self.key = PermutationStateKey(self.non_zero, tuple(self.classes))
            self.total_states = self._fac_rat(self.counts)

        @staticmethod
        def _fac_rat(counts):
            import math

            subfac = np.prod([math.factorial(x) for x in counts])
            ndim_fac = math.factorial(np.sum(counts))

            return ndim_fac // subfac

        @staticmethod
        def _subtree_counts(total, ndim, counts, where):
            """
            Computes the number of states in the tree built from decrementing counts[where] by 1
            Is it trivially simple? Yes
            But there's power to having it be named.
            :param total:
            :type total:
            :param ndim:
            :type ndim:
            :param counts:
            :type counts:
            :param where:
            :type where:
            :return:
            :rtype:
            """
            mprod = total * counts[where]
            if mprod % ndim != 0:
                raise ValueError("subtree counts {} don't comport with dimension {}".format(
                    mprod, ndim
                ))
            return mprod // ndim

        @staticmethod
        def _unique_permutations(elements):
            """
            From StackOverflow, an efficient enough
            method to get unique permutations
            :param perms:
            :type perms:
            :return:
            :rtype:
            """

            class unique_element:
                def __init__(self, value, occurrences):
                    self.value = value
                    self.occurrences = occurrences

            def perm_unique_helper(listunique, result_list, d):
                if d < 0:
                    yield tuple(result_list)
                else:
                    for i in listunique:
                        if i.occurrences > 0:
                            result_list[d] = i.value
                            i.occurrences -= 1
                            for g in perm_unique_helper(listunique, result_list, d - 1):
                                yield g
                            i.occurrences += 1

            if not hasattr(elements, 'count'):
                elements = list(elements)
            eset = set(elements)
            listunique = [unique_element(i, elements.count(i)) for i in eset]
            u = len(elements)
            return list(sorted(list(perm_unique_helper(listunique, [0] * u, u - 1)), reverse=True))

        def get_perm_indices(self, states, assume_sorted=True):
            """
            Gets the indices for a set of states.
            Does this by looping through the state, decrementing the appropriate character class,
            computing the number of nodes in the child tree (done based on the initial total states),
            and then adding up all those terms.
            We make use of the assumption that states are sorted to avoid doing more work than necessary
            by reusing stuff from the previous state

            :param state:
            :type state:
            :return:
            :rtype:
            """
            if not assume_sorted:
                raise NotImplementedError("need to sort")

            num_before = 0
            cur_total = self.total_states
            ndim = self.ndim
            cur_dim = self.ndim
            cur_classes = np.copy(self.classes)
            cur_counts = np.copy(self.counts)
            stack = collections.deque() # stack of current data to reuse
            # determine where each successive state differs so we can know how much to reuse
            diffs = np.diff(states, axis=0)
            inds = np.full((len(states),), -1)
            for sn,state in enumerate(states):
                if sn > 0:
                    # we reuse as much work as we can by only popping a few elements off of the class/counts stacks
                    num_diff = ndim - np.where(diffs[sn-1] != 0)[0][0] # number of differing states
                    if num_diff == 0: # same state so just reuse the previous value
                        if inds[sn - 1] == -1:
                            raise ValueError("state {} tried to reused bad value from state {}".format(
                                states[sn], states[sn-1]
                            ))
                        inds[sn] = inds[sn - 1]
                        continue
                    # we pop until we know the states agree once more
                    for n in range(num_diff):
                        num_before, cur_total, cur_classes, cur_counts = stack.pop()
                    cur_dim = num_diff
                    # print("  ::>", "{:>2}".format(sn), len(stack), state)
                    state = state[-num_diff:]
                    # print("    + ", state)
                    # print("    +", "{:>2}".format(num_before))
                # tree traversal, counting leaves in the subtrees
                for i, el in enumerate(state):
                    cur_num = num_before
                    for j, v in enumerate(cur_classes):
                        subtotal = self._subtree_counts(cur_total, cur_dim, cur_counts, j)
                        if v == el:
                            stack.append((cur_num, cur_total, cur_classes, cur_counts.copy()))
                            cur_total = subtotal
                            cur_dim -= 1
                            cur_counts[j] -= 1
                            if cur_counts[j] <= 0: # just to be safe because why not
                                cur_classes = np.delete(cur_classes, j)
                                cur_counts = np.delete(cur_counts, j)
                            break
                        else:
                            num_before += subtotal
                    # short circuit if we've gotten down to a terminal node where
                    # there is just one unique element
                    if len(cur_classes) == 1:
                        # print("    +", "{:>2}".format(num_before), i, j, cur_total)
                        tup = (cur_num, cur_total, cur_classes, cur_counts)
                        for x in range(len(state) - (i+1)):
                            stack.append(tup)
                        break
                inds[sn] = num_before
                # print("    =", "{:>2}".format(num_before))

            return inds

        def from_perm_indices(self, inds, assume_sorted=False):
            """
            Just loops through the unique permutations
            and returns the appropriate ones for inds.
            Done all in one call for efficiency reasons
            :param inds:
            :type inds: np.ndarray
            :return: permutation array
            :rtype: np.ndarray
            """

            if len(inds) == 0:
                return np.array([], dtype='int8')

            if not assume_sorted:
                sorting = np.argsort(inds)
                inds = inds[sorting]
            else:
                sorting = None

            perms = []
            for n, p in enumerate(self._unique_permutations(self.partition)):
                while n == inds[0]:
                    perms.append(p)
                    inds = inds[1:]
                    if len(inds) == 0:
                        break
                if len(inds) == 0:
                    break
            if len(inds) > 0:
                raise ValueError("indices {} are beyond the number of permutations supported by {}".format(
                    inds,
                    self
                ))

            perms = np.array(perms, dtype='int8')
            if sorting is not None:
                perms = perms[np.argsort(sorting)]

            return perms

        def __repr__(self):
            return "{}({}, ndim={}, states={})".format(
                type(self).__name__,
                self.partition[np.where(self.partition != 0)],
                self.ndim,
                self.total_states
            )
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
            indexers = [ self.PartitionPermutationIndexer(x) for x in partitions[sorting,] ]
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
        Currently pretty clumsy.
        We just search and sum for each state.

        We're going to move to a more intelligent process
        where we initially sort and group the states by their associated integer partition
        and then reverse lexicographically so that we can do blocks of states at once.
        This especially allows us to be faster inside the actual permutation indexer, since
        there when we descend the indexing tree we usually don't need to go all the way back
        up for every successive state

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
                # max_part = max(len(x.classes) for x in partitions)
                sort_classes = np.array([
                    np.concatenate([-np.flip(s), np.flip(z), [c[0]]])
                    for c, s, z in zip(state_dat, states, np.sort(states, axis=-1))
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
                        block_inds = idxer.get_perm_indices(idxer_block)
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


