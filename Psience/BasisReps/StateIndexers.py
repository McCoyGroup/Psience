"""
Provides utilities for converting state vectors to indices.
Used implicitly by various conversion utilities
"""

import abc, numpy as np

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
    def to_indices(self, states):
        return np.ravel_multi_index(states, self.dims)
    def from_indices(self, indices):
        return np.unravel_index(indices, self.dims)

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
            self.total_states = self._fac_rat(self.counts)

        @staticmethod
        def _fac_rat(counts):
            import math

            subfac = np.prod([math.factorial(x) for x in counts])
            ndim_fac = math.factorial(np.sum(counts))

            return ndim_fac / subfac

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
            return total * (counts[where] / ndim)

        def get_perm_idx(self, state):
            """
            Gets the index for a single state.
            Does this by looping through the state, decrementing the appropriate character class,
            computing the number of nodes in the child tree (done based on the initial total states),
            and then adding up all those terms

            :param state:
            :type state:
            :return:
            :rtype:
            """

            num_before=0
            cur_dim = self.ndim
            cur_total = self.total_states
            cur_classes = np.copy(self.classes)
            cur_counts = np.copy(self.counts)
            # print(state)
            # print(cur_classes, cur_counts, cur_total, cur_dim)
            for i, el in enumerate(state):
                for j, v in enumerate(cur_classes):
                    subtotal = self._subtree_counts(cur_total, cur_dim, cur_counts, j)
                    # print(i, el, subtotal)
                    if v == el:
                        cur_total = subtotal
                        cur_dim -= 1
                        cur_counts[j] -= 1
                        if cur_counts[j] <= 0: # just to be safe because why not
                            cur_classes = np.delete(cur_classes, j)
                            cur_counts = np.delete(cur_counts, j)
                        break
                    else:
                        num_before += subtotal
                if len(cur_classes) == 1:
                    break
            return num_before

    def integer_partitions(self, num):
        """
        Gives basically the second sort criterion by calculating integer partitions
        (i.e. number of way to split up num quanta across modes)
        Sorted by default, which is our saving grace
        :param num:
        :type num:
        :return:
        :rtype:
        """
        if num not in self.partitioning_cache:
            indexers = list(reversed([
                self.PartitionPermutationIndexer(
                    np.array(list(reversed(x)) + [0]*(self.ndim - len(x)))
                ) for x in self._accel_asc(num)
            ]))
            total_dim = np.sum([x.total_states for x in indexers])
            if num > 0:
                prev_dim = self.integer_partitions(num-1)[1]
                total_dim += prev_dim # recursively we get it all
            else:
                prev_dim = 0
            self.partitioning_cache[num] = (prev_dim, total_dim, indexers)

        return self.partitioning_cache[num]

    def _get_inds_by_quanta(self, states, quanta):
        """
        Currently pretty clumsy.
        We just search and sum for each state.
        :param states:
        :type states:
        :param quanta:
        :type quanta:
        :return:
        :rtype:
        """
        # sorting metric
        baseline, _, partitions = self.integer_partitions(quanta)
        inds = np.full(len(states), -1)
        for i, s in enumerate(states):
            idx = baseline
            cls, cnt = np.unique(s, return_counts=True)
            cls = np.flip(cls)
            cnt = np.flip(cnt)
            for indexer in partitions:
                icts = indexer.counts #type: np.ndarray
                icls = indexer.classes  # type: np.ndarray
                if (
                        len(icts) == len(cnt)
                        and (icts == cnt).all()
                        and (icls == cls).all()
                ): # it's a match!
                    idx += indexer.get_perm_idx(s)
                    break
                else:
                    # add on number of skipped tree nodes
                    idx += indexer.total_states
            inds[i] = idx
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
        nquants = np.sum(states, axis=1)
        inds = np.full_like(nquants, -1) # set up return value
        for q in np.unique(nquants):
            where_are_we = np.where(nquants == q)[0]
            these_states = states[where_are_we]
            inds[where_are_we] = self._get_inds_by_quanta(these_states, q)
        if smol:
            inds = inds[0]
        return inds

    def from_indices(self, indices):
        """
        Inverts the index calculation...not implemented yet
        :param indices:
        :type indices:
        :return:
        :rtype:
        """
        raise NotImplementedError("womp")




