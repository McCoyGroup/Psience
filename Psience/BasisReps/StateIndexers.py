"""
Provides utilities for converting state vectors to indices.
Used implicitly by various conversion utilities
"""

import abc, numpy as np, collections
from McUtils.Combinatorics import SymmetricGroupGenerator

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
        self.symmetric_group = SymmetricGroupGenerator(self.ndim)
        self.partitioning_cache = {}

    def to_state(self, serializer=None):
        return {
            'ndim': self.ndim
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(data['ndim'])

    def to_indices(self, states):
        """
        Finds the appropriate integer partitioning for each state
        :param states: 2D array of states as excitations
        :type states: np.ndarray
        :return:
        :rtype:
        """

        return self.symmetric_group.to_indices(states)

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

        return self.symmetric_group.from_indices(indices)


