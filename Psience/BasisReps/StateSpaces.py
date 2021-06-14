"""
Provides a relatively haphazard set of simple classes to keep track of state information.
By providing a single interface here, we can avoid recomputing information over and over.
"""

import numpy as np, itertools as ip, enum, scipy.sparse as sp
import abc

from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Scaffolding import MaxSizeCache
from McUtils.Combinatorics import SymmetricGroupGenerator, UniquePermutations

__all__ = [
    "AbstractStateSpace",
    "BasisStateSpace",
    "BasisMultiStateSpace",
    "SelectionRuleStateSpace",
    "BraKetSpace",
    "StateSpaceMatrix"
]

class AbstractStateSpace(metaclass=abc.ABCMeta):
    """
    Represents a generalized state space which will provide core
    methods to index into a basis and generate representations
    """

    class StateSpaceSpec(enum.Enum):
        Excitations = "excitations"
        Indices = "indices"

    # for flexibility later.
    # we make use of explicitly narrowed dtypes
    # to cut down on the memory footprint
    excitations_dtype = np.dtype('int8')
    indices_dtype = np.dtype('int64')

    def __init__(self, basis):
        """
        :param basis:
        :type basis: RepresentationBasis
        """
        self.basis = basis
        self._indices = None
        self._excitations = None
        self._indexer = None
        self._uindexer = None
        self._exc_indexer = None
        self._uexc_indexer = None
        self._uinds = None
        self._sort_uinds = None

    @property
    def ndim(self):
        return self.basis.ndim

    def _pull_exc(self):
        if self._excitations is None:
            self._excitations = self.as_excitations()
        return self._excitations
    @property
    def excitations(self):
        return self._pull_exc()
    @excitations.setter
    def excitations(self, exc):
        self._excitations = exc

    @property
    def mode(self):
        return self.get_mode()
    @abc.abstractmethod
    def get_mode(self):
        """
        Returns the mode (indices or excitations) for the held states

        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")
    @property
    def has_indices(self):
        return self._indices is not None
    @property
    def has_excitations(self):
        return self._excitations is not None

    def _pull_inds(self):
        if self._indices is None:
            self._indices = self.as_indices()
        return self._indices
    @property
    def indices(self):
        return self._pull_inds()
    @indices.setter
    def indices(self, inds):
        self._indices = inds

    @property
    def indexer(self):
        if self._indexer is None:
            self._indexer = nput.argsort(self.indices)
        return self._indexer
    @indexer.setter
    def indexer(self, idxer):
        self._indexer = idxer

    @indexer.setter
    def indexer(self, idxer):
        self._indexer = idxer

    @property
    def exc_indexer(self):
        if self._exc_indexer is None:
            self._exc_indexer = nput.argsort(self.excitations)
        return self._exc_indexer
    @exc_indexer.setter
    def exc_indexer(self, idxer):
        self._exc_indexer = idxer



    def find(self, to_search, check=True):
        """
        Finds the indices of a set of indices inside the space

        :param to_search: array of ints
        :type to_search: np.ndarray | AbstractStateSpace
        :return:
        :rtype:
        """
        if not isinstance(to_search, np.ndarray):
            if isinstance(to_search, AbstractStateSpace) or hasattr(to_search, 'indices'):
                to_search = to_search.indices
            else:
                to_search = np.array(to_search)
        if to_search.ndim > 1:
            raise ValueError("currently only accept subspaces as indices or AbstractStateSpaces")
        vals = np.searchsorted(self.indices, to_search, sorter=self.indexer)
        if isinstance(vals, (np.integer, int)):
            vals = np.array([vals])
        # we have the ordering according to the _sorted_ version of `indices`
        # so now we need to invert that back to the unsorted version
        if len(self.indexer) > 0:
            big_vals = vals == len(self.indexer)
            vals[big_vals] = -1
            vals = self.indexer[vals]
            # now because of how searchsorted works, we need to check if the found values
            # truly agree with what we asked for...although I could maybe be smarter on this...
            bad_vals = self.indices[vals] != to_search
            if vals.shape == ():
                if bad_vals:
                    vals = -1
            else:
                vals[bad_vals] = -1
        else:
            bad_vals = np.full_like(to_search, True)
            vals = np.full_like(vals, -1)
        if check and bad_vals.any():
            raise IndexError("{} not in {}".format(to_search[bad_vals], self))
        return vals

    def find(self, to_search, check=True):
        """
        Finds the indices of a set of indices inside the space

        :param to_search: array of ints
        :type to_search: np.ndarray | AbstractStateSpace
        :return:
        :rtype:
        """
        if not isinstance(to_search, np.ndarray):
            if isinstance(to_search, AbstractStateSpace) or hasattr(to_search, 'indices'):
                to_search = to_search.indices
            else:
                to_search = np.asanyarray(to_search)
        # if to_search.ndim > 1:
        #     raise ValueError("currently only accept subspaces as indices or AbstractStateSpaces")

        if to_search.ndim == 0:
            to_search = np.array([to_search], dtype=to_search.dtype)
        if to_search.ndim == 1:
            vals, _ = nput.find(self.indices, to_search, sorting=self.indexer, check=check)
        else:
            vals, self._exc_indexer = nput.find(self.excitations, to_search, sorting=self._exc_indexer, check=check)

        return vals

    def __len__(self):
        if self._indices is not None:
            return len(self.indices)
        else:
            return len(self.excitations)

    @property
    def unique_len(self):
        if self._indices is not None:
            return len(self.unique_indices)
        else:
            return len(self.unique_excitations)

    @property
    def unique_indices(self):
        """
        Returns the unique indices
        :return:
        :rtype:
        """
        return self.as_unique_indices()

    @property
    def unique_excitations(self):
        """
        Returns the unique excitations
        :return:
        :rtype:
        """
        return self.as_unique_excitations()

    @abc.abstractmethod
    def as_indices(self):
        """
        Returns the index version of the stored states
        :return:
        :rtype: np.ndarray
        """
        return NotImplementedError("abstract base class")

    def as_unique_indices(self, sort=False):
        """
        Returns unique indices

        :return:
        :rtype:
        """
        if len(self.indices) == 0:
            return self.indices
        inds = self.indices
        if sort and self._sort_uinds is None:
            _, self._indexer, sort_uinds = nput.unique(inds, sorting=self._indexer, return_index=True)
            self._sort_uinds = sort_uinds
            self._uinds = np.sort(sort_uinds)
        elif not sort and self._uinds is None:
            _, self._indexer, sort_uinds = nput.unique(inds, sorting=self._indexer, return_index=True)
            self._sort_uinds = sort_uinds
            self._uinds = np.sort(sort_uinds)

        if sort:
            return inds[self._sort_uinds,]
        else:
            return inds[self._uinds,]

    @abc.abstractmethod
    def as_excitations(self):
        """
        Returns the excitation version of the stored states
        :return:
        :rtype: np.ndarray
        """
        return NotImplementedError("abstract base class")

    def as_unique_excitations(self, sort=False):
        """
        Returns unique excitations
        :return:
        :rtype:
        """
        if len(self.excitations) == 0:
            return self.excitations
        exc = self.excitations
        if sort and self._sort_uinds is None:
            inds = self.as_unique_indices(sort=sort)
        elif not sort and self._uinds is None:
            _, self._exc_indexer, uinds = nput.unique(exc, axis=0, sorting=self._exc_indexer, return_index=True)
            self._uinds = np.sort(uinds)

        if sort:
            return exc[self._sort_uinds,]
        else:
            return exc[self._uinds,]

    @abc.abstractmethod
    def get_representation_indices(self,
                                   other=None,
                                   selection_rules=None,
                                   freqs=None,
                                   freq_threshold=None
                                   ):
        """
        Returns bra and ket indices that can be used as indices to generate representations

        :param other:
        :type other:
        :param selection_rules:
        :type selection_rules:
        :param freqs:
        :type freqs:
        :param freq_threshold:
        :type freq_threshold:
        :return:
        :rtype: (np.ndarray, np.ndarray)
        """
        return NotImplementedError("abstract base class")

    @abc.abstractmethod
    def get_representation_brakets(self,
                                   other=None,
                                   selection_rules=None,
                                   freqs=None,
                                   freq_threshold=None
                                   ):
        """
        Returns a BraKetSpace that can be used as generate representations

        :param other:
        :type other:
        :param selection_rules:
        :type selection_rules:
        :param freqs:
        :type freqs:
        :param freq_threshold:
        :type freq_threshold:
        :return:
        :rtype: BraKetSpace
        """
        return NotImplementedError("abstract base class")

    @abc.abstractmethod
    def take_states(self, states):
        """
        Takes the intersection of self and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")
    @abc.abstractmethod
    def take_subspace(self, sel):
        """
        Takes a subset of the states

        :param sel:
        :type sel:
        :return:
        :rtype: AbstractStateSpace
        """
        raise NotImplementedError("abstract base class")
    @abc.abstractmethod
    def take_subdimensions(self, inds):
        """
        Takes a subset of the state dimensions

        :param sel:
        :type sel:
        :return:
        :rtype: AbstractStateSpace
        """
        raise NotImplementedError("abstract base class")

    @abc.abstractmethod
    def drop_states(self, states):
        """
        Takes the difference of self and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")
    @abc.abstractmethod
    def drop_subspace(self, sel):
        """
        Drops a subset of the states

        :param sel:
        :type sel:
        :return:
        :rtype: AbstractStateSpace
        """
        raise NotImplementedError("abstract base class")
    @abc.abstractmethod
    def drop_subdimensions(self, inds):
        """
        Drops a subset of the state dimensions

        :param sel:
        :type sel:
        :return:
        :rtype: AbstractStateSpace
        """
        raise NotImplementedError("abstract base class")

    @classmethod
    def get_states_with_quanta(cls, n, ndim):
        """
        Returns the states with number of quanta equal to n

        :param quanta:
        :type quanta:
        :return:
        :rtype:
        """

        # generate basic padded partitions of `n`
        partitions = SymmetricGroupGenerator(ndim).get_terms(n, flatten=True)

        # then concatenate unique permutations for each and cast to numpy
        return np.array(
            partitions,
            dtype=cls.excitations_dtype
        )

    @abc.abstractmethod
    def to_single(self):
        """
        Flattens any complicated state space structure into a
        single space like a `BasisStateSpace`

        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")

    def split(self, chunksize):
        """
        Subclass overridable function to allow for spaces to be
        split up into chunks
        :param chunksize:
        :type chunksize:
        :return:
        :rtype:
        """
        raise TypeError("{} can't be split".format(type(self).__name__))

class BasisStateSpace(AbstractStateSpace):
    """
    Represents a subspace of states inside a representation basis.
    Useful largely to provide consistent, unambiguous representations of multiple states across
    the different representation-generating methods in the code base.
    """

    def __init__(self, basis, states, mode=None):
        """
        :param basis:
        :type basis: RepresentationBasis
        :param states:
        :type states: Iterable[int]
        :param mode: whether the states were supplied as indices or as excitations
        :type mode: None | str | StateSpaceSpec
        """

        super().__init__(basis)

        self._init_states = np.asanyarray(states)#, dtype=int)
        if mode is not None and not isinstance(mode, self.StateSpaceSpec):
            mode = self.StateSpaceSpec(mode)
        self._init_state_types = mode
        self._indices = None
        self._excitations = None
        if len(self._init_states) > 0:
            if self.infer_state_inds_type() == self.StateSpaceSpec.Indices:
                self._indices = self._init_states.astype(int)
            else:
                self._excitations = self._init_states.astype(self.excitations_dtype)
        else:
            if mode is None:
                self._init_state_types = self.StateSpaceSpec.Indices

        self._sorted = None
        self._indexer = None

    def to_state(self, serializer=None):
        return {
            'basis':self.basis,
            'indices':self.indices
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(
            serializer.deserialize(data['basis']),
            serializer.deserialize(data['indices']),
            mode=cls.StateSpaceSpec.Indices)

    @property
    def indices(self):
        """
        Returns held indices

        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        return self._pull_inds()
    @indices.setter
    def indices(self, inds):
        """
        Sets indices
        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        self._indices = np.asanyarray(inds).astype(self.indices_dtype)

    @property
    def excitations(self):
        """
        Returns held excitations
        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        return self._pull_exc()
    @excitations.setter
    def excitations(self, exc):
        """
        Sets held excitations
        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        self._excitations = np.asanyarray(exc).astype(self.excitations_dtype)

    @classmethod
    def from_quanta(cls, basis, quants):
        """
        Returns states with `quants` quanta of excitation
        using the basis `basis`

        :param basis:
        :type basis: RepresentationBasis
        :param quants: set of integers
        :type quants: int | Iterable[int]
        :return: BasisStateSpace
        :rtype:
        """

        if isinstance(quants, int):
            quants = [quants]

        states = np.concatenate(
            [cls.get_states_with_quanta(n, basis.ndim) for n in quants]
        )

        return cls(basis, states, mode=cls.StateSpaceSpec.Excitations)

    def get_mode(self):
        return self.StateSpaceSpec.Indices if self.has_indices else self.StateSpaceSpec.Excitations

    def infer_state_inds_type(self):
        if self._init_state_types is not None:
            return self._init_state_types
        else:
            end_shape = self._init_states.shape[-1]
            ndim = self.basis.ndim
            if end_shape == ndim:
                return self.StateSpaceSpec.Excitations
            else:
                return self.StateSpaceSpec.Indices

    def as_excitations(self):
        """
        Returns states as sets of excitations, rather than indices indo the basis functions.
        For 1D, this just means converting a list of states into tuples of length 1.

        :param states:
        :type states:
        :return:
        :rtype:
        """

        states = self._init_states
        states_type = self.infer_state_inds_type()
        if len(states) == 0:
            return np.array([], dtype=self.excitations_dtype)

        if states_type is self.StateSpaceSpec.Excitations:
            return np.reshape(states, (-1, self.ndim))
        elif states_type is self.StateSpaceSpec.Indices:
            states, self._indexer, uinds, inv = nput.unique(states, sorting=self._indexer, return_index=True, return_inverse=True)
            self._sort_uinds = uinds
            self._uinds = np.sort(uinds)
            raw_exc = self.basis.unravel_state_inds(states.flatten())
            return raw_exc[inv]
        else:
            raise ValueError("don't know what to do with state spec {}".format(
                states_type
            ))

    def as_indices(self):
        """
        Returns states as sets of excitations, rather than indices indo the basis functions.
        For 1D, this just means converting a list of states into tuples of length 1.

        :param states:
        :type states:
        :return:
        :rtype:
        """

        states = self._init_states
        states_type = self.infer_state_inds_type()
        if states_type is self.StateSpaceSpec.Excitations:
            states, self._exc_indexer, uinds, inv = nput.unique(states,
                                                                sorting=self._exc_indexer,
                                                                return_index=True,
                                                                axis=0,
                                                                return_inverse=True
                                                                )
            self._uinds = np.sort(uinds)
            # if len(states) > 1000:
            #     raise Exception('why')
            raw_inds = self.basis.ravel_state_inds(np.reshape(states, (-1, self.ndim)))
            return raw_inds[inv]
        elif states_type is self.StateSpaceSpec.Indices:
            return states.flatten()
        else:
            raise ValueError("don't know what to do with state spec {}".format(
                states_type
            ))

    def to_single(self):
        """
        Basically a no-op
        :return:
        :rtype:
        """
        return self

    def is_unique(self):
        """
        Returns `True` if the number of states is equal to number of unique states
        :return:
        :rtype:
        """
        if self.mode == self.StateSpaceSpec.Indices:
            return len(self.unique_indices) == len(self.indices)
        else:
            return len(self.unique_excitations) == len(self.excitations)

    def is_sorted(self, allow_indeterminate=True):
        """
        Checks and then sets a flag
        :return:
        :rtype:
        """
        if self._sorted is None:
            if allow_indeterminate and self._indexer is None:
                return None
            else:
                self._sorted = (self.indexer == np.arange(len(self.indexer))).all()
        return self._sorted

    def take_unique(self, sort=False, use_indices=False):
        """
        Returns only the unique states, but preserves
        ordering and all of that unless explicitly allowed not to
        :return:
        :rtype:
        """
        if self.is_unique():
            if not sort:
                return self
            elif self.is_sorted():
                return self

        if use_indices or sort or (self.mode == self.StateSpaceSpec.Indices):
            states = self.as_unique_indices(sort=sort)
            spec = self.StateSpaceSpec.Indices
            new = type(self)(self.basis, states, mode=spec)
            if self._excitations is not None:
                new.excitations = self.as_unique_excitations(sort=sort)
            if sort or self.is_sorted(allow_indeterminate=True):
                # we have strict ordering relationships we can use since we know the
                # new states are also sorted
                new._indexer = np.arange(len(states))
                new._sorted = True
        else:
            states = self.unique_excitations
            spec = self.StateSpaceSpec.Excitations
            new = type(self)(self.basis, states, mode=spec)
            if (self.mode == self.StateSpaceSpec.Indices):
                new.indices = self.unique_indices
        return new

    def as_sorted(self):
        """
        Returns only the unique states, but preserves
        ordering and all of that unless explicitly allowed not to
        :return:
        :rtype:
        """

        indexer = self.indexer
        spec = self.StateSpaceSpec.Indices
        states = self.indices[indexer]
        new = type(self)(self.basis, states, mode=spec)
        new._sorted = True
        new._indexer = np.arange(len(states))

        if self._uinds is not None:
            # we relate the old positions where the unique elements happened to the
            # new ones by argsorting the current indexer so we can know "where" something
            # went and then pulling the uinds from that
            where_map = np.argsort(indexer)
            if self._sort_uinds is not None:
                new._uinds = where_map[new._uinds]

        if self.has_excitations():
            new.excitations = self.excitations[indexer,]
            if self._exc_indexer is not None:
                # we relate the old sorting to the new sorting by
                # converting reordering the old sorting and then
                # sorting that
                new._exc_indexer = np.argsort(self._exc_indexer[indexer])

        return new

    def apply_selection_rules(self, selection_rules, filter_space=None, iterations=1):
        """
        Generates a new state space from the application of `selection_rules` to the state space.
        Returns a `BasisMultiStateSpace` where each state tracks the effect of the application of the selection rules
        up to the number of iteration specified.

        :param basis:
        :type basis:
        :param selection_rules:
        :type selection_rules:
        :param states:
        :type states:
        :param iterations:
        :type iterations:
        :param filter_space:
        :type filter_space:
        :return:
        :rtype: SelectionRuleStateSpace
        """

        return SelectionRuleStateSpace.from_rules(self, selection_rules, filter_space=filter_space, iterations=iterations)

    def get_representation_indices(self,
                                   other=None,
                                   selection_rules=None,
                                   freqs=None,
                                   freq_threshold=None,
                                   filter=None,
                                   return_filter=False
                                   ):
        """
        Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
        in this state space.
        Basically just takes all pairs of indices.
        Only returns the upper-triangle indices

        :return:
        :rtype:
        """

        # if other is None:
        #     other = self

        filter = None
        if freq_threshold is None:
            if selection_rules is None:
                # TODO: not sure I want to be happening here
                l_inds = self.indices
                r_inds = other.indices
                # TODO: make this less slow... ip.product can be brutal
                pairs = np.array(list(ip.product(l_inds, r_inds)))
                _, upos = np.unique(pairs, axis=0, return_index=True)
                m_pairs = pairs[np.sort(upos)].T
            else:
                # Get the representation indices that can be coupled under the supplied set of selection rules
                # Currently this is clumsy.
                # We do this by computing transformed states finding where this intersects with the other space
                if filter is None:
                    filter = other
                transf = self.apply_selection_rules(selection_rules, filter_space=filter)
                if not isinstance(transf, AbstractStateSpace):
                    transf, filter = transf
                m_pairs = transf.get_representation_indices()

        else:
            raise NotImplementedError("Changed up how comb will apply, but haven't finished the implementation")
            exc = self.excitations
            m = np.array(list(ip.combinations(m, 2)))
            # now if we have a frequency comb, we apply it

            if freqs is None:
                raise ValueError("to apply frequency difference threshold, need harmonic frequencies")
            # state_freq = np.sum(freqs[np.newaxis, :]*(states + 1/2), axis=1)

            for i, c in enumerate(h1_couplings):
                s, h1_terms = c
                freq = np.sum(freqs * (s + 1 / 2), axis=1)
                h1_freq = np.sum(freqs[np.newaxis, :] * (h1_terms + 1 / 2), axis=1)
                diffs = h1_freq - freq
                thresh = diffs < freq_threshold
                # for d in h1_freq_diffs[1:]:
                #     # if any of the coupled states is below the threshold for any of the target states we include it
                #     # in our basis
                #     h1_thresh = np.logical_or(h1_thresh, d < freq_threshold)
                h1_sel = np.where(thresh)
                h1_terms = h1_terms[h1_sel]
                h1_couplings[i] = (s, h1_terms)

        if return_filter:
            return m_pairs, filter
        else:
            return m_pairs

    def get_representation_brakets(self,
                                   other=None,
                                   selection_rules=None,
                                   freqs=None,
                                   freq_threshold=None,
                                   filter=None,
                                   return_filter=False
                                   ):
        """
        Generates a `BraKetSpace` that can be fed into a `Representation`
        Only returns the upper-triangle pairs because we assume symmetry

        :return:
        :rtype:
        """

        #TODO: this doesn't _really_ need the indices...
        #       and it'd be faster if I can do this without them

        # if (self.has_indices and (other is None or other.has_indices) or (
        #         self.has_indices
        #         and not (self.has_excitations and other.has_excitations)
        #         and len(self) > len(other)
        # ):
        inds = self.get_representation_indices(other=other,
                                               selection_rules=selection_rules,
                                               freqs=freqs,
                                               freq_threshold=freq_threshold,
                                               filter=filter,
                                               return_filter=return_filter
                                               )

        if return_filter:
            inds, filter = inds
        else:
            filter = None

        if len(inds) > 0:
            row_space = BasisStateSpace(self.basis, inds[0], mode=self.StateSpaceSpec.Indices)
            col_space = BasisStateSpace(self.basis, inds[1], mode=self.StateSpaceSpec.Indices)
            bras = self.to_single().take_states(row_space)
            if other is not None:
                kets = other.to_single().take_states(col_space)
            else:
                kets = self.to_single().take_states(col_space)
        else:
            bras = BasisStateSpace(self.basis, [])
            kets = BasisStateSpace(self.basis, [])

        if return_filter:
            # try:
            return BraKetSpace(bras, kets), filter
            # except:
            #     coops = nput.contained(col_space.indices, other.indices, invert=True)[0]
            #     meh = inds[0][coops]
            #     raise Exception(
            #         nput.difference(kets.excitations, other.excitations)[0],
            #         meh,
            #         self.basis.indexer.symmetric_group.from_indices(meh[:20]),
            #         other.indices[:20]
            #     )
        else:
            return BraKetSpace(bras, kets)

    def take_subspace(self, sel, assume_sorted=False):
        """
        Returns a subsample of the space.
        Intended to be a cheap operation, so samples
        along either the indices or the excitations, depending
        on which we have
        If we know the subsample is sorted then we can actually reuse more information
        and so we make use of that

        :param sel:
        :type sel:
        :return:
        :rtype:
        """

        if assume_sorted and not self.is_sorted():
            return self.as_sorted().take_subspace(sel, assume_sorted=True)

        if self.has_excitations:
            subspace = type(self)(
                self.basis,
                self.excitations[sel,],
                mode=self.StateSpaceSpec.Excitations
            )
            if self.has_indices:
                subspace.indices = self.indices[sel,]
            if assume_sorted: #from earlier directly implies `is_sorted`
                raise NotImplementedError('need to account for the implications of this...')
        else:
            subspace = type(self)(
                self.basis,
                self.indices[sel,],
                mode=self.StateSpaceSpec.Indices
            )

        return subspace
    def take_subdimensions(self, inds):
        """
        Returns a subsample of the space with some dimensions
        dropped
        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        return type(self)(
            self.basis.take_subdimensions(inds),
            self._excitations[:, inds],
            mode=self.StateSpaceSpec.Excitations
        )
    def take_states(self, states, sort=False, assume_sorted=False):
        """
        Takes the set of specified states from the space.
        A lot like take_subspace, but operates on states, not indices
        :param states:
        :type states:
        :return:
        :rtype:
        """
        found = self.find(states, check=False)
        sel = found[found >= 0]
        if sort and not assume_sorted:
            sel = np.sort(sel)
            assume_sorted = True
        return self.take_subspace(sel, assume_sorted=assume_sorted)

    def drop_subspace(self, sel):
        """
        Returns a subsample of the space.
        Intended to be a cheap operation, so samples
        along either the indices or the excitations, depending
        on which we have

        :param sel:
        :type sel:
        :return:
        :rtype:
        """

        if self._excitations is not None:
            diff = np.setdiff1d(np.arange(len(self.excitations)), sel)
        else:
            diff = np.setdiff1d(np.arange(len(self.indices)), sel)
        return self.take_subspace(diff)
    def drop_subdimensions(self, inds):
        """
        Returns a subsample of the space with some dimensions
        dropped
        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        diff = np.setdiff1d(np.arange(self.basis.ndim), inds)
        return self.take_subdimensions(diff)
    def drop_states(self, states):
        """
        Takes the set of specified states from the space.
        A lot like take_subspace, but operates on states, not indices
        :param states:
        :type states:
        :return:
        :rtype:
        """
        found = self.find(states, check=False)
        sel = found[found >= 0]
        return self.drop_subspace(sel)

    def split(self, chunksize):
        """
        Splits the space up into chunks of at max chunksize
        :param chunksize:
        :type chunksize: int
        :return:
        :rtype: Iterable[BasisStateSpace]
        """

        nblocks = np.ceil(len(self) / chunksize)

        if self._indices is not None:
            ind_chunks = np.array_split(self._indices, nblocks, axis=0)
            if self._excitations is not None:
                exc_chunks = np.array_split(self._excitations, nblocks, axis=0)
            else:
                exc_chunks = [None] * len(ind_chunks)
        else:
            exc_chunks = np.array_split(self._excitations, nblocks, axis=0)
            ind_chunks = [None] * len(exc_chunks)

        spaces = []
        for i,e in zip(ind_chunks, exc_chunks):
            if i is None:
                new = type(self)(self.basis, e, mode = self.StateSpaceSpec.Excitations)
            else:
                new = type(self)(self.basis, i, mode=self.StateSpaceSpec.Indices)
                if e is not None:
                    new.excitations = e
            spaces.append(new)
        return spaces

    def union(self, other, sort=False, use_indices=True):
        """
        Returns a merged version of self and other, making
        use of as much of the information inherent in both as is possible

        :param other:
        :type other: BasisStateSpace
        :return:
        :rtype:
        """

        if self.basis is not other.basis:
            raise ValueError("can't merge state spaces over different bases ({} and {})".format(
                self.basis,
                other.basis
            ))

        if use_indices or sort or (self.has_indices and other.has_indices) or(
                self.has_indices
                and not (self.has_excitations and other.has_excitations)
                and len(self) > len(other)
        ): # no need to be wasteful and recalc stuff, right?
            # create merge based on indices and then
            # secondarily based on excitations if possible
            self_inds = self.unique_indices
            other_inds = other.unique_indices

            if len(self_inds) == 0:
                return other.take_unique(sort=sort)
                # new_inds = other_inds
            elif len(other_inds) == 0:
                return self.take_unique(sort=sort)
                # new_inds = self_inds
            else:
                new_inds = np.concatenate([self_inds, other_inds], axis=0)
            _, _, uinds = nput.unique(new_inds, return_index=True)
            # if the number of unique states hasn't increased,
            # one of the spaces is contained in the other
            if len(uinds) == len(self.unique_indices):
                return self.take_unique(sort=sort)
            elif len(uinds) == len(other.unique_indices):
                return other.take_unique(sort=sort)
            else:
                if not sort:
                    sorting = nput.argsort(uinds)
                    indexer = nput.argsort(sorting)
                    uinds = uinds[sorting]
                else:
                    indexer = np.arange(len(uinds))

                new_inds = new_inds[uinds,]
                new = BasisStateSpace(self.basis, new_inds, mode=self.StateSpaceSpec.Indices)
                new._indexer = indexer
                new._uinds = np.arange(len(uinds))

                if self._excitations is not None or other._excitations is not None:
                    self_exc = self.unique_excitations
                    other_exc = other.unique_excitations
                    if len(self_exc) == 0:
                        new_exc = other_exc
                    elif len(other_exc) == 0:
                        new_exc = self_exc
                    else:
                        new_exc = np.concatenate([self_exc, other_exc], axis=0)
                    new._excitations = new_exc[uinds,]
        else:
            # create merge based on excitations and then
            # secondarily based on indices if possible
            self_exc = self.unique_excitations
            other_exc = other.unique_excitations
            if len(self_exc) == 0:
                return other.take_unique(sort=sort)
                # new_exc = other_exc
            elif len(other_exc) == 0:
                return self.take_unique(sort=sort)
                # new_exc = self_exc
            else:
                new_exc = np.concatenate([self_exc, other_exc], axis=0)

            _, uinds = np.unique(new_exc, axis=0, return_index=True)

            # if the number of unique states hasn't increased,
            # one of the spaces is contained in the other
            if len(uinds) == len(self.unique_excitations):
                new = self.take_unique(sort=sort)
            elif len(uinds) == len(other.unique_excitations):
                new = other.take_unique(sort=sort)
            else:
                sorting = np.argsort(uinds)
                uinds = uinds[sorting]
                new_exc = new_exc[uinds,]

                new = BasisStateSpace(self.basis, new_exc, mode=self.StateSpaceSpec.Excitations)
                new._exc_indexer = nput.argsort(sorting)
                new._uinds = np.arange(len(uinds))

                if other._indices is not None:
                    self_inds = self.indices
                    other_inds = other.indices
                    new_inds = np.concatenate([self_inds, other_inds], axis=0)
                    new._indices = new_inds[uinds,]

        return new

    def intersection(self, other, sort=False, use_indices=False):
        """
        Returns an intersected self and other

        :param other:
        :type other: BasisStateSpace
        :return:
        :rtype:
        """

        if self.basis is not other.basis:
            raise ValueError("can't merge state spaces over different bases ({} and {})".format(
                self.basis,
                other.basis
            ))

        # create intersection based on indices and then
        # make use of this subselection to resample the basis

        if use_indices or sort or (self.has_indices and other.has_indices) or (
                self.has_indices
                and not (self.has_excitations and other.has_excitations)
                and len(self) > len(other)
        ): # no need to be wasteful and recalc stuff, right?
            # unfortunately no way to only use Excitation data here...
            self_inds = self.unique_indices
            other_inds = other.unique_indices

            if sort:
                new_inds, sortings, _ = nput.intersection(
                    self_inds, other_inds,
                    sortings=(self._uindexer, other._uindexer),
                    assume_unique=True
                )
            else:
                new_inds, sortings, _,  x_inds, _ = nput.intersection(
                    self_inds, other_inds,
                    sortings=(self._uindexer, other._uindexer),
                    assume_unique=True, return_indices=True
                )
            self._uindexer, other._uindexer = sortings
            # now check to make sure we're not being wasteful and destroying an object
            # that doesn't need to be destroyed
            if len(new_inds) == len(self_inds):
                if self.is_unique:
                    return self
                else:
                    return self.take_unique(sort=sort)
            elif len(new_inds) == len(other_inds):
                if other.is_unique:
                    return other
                else:
                    return other.take_unique(sort=sort)
            else:
                if sort:
                    new = self.take_states(new_inds)
                    new._indexer = np.arange(len(new_inds))
                    new._uinds = np.arange(len(new_inds))
                    new._sort_uinds = np.arange(len(new_inds))
                    return new
                else:
                    new = self.take_subspace(x_inds)
                    new._uinds = np.arange(len(new_inds))
                    return new
        else:
            self_exc = self.unique_excitations
            other_exc = other.unique_excitations

            _, sortings, _, new_inds, _ = nput.intersection(
                self_exc, other_exc,
                sortings=(self._uexc_indexer, other._uexc_indexer),
                assume_unique=True, return_indices=True
            )
            self._uexc_indexer, other._uexc_indexer = sortings

            # now check to make sure we're not being wasteful and destroying an object
            # that doesn't need to be destroyed
            if len(new_inds) == len(self_exc):
                if self.is_unique:
                    return self
                else:
                    return self.take_unique()
            elif len(new_inds) == len(other_exc):
                if other.is_unique:
                    return other
                else:
                    return other.take_unique()
            else:
                new_inds = np.sort(new_inds)
                new = self.take_subspace(new_inds, assume_sorted=True)
                new._uinds = np.arange(len(new_inds))
                return new

    def difference(self, other, sort=False, use_indices=False):
        """
        Returns an diff'ed self and other

        :param other:
        :type other: BasisStateSpace
        :return:
        :rtype:
        """

        if self.basis is not other.basis:
            raise ValueError("can't take a difference of state spaces over different bases ({} and {})".format(
                self.basis,
                other.basis
            ))

        # create intersection difference on indices and then
        # make use of this subselection to resample the basis
        if use_indices or sort or (self.has_indices and other.has_indices) or (
                self.has_indices and not (
                    self.has_excitations and other.has_excitations
                ) and len(self) > len(other)
        ):  # no need to be wasteful and recalc stuff, right?
            self_inds = self.unique_indices
            other_inds = other.unique_indices

            new_inds, sortings, merge_sorting = nput.difference(
                self_inds, other_inds,
                sortings=(self._uindexer, other._uindexer),
                assume_unique=True
            )
            self._uindexer, other._uindexer = sortings

            # now we check that we're not destroying an object that can be
            # reused
            if len(new_inds) == len(self_inds):
                if self.is_unique:
                    return self
                else:
                    return self.take_unique(sort=sort)
            else:
                if sort:
                    new = self.take_states(new_inds)
                    new._indexer = np.arange(len(new_inds))
                    new._uinds = np.arange(len(new_inds))
                    new._sort_uinds = np.arange(len(new_inds))
                    return new
                else:
                    _, _, _, found_inds, _ = nput.intersection(
                        self_inds, new_inds,
                        sortings=(self._uindexer, None),
                        assume_unique=True, return_indices=True
                    )
                    found_inds = np.sort(found_inds)
                    new = self.take_subspace(found_inds)
                    new._uinds = np.arange(len(found_inds))
                    return new
        else:
            self_exc = self.unique_excitations
            other_exc = other.unique_excitations

            self_inds, dtype, _, _ = nput.coerce_dtype(self_exc)
            other_inds, _, _, _ = nput.coerce_dtype(other_exc, dtype=dtype)

            new_inds, sortings = nput.difference(
                self_inds, other_inds,
                sortings=(self._uexc_indexer, other._uexc_indexer),
                assume_unique=True
            )
            self._uexc_indexer, other._uexc_indexer = sortings
            # now check to make sure we're not being wasteful and destroying an object
            # that doesn't need to be destroyed
            if len(new_inds) == len(self_inds):
                if self.is_unique:
                    return self
                else:
                    return self.take_unique()
            else:
                _, _, _, found_inds, _ = nput.intersection(
                    self_inds, new_inds,
                    sortings=(self._uexc_indexer, None),
                    assume_unique=True, return_indices=True
                )
                found_inds = np.sort(found_inds)
                new = self.take_subspace(found_inds)
                new._uinds = np.arange(len(found_inds))
                return new

    def __repr__(self):
        return "{}(nstates={}, basis={})".format(
            type(self).__name__,
            len(self),
            self.basis
        )

    # def __getitem__(self, item):
    #     """
    #     Pulls a sample of the held states
    #     :param item:
    #     :type item:
    #     :return:
    #     :rtype:
    #     """
    #     states = self._init_states[item]
    #     return self.take_states(states)

    def __eq__(self, other):
        """

        :param other:
        :type other:
        :return:
        :rtype:
        """
        if self.basis is not other.basis:
            raise ValueError("can't compare {} and {} with different bases ({} and {})".format(
                type(self).__name__,
                type(other).__name__,
                self.basis,
                other.basis
            ))
        return (
            len(self) == len(other)
            and (self.indices == other.indices).all()
        )

class BasisMultiStateSpace(AbstractStateSpace):
    """
    Represents a collection of `BasisStateSpace` objects.
    This is commonly generated by application of selection rules to a standard `BasisStateSpace`.
    Each of these state spaces is nominally independent of the rest, allowing for combinatorial
    efficiency later down the line.
    """

    def __init__(self, spaces):
        """
        :param spaces: array of `BasisStateSpace` objects
        :type spaces: np.ndarray
        :param selection_rules: array of rules used to generate the subspace
        :type selection_rules: np.ndarray
        """
        self.spaces = np.asanyarray(spaces, dtype=object)
        super().__init__(self.basis)

    def to_state(self, serializer=None):
        return {
            'spaces':self.spaces
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(serializer.deserialize(data['spaces']))

    @property
    def representative_space(self):
        return self.spaces.flatten()[0]

    @property
    def basis(self):
        return self.representative_space.basis
    @basis.setter
    def basis(self, b):
        if b is not self.basis:
            raise NotImplementedError("can't change basis after construction")

    @property
    def ndim(self):
        return self.representative_space.ndim

    @property
    def nstates(self):
        return int(np.product(self.spaces.shape))

    def __iter__(self):
        return iter(self.spaces)
    @property
    def flat(self):
        return self.spaces.flat

    def get_mode(self):
        if all(s.has_indices for s in self.spaces.flat):
            return self.StateSpaceSpec.Indices
        elif all(s.has_excitations for s in self.spaces.flat):
            return self.StateSpaceSpec.Excitations
        else:
            raise ValueError("not sure what to do with mixed-mode state spaces...")

    def as_indices(self):
        """
        Pulls the full set of indices out of all of the
        held spaces and returns them as a flat vector
        :return:
        :rtype:
        """

        # we figure out which held spaces still need indices
        needs_inds = [x for x in self.spaces.flat if (not x.has_indices and len(x) > 0)]
        if len(needs_inds) > 0:
            # then we join these up into one big space and use that to efficiently calculate indices
            single_space = BasisStateSpace(self.basis, np.concatenate([x.excitations for x in needs_inds], axis=0), mode=self.StateSpaceSpec.Excitations)
            full_inds = single_space.indices
            # and finally we reassociate these with the appropriate OG spaces
            s = 0
            for space in needs_inds:
                l = len(space)
                space.indices = full_inds[s:s+l]
                s += l

        ind_arrays = [space.indices for space in self.spaces.flat if len(space) > 0]
        if len(ind_arrays) == 0:
            inds = np.array([], dtype=int)
        else:
            inds = np.concatenate(ind_arrays)
        return inds

    def as_excitations(self):
        """
        Pulls the full set excitations out of all of the
        held spaces and returns them as a flat vector
        :return:
        :rtype:
        """

        # we figure out which held spaces still need excitations
        needs_exc = [x for x in self.spaces.flat if (not x.has_excitations and len(x) > 0)]
        if len(needs_exc) > 0:
            # then we join these up into one big space and use that to efficiently calculate excitations
            single_space = BasisStateSpace(self.basis, np.concatenate([x.indices for x in needs_exc], axis=0),
                                           mode=self.StateSpaceSpec.Indices)
            full_exc = single_space.excitations # TODO: think about how to also capture the unique pos. calc'd inside...
            # and finally we reassociate these with the appropriate OG spaces
            s = 0
            for space in needs_exc:
                l = len(space)
                space.excitations = full_exc[s:s + l]
                s += l

        exc_arrays = [space.excitations for space in self.spaces.flat if len(space) > 0]
        if len(exc_arrays) == 0:
            exc = np.array([], dtype=self.excitations_dtype)
        else:
            exc = np.concatenate(exc_arrays, axis=0)

        return exc

    def __getitem__(self, item):
        it = self.spaces[item]
        if isinstance(it, np.ndarray):
            # multispace
            it = type(self)(it)
        return it

    def get_representation_indices(self,
                                   freqs=None,
                                   freq_threshold=None,
                                   other=None,
                                   selection_rules=None,
                                   return_filter=False
                                   ):
        """
        Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
        in this state space.
        Basically just takes all pairs of indices.

        :return:
        :rtype:
        """

        raise NotImplementedError('unoptimized code path')

        if other is not None:
            raise ValueError("haven't implemented getting indices to a set of 'other' states...")
        if selection_rules is not None:
            raise ValueError("applying selection_rules in the index generation isn't well defined...")

        m_pairs = None
        for space in self.spaces.flat:
            m = space.get_representation_indices(
                freqs=None,
                freq_threshold=None
            ).T
            # if self.selection_rules is not None:
            #     q_changes = np.unique([sum(np.abs(x)) for x in self.selection_rules])
            #     m = self.filter_representation_inds(m, q_changes).T
            # else:
            #     m = m.T
            if m_pairs is None:
                m_pairs = m
            else:
                m_pairs = np.concatenate([m_pairs, m], axis=0)
                _, upos = np.unique(m_pairs, return_index=True, axis=0)
                m_pairs = m_pairs[np.sort(upos)]

        m_pairs = m_pairs.T

        return m_pairs

    def get_representation_brakets(self,
                                   freqs=None,
                                   freq_threshold=None,
                                   other=None,
                                   selection_rules=None,
                                   filter=None,
                                   return_filter=False
                                   ):

        return BasisStateSpace.get_representation_brakets(self,
                                                          other=other,
                                                          freqs=freqs,
                                                          freq_threshold=freq_threshold,
                                                          filter=filter,
                                                          return_filter=return_filter
                                                          )

    def to_single(self):
        """
        Condenses the multi state space down to
        a single BasisStateSpace

        :return:
        :rtype:
        """

        if self.mode == self.StateSpaceSpec.Indices:
            states = BasisStateSpace(
                self.basis,
                self.indices,
                mode=BasisStateSpace.StateSpaceSpec.Indices
            )
            if self.spaces[0]._excitations is not None:
                states.excitations = self.excitations
        else:
            states = BasisStateSpace(
                self.basis,
                self.excitations,
                mode=BasisStateSpace.StateSpaceSpec.Excitations
            )
        return states

    def take_states(self, states):
        """
        Takes the intersection of each held space and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """
        def take_inter(space, states=states):
            try:
                return space.take_states(states)
            except:
                raise ValueError(space, len(self.spaces[0]), self.spaces[0], self.spaces.shape)
        if self.spaces.ndim == 1:
            new_spaces = np.array([s.take_states(states) for s in self.spaces])
        else:
            new_spaces = np.apply_along_axis(take_inter, -1, self.spaces)
        return type(self)(new_spaces)
    def take_subspace(self, sel):
        """
        Takes the specified states, making sure each held space
        only contains states in `sel`
        :param sel:
        :type sel:
        :return:
        :rtype:
        """

        subsel = self.indices[sel,]
        return self.take_states(subsel)
    def take_subdimensions(self, inds):
        """
        Takes the subdimensions from each space
        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        def take(space, inds=inds):
            return space.take_subdimensions(inds)
        new_spaces = np.apply_along_axis(take, -1, self.spaces)
        return type(self)(new_spaces)

    def drop_states(self, states):
        """
        Take the difference of each held space and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """
        def take_diff(space, states=states):
            try:
                return space.drop_states(states)
            except:
                raise ValueError(space, len(self.spaces[0]), self.spaces[0], self.spaces.shape)
        if self.spaces.ndim == 1:
            new_spaces = np.array([s.take_states(states) for s in self.spaces])
        else:
            new_spaces = np.apply_along_axis(take_diff, -1, self.spaces)
        return type(self)(new_spaces)
    def drop_subspace(self, sel):
        """
        Takes the specified states, making sure each held space
        only contains states in `sel`
        :param sel:
        :type sel:
        :return:
        :rtype:
        """

        subsel = self.indices[sel,]
        return self.drop_subspace(subsel)
    def drop_subdimensions(self, inds):
        """
        Takes the subdimensions from each space
        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        def take(space, inds=inds):
            return space.drop_subdimensions(inds)
        new_spaces = np.apply_along_axis(take, -1, self.spaces)
        return type(self)(new_spaces)

    def __repr__(self):
        return "{}(nstates={}, shape={}, basis={})".format(
            type(self).__name__,
            len(self),
            self.spaces.shape,
            self.basis
        )

class SelectionRuleStateSpace(BasisMultiStateSpace):
    """
    A `BasisMultiStateSpace` subclass that is only built from applying selection rules to an initial space
    """
    def __init__(self, init_space, excitations, selection_rules=None):
        """
        :param init_space:
        :type init_space:
        :param excitations:
        :type excitations:
        :param selection_rules:
        :type selection_rules:
        """

        if isinstance(excitations, np.ndarray) and excitations.dtype == int:
            excitations = [BasisStateSpace(init_space.basis, x) for x in excitations]
        self._base_space = init_space
        self.sel_rules = selection_rules
        super().__init__(excitations)

    def to_state(self, serializer=None):
        return {
            'base_space':self._base_space,
            'spaces':self.spaces,
            'selection_rules':self.sel_rules
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(
            serializer.deserialize(data['base_space']),
            serializer.deserialize(data['spaces']),
            selection_rules=serializer.deserialize(data['selection_rules'])
            )

    def as_indices(self):
        """
        Pulls the full set indices out of all of the
        held spaces and returns them as a flat vector
        :return:
        :rtype:
        """
        sups = super().as_indices()
        base_inds = self._base_space.as_indices()
        if len(sups) == 0:
            return base_inds
        return np.concatenate([base_inds, sups])
    def as_excitations(self):
        """
        Pulls the full set excitations out of all of the
        held spaces and returns them as a flat vector
        :return:
        :rtype:
        """
        sups = super().as_excitations()
        if len(sups) == 0:
            return self._base_space.as_excitations()
        else:
            return np.concatenate([self._base_space.as_excitations(), sups], axis=0)

    @property
    def representative_space(self):
        return self._base_space

    def take_states(self, states):
        """
        Takes the intersection of each held space and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """

        def take_inter(space, states=states):
            try:
                return space.take_states(states)
            except:
                raise ValueError(space, len(self.spaces[0]), self.spaces[0], self.spaces.shape)

        if self.spaces.ndim == 1:
            new_spaces = np.array([s.take_states(states) for s in self.spaces])
        else:
            new_spaces = np.apply_along_axis(take_inter, -1, self.spaces)

        return type(self)(self._base_space, new_spaces)
    def take_subspace(self, states):
        """
        Takes the intersection of each held space and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """

        def take_inter(space, states=states):
            try:
                return space.take_subspace(states)
            except:
                raise ValueError(space, len(self.spaces[0]), self.spaces[0], self.spaces.shape)

        if self.spaces.ndim == 1:
            new_spaces = np.array([s.take_subspace(states) for s in self.spaces])
        else:
            new_spaces = np.apply_along_axis(take_inter, -1, self.spaces)

        return type(self)(self._base_space, new_spaces)
    def take_subdimensions(self, inds):
        """
        Takes the subdimensions from each space
        :param inds:
        :type inds:
        :return:
        :rtype:
        """

        def take(space, inds=inds):
            return space.take_subdimensions(inds)
        new_spaces = np.apply_along_axis(take, -1, self.spaces)
        return type(self)(self._base_space.take_subdimensions(inds), new_spaces)

    def drop_states(self, states):
        """
        Takes the intersection of each held space and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """

        def take_inter(space, states=states):
            try:
                return space.drop_states(states)
            except:
                raise ValueError(space, len(self.spaces[0]), self.spaces[0], self.spaces.shape)

        if self.spaces.ndim == 1:
            new_spaces = np.array([s.drop_states(states) for s in self.spaces])
        else:
            new_spaces = np.apply_along_axis(take_inter, -1, self.spaces)

        return type(self)(self._base_space, new_spaces)
    def drop_subspace(self, inds):
        """
        Takes the intersection of each held space and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """

        def take_inter(space, states=inds):
            try:
                return space.drop_subspace(states)
            except:
                raise ValueError(space, len(self.spaces[0]), self.spaces[0], self.spaces.shape)

        if self.spaces.ndim == 1:
            new_spaces = np.array([s.drop_subspace(inds) for s in self.spaces])
        else:
            new_spaces = np.apply_along_axis(take_inter, -1, self.spaces)

        return type(self)(self._base_space, new_spaces)
    def drop_subdimensions(self, inds):
        """
        Takes the subdimensions from each space
        :param inds:
        :type inds:
        :return:
        :rtype:
        """

        def take(space, inds=inds):
            return space.drop_subdimension(inds)
        new_spaces = np.apply_along_axis(take, -1, self.spaces)
        return type(self)(self._base_space.drop_subdimension(inds), new_spaces)

    def get_representation_indices(self,
                                   other=None,
                                   freqs=None,
                                   freq_threshold=None,
                                   selection_rules=None,
                                   filter=None,
                                   return_filter=False
                                   ):
        """
        This is where this pays dividends, as we know that only the init_space and the held excitations can couple
        which reduces the combinatoric work by a factor of like 2.
        :return:
        :rtype:
        """

        if other is not None:
            raise ValueError("haven't implemented getting indices to a set of 'other' states...")
        if selection_rules is not None:
            raise ValueError("applying selection_rules in the index generation isn't well defined...")
        if freq_threshold is not None:
            raise ValueError("Haven't implemented freq. threshold yet...")

        inds_base = self.representative_space.indices
        inds_l = []
        inds_r = []
        if filter is None:
            if other is not None:
                filter = other
            else:
                filter = None
        for i, s in zip(inds_base, self.flat):
            if isinstance(s, SelectionRuleStateSpace):
                j = s.representative_space.indices
                inds_l.append(np.full(len(j), i, dtype=int))
                inds_r.append(j)
                if other is not None:
                    others, filter = s.get_representation_indices(
                        freqs=freqs, freq_threshold=freq_threshold, other=other,
                        selection_rules=selection_rules,
                        filter=filter,
                        return_filter=True
                    )
                else:
                    others = s.get_representation_indices(
                        freqs=freqs, freq_threshold=freq_threshold, other=other,
                        selection_rules=selection_rules,
                        filter=filter,
                        return_filter=False
                    )
                inds_l.append(others[0])
                inds_r.append(others[1])
            else:
                j = s.unique_indices
                inds_l.append(np.full(len(j), i, dtype=int))
                inds_r.append(j)

        inds_l = np.concatenate(inds_l)
        inds_r = np.concatenate(inds_r)

        upairs = np.array([inds_l, inds_r])
        # _, upos = np.unique(pairs, axis=0, return_index=True)
        # upairs = pairs[np.sort(upos), ...]
        # raise Exception(len(upos), len(pairs), len(inds_l), len(upairs))

        if return_filter:
            return upairs, filter
        else:
            return upairs

    def filter_representation_inds(self, ind_pairs, q_changes):
        """
        Filters representation indices by the allowed #quantum changes.
        Not sure I'll even need this, if `get_representation_indices` is tight enough.

        :param ind_pairs:
        :type ind_pairs:
        :param q_changes:
        :type q_changes:
        :return:
        :rtype:
        """

        raise NotImplementedError("unoptimized code path")

        b1 = self.to_single().take_states(ind_pairs[0])
        b2 = self.to_single().take_states(ind_pairs[1])

        e1 = b1.excitations
        e2 = b2.excitations

        diffs = np.sum(np.abs(e2 - e1), axis=1)
        good_doggo = np.full((len(diffs),), False, dtype=bool)
        for q in q_changes:
            good_doggo = np.logical_or(good_doggo, diffs == q)

        new_stuff = np.array([d[good_doggo] for d in ind_pairs])
        return new_stuff

    @classmethod
    def _from_permutations(cls, space, permutations, filter_space, selection_rules):
        """
        Applies a set of permutations to an initial space and returns a new state space.

        :param space:
        :type space:
        :param permutations:
        :type permutations:
        :param filter_space:
        :type filter_space: BasisStateSpace
        :return:
        :rtype:
        """

        og = space.excitations.astype(cls.excitations_dtype)


        excitations = np.full((len(og),), None, dtype=object)

        use_sparse = isinstance(permutations, SparseArray)
        if use_sparse:
            raise NotImplementedError('needs fixing to work')
            # we drop down to the scipy wrappers until I can improve broadcasting in SparseArray
            og = sp.csc_matrix(og, dtype='int8')
            permutations = permutations.data

        if filter_space is not None:
            filter_exc = filter_space.unique_excitations.astype(cls.excitations_dtype)
            # two quick filters
            filter_min = np.min(filter_exc)
            filter_max = np.max(filter_exc)
            # and the sieve for nd stuff

            nrows, ncols = filter_exc.shape
            dtype = {'names': ['f{}'.format(i) for i in range(ncols)],
                     'formats': ncols * [filter_exc.dtype]}
            filter_inds = filter_exc.view(dtype)
        else:
            filter_min = 0
            filter_max = None
            filter_inds = None

        # and now we add these onto each of our states to build a new set of states
        for i, o in enumerate(og):
            if use_sparse:
                raise NotImplementedError('needs fixing to work')
                new_states = sp.vstack([o]*(permutations.shape[0])) + permutations
                mins = np.min(new_states, axis=1).toarray()
                well_behaved = np.where(mins >= filter_min)
                new_states = new_states[well_behaved]
                if filter_space is not None:
                    maxes = np.max(new_states, axis=1).toarray()
                    well_behaved = np.where(maxes <= filter_max)
                    new_states = new_states[well_behaved]
                _, pos = np.unique(new_states, axis=0, return_index=True)
                new_states = new_states[np.sort(pos)]
            else:
                new_states = o[np.newaxis, :] + permutations
                mins = np.min(new_states, axis=1)
                well_behaved = np.where(mins >= filter_min)
                new_states = new_states[well_behaved]
                if filter_space is not None:
                    maxes = np.max(new_states, axis=1)
                    well_behaved = np.where(maxes <= filter_max)
                    new_states = new_states[well_behaved]
                _, pos = np.unique(new_states, axis=0, return_index=True)
                new_states = new_states[np.sort(pos)]

            # finally, if we have a filter space, we apply it
            if filter_space is not None:
                as_inds = new_states.view(dtype) # convert to indices
                new_states = np.intersect1d(filter_inds, as_inds)
                new_states = new_states.view(filter_exc.dtype).reshape(-1, ncols)

                # dropped = np.setdiff1d(as_inds, filter_inds, assume_unique=True) # the indices that _aren't_ in the filter set
                # new_states = np.setdiff1d(as_inds, dropped, assume_unique=True) # the indices that _aren't_ _not_ in the filter set
                excitations[i] = BasisStateSpace(space.basis, new_states, mode=BasisStateSpace.StateSpaceSpec.Excitations)
            else:
                excitations[i] = BasisStateSpace(space.basis, new_states, mode=BasisStateSpace.StateSpaceSpec.Excitations)

        return cls(space, excitations, selection_rules)

    @classmethod
    def _generate_selection_rule_permutations(cls, space, selection_rules):
        """
        Turns a set of selection rule specs into permutations

        :param space:
        :type space:
        :param selection_rules:
        :type selection_rules:
        :return:
        :rtype:
        """

        nmodes = space.ndim
        # group selection rules by how many modes they touch
        sel_rule_groups = {}
        for s in selection_rules:
            k = len(s)
            if k in sel_rule_groups:
                sel_rule_groups[k].append(s)
            else:
                sel_rule_groups[k] = [s]

        # now generate unique permutations based on these
        perms = []
        for n,v in sel_rule_groups.items():
            if n <= nmodes:
                padding = [0] * (nmodes - n)
                for p in v:
                    p = np.array(list(p) + padding, dtype=cls.excitations_dtype)
                    perms.append(UniquePermutations(p).permutations())
        permutations = np.concatenate(perms)

        return permutations

    @classmethod
    def _apply_rules_recursive(cls, space, permutations, filter_space, selection_rules, iterations=1):
        new = cls._from_permutations(space, permutations, filter_space, selection_rules)
        if iterations > 1:
            for n, s in enumerate(new):
                new[n] = cls._apply_rules_recursive(s,  permutations, filter_space, selection_rules, iterations=iterations-1)
        return new

    @classmethod
    def from_rules(cls, space, selection_rules, filter_space=None, iterations=1, method='new'):
        """
        :param space: initial space to which to apply the transformations
        :type space: BasisStateSpace | BasisMultiStateSpace
        :param selection_rules: different possible transformations
        :type selection_rules: Iterable[Iterable[int]]
        :param iterations: number of times to apply the transformations
        :type iterations: int
        :param filter_space: a space within which all generated `BasisStateSpace` objects must be contained
        :type filter_space: BasisStateSpace | None
        :return:
        :rtype: SelectionRuleStateSpace
        """

        if iterations > 1:
            raise NotImplementedError(
                "things have changed and higher iterations aren't currently supported but could be supported in the future by being smart with the selection rules"
            )

        if method == 'legacy':
            if filter_space is None:
                permutations = cls._generate_selection_rule_permutations(space, selection_rules)
                new = cls._apply_rules_recursive(space, permutations, filter_space, selection_rules,
                                                 iterations=iterations)
            else:
                # TODO: add a check to see if it's faster to do a quadratic-ish time filter
                #       or a generate and intersect
                #       NOTE: this results in a different ordering than the other approach would...
                #             but that shouldn't be an explicit issue?

                permutations = cls._generate_selection_rule_permutations(space, selection_rules)
                new = cls._find_space_intersections(space, filter_space, permutations, selection_rules)
                # raise Exception(new.excitations)
            return new
        else:
            exc = space.excitations
            symmetric_group_inds = hasattr(space.basis.indexer, 'symmetric_group')
            if symmetric_group_inds:
                symm_grp = space.basis.indexer.symmetric_group #type: SymmetricGroupGenerator
            else:
                symm_grp = SymmetricGroupGenerator(exc.shape[-1])
            if filter_space is None:
                new_exc, new_inds = symm_grp.take_permutation_rule_direct_sum(exc, selection_rules,
                                                                              return_indices=True, split_results=True)
            else:
                if isinstance(filter_space, BasisStateSpace):
                    filter_space = (filter_space.excitations, filter_space.indices)

                new_exc, new_inds, filter = symm_grp.take_permutation_rule_direct_sum(exc, selection_rules, filter_perms=filter_space,
                                                                              return_filter=True,
                                                                              return_indices=True, split_results=True)


                # print("?"*100)
                # print(len(new_inds), len(new_exc))
                # assert len(new_inds) == len(new_exc)
                # assert [len(x) for x in new_inds] == [len(x) for x in new_exc]
                # assert symm_grp.to_indices(np.concatenate(new_exc, axis=0)).tolist() == np.concatenate(new_inds).tolist()
                # print(
                #     [(i, j, exc[i]) for i,es in enumerate(new_exc) for j,e in enumerate(es) if len(es) > 0 and sum(e) == 5 and max(e) == 1 ],
                #     symm_grp.to_indices(np.concatenate(new_exc, axis=0)).tolist() == np.concatenate(new_inds).tolist()
                # )
                # # print(exc[:20], selection_rules, (filter_space[0][:10], filter_space[2][:10]))
                # print("="*100)

            new = []
            for e,i in zip(new_exc, new_inds):
                new_space = BasisStateSpace(space.basis, e, mode=BasisStateSpace.StateSpaceSpec.Excitations)
                new_space.indices = i
                new.append(new_space)
            new = np.array(new, dtype=object)

            if filter_space is None:
                return cls(space, new, selection_rules)
            else:
                return cls(space, new, selection_rules), filter

    @classmethod
    def _find_space_intersections(cls, space, filter_space, perms, selection_rules):
        """
        :param space: set of initial states to generate connections off of
        :type space: BasisStateSpace
        :param filter_space: set of connected states to test against
        :type filter_space: BasisStateSpace
        :param rules: set of possible selection changes that would work
        :type rules: np.ndarry
        :return:
        :rtype:
        """

        # we filter on the filter_space by first computing the
        exc_1 = space.excitations
        filter_space = filter_space.take_unique()
        exc_2 = filter_space.excitations

        nrows, ncols = perms.shape
        dtype = {'names': ['f{}'.format(i) for i in range(ncols)],
                 'formats': ncols * [perms.dtype]}
        rules_view = perms.view(dtype) # convert to indices


        # calculate _all_ exc_2 - exc_1 diffs to test agains
        # import time
        # start = time.time()

        nq_changes = np.sum(perms, axis=1)
        excitations = np.full((len(exc_1),), None, dtype=object)
        for i, state in enumerate(exc_1):
            # start=time.time()
            diffs = exc_2 - state[np.newaxis, :]
            quanta_diffs = np.sum(np.abs(diffs), axis=1)
            mask = np.full(len(quanta_diffs), False)
            for nq in nq_changes:
                where_q = np.where(quanta_diffs == nq)
                mask[where_q] = True

            proper_inds = np.where(mask)
            if len(proper_inds) == 0:
                excitations[i] = BasisStateSpace(filter_space.basis, [], mode=filter_space.StateSpaceSpec.Indices)
                continue

            proper_inds = proper_inds[0]
            test_diffs = diffs[proper_inds]
            # raise Exception(dtype, test_diffs.dtype)
            _, inter_inds, _ = np.intersect1d(test_diffs.view(dtype), rules_view, return_indices=True)
            clean_inds = proper_inds[inter_inds]
            if len(clean_inds) > 0:
                # clean_inds = np.sort(clean_inds)
                excitations[i] = filter_space.take_subspace(clean_inds)
            else:
                excitations[i] = BasisStateSpace(filter_space.basis, [], mode=filter_space.StateSpaceSpec.Indices)

        return cls(space, excitations, selection_rules)

    @staticmethod
    def _pick(indices, to_search, indexer):
        if not isinstance(to_search, np.ndarray):
            if isinstance(to_search, AbstractStateSpace) or hasattr(to_search, 'indices'):
                to_search = to_search.indices
            else:
                to_search = np.array(to_search)
        vals = np.searchsorted(indices, to_search, sorter=indexer)
        if isinstance(vals, (np.integer, int)):
            vals = np.array([vals])
        # we have the ordering according to the _sorted_ version of `indices`
        # so now we need to invert that back to the unsorted version
        if len(indexer) > 0:
            big_vals = vals == len(indexer)
            vals[big_vals] = -1
            vals = indexer[vals]
            # now because of how searchsorted works, we need to check if the found values
            # truly agree with what we asked for
            bad_vals = indices[vals] != to_search
            if vals.shape == ():
                if bad_vals:
                    vals = -1
            else:
                vals[bad_vals] = -1
        else:
            vals = np.full_like(vals, -1)

        return vals

    def union(self, other, handle_subspaces=True):
        """
        Returns a merged version of self and other, adding
        any states in other to self and merging where they intersect

        :param other:
        :type other: SelectionRuleStateSpace
        :return:
        :rtype:
        """

        if not isinstance(other, SelectionRuleStateSpace):
            raise TypeError("union with {} only defined over subclasses of {}".format(
                type(self).__name__,
                SelectionRuleStateSpace.__name__
            ))

        if self.basis is not other.basis:
            raise ValueError("can't merge state spaces over different bases ({} and {})".format(
                self.basis,
                other.basis
            ))

        excitation_mode = (
                self.representative_space.has_excitations
                and other.representative_space.has_excitations
        )
        if excitation_mode: # special case I guess?
            self_exc = self.representative_space.excitations
            other_exc = other.representative_space.excitations

            nrows, ncols = self_exc.shape
            dtype = {'names': ['f{}'.format(i) for i in range(ncols)],
                     'formats': ncols * [self_exc.dtype]}

            self_inds = self_exc.view(dtype)
            other_inds = other_exc.view(dtype)
            other_exclusions = np.setdiff1d(other_inds, self_inds)

            _, where_inds, _ = np.intersect1d(other_inds, other_exclusions, return_indices=True)
            where_inds = np.sort(where_inds)

        else:
            self_inds = self.representative_space.indices
            other_inds = other.representative_space.indices
            other_exclusions = np.setdiff1d(other_inds, self_inds)

            _, where_inds, _ = np.intersect1d(other_inds, other_exclusions, return_indices=True)
            where_inds = np.sort(where_inds)

        new_rep = self.representative_space.union(other.representative_space)

        new_spaces = np.concatenate(
            [
                self.spaces,
                other.spaces[where_inds,]
            ],
            axis=0
        )

        if handle_subspaces:
            if excitation_mode:
                new_inds = new_rep.excitations.view(dtype)
            else:
                new_inds = new_rep.indices
            inter_inds = np.intersect1d(self_inds, other_inds)
            # try:
            _, where_new, _ = np.intersect1d(new_inds, inter_inds, return_indices=True)
            # except:
            #     raise Exception(new_inds, inter_inds)
            _, where_old, _ = np.intersect1d(other_inds, inter_inds, return_indices=True)
            # try:
            #     deb = self.debug_flag
            # except AttributeError:
            #     pass
            # else:
            #     if deb:
            #         raise Exception(
            #             inter,
            #             new_inds[where_new],
            #             other_inds[where_old]
            #         )
            for i_new, i_old in zip(where_new, where_old):
                new_spaces[i_new] = new_spaces[i_new].union(other[i_old])

        return type(self)(new_rep, new_spaces)

    def intersection(self, other, handle_subspaces=True, use_indices=False):
        """
        Returns an intersected self and other

        :param other:
        :type other: SelectionRuleStateSpace
        :return:
        :rtype:
        """

        if not isinstance(other, SelectionRuleStateSpace):
            raise TypeError("intersection with {} only defined over subclasses of {}".format(
                type(self).__name__,
                SelectionRuleStateSpace.__name__
            ))

        if self.basis is not other.basis:
            raise ValueError("can't merge state spaces over different bases ({} and {})".format(
                self.basis,
                other.basis
            ))

        # create intersection based on indices and then
        # make use of this subselection to resample the basis


        if not use_indices and (
                self.representative_space.has_excitations
                and other.representative_space.has_excitations
        ): # special case I guess?
            self_exc = self.representative_space.excitations
            other_exc = other.representative_space.excitations

            nrows, ncols = self_exc.shape
            dtype = {'names': ['f{}'.format(i) for i in range(ncols)],
                     'formats': ncols * [self_exc.dtype]}

            self_inds = self_exc.view(dtype)
            other_inds = other_exc.view(dtype)
            _, where_inds, other_where = np.intersect1d(self_inds, other_inds, return_indices=True)

            where_inds = np.sort(where_inds)

        else:
            self_inds = self.representative_space.indices
            other_inds = other.representative_space.indices

            _, where_inds, other_where = np.intersect1d(self_inds, other_inds, return_indices=True)
            where_inds = np.sort(where_inds)

        new_spaces = self.spaces[where_inds,]
        if handle_subspaces:
            for n,i in enumerate(other_where):
                new_spaces[n] = new_spaces[n].intersection(other[n])

        new_rep = self.representative_space.take_subspace(where_inds)#(other.representative_space)

        return type(self)(new_rep, new_spaces)
    def difference(self, other, handle_subspaces=True):
        """
        Returns an diff'ed self and other.
        We get fundamentally different behaviour for `handle_subspaces` than without it.
        If we have it _on_ then differences are computed for each states in the intersection of
          the primary (key) states.
        If we have it off then the difference in the key states is computed and nothing more is
        done.

        :param other:
        :type other: SelectionRuleStateSpace
        :return:
        :rtype: SelectionRuleStateSpace
        """

        if not isinstance(other, SelectionRuleStateSpace):
            raise TypeError("difference with {} only defined over subclasses of {} (not {})".format(
                type(self).__name__,
                SelectionRuleStateSpace.__name__,
                type(other).__name__
            ))

        if self.basis is not other.basis:
            raise ValueError("can't take a difference of state spaces over different bases ({} and {})".format(
                self.basis,
                other.basis
            ))

        self_inds = self.representative_space.indices
        other_inds = other.representative_space.indices

        if handle_subspaces:
            new_rep = self.representative_space.intersection(other.representative_space)
            intersected_inds = new_rep.indices
            indexer = np.argsort(self_inds)
            where_inds = self._pick(self_inds, intersected_inds, indexer)
            new_spaces = self.spaces[where_inds,]
            sub_indexer = np.argsort(other_inds)
            where_inds = self._pick(other_inds, intersected_inds, sub_indexer)
            for n, i in enumerate(where_inds):
                new_spaces[n] = new_spaces[n].difference(other[n])
        else:
            new_rep = self.representative_space.difference(other.representative_space)
            diff_inds = new_rep.indices
            indexer = np.argsort(self_inds)
            where_inds = self._pick(self_inds, diff_inds, indexer)
            new_spaces = self.spaces[where_inds,]

        return type(self)(new_rep, new_spaces)

    def __getitem__(self, item):
        it = self.spaces[item]
        if isinstance(it, np.ndarray):
            # multispace
            rinds = self.representative_space.indices
            init = self.representative_space.take_states(rinds[item])
            it = type(self)(init, it)
        return it
    def __setitem__(self, item, vals):
        self.spaces[item] = vals
        self._indices = None
        self._excitations = None
        self._indexer = None
        self._uinds = None

    def __repr__(self):
        return "{}(ogstates={}, nstates={}, basis={})".format(
            type(self).__name__,
            len(self.spaces),
            len(self),
            self.basis
        )

class BraKetSpace:
    """
    Represents a set of pairs of states that can be fed into a `Representation` or `Operator`
    to efficiently tell it what terms it need to calculate.
    This basically just implements a bunch of stuff for generating a Graph defining
    the connections between states.
    """
    def __init__(self,
                 bra_space,
                 ket_space
                 ):
        """
        :param bra_space:
        :type bra_space: BasisStateSpace
        :param ket_space:
        :type ket_space: BasisStateSpace
        """
        self.bras = bra_space
        self.kets = ket_space
        self.ndim = self.bras.ndim
        self._orthogs = None
        self._preind_trie = None
        if len(bra_space) != len(ket_space) or (bra_space.ndim != ket_space.ndim):
            raise ValueError("Bras {} and kets {} have different dimension".format(bra_space, ket_space))
        self.state_pairs = (
            self.bras.excitations.T,
            self.kets.excitations.T
        )


    @classmethod
    def from_indices(cls, inds, basis=None, quanta=None):
        if basis is None:
            if quanta is None:
                raise ValueError("{}.{}: either basis of number of quanta (assumes harmonic basis) is required".format(
                    cls.__name__,
                    'from_indices'
                ))
            from .HarmonicOscillator import HarmonicOscillatorProductBasis
            basis = HarmonicOscillatorProductBasis(quanta)

        # we need to coerce the indices into a set of bra states and ket states...
        inds = np.asarray(inds, dtype=int) # easier for us to work with @ the cost of a bit of perf.

        if inds.shape[1] == 2: # we assume we have a (NDim X (Bra, Ket), X States) array
            bra_states = inds[:, 0, ...]
            ket_states = inds[:, 1, ...]
        elif inds.shape[0] == 2:
            bra_states = inds[0]
            ket_states = inds[1]
        else:
            raise ValueError("don't know what to do with indices array of shape {}".format(inds.shape))

        return BraKetSpace(
            BasisStateSpace(basis, bra_states, mode=BasisStateSpace.StateSpaceSpec.Indices),
            BasisStateSpace(basis, ket_states, mode=BasisStateSpace.StateSpaceSpec.Indices)
        )

    def __len__(self):
        return len(self.state_pairs[0][0])

    def __repr__(self):
        return "{}(nstates={})".format(type(self).__name__, len(self))

    def load_non_orthog(self,
                        use_aggressive_caching=None,
                        use_preindex_trie=None,
                        preindex_trie_depth=None
                        ):
        if self._orthogs is None:

            exc_l, exc_r = self.state_pairs
            exc_l = np.asanyarray(exc_l)#, dtype='int8')
            exc_r = np.asanyarray(exc_r)#, dtype='int8')
            womp = np.equal(exc_l, exc_r)

            if use_preindex_trie:
                trie = self.OrthogoIndexerTrie(womp, max_depth=preindex_trie_depth)
            else:
                trie = None

            # caching version tries very hard to be clever, non-caching version
            # only caches up to the specified trie depth
            if use_aggressive_caching:
                self._orthogs = self.CachingOrthogonalIndexCalculator(womp, trie)
            else:
                self._orthogs = self.OrthogonalIndexCalculator(womp, trie)

    class OrthogoIndexerTrie:
        """
        A trie that makes it faster to determine which states are
        non-orthogonal excluding some subset of indices
        """
        default_max_depth=8 # more depth means faster evals but more memory; worth getting a trade-off here
        def __init__(self, base_orthogs, max_depth=None):
            self.orthogs = base_orthogs
            if max_depth is None:
                max_depth = self.default_max_depth
            self.max_depth = max_depth
            self.trie = {}#{'idx':np.arange(len(self.orthogs[0]))}

        def get_idx_terms(self, idx):
            """
            idx is a assumed sorted and since so many
            of the request idx will start with the same
            first few indices, we'll get a significant benefit at the end
            of the day
            :param idx:
            :type idx:
            :return:
            :rtype:
            """
            trie = self.trie
            tests = self.orthogs
            if len(idx) > 1:
                # we explicilty calculate paths of length 1 because of the memory
                # required to store them
                for n,i in enumerate(idx[:self.max_depth]):
                    if i not in trie:
                        if n == 0:
                            trie[i] = {}
                        elif n == 1:
                            cur_idx = np.where(tests[idx[0]])
                            if len(cur_idx) > 0:
                                cur_idx = cur_idx[0]
                            trie[i] = {'idx': cur_idx[tests[i, cur_idx]]}
                        else:
                            cur_idx = trie['idx']
                            trie[i] = {'idx': cur_idx[tests[i, cur_idx]]}
                    trie = trie[i]

                return trie['idx'], idx[self.max_depth:]
            elif len(idx) == 1:
                vals = np.where(tests[idx[0]])
                if len(vals) > 0:
                    vals = vals[0]
                return vals, ()
            else:
                return np.array([], dtype=BasisStateSpace.excitations_dtype), ()
                # idx_1 = np.arange(len(self.orthogs[0]))
                # for i in idx[:self.max_depth]:
                #     if i not in trie:
                #         cur_idx = trie['idx']
                #         trie[i] = {'idx':cur_idx[tests[i, cur_idx]]}
                #     trie = trie[i]

    class CachingOrthogonalIndexCalculator:
        """
        Provides a way to use Boolean algebra to construct new values from old values on
        the cheap
        """

        def __init__(self, eq_tests, trie):
            self.tests = eq_tests
            self.pretest = eq_tests.flatten().all()
            self.cache = {}#MaxSizeCache()
            self.trie = trie

        def get_idx_comp(self, item):
            """
            Pulls a set of indices either by direct computation
            or by using some Boolean algebra to refine the number of
            cases under which we actually need to do a full computation.
            Allows us to reuse old info.
            :param item:
            :type item:
            :return:
            :rtype:
            """
            if self.pretest:
                return np.arange(len(self.tests[0]))
            # skip doing any work if we're in this special case
            rest_inds = np.setdiff1d(np.arange(len(self.tests)), item)
            if len(rest_inds) == 0:
                return np.arange(len(self.tests[0]))



            # we first determine which sets of indices we've already
            # computed so that we can reuse that information
            prev_cache = None
            cache = self.cache
            for n, i in enumerate(item):
                if i in cache:
                    cache = cache[i]
                else:
                    # we check if we're only one element away from the end
                    # i.e. if we can use the current cache info to efficiently calculate the
                    # next cache bit
                    if n == len(item) - 1 and 'vals' in cache:
                        cache[i] = {
                            'vals':self._get_nonorthog_from_prev(cache, item)
                        }
                        cache = cache[i]
                    else:
                        # otherwise we just fall through
                        cache = self.cache
                        for j in item:
                            if j not in cache:
                                cache[j] = {}
                            cache = cache[j]
                        if self.trie is None:
                            cache['vals'] = self._pull_nonorthog(item)
                        else:
                            cache['vals'] = self._pull_nonorthog_trie(item, self.trie)
                    break
            else:
                # this means we made it all the way through,
                # so now we either check to see if we've already computed the indices
                if 'vals' not in cache:
                    # we see if we can compute this stuff from a later term
                    for next_key in cache.keys():
                        if 'vals' in cache[next_key]:
                            break
                    else:
                        next_key = None
                    if next_key is not None:
                        cache['vals'] = self._get_nonorthog_from_next(cache, next_key)
                    elif prev_cache is not None and 'vals' in prev_cache:
                        # we use prior info to be somewhat more efficient
                        cache['vals'] = self._get_nonorthog_from_prev(prev_cache, item)
                    else:
                        # directly compute as a fallback
                        if self.trie is None:
                            cache['vals'] = self._pull_nonorthog(item)
                        else:
                            cache['vals'] = self._pull_nonorthog_trie(item, self.trie)

            return cache['vals']

        def _get_nonorthog_from_next(self, cache, next_key):
            """
            Adds the restriction from next_key to the values already in cache
            :param cache:
            :type cache:
            :param next_key:
            :type next_key:
            :return:
            :rtype:
            """
            cur_pos = cache[next_key]['vals']
            return cur_pos[self.tests[next_key][cur_pos]]

        def _get_nonorthog_from_prev(self, cache, inds):
            """
            Calculates the terms associated with removing the restriction
            on term i from the vals in cache

            :param cache:
            :type cache:
            :param next:
            :type next:
            :return:
            :rtype:
            """
            cur_inds = cache['vals']
            next_tests = self.tests[inds[-1]]

            if self.trie is not None:
                # we use the trie to do an initial index filtering which allows
                # us to do (potentially) even less work than before
                unused = np.setdiff1d(np.arange(len(self.tests)), inds)
                init_inds, rest = self.trie.get_idx_terms(unused)
                init_inds = np.setdiff1d(init_inds, cur_inds)
                init_inds = init_inds[np.logical_not(next_tests[init_inds])]

                for e in rest:
                    init_inds = init_inds[self.tests[e, init_inds]]

                recalcs = init_inds

            else:
                # we can reuse the entirety of cur_inds, since the next case is just a relaxation on those
                # the new positions we need to calculate are just the intersection of the complement of cur_inds and where
                # next_tests fails...which is just the symmetric difference of where next_tests fails and cur_inds
                recalc_pos = np.setdiff1d(np.where(np.logical_not(next_tests))[0], cur_inds)

                if self.trie is None:
                    recalcs = self._pull_nonorthog(inds, subinds=recalc_pos)
                else:
                    recalcs = self._pull_nonorthog_trie(inds, self.trie, subinds=recalc_pos)

            return np.sort(np.concatenate([cur_inds, recalcs]), kind='mergesort') # this sort might kill any benefit...

        def _pull_nonorthog(self, inds, subinds=None):
            """
            Directly computes orthogonality relations
            :param inds:
            :type inds:
            :return:
            :rtype:
            """
            orthos = self.tests
            unused = np.delete(np.arange(len(orthos)), inds)
            if subinds is None:
                cur_inds = np.where(self.tests[unused[0]])[0]
                unused = unused[1:]
                for e in unused:
                    cur_inds = cur_inds[self.tests[e, cur_inds]]
                return cur_inds
            else:
                cur_inds = subinds
                for e in unused:
                    cur_inds = cur_inds[self.tests[e, cur_inds]]
                return cur_inds

        def _pull_nonorthog_trie(self, inds, preindex_trie, subinds=None):
            orthos = self.tests
            unused = np.delete(np.arange(len(orthos)), inds)
            init_inds, rest = preindex_trie.get_idx_terms(unused)
            if subinds is not None:
                init_inds = np.intersect1d(init_inds, subinds)
            for e in rest:
                init_inds = init_inds[self.tests[e, init_inds]]
            return init_inds

    # def load_preindex_trie(self, tests, max_depth=None):
    #     """
    #     Loads the preindex trie we'll use to speed up non-orthogonality calcs.
    #     If only a single non-orthog is needed, this probably isn't the way to go...
    #
    #     :param max_comp:
    #     :type max_comp:
    #     :param max_depth:
    #     :type max_depth:
    #     :return:
    #     :rtype:
    #     """
    #     if self._preind_trie is None:
    #         self._preind_trie = self.OrthogoIndexerTrie(tests, max_depth=max_depth)
    #     return self._preind_trie

    class OrthogonalIndexSparseCalculator:
        def __init__(self, orthogs):
            self.tests = sp.csc_matrix(orthogs.T)
            # if len(np.where(orthogs == 0)[0]) > 0:
            #     raise Exception(self.tests.nnz, np.prod(self.tests.shape))
        def get_idx_comp(self, inds):
            nterms = self.tests.shape[1]
            ninds = len(inds)
            sampling_vector = np.ones(nterms)
            sampling_vector[inds,] = 0

            out = self.tests.dot(sampling_vector)
            targ_val = nterms-ninds
            return np.where(out==targ_val)[0]

    def clear_cache(self):
        self._orthogs = None

    class OrthogonalIndexCalculator:
        def __init__(self, tests, trie):
            self.tests = tests # a bunch of t/f statements
            self.pretest = tests.flatten().all()
            self.trie = trie # cached prefiltering
        def get_idx_comp(self, item):
            """
            :param item:
            :type item:
            :return:
            :rtype:
            """

            if self.pretest:
                return np.arange(len(self.tests[0]))

            rest_inds = np.setdiff1d(np.arange(len(self.tests)), item)
            if len(rest_inds) == 0:
                return np.arange(len(self.tests[0]))

            if self.trie is None:
                cur_inds = np.where(self.tests[rest_inds[0]])[0]
                rest_inds = rest_inds[1:]
            else:
                cur_inds, rest_inds = self.trie.get_idx_terms(rest_inds)

            for e in rest_inds:
                cur_inds = cur_inds[self.tests[e, cur_inds]]
            return cur_inds

    aggressive_caching_enabled=True
    preindex_trie_enabled=True
    def get_non_orthog(self,
                       inds,
                       assume_unique=False,
                       use_aggressive_caching=None,
                       use_preindex_trie=None,
                       preindex_trie_depth=None
                       ):
        """
        Returns whether the states are non-orthogonal under the set of indices.

        :param inds:
        :type inds:
        :return:
        :rtype:
        """

        if use_aggressive_caching is None:
            use_aggressive_caching = self.aggressive_caching_enabled
        if use_preindex_trie is None:
            use_preindex_trie = self.preindex_trie_enabled

        if not assume_unique:
            inds = np.unique(inds)

        self.load_non_orthog(
            use_aggressive_caching=use_aggressive_caching,
            use_preindex_trie=use_preindex_trie,
            preindex_trie_depth=preindex_trie_depth
        )
        inds = tuple(np.sort(inds))

        return self._orthogs.get_idx_comp(inds)

    def get_sel_rules_from1d(self, inds, rules):
        from collections import OrderedDict

        uinds = OrderedDict((k, None) for k in inds)
        uinds = np.array(list(uinds.keys()))
        # from McUtils.Scaffolding import Logger
        # logger = Logger.lookup("debug")
        # logger.log_print("uinds: {u}", u=uinds)
        mm = {k: i for i, k in enumerate(uinds)}
        subdim = len(uinds)
        sel_bits = [None] * subdim
        for s, i in zip(rules, inds):
            n = mm[i]  # makes sure that we fill in in the same order as uinds
            if sel_bits[n] is None:
                sel_bits[n] = s
            else:
                new_rules = np.add.outer(s, sel_bits[n])
                _, pos = np.unique(new_rules, axis=0, return_index=True)
                sel_bits[n] = new_rules[np.sort(pos)].flatten()  # make the next set of rules

        return sel_bits

    def get_sel_rule_filter(self, rules):

        import functools as fp

        bra, ket = self.state_pairs

        # we have a set of rules for every dimension in states...
        # so we define a function that will tell us where any of the possible selection rules are satisfied
        # and then we apply that in tandem to the quanta changes between bra and ket and the rules
        rules = [np.array(r) for r in rules]

        # from McUtils.Scaffolding import Logger
        # logger = Logger.lookup("debug")
        # logger.log_print("oh! {r}", r=rules)
        apply_rules = lambda diff, rule: fp.reduce(
            lambda d, i: np.logical_or(d, diff == i),
            rule[1:],
            diff == rule[0]
        )
        sels = [
            apply_rules(k - b, r) for b, k, r in zip(bra, ket, rules)
        ]
        # then we figure out which states by satisfied all the rules
        return fp.reduce(np.logical_and, sels[1:], sels[0])

    def take_subspace(self, sel):
        return type(self)(
            self.bras.take_subspace(sel),
            self.kets.take_subspace(sel)
        )

    def take_subdimensions(self, inds):
        return type(self)(
            self.bras.take_subdimensions(inds),
            self.kets.take_subdimensions(inds)
        )

    def apply_non_orthogonality(self,
                                inds,
                                use_aggressive_caching=None,
                                use_preindex_trie=None,
                                preindex_trie_depth=None,
                                assume_unique=False
                                ):
        """
        Takes the bra-ket pairs that are non-orthogonal under the indices `inds`

        :param inds:
        :type inds:
        :param assume_unique:
        :type assume_unique:
        :return:
        :rtype:
        """
        if use_aggressive_caching is None:
            use_aggressive_caching = self.aggressive_caching_enabled
        if use_preindex_trie is None:
            use_preindex_trie = self.preindex_trie_enabled
        non_orthog = self.get_non_orthog(inds,
                                         use_aggressive_caching=use_aggressive_caching,
                                         assume_unique=assume_unique,
                                         use_preindex_trie=use_preindex_trie,
                                         preindex_trie_depth=preindex_trie_depth
                                         )
        return self.take_subspace(non_orthog), non_orthog

    def apply_sel_rules(self, rules):
        """
        Applies selections rules
        :param rules:
        :type rules:
        :return:
        :rtype:
        """

        all_sels = self.get_sel_rule_filter(rules)
        return self.take_subspace(all_sels), all_sels

    def adjacency_matrix(self, total_space=None):
        """
        Generates the (sparse) unweighted adjacency matrix for the bras & kets
        :return:
        :rtype:
        """

        base_inds = np.concatenate([self.bras.unique_indices, self.kets.unique_indices])
        _, upos = np.unique(base_inds, return_index=True)
        inds = base_inds[np.sort(upos)]
        if total_space is None:
            total_space = BasisStateSpace(self.bras.basis, inds, mode=BasisStateSpace.StateSpaceSpec.Indices)
        bra_inds = total_space.find(self.bras.indices)
        ket_inds = total_space.find(self.kets.indices)

        # raise Exception(len(ket_inds), len(bra_inds))

        utri = np.array([bra_inds, ket_inds]).T
        ltri = np.array([ket_inds, bra_inds]).T
        indies = np.unique(np.concatenate([utri, ltri]), axis=0)

        return SparseArray(
            (
                np.ones(len(indies)),
                indies.T
            ),
            shape=(len(total_space), len(total_space))
        )

    def split(self, chunksize):
        """
        splits the brakets into blocks of at max chunksize
        :param chunksize:
        :type chunksize: int
        :return:
        :rtype: Iterable[BraKetSpace]
        """

        bra_splits = self.bras.split(chunksize)
        ket_splits = self.kets.split(chunksize)

        # raise Exception(ket_splits)

        return [type(self)(b, k) for b,k in zip(bra_splits, ket_splits)]

# # We might need this in the future, but it feels premature to set it up just yet...
# class SharedBasisStateSpace(BasisStateSpace):
#     """
#     A mutable basis state space that can register "listeners"
#     """

class StateSpaceMatrix:
    """
    A `SparseArray` that holds onto a `BasisStateSpace` that keeps track of the
    total set of states involved.
    By default is assumed real-symmetric. This can be relaxed in the future.

    TODO: The idea is good, but calculating what is "in" the array and what is "out"
            every single time this is applied could be slow...
          We'll need to test to see how slow
    """

    def __init__(self, initial_basis, initial_vals=None, column_space=None, symmetric=True):
        """

        :param initial_basis:
        :type initial_basis: BasisStateSpace | RepresentationBasis
        """

        if isinstance(initial_basis, BraKetSpace):
            if initial_vals is not None:
                bras = initial_basis.bras
                kets = initial_basis.kets
                initial_basis = bras.union(kets).take_unique(sort=True)
                row_inds = initial_basis.find(bras)
                col_inds = initial_basis.find(kets)
                if symmetric:
                    # this is all done over the total state space rather than the column
                    # space to start so that the symmetry can be made manifest
                    # but then it needs to be reduced over the column space at the end
                    # a basic dupe of this appears later
                    upper_triangle = np.where(row_inds <= col_inds)[0]
                    initial_vals = initial_vals[upper_triangle]
                    col_inds = col_inds[upper_triangle]
                    row_inds = row_inds[upper_triangle]
                    flippers = np.where(col_inds != row_inds)[0]
                    initial_vals = np.concatenate([initial_vals, initial_vals[flippers]])
                    new_row_inds = np.concatenate([row_inds, col_inds[flippers]])
                    new_col_inds = np.concatenate([col_inds, row_inds[flippers]])
                    row_inds = new_row_inds
                    col_inds = new_col_inds
                if column_space is not None:
                    column_space = column_space.take_unique(sort=True)
                    col_inds = column_space.find(initial_basis.take_subspace(col_inds))
                initial_vals = (initial_vals, (row_inds, col_inds))
            else:
                initial_basis = initial_basis.bras.union(initial_basis.kets).take_unique(sort=True)
                if column_space is not None:
                    column_space = column_space.take_unique(sort=True)
        else:
            if initial_vals is not None and not isinstance(initial_vals, SparseArray):
                raise ValueError("can only initialize a {} with values if a {} is passed too".format(
                    type(self).__name__,
                    BraKetSpace.__name__
                ))
            if not isinstance(initial_basis, BasisStateSpace):
                initial_basis = BasisStateSpace(initial_basis, [], mode=BasisStateSpace.StateSpaceSpec.Indices)

        self.symmetric = symmetric
        self.states = initial_basis
        self.column_space = column_space
        if column_space is not None and len(column_space) > len(initial_basis):
            raise ValueError("can't have column space {} be larger than total space {}".format(column_space, initial_basis))
        self._brakets = None
        if initial_vals is None:
            self.array = SparseArray.empty(self.shape)
        elif isinstance(initial_vals, SparseArray):
            self.array = initial_vals
        else:
            self.array = SparseArray.from_data(initial_vals, shape=self.shape)

    @property
    def shape(self):
        if self.column_space is None:
            return (len(self.states), len(self.states))
        else:
            return (len(self.states), len(self.column_space))
    @property
    def basis(self):
        """
        Returns the basis for the matrix rep

        :return:
        :rtype:
        """
        return self.states.basis

    @property
    def brakets(self):
        """
        Returns the BraKetSpace for the held indices
        :return:
        :rtype:
        """
        if self._brakets is None:
            self._brakets = self._get_brakets()
        return self._brakets

    @classmethod
    def identity_from_space(cls, space, column_space=None):
        """
        Returns a StateSpaceMatrix where the diagonal is filled with 1s
        :param space:
        :type space:
        :param column_space:
        :type column_space:
        :return:
        :rtype:
        """
        if column_space is not None:
            vals = np.ones(len(column_space))
            brakets = BraKetSpace(space, column_space)
        else:
            vals = np.ones(len(space))
            brakets = BraKetSpace(space, space)
        return cls(brakets, initial_vals=vals, column_space=column_space)

    def _get_brakets(self):
        """
        Builds a BraKet space for the current indices
        :return:
        :rtype:
        """
        if len(self.states) == 0:
            return BraKetSpace(
                BasisStateSpace(self.basis, [], mode=BasisStateSpace.StateSpaceSpec.Indices),
                BasisStateSpace(self.basis, [], mode=BasisStateSpace.StateSpaceSpec.Indices)
            )
        arr = self.array #type: SparseArray
        inds = arr.block_inds[1]
        row_inds, col_inds = inds
        bras = self.states.take_subspace(row_inds)
        kets = self.states.take_subspace(col_inds)

        return BraKetSpace(bras, kets)

    def extend_basis(self, states, extend_columns=True):
        """
        Extends the held state space and resizes the held array if need be

        :param states:
        :type states: BasisStateSpace
        :return:
        :rtype:
        """

        # we maintain sorting of the states so that
        # later operations do not require any sorts
        new_states = self.states.union(states, sort=True)
        if new_states is not self.states:
            # this requires us to first find the proper
            cur_vals, _ = self.array.block_data
            row_inds = new_states.find(self.brakets.bras)
            if self.column_space is not None:
                col_inds = self.column_space.find(self.brakets.kets)
            elif extend_columns:
                col_inds = new_states.find(self.brakets.kets)
            else:
                self.column_space = self.states
                col_inds = self.column_space.find(self.brakets.kets)
            self.states = new_states

            new_array = SparseArray.from_data(
                (
                    cur_vals,
                    (
                        row_inds,
                        col_inds
                    )
                ),
                shape=self.shape
                # layout=self.array.fmt,
                # dtype=self.array.dtype
            )
            self.array = new_array

    def _compute_uncached_values(self, func, brakets):
        """

        :param func: A function that can take a braket spec and compute values
        :type func:
        :param brakets: A set of brakets to compute values for
        :type brakets:
        :return:
        :rtype:
        """

        new_terms = self._get_uncached_states(brakets)

        if len(new_terms) > 0:
            # we then calculate values for the new_pos
            row_inds = self.states.find(new_terms.bras)
            col_inds = self.states.find(new_terms.kets)
            # raise Exception(row_inds[:5], col_inds[:5], new_terms.bras.excitations[:5], new_terms.kets.excitations[:5])
            if self.symmetric:
                upper_triangle = np.where(row_inds <= col_inds)[0]
                new_terms = new_terms.take_subspace(upper_triangle)
                row_inds = row_inds[upper_triangle]
                col_inds = col_inds[upper_triangle]
            new_vals = func(new_terms)
            if self.symmetric:
                # this is all done over the total state space rather than the column
                # space to start so that the symmetry can be made manifest
                # but then it needs to be reduced over the column space at the end
                # a basic dupe of this appears later
                flippers = np.where(col_inds != row_inds)[0]
                new_vals = np.concatenate([new_vals, new_vals[flippers]])
                new_row_inds = np.concatenate([row_inds, col_inds[flippers]])
                new_col_inds = np.concatenate([col_inds, row_inds[flippers]])
                row_inds = new_row_inds
                col_inds = new_col_inds

            if self.column_space is not None:
                col_inds = self.column_space.find(self.states.take_subspace(col_inds))

            self.array[row_inds, col_inds] = new_vals

            self._brakets = None # easier than a merge at this point in time...

    def compute_values(self, func, brakets):
        """
        Computes new values into the held `SparseArray` based on the function and brakets provided
        and returns the entire array of values

        :param func: A function that can take a braket spec and compute values
        :type func:
        :param brakets: A set of brakets to compute values for
        :type brakets:
        :return:
        :rtype:
        """

        self._compute_uncached_values(func, brakets)

        row_inds = self.states.find(brakets.bras)
        col_inds = self.states.find(brakets.kets)
        return self.array[row_inds, col_inds]

    def _get_uncached_states(self, brakets):
        """
        Determines which brakets don't have values in the `SparseArray`

        :param brakets:
        :type brakets: BraKetSpace
        :return:
        :rtype:
        """

        # First check for basic compatibility
        if brakets.bras.basis is not self.states.basis:
            raise ValueError("{} with basis {} is incompatible with held basis {}".format(
                type(brakets).__name__,
                brakets.bras.basis,
                self.states.basis
            ))

        # Next find the complement of the bras/kets with the stored states and extend the
        # stored stuff if needed
        comp = brakets.bras.difference(self.states).union(brakets.kets.difference(self.states))
        self.extend_basis(comp)

        # Then pull the row/column indices for the test states and linearize these
        row_inds = self.states.find(brakets.bras)
        col_inds = self.states.find(brakets.kets)
        # and wrap all of this up into a linear index
        flat_idx = np.ravel_multi_index((row_inds, col_inds), self.array.shape)
        # then pull the SparseArray flat idx
        flat_inds, _ = self.array.block_inds
        new_stuff = np.where(np.logical_not(np.isin(flat_idx, flat_inds)))[0]
        return brakets.take_subspace(new_stuff)

    def dot(self, other):
        """
        Performs a dot product between the held SparseArray and another
        StateSpaceMatrix

        :param other: other matrix
        :type other: StateSpaceMatrix
        :return:
        :rtype:
        """
        import copy

        # First check for basic compatibility
        if other.basis is not self.basis:
            raise ValueError("{} with basis {} is incompatible with held basis {}".format(
                type(other).__name__,
                other.basis,
                self.basis
            ))

        # we assume the state spaces remain sorted
        if self.states != other.states:
            # like this step might be quite slow...
            new_states_other = self.states.difference(other.states)
            new_states_self = other.states.difference(self.states)
            if len(new_states_self) > 0:
                self.extend_basis(new_states_self)

            if len(new_states_other) > 0:
                other = copy.copy(other)
                other.extend_basis(new_states_other, extend_columns=False)

        new_array = self.array.dot(other.array)

        return type(self)(self.states, initial_vals=new_array, column_space=other.column_space, symmetric=False)

    def __getitem__(self, item):
        if not isinstance(item, BraKetSpace):
            raise TypeError("{} can only be indexed by a {}".format(
                type(self).__name__,
                BraKetSpace.__name__
            ))

        row_inds = self.states.find(item.bras)
        col_inds = self.states.find(item.kets)
        return self.array[row_inds, col_inds]

    def __setitem__(self, item, vals):
        if not isinstance(item, BraKetSpace):
            raise TypeError("{} can only be indexed by a {}".format(
                type(self).__name__,
                BraKetSpace.__name__
            ))

        row_inds = self.states.find(item.bras)
        col_inds = self.states.find(item.kets)
        self.array[row_inds, col_inds] = vals

    def __repr__(self):
        return "{}({}, basis={})".format(type(self).__name__, self.array, self.states.basis)