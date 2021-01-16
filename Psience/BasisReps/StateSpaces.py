"""
Provides a relatively haphazard set of simple classes to keep track of state information.
By providing a single interface here, we can avoid recomputing information over and over.
"""

import numpy as np, itertools as ip, enum, scipy.sparse as sp
import abc

from McUtils.Numputils import SparseArray

__all__ = [
    "BasisStateSpace",
    "BasisMultiStateSpace",
    "SelectionRuleStateSpace",
    "BraKetSpace"
]

class AbstractStateSpace(metaclass=abc.ABCMeta):
    """
    Represents a generalized state space which will provide core
    methods to index into a basis and generate representations
    """

    class StateSpaceSpec(enum.Enum):
        Excitations = "excitations"
        Indices = "indices"

    def __init__(self, basis):
        """
        :param basis:
        :type basis: RepresentationBasis
        """
        self.basis = basis
        self._indices = None
        self._excitations = None
        self._indexer = None
        self._uinds = None

    @property
    def ndim(self):
        return self.basis.ndim

    @property
    def excitations(self):
        if self._excitations is None:
            self._excitations = self.as_excitations()
        return self._excitations
    @excitations.setter
    def excitations(self, exc):
        self._excitations = exc

    @property
    def indices(self):
        if self._indices is None:
            self._indices = self.as_indices()
        return self._indices
    @indices.setter
    def indices(self, inds):
        self._indices = inds

    @property
    def indexer(self):
        if self._indexer is None:
            self._indexer = np.argsort(self.indices)
        return self._indexer
    @indexer.setter
    def indexer(self, idxer):
        self._indexer = idxer

    def find(self, to_search):
        """
        Finds the indices of a set of indices inside the space

        :param to_search: array of ints
        :type to_search: np.ndarray
        :return:
        :rtype:
        """
        if not isinstance(to_search, np.ndarray):
            if hasattr(to_search, 'indices'):
                to_search = to_search.indices
            else:
                to_search = np.array(to_search)
        vals = np.searchsorted(self.indices, to_search, sorter=self.indexer)
        # we have the ordering according to the _sorted_ version of `indices`
        # so now we need to invert that back to the unsorted version
        vals = self.indexer[vals]
        # now because of how searchsorted works, we need to check if the found values
        # truly agree with what we asked for
        bad_vals = self.indices[vals] != to_search
        vals[bad_vals] = -1
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

    def as_unique_indices(self):
        """
        Returns unique indices
        :return:
        :rtype:
        """
        inds = self.indices
        if self._uinds is None:
            _, uinds = np.unique(inds, return_index=True)
            self._uinds = np.sort(uinds)
        return inds[self._uinds,]

    @abc.abstractmethod
    def as_excitations(self):
        """
        Returns the excitation version of the stored states
        :return:
        :rtype: np.ndarray
        """
        return NotImplementedError("abstract base class")

    def as_unique_excitations(self):
        """
        Returns unique excitations
        :return:
        :rtype:
        """
        exc = self.excitations
        if self._uinds is None:
            _, uinds = np.unique(exc, axis=0, return_index=True)
            self._uinds = np.sort(uinds)
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

        eset = set(elements)
        listunique = [unique_element(i, elements.count(i)) for i in eset]
        u = len(elements)
        return list(sorted(list(perm_unique_helper(listunique, [0] * u, u - 1)), reverse=True))

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
        partitions = np.array(
            [list(reversed(x)) + [0] * (ndim - len(x)) for x in cls._accel_asc(n) if len(x) <= ndim]
        )
        sorting = np.lexsort(partitions.T)

        # then concatenate unique permutations for each and cast to numpy
        return np.array(
            sum((cls._unique_permutations(list(x)) for x in partitions[sorting,]), []),
            dtype=int
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

        self._init_states = np.asarray(states, dtype=int)
        if mode is not None and not isinstance(mode, self.StateSpaceSpec):
            mode = self.StateSpaceSpec(mode)
        self._init_state_types = mode
        self._indices = None
        self._excitations = None
        if self.infer_state_inds_type() == self.StateSpaceSpec.Indices:
            self._indices = self._init_states
        else:
            self._excitations = self._init_states
        self._indexer = None

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

        if states_type is self.StateSpaceSpec.Excitations:
            return np.reshape(states, (-1, self.ndim))
        elif states_type is self.StateSpaceSpec.Indices:
            return self.basis.unravel_state_inds(states.flatten())
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
        # print("oops")
        states = self._init_states
        states_type = self.infer_state_inds_type()
        # print("???", states_type)
        if states_type is self.StateSpaceSpec.Excitations:
            return self.basis.ravel_state_inds(np.reshape(states, (-1, self.ndim)))
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

    def take_unique(self):
        """
        Returns only the unique states, but preserves
        ordering and all of that
        :return:
        :rtype:
        """
        if self._indices is not None:
            states = self.unique_indices
            spec = self.StateSpaceSpec.Indices
            new = type(self)(self.basis, states, mode=spec)
            if self._excitations is not None:
                new.excitations = self.unique_excitations
        else:
            states = self.unique_excitations
            spec = self.StateSpaceSpec.Excitations
            new = type(self)(self.basis, states, mode=spec)
            if self._indices is not None:
                new.indices = self.unique_indices
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
                                   freq_threshold=None
                                   ):
        """
        Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
        in this state space.
        Basically just takes all pairs of indices.
        Only returns the upper-triangle indices

        :return:
        :rtype:
        """

        if other is None:
            other = self

        if freq_threshold is None:
            if selection_rules is None:
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
                transf = self.apply_selection_rules(selection_rules, filter_space=other)
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

        return m_pairs

    def get_representation_brakets(self,
                                   other=None,
                                   selection_rules=None,
                                   freqs=None,
                                   freq_threshold=None
                                   ):
        """
        Generates a `BraKetSpace` that can be fed into a `Representation`
        Basically just takes all pairs of indices.
        Only returns the upper-triangle indices

        :return:
        :rtype:
        """

        inds = self.get_representation_indices(other=other,
                                               selection_rules=selection_rules,
                                               freqs=freqs,
                                               freq_threshold=freq_threshold
                                               )
        bras = self.to_single().take_states(inds[0])
        kets = self.to_single().take_states(inds[1])
        return BraKetSpace(bras, kets)

    def take_subspace(self, sel):
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
            subspace = type(self)(
                self.basis,
                self.excitations[sel,],
                mode=self.StateSpaceSpec.Excitations
            )
            if self._indices is not None:
                subspace.indices = self.indices[sel,]
        else:
            subspace = type(self)(
                self.basis,
                self.indices[sel,],
                mode=self.StateSpaceSpec.Indices
            )
            if self._excitations is not None:
                subspace.excitations = self.excitations[sel,]
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
    def take_states(self, states):
        """
        Takes the set of specified states from the space.
        A lot like take_subspace, but operates on states, not indices
        :param states:
        :type states:
        :return:
        :rtype:
        """
        found = self.find(states)
        sel = found[found >= 0]
        return self.take_subspace(sel)

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
        self.spaces = np.asarray(spaces, dtype=object)
        super().__init__(self.basis)

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

    def as_indices(self):
        """
        Pulls the full set indices out of all of the
        held spaces and returns them as a flat vector
        :return:
        :rtype:
        """
        inds = None
        for space in self.spaces.flat:
            new_inds = space.indices
            if inds is None:
                inds = new_inds
            else:
                inds = np.concatenate([inds, new_inds])
                # _, pos = np.unique(cat, return_index=True)
                # inds = cat[np.sort(pos)]
        return inds

    def as_excitations(self):
        """
        Pulls the full set excitations out of all of the
        held spaces and returns them as a flat vector
        :return:
        :rtype:
        """
        exc = None
        for space in self.spaces.flat:
            new_exc = space.excitations
            if exc is None:
                exc = new_exc
            else:
                exc = np.concatenate([exc, new_exc], axis=0)
                # _, pos = np.unique(cat, axis=0, return_index=True)
                # exc = cat[np.sort(pos)]

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
                                   selection_rules=None
                                   ):
        """
        Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
        in this state space.
        Basically just takes all pairs of indices.

        :return:
        :rtype:
        """

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
                                   selection_rules=None
                                   ):
        return BasisStateSpace.get_representation_brakets(self, freqs=freqs, freq_threshold=freq_threshold)

    def to_single(self):
        """
        Condenses the multi state space down to
        a single BasisStateSpace

        :return:
        :rtype:
        """

        states = BasisStateSpace(
            self.basis,
            self.indices,
            mode=BasisStateSpace.StateSpaceSpec.Indices
        )
        states.excitations = self.excitations
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

    def get_representation_indices(self,
                                   other=None,
                                   freqs=None,
                                   freq_threshold=None,
                                   selection_rules=None
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
        for i, s in zip(inds_base, self.flat):
            j = s.indices
            inds_l.append(np.full(len(j), i, dtype=int))
            inds_r.append(j)

        inds_l = np.concatenate(inds_l)
        inds_r = np.concatenate(inds_r)

        # print(inds_base, self.spaces)

        pairs = np.array([inds_l, inds_r]).T
        _, upos = np.unique(pairs, axis=0, return_index=True)
        upairs = pairs[np.sort(upos), ...]
        # raise Exception(len(upos), len(pairs), len(inds_l), len(upairs))
        return upairs.T

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

        og = space.excitations
        excitations = np.full((len(og),), None, dtype=object)

        use_sparse = isinstance(permutations, SparseArray)
        if use_sparse:
            # we drop down to the scipy wrappers until I can improve broadcasting in SparseArray
            og = sp.csc_matrix(og, dtype=int)
            permutations = permutations.data

        if filter_space is not None:
            filter_exc = filter_space.excitations
            # two quick filters
            filter_min = np.min(filter_exc)
            filter_max = np.max(filter_exc)
            # and the slow sieve
            filter_inds = filter_space.indices
        else:
            filter_min = 0
            filter_max = None
            filter_inds = None


        # and now we add these onto each of our states to build a new set of states
        for i, o in enumerate(og):
            if use_sparse:
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
                as_inds = space.basis.ravel_state_inds(new_states) # convert to indices
                dropped = np.setdiff1d(as_inds, filter_inds, assume_unique=True) # the indices that _aren't_ in the filter set
                new_states = np.setdiff1d(as_inds, dropped, assume_unique=True) # the indices that _aren't_ _not_ in the filter set
                excitations[i] = BasisStateSpace(space.basis, new_states, mode=BasisStateSpace.StateSpaceSpec.Indices)
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
                    p = list(p) + padding
                    perms.extend(cls._unique_permutations(p))
        permutations = np.array(perms, dtype=int)

        return permutations

    @classmethod
    def from_rules(cls, space, selection_rules, filter_space=None, iterations=1):
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

        permutations = cls._generate_selection_rule_permutations(space, selection_rules)

        # print(len(permutations))

        # we set up storage for our excitations
        new = space
        for i in range(iterations):
            new = cls._from_permutations(new, permutations, filter_space, selection_rules)

        return new

    def __getitem__(self, item):
        it = self.spaces[item]
        if isinstance(it, np.ndarray):
            # multispace
            init = self.base_space[it]
            it = type(self)(init, it)
        return it

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
        :type bra_space: AbstractStateSpace
        :param ket_space:
        :type ket_space: AbstractStateSpace
        """
        self.bras = bra_space
        self.kets = ket_space
        self.ndim = self.bras.ndim
        self._orthogs = {}
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

    def load_non_orthog(self):
        if 'base' not in self._orthogs:
            exc_l, exc_r = self.state_pairs
            exc_l = np.asarray(exc_l, dtype=int)
            exc_r = np.asarray(exc_r, dtype=int)
            # print(type(exc_l), type(exc_r), exc_l, exc_r)
            womp = np.equal(exc_l, exc_r)
            # print(womp, len(exc_l), len(exc_r))
            self._orthogs['base'] = womp

    def get_non_orthog(self, inds, assume_unique=False):
        """
        Returns whether the states are non-orthogonal under the set of indices.

        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        if not assume_unique:
            inds = np.unique(inds)
        inds = tuple(np.sort(inds))
        if inds not in self._orthogs:
            self.load_non_orthog()
            orthos = self._orthogs['base']
            unused = np.delete(np.arange(len(orthos)), inds)
            self._orthogs[inds] = np.where(np.prod(orthos[unused], axis=0) == 1)[0]
        return self._orthogs[inds]

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

    def apply_non_orthogonality(self, inds, assume_unique=False):
        """
        Takes the bra-ket pairs that are non-orthogonal under the
        indices `inds`
        :param inds:
        :type inds:
        :param assume_unique:
        :type assume_unique:
        :return:
        :rtype:
        """
        non_orthog = self.get_non_orthog(inds, assume_unique=assume_unique)
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
