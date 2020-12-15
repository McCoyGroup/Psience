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

    @property
    def ndim(self):
        return self.basis.ndim

    @property
    def excitations(self):
        if self._excitations is None:
            self._excitations = self.as_excitations()
        return self._excitations

    @property
    def indices(self):
        if self._indices is None:
            self._indices = self.as_indices()
        return self._indices

    @property
    def indexer(self):
        if self._indexer is None:
            self._indexer = np.argsort(self.indices)
        return self._indexer

    def find(self, to_search):
        """
        Finds the indices of a set of indices inside the space

        :param to_search: array of ints
        :type to_search: np.ndarray
        :return:
        :rtype:
        """
        if not isinstance(to_search, np.ndarray) and hasattr(to_search, 'indices'):
            to_search = to_search.indices
        return np.searchsorted(self.indices, to_search, sorter=self.indexer)

    def __len__(self):
        return len(self.indices)

    @abc.abstractmethod
    def as_indices(self):
        """
        Returns the index version of the stored states
        :return:
        :rtype: np.ndarray
        """
        return NotImplementedError("abstract base class")

    @abc.abstractmethod
    def as_excitations(self):
        """
        Returns the excitation version of the stored states
        :return:
        :rtype: np.ndarray
        """
        return NotImplementedError("abstract base class")

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
    def take_intersection(self, states):
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
        self._indexer = None

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
                m_pairs = np.unique(np.array(list(ip.product(l_inds, r_inds))), axis=0).T
            else:
                # Get the representation indices that can be coupled under the supplied set of selection rules
                # Currently this is clumsy.
                # We do this by computing transformed states for the whichever space is smaller and finding where this intersects with the larger space
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
        return BraKetSpace(
            BasisStateSpace(self.basis, inds[0], mode=BasisStateSpace.StateSpaceSpec.Indices),
            BasisStateSpace(self.basis, inds[1], mode=BasisStateSpace.StateSpaceSpec.Indices)
        )

    # def __getitem__(self, item):
    #     inds = self.indices[item]
    #     if isinstance(inds, np.ndarray):
    #         return type(self)(self.basis, inds, mode=self.StateSpaceSpec.Indices)
    #     else:
    #         return inds

    def take_intersection(self, states):
        """
        Takes the intersection of self and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """
        found = self.find(states)
        sel = found[found >= 0]
        return self.take_subspace(sel)
    def take_subspace(self, sel):
        return type(self)(self.basis,
                          self.indices[sel],
                          mode=self.StateSpaceSpec.Indices)
    def take_subdimensions(self, inds):
        return type(self)(self.basis.take_subdimensions(inds),
                          self._excitations[:, inds],
                          mode=self.StateSpaceSpec.Excitations)

    def __repr__(self):
        return "{}(nstates={}, basis={})".format(
            type(self).__name__,
            len(self),
            self.basis
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

    def __len__(self):
        return len(np.unique(self.indices))

    @property
    def nstates(self):
        return int(np.product(self.spaces.shape))

    def __iter__(self):
        return iter(self.spaces)
    @property
    def flat(self):
        return self.spaces.flat


    def as_indices(self):
        inds = None
        for space in self.spaces.flat:
            new_inds = space.indices
            if inds is None:
                inds = new_inds
            else:
                inds = np.unique(np.concatenate([inds, new_inds]))
        return inds

    def as_excitations(self):

        inds = None
        for space in self.spaces.flat:
            new_inds = space.excitations
            if inds is None:
                inds = new_inds
            else:
                inds = np.unique(np.concatenate([inds, new_inds], axis=0), axis=0)

        return inds

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
                m_pairs = np.unique(
                    np.concatenate([m_pairs, m], axis=0),
                    axis=0
                )

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
        return BasisStateSpace(self.basis,
                               self.indices,
                               mode=BasisStateSpace.StateSpaceSpec.Indices
                               )

    def take_intersection(self, states):
        """
        Takes the intersection of each held space and the specified states
        :param states:
        :type states:
        :return:
        :rtype:
        """
        def take_inter(space, states=states):
            return space.take_intersection(states)
        new_spaces = np.apply_along_axis(
            take_inter,
            -1,
            self.spaces
        )
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
        return self.take_intersection(subsel)
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

        return np.unique(
            np.array([inds_l, inds_r]).T,
            axis=0
        ).T

    def get_representation_brakets(self,
                                   freqs=None,
                                   freq_threshold=None,
                                   other=None,
                                   selection_rules=None
                                   ):
        return BasisStateSpace.get_representation_brakets(self, freqs=freqs, freq_threshold=freq_threshold)

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

        b1 = BasisStateSpace(self.basis, ind_pairs[0], mode=BasisStateSpace.StateSpaceSpec.Indices)
        b2 = BasisStateSpace(self.basis, ind_pairs[1], mode=BasisStateSpace.StateSpaceSpec.Indices)

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
                new_states = np.unique(new_states, axis=0)
            else:
                new_states = o[np.newaxis, :] + permutations
                mins = np.min(new_states, axis=1)
                well_behaved = np.where(mins >= filter_min)
                new_states = new_states[well_behaved]
                if filter_space is not None:
                    maxes = np.max(new_states, axis=1)
                    well_behaved = np.where(maxes <= filter_max)
                    new_states = new_states[well_behaved]
                new_states = np.unique(new_states, axis=0)
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

        nmodes = space.ndim
        # group selection rules by how many modes they touch
        sel_rule_groups = {}
        for s in selection_rules:
            k = len(s)
            if k in sel_rule_groups:
                sel_rule_groups[k].append(s)
            else:
                sel_rule_groups[k] = [s]
        # this determines how many permutations we have
        num_perms = sum(len(sel_rule_groups[k]) * (nmodes ** k) for k in sel_rule_groups)

        # then we can set up storage for each of these
        permutations = np.zeros((num_perms, nmodes), dtype=int)

        # loop through the numbers of modes and create the appropriate permutations
        prev = []
        for k, g in sel_rule_groups.items():
            g = np.array(g, dtype=int)
            # where to start filling into perms
            i_start = sum(pl * (nmodes ** pk) for pl, pk in prev)
            l = len(g)
            prev.append([l, k])
            if k != 0:  # special case for the empty rule
                inds = np.indices((nmodes,) * k)
                inds = inds.transpose(tuple(range(1, k + 1)) + (0,))
                inds = np.reshape(inds, (-1, k))
                for i, perm in enumerate(inds):
                    uinds = np.unique(perm)
                    if len(uinds) < k:  # must act on totally different modes
                        continue

                    sind = i_start + l * i
                    eind = i_start + l * (i + 1)
                    permutations[sind:eind, perm] = g

        # TODO: get sparsity working, rather than just slowing shit down
        # use_sparse = nmodes > 6
        # use_sparse = use_sparse and max(
        #     sel_rule_groups.keys()) < nmodes // 2  # means we leave the majority of modes untouched
        # if use_sparse:
        #     permutations = SparseArray(sp.csc_matrix(permutations), shape=permutations.shape)

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
        return len(self.bras)

    def __repr__(self):
        return "{}(nstates={})".format(type(self).__name__, len(self))

    def load_non_orthog(self):
        if 'base' not in self._orthogs:
            exc_l, exc_r = self.state_pairs
            self._orthogs['base'] = exc_l == exc_r

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
                sel_bits[n] = np.unique(np.add.outer(s, sel_bits[n])).flatten()  # make the next set of rules

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
        non_orthog = self.get_non_orthog(inds, assume_unique=assume_unique)
        return self.take_subspace(non_orthog), non_orthog

    def apply_sel_rules(self, rules):
        """
        Applies selections rules, assuming we have
        :param rules:
        :type rules:
        :return:
        :rtype:
        """

        all_sels = self.get_sel_rule_filter(rules)
        return self.take_subspace(all_sels), all_sels

    def adjacency_matrix(self):
        """
        Generates the (sparse) unweighted adjacency matrix for the bras & kets
        :return:
        :rtype:
        """

        total_space = BasisStateSpace(
            self.bras.basis,
            np.unique(np.concatenate([self.bras.indices, self.kets.indices])),
            mode=BasisStateSpace.StateSpaceSpec.Indices
        )
        bra_inds = total_space.find(self.bras)
        ket_inds = total_space.find(self.kets)

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
