"""
Provides a very general specification for a RepresentationBasis object that can be
used to define matrix representations
"""
import abc, numpy as np, scipy.sparse as sp, itertools as ip, enum

from McUtils.Numputils import SparseArray

from .Terms import Representation
from .Operators import Operator, ContractedOperator

__all__ = [
    "RepresentationBasis",
    "SimpleProductBasis",
    "BasisStateSpace",
    "BasisMultiStateSpace",
    "SelectionRuleStateSpace"
]

class RepresentationBasis(metaclass=abc.ABCMeta):
    """
    Metaclass for representations.
    Requires concrete implementations of the position and momentum operators.
    """
    name = "Basis"
    def __init__(self, function_generator, n_quanta):
        """

        :param function_generator:
        :type function_generator:
        :param n_quanta: numbers of quanta
        :type n_quanta: int
        """
        self.quanta = n_quanta
        self.generator = function_generator

    @property
    def dimensions(self):
        return (self.quanta,)

    @property
    def ndim(self):
        return 1

    def ravel_state_inds(self, idx):
        """
        Converts state indices from an array of quanta to an array of indices...except in 1D this really isn't doing anything

        :param idx: indices
        :type idx: Iterable[Iterable[int]]
        :return: array of state indices in the basis
        :rtype: tuple[int]
        """

        idx = np.asarray(idx, dtype=int)
        last_dim = idx.shape[-1]
        if last_dim == 1:
            return idx.reshape(idx.shape[:-1])
        else:
            return idx

    def unravel_state_inds(self, idx):
        """
        Converts state indices from an array of ints to an array of quanta...except in 1D this really isn't doing anything

        :param idx: indices
        :type idx: Iterable[int]
        :return: array of state tuples in the basis
        :rtype: tuple[tuple[int]]
        """

        idx = np.asarray(idx, dtype=int)
        last_dim = idx.shape[-1]
        if last_dim == 1:
            return idx
        else:
            return np.expand_dims(idx, axis=-1)

    def __getitem__(self, item):
        if self.generator is None:
            raise ValueError("basis function generator is None (i.e. explicit functions can't be returned)")
        return self.generator(item)
    def __repr__(self):
        unnamed = self.name == "Basis"
        return "{}({}dimension={})".format(
            type(self).__name__,
            "'{}',".format(self.name) if not unnamed else "",
            self.dimensions#self.generator
        )

    @abc.abstractmethod
    def p(self, n):
        """
        Generates the momentum matrix up to n-quanta

        There's one big subtlety to what we're doing here, which is that
          for efficiency reasons we return an entirely real matrix
        The reason for that is we assumed people would mostly use it in the context
          of stuff like pp, pQp, or pQQp, in which case the imaginary part pulls out
          and becomes a negative sign
        We actually use this assumption across _all_ of our representations
        :param n:
        :type n:
        :return:
        :rtype:
        """
        raise NotImplemented

    @abc.abstractmethod
    def x(self, n):
        """
        Generates the coordinate matrix up to n-quanta

        :param n:
        :type n:
        :return:
        :rtype:
        """
        raise NotImplemented

    def I(self, n):
        """
        Generates the identity matrix up to n-quanta

        :param n:
        :type n:
        :return:
        :rtype:
        """
        return sp.csc_matrix(sp.eye(n))

    @property
    def operator_mapping(self):
        return {'x':self.x, 'p':self.p, 'I':self.I}

    selection_rules_mapping = {'x': None, 'p': None, 'I': [0]}
    def operator(self, *terms):
        funcs = [self.operator_mapping[f] if isinstance(f, str) else f for f in terms]
        q = (self.quanta,)
        op = Operator(funcs, q)
        return op
    def representation(self, *terms, logger=None):
        """
        Provides a representation of a product operator specified by 'terms'
        :param terms:
        :type terms:
        :return:
        :rtype:
        """

        q=self.quanta
        return Representation(self.operator(*terms), self, logger=logger)

    @classmethod
    def _sel_rules_from_rules(cls, rules):
        """
        We take the combo of the specified rules, where we take successive products of 1D rules with the
        current set of rules following the pattern that
            1. a 1D change can apply to any index in an existing rule
            2. a 1D change can be appended to an existing rule

        We ensure at each step that the rules remain sorted & duplicates are removed so as to keep the rule sets compact.
        This is done in simple python loops, because doing it with arrayops seemed harder & not worth it for a relatively cheap operation.

        :param rules:
        :type rules:
        :return:
        :rtype:
        """
        from collections import OrderedDict

        cur = set()
        ndim = len(rules)
        # print([tuple(range(ndim)) for _ in range(ndim)])
        for p in ip.product(*(tuple(range(ndim)) for _ in range(ndim))): # loop over different permutations of the operators
            # figure out what inidces are touched by which bits of p
            opgroups = OrderedDict()
            for i in p:
                if i in opgroups:
                    opgroups[i].append(rules[i])
                else:
                    opgroups[i] = [rules[i]]

            # now generate the set of 1D changes from each opgroup
            ichng = []
            for subrules in opgroups.values():
                # print("oooo>", subrules)
                base = np.array(subrules[0])
                for new in subrules[1:]:
                    base = np.unique(np.add.outer(new, base).flatten())
                ichng.append(base)

            # we do a check for the zero case before the product then drop all zeros
            # for speed reasons
            if all(0 in r for r in ichng):
                cur.add(())

            for i, r in enumerate(ichng):
                ichng[i] = r[r != 0]

            # print(">>>>", ichng)
            for t in ip.product(*ichng):
                # print(p)
                t = tuple(sorted(t))
                cur.add(t)

        # for dim, rule_list in enumerate(rules[1:]):
        #     new = set()
        #     # loop over the indices that can be touched by the rule
        #     for p in ip.permutations():
        #         ...
        #
        #     # # loop over the ints held in the rule
        #     # for r in rule_list:
        #     #     # transform the existing rules
        #     #     print(r, cur)
        #     #     for x in cur:
        #     #         for i in range(len(x)):
        #     #             test = tuple(
        #     #                 sorted(
        #     #                     y + r if j == i else y for j, y in enumerate(x)
        #     #                     if y + r != 0
        #     #                 )
        #     #             )
        #     #             if len(test) > 0:
        #     #                 new.add(test)
        #     #         test = tuple(sorted(x + (r,)))
        #     #         new.add(test)
        #     cur = new
        #     # cur.update(new)

        return list(sorted(cur, key=lambda l:len(l)*100 + sum(l)))


    def selection_rules(self, *terms):
        """
        Generates the full set of possible selection rules for terms

        :param terms:
        :type terms:
        :return:
        :rtype:
        """

        # first we get the rules the basis knows about
        rules = [None] * len(terms)
        for i, t in enumerate(terms):
            if t in self.selection_rules_mapping:
                rules[i] = self.selection_rules_mapping[t]
            else:
                raise ValueError("don't know how to infer selection rules from {}".format(
                    t
                ))

        return self._sel_rules_from_rules(rules)
        # TODO: add in to PerturbationTheory stuff + need to add initial_states & final_states for transition moment code... (otherwise _waaay_ too slow)
        transitions_h1 = [
            [-1],
            [1],
            [-3],
            [3],
            [-1, -1, -1],
            [-1, -1, 1],
            [-1, 1, 1],
            [1, 1, 1],
            [1, 2],
            [-1, 2],
            [1, -2],
            [-1, -2]
        ]
    # def __repr__(self):
    #     return "{}('{}')".format(
    #         type(self).__name__,
    #         self.name
    #     )

class SimpleProductBasis(RepresentationBasis):
    """
    Defines a direct product basis from a 1D basis.
    Mixed product bases aren't currently supported.
    """
    def __init__(self, basis_type, n_quanta):
        """
        :param basis_type: the type of basis to do a product over
        :type basis_type: type
        :param n_quanta: the number of quanta for the representations
        :type n_quanta: Iterable[int]
        """
        self.basis_type = basis_type
        self.bases = tuple(basis_type(n) for n in n_quanta)
        super().__init__(self.get_function, None)

    @property
    def ndim(self):
        return len(self.bases)

    @property
    def dimensions(self):
        return self.quanta

    @property
    def quanta(self):
        return tuple(b.quanta for b in self.bases)
    @quanta.setter
    def quanta(self, n):
        if n is not None:
            raise ValueError("{}: '{}' can't be set directly".format(
                type(self).__name__,
                'quanta'
            ))

    @property
    def selection_rules_mapping(self):
        return self.bases[0].selection_rules_mapping

    def ravel_state_inds(self, idx):
        """
        Converts state indices from an array of quanta to an array of indices

        :param idx: indices
        :type idx: Iterable[Iterable[int]]
        :return: array of state indices in the basis
        :rtype: tuple[int]
        """
        idx = np.asarray(idx, dtype=int)
        if idx.ndim == 1:
            idx = idx[np.newaxis]
        return np.ravel_multi_index(idx.T, self.quanta)

    def unravel_state_inds(self, idx):
        """
        Converts state indices from an array of ints to an array of quanta

        :param idx: indices
        :type idx: Iterable[int]
        :return: array of state tuples in the basis
        :rtype: tuple[tuple[int]]
        """

        return np.array(np.unravel_index(idx, self.quanta)).T

    def get_function(self, idx):
        fs = tuple(b[n] for b, n in zip(self.bases, idx))
        return lambda *r, _fs=fs, **kw: np.prod(f(*r, **kw) for f in _fs)

    def operator(self, *terms, coeffs=None, axes=None):
        """
        Builds an operator based on supplied terms, remapping names where possible.
        If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.

        :param terms:
        :type terms:
        :param coeffs:
        :type coeffs:
        :param axes:
        :type axes:
        :return:
        :rtype:
        """

        funcs = [self.bases[0].operator_mapping[f] if isinstance(f, str) else f for f in terms]
        sel_rules = [self.bases[0].selection_rules_mapping[f] if isinstance(f, str) else None for f in terms]
        if any(s is None for s in sel_rules):
            sel_rules = None
        else:
            sel_rules = [np.asarray(s) for s in sel_rules]
        q = self.quanta
        # determine the symmetries up front to make stuff faster
        ids = [hash(f) for f in terms]
        mapping = {k:i for i,k in enumerate(ids)}
        labels = [mapping[k] for k in ids]
        if coeffs is not None:
            op = ContractedOperator(coeffs, funcs, q, axes=axes, symmetries=labels, selection_rules=sel_rules)
        else:
            op = Operator(funcs, q, symmetries=labels, selection_rules=sel_rules)
        return op
    def representation(self, *terms, coeffs=None, axes=None, logger=None):
        """
        Provides a representation of a product operator specified by _terms_.
        If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.

        :param terms:
        :type terms:
        :return:
        :rtype:
        """
        return Representation(self.operator(*terms, coeffs=coeffs, axes=axes), self, logger=logger)
    def x(self, n):
        """
        Returns the representation of x in the multi-dimensional basis with every term evaluated up to n quanta
        Whether this is what we want or not is still TBD
        :param n:
        :type n:
        :return:
        :rtype:
        """
        return self.representation(self.bases[0].x)[:n, :n]
    def p(self, n):
        """
        Returns the representation of p in the multi-dimensional basis with every term evaluated up to n quanta
        Whether this is what we want or not is still TBD
        :param n:
        :type n:
        :return:
        :rtype:
        """
        return self.representation(self.bases[0].p)[:n, :n]

class BasisStateSpace:
    """
    Represents a subspace of states inside a representation basis.
    Useful largely to provide consistent, unambiguous representations of multiple states across
    the different representation-generating methods in the code base.
    """

    class StateSpaceSpec(enum.Enum):
        Excitations = "excitations"
        Indices = "indices"

    def __init__(self, basis, states, mode=None):
        """
        :param basis:
        :type basis: RepresentationBasis
        :param states:
        :type states: Iterable[int]
        :param mode: whether the states were supplied as indices or as excitations
        :type mode: None | str | StateSpaceSpec
        """

        self.basis = basis
        self._init_states = np.asarray(states, dtype=int)
        if mode is not None and not isinstance(mode, self.StateSpaceSpec):
            mode = self.StateSpaceSpec(mode)
        self._init_state_types = mode
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
        return np.searchsorted(self.indices, to_search, sorter=self.indexer)

    def __len__(self):
        return len(self.indices)

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

    # def __getitem__(self, item):
    #     inds = self.indices[item]
    #     if isinstance(inds, np.ndarray):
    #         return type(self)(self.basis, inds, mode=self.StateSpaceSpec.Indices)
    #     else:
    #         return inds

    def __repr__(self):
        return "{}(nstates={}, basis={})".format(
            type(self).__name__,
            len(self),
            self.basis
        )

class BasisMultiStateSpace:
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
        self._indexer = None
        self._indices = None
        self._excitations = None

    @property
    def representative_space(self):
        return self.spaces.flatten()[0]

    @property
    def basis(self):
        return self.representative_space.basis

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

    @property
    def indices(self):
        """
        Returns all of the indices inside all of the held state spaces

        :return:
        :rtype:
        """

        if self._indices is None:
            self._indices = self._get_inds()
        return self._indices

    def _get_inds(self):
        inds = None
        for space in self.spaces.flat:
            new_inds = space.indices
            if inds is None:
                inds = new_inds
            else:
                inds = np.unique(np.concatenate([inds, new_inds]))
        return inds

    @property
    def excitations(self):
        """
        Returns all of the indices inside all of the held state spaces

        :return:
        :rtype:
        """

        if self._excitations is None:
            self._excitations = self._get_exc()
        return self._excitations

    def _get_exc(self):

        inds = None
        for space in self.spaces.flat:
            new_inds = space.excitations
            if inds is None:
                inds = new_inds
            else:
                inds = np.unique(np.concatenate([inds, new_inds], axis=0), axis=0)

        return inds

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
        return np.searchsorted(self.indices, to_search, sorter=self.indexer)

    def __getitem__(self, item):
        it = self.spaces[item]
        if isinstance(it, np.ndarray):
            # multispace
            it = type(self)(it)
        return it

    def get_representation_indices(self,
                                   freqs=None,
                                   freq_threshold=None
                                   ):
        """
        Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
        in this state space.
        Basically just takes all pairs of indices.

        :return:
        :rtype:
        """

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

    def to_single(self):
        return BasisStateSpace(self.basis,
                               self.indices,
                               mode=BasisStateSpace.StateSpaceSpec.Indices
                               )

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

        super().__init__(excitations)
        self.base_space = init_space
        self.sel_rules = selection_rules

    @property
    def basis(self):
        return self.base_space.basis
    @property
    def ndim(self):
        return self.base_space.ndim

    def get_representation_indices(self,
                                   freqs=None,
                                   freq_threshold=None
                                   ):
        """
        This is where this pays dividends, as we know that only the init_space and the held excitations can couple
        which reduces the combinatoric work by a factor of like 2.
        :return:
        :rtype:
        """

        if freq_threshold is not None:
            raise ValueError("Haven't implemented freq. threshold yet...")

        inds_base = self.base_space.indices
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