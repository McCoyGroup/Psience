"""
Provides a very general specification for a RepresentationBasis object that can be
used to define matrix representations
"""

import abc, numpy as np, scipy.sparse as sp, itertools as ip
from .StateIndexers import BaseStateIndexer, ArrayStateIndexer, PermutationStateIndexer
from McUtils.Scaffolding import MaxSizeCache

__all__ = [
    "RepresentationBasis",
    "SimpleProductBasis"
]

__reload_hook__ = ['.StateIndexers', '.Terms', '.Operators']

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
        from .Operators import Operator
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
        from .Terms import Representation
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
        # transitions_h1 = [
        #     [-1],
        #     [1],
        #     [-3],
        #     [3],
        #     [-1, -1, -1],
        #     [-1, -1, 1],
        #     [-1, 1, 1],
        #     [1, 1, 1],
        #     [1, 2],
        #     [-1, 2],
        #     [1, -2],
        #     [-1, -2]
        # ]
    # def __repr__(self):
    #     return "{}('{}')".format(
    #         type(self).__name__,
    #         self.name
    #     )

class SimpleProductBasis(RepresentationBasis):
    """
    Defines a direct product basis from a 1D basis.
    Mixed product bases aren't currently supported, but this provides
    at least a sample for how that kind of things could be
    generated.
    """

    _indexer_cache = MaxSizeCache()
    array_indexer_cutoff = 6
    def __init__(self, basis_type, n_quanta, indexer=None):
        """
        :param basis_type: the type of basis to do a product over
        :type basis_type: type
        :param n_quanta: the number of quanta for the representations
        :type n_quanta: Iterable[int]
        :param indexer: an object that can turn state specs into indices and indices into state specs
        :type indexer: BaseStateIndexer
        """
        self.basis_type = basis_type
        self.bases = tuple(basis_type(n) for n in n_quanta)
        if indexer is None:
            if self.ndim > self.array_indexer_cutoff: # use the faster version when there are few permutations
                if self.ndim not in self._indexer_cache:
                    indexer = PermutationStateIndexer(self.ndim)
                    self._indexer_cache[self.ndim] = indexer
                else:
                    indexer = self._indexer_cache[self.ndim]
            else:
                dims = self.dimensions
                if dims not in self._indexer_cache:
                    indexer = ArrayStateIndexer(dims)
                    self._indexer_cache[dims] = indexer
                else:
                    indexer = self._indexer_cache[dims]
        self.indexer = indexer
        super().__init__(self.get_function, None)

    def to_state(self, serializer=None):
        return {
            'basis': self.basis_type,
            'quanta': self.quanta,
            'indexer': self.indexer
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(serializer.deserializer(data['basis']),
                   serializer.deserializer(data['quanta']),
                   indexer=serializer.deserializer(data['indexer']))

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

        return self.indexer.to_indices(idx)

    def unravel_state_inds(self, idx):
        """
        Converts state indices from an array of ints to an array of quanta

        :param idx: indices
        :type idx: Iterable[int]
        :return: array of state tuples in the basis
        :rtype: tuple[tuple[int]]
        """

        return self.indexer.from_indices(idx)

        # return np.array(np.unravel_index(idx, self.quanta)).T

    def get_function(self, idx):
        fs = tuple(b[n] for b, n in zip(self.bases, idx))
        return lambda *r, _fs=fs, **kw: np.prod(f(*r, **kw) for f in _fs)

    def operator(self, *terms, coeffs=None, axes=None, parallelizer=None, logger=None):
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
        from .Operators import Operator, ContractedOperator

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
            op = ContractedOperator(coeffs, funcs, q, axes=axes, symmetries=labels,
                                    selection_rules=sel_rules,
                                    parallelizer=parallelizer,
                                    logger=logger
                                    )
        else:
            op = Operator(funcs, q, symmetries=labels, selection_rules=sel_rules, parallelizer=parallelizer, logger=logger)
        return op
    def representation(self, *terms, coeffs=None, axes=None, logger=None, parallelizer=None):
        """
        Provides a representation of a product operator specified by _terms_.
        If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.

        :param terms:
        :type terms:
        :return:
        :rtype:
        """
        from .Terms import Representation
        return Representation(
            self.operator(*terms, coeffs=coeffs, axes=axes, parallelizer=parallelizer, logger=logger), self, logger=logger)
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

    def take_subdimensions(self, dims):
        """
        Casts down to lower dimensional space

        :param dims:
        :type dims:
        :return:
        :rtype:
        """
        qq = self.quanta
        return type(self)(self.basis_type, [qq[d] for d in dims])

    def get_state_space(self, quanta):
        from .StateSpaces import BasisStateSpace
        return BasisStateSpace.from_quanta(self, quanta)
