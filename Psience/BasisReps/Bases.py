"""
Provides a very general specification for a RepresentationBasis object that can be
used to define matrix representations
"""

import abc, numpy as np, scipy.sparse as sp, itertools as ip
from .StateIndexers import BaseStateIndexer, ArrayStateIndexer, PermutationStateIndexer
from McUtils.Scaffolding import MaxSizeCache
from McUtils.Combinatorics import SymmetricGroupGenerator

__all__ = [
    "RepresentationBasis",
    "SimpleProductBasis"
]

__reload_hook__ = ['.StateIndexers', '.Representations', '.Operators']

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
        :param n_quanta: numbers of quanta (hold over from initial implementation)
        :type n_quanta: int
        """
        self.quanta = n_quanta
        self.generator = function_generator

    @abc.abstractmethod
    def __eq__(self, other):
        raise NotImplementedError("abstract equality not implemented")

    @property
    def dimensions(self):
        """
        Returns the dimensions of the basis

        :return:
        :rtype:
        """
        return (self.quanta,)

    @property
    def ndim(self):
        """
        Returns the number of dimensions of the basis

        :return:
        :rtype:
        """
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
        Generates the momentum matrix up to n-quanta.
        There's one big subtlety to what we're doing here, which is that
          for efficiency reasons we return an entirely real matrix
        The reason for that is we assumed people would mostly use it in the context
          of stuff like pp, pQp, or pQQp, in which case the imaginary part pulls out
          and becomes a negative sign
        We actually use this assumption across _all_ of our representations.

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
    def operator(self, *terms, logger=None, parallelizer=None, chunk_size=None):
        """
        Provides an `Operator` to handle the given terms

        :param terms:
        :type terms:
        :param logger:
        :type logger:
        :param parallelizer:
        :type parallelizer:
        :param chunk_size:
        :type chunk_size:
        :return:
        :rtype:
        """
        from .Operators import Operator
        funcs = [self.operator_mapping[f] if isinstance(f, str) else f for f in terms]
        q = (self.quanta,)
        op = Operator(funcs, q, logger=logger, parallelizer=parallelizer, chunk_size=chunk_size)
        return op
    def representation(self, *terms, logger=None, name=None, parallelizer=None, chunk_size=None, memory_constrained=False):
        """
        Provides a representation of a product operator specified by `terms`

        :param terms:
        :type terms:
        :return:
        :rtype:
        """
        from .Representations import Representation
        q=self.quanta
        return Representation(self.operator(*terms,
                                            logger=logger, parallelizer=parallelizer, chunk_size=chunk_size),
                              self, name=name, memory_constrained=memory_constrained
                              )

    @classmethod
    def _sel_rules_from_rules(cls, rules, max_len=None):
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

        rules = [np.sort(x) for x in rules]
        ndim = len(rules)
        if max_len is None:
            max_len = ndim

        cur_rules = {(0,) * max_len}
        for r in rules:
            new_rules = set()
            for e in cur_rules:
                for s in r:
                    for i in range(max_len):
                        shift = e[i] + s
                        new = e[:i] + (shift,) + e[i+1:]
                        new = tuple(sorted(new, key=lambda l: -abs(l)*10-(1 if l> 0 else 0)))
                        new_rules.add(new)
                        continue


                        # keep things sorted by absolute value
                        # to cut down on the amount of sorting we have to do
                        # later
                        break_time = e[i] == 0
                        abs_shift = abs(shift)
                        if i > 0 and abs_shift > abs(e[i-1]): # shift this up
                            j = i-1
                            while j > 0 and abs_shift > abs(e[j]):
                                j-=1
                            while j > 0 and abs_shift == abs(e[j]) and shift < 0 and e[j] > 0:
                                j -= 1
                            new = e[:j] + (shift,) + e[j:i] + e[i+1:]
                        elif i < max_len - 1 and abs_shift < abs(e[i+1]):
                            j = i+1
                            while j < max_len-1 and abs_shift < abs(e[j]):
                                j+=1
                            while j < max_len-1 and abs_shift == abs(e[j]) and shift > 0 and e[j] < 0:
                                j += 1
                            new = e[:i] + e[i+1:j] + (shift,) + e[j:]
                        else: # doesn't need to move
                            j = None
                            new = e[:i] + (shift,) + e[i+1:]
                        new_rules.add(new)
                        if break_time:
                            break
            # print(cur_rules)
            cur_rules = new_rules

        cur_rules = np.array(list(cur_rules))
        # print(np.where(cur_rules[:, -1] == -1))
        new_rules = []
        for r in cur_rules:
            w = np.where(r==0)
            if len(w) > 0:
                w = w[0]
            if len(w) == 0:
                new_rules.append(tuple(r))
            else:
                new_rules.append(tuple(r[:w[0]]))

        new_rules = list(sorted(new_rules, key=lambda l:len(l)*100 + sum(l)))
        # print(new_rules)

        return new_rules

    def selection_rule_steps(self, *terms):
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

        return rules

    def selection_rules(self, *terms):
        """
        Generates the full set of possible selection rules for terms

        :param terms:
        :type terms:
        :return:
        :rtype:
        """

        # first we get the rules the basis knows about
        rules = self.selection_rule_steps(*terms)

        # print(rules)
        return self._sel_rules_from_rules(rules)

class SimpleProductBasis(RepresentationBasis):
    """
    Defines a direct product basis from a 1D basis.
    Mixed product bases aren't currently supported, but this provides
    at least a sample for how that kind of things could be
    generated.
    """

    _indexer_cache = MaxSizeCache()
    array_indexer_cutoff = 0
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
        """
        Provides the number of dimensions of the basis

        :return:
        :rtype:
        """
        return len(self.bases)

    @property
    def dimensions(self):
        """
        Provides the dimensions of the basis

        :return:
        :rtype:
        """
        return self.quanta

    def __eq__(self, other):
        """
        :param other:
        :type other: SimpleProductBasis
        :return:
        :rtype:
        """
        return all(s==o for s,o in zip(self.bases, other.bases))

    @property
    def quanta(self):
        """
        Provides the quanta in each dimension of the basis

        :return:
        :rtype:
        """
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

    def operator(self, *terms, coeffs=None, axes=None,
                 parallelizer=None, logger=None, chunk_size=None):
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
                                    logger=logger,
                                    chunk_size=chunk_size
                                    )
        else:
            op = Operator(funcs, q, symmetries=labels, selection_rules=sel_rules,
                          parallelizer=parallelizer, logger=logger, chunk_size=chunk_size)
        return op
    def representation(self, *terms, coeffs=None, name=None, axes=None,
                       logger=None, parallelizer=None, chunk_size=None, memory_constrained=False):
        """
        Provides a representation of a product operator specified by _terms_.
        If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.

        :param terms:
        :type terms:
        :return:
        :rtype:
        """
        from .Representations import Representation
        return Representation(
            self.operator(*terms, coeffs=coeffs, axes=axes,
                          parallelizer=parallelizer, logger=logger, chunk_size=chunk_size),
            self, logger=logger, name=name, memory_constrained=memory_constrained
        )
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
