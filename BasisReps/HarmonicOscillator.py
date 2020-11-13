"""
Provides representations based on starting from a harmonic oscillator basis
"""

__all__ = [
    "HarmonicOscillatorBasis",
    "HarmonicOscillatorProductBasis"
]

import scipy.sparse as sp, numpy as np, functools, itertools

from .Bases import *
from .Operators import Operator, ContractedOperator
from McUtils.Data import WavefunctionData

class HarmonicOscillatorBasis(RepresentationBasis):
    """
    Provides a concrete implementation of RepresentationBasis using the H.O.
    Need to make it handle units a bit better.
    Currently 1D, need to make multidimensional in the future.
    """
    name='HarmonicOscillator'
    def __init__(self, n_quanta, m=None, re=None, dimensionless=True):
        self.dimensionless = dimensionless
        if not dimensionless and any(x is None for x in (m, re)):
            raise ValueError("if not dimensionless, parameters cannot be 'None'")
        self.data = WavefunctionData['HarmonicOscillator']
        # raise Exception(self.data)
        super().__init__(self.data, n_quanta)

    selection_rules_mapping = {
        'x':[-1, 1],
        'p':[-1, 1],
        "I":[0]
    }

    @functools.lru_cache(maxsize=128)
    def p(self, n):
        mat = self.pmatrix_ho(n)
        if not self.dimensionless:
            raise NotImplementedError
        return mat
    @staticmethod
    def pmatrix_ho(n):
        """
        There's one big subtlety to what we're doing here, which is that
          for efficiency reasons we return an entirely real matrix
        The reason for that is we assumed people would mostly use it in the context
          of stuff like pp, pQp, or pQQp, in which case the imaginary part pulls out
          and becomes a negative sign
        We actually use this assumption across _all_ of our representations
        :param n:
        :type n:
        :return:
        :rtype: sp.csr_matrix
        """
        # the imaginary terms pull out and just become a negative sign
        ar = 1 / np.sqrt(2) * np.sqrt(np.arange(1, n))
        bands = [
            [  ar,  1],
            [ -ar, -1]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

    @functools.lru_cache(maxsize=128)
    def x(self, n):
        return self.qmatrix_ho(n)
    @staticmethod
    def qmatrix_ho(n):
        """

        :param n:
        :type n:
        :return:
        :rtype: sp.csr_matrix
        """

        ar = 1 / np.sqrt(2) * np.sqrt(np.arange(1, n))
        bands = [
            [ar,  1],
            [ar, -1]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

    def operator(self, *terms):
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

        if all(isinstance(t, str) for t in terms):
            term_eval = HarmonicProductOperatorTermEvaluator(terms).diag
            def computer(idx, term_eval=term_eval): # for because of how Operator works...
                if len(idx) > 1:
                    raise ValueError("multidimensional index by one-dimensional basis ?_?")
                return term_eval, term_eval.selection_rules
            return Operator(computer, (self.quanta,), prod_dim=1)
        else:
            return super().operator(*terms)

        # funcs = [self.bases[0].operator_mapping[f] if isinstance(f, str) else f for f in terms]
        # q = self.quanta
        # # determine the symmetries up front to make stuff faster
        # ids = [hash(f) for f in terms]
        # mapping = {k:i for i,k in enumerate(ids)}
        # labels = [mapping[k] for k in ids]
        # if coeffs is not None:
        #     op = ContractedOperator(coeffs, funcs, q, axes=axes, symmetries=labels)
        # else:
        #     op = Operator(funcs, q, symmetries=labels)
        # return op

    # def representation(self, *terms):
    #     """
    #     Provides a representation of a product operator specified by 'terms'
    #     :param terms:
    #     :type terms:
    #     :return:
    #     :rtype:
    #     """
    #
    #     q=self.quanta
    #     return TermComputer(self.operator(*terms), q)

class HarmonicOscillatorProductBasis(SimpleProductBasis):

    """
    Tiny, tiny layer on `SimpleProductBasis` that makes use of some analytic work done
    to support representations of `p` and `x`.
    """

    def __init__(self, n_quanta):
        super().__init__(HarmonicOscillatorBasis, n_quanta)

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

        if all(isinstance(t, str) for t in terms):
            computer = HarmonicProductOperatorTermEvaluator(terms)
            if coeffs is None:
                op = Operator(computer, self.quanta, prod_dim=len(terms))
            else:
                op = ContractedOperator(coeffs, computer, self.quanta, prod_dim=len(terms))
            return op
        else:
            return super().operator(*terms, coeffs=coeffs, axes=axes)


class HarmonicProductOperatorTermEvaluator:
    """
    A simple class that can be used to directly evaluate any operator built as a product of `p` and `x` terms.
    Automatically dispatches to provide different permutations.
    """
    _cache = {}
    def __init__(self, terms):
        """
        :param terms: list of 'x' and 'p' terms
        :type terms:
        """
        self.terms = terms

        p_pos = []
        for i, o in enumerate(terms):
            if o not in {"q", "p", "x"}:
                raise ValueError("Don't know what do with operator {}".format(o))
            elif o == "p":
                p_pos.append(i)

        self.p_pos = p_pos

        self.i_phase = (-1) ** (len(p_pos) // 2)  # the -1 terms coming from the complex part of the phase
        # the residual complex part, in general we need this to be zero
        self.is_complex = len(p_pos) % 2 == 1
        if self.is_complex:
            raise ValueError("For efficiency, complex-valued operators not supported right now")

    def take_suboperator(self, inds):
        """
        Subsamples terms and returns an evaluator for the corresponding operator.
        Basically just determines how `inds` gets broken down into its component pieces.
        These pieces get strung into 1D evaluators and then pushed into a `ProdOp`.


        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        from collections import OrderedDict

        if len(inds) != len(self.terms):
            raise ValueError("number of indices specified {} and number of terms {} must agree (from {} into {})".format(
                len(inds),
                len(self.terms),
                inds,
                self.terms
            ))

        ind_map = OrderedDict()
        for i, t in zip(inds, self.terms):
            if i in ind_map:
                ind_map[i].append(t)
            else:
                ind_map[i] = [t]

        term_spec = tuple(tuple(t) for t in ind_map.values())

        if term_spec in self._cache:
            prop = self._cache[term_spec]
        else:
            terms_1D = [self.TermEvaluator1D.load_cached(terms) for terms in term_spec]
            uinds = list(ind_map.keys())
            prop = self.ProdOp(terms_1D, uinds, self.i_phase, self.is_complex)
            self._cache[term_spec] = prop

        return prop
    @property
    def diag(self):
        return self.take_suboperator([0]*len(self.terms))

    def __getitem__(self, inds):
        """
        Syntactic sugar to delegate to `take_suboperator`

        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        return self.take_suboperator(inds)

    def __call__(self, inds):
        """
        Syntactic sugar to play nicely with `Operator`

        :param inds:
        :type inds:
        :return:
        :rtype:
        """
        op = self[inds]
        return op, op.selection_rules

    class ProdOp:
        """
        Simple subclass to manage product operator evaluation
        """
        def __init__(self, ops, uinds, i_phase, is_complex):
            self.ops = ops
            self.uinds = uinds
            self.i_phase = i_phase
            self.is_complex = is_complex
        def __repr__(self):
            return "{}({})".format(
                type(self).__name__,
                "*".join("({})".format(",".join(o.terms)) for o in self.ops)
            )
        def eval_states(self, states):

            # print(states)

            chunk = None
            for s, op in zip(states, self.ops):
                # print(op)
                blob = np.asarray(op(s))
                # print(states)
                if chunk is None:
                    if isinstance(blob, (int, np.integer, float, np.floating)):
                        chunk = np.array([blob])
                    else:
                        chunk = np.asarray(blob)
                else:
                    chunk = chunk * blob

            chunk = self.i_phase * chunk
            # if complex...do what?
            return chunk

        def __call__(self, states):
            # print("?????????")
            return self.eval_states(states)

        @property
        def selection_rules(self):
            return [o.selection_rules for o in self.ops]

    class TermEvaluator1D:
        """
        1D evaluator for terms looking like `x`, `p`, `q`, etc.
        All of the overall `(-i)^N` info is in the `ProdOp` class that's expected to hold this.
        Only maintains phase info & calculates elements.
        """
        _evaluators = {}
        def __init__(self, terms):
            self.terms = terms
            self.N = len(terms)
            self.generators = {}
            p_pos = []
            for i, o in enumerate(terms):
                if o == "p":
                    p_pos.append(i)
            self.p_pos = tuple(p_pos)

        def __repr__(self):
            return "{}({})".format(
                type(self).__name__,
                ", ".join(self.terms)
            )

        @classmethod
        def load_cached(cls, terms):
            if terms in cls._evaluators:
                eval = cls._evaluators[terms]
            else:
                eval = cls(terms)
                cls._evaluators[terms] = eval

            return eval

        @property
        def selection_rules(self):
            return list(range(-self.N, self.N+1, 2))

        # @classmethod
        # def rho_element_generator(cls, op_seq, p_pos):
        #     """
        #     A way to represent the elements of an operator composed of a sequence of `p` and `x` terms.
        #     Returns a callable so that key parts (esp. type-checking stuff) can be reused.
        #
        #     :param op_seq: sequence of position and momentum operators
        #     :type op_seq:
        #     :return: function that makes a function to get the rho coefficients...
        #     :rtype:
        #     """
        #
        #     i_phase = (-1) ** (len(p_pos) // 2)  # the -1 terms coming from the complex part of the phase
        #     # the residual complex part
        #     is_complex = len(p_pos) % 2 == 1
        #     def get_terms(a, N=len(op_seq), i_phase=i_phase, is_complex=is_complex, sel=p_pos):
        #         return cls._rho_term_generator(a, N, i_phase, is_complex, sel)
        #
        #     return get_terms

        def __call__(self, states):
            return self.evaluate_state_terms(states)

        def evaluate_state_terms(self, states):
            """
            Evaluates terms coming from different state excitations.
            Doesn't do any pre-filtering, since that's expected to be in the caller.

            :param states:
            :type states:
            :return:
            :rtype:
            """

            deltas = states[1] - states[0]
            delta_vals = np.unique(deltas)
            delta_sels = [np.where(deltas==a) for a in delta_vals]

            # print(">>>>", delta_vals)

            biggo = np.zeros(states[0].shape)
            for a, s in zip(delta_vals, delta_sels):
                gen = self.load_generator(a)
                og = states[0][s]
                # print("  ", a, s, og, states[1][s])
                biggo[s] = gen(og)

            return biggo

        def load_generator(self, a):
            # print(a)
            if a not in self.generators:
                gen = self.rho_term_generator(a, self.N, self.p_pos)
                self.generators[a] = gen
            else:
                gen = self.generators[a]
            return gen

        @classmethod
        def rho_term_generator(cls, a, N, sel):
            """
            Returns a function to be called on a quantum number to get the coefficient associated with exciting that mode by `a` quanta over
            `N` steps w/ phase info coming from where the momenta-terms are.

            :param a:
            :type a:
            :param N:
            :type N:
            :param i_phase:
            :type i_phase:
            :param is_complex:
            :type is_complex:
            :param sel:
            :type sel:
            :return:
            :rtype:
            """

            if abs(a) > N or (N - a) % 2 != 0:
                def terms(ni):
                    if isinstance(ni, (int, np.integer)):
                        return 0
                    else:
                        return np.zeros(ni.shape)
                return terms

            if a < 0:
                up_steps = (N - abs(a)) // 2
                down_steps = (N + abs(a)) // 2
            else:
                down_steps = (N - abs(a)) // 2
                up_steps = (N + abs(a)) // 2

            # this can be sped up considerably by using cleverer integer partitioning stuff and whatno
            partitions = np.unique(
                np.array(list(
                    itertools.permutations(
                        [-1] * down_steps
                        + [1] * up_steps
                    )
                )),
                axis=0
            )

            # print(a, down_steps, up_steps, partitions)
            # determine the 'phase' of each 'path'
            phases = np.prod(partitions[:, sel], axis=1) / np.sqrt(2) ** N
            # print("   > p", sel, phases)
            # then compute the 'paths' themselves
            unitized = np.abs(np.sign(partitions - 1))
            paths = np.cumsum(partitions, axis=1) + unitized

            def terms(ni, phases=phases, paths=paths):
                return cls.rho(phases, paths, ni)

            return terms

        @classmethod
        def rho(cls, phases, paths, ni):
            """

            :param phases:
            :type phases:
            :param paths:
            :type paths:
            :param ni:
            :type ni:
            :return:
            :rtype:
            """
            # finally define a function that will apply the the "path" to a quantum number and add everything up
            if not isinstance(ni, (int, np.integer)):
                # print("   ??", ni)
                ni = np.asarray(ni)
                # we need to make the broadcasting work
                # requires basically that we copy the shape from ni to the beginning for phases and paths
                # so basically we figure out how many (1) we need to insert at the end of phases and paths
                # and add a (1) to the start of ni
                nidim = ni.ndim
                for i in range(nidim):
                    phases = np.expand_dims(phases, axis=-1)
                    paths = np.expand_dims(paths, axis=-1)
                #
                # phases = np.broadcast_to(phases, phases.shape + (1,)*nidim)
                # paths = np.broadcast_to(paths, paths.shape + (1,)*nidim)
                ni = np.expand_dims(ni, 0)

            path_displacements = paths + ni
            path_displacements[path_displacements < 0] = 0
            path_terms = np.sqrt(np.prod(path_displacements, axis=1))

            return np.sum(phases * path_terms, axis=0)