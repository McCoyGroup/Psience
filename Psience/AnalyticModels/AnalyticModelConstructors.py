
import numpy as np
from .Helpers import sym, AnalyticModelBase
from ..Data import KEData

__all__ = [
    'AnalyticPotentialConstructor',
    'AnalyticKineticEnergyConstructor'
]

class AnalyticPotentialConstructor(AnalyticModelBase):
    """

    """
    @classmethod
    def symbolic_morse(cls, instance=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        return cls.morse(*cls.symbol_list(["De", "a", "r", "re"], instance=instance))
    @staticmethod
    def morse(De, a, r, re):
        """
        :param De:
        :type De:
        :param a:
        :type a:
        :param r:
        :type r:
        :return:
        :rtype: sym.expr
        """
        return De * (1 - sym.exp(-a * (r-re)))**2
    @staticmethod
    def harm(k, x, x_e):
        """

        :param De:
        :type De:
        :param a:
        :type a:
        :param r:
        :type r:
        :return:
        :rtype: sym.expr
        """
        return k*(x-x_e)**2
    @classmethod
    def symbolic_harmonic(cls, instance=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        return cls.harm(*cls.symbol_list(["k", "q", "qe"], instance=instance))

class AnalyticKineticEnergyConstructor(AnalyticModelBase):
    """
    Provides G and V' elements from Frederick and Woywood
    """

    @classmethod
    def _get_coord_key(cls, inds1:'Iterable[int]', inds2:'Iterable[int]', coord_types=None):
        # we infer which type of element to use and return that
        # with appropriate substitutions
        if coord_types is None:
            coord_types = [None, None]
        coord_types = list(coord_types)
        if coord_types[0] is None:
            coord_types[0] = cls.infer_coord_type(inds1)
        if coord_types[1] is None:
            coord_types[1] = cls.infer_coord_type(inds2)
        type1, type2 = coord_types
        sorting = ['r', 'a', 't', 'y']
        if sorting.index(coord_types[0]) > sorting.index(coord_types[1]):
            inds1, inds2 = inds2, inds1
            type1, type2 = type2, type1
        mapping = {i:k for i,k in zip(inds1, tuple(range(1, len(inds1) + 1)))}
        n = len(inds1) + 1
        for i in inds2:
            if i not in inds1:
                mapping[i] = n
                n += 1

        return ((type1, type2), tuple(mapping[i] for i in inds1), tuple(mapping[i] for i in inds2)), mapping

    @classmethod
    def g(cls, inds1:'Iterable[int]', inds2:'Iterable[int]', coord_types=None):
        key, mapping = cls._get_coord_key(inds1, inds2, coord_types=coord_types)
        try:
            expr = KEData[(key,)][0]
        except KeyError:
            expr = 0
        if not isinstance(expr, (int, np.integer)):
            subs = tuple((s, cls.reindex_symbol(s, mapping)) for s in expr.free_symbols)
            expr = expr.subs(subs)
        return expr

    @classmethod
    def vp(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None):
        key, mapping = cls._get_coord_key(inds1, inds2, coord_types=coord_types)
        try:
            expr = KEData[(key,)][1]
        except KeyError:
            expr = 0
        if not isinstance(expr, (int, np.integer)):
            subs = tuple((s, cls.reindex_symbol(s, mapping)) for s in expr.free_symbols)
            expr = expr.subs(subs)
        return expr

    @classmethod
    def infer_coord_type(cls, inds):
        if len(inds) == 2:
            return "r"
        elif len(inds) == 3:
            return "a"
        elif len(inds) == 4:
            return "t"
        else:
            raise ValueError("too many indices")

    @classmethod
    def _g(cls, inds1:'Iterable[int]', inds2:'Iterable[int]'):
        inds1 = list(inds1)
        inds2 = list(inds2)
        if len(inds1) < len(inds2):
            iter_inds = inds1
            rem_inds = inds2
        else:
            iter_inds = inds2
            rem_inds = inds1
        d1 = cls.displacement_vectors(iter_inds)
        d2 = cls.displacement_vectors(rem_inds)
        g_contribs = []
        for n,i in enumerate(iter_inds):
            try:
                m = rem_inds.index(i)
            except ValueError:
                pass
            else:
                # raise Exception(d1[n], d2[m])
                a = d1[n]
                b = d2[m]
                # print(a, b)
                g_contribs.append(1/cls.symbolic_m(i)*a.dot(b))

        g = sum(g_contribs)
        if isinstance(g, sym.Expr):
            g = g.expand().simplify().expand()#.subs(sym.Abs, sym.Id).simplify()
        return g

class AnalyticModel:
    """
    Provides a symbolic representation of an analytically evaluatable Hamiltonian
    which can be used to get derived expressions to evaluate.
    """

    def __init__(self, coordinates, potential):
        self.coords = coordinates
        self._syms = None
        self._base_g = None
        self._g = None
        self.pot = potential
    @property
    def internal_coordinates(self):
        if self._syms is None:
            self._load_symbols()
        return self._syms[0]

    def _load_symbols(self):
        sym_list = []
        sym_set = set()
        for x in self.coords:
            for s in x.free_symbols:
                if s not in sym_set:
                    sym_set.add(s)
                    sym_list.append(s)
        self._syms = (tuple(sym_list), tuple(self._parse_symbol(s) for s in sym_list))

    def _parse_symbol(self, sym):
        name = sym.name
        t, a = name.split("[")
        a = a.strip("]").split(",")
        return (t, tuple(int(x) for x in a))

    def jacobian(self):
        ics = self.internal_coordinates
        return [AnalyticModelBase.take_derivs(c, ics) for c in self.coords]
    def _base_gmat(self):
        if self._base_g is None:
            if self._syms is None:
                self._load_symbols()
            symlist = self._syms[1]
            return [[AnalyticKineticEnergyConstructor.g(a[1], b[1], coord_types=[a[0], b[0]]) for b in symlist] for a in symlist]
        return self._base_g
    def g(self, order=0):
        # Gmatrix elements will basically involve taking some
        # kind of direct product of coordinate reps
        if self._g is None:
            J_t = self.jacobian()
            G = self._base_gmat()
            def transpose(A):
                n = len(A)
                m = len(A[0])
                return [[A[i][k] for i in range(n)] for k in range(m)]
            def dot(a, b):
                if isinstance(a[0], list): # matmul
                    b = transpose(b)
                    return [[dot(a[i], b[j]) for j in range(len(b))]  for i in range(len(a))]
                else: # vector dot
                    return sum(a*b for a,b in zip(a,b))
            J = transpose(J_t)
            self._g = dot(dot(J_t, G), J)
        if order > 0:
            raise NotImplementedError("need to handle coordinate transforms")
            ics = self.internal_coordinates
            for i in range(order):
                ...


    def v(self, order=2):
        # we provide a Taylor series expansion of the potential
        raise NotImplementedError("need to handle composite coordinates")
    def vp(self, order=0):
        raise NotImplementedError("need to handle composite coordinates")

    def r(self, i, j):
        return AnalyticModelBase.symbolic_r(i, j)
    def a(self, i, j, k):
        return AnalyticModelBase.symbolic_a(i, j, k)
    def t(self, i, j, k, l):
        return AnalyticModelBase.symbolic_t(i, j, k, l)
    def y(self, i, j, k, l):
        return AnalyticModelBase.symbolic_t(i, j, k, l)