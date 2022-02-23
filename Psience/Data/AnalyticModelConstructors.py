
import numpy as np
from .Helpers import sym, AnaylticModelBase
from .KEData import KEData

__all__ = [
    'AnalyticPotentialConstructor',
    'AnalyticGMatrixConstructor',
    "AnalyticPseudopotentialConstructor"
]

class AnalyticPotentialConstructor(AnaylticModelBase):
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


class AnalyticKineticEnergyConstructor(AnaylticModelBase):
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