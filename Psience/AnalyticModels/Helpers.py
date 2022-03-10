
__all__ = [
    "sym",
    "SymbolicCaller",
    "AnalyticModelBase"
]

import numpy as np

class SympyShim:
    """
    Provides a loader that can either load sympy when requested
    or throw an error if it can't be loaded
    """
    sym = None
    @classmethod
    def _load_sympy(self):
        if self.sym is None:
            import sympy as sym
            self.sym = sym
        return self.sym
    def __getattr__(self, item):
        return getattr(self._load_sympy(), item)
sym = SympyShim()

class SymbolicCaller:
    """
    Delegates to `__call__` through `__getitem__` for symbolics
    """
    def __init__(self, sym):
        self.sym = sym
    def __getitem__(self, item):
        if isinstance(item, int):
            return self.sym(item)
        else:
            return self.sym(*item)

class AnalyticModelBase:
    """
    Provides a base class for analytic models
    """
    @classmethod
    def take_derivs(cls, expr, vars):
        """
        Takes derivatives of `expr` with respect to `vars` even if `expr` is an array

        :param expr:
        :type expr:
        :param vars:
        :type vars:
        :return:
        :rtype:
        """
        if isinstance(expr, list):
            return [cls.take_derivs(e, vars) for e in expr]
        else:
            return [sym.diff(expr, v) for v in vars]

    @classmethod
    def eval_exprs(cls, expr, subs):
        """
        Evaluates `expr` with the given substitutions

        :param expr:
        :type expr:
        :param subs:
        :type subs:
        :return:
        :rtype:
        """
        if isinstance(expr, list):
            return [cls.eval_exprs(e, subs) for e in expr]
        elif isinstance(expr, (int, float, np.integer, np.floating)):
            return expr
        else:
            return expr.subs(subs)

    @classmethod
    def symbol_list(cls, names, instance=None):
        """
        Gets a list of symbols for `names` with a given instance number

        :param names:
        :type names:
        :param instance:
        :type instance:
        :return:
        :rtype:
        """
        instance_tag = "" if instance is None else str(instance)
        return sym.symbols(" ".join(n+instance_tag for n in names))
    @staticmethod
    def symbolic_m(i):
        """
        Provides a symbolic representation of a mass

        :param i:
        :type i:
        :return:
        :rtype:
        """
        return sym.Symbol('m[{i}]'.format(i=i), positive=True)
    @staticmethod
    def symbolic_r(i, j):
        """
        Provides a symbolic representation of a bond length

        :param i:
        :type i:
        :param j:
        :type j:
        :return:
        :rtype:
        """
        if j < i:
            i, j = j, i
        return sym.Symbol('r[{i},{j}]'.format(i=i, j=j), positive=True)
    @staticmethod
    def symbolic_a(i, j, k):
        """
        Provides a symbolic representation of a bond angle

        :param i:
        :type i:
        :param j:
        :type j:
        :param k:
        :type k:
        :return:
        :rtype:
        """
        if k < i:
            i, k = k, i
        return sym.Symbol('a[{i},{j},{k}]'.format(i=i, j=j, k=k), real=True)
    @staticmethod
    def symbolic_t(i, j, k, l):
        """
        Provides a symbolic representation of a dihedral angle

        :param i:
        :type i:
        :param j:
        :type j:
        :param k:
        :type k:
        :param l:
        :type l:
        :return:
        :rtype:
        """
        if l < i:
            i, j, k, l = l, k, j, i
        return sym.Symbol('t[{i},{j},{k},{l}]'.format(i=i, j=j, k=k, l=l), real=True)
    @staticmethod
    def symbolic_y(i, j, k, l):
        """
        Provides a symbolic representation of a book angle

        :param i:
        :type i:
        :param j:
        :type j:
        :param k:
        :type k:
        :param l:
        :type l:
        :return:
        :rtype:
        """
        return sym.Symbol('y[{i},{j},{k},{l}]'.format(i=i, j=j, k=k, l=l), real=True)

    @staticmethod
    def reindex_symbol(symbol, mapping):
        """
        Changes the indices on symbols using the given mapping

        :param symbol:
        :type symbol:
        :param mapping:
        :type mapping:
        :return:
        :rtype:
        """
        name = symbol.name
        for k,v in mapping.items():
            name = name.replace(str(k), str(k)+"$")
        for k,v in mapping.items():
            name = name.replace(str(k)+"$", str(v))
        return sym.Symbol(name, positive=symbol.is_positive, real=symbol.is_real)

    @classmethod
    def lam(cls, i, j, k):
        """
        Provides the `lambda` expression from Frederick and Woywood

        :param i:
        :type i:
        :param j:
        :type j:
        :param k:
        :type k:
        :return:
        :rtype:
        """
        r12 = cls.symbolic_r(i, j)
        r23 = cls.symbolic_r(j, k)
        a = cls.symbolic_a(i, j, k)
        return 1/sym.sin(a)*(1/r12 - sym.cos(a)/r23)