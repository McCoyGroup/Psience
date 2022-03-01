
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
    def __init__(self, sym):
        self.sym = sym
    def __getitem__(self, item):
        if isinstance(item, int):
            return self.sym(item)
        else:
            return self.sym(*item)

class AnalyticModelBase:
    @classmethod
    def take_derivs(cls, expr, vars):
        if isinstance(expr, list):
            return [cls.take_derivs(e, vars) for e in expr]
        else:
            return [sym.diff(expr, v) for v in vars]

    @classmethod
    def eval_exprs(cls, expr, subs):
        if isinstance(expr, list):
            return [cls.eval_exprs(e, subs) for e in expr]
        elif isinstance(expr, (int, float, np.integer, np.floating)):
            return expr
        else:
            return expr.subs(subs)

    @classmethod
    def symbol_list(cls, names, instance=None):
        instance_tag = "" if instance is None else str(instance)
        return sym.symbols(" ".join(n+instance_tag for n in names))
    @staticmethod
    def symbolic_m(i):
        return sym.Symbol('m[{i}]'.format(i=i), positive=True)
    @staticmethod
    def symbolic_r(i, j):
        if j < i:
            i, j = j, i
        return sym.Symbol('r[{i},{j}]'.format(i=i, j=j), positive=True)
    @staticmethod
    def symbolic_a(i, j, k):
        if k < i:
            i, k = k, i
        return sym.Symbol('a[{i},{j},{k}]'.format(i=i, j=j, k=k), real=True)
    @staticmethod
    def symbolic_t(i, j, k, l):
        if l < i:
            i, j, k, l = l, k, j, i
        return sym.Symbol('t[{i},{j},{k},{l}]'.format(i=i, j=j, k=k, l=l), real=True)
    @staticmethod
    def symbolic_y(i, j, k, l):
        return sym.Symbol('y[{i},{j},{k},{l}]'.format(i=i, j=j, k=k, l=l), real=True)

    @staticmethod
    def reindex_symbol(symbol, mapping):
        name = symbol.name
        for k,v in mapping.items():
            name = name.replace(str(k), str(k)+"$")
        for k,v in mapping.items():
            name = name.replace(str(k)+"$", str(v))
        return sym.Symbol(name, positive=symbol.is_positive, real=symbol.is_real)

    @classmethod
    def lam(cls, i, j, k):
        r12 = cls.symbolic_r(i, j)
        r23 = cls.symbolic_r(j, k)
        a = cls.symbolic_a(i, j, k)
        return 1/sym.sin(a)*(1/r12 - sym.cos(a)/r23)