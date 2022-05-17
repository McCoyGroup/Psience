
__all__ = [
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
    # @classmethod
    # def _load_attr(self, item):
    #     return getattr(self._load_sympy(), item)
    # def tensorcontraction(self, ):
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
    sym = sym

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
            expr = sym.Array(expr)
        return sym.Array([sym.diff(expr, v) for v in vars])

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
        # if isinstance(expr, sym.Array):
        #     return [cls.eval_exprs(e, subs) for e in expr]
        # elif isinstance(expr, (int, float, np.integer, np.floating)):
        #     return expr
        # else:
        if isinstance(expr, list):
            try:
                expr = sym.Array(expr)
            except ValueError:
                pass
        elif isinstance(expr, np.ndarray) and expr.dtype == object:
            expr = expr.tolist()
            try:
                expr = sym.Array(expr)
            except ValueError:
                pass
        if isinstance(expr, (int, float, np.integer, np.floating, np.ndarray)):
            return expr
        elif not isinstance(expr, list):
            return expr.subs(subs)
        else:
            return [cls.eval_exprs(e, subs) for e in expr]

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
    def symbol(base, *args, **kwargs):
        template = "{}["+",".join(["{}"]*len(args))+"]"
        return sym.Symbol(template.format(base, *args), **kwargs)
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
    @classmethod
    def var(cls, *args):
        if len(args) == 2:
            return cls.symbolic_r(*args)
        elif len(args) == 3:
            return cls.symbolic_a(*args)
        elif len(args) == 4:
            return cls.symbolic_t(*args)
        else:
            raise NotImplementedError("huh")

    @classmethod
    def reindex_symbol(cls, symbol, mapping, target_symbols=None):
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
        qual_name = name
        for k,v in mapping.items():
            name = name.replace(str(k)+"$", str(v))
        if target_symbols is not None and name not in target_symbols:
            t, a = name.split("[")
            a = a.strip("]").split(",")
            inds = tuple(int(x) for x in a)
            if t == "r":
                test = cls.symbolic_r(inds[1], inds[0])
                if test.name in target_symbols:
                    name = test.name
            elif t == "a":
                test = cls.symbolic_a(inds[2], inds[1], inds[0])
                if test.name in target_symbols:
                    name = test.name
            elif t == "t":
                test = cls.symbolic_t(inds[3], inds[2], inds[1], inds[0])
                if test.name in target_symbols:
                    name = test.name

        return [sym.Symbol(qual_name, positive=symbol.is_positive, real=symbol.is_real), sym.Symbol(name, positive=symbol.is_positive, real=symbol.is_real)]

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

    @classmethod
    def is_identity(cls, A):
        N = len(A)
        for i in range(N):
            m = len(A[i])
            if m != N:
                return False
            for j in range(m):
                if (j != i and A[i][j] != 0) or (j == i and A[i][j] != 1):
                    return False
        return True
    @classmethod
    def transpose(cls, A):
        if isinstance(A, list):
            A = sym.Array(A)
        return sym.transpose(A)
    @classmethod
    def dot(cls, a, b, axes=None):
        if isinstance(a, (list, sym.Matrix)):
            a = sym.Array(a)
        if isinstance(b, (list, sym.Matrix)):
            b = sym.Array(b)
        if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
            if axes is None:
                return np.dot(a, b)
            else:
                return np.tensordot(a, b, axes=axes)
        else:
            if isinstance(a, np.ndarray):
                a = sym.Array(a)
            if isinstance(b, np.ndarray):
                b = sym.Array(b)
            if axes is None:
                axes = (a.rank()-1, a.rank())
            else:
                axa, axb = axes
                if isinstance(axb, int):
                    axb += a.rank()-1
                else:
                    ar = a.rank()-1
                    axb = [ar + x for x in axb]
                axes = [axa, axb]
            return sym.tensorcontraction(sym.tensorproduct(a, b), axes)
    @classmethod
    def contract(cls, a, axes):
        if isinstance(a, list):
            a = sym.Array(a)
        return sym.tensorcontraction(a, axes)