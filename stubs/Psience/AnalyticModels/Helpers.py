__all__ = ['SymbolicCaller', 'AnalyticModelBase']
import numpy as np

class SympyShim:
    """
    Provides a loader that can either load sympy when requested
    or throw an error if it can't be loaded
    """
    sym = None

    @classmethod
    def _load_sympy(self):
        ...

    def __getattr__(self, item):
        ...
sym = SympyShim()

class SymbolicCaller:
    """
    Delegates to `__call__` through `__getitem__` for symbolics
    """

    def __init__(self, sym):
        ...

    def __getitem__(self, item):
        ...

class AnalyticModelBase:
    """
    Provides a base class for analytic models
    """
    sym = sym

    @classmethod
    def get_numeric_types(cls):
        ...

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
        ...

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
        ...

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
        ...

    @staticmethod
    def symbolic_x(i):
        """
        Provides a symbolic representation of a position

        :param i:
        :type i:
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def symbolic_n(i, j, k):
        """
        Provides a symbolic representation of a normal to a plane

        :param i:
        :type i:
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def symbolic_m(i):
        """
        Provides a symbolic representation of a mass

        :param i:
        :type i:
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def symbol(base, *args, **kwargs):
        ...

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
        ...

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
        ...

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
        ...

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
        ...

    @classmethod
    def infer_coord_type(cls, inds):
        ...

    @classmethod
    def var(cls, *args, coord_type=None):
        ...

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
        ...

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
        ...

    @classmethod
    def is_identity(cls, A):
        ...

    @classmethod
    def transpose(cls, A):
        ...

    @classmethod
    def dot(cls, a, b, axes=None):
        ...

    @classmethod
    def contract(cls, a, axes):
        ...

    @classmethod
    def transform_coordinates(cls, rotation, coord_vec=None, coord_name_fmt='q{id}[{num}]'):
        ...