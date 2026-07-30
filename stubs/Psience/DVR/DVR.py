import numpy as np
from McUtils.Scaffolding import ParameterManager
from .ColbertMiller import PolarDVR, RingDVR, CartesianDVR
from .DirectProduct import DirectProductDVR
from .Extensions import SelfConsistentDVR, PotentialOptimizedDVR
__all__ = ['DVR']

class DVRConstructor:
    _domain_map = None

    @classmethod
    def load_domain_map(cls):
        """
        **LLM Docstring**

        Build the default mapping from special-cased coordinate domains to the appropriate specialized 1D DVR class: `(0, pi)` for a polar/angular coordinate, `(0, 2*pi)` for a periodic (ring) coordinate, and `None` as the catch-all default for ordinary Cartesian-like coordinates.

        :return: the domain-to-DVR-class mapping
        :rtype: dict
        """
        ...

    @classmethod
    def infer_DVR_type(cls, domain):
        """
        **LLM Docstring**

        Infer which specialized 1D DVR class to use for a given coordinate domain, matching it (via `np.allclose`) against the special-cased domains in the (lazily loaded and cached) domain map, and falling back to the default (Cartesian) class if none match.

        :param domain: the `(min, max)` coordinate domain to infer a DVR class for
        :type domain: tuple
        :return: the inferred DVR class
        :rtype: type
        """
        ...

    @classmethod
    def construct(cls, domain=None, divs=None, potential_function=None, g=None, g_deriv=None, mass=None, po_divs=25, classes=None, scf=False, potential_optimize=False, logger=True, **base_opts):
        """
        **LLM Docstring**

        Dispatch to build an appropriate DVR object from a domain/division specification, either a single 1D DVR (of an inferred or explicitly given class) or, for multiple dimensions, a `DirectProductDVR` combining per-dimension 1D DVRs -- optionally wrapped in a `SelfConsistentDVR` and/or `PotentialOptimizedDVR` if `scf`/`potential_optimize` are requested.

        :param domain: the coordinate domain(s): a single `(min, max)` pair for 1D, or a list of such pairs for multiple dimensions
        :type domain: tuple | list[tuple]
        :param divs: the number of grid points per dimension
        :type divs: int | list[int]
        :param potential_function: the potential-energy function to evaluate on the grid
        :type potential_function: callable | None
        :param g: the kinetic-coupling (G-matrix) function/matrix; if callable, treated as a single-dimension `g` function, otherwise its diagonal is extracted per dimension
        :type g: callable | np.ndarray | None
        :param g_deriv: the derivative of `g` with respect to each coordinate
        :type g_deriv: callable | list | None
        :param mass: the mass (or per-dimension masses) to use if `g` isn't given
        :type mass: float | list[float] | None
        :param po_divs: the number of grid points to use for the potential-optimized basis, per dimension
        :type po_divs: int | list[int]
        :param classes: explicit DVR class(es) to use per dimension, bypassing automatic inference
        :type classes: type | list[type] | None
        :param scf: whether to wrap the resulting multi-dimensional DVR in a `SelfConsistentDVR`
        :type scf: bool
        :param potential_optimize: whether to wrap the resulting DVR in a `PotentialOptimizedDVR` (built from the minimum directly, or from the SCF result if `scf` is also set)
        :type potential_optimize: bool
        :param logger: logger (or logging flag) to pass through to the constructed DVR(s)
        :type logger: object | bool
        :param base_opts: extra options forwarded to the underlying DVR constructor(s), filtered per target class where relevant
        :type base_opts: dict
        :return: the constructed DVR object
        :rtype: object
        :raises ValueError: if `domain` or `divs` is `None`
        """
        ...

def DVR(domain=None, divs=None, classes=None, potential_function=None, g=None, g_deriv=None, scf=False, potential_optimize=False, **base_opts):
    """
    Constructs a DVR object

    :param domain:
    :type domain:
    :param divs:
    :type divs:
    :param classes:
    :type classes:
    :param potential_function:
    :type potential_function:
    :param g:
    :type g:
    :param g_deriv:
    :type g_deriv:
    :param base_opts:
    :type base_opts:
    :return:
    :rtype:
    """
    ...