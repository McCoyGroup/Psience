import numpy as np
from .BaseDVR import BaseDVR
import typing

class InitialBasis(typing.Protocol):

    @property
    def dimensions(self) -> 'Iterable[int]':
        """
        **LLM Docstring**

        Protocol property for the number of basis functions along each dimension of the underlying finite basis. Not implemented; concrete basis types satisfying this protocol must provide their own.

        :return: the per-dimension basis sizes
        :rtype: Iterable[int]
        """
        ...

    def x(self, n: int) -> np.ndarray:
        """
        **LLM Docstring**

        Protocol stub for the position-operator matrix representation in the first `n` basis functions. Not implemented; concrete basis types must provide their own.

        :param n: the number of basis functions to represent the operator in
        :type n: int
        :return: the position-operator matrix
        :rtype: np.ndarray
        """
        ...

    def p2(self, n: int) -> np.ndarray:
        """
        **LLM Docstring**

        Protocol stub for the momentum-squared operator matrix representation in the first `n` basis functions. Not implemented; concrete basis types must provide their own.

        :param n: the number of basis functions to represent the operator in
        :type n: int
        :return: the momentum-squared operator matrix
        :rtype: np.ndarray
        """
        ...

    def p(self, n: int) -> np.ndarray:
        """
        **LLM Docstring**

        Protocol stub for the momentum-operator matrix representation in the first `n` basis functions. Not implemented; concrete basis types must provide their own.

        :param n: the number of basis functions to represent the operator in
        :type n: int
        :return: the momentum-operator matrix
        :rtype: np.ndarray
        """
        ...

class FiniteBasisDVR(BaseDVR):

    def __init__(self, basis: InitialBasis, domain=None, divs=None, **opts):
        """
        **LLM Docstring**

        Set up a DVR built from an arbitrary finite (non-grid) initial basis, using the basis's own dimension count as the number of DVR points.

        :param basis: the initial finite basis to build the DVR from
        :type basis: InitialBasis
        :param domain: accepted for interface consistency with `BaseDVR` but not used (always passed through as `None`)
        :type domain: tuple | None
        :param divs: accepted for interface consistency but overridden by `basis.dimensions[0]`
        :type divs: int | None
        :param opts: extra options forwarded to the base `BaseDVR.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def get_grid(self, domain=None, divs=None, **kwargs):
        """
        **LLM Docstring**

        Build the DVR grid by diagonalizing the position-operator matrix of the underlying finite basis.

        :param domain: accepted for interface consistency but not used
        :type domain: tuple | None
        :param divs: accepted for interface consistency but not used (the basis's own dimension count is used instead)
        :type divs: int | None
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: `(grid_points, transformation)` -- the position eigenvalues and the eigenvector matrix transforming from the initial basis to the DVR grid basis
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def real_momentum(self, grid=None, mass=None, hb=1, **kwargs):
        """
        **LLM Docstring**

        Build the momentum-operator matrix in the DVR grid basis, by transforming the initial basis's momentum-operator matrix through the position-eigenvector transformation.

        :param grid: the `(grid_points, transformation)` pair from `get_grid`
        :type grid: tuple
        :param mass: accepted for interface consistency but not used
        :type mass: float | None
        :param hb: accepted for interface consistency but not used
        :type hb: float
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the momentum-operator matrix in the DVR basis
        :rtype: np.ndarray
        """
        ...

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        """
        **LLM Docstring**

        Build the kinetic-energy matrix in the DVR grid basis, using the initial basis's momentum-squared operator (falling back to squaring its momentum-operator matrix if `p2` isn't implemented) transformed through the position-eigenvector transformation.

        :param grid: the `(grid_points, transformation)` pair from `get_grid`
        :type grid: tuple
        :param mass: the particle mass
        :type mass: float
        :param hb: the value of hbar to use
        :type hb: float
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the kinetic-energy matrix in the DVR basis
        :rtype: np.ndarray
        """
        ...

    def potential_energy(self, grid=None, potential_function=None, potential_values=None, potential_grid=None, logger=None, **pars):
        """
        **LLM Docstring**

        Compute the potential-energy matrix using just the grid points (rather than the full `(grid_points, transformation)` pair), delegating to the base `BaseDVR.potential_energy`.

        :param grid: the `(grid_points, transformation)` pair from `get_grid`, of which only the grid points are used
        :type grid: tuple
        :param potential_function: the potential-energy function to evaluate on the grid
        :type potential_function: callable | None
        :param potential_values: precomputed potential values at the grid points
        :type potential_values: np.ndarray | None
        :param potential_grid: an explicit `(points, values)` pair to interpolate the potential from
        :type potential_grid: tuple | None
        :param logger: logger for diagnostics
        :type logger: Logger | None
        :param pars: extra options forwarded to the base `potential_energy`
        :type pars: dict
        :return: the potential-energy matrix
        :rtype: np.ndarray
        """
        ...

class HarmonicDVR(FiniteBasisDVR):

    def __init__(self, divs=None, **opts):
        """
        **LLM Docstring**

        Set up a DVR built from a harmonic-oscillator basis of the requested size.

        :param divs: the number of harmonic-oscillator basis functions to use
        :type divs: int | None
        :param opts: extra options forwarded to `FiniteBasisDVR.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

class WavefunctionBasisDVR(FiniteBasisDVR):

    def __init__(self, wfns=None, **opts):
        """
        **LLM Docstring**

        Set up a DVR built from an existing set of wavefunctions used as the initial basis.

        :param wfns: the wavefunctions to use as the initial basis
        :type wfns: object | None
        :param opts: extra options forwarded to `FiniteBasisDVR.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...