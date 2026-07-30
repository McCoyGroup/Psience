"""
Provides a direct-product extension to a system of 1D DVRs
"""
import numpy as np, scipy.sparse as sp
from functools import reduce
from .BaseDVR import BaseDVR
from .ColbertMiller import *
__all__ = ['DirectProductDVR', 'CartesianNDDVR', 'RingNDDVR', 'SphericalDVR']

class DirectProductDVR(BaseDVR):

    def __init__(self, dvrs_1D, zero_threshold=1e-14, **base_opts):
        """
        :param dvrs_1D: a series of 1D DVRs that can provide the inputs we'll product together
        :type dvrs_1D: Iterable[AbstractDVR]
        :param base_opts:
        :type base_opts:
        """
        ...

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation showing the class name, the component 1D DVRs, and the potential function.

        :return: string of the form `ClassName(dvrs, pot=potential_function)`
        :rtype: str
        """
        ...

    def get_grid(self, domain=None, divs=None, **kwargs):
        """
        **LLM Docstring**

        Build the full multi-dimensional DVR grid as the Cartesian product of each component 1D DVR's own grid, handling component DVRs (like `FiniteBasisDVR`) that additionally return a basis transformation alongside their grid points.

        :param domain: per-dimension domains to build fresh 1D grids with, instead of each component DVR's own stored grid
        :type domain: list[tuple] | None
        :param divs: per-dimension division counts, paired with `domain`
        :type divs: list[int] | None
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the multi-dimensional grid, or `(grid, transformations)` if any component DVR returned a basis transformation
        :rtype: np.ndarray | tuple[np.ndarray, list]
        """
        ...

    def grid(self, domain=None, divs=None, **kwargs):
        """
        **LLM Docstring**

        Build (or retrieve) the multi-dimensional DVR grid, falling back to this DVR's stored `domain`/`divs` if not given explicitly, via `get_grid`.

        :param domain: per-dimension domains; defaults to `self.domain`
        :type domain: list[tuple] | None
        :param divs: per-dimension division counts; defaults to `self.divs`
        :type divs: list[int] | None
        :param kwargs: extra options forwarded to `get_grid`
        :type kwargs: dict
        :return: the multi-dimensional grid, or `(grid, transformations)`
        :rtype: np.ndarray | tuple
        """
        ...

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, include_kinetic_coupling=True, **kwargs):
        """
        **LLM Docstring**

        Assemble the full multi-dimensional kinetic-energy operator as a sparse matrix, either as a simple Kronecker sum of the per-dimension 1D kinetic-energy operators (constant-mass case), or, when a (possibly coordinate-dependent) kinetic-coupling tensor `g`/`g_deriv` is supplied, by explicitly building the diagonal G-matrix-weighted kinetic terms, the pseudopotential-like `g_deriv` correction, and the off-diagonal momentum-momentum coupling terms between each pair of dimensions.

        :param grid: the multi-dimensional grid (or `(grid, transformations)` pair) to build the operator on; each dimension's subgrid is extracted from it
        :type grid: np.ndarray | tuple
        :param mass: the mass (or per-dimension masses) to use if `g` isn't given
        :type mass: float | list[float] | None
        :param hb: the value of hbar (or per-dimension values) to use
        :type hb: float | list[float]
        :param g: the (possibly coordinate-dependent) kinetic-coupling tensor, `g[i][j]` giving the coupling between dimensions `i` and `j` as a constant or a callable of the flattened grid
        :type g: list[list] | None
        :param g_deriv: the second-derivative correction terms for the diagonal `g` entries, per dimension
        :type g_deriv: list | None
        :param logger: logger used to report progress when kinetic coupling is being evaluated
        :type logger: Logger | None
        :param include_kinetic_coupling: whether to include the off-diagonal momentum-coupling terms (only relevant when `g` is given and has nonzero off-diagonal entries)
        :type include_kinetic_coupling: bool
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the assembled kinetic-energy matrix, as a sparse matrix (or dense array if it ends up sufficiently non-sparse)
        :rtype: scipy.sparse.spmatrix | np.ndarray
        :raises ValueError: if `g` is given without a corresponding `g_deriv`
        """
        ...

    def kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, **kwargs):
        """
        Computes the N-dimensional kinetic energy
        :param grid:
        :type grid:
        :param mass:
        :type mass:
        :param hb:
        :type hb:
        :param g:
        :type g:
        :param g_deriv:
        :type g_deriv:
        :return:
        :rtype:
        """
        ...

class CartesianNDDVR(DirectProductDVR):
    """
    Provides an ND-DVR over different domains
    """

    def __init__(self, domains, **base_opts):
        """
        **LLM Docstring**

        Set up an N-dimensional Cartesian DVR by building an independent `CartesianDVR` for each requested `(min, max, divs)` domain specification and combining them via `DirectProductDVR`.

        :param domains: the per-dimension `(min, max, divs)` domain specifications
        :type domains: Iterable[tuple]
        :param base_opts: extra options forwarded to the base `DirectProductDVR.__init__`
        :type base_opts: dict
        :return: None
        :rtype: None
        """
        ...

class RingNDDVR(DirectProductDVR):
    """
    Provides an ND-DVR for products of periodic (0, 2Pi) ranges
    """

    def __init__(self, divs, **base_opts):
        """
        **LLM Docstring**

        Set up an N-dimensional DVR over a product of periodic `(0, 2*pi)` ring domains by building an independent `RingDVR` for each requested division count and combining them via `DirectProductDVR`.

        :param divs: the per-dimension number of grid points
        :type divs: Iterable[int]
        :param base_opts: extra options forwarded to the base `DirectProductDVR.__init__`
        :type base_opts: dict
        :return: None
        :rtype: None
        """
        ...

class SphericalDVR(DirectProductDVR):

    def __init__(self, r_max, divs, **base_opts):
        """
        **LLM Docstring**

        Set up a 3D spherical-coordinate DVR by combining a `RadialDVR` (over `(0, r_max)`), a `RingDVR`, and a `PolarDVR` via `DirectProductDVR`. Note the `PolarDVR` component is constructed with `domain=(0, 2*np.pi)`, matching the `RingDVR` above it rather than the `(0, pi)` domain `PolarDVR` normally defaults to and is designed for -- this looks like a copy-paste bug.

        :param r_max: the maximum radial extent
        :type r_max: float
        :param divs: the `(radial_divs, angular_divs)` number of grid points for the radial and angular dimensions (the same `angular_divs` is used for both the ring and polar dimensions)
        :type divs: tuple[int, int]
        :param base_opts: extra options forwarded to the base `DirectProductDVR.__init__`
        :type base_opts: dict
        :return: None
        :rtype: None
        """
        ...