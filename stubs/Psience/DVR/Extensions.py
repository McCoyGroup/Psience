import functools
import numpy as np
from ..VSCF import GridSCF
from McUtils.Scaffolding import ParameterManager
from .BaseDVR import BaseDVR
from .Wavefunctions import DVRWavefunctions
from .DirectProduct import DirectProductDVR
from .FiniteBasisDVR import WavefunctionBasisDVR
__all__ = ['SelfConsistentDVR', 'PotentialOptimizedDVR']

class SCFWavefunctionGenerator:

    def __init__(self, dvr_1D: BaseDVR):
        """
        **LLM Docstring**

        Wrap a 1D DVR so it can be re-solved repeatedly with updated potential values (reusing the cached grid/kinetic-energy operator after the first solve), for use as a per-dimension wavefunction generator in a self-consistent-field iteration.

        :param dvr_1D: the 1D DVR to wrap
        :type dvr_1D: BaseDVR
        :return: None
        :rtype: None
        """
        ...

    def __call__(self, pot, **kwargs):
        """
        **LLM Docstring**

        Solve the wrapped 1D DVR for the given potential values, reusing the cached grid and kinetic-energy operator from the previous call (if any) to avoid recomputing them each SCF iteration.

        :param pot: the potential values (on this dimension's DVR grid) to solve with
        :type pot: np.ndarray
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the resulting 1D wavefunctions
        :rtype: DVRWavefunctions
        """
        ...

class SelfConsistentDVR(GridSCF):

    def __init__(self, base_dvr: 'DirectProductDVR', **opts):
        """
        **LLM Docstring**

        Set up a self-consistent-field (SCF) treatment of a multi-dimensional `DirectProductDVR`, precomputing the grid and potential-energy values and building the per-dimension SCF wavefunction generators before delegating to `GridSCF.__init__`.

        :param base_dvr: the multi-dimensional direct-product DVR to run the SCF procedure on
        :type base_dvr: DirectProductDVR
        :param opts: extra options, filtered and forwarded to the base `GridSCF.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def create_grid_vals(self):
        """
        **LLM Docstring**

        Compute the full grid and potential-energy values for the underlying multi-dimensional DVR, used to initialize the SCF procedure.

        :return: `(grid, pe)` -- the DVR grid and the potential-energy values reshaped to match the grid's spatial shape
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def create_solvers(self, grid, pe):
        """
        **LLM Docstring**

        Build the per-dimension SCF wavefunction generators, first rebinding each 1D DVR's `g`/`g_deriv` kinetic-coupling functions (if callable) so they can be evaluated at the fixed "other-dimension" grid point (from the initial guess) while varying only their own coordinate.

        :param grid: the full multi-dimensional grid
        :type grid: np.ndarray
        :param pe: the potential-energy values on the grid
        :type pe: np.ndarray
        :return: the list of per-dimension `SCFWavefunctionGenerator` objects
        :rtype: list[SCFWavefunctionGenerator]
        """
        ...

    @staticmethod
    def _wrap_g(g, gp, i):
        """
        **LLM Docstring**

        Wrap a per-dimension kinetic-coupling function `g` so that, when evaluated during the SCF procedure for dimension `i`, the other dimensions' coordinates are held fixed at the reference point `gp` while only dimension `i`'s coordinate is allowed to vary over the passed-in 1D grid.

        :param g: the kinetic-coupling function to wrap; passed through unchanged if not callable
        :type g: callable | object
        :param gp: the fixed reference point (coordinates for every dimension) to hold the other dimensions at
        :type gp: np.ndarray
        :param i: the dimension index that should vary
        :type i: int
        :return: the wrapped function, or `g` unchanged if it wasn't callable
        :rtype: callable | object
        """
        ...

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation showing the class name and the wrapped base DVR.

        :return: string of the form `ClassName(base_dvr)`
        :rtype: str
        """
        ...

class PotentialOptimizedDVR(DirectProductDVR):

    def __init__(self, wfns_1D: 'Iterable[DVRWavefunctions]', **base_opts):
        """
        **LLM Docstring**

        Build a multi-dimensional DVR whose per-dimension bases are potential-optimized (built from a given set of 1D wavefunctions rather than a fixed grid), via `WavefunctionBasisDVR` and `DirectProductDVR`.

        :param wfns_1D: the per-dimension 1D wavefunctions to use as the optimized basis
        :type wfns_1D: Iterable[DVRWavefunctions]
        :param base_opts: extra options forwarded to the base `DirectProductDVR.__init__`
        :type base_opts: dict
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def from_minimum(cls, base_dvr: 'DirectProductDVR|SelfConsistentDVR', **opts):
        """
        **LLM Docstring**

        Build a `PotentialOptimizedDVR` using the SCF wavefunctions computed from an initial (unconverged, single-iteration) guess at the potential minimum, wrapping `base_dvr` in a `SelfConsistentDVR` first if it isn't already one.

        :param base_dvr: the base multi-dimensional DVR (or an already-built `SelfConsistentDVR`) to optimize the basis from
        :type base_dvr: DirectProductDVR | SelfConsistentDVR
        :param opts: extra options, overriding the base DVR's own stored options, forwarded to the constructor
        :type opts: dict
        :return: the potential-optimized DVR
        :rtype: PotentialOptimizedDVR
        """
        ...

    @classmethod
    def from_scf(cls, scf_dvr: 'DirectProductDVR|SelfConsistentDVR', wfns=None, **opts):
        """
        **LLM Docstring**

        Build a `PotentialOptimizedDVR` using the wavefunctions from a fully-converged SCF run (or an explicitly supplied set of wavefunctions), wrapping `scf_dvr` in a `SelfConsistentDVR` first if it isn't already one.

        :param scf_dvr: the base multi-dimensional DVR (or an already-built `SelfConsistentDVR`) to run/use for the optimized basis
        :type scf_dvr: DirectProductDVR | SelfConsistentDVR
        :param wfns: explicit wavefunctions to use instead of running the SCF procedure
        :type wfns: Iterable[DVRWavefunctions] | None
        :param opts: extra options, overriding the base DVR's own stored options, forwarded to the constructor
        :type opts: dict
        :return: the potential-optimized DVR
        :rtype: PotentialOptimizedDVR
        """
        ...