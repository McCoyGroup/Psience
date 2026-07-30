"""
Redoes what was originally PyDVR but in the _right_ way using proper subclassing and abstract properties
"""
import abc, numpy as np, scipy.sparse as sp, scipy.interpolate as interp
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger
__all__ = ['BaseDVR', 'DVRResults', 'DVRException']

class BaseDVR(metaclass=abc.ABCMeta):
    """
    Provides the abstract interface for creating a
    convenient runnable DVR that can be cleanly subclassed to provide
    extensions
    """

    def __init__(self, domain=None, divs=None, potential_function=None, logger=None, **base_opts):
        """
        :param base_opts: base opts to use when running
        :type base_opts:
        """
        ...

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation showing the class name, domain, number of grid points, and potential function.

        :return: string of the form `ClassName(domain, pts=divs, pot=potential_function)`
        :rtype: str
        """
        ...

    def _logger(self, logger):
        """
        **LLM Docstring**

        Resolve which logger to use for an operation: the given `logger` if not `None`, otherwise this DVR's own stored logger.

        :param logger: an explicit logger to use, or `None` to fall back to `self.logger`
        :type logger: Logger | None
        :return: the resolved logger
        :rtype: Logger
        """
        ...

    @abc.abstractmethod
    def get_grid(self, domain=None, divs=None, **kwargs):
        """
        **LLM Docstring**

        Abstract hook for building this DVR's 1D grid over the given domain/division count. Concrete DVR subclasses must implement this.

        :param domain: the coordinate domain to build the grid over
        :type domain: tuple | None
        :param divs: the number of grid points
        :type divs: int | None
        :param kwargs: extra representation-specific options
        :type kwargs: dict
        :return: never returns on the base class
        :rtype: np.ndarray
        :raises NotImplementedError: always, on the base class
        """
        ...

    def grid(self, domain=None, divs=None, **kwargs):
        """
        **LLM Docstring**

        Build (or retrieve) this DVR's grid, falling back to the stored `domain`/`divs` if not given explicitly, via `get_grid`.

        :param domain: the coordinate domain to build the grid over; defaults to `self.domain`
        :type domain: tuple | None
        :param divs: the number of grid points; defaults to `self.divs`
        :type divs: int | None
        :param kwargs: extra options forwarded to `get_grid`
        :type kwargs: dict
        :return: the DVR grid
        :rtype: np.ndarray
        :raises ValueError: if neither `domain` nor `self.domain`, or neither `divs` nor `self.divs`, is available
        """
        ...

    @abc.abstractmethod
    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        """
        **LLM Docstring**

        Abstract hook for building this DVR's 1D kinetic-energy operator matrix on the given grid. Concrete DVR subclasses must implement this.

        :param grid: the DVR grid to build the operator on
        :type grid: np.ndarray | None
        :param mass: the particle mass
        :type mass: float | None
        :param hb: the value of hbar to use
        :type hb: float
        :param kwargs: extra representation-specific options
        :type kwargs: dict
        :return: never returns on the base class
        :rtype: np.ndarray
        :raises NotImplementedError: always, on the base class
        """
        ...

    def handle_kinetic_coupling(self, grid, ke_1D, g, g_deriv, hb=1, logger=None, **kwargs):
        """
        **LLM Docstring**

        Apply a (possibly coordinate-dependent) kinetic-coupling correction to a 1D kinetic-energy matrix: multiplies each off-diagonal element by the averaged `g`-value of its two grid points and adds a diagonal correction from `g_deriv`, following the standard variable-mass DVR kinetic-energy formula.

        :param grid: the DVR grid the kinetic-energy matrix is defined on
        :type grid: np.ndarray
        :param ke_1D: the base (unit-`g`) 1D kinetic-energy matrix to correct
        :type ke_1D: np.ndarray
        :param g: the kinetic-coupling function/value; if `None`, `ke_1D` is returned unchanged
        :type g: callable | float | np.ndarray | None
        :param g_deriv: the second-derivative correction term for `g`; required if `g` is given
        :type g_deriv: callable | float | np.ndarray | None
        :param hb: the value of hbar to use
        :type hb: float
        :param logger: logger for diagnostics
        :type logger: Logger | None
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the (possibly `g`-corrected) kinetic-energy matrix
        :rtype: np.ndarray
        :raises ValueError: if `g` is given without a corresponding `g_deriv`
        """
        ...

    def kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs):
        """
        **LLM Docstring**

        Build the full kinetic-energy operator, computing the base 1D operator (via `get_kinetic_energy`) and then applying any kinetic-coupling correction (via `handle_kinetic_coupling`); when a kinetic-coupling function `g` is supplied, the mass is fixed to `1` since `g` already encodes the effective mass.

        :param grid: the DVR grid to build the operator on; defaults to `self.grid()`
        :type grid: np.ndarray | None
        :param mass: the particle mass; required unless `g` is given
        :type mass: float | None
        :param hb: the value of hbar to use
        :type hb: float
        :param g: the kinetic-coupling function/value
        :type g: callable | float | np.ndarray | None
        :param g_deriv: the second-derivative correction term for `g`
        :type g_deriv: callable | float | np.ndarray | None
        :param kwargs: extra options forwarded to `get_kinetic_energy`/`handle_kinetic_coupling`
        :type kwargs: dict
        :return: the kinetic-energy operator matrix
        :rtype: np.ndarray
        :raises ValueError: if no mass is available and `g` isn't given
        """
        ...

    def real_momentum(self, grid=None, mass=None, hb=1, **kwargs):
        """
        **LLM Docstring**

        Abstract hook for the real part of the momentum-operator matrix on the given grid. Not implemented on the base class; concrete DVR subclasses that support it must override this.

        :param grid: the DVR grid to build the operator on
        :type grid: np.ndarray | None
        :param mass: the particle mass
        :type mass: float | None
        :param hb: the value of hbar to use
        :type hb: float
        :param kwargs: extra representation-specific options
        :type kwargs: dict
        :return: never returns on the base class
        :rtype: np.ndarray
        :raises NotImplementedError: always, on the base class
        """
        ...

    def potential_energy(self, grid=None, potential_function=None, potential_values=None, potential_grid=None, logger=None, **pars):
        """
        Calculates the potential energy at the grid points based
        on dispatching on the input form of the potential

        :param grid: the grid of points built earlier in the DVR
        :type grid:
        :param potential_function: a function to evaluate the potential energy at the points
        :type potential_function:
        :param potential_values: the values of the potential at the DVR points
        :type potential_values:
        :param potential_grid: a grid of points and values to be interpolated
        :type potential_grid:
        :param pars: ignored keyword arguments
        :type pars:
        :return:
        :rtype:
        """
        ...

    def hamiltonian(self, kinetic_energy=None, potential_energy=None, potential_threshold=None, **pars):
        """
        Calculates the total Hamiltonian from the kinetic and potential matrices

        :param kinetic_energy:
        :type kinetic_energy:
        :param potential_energy:
        :type potential_energy: np.ndarray | sp.spmatrix
        :param potential_threshold:
        :type potential_threshold:
        :param pars:
        :type pars:
        :return:
        :rtype:
        """
        ...

    def wavefunctions(self, hamiltonian=None, num_wfns=25, nodeless_ground_state=False, diag_mode=None, logger=None, **pars):
        """
        Calculates the wavefunctions for the given Hamiltonian.
        Doesn't support any kind of pruning based on potential values although that might be a good feature
        to support explicitly in the future

        :param hamiltonian:
        :type hamiltonian:
        :param num_wfns:
        :type num_wfns:
        :param nodeless_ground_state:
        :type nodeless_ground_state:
        :param diag_mode:
        :type diag_mode:
        :param pars:
        :type pars:
        :return:
        :rtype:
        """
        ...

    def run(self, result='wavefunctions', logger=None, grid=None, potential_energy=None, kinetic_energy=None, hamiltonian=None, **opts):
        """
        :return:
        :rtype: DVRResults
        """
        ...

class DVRResults:
    """
    A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension
    """

    def __init__(self, grid=None, kinetic_energy=None, potential_energy=None, hamiltonian=None, wavefunctions=None, parent=None, **opts):
        """
        **LLM Docstring**

        Store the full set of intermediate and final results from a DVR run (grid, kinetic/potential-energy operators, Hamiltonian, wavefunctions) alongside the DVR object that produced them and any extra run options.

        :param grid: the DVR grid used
        :type grid: np.ndarray | None
        :param kinetic_energy: the kinetic-energy operator matrix
        :type kinetic_energy: np.ndarray | None
        :param potential_energy: the potential-energy operator matrix
        :type potential_energy: np.ndarray | None
        :param hamiltonian: the full Hamiltonian matrix
        :type hamiltonian: np.ndarray | None
        :param wavefunctions: the resulting wavefunctions
        :type wavefunctions: DVRWavefunctions | None
        :param parent: the `BaseDVR` object that produced these results
        :type parent: BaseDVR | None
        :param opts: extra run options/metadata to store
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def dimension(self):
        """
        **LLM Docstring**

        The number of spatial dimensions of the underlying DVR grid.

        :return: the dimensionality
        :rtype: int
        """
        ...

    def plot_potential(self, plot_class=None, figure=None, plot_units=None, energy_threshold=None, zero_shift=False, **opts):
        """
        Simple plotting function for the potential.
        Should be updated to deal with higher dimensional cases

        :param plot_class: the graphics class to use for the plot
        :type plot_class: McUtils.Plots.Graphics
        :param opts: plot styling options
        :type opts:
        :return:
        :rtype: McUtils.Plots.Graphics
        """
        ...

class DVRException(Exception):
    """
    Base exception class for working with DVRs
    """