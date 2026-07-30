"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""
import numpy as np
from McUtils.Plots import Graphics
from McUtils.Zachary import Interpolator
from Psience.Wavefun import Wavefunction, Wavefunctions
from .BaseDVR import DVRResults
__all__ = ['DVRWavefunctions', 'DVRWavefunction']

class DVRWavefunction(Wavefunction):

    def __init__(self, energy, data, parent=None, grid=None, index=None, **opts):
        """
        **LLM Docstring**

        Build a single DVR wavefunction from its value at the DVR grid points, falling back to the parent collection's grid if none is given directly.

        :param energy: the wavefunction's energy
        :type energy: float
        :param data: the wavefunction's values at the grid points
        :type data: np.ndarray
        :param parent: the `DVRWavefunctions` collection this wavefunction belongs to
        :type parent: DVRWavefunctions | None
        :param grid: the DVR grid this wavefunction is defined on; defaults to `parent.grid`
        :type grid: np.ndarray | None
        :param index: this wavefunction's index within its parent collection
        :type index: object | None
        :param opts: extra options forwarded to the base `Wavefunction.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def get_dimension(self):
        """
        **LLM Docstring**

        The number of degrees of freedom this wavefunction is defined over, inferred from the trailing dimension of its grid.

        :return: the dimensionality
        :rtype: int
        """
        ...

    def plot(self, figure=None, grid=None, **opts):
        """
        **LLM Docstring**

        Plot the wavefunction using its own DVR grid and values, delegating to the base `Wavefunction.plot`.

        :param figure: an existing figure to draw into
        :type figure: object | None
        :param grid: the grid to plot over; defaults to `self.grid`
        :type grid: np.ndarray | None
        :param opts: extra options forwarded to the base `plot`
        :type opts: dict
        :return: the resulting figure
        :rtype: object
        """
        ...

    def expectation(self, op, other=None):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """
        ...

    @property
    def interp(self):
        """
        **LLM Docstring**

        A (lazily built and cached) continuous interpolant of the wavefunction's grid values, used by `evaluate` to evaluate the wavefunction off-grid.

        :return: the cached (or newly built) interpolator
        :rtype: Interpolator
        """
        ...

    def evaluate(self, points):
        """
        Evaluates the functions at the given points

        :return:
        :rtype:
        """
        ...

    def marginalize_out(self, dofs):
        """
        Computes the projection of the current wavefunction onto a set of degrees
        of freedom

        :return:
        :rtype:
        """
        ...

class DVRWavefunctions(Wavefunctions):
    wavefunction_class = DVRWavefunction

    def __init__(self, energies=None, wavefunctions=None, grid=None, results: DVRResults=None, **opts):
        """
        **LLM Docstring**

        Build a collection of DVR wavefunctions sharing a common grid and the `DVRResults` object they were solved from.

        :param energies: the energies of each wavefunction in the collection
        :type energies: np.ndarray | None
        :param wavefunctions: the matrix of wavefunction values at the grid points, one column per state
        :type wavefunctions: np.ndarray | None
        :param grid: the DVR grid the wavefunctions are defined on
        :type grid: np.ndarray | None
        :param results: the `DVRResults` object (grid, kinetic energy, potential energy, etc.) the wavefunctions were solved from
        :type results: DVRResults | None
        :param opts: extra options forwarded to the base `Wavefunctions.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation showing the class name, the number of wavefunctions, and the DVR object they were solved from.

        :return: string of the form `ClassName(num=N, DVR=dvr)`
        :rtype: str
        """
        ...

    def plot(self, figure=None, **opts):
        """
        Plots the held wavefunctions

        :param figure:
        :type figure:
        :param graphics_class:
        :type graphics_class:
        :param plot_style:
        :type plot_style:
        :param scaling:
        :type scaling:
        :param shift:
        :type shift:
        :param opts:
        :type opts:
        :return:
        :rtype: Graphics
        """
        ...

    def expectation(self, op, other=None, multiplicative=True):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: DVRWavefunctions | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """
        ...

    def transform_operator(self, M):
        """
        **LLM Docstring**

        Transform an operator matrix given in the DVR grid-point basis into the basis of these wavefunctions, by sandwiching it between the wavefunction coefficient matrix and its transpose.

        :param M: the operator matrix, in the DVR grid basis (dense or sparse)
        :type M: np.ndarray | scipy.sparse.spmatrix
        :return: the operator matrix in the wavefunction basis
        :rtype: np.ndarray
        """
        ...

    def coordinate(self):
        """
        **LLM Docstring**

        The position-operator matrix in the wavefunction basis, computed as the expectation value of the DVR grid points themselves.

        :return: the coordinate-operator matrix
        :rtype: np.ndarray
        """
        ...

    def momentum(self):
        """
        **LLM Docstring**

        The real part of the momentum-operator matrix in the wavefunction basis, computed from the underlying DVR's `real_momentum` operator and transformed into the wavefunction basis.

        :return: the momentum-operator matrix
        :rtype: np.ndarray
        """
        ...

    def laplacian(self):
        """
        **LLM Docstring**

        The Laplacian operator matrix in the wavefunction basis, derived from a fresh (unit-mass, zero-potential, uncoupled) kinetic-energy calculation on the underlying DVR and transformed into the wavefunction basis.

        :return: the Laplacian operator matrix
        :rtype: np.ndarray
        """
        ...

    def kinetic_energy(self):
        """
        **LLM Docstring**

        The kinetic-energy operator matrix in the wavefunction basis, transformed from the stored `DVRResults.kinetic_energy` (in the grid basis).

        :return: the kinetic-energy operator matrix
        :rtype: np.ndarray
        """
        ...

    def potential_energy(self):
        """
        **LLM Docstring**

        The potential-energy operator matrix in the wavefunction basis, transformed from the stored `DVRResults.kinetic_energy` (in the grid basis).

        :return: the (mistakenly computed) kinetic-energy operator matrix
        :rtype: np.ndarray
        """
        ...