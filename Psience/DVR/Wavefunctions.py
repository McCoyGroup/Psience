"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""

import numpy as np

from McUtils.Plots import Graphics
from McUtils.Zachary import Interpolator

from Psience.Wavefun import Wavefunction, Wavefunctions
from .BaseDVR import DVRResults

__all__ = ["DVRWavefunctions", "DVRWavefunction"]

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
        super().__init__(energy, data, parent=parent, index=index, **opts)
        if grid is None:
            grid = self.parent.grid
        self.grid = grid
        self._interp = None
    def get_dimension(self):
        """
        **LLM Docstring**

        The number of degrees of freedom this wavefunction is defined over, inferred from the trailing dimension of its grid.

        :return: the dimensionality
        :rtype: int
        """
        return self.grid.shape[-1]

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
        if grid is None:
            grid = self.grid
        return super().plot(figure=figure, grid=grid, values=self.data, **opts)

    def expectation(self, op, other=None):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """
        wf = op(self.data)
        if not isinstance(other, np.ndarray):
            other = other.data
        return np.dot(wf, other)

    @property
    def interp(self):
        """
        **LLM Docstring**

        A (lazily built and cached) continuous interpolant of the wavefunction's grid values, used by `evaluate` to evaluate the wavefunction off-grid.

        :return: the cached (or newly built) interpolator
        :rtype: Interpolator
        """
        if self._interp is None:
            self._interp = Interpolator(self.grid, self.data)
        return self._interp
    def evaluate(self, points):
        """
        Evaluates the functions at the given points

        :return:
        :rtype:
        """
        return self.interp(points)

    def marginalize_out(self, dofs):
        """
        Computes the projection of the current wavefunction onto a set of degrees
        of freedom

        :return:
        :rtype:
        """

        # if isinstance(dofs, (int, np.integer)):
        #     dofs = [dofs]
        # dofs = np.flip(np.sort(dofs))

        raise NotImplementedError("DVR projections not yet implemented")

class DVRWavefunctions(Wavefunctions):
    # most evaluations are most efficient done in batch for DVR wavefunctions so we focus on the batch object
    wavefunction_class = DVRWavefunction
    def __init__(self, energies=None, wavefunctions=None, grid=None, results:DVRResults=None, **opts):
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
        super().__init__(energies=energies, wavefunctions=wavefunctions, results=results, grid=grid, **opts) # add all opts
        self.results = results
        self.grid = grid
    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation showing the class name, the number of wavefunctions, and the DVR object they were solved from.

        :return: string of the form `ClassName(num=N, DVR=dvr)`
        :rtype: str
        """
        return "{}(num={}, DVR={})".format(
            type(self).__name__,
            len(self),
            self.results.parent
        )

    # def __len__(self):
    #     return len(self.energies)
    # def __iter__(self):
    #     for i in range(len(self)):
    #         yield self.__getitem__(i)
    # def __getitem__(self, item):
    #     """
    #     Provides a single `DVRWavefunction` or slice of `DVRWavefunctions`
    #     :param item:
    #     :type item:
    #     :return:
    #     :rtype: DVRWavefunction | DVRWavefunctions
    #     """
    #     if not isinstance(item, (int, np.integer)):
    #         return type(self)(
    #             energies=self.energies[item],
    #             wavefunctions=self.wavefunctions[:, item],#.reshape((len(self.wavefunctions), -1)), # WTF? Was I transposing here?...?
    #             wavefunction_class=self.wavefunction_class,
    #             results=self.results,
    #             grid=self.grid,
    #             **self.opts
    #         )
    #     else:
    #         return self.wavefunction_class(
    #             self.energies[item],
    #             self.wavefunctions[:, item],
    #             parent=self,
    #             **self.opts
    #         )

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

        grid = self.grid

        dim = len(grid.shape) # mesh grid...
        if dim > 1 and grid.shape[-1] == dim-1: # check whether we have a mesh of points that we need to reshape
            grid = np.moveaxis(grid, grid.ndim, 0)

        return super().plot(
            figure=figure,
            grid=grid,
            **opts
        )

    def expectation(self, op, other=None, multiplicative=True):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: DVRWavefunctions | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """
        if other is None:
            other = self
        if not multiplicative:
            raise ValueError("don't have non-multiplicative operators supported yet...")
        else:
            if not isinstance(op, np.ndarray):
                op = op(self.grid.reshape((-1,) + self.grid.shape[2:]))
            wfs = self.wavefunctions
            for _ in range(op.ndim-1):
                wfs = np.expand_dims(wfs, -1)
            wfs = np.expand_dims(op, 1) * wfs
            # print(self.wavefunctions.shape, wfs.shape)
        if not isinstance(other, np.ndarray):
            other = other.wavefunctions
        ev = np.tensordot(other, wfs, axes=[0, 0])
        ev = ev.transpose([1, 0] + list(range(2, ev.ndim)))
        # print("--->", wfs.shape, other.shape, ev.shape)
        return ev

    def transform_operator(self, M):
        """
        **LLM Docstring**

        Transform an operator matrix given in the DVR grid-point basis into the basis of these wavefunctions, by sandwiching it between the wavefunction coefficient matrix and its transpose.

        :param M: the operator matrix, in the DVR grid basis (dense or sparse)
        :type M: np.ndarray | scipy.sparse.spmatrix
        :return: the operator matrix in the wavefunction basis
        :rtype: np.ndarray
        """
        if hasattr(M, 'toarray'):
            M = M.toarray()
        return np.dot(np.dot(self.wavefunctions.T, M), self.wavefunctions)

    def coordinate(self):
        """
        **LLM Docstring**

        The position-operator matrix in the wavefunction basis, computed as the expectation value of the DVR grid points themselves.

        :return: the coordinate-operator matrix
        :rtype: np.ndarray
        """
        return self.expectation(self.results.grid)
    def momentum(self):
        """
        **LLM Docstring**

        The real part of the momentum-operator matrix in the wavefunction basis, computed from the underlying DVR's `real_momentum` operator and transformed into the wavefunction basis.

        :return: the momentum-operator matrix
        :rtype: np.ndarray
        """
        dvr = self.results.parent
        p = dvr.real_momentum(grid=self.results.grid, **dvr.opts)
        return self.transform_operator(p)
    def laplacian(self):
        """
        **LLM Docstring**

        The Laplacian operator matrix in the wavefunction basis, derived from a fresh (unit-mass, zero-potential, uncoupled) kinetic-energy calculation on the underlying DVR and transformed into the wavefunction basis.

        :return: the Laplacian operator matrix
        :rtype: np.ndarray
        """
        dvr = self.results.parent
        res = dvr.run(mass=1, g=None, hb=1, potential_function=lambda g:np.zeros(len(g)), result='kinetic_energy')
        p2 = -2*res.kinetic_energy
        return self.transform_operator(p2)
    def kinetic_energy(self):
        """
        **LLM Docstring**

        The kinetic-energy operator matrix in the wavefunction basis, transformed from the stored `DVRResults.kinetic_energy` (in the grid basis).

        :return: the kinetic-energy operator matrix
        :rtype: np.ndarray
        """
        # import McUtils.Plots as plt
        # plt.ArrayPlot(self.results.kinetic_energy)
        # plt.ArrayPlot(self.transform_operator(self.results.kinetic_energy)).show()
        # print(self.results.kinetic_energy)
        return self.transform_operator(self.results.kinetic_energy)
    def potential_energy(self):
        """
        **LLM Docstring**

        The potential-energy operator matrix in the wavefunction basis, transformed from the stored `DVRResults.kinetic_energy` (in the grid basis).

        :return: the (mistakenly computed) kinetic-energy operator matrix
        :rtype: np.ndarray
        """
        return self.transform_operator(self.results.potential_energy)