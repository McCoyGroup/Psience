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
        super().__init__(energy, data, parent=parent, index=index, **opts)
        if grid is None:
            grid = self.parent.grid
        self.grid = grid
        self._interp = None
    def get_dimension(self):
        return self.grid.shape[-1]

    def plot(self, figure=None, grid=None, **opts):
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
        super().__init__(energies=energies, wavefunctions=wavefunctions, results=results, grid=grid, **opts) # add all opts
        self.results = results
        self.grid = grid
    def __repr__(self):
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
        if hasattr(M, 'toarray'):
            M = M.toarray()
        return np.dot(np.dot(self.wavefunctions.T, M), self.wavefunctions)

    def coordinate(self):
        return self.expectation(self.results.grid)
    def momentum(self):
        dvr = self.results.parent
        p = dvr.real_momentum(grid=self.results.grid, **dvr.opts)
        return self.transform_operator(p)
    def laplacian(self):
        dvr = self.results.parent
        res = dvr.run(mass=1, g=None, hb=1, potential_function=lambda g:np.zeros(len(g)), result='kinetic_energy')
        p2 = -2*res.kinetic_energy
        return self.transform_operator(p2)
    def kinetic_energy(self):
        # import McUtils.Plots as plt
        # plt.ArrayPlot(self.results.kinetic_energy)
        # plt.ArrayPlot(self.transform_operator(self.results.kinetic_energy)).show()
        # print(self.results.kinetic_energy)
        return self.transform_operator(self.results.kinetic_energy)
    def potential_energy(self):
        return self.transform_operator(self.results.kinetic_energy)