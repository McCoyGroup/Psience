"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""

import numpy as np

from McUtils.Plots import Graphics, Plot, Plot3D

from Psience.Wavefun import Wavefunction, Wavefunctions
from .BaseDVR import DVRResults

__all__ = ["DVRWavefunctions", "DVRWavefunction"]

class DVRWavefunction(Wavefunction):
    def plot(self, figure=None, grid=None, index=0, scaling=1, shift=0, **opts):
        """
        Plots a single wave function on the grid

        :param figure:
        :type figure:
        :param grid:
        :type grid:
        :param index:
        :type index:
        :param scaling:
        :type scaling:
        :param shift:
        :type shift:
        :param opts:
        :type opts:
        :return:
        :rtype:
        """

        if grid is None:
            grid = self.parent.results

        dim = len(grid.shape)
        if dim > 1 and grid.shape[-1] == dim-1: # check whether we have a mesh of points that we need to reshape
            unroll = np.roll(np.arange(len(grid.shape)), 1)
            grid = grid.transpose(unroll)

        if not isinstance(scaling, (int, float, np.integer, np.floating)):
            scaling = scaling[index]
        if not isinstance(shift, (int, float, np.integer, np.floating)):
            shift = shift[index]

        if dim == 1:
            return Plot(grid, self.data*scaling+shift, figure=figure, **opts)
        else:
            return Plot3D(*grid, self.data.reshape(grid[0].shape), figure=figure, **opts)

    def expectation(self, op, other):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """
        import numpy as np

        wf = op(self.data)
        if not isinstance(other, np.ndarray):
            other = other.data
        return np.dot(wf, other)

    def probability_density(self):
        """Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        import numpy as np

        return np.power(self.data, 2)

class DVRWavefunctions(Wavefunctions):
    # most evaluations are most efficient done in batch for DVR wavefunctions so we focus on the batch object
    def __init__(self, energies=None, wavefunctions=None, grid=None, wavefunction_class=DVRWavefunction, results:DVRResults=None, **opts):
        super().__init__(energies=energies, wavefunctions=wavefunctions, wavefunction_class=wavefunction_class, **opts)
        self.results = results
        self.grid = grid
    def __repr__(self):
        return "{}(num={}, DVR={})".format(
            type(self).__name__,
            len(self),
            self.results.parent
        )

    def __len__(self):
        return len(self.energies)
    def __iter__(self):
        for i in range(len(self)):
            yield self.__getitem__(i)
    def __getitem__(self, item):
        """
        Provides a single `DVRWavefunction` or slice of `DVRWavefunctions`
        :param item:
        :type item:
        :return:
        :rtype: DVRWavefunction | DVRWavefunctions
        """
        if not isinstance(item, (int, np.integer)):
            return type(self)(
                energies=self.energies[item],
                wavefunctions=self.wavefunctions[:, item].reshape((len(self.wavefunctions), -1)),
                wavefunction_class=self.wavefunction_class,
                results=self.results,
                grid=self.grid,
                **self.opts
            )
        else:
            return self.wavefunction_class(
                self.energies[item],
                self.wavefunctions[:, item],
                parent=self,
                **self.opts
            )

    def plot(self, figure=None, graphics_class=None, plot_style=None, scaling=1, shift=0, **opts):
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

        dim = len(grid.shape)
        if dim > 1 and grid.shape[-1] == dim-1: # check whether we have a mesh of points that we need to reshape
            grid = np.moveaxis(grid, grid.ndim, 0)

        return super().plot(
            figure=figure,
            graphics_class=graphics_class,
            plot_style=plot_style,
            grid=grid,
            scaling=scaling,
            shift=shift,
            **opts
        )

    def expectation(self, op, other=None):
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
        if isinstance(op, np.ndarray):
            wfs = self.wavefunctions
            for _ in range(op.ndim-1):
                wfs = np.expand_dims(wfs, -1)
            # print(np.expand_dims(op, 1).shape, wfs.shape)
            wfs = np.expand_dims(op, 1) * wfs
            # print(self.wavefunctions.shape, wfs.shape)
        else:
            wfs = op(self.wavefunctions)
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

    def probability_density(self):
        """Computes the probability density of the set of wavefunctions

        :return:
        :rtype:
        """
        return np.power(self.wavefunctions, 2)

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