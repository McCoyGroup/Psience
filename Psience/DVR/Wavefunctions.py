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
            if figure is None:
                return Plot(grid, self.data*scaling+shift, **opts)
            else:
                return figure.plot(grid, self.data*scaling+shift, **opts)
        else:
            if figure is None:
                return Plot3D(*grid, self.data.reshape(grid[0].shape), **opts)
            else:
                return figure.plot(*grid, self.data.reshape(grid[0].shape), **opts)

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
    def __init__(self, energies=None, wavefunctions=None, wavefunction_class=DVRWavefunction, results:DVRResults=None, **opts):
        super().__init__(energies=energies, wavefunctions=wavefunctions, wavefunction_class=wavefunction_class, **opts)
        self.results = results
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

        grid = self.results.grid

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

    def expectation(self, op, other):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: DVRWavefunctions | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """

        wfs = op(self.wavefunctions)
        if not isinstance(other, np.ndarray):
            other = other.wavefunctions
        return np.dot(wfs, other)

    def probability_density(self):
        """Computes the probability density of the set of wavefunctions

        :return:
        :rtype:
        """

        return np.power(self.wavefunctions, 2)