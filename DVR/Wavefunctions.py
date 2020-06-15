"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""

from Psience.Wavefun import Wavefunction, Wavefunctions

class DVRWavefunction(Wavefunction):
    def plot(self, figure = None, grid = None, **opts):
        import numpy as np

        if grid is None:
            grid = self.opts['grid']

        dim = len(grid.shape)
        if dim > 1 and grid.shape[-1] == dim-1: # check whether we have a mesh of points that we need to reshape
            unroll = np.roll(np.arange(len(grid.shape)), 1)
            grid = grid.transpose(unroll)

        if dim == 1:
            if figure is None:
                from McUtils.Plots import Plot
                return Plot(grid, self.data, **opts)
            else:
                return figure.plot(grid, self.data, **opts)
        else:
            if figure is None:
                from McUtils.Plots import Plot3D
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
    def plot(self, figure = None, graphics_class = None, plot_style = None, **opts):
        import numpy as np

        grid = self.opts['grid']

        dim = len(grid.shape)
        if dim > 1 and grid.shape[-1] == dim-1: # check whether we have a mesh of points that we need to reshape
            unroll = np.roll(np.arange(len(grid.shape)), 1)
            grid = grid.transpose(unroll)

        super().plot(
            figure=figure,
            graphics_class=graphics_class,
            plot_style=plot_style,
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
        import numpy as np

        wfs = op(self.wavefunctions)
        if not isinstance(other, np.ndarray):
            other = other.wavefunctions
        return np.dot(wfs, other)

    def probability_density(self):
        """Computes the probability density of the set of wavefunctions

        :return:
        :rtype:
        """
        import numpy as np

        return np.power(self.wavefunctions, 2)