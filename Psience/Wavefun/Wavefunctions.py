"""
Provides very general support for an abstract wavefunction object
Allows different methods to provide their own concrete implementation details
"""
import abc
from abc import *
import numpy as np

from McUtils.Plots import Graphics, Plot, TriContourPlot
from ..Spectra import DiscreteSpectrum

__all__ = [
    "Wavefunction",
    "Wavefunctions",
    "WavefunctionException"
]

class WavefunctionException(Exception):
    pass

class Wavefunction:
    """Represents a single wavefunction object"""
    def __init__(self, energy, data, parent=None, index=None, **opts):
        self.energy = energy
        self.data   = data
        self.parent = parent
        self.index = index
        self.opts   = opts

    @abc.abstractmethod
    def get_dimension(self):
        raise NotImplementedError("abstract base method")
    @property
    def ndim(self):
        return self.get_dimension()

    def plot(self,
             figure=None, domain=None, grid=None, values=None, plot_points=100,
             index=0, scaling=1, shift=0, plotter=None, plot_density=False,
             zero_tol=1e-8, contour_levels=None,
             **opts
             ):
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

        if grid is None and domain is None:
            raise ValueError("can't plot a wave function without a specified domain")

        if grid is None:
            if isinstance(domain[0], (int, np.integer, float, np.floating)):
                domain = [domain]
            if isinstance(plot_points, (int, np.integer)):
                plot_points = [plot_points] * len(domain)

            grids = []
            for dom, pts in zip(domain, plot_points):
                grids.append(np.linspace(*dom, pts))
            grid = np.moveaxis(np.array(np.meshgrid(*grids)), 0, -1).reshape(-1, len(domain)) # vector of points

        grid = np.asanyarray(grid)
        if grid.ndim == 1:
            grid = grid[:, np.newaxis]
        elif grid.ndim > 2:
            grid = grid.reshape(-1, grid.shape[-1])
        dim = grid.shape[-1]

        if dim > 2 and plotter is None: # if people want to try, let 'em
            raise ValueError("can't plot data with dimension higher than 2, take a projection first")

        # allows us to scale wave functions independently
        if not isinstance(scaling, (int, float, np.integer, np.floating)):
            scaling = scaling[index]
        if not isinstance(shift, (int, float, np.integer, np.floating)):
            shift = shift[index]

        if values is None:
            if plot_density:
                values = self.probability_density(grid)
            else:
                values = self.evaluate(grid)
        values[np.abs(values) < zero_tol] = 0.

        values = values * scaling + shift

        if contour_levels is not None and 'levels' not in opts:
            max_val = np.max(np.abs(values))
            opts['levels'] = np.linspace(-max_val, max_val, contour_levels)

        if plotter is None:
            if dim == 1:
                plotter = Plot
            else:
                plotter = TriContourPlot

        return plotter(*grid.T, values, figure=figure, **opts)

    def projection_plot(self,
                        coords,
                        figure=None,
                        **plot_options
                        ):
        """
        A convenience function to plot multiple projections
        on the same set of axes

        :param coords:
        :type coords:
        :param figure:
        :type figure:
        :param plot_options:
        :type plot_options:
        :return:
        :rtype:
        """
        if isinstance(coords, (int, np.integer)):
            coords = [[coords]]
        elif isinstance(coords[0], (int, np.integer)):
            coords = [[c] for c in coords]

        for proj_inds in coords:
            fig = self.project(proj_inds).plot(
                figure=figure,
                **plot_options
            )
            if figure is None:
                figure = fig

        return figure

    @abstractmethod
    def expectation(self, op, other=None):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction
        :param op:
        :type op:
        :return:
        :rtype:
        """
        pass
    def overlap(self, other):
        return self.expectation(lambda w:w, other=other)
    @abstractmethod
    def evaluate(self, points):
        """
        Evaluates the current wavefunction

        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base method")
    @property
    def probability_density(self):
        """
        Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        return lambda pts:self.evaluate(pts)**2 # we're assuming real I guess...
    @abstractmethod
    def marginalize_out(self, dofs):
        """
        Integrates out the contributions from the degrees of freedom `dofs`

        :return:
        :rtype: Wavefunction
        """
        raise NotImplementedError("abstract base method")
    def project(self, dofs):
        """
        Computes the projection of the current wavefunction onto a set of degrees
        of freedom, returning a projected wave function object

        :return:
        :rtype: Wavefunction
        """

        dof_complement = np.setdiff1d(np.arange(self.ndim), dofs)
        return self.marginalize_out(dof_complement)

class Wavefunctions:
    """
    An object representing a set of wavefunctions.
    Provides concrete, but potentially inefficient methods for doing all the wavefunction ops.

    """
    wavefunction_class = Wavefunction
    def __init__(self,
                 energies=None, wavefunctions=None,
                 indices=None, wavefunction_class=None, **opts):
        self.wavefunctions = wavefunctions
        self.energies = energies
        self.wavefunction_class = self.wavefunction_class if wavefunction_class is None else wavefunction_class
        self.indices = indices
        self.opts = opts

    def get_wavefunctions(self, which):
        inds = self.indices
        if inds is None:
            inds = np.arange(len(self.wavefunctions))
        if not isinstance(which, (int, np.integer)):
            return type(self)(
                energies=self.energies[which],
                wavefunctions=self.wavefunctions[:, which],
                wavefunction_class=self.wavefunction_class,
                indices=inds[which],
                **self.opts
            )
        else:
            return self.wavefunction_class(
                self.energies[which],
                self.wavefunctions[:, which],
                parent=self,
                index=inds[which],
                **self.opts
            )
    def __getitem__(self, item):
        """Returns a single Wavefunction object"""
        # iter comes for free with this
        return self.get_wavefunctions(item)
    def __len__(self):
        return len(self.energies)
    def __iter__(self):
        for i in range(len(self)):
            yield self.__getitem__(i)

    def frequencies(self, start_at = 0):
        return np.concatenate([self.energies[:start_at], self.energies[1+start_at:]]) - self.energies[start_at]

    def get_spectrum(self,
                     dipole_function,
                     start_at=0,
                     **options
                     ):
        freqs = self.frequencies(start_at=start_at)
        transition_moments = self.expectation(dipole_function,
                                              **options,
                                              )[start_at]
        transition_moments = np.concatenate([
            transition_moments[:start_at],
            transition_moments[start_at+1:]
        ])
        return DiscreteSpectrum.from_transition_moments(
            freqs,
            transition_moments
        )

    def plot(self, figure=None, graphics_class=None, **opts):
        """Plots all of the wavefunctions on one set of axes

        :param opts:
        :type opts:
        :return:
        :rtype:
        """

        k = "plot_defaults"
        opts = dict(self.opts[k] if k in self.opts else (), **opts)

        if figure == None:
            dim = self.opts['dimension'] if 'dimension' in self.opts else 1
            if graphics_class is None:
                if dim ==1:
                    graphics_class = Graphics
                elif dim == 2:
                    graphics_class = Graphics#Graphics3D
                else:
                    raise WavefunctionException(
                        "{}.{}: don't know how to plot wavefunctions of dimension {}".format(
                            type(self).__name__, 'plot', dim
                        )
                    )
            figure = graphics_class(strict=False, **opts)

        for i, wfn in enumerate(self):
            ind = wfn.index
            if ind is None:
                ind = i
            wfn.plot(figure, index=ind, **opts)

        return figure

    def expectation(self, op, other=None):
        """
        Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunctions
        :param op:
        :type op:
        :return:
        :rtype:
        """
        if other is None:
            other = self

        res = []
        for wfn in self:
            subres = []
            for ofn in other:
                subres.append(wfn.expectation(op, ofn))
            res.append(subres)
        return np.array(res)
    def overlap(self, other):
        return self.expectation(lambda w:w, other=other)

    def coordinate(self):
        """
        Provides the coordinate operator in the wavefunction basis

        :return:
        :rtype:
        """
        raise NotImplementedError("no coordinate rep implemented for {}".format(self))
    def momentum(self):
        """
        Provides the real part of the representation of the momentum operator in the wavefunction basis

        :return:
        :rtype:
        """
        raise NotImplementedError("no momentum implemented for {}".format(self))
    def laplacian(self):
        """
        Provides the representation of the laplacian in the wavefunction basis

        :return:
        :rtype:
        """
        raise NotImplementedError("no momentum implemented for {}".format(self))
    def kinetic_energy(self):
        """
        Provides the representation of the KE in the wavefunction basis

        :return:
        :rtype:
        """
        raise NotImplementedError("no KE implemented for {}".format(self))