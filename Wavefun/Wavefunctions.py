from abc import *

__all__ = [
    "Wavefunction",
    "Wavefunctions",
    "WavefunctionException"
]

class WavefunctionException(Exception):
    pass

class Wavefunction:
    """Represents a single wavefunction object"""
    def __init__(self, energy, data, parent = None, **opts):
        self.energy = energy
        self.data   = data
        self.parent = parent
        self.opts   = opts
    @abstractmethod
    def plot(self, figure = None, **opts):
        """Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)

        :param figure:
        :type figure: Graphics | Graphics3D
        :return:
        :rtype:
        """
        pass
    @abstractmethod
    def expectation(self, op, other):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction
        :param op:
        :type op:
        :return:
        :rtype:
        """
        pass
    @abstractmethod
    def probability_density(self):
        """Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        return self.expectation(lambda a:a, self)

class Wavefunctions:
    """An object representing a set of wavefunctions.
    Provides concrete, but potentially inefficient methods for doing all the wavefunction ops.

    """

    def __init__(self, energies = None, wavefunctions = None, wavefunction_class = None, **opts):
        self.wavefunctions = wavefunctions
        self.energies = energies
        self.wavefunction_class = wavefunction_class
        self.opts = opts

    def get_wavefunctions(self, which):
        if isinstance(which, slice):
            return type(self)(
                energies = self.energies[which],
                wavefunctions = self.wavefunctions[which],
                wavefunction_class = self.wavefunction_class,
                **self.opts
            )
        else:
            return self.wavefunction_class(self.energies[which], self.wavefunctions[which], parent = self, **self.opts)
    def __getitem__(self, item):
        """Returns a single Wavefunction object"""
        # iter comes for free with this
        return self.get_wavefunctions(item)
    def __iter__(self):
        for eng,wfn in zip(self.energies, self.wavefunctions):
            yield self.wavefunction_class(eng, wfn, parent = self, **self.opts)

    def frequencies(self, start_at = 0):
        return self.energies[1+start_at:] - self.energies[start_at]

    def plot(self, figure = None, graphics_class = None, plot_style = None, **opts):
        """Plots all of the wavefunctions on one set of axes

        :param opts:
        :type opts:
        :return:
        :rtype:
        """
        from McUtils.Plots import Graphics, Graphics3D

        k = "plot_defaults"
        opts = dict(self.opts[k] if k in self.opts else (), **opts)

        if figure == None:
            dim = self.opts['dimension'] if 'dimension' in self.opts else 1
            if graphics_class is None:
                if dim ==1:
                    graphics_class = Graphics
                elif dim == 2:
                    graphics_class = Graphics3D
                else:
                    raise WavefunctionException(
                        "{}.{}: don't know how to plot wavefunctions of dimension {}".format(
                            type(self).__name__, 'plot', dim
                        )
                    )
            figure = graphics_class(**opts)

        if plot_style is None:
            plot_style = {}
        for wfn in self:
            wfn.plot(figure, **plot_style)

        return figure
