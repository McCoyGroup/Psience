from abc import *
from McUtils.Plots import Graphics

class Wavefunction:
    """Represents a single wavefunction object"""
    def __init__(self, energy, data, parent = None, **opts):
        self.energy = energy
        self.data   = data
        self.parent = parent
        self.opts   = opts
    @abstractmethod
    def plot(self, figure = None):
        """Uses matplotlib to plot the wavefunction on the passed figure (makes a new one if none)

        :param figure:
        :type figure: Graphics
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

    def __getitem__(self, item):
        """Returns a single Wavefunction object"""
        # iter comes for free with this
        if isinstance(item, slice):
            return type(self)(
                energies = self.energies[item],
                wavefunctions = self.wavefunctions[item],
                wavefunction_class = self.wavefunction_class,
                **self.opts
            )
        else:
            self.wavefunction_class(self.energies[item], self.wavefunctions[item], parent = self, **self.opts)
    def __iter__(self):
        for eng,wfn in zip(self.energies, self.wavefunctions):
            yield self.wavefunction_class(eng, wfn, parent = self, **self.opts)

    def plot(self, **opts):
        """Plots all of the wavefunctions on one set of axes

        :param opts:
        :type opts:
        :return:
        :rtype:
        """

        k = "plot_defaults"
        opts = dict(self.opts[k] if k in self.opts else (), **opts)

        p = Graphics(**opts)
        for wfn in self:
            wfn.plot(p)

        return p