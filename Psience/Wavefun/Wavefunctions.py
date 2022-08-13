"""
Provides very general support for an abstract wavefunction object
Allows different methods to provide their own concrete implementation details
"""
from abc import *
import numpy as np

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
    @abstractmethod
    def plot(self, figure = None, index = None, **opts):
        """Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)

        :param figure:
        :type figure: Graphics | Graphics3D
        :return:
        :rtype:
        """
        pass
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
        if isinstance(which, slice):
            return type(self)(
                energies=self.energies[which],
                wavefunctions=self.wavefunctions[which],
                wavefunction_class=self.wavefunction_class,
                indices=inds[which],
                **self.opts
            )
        else:
            return self.wavefunction_class(self.energies[which], self.wavefunctions[which], parent=self, index=inds[which], **self.opts)
    def __getitem__(self, item):
        """Returns a single Wavefunction object"""
        # iter comes for free with this
        return self.get_wavefunctions(item)
    def __iter__(self):
        eng_list = list(self.energies)
        wf_list = list(self.wavefunctions) # is this right?
        for eng, wfn in zip(eng_list, wf_list):
            yield self.wavefunction_class(eng, wfn, parent = self, **self.opts)

    def frequencies(self, start_at = 0):
        return self.energies[1+start_at:] - self.energies[start_at]

    def plot(self, figure=None, graphics_class=None, **opts):
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
