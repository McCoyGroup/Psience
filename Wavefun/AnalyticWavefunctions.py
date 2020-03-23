from .Wavefunctions import *
import numpy as np, scipy.special as orthogs
from McUtils.Plots import *


class AnalyticWavefunction(Wavefunction):
    """
    Mini-wrapper so that we can do stuff like evaluate our wavefunctions...
    """
    def evaluate(self, *args, **kwargs):
        return self.data(args, **kwargs)

    def plot(self, figure=None, plot_class = None, domain = (-5, 5), **opts):
        """Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)

        :param figure:
        :type figure: Graphics | Graphics3D
        :return:
        :rtype:
        """

        discrete = np.linspace(*domain, 100)
        data = self.evaluate(discrete, **opts)
        if plot_class is None:
            plot_class = Plot # this need to be more general
        return plot_class(discrete, data, figure=figure, **opts)

    def expect(self, operator, domain = (-5, 5), **opts):
        # We should add a Mesh.from_domain method and a Integrator class to Zachart to do all this
        # raise WavefunctionException("AnalyticWavefunctions don't have a general way to get expectation values..")
        discrete = np.linspace(domain, 100)
        data = self.evaluate(discrete)
        raise WavefunctionException("This will require update to McUtils.Zachary first")

    def expectation(self, op, other):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction
        :param op:
        :type op:
        :return:
        :rtype:
        """
        raise WavefunctionException("This will require update to McUtils.Zachary first")

    def probability_density(self):
        """Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        return self.data

class AnalyticWavefunctions(Wavefunctions):
    """
    Represents wavefunctions with analytic forms, most commonly harmonic oscillators
    """
    def __init__(self, energy_expr, wfn_expr, **ops):
        self._energy_expr = energy_expr
        self._wfn_expr = wfn_expr
        if 'wavefunction_class' not in ops:
            ops['wavefunction_class'] = AnalyticWavefunction
        super().__init__(**ops)
    def get_wavefunctions(self, which):
        energy = self._energy_expr(which, **self.opts)
        wfn = self._wfn_expr(which, **self.opts)
        return self.wavefunction_class(energy, wfn)

class HarmonicOscillatorWavefunctions(AnalyticWavefunctions):
    def __init__(self, frequency = 1, hb = 1):
        energy = lambda n: hb*frequency*(n+1/2)
        wfn = lambda n: lambda x, H=orthogs.hermite(n): H(x)*np.exp(-(x**2))
        super().__init__(energy, wfn)

# need a clean way to express multi-quanta H.O.s for the VPT stuff


