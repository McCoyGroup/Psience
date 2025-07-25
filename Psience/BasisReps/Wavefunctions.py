"""
Provides Wavefunctions to use with the bases we define
"""

import numpy as np
from McUtils.Plots import *

from ..Wavefun import *

from .Bases import RepresentationBasis
from .Representations import Representation
from .Operators import Operator

__all__ = [
    "AnalyticWavefunctions",
    "AnalyticWavefunction",
    "WavefunctionBasis"
]

class WavefunctionBasis(RepresentationBasis):
    def __init__(self, wavefunctions:Wavefunctions):
        self.wfns = wavefunctions
        super().__init__(
            self.wfns.__getitem__,
            len(self.wfns)
        )

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.wfns == other.wfns
    def x(self, n):
        return self.wfns.coordinate()[:n, :n]
    def p(self, n):
        return self.wfns.momentum()[:n, :n]
    def p2(self, n):
        return self.wfns.laplacian()[:n, :n]
    def kinetic_energy(self, n):
        return self.wfns.kinetic_energy()[:n, :n]

class AnalyticWavefunction(Wavefunction):
    """
    Little extension to RepresentationBasis so that we can use p and x and stuff
    to evaluate out matrix elements and stuff.
    This will be more progressively more tightly meshed with RepresentationBasis in the future,
    but for now I just want to provide the scaffolding so that I can build off of it.
    """
    def __init__(self, energy, data, **opts):
        super().__init__(energy, data, **opts)

    def evaluate(self, *args, **kwargs):
        return self.data(args, **kwargs)

    def plot(self, figure=None, plot_class = None, domain = (-5, 5), **opts):
        """Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)

        :param figure:
        :type figure: Graphics | Graphics3D
        :return:
        :rtype:
        """

        # I really need to add support for Plot(f, domain)
        # and ContourPlot(f, domain) where domain can have fixed values for some of the
        # elements...
        discrete = np.linspace(*domain, 100)
        data = self.evaluate(discrete, **opts)
        if plot_class is None:
            plot_class = Plot
            # this need to be more general, but that should really be done on the Plots side
        return plot_class(discrete, data, figure=figure, **opts)

    def expect(self, operator):
        """
        Provides expectation values of operators, but the operators have to be Operator objects...
          basically all the logic is inside the operator, but this is worth it for use in ExpansionWavefunction
        We can also potentially add support for ExpansionOperators or SymbolicOperators in the future that are
          able to very cleanly reuse stuff like the `p` matrix that a RepresentationBasis defines

        :param operator: the operator to take the expectation of
        :type operator: Operator
        """
        return operator[self.index, self.index]

    def expectation(self, op, other):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other: the other wavefunction
        :type other: AnalyticWavefunction
        :param op: the operator to take the matrix element of
        :type op: Operator
        :return:
        :rtype:
        """
        raise NotImplementedError("...woof")
        o = Representation(op, ...)
        # What kind of quanta info do I have?
        # These things clearly need to know about multi-dimensional systems
        return o[self.index, other.index]

    def probability_density(self):
        """Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        return self.data

    def marginalize_out(self, dofs):
        """
        Computes the projection of the current wavefunction onto a set of degrees
        of freedom

        :return:
        :rtype:
        """
        raise NotImplementedError("analytic wave function projections not yet implemented")

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
        return self.wavefunction_class(energy, wfn, index=which)