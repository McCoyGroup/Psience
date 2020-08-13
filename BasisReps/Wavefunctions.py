"""
Provides Wavefunctions to use with the bases we define
"""

import numpy as np
from McUtils.Plots import *

from ..Wavefun import *

from .Terms import TermComputer
from .Operators import Operator

class AnalyticWavefunction(Wavefunction):
    """
    Little extension to RepresentationBasis so that we can use p and x and stuff
    to evaluate out matrix elements and stuff
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
        o = TermComputer(op, ...)
        # What kind of quanta info do I have?
        # These things clearly need to know about multi-dimensional systems
        return o[self.index, other.index]

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
        return self.wavefunction_class(energy, wfn, index=which)

class ExpansionWavefunction(AnalyticWavefunction):
    """
    Wrapper that takes an expansion alongside its basis
    """

    def __init__(self, energy, coefficients, analytic_wavefunctions):
        super().__init__(energy, {
            'coeffs':coefficients,
            'basis':analytic_wavefunctions
        })

    def evaluate(self, *args, **kwargs):
        return np.dot(self.data['coeffs'], np.array([f(args, **kwargs) for f in self.data['basis']]))

    def expect(self, operator):
        """

        :param operator:
        :type operator:
        :return:
        :rtype:
        """
        op_vector = operator[
                tuple(x.index for x in self.data['basis']),
                tuple(x.index for x in self.data['basis'])
            ]
        # For multidimensional bases this will probably not work without modification
        # In that case, the operator returns a slice of a high-dimensional object,
        #  and we'll need multidimensional index information for this to work right
        return np.dot(self.data['coeffs'], op_vector)

    def expectation(self, operator, other):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction
        :param op:
        :type op:
        :return:
        :rtype:
        """
        op_matrix = operator[
                tuple(x.index for x in self.data['basis']),
                tuple(other.index for x in self.data['basis'])
            ]
        if not isinstance(G, int):
            # takes an (e.g.) 5-dimensional SparseTensor and turns it into a contracted 2D one
            subKE = pp[inds]
            if isinstance(subKE, np.ndarray):
                ke = np.tensordot(subKE.squeeze(), G, axes=[[0, 1], [0, 1]])
            else:
                ke = subKE.tensordot(G, axes=[[0, 1], [0, 1]]).squeeze()
        else:
            ke = 0
        # See the note about needing to handle multidimensional cases better in `expect`
        return

    def probability_density(self):
        """Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        return self.data


class ExpansionWavefunctions(AnalyticWavefunctions):
    """
    Represents wavefunctions with analytic forms, most commonly harmonic oscillators
    """
    def __init__(self, coeffs, basis_wavefunctions, **ops):
        self._coeffs = coeffs
        self._basis = basis_wavefunctions
        self._energy_expr = self.get_energies()
        self._wfn_expr = wfn_expr
        if 'wavefunction_class' not in ops:
            ops['wavefunction_class'] = AnalyticWavefunction
        super().__init__(**ops)
    def get_wavefunctions(self, which):
        energy = self._energy_expr(which, **self.opts)
        wfn = self._wfn_expr(which, **self.opts)
        return self.wavefunction_class(energy, wfn)


