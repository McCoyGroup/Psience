"""
Provides Wavefunctions to use with the bases we define
"""

import numpy as np
from McUtils.Plots import *

from ..Wavefun import *

from .Representations import Representation
from .Operators import Operator

__all__ = [
    "AnalyticWavefunctions",
    "AnalyticWavefunction",
    "ExpansionWavefunctions",
    "ExpansionWavefunction"
]

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
    Simple wave function that takes a set of expansion coefficients alongside its basis.
    Technically this should be called a _linear expansion wave function_, but
    that was too long for my taste.
    """
    def __init__(self, energy, coefficients, basis_wfns):
        """
        :param energy: energy of the wavefunction
        :type energy: float
        :param coefficients: expansion coefficients
        :type coefficients: Iterable[float]
        :param basis_wfns: basis functions for the expansion
        :type basis_wfns: Wavefunctions
        """
        super().__init__(
            energy,
            {
                'coeffs':coefficients,
                'basis':basis_wfns
            }
        )

    @property
    def coeffs(self):
        return self.data['coeffs']
    @property
    def basis(self):
        return self.data['basis']

    def evaluate(self, *args, **kwargs):
        """
        Evaluates the wavecfunction as any other linear expansion.

        :param args: coordinates + any other args the basis takes
        :type args:
        :param kwargs: any keyword arguments the basis takes
        :type kwargs:
        :return: values of the wavefunction
        :rtype:
        """
        return np.dot(self.data['coeffs'], np.array([f(args, **kwargs) for f in self.data['basis']]))

    def expect(self, operator):
        """
        Provides the expectation value of the operator `op`.
        Uses the basis to compute the reps and then expands with the expansion coeffs.

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

    def expectation(self, op, other):
        """
        Computes the expectation value of operator `op` over the wavefunction `other` and `self`.
        **Note**: _the basis of `other`, `self`, and `op` are assumed to be the same_.

        :param op: an operator represented in the basis of the expansion
        :type op: Operator
        :param other: the other wavefunction to expand over
        :type other: ExpansionWavefunction
        :return:
        :rtype:
        """
        op_matrix = op[
                tuple(x.index for x in self.data['basis']),
                tuple(o.index for o in other.basis)
            ]
        # See the note about needing to handle multidimensional cases better in `expect`
        return np.dot(self.data('coeffs'), np.dot(op_matrix, other.coeffs))

    def probability_density(self):
        """Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        raise NotImplementedError

class ExpansionWavefunctions(Wavefunctions):
    """
    Simple expansion wave function rep that takes multiple sets of coefficients.
    As with all objects deriving from `Wavefunctions`, can be iterated through to
    provide a manifold of standalone `ExpansionWavefunction` objects.
    Currently there are major conceptual issues, as I need this to _both_ support `AnalyticWavefunctions`
    and `RepresentationBasis` as the basis...
    which means `AnalyticWavefunctions` needs to track a basis...
    but `AnalyticWavefunctions` wasn't really designed for that, so I need to go back and figure out how
    that binding is going to be managed.
    """
    def __init__(self, energies, coefficients, basis_wfns, **ops):
        """
        :param energies: energies for the stored wavefunctions
        :type energies: Iterable[float]
        :param coefficients: expansion coefficients
        :type coefficients: Iterable[Iterable[float]]
        :param basis_wfns: wavefunctions to use as the basis for the expansion
        :type basis_wfns: Wavefunctions
        :param ops: extra options for feeding through to `Wavefunctions`
        :type ops:
        """
        # self._coeffs = coefficients
        # self._energies = energies
        self._basis = basis_wfns
        if 'wavefunction_class' not in ops:
            ops['wavefunction_class'] = ExpansionWavefunction
        super().__init__(energies, coefficients, **ops)

    # def expect(self, op, other):
    #     """
    #     Provides expectation values of the wavefunctions o
    #     :param op:
    #     :type op:
    #     :param other:
    #     :type other:
    #     :return:
    #     :rtype:
    #     """
    #     return NotImplemented


    # def get_wavefunctions(self, which):
    #     energy = self.energies[which]
    #     wfn = self.wavefunctions[which]
    #     return self.wavefunction_class(energy, wfn, self._basis)