"""
Provides classes to support wave functions coming from VPT calculations
"""

import numpy as np, itertools as ip

from ..BasisReps import SimpleProductBasis, ExpansionWavefunctions, ExpansionWavefunction

from .Terms import DipoleTerms
from .Hamiltonian import PerturbationTheoryCorrections

__all__ = [
    'PerturbationTheoryWavefunction',
    'PerturbationTheoryWavefunctions',
]

class PerturbationTheoryWavefunction(ExpansionWavefunction):
    """
    These things are fed the first and second order corrections
    """
    def __init__(self, mol, basis, corrections):
        """
        :param mol: the molecule the wavefunction is for
        :type mol: Molecule
        :param basis: the basis the expansion is being done in
        :type basis: RepresentationBasis
        :param corrections: the corrections to the terms
        :type corrections: PerturbationTheoryCorrections
        """
        self.mol = mol
        self.corrs = corrections
        self.rep_basis = basis
        super().__init__(self.corrs.energies, self.corrs.wavefunctions, None)
    @property
    def order(self):
        return self.corrs.order
    def expectation(self, operator, other):
        return NotImplemented
    @property
    def zero_order_energy(self):
        return self.corrs.energies[0]

class PerturbationTheoryWavefunctions(ExpansionWavefunctions):
    """
    These things are fed the first and second order corrections
    """
    def __init__(self, mol, basis, corrections):
        """
        :param mol: the molecule the wavefunction is for
        :type mol: Molecule
        :param basis: the basis the expansion is being done in
        :type basis: SimpleProductBasis
        :param corrections: the corrections to the terms
        :type corrections: PerturbationTheoryCorrections
        """
        self.mol = mol
        self.corrs = corrections
        self.rep_basis = basis # temporary hack until I decided how to merge the idea of
                               # AnalyticWavefunctions with a RepresentationBasis
        self._tm_dat = None
        super().__init__(
            self.corrs.energies,
            self.corrs.wfn_corrections,
            None
        )

    @property
    def order(self):
        return self.corrs.order
    def expectation(self, operator, other):
        return NotImplemented
    @property
    def zero_order_energies(self):
        return self.corrs.energy_corrs[:, 0]
    def _transition_moments(self, mu_x, mu_y, mu_z):
        """
        Calculates the x, y, and z components of the
        transition moment between the wavefunctions stored

        :param mu_x: dipole x components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_x: Iterable[np.ndarray]
        :param mu_y: dipole y components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_y: Iterable[np.ndarray]
        :param mu_z: dipole z components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_z: Iterable[np.ndarray]
        :return:
        :rtype:
        """

        M = self.corrs.states['coupled_states'] # the coupled subspace space we're working in
        corr_vecs = self.corrs.corrections['wavefunctions'][..., M]
        transition_moment_components = np.zeros((3, 3)).tolist() # x, y, and z components of the 0th, 1st, and 2nd order stuff
        # raise Exception([mu_1[0].shape)
        mu = [mu_x, mu_y, mu_z]
        for a in range(3): #x, y, and z
            # ...I might need to add the constant term into this?
            mu_1, mu_2, mu_3 = mu[a]
            m1 = self.rep_basis.representation("x", coeffs=mu_1)
            m1 = m1[np.ix_(M, M)]
            m2 = 1/2*self.rep_basis.representation("x", "x", coeffs=mu_2)
            m2 = m2[np.ix_(M, M)]
            m3 = 1/6*self.rep_basis.representation("x", "x", "x", coeffs=mu_3)
            m3 = m3[np.ix_(M, M)]

            mu_terms = [m1, m2, m3]
            corr_terms = [corr_vecs[:, 0], corr_vecs[:, 1], corr_vecs[:, 2]]

            for q in range(3): # total quanta
                transition_moment_components[q][a] = [
                    np.dot(np.dot(corr_terms[i], mu_terms[k]), corr_terms[j].T)
                    for i, j, k in ip.product(range(q+1), range(q+1), range(q+1)) if i+j+k==q
                ]

        # we calculate it explicitly like this up front in case we want to use it later since the shape
        # can be a bit confusing ([0-order, 1-ord, 2-ord], [x, y, z])
        tmom = [
            sum(
                sum(ip.chain(transition_moment_components[i][j]))
                for i in range(3) # correction order
            ) for j in range(3) # xyz
        ]
        # raise Exception(tmom)
        return tmom, transition_moment_components

    @property
    def transition_moments(self):
        """
        Computes the transition moments between wavefunctions stored in the object

        :return:
        :rtype:
        """

        if self._tm_dat is None:
           self._tm_dat = self._transition_moments(*DipoleTerms(self.mol).get_terms())
        return self._tm_dat[0]

    @property
    def transition_moment_corrections(self):
        """
        Computes the transition moment corrections between wavefunctions stored in the object

        :return:
        :rtype:
        """

        if self._tm_dat is None:
            self._tm_dat = self._transition_moments(*DipoleTerms(self.mol).get_terms())
        return self._tm_dat[1]

    @property
    def oscillator_strengths(self):
        """
        Computes the oscillator strengths for transitions from the ground state to the other states

        :return:
        :rtype:
        """

        tms = self.transition_moments
        gs_tms = np.array([tms[i][0] for i in range(3)]).T
        osc = np.linalg.norm(gs_tms, axis=1) ** 2
        return osc

    @property
    def intensities(self):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        eng = self.energies
        units = 3.554206329390961e6
        return units * (eng - eng[0]) * self.oscillator_strengths

    @property
    def zero_order_intensities(self):
        """
        Computes the harmonic intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        eng = self.zero_order_energies
        tm = np.array([self.transition_moment_corrections[0][x][0][0] for x in range(3)]).T
        osc = np.linalg.norm(tm, axis=1)**2
        units = 3.554206329390961e6
        return units*(eng - eng[0]) * osc


# class PerturbationTheoryCorrections:
#     def __init__(self, corrections, energies, hamiltonian):
#         ...
#
#     def _fmt_corr_matrix(self, states, coupled_states, H1_blocks, coeffs):
#         """A useful debug function"""
#         h1_corr = UnitsData.convert("Hartrees", "Wavenumbers") * H1_blocks * coeffs
#         return "\n".join(
#             [
#                 "Corr:",
#                 "    " + (" ".join(" {2}{1}{0}".format(*x) for x in self.get_state_quantum_numbers(coupled_states))),
#                 *("{2}{1}{0} ".format(*s) + (" ".join("{:<+4.0f}".format(x) for x in h)) for s, h in
#                   zip(self.get_state_quantum_numbers(states), h1_corr))
#                 ]
#         )
#
#     def _fmt_corr2_matrix(self, states, coupled_states, H1_blocks, coeffs):
#         """A useful debug function"""
#         h1_corr = UnitsData.convert("Hartrees", "Wavenumbers") * H1_blocks
#         return "\n".join(
#             [
#                 "Corr:",
#                 "    " + (" ".join(" {2}{1}{0}".format(*x) for x in self.get_state_quantum_numbers(coupled_states))),
#                 *("{2}{1}{0} ".format(*s) + (" ".join("{:<+4.0f}".format(x) for x in h)) for s, h in
#                   zip(self.get_state_quantum_numbers(states), h1_corr))
#                 ]
#         )

# class PerturbationTheoryWavefunctions(Wavefunctions):
#     """
#     Represents a set of wavefunctions coming out of a VPT2 calculation.
#     Mostly just a wrapper on a PerturbationTheoryHamiltonian
#     """
#     def __init__(self, states, hamiltonian):
#         """
#         :param energies:
#         :type energies:
#         :param coeffs:
#         :type coeffs:
#         :param hamiltonian:
#         :type hamiltonian:
#         """
#         self.hamiltonian = hamiltonian
#         self.states = states
#         self._energies = None
#         super().__init__(
#             energies=None,
#             wavefunctions=None,
#             wavefunction_class=PerturbationTheoryWavefunction
#         )
#     @property
#     def energies(self):
#         return self.get_energies()
#     def get_energies(self):
#         ...