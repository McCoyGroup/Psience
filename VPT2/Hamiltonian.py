"""
Provides support for build perturbation theory Hamiltonians
"""

import numpy as np, itertools, enum

from McUtils.Numputils import SparseArray, vec_outer
from McUtils.Misc import Logger

from ..Molecools import Molecule
from ..BasisReps import HarmonicOscillatorBasis, SimpleProductBasis

from .Common import PerturbationTheoryException
from .Terms import PotentialTerms, KineticTerms, CoriolisTerm, PotentialLikeTerm

__all__ = [
    'PerturbationTheoryHamiltonian',
    'PerturbationTheoryCorrections'
]


class PerturbationTheoryCorrections:
    """
    Represents a set of corrections from perturbation theory.
    Can be used to correct other operators in the basis of the original calculation.

    """
    def __init__(self,
                 states,
                 corrections,
                 hamiltonians
                 ):
        """
        :param states: a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
        :type states: dict
        :param corrections: the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT
        :type corrections: dict
        :param hamiltonians: the set of Hamiltonian matrices used as an expansion
        :type hamiltonians: Iterable[np.ndarray]
        """
        self.states = states['states']
        self.coupled_states = states['coupled_states']
        self.total_basis = states['total_states']
        self.energy_corrs = corrections['energies']
        self.wfn_corrections = corrections['wavefunctions']
        if 'degenerate_states' in states:
            self.degenerate_states = states['degenerate_states']
        else:
            self.degenerate_states = None

        if 'degenerate_transformation' in corrections:
            self.degenerate_transf = corrections['degenerate_transformation']
        else:
            self.degenerate_transf = None

        if 'degenerate_energies' in corrections:
            self.degenerate_energies = corrections['degenerate_energies']
        else:
            self.degenerate_energies = None
        self.hams = hamiltonians

    @property
    def degenerate(self):
        return self.degenerate_transf is not None

    @property
    def energies(self):
        if self.degenerate:
            return self.degenerate_energies
        else:
            return np.sum(self.energy_corrs, axis=1)

    @property
    def order(self):
        return len(self.energy_corrs[0])

    def operator_representation(self, operator_expansion, order=None, subspace=None):
        """
        Generates the representation of the operator in the basis of stored states

        :param operator_expansion: the expansion of the operator
        :type operator_expansion: Iterable[float] | Iterable[np.ndarray]
        :param order: the order of correction to go up to
        :type order: Iterable[float] | Iterable[np.ndarray]
        :param subspace: the subspace of terms in which the operator expansion is defined
        :type subspace: Iterable[int]
        :return: the set of representation matrices for this operator
        :rtype: Iterable[np.ndarray]
        """

        mordor = self.order - 1
        if order is None:
            order = mordor
        if order > mordor:
            raise PerturbationTheoryException("{}: can't correct up to order {} when zero-order states were only corrected up to order {}".format(
                type(self).__name__,
                order,
                mordor
            ))
        order = order + 1 # so that we actually do get up to the request order after accounting for the zeros...
        if len(operator_expansion) < order:
            operator_expansion = list(operator_expansion) + [0]*(order - len(operator_expansion))

        wfn_corrs = [
            self.wfn_corrections[:, k, subspace] if subspace is not None else self.wfn_corrections[:, k, subspace]
            for k in range(order)
        ]

        reps = [np.zeros(1)] * order
        for k in range(order):
            op = None
            for a in range(k+1):
                for b in range(k-a+1):
                    c = k - (a + b)
                    rop = operator_expansion[c]
                    if isinstance(rop, (int, float, np.integer, np.floating)): # constant reps...
                        if rop != 0: # cheap easy check
                            subrep = rop * np.dot(wfn_corrs[a], wfn_corrs[b].T)
                            if op is None:
                                op = subrep
                            else:
                                op += subrep
                    else:
                        subrep = np.dot(np.dot(wfn_corrs[a], rop), wfn_corrs[b].T)
                        if op is None:
                            op = subrep
                        else:
                            op += subrep
            reps[k] = op

        return reps

    def savez(self, file):
        keys = dict(
            states=self.states,
            coupled_states=self.coupled_states,
            total_states=self.total_basis,
            energies=self.energy_corrs,
            wavefunctions=self.wfn_corrections,
            hamiltonians=self.hams
        )
        if self.degenerate_states is not None:
            keys['degenerate_states'] = self.degenerate_states
        if self.degenerate_transf is not None:
            keys['degenerate_transformation'] = self.degenerate_transf
        if self.degenerate_energies is not None:
            keys['degenerate_energies'] = self.degenerate_energies
        np.savez(file, **keys)

    @classmethod
    def loadz(cls, file):
        keys = np.load(file)
        return cls(
            {
                "states":keys['states'],
                 "coupled_states":keys['coupled_states'],
                 "total_states":keys['total_states'],
                 "degenerate_states":keys['degenerate_states'] if 'degenerate_states' in keys else None
             },
            {
                "energies":keys['energies'],
                "wavefunctions":keys['wavefunctions'],
                "degenerate_transformation": keys['degenerate_transformation'] if 'degenerate_transformation' in keys else None,
                "degenerate_energies": keys['degenerate_energies'] if 'degenerate_energies' in keys else None
            },
            keys['hamiltonians']
        )


# PerturbationTheoryCorrections = collections.namedtuple(
#     "PerturbationTheoryCorrections",
#     [
#         "energies",
#         "wavefunctions",
#         "states",
#         "corrections",
#         "hamiltonians"
#     ]
# )

class PerturbationTheoryHamiltonian:
    """
    Represents the main Hamiltonian used in the perturbation theory calculation.
    Uses a harmonic oscillator basis for representing H0, H1, and H2 (and only goes up to H2 for now)
    """

    def __init__(self,
                 molecule=None,
                 n_quanta=None,
                 modes=None,
                 mode_selection=None,
                 coriolis_coupling = True,
                 log = None
                 ):
        """
        :param molecule: the molecule on which we're doing perturbation theory
        :type molecule:  Molecule
        :param n_quanta: the numbers of quanta to use when representing the entire state space
        :type n_quanta: int | None
        :param modes: the set of modes to use as the basis
        :type modes: None | MolecularNormalModes
        :param mode_selection: the subset of modes to use when doing expansions
        :type mode_selection: None | Iterable[int]
        :param coriolis_coupling: whether to add coriolis coupling if not in internals
        :type coriolis_coupling: bool
        """

        if molecule is None:
            raise PerturbationTheoryException("{} requires a Molecule to do its dirty-work")
        # molecule = molecule.get_embedded_molecule()
        self.molecule = molecule
        if modes is None:
            modes = molecule.normal_modes
        mode_n = modes.basis.matrix.shape[1] if mode_selection is None else len(mode_selection)
        self.mode_n = mode_n
        if n_quanta is None:
            n_quanta = 10 # dunno yet how I want to handle this since it should really be defined by the order of state requested...
        self.n_quanta = np.full((mode_n,), n_quanta) if isinstance(n_quanta, (int, np.int)) else tuple(n_quanta)
        self.modes = modes
        self.V_terms = PotentialTerms(self.molecule, modes=modes, mode_selection=mode_selection)
        self.G_terms = KineticTerms(self.molecule, modes=modes, mode_selection=mode_selection)
        if coriolis_coupling and (self.molecule.internal_coordinates is None):
            self.coriolis_terms = CoriolisTerm(self.molecule, modes=modes, mode_selection=mode_selection)
        else:
            self.coriolis_terms = None
        self.watson_term = PotentialLikeTerm(self.molecule, modes=modes, mode_selection=mode_selection)

        self._h0 = self._h1 = self._h2 = None

        self.basis = SimpleProductBasis(HarmonicOscillatorBasis, self.n_quanta)

        if log is None or isinstance(log, Logger):
            self.logger = log
        elif log is True:
            self.logger = Logger()
        elif log is False:
            self.logger = None
        else:
            self.logger = Logger(log)

    @classmethod
    def from_fchk(cls, file,
                  internals=None,
                  mode_selection=None,
                  **kw
                  ):
        """
        :param file: fchk file to load from
        :type file: str
        :param internals: internal coordinate specification as a Z-matrix ordering
        :type internals: Iterable[Iterable[int]]
        :param n_quanta:
        :type n_quanta: int | Iterable[int]
        :return:
        :rtype:
        """

        molecule = Molecule.from_file(file, zmatrix=internals, mode='fchk')
        return cls(molecule=molecule, mode_selection=mode_selection, **kw)

    @property
    def H0(self):
        """
        Provides the representation for H0 in this basis
        """
        if self._h0 is None:
            self._h0 = (
                -1/2 * self.basis.representation('p', 'p', coeffs=self.G_terms[0])
                 + 1/2 * self.basis.representation('x', 'x', coeffs=self.V_terms[0])
            )

        return self._h0

    @property
    def H1(self):
        """
        Provides the representation for H1 in this basis
        """
        if self._h1 is None:
            self._h1 = (
                -1/2 * self.basis.representation('p', 'x', 'p', coeffs=self.G_terms[1], axes=[[0, 1, 2], [1, 0, 2]])
                + 1/6 * self.basis.representation('x', 'x', 'x', coeffs=self.V_terms[1])
            )
        return self._h1

    @property
    def H2(self):
        """
        Provides the representation for H2 in this basis
        """
        if self._h2 is None:
            self._h2 = (
                -1/4 * self.basis.representation('p', 'x', 'x', 'p',
                                                coeffs=self.G_terms[2], axes=[[0, 1, 2, 3], [2, 0, 1, 3]])
                + 1/24 * self.basis.representation('x', 'x', 'x', 'x', coeffs=self.V_terms[2])
            )
            if self.coriolis_terms is not None:
                total_cor = self.coriolis_terms[0]+self.coriolis_terms[1]+self.coriolis_terms[2]
                # import McUtils.Plots as plt
                # plt.ArrayPlot(
                #     total_cor.reshape(total_cor.shape[0] ** 2, total_cor.shape[0] ** 2)
                # ).show()
                self._h2 += -1 * self.basis.representation('x', 'p', 'x', 'p', coeffs=total_cor)
            else:
                self._h2 += 0 * self.basis.representation(coeffs=0)

            self._h2 += 1/8 * self.basis.representation(coeffs=self.watson_term[0])

        return self._h2

    @property
    def perturbations(self):
        return (self.H0, self.H1, self.H2)

    @staticmethod
    def _Nielsen_xss(s, w, v3, v4, zeta, Be, ndim):
        # actually pulled from the Stanton VPT4 paper since they had
        # the same units as I do...
        # we split this up into 3rd derivative, 4th derivative, and coriolis terms

        xss_4 = 1 / 16 * v4[s, s, s, s]
        xss_3 = -(
                5/48 * (v3[s, s, s] ** 2 / w[s])
                + 1/16 * sum((
                              (v3[s, s, t] ** 2) / w[t]
                              * (8 * (w[s] ** 2) - 3 * (w[t] ** 2))
                              / (4 * (w[s] ** 2) - (w[t] ** 2))
                      ) for t in range(ndim) if t != s
                      )
        )
        xss_cor = 0.
        return [xss_3, xss_4, xss_cor]

    @staticmethod
    def _Nielsen_xst(s, t, w, v3, v4, zeta, Be, ndim):
        # actually pulled from Stanton VPT4 paper
        xst_4 = 1 / 4 * v4[s, s, t, t]
        xst_3 = - 1 / 2 * (
                v3[s, s, t] ** 2 * w[s] / (4 * w[s] ** 2 - w[t] ** 2)
                + v3[s, t, t] ** 2 * w[t] / (4 * w[t] ** 2 - w[s] ** 2)
                + v3[s, s, s] * v3[s, s, t] / (2 * w[s])
                + v3[t, t, t] * v3[t, t, s] / (2 * w[t])
                - sum((
                              (
                                      (v3[s, t, r] ** 2) * w[r] * (w[s] ** 2 + w[t] ** 2 - w[r] ** 2)
                                      / (
                                          # This fucking delta_ijk term I don't know what it should be
                                          # because no one has it in my units
                                          # and none of the force-field definitions are consistent
                                              w[s] ** 4 + w[t] ** 4 + w[r] ** 4
                                              - 2 * ((w[s] * w[t]) ** 2 + (w[s] * w[r]) ** 2 + (w[t] * w[r]) ** 2)
                                      )
                              )
                              - v3[s, s, r] * v3[t, t, r] / (2 * w[r])
                      ) for r in range(ndim) if r != s and r != t
                      )
        )
        xst_cor = sum((
                                  Be[a] * (zeta[a, s, t] ** 2) * (w[t] / w[s] + w[t] / w[s])
                          ) for a in range(3))

        return [xst_3, xst_4, xst_cor]

    @staticmethod
    def _Barone_xss(s, w, v3, v4, zeta, Be, ndim):
        # unclear what the units should be from Barone's paper...
        l = w ** 2

        x_ss_coeff = 1 / (16 * l[s])
        x_ss_4 = x_ss_coeff * v4[s, s, s, s]
        x_ss_3 = -x_ss_coeff * (
                (5 / 3) * (v3[s, s, s] ** 2) / l[s]
                + sum((
                              (v3[s, s, t] ** 2) / l[t] * (8 * l[s] - 3 * l[t]) / (4 * l[s] - l[t])
                      ) for t in range(ndim) if t != s
                      )
        )
        x_ss_cor = 0.
        return  [x_ss_3, x_ss_4, x_ss_cor]

    @staticmethod
    def _Barone_xst(s, t, w, v3, v4, zeta, Be, ndim):
        from McUtils.Data import UnitsData
        # w2h = UnitsData.convert("Wavenumbers", "Hartrees")
        l = w**2
        x_coeff = 1 / np.sqrt(16 * w[s] * w[t])
        x4 = x_coeff * v4[s, s, t, t]
        x3 = - x_coeff * (
                + 2 * (v3[s, s, t] ** 2) / (4 * l[s] - l[t])
                + 2 * (v3[s, t, t] ** 2) / (4 * l[t] - l[s])
                + v3[s, s, s] * v3[s, t, t] / l[s]
                + v3[t, t, t] * v3[s, s, t] / l[t]
                - sum((
                              2 * (l[s] + l[t] - l[u]) * v3[s, t, u]
                              / (8 * w[s] * w[s] * w[u]) * (
                                      1 / (w[s] + w[t] + w[u])
                                      - 1 / (w[s] + w[t] - w[u])
                                      - 1 / (w[s] + w[u] - w[t])
                                      - 1 / (w[t] + w[u] - w[s])
                              )
                              - v3[s, s, u] * v3[t, t, u] / l[u]
                      ) for u in range(ndim) if u != t and u != s)
        )
        xcor = x_coeff * (
                4 * (l[s] + l[t]) * sum(( Be[a] * (zeta[a, s, t] ** 2) ) for a in range(3))
        )

        return  [x3, x4, xcor]

    @classmethod
    def _get_Nielsen_energies(cls, states, freqs, v3, v4, zeta, Be):
        """
        Returns energies using Harald Nielsen's formulae up to second order. Assumes no degeneracies.
        If implemented smarter, would be much, much faster than doing full-out perturbation theory, but less flexible.
        Good for validation, too.


        :param states: states to get energies for as lists of quanta in degrees of freedom
        :type states: Iterable[Iterable[int]]
        :param freqs: Harmonic frequencies
        :type freqs: np.ndarray
        :param v3: Cubic force constants
        :type v3: np.ndarray
        :param v4: Quartic force constants
        :type v4: np.ndarray
        :param zeta: Coriolis couplings
        :type zeta: np.ndarray
        :param Be: Moments of inertia
        :type Be: np.ndarray
        :return:
        :rtype:
        """

        ndim = len(freqs)
        x_mat_linear = np.array([
                cls._Nielsen_xss(s, freqs, v3, v4, zeta, Be, ndim) if s==t else
                cls._Nielsen_xst(s, t, freqs, v3, v4, zeta, Be, ndim)
                for s in range(ndim) for t in range(s, ndim)
          ]).T
        x_mat = np.zeros((3, ndim, ndim))
        ri, ci = np.triu_indices(ndim)
        x_mat[:, ri, ci] = x_mat_linear
        x_mat[:, ci, ri] = x_mat_linear

        from McUtils.Data import UnitsData
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        raise Exception(x_mat * h2w)

        states = np.array(states) + 1/2 # n+1/2 for harmonic vibrations

        x_mat = np.sum(x_mat, axis=0)
        e_harm = np.tensordot(freqs, states, axes=[0, 1])
        outer_states = vec_outer(states, states)
        # raise Exception(states, outer_states)
        e_anharm = np.tensordot(x_mat, outer_states, axes=[[0, 1], [1, 2]])

        return e_harm, e_anharm

    def get_Nielsen_energies(self, states):
        """

        :param states:
        :type states:
        :return:
        :rtype:
        """


        from McUtils.Data import UnitsData
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # TODO: figure out WTF the units on this have to be...

        freqs = self.modes.freqs
        v3 = self.V_terms[1]
        v4 = self.V_terms[2]

        # raise Exception(np.round( 6 * v3 * h2w))

        zeta, Be = self.coriolis_terms.get_zetas_and_momi()

        harm, anharm = self._get_Nielsen_energies(states, freqs, v3, v4, zeta, Be)

        # harm = harm / h2w
        anharm = anharm

        return harm, anharm

    @staticmethod
    def get_coupled_space(states, order, freqs=None, freq_threshold=None):
        """
        Returns the set of states that couple the given states up to the given order at each level of perturbation (beyond zero order)

        :param state: the states of interest
        :type state: Iterable[int]
        :param order: the order of perturbation theory we're doing
        :type order: int
        :param freqs: the zero-order frequencies in each vibrational mode being coupled
        :type freqs: Iterable[float]
        :param freq_threshold: the threshold for the maximum frequency difference between states to be considered
        :type freq_threshold: None | float
        :return: the sets of coupled states
        :rtype: tuple[tuple[int]]
        """

        states = np.array(states, dtype=int)
        nmodes = states.shape[-1]
        possible_h1_states = np.array([]) # the states that can be coupled through H1
        # first we generate the possible transformations +-1, +-3, +1+2, +1-2 -1+2, -1-2
        transitions_h1 = [
            [-1], [1], [-3], [3],
            [1, 2], [-1, 2], [1, -2], [-1, -2]
        ]
        state_rules_h1 = [(list(x) + [0], abs(sum(x))) for x in transitions_h1]
        permutations = np.array(
            sum((
                [
                    p for p in itertools.product(*([s] * nmodes))
                    if abs(sum(p)) == k
                ] for s, k in state_rules_h1),
                []
            ))
        initial_states = states
        for i in range(1, order, 2): # H1 only applies extra corrections odd orders
            new_states = np.concatenate([s[np.newaxis, :] + permutations for s in initial_states], axis=0)
            well_behaved = np.where(np.min(new_states, axis=1) >= 0)
            new_states = np.unique(new_states[well_behaved], axis=0)
            if len(possible_h1_states) > 0:
                possible_h1_states = np.unique(np.concatenate([possible_h1_states, new_states], axis=0), axis=0)
            else:
                possible_h1_states = new_states
            initial_states = possible_h1_states

        # from second order corrections
        possible_h2_states = states # the states that can be coupled to state through H2
        # first we generate the possible transformations +-2, +-4, +2-2, +2+2 -2+2, -2-2, -1+3, -1-3, +1+3, +1-3
        transitions_h2 = [
            [], [-2], [2], [-4], [4],
            [-2, -2], [-2, 2], [2, 2],
            [-1, -3], [-1, 3], [1, -3], [-3, 1]
        ]
        state_rules_h2 = [(list(x) + [0], abs(sum(x))) for x in transitions_h2]
        permutations = np.array(
            sum((
                [
                    p for p in itertools.product(*([s] * nmodes))
                    if abs(sum(p)) == k
                ] for s, k in state_rules_h2),
                []
            ))
        initial_states = states
        for i in range(2, order+1, 2): # H2 only applies extra corrections at even orders
            new_states = np.concatenate([s[np.newaxis, :] + permutations for s in initial_states], axis=0)
            well_behaved = np.where(np.min(new_states, axis=1) >= 0)
            new_states = np.unique(new_states[well_behaved], axis=0)
            if len(possible_h2_states) > 0:
                possible_h2_states = np.unique(np.concatenate([possible_h2_states, new_states], axis=0), axis=0)
            else:
                possible_h2_states = new_states
            initial_states = possible_h2_states

        #now if we have a frequency comb, we apply it
        if freq_threshold is not None:
            if freqs is None:
                raise ValueError("to apply frequency difference threshold, need harmonic frequencies")
            state_freq = np.sum(freqs[np.newaxis, :]*(states + 1/2), axis=1)

            h1_freq = np.sum(freqs[np.newaxis, :] * (possible_h1_states + 1/2), axis=1)
            h1_freq_diffs = np.abs(np.subtract.outer(state_freq, h1_freq))
            h1_thresh = h1_freq_diffs[0] < freq_threshold
            for d in h1_freq_diffs[1:]:
                # if any of the coupled states is below the threshold for any of the target states we include it
                # in our basis
                h1_thresh = np.logical_or(h1_thresh, d < freq_threshold)
            h1_sel = np.where(h1_thresh)
            h1_states = possible_h1_states[h1_sel]

            h2_freq = np.sum(freqs[np.newaxis, :] * (possible_h2_states + 1/2), axis=1)
            h2_freq_diffs = np.abs(np.subtract.outer(state_freq, h2_freq))
            h2_thresh = h2_freq_diffs[0] < freq_threshold
            for d in h2_freq_diffs[1:]:
                # if any of the coupled states is below the threshold for any of the target states we include it
                # in our basis
                h2_thresh = np.logical_or(h2_thresh, d < freq_threshold)
            h2_sel = np.where(h2_thresh)
            # print(h2_sel)
            h2_states = possible_h2_states[h2_sel]

        else:
            h1_states = possible_h1_states
            h2_states = possible_h2_states

        return h1_states, h2_states

    @staticmethod
    def _get_VPT_representations(
            h_reps,
            states,
            coupled_states
    ):
        """
        Gets the sparse representations of h_reps inside the basis of coupled states

        :param h_reps:
        :type h_reps:
        :param states:
        :type states:
        :param coupled_states:
        :type coupled_states:
        :param total_states:
        :type total_states:
        :return:
        :rtype:
        """

        if len(coupled_states) != len(h_reps) - 1:
            raise ValueError("coupled states must be specified for all perturbations (got {}, expected {})".format(
                len(coupled_states),
                len(h_reps) - 1
            ))

        # determine the total coupled space
        coupled_spaces = []
        states = np.array(states).astype(int)
        coupled_spaces.append(states)
        for m in coupled_states:
            m = np.array(m).astype(int)  # for safety
            coupled_spaces.append(m)
        total_coupled_space = np.unique(np.concatenate(coupled_spaces))

        # determine indices of subspaces within this total space
        coupled_space_inds = []
        sorter = np.argsort(total_coupled_space)
        for m in coupled_spaces:
            coupled_space_inds.append(np.searchsorted(total_coupled_space, m, sorter=sorter))

        # get explicit matrix reps inside the separate coupled subspaces
        N = len(total_coupled_space)
        # I should try to walk away from using scipy.sparse here and instead
        # shift to SparseArray, since it'll support swapping out the back end better...
        H = [np.zeros(1)] * len(h_reps)
        diag = h_reps[0][total_coupled_space, total_coupled_space] # this is just diagonal
        H[0] = SparseArray.from_diag(diag)
        for i,h in enumerate(h_reps[1:]):
            # calculate matrix elements in the coupled subspace
            m = coupled_spaces[i+1]
            if len(m) > 0:
                m_pairs = np.array(list(itertools.combinations_with_replacement(m, 2))).T
                # print(m_pairs)
                sub = h[m_pairs[0], m_pairs[1]]
            else:
                sub = 0
            if isinstance(sub, (int, np.integer, np.floating, float)):
                if sub == 0:
                    sub = SparseArray.empty((N, N), dtype=float)
                else:
                    raise ValueError("Using a constant shift of {} will force Hamiltonians to be dense...".format(sub))
                    sub = np.full((N, N), sub)
            else:
                # figure out the appropriate inds for this data in a sparse representation
                inds = coupled_space_inds[i+1]
                # upper triangle of indices
                up_tri = np.array(tuple(itertools.combinations_with_replacement(inds, 2)))
                # lower triangle is made by transposition
                low_tri = np.array([up_tri[:, 1], up_tri[:, 0]]).T
                # but npw we need to remove the duplicates, because many sparse matrix implementations
                # will sum up any repeated elements
                full_inds = np.concatenate([up_tri, low_tri])
                full_dat = np.concatenate([sub, sub])
                full_inds, idx = np.unique(full_inds, axis=0, return_index=True)
                full_dat = full_dat[idx]
                sub = SparseArray((full_dat, full_inds.T), shape=(N, N))

            H[i+1] = sub #type: np.ndarray

        return H, coupled_space_inds, total_coupled_space

    @staticmethod
    def _martin_test(h_reps, states, threshold):
        """
        Applies the Martin Test to a set of states and perturbations to determine which resonances need to be
        treated variationally. Everything is done within the set of indices for the representations.

        :param h_reps: The representation matrices of the perturbations we're applying.
        :type h_reps: Iterable[np.ndarray | SparseArray]
        :param states: The indices of the states to which we're going apply to the Martin test.
        :type states: Iterable[int]
        :param threshold: The threshold for what should be treated variationally (in the same energy units as the Hamiltonians)
        :type threshold: float
        :return: Pairs of coupled states
        :rtype: tuple[Iterable[int], Iterable[int]]
        """

        H0 = h_reps[0]
        H1 = h_reps[1]
        energies = np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag
        state_energies = energies[states]
        diffs = np.abs(np.array([e - energies for e in state_energies]))
        # raise Exception(diffs)
        for n, s in enumerate(states):
            diffs[n, s] = 1

        H1_blocks = H1[states, :]
        if isinstance(H1_blocks, SparseArray):
            H1_blocks = H1_blocks.toarray()
        corr_mat = (H1_blocks ** 4) / (diffs ** 3)

        # raise Exception(H1_blocks.shape, diffs.shape)

        deg_states = []
        for s, block in zip(states, corr_mat):
            big = np.where(np.abs(block) > threshold)[0]
            if len(big) > 0:
                deg_states.extend((s, d) for d in big)

        if len(deg_states) == 0:
            return None
        else:
            return np.array(deg_states).T

    @classmethod
    def _get_VPT_corrections(
            cls,
            h_reps,
            states,
            coupled_states,
            total_states,
            order,
            degenerate_states=None,
            logger=None
    ):
        """
        Applies perturbation theory to the constructed representations of H0, H1, etc.

        :param h_reps: series of perturbations as indexable objects
        :type h_reps: Iterable[np.ndarray]
        :param states: index of the states to get corrections for
        :type states: Iterable[int]
        :param coupled_states: indices of states to couple for each level of perturbation
        :type coupled_states: Iterable[Iterable[int]]
        :param total_states: the full number of state indices
        :type total_states: int
        :param order: the order of perturbation theory to apply
        :type order: int
        :param degenerate_states: the pairs of degeneracies to account for
        :type degenerate_states: None | (Iterable[int], Iterable[int])
        :return: the vibrational perturbation theory corrections for a single target state
        :rtype: PerturbationTheoryCorrections
        """

        # We use the iterative equations
        #            En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
        #     <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
        #           |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=1...k-1) + <n^(0)|n^(k)> |n^(0)>
        #  where Pi_n is the perturbation operator [1/(E_m-E_n) for m!=n]

        if logger is not None:
            logger.log_print(
                "\n    ".join([
                    "Computing PT corrections:",
                    "perturbations: {pert_num}",
                    "order: {ord}",
                    "states: {state_num}",
                    "basis sizes {basis_size}"
                ]),
                pert_num=len(h_reps) - 1,
                ord=order,
                state_num=len(states),
                basis_size=[len(b) for b in coupled_states]
            )
        H, coupled_inds, total_coupled_space = cls._get_VPT_representations(h_reps, states, coupled_states)
        del coupled_states # for safety as I debug
        N = len(total_coupled_space)
        state_inds = coupled_inds[0]
        if degenerate_states is not None:
            # we check this twice because the Martin test can return None
            if isinstance(degenerate_states, (int, np.integer, np.floating, float)):
                thresh = degenerate_states
                if logger is not None:
                    logger.log_print(
                        "    applying Martin test with threshold {}",
                        thresh
                    )
                degenerate_states = cls._martin_test(
                    H,
                    state_inds, # state indices in the coupled_states
                    thresh
                )
            elif all(isinstance(x, (int, np.integer)) for x in degenerate_states):
                # means we just got the degenerate subspace to work with
                pairs = [
                    l for l in itertools.product(degenerate_states, degenerate_states) if l[0] < l[1]
                    ]
                degenerate_states = np.array(pairs).T
            elif degenerate_states is not None:
                try:
                    degenerate_states = degenerate_states(H, coupled_states)
                except (TypeError, ValueError):
                    pass

        if degenerate_states is None:
            deg_vals = None
        else:
            if logger is not None:
                logger.log_print(
                    "    got {} degenerate states",
                    len(degenerate_states[0])
                )
            deg_i, deg_j = degenerate_states
            deg_vals = np.array([h[deg_i, deg_j] for h in H])

            deg_i, deg_j = degenerate_states
            H_non_deg = [np.zeros(1)] * len(H)
            for i,h in enumerate(H):
                h = h.copy()
                h[deg_i, deg_j] = 0.
                h[deg_j, deg_i] = 0.
                H_non_deg[i] = h
            H = H_non_deg


        all_energies = np.zeros((len(states), order + 1))
        all_overlaps = np.zeros((len(states), order + 1))
        all_corrs = np.zeros((len(states), order + 1, N))
        all_wfns = np.zeros((len(states), order + 1, total_states))

        H0 = H[0]
        e_vec_full = np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag
        if isinstance(e_vec_full, SparseArray):
            e_vec_full = e_vec_full.toarray()
            raise Exception(e_vec_full)
        for n, energies, overlaps, corrs, wfns in zip(
                states, all_energies, all_overlaps, all_corrs, all_wfns
        ):
            # taking advantage of mutability of the arrays here...

            # find the state index in the coupled subspace
            n_ind = np.where(total_coupled_space==n)[0][0]
            # generate the perturbation operator
            E0 = e_vec_full[n_ind]
            e_vec = e_vec_full - E0
            e_vec[n_ind] = 1
            pi = 1/e_vec
            pi[n_ind] = 0
            pi = SparseArray.from_diag(pi)

            energies[0] = E0
            overlaps[0] = 1
            corrs[0, n_ind] = 1

            def dot(a, b):
                if isinstance(a, np.ndarray):
                    return np.dot(a, b)
                else:
                    return a.dot(b)

            for k in range(1, order+1): # to actually go up to k
                #         En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
                Ek = (
                        H[k][n_ind, n_ind]
                        + sum(
                             dot(H[k-i][n_ind], corrs[i]) - energies[k-i]*overlaps[i]
                            for i in range(1, k)
                        )
                )
                energies[k] = Ek
                #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
                ok = -1/2 * sum(dot(corrs[i], corrs[k-i]) for i in range(1, k))
                overlaps[k] = ok
                #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=0...k-1) + <n^(0)|n^(k)> |n^(0)>
                corrs[k] = sum(
                    dot(pi, energies[k-i]*corrs[i] - dot(H[k-i], corrs[i]))
                    for i in range(0, k)
                )
                corrs[k][n_ind] = ok # pi (the perturbation operator) ensures it's zero before this

            # now we broadcast the corrections back up so that they're good for the _entire_ population of states...
            # this should be modified to work with sparse arrays in the future
            for wfn, cor in zip(wfns, corrs):
                if isinstance(total_coupled_space[0], (int, np.integer)):
                    cs = (total_coupled_space,)
                else:
                    cs = total_coupled_space
                wfn[cs] = cor

        corrs = PerturbationTheoryCorrections(
            {
                "states":states,
                 "coupled_states":total_coupled_space,
                 "total_states":total_states,
                 "degenerate_states":degenerate_states
             },
            {
                "energies":all_energies,
                "wavefunctions":all_wfns,
                "degenerate_transformation":None,

                "degenerate_energies": None
            },
            H
        )

        if degenerate_states is not None:
            # means we need to do the second level of corrections
            # we're going to need to do a full diagonalization, but we handle this
            # by pulling only the coupled elements down to a dense matrix and
            # diagonalizing there before returning the degenerate rotation as a sparse
            # matrix

            # first we figure out what the space of degenerate states is
            # has to be an NxN, so we figure out what the total space of potential degeneracies is
            degenerate_space = np.unique(np.concatenate(degenerate_states))
            deg_dim = len(degenerate_space)
            # then we have to work out how indices in the larger space map onto those in the smaller one
            remap = {i:k for k,i in enumerate(degenerate_space)}
            mapped_inds = (
                np.array([remap[i] for i in degenerate_states[0]]),
                np.array([remap[j] for j in degenerate_states[1]])
            )
            # at this point, we can fill H_deg in the basis of degenerate states
            H_deg = [np.zeros(1)] * len(deg_vals[1:])
            for i, v in enumerate(deg_vals[1:]):
                H_deg[i] = np.zeros((deg_dim, deg_dim))
                H_deg[i][mapped_inds[0], mapped_inds[1]] = v
                H_deg[i][mapped_inds[1], mapped_inds[0]] = v
            # now we need to transform from the basis of zero-order states to the basis of non-degenerate states
            # which we do by using our previously built PerturbationTheoryCorrections
            if logger is not None:
                logger.log_print(
                    "    generating representation of resonance terms in Hamiltonian"
                )
            H_deg = corrs.operator_representation(H_deg, subspace=degenerate_space)
            # now that we've corrected those elements, we add on the diagonal terms,
            # add things up, and diagonalize
            H_to_diag = np.sum(H_deg, axis=0)
            H_to_diag[np.diag_indices_from(H_to_diag)] = corrs.energies
            deg_engs, deg_transf = np.linalg.eigh(H_to_diag)

            # finally we re-sort so that the new wavefunctions look maximally like the non-degenerate ones
            deg_transf = deg_transf.T # easier for the moment...
            state_set = set(np.arange(len(deg_transf)))
            sort_transf = deg_transf.copy()

            for i in range(len(deg_transf)):
                max_ov = np.max(deg_transf[:, i]**2)
                ov_thresh = .5
                if max_ov < ov_thresh: # there must be a single mode that has more than 50% of the initial state character?
                    if logger is not None:
                        logger.log_print(
                             "    state {} is more than 50% mixed",
                            i
                        )
                #     raise PerturbationTheoryException("mode {} is has no contribution of greater than {}".format(
                #         i, ov_thresh
                #     ))
            sorting = [ -1 ] * len(deg_transf)
            for i in range(len(deg_transf)):
                o = np.argmax(abs(sort_transf[:, i]))
                sorting[i] = o
                sort_transf[o] = np.zeros(len(sort_transf))
            # sorting = [ np.argmax(abs(sort_transf[:, i])) for i in range(len(deg_transf)) ]
            # if len(sorting) != len(np.unique(sorting)):
            #     raise PerturbationTheoryException("After diagonalizing can't distinguish modes...")
            deg_engs = deg_engs[sorting,]
            deg_transf = deg_transf[sorting, :]

            corrs.degenerate_energies = deg_engs
            corrs.degenerate_transf = deg_transf

        return corrs

    def get_representations(self, coupled_states):
        """
        Returns the representations of the perturbations in the basis of coupled states

        :param coupled_states:
        :type coupled_states:
        :return:
        :rtype:
        """

        raise NotImplementedError("Implementation has changed and this method needs to catch up")
        total_states = int(np.prod(self.basis.quanta))
        if coupled_states is None:
            coupled_states = np.arange(total_states)  # expensive but...I guess what are you gonna do?
        elif not isinstance(coupled_states[0], int):
            coupled_states = self.basis.ravel_state_inds(coupled_states)

        m = coupled_states
        N = len(m)
        h_reps = self.perturbations
        H = [None] * len(h_reps)
        for i, h in enumerate(h_reps):
            sub = h[np.ix_(m, m)]
            if isinstance(sub, (int, np.integer, np.floating, float)):
                sub = np.full((N, N), sub)
            H[i] = sub  # type: np.ndarray
        return H

    def get_wavefunctions(self,
                          states,
                          coupled_states=None,
                          degeneracies=None,
                          order=2
                          ):
        """
        Gets a set of `PerturbationTheoryWavefunctions` from the perturbations defined by the Hamiltonian

        :param states: the states to get the index for, given either as indices or as a numbers of quanta
        :type states: Iterable[int] | Iterable[Iterable[int]]
        :param coupled_states: the list of states to explicitly allow to couple in
        :type coupled_states: Iterable[int] | Iterable[Iterable[int]]
        :param degeneracies: the pairs of states to be treated via degenerate perturbation theory
        :type degeneracies: (Iterable[int], Iterable[int])  | (Iterable[Iterable[int]], Iterable[Iterable[int]])
        :return: generated wave functions
        :rtype: PerturbationTheoryWavefunctions
        """

        from .Wavefunctions import PerturbationTheoryWavefunction, PerturbationTheoryWavefunctions

        if not isinstance(states[0], (int, np.integer)):
            states = self.basis.ravel_state_inds(states)
        total_states = int(np.prod(self.basis.quanta))
        if coupled_states is None or isinstance(coupled_states, (int, np.integer, float, np.floating)):
            state_modes = self.basis.unravel_state_inds(states)
            # pull the states that we really want to couple
            coupled_states = self.get_coupled_space(state_modes, order, freqs=self.modes.freqs, freq_threshold=coupled_states)
            coupled_states = [self.basis.ravel_state_inds(c) for c in coupled_states ]
        elif isinstance(coupled_states[0], int): # we use the same states in H1 and H2
            coupled_states = [coupled_states, coupled_states]
        elif len(coupled_states) > 2 and isinstance(coupled_states[0][0], int): # we use the same states in H1 and H2, but got them as proper modes
            coupled_states = self.basis.ravel_state_inds(coupled_states)
            coupled_states = [coupled_states, coupled_states]
        elif len(coupled_states) == 2 and not isinstance(coupled_states[0][0], int): # we got different states for H1 and H2, but they're as modes
            coupled_states = [ self.basis.ravel_state_inds(c) for c in coupled_states ]

        if (
                degeneracies is not None
                and not isinstance(degeneracies, (int, float, np.integer, np.floating))
                and not isinstance(degeneracies[0], (int, np.integer))
                and not isinstance(degeneracies[0][0], (int, np.integer))
        ):
            degeneracies = (
                self.basis.ravel_state_inds(degeneracies[0]),
                self.basis.ravel_state_inds(degeneracies[1])
            )

        corrs = self._get_VPT_corrections(
            self.perturbations,
            states,
            coupled_states,
            total_states,
            order,
            degenerate_states=degeneracies,
            logger=self.logger
            )

        return PerturbationTheoryWavefunctions(self.molecule, self.basis, corrs)