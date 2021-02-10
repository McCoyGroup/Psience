"""
Provides support for build perturbation theory Hamiltonians
"""

import numpy as np, itertools, time

from McUtils.Numputils import SparseArray, vec_outer
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer
from McUtils.Parallelizers import Parallelizer
from McUtils.Data import UnitsData

from ..Molecools import Molecule
from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, BraKetSpace, HarmonicOscillatorProductBasis

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
                 hamiltonians,
                 states,
                 coupled_states,
                 total_basis,
                 energy_corrs,
                 wfn_corrections,
                 degenerate_states=None,
                 degenerate_transformation=None,
                 degenerate_energies=None
                 ):
        """
        :param hamiltonians:
        :type hamiltonians: Iterable[SparseArray]
        :param states:
        :type states: BasisStateSpace
        :param coupled_states:
        :type coupled_states: BasisMultiStateSpace
        :param total_basis:
        :type total_basis: BasisMultiStateSpace
        :param energy_corrs:
        :type energy_corrs: np.ndarray
        :param wfn_corrections:
        :type wfn_corrections: Iterable[SparseArray]
        :param degenerate_states:
        :type degenerate_states: None | np.ndarray
        :param degenerate_transformation:
        :type degenerate_transformation: None | np.ndarray
        :param degenerate_energies:
        :type degenerate_energies: None | np.ndarray
        """
        self.hams = hamiltonians
        self.states = states
        self.coupled_states = coupled_states
        self.total_basis = total_basis
        self.energy_corrs = energy_corrs
        self.wfn_corrections = wfn_corrections
        self.degenerate_states = degenerate_states
        self.degenerate_transf = degenerate_transformation
        self.degenerate_energies = degenerate_energies

    @classmethod
    def from_dicts(cls,
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
        state_space = states['states']
        coupled_states = states['coupled_states']
        total_basis = states['total_states']
        energy_corrs = corrections['energies']
        wfn_corrections = corrections['wavefunctions']
        if 'degenerate_states' in states:
            degenerate_states = states['degenerate_states']
        else:
            degenerate_states = None

        if 'degenerate_transformation' in corrections:
            degenerate_transf = corrections['degenerate_transformation']
        else:
            degenerate_transf = None

        if 'degenerate_energies' in corrections:
            degenerate_energies = corrections['degenerate_energies']
        else:
            degenerate_energies = None

        return cls(
            hamiltonians,
            state_space,
            coupled_states,
            total_basis,
            energy_corrs,
            wfn_corrections,
            degenerate_states=degenerate_states,
            degenerate_transformation=degenerate_transf,
            degenerate_energies=degenerate_energies
        )

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

    def take_subspace(self, space):
        """
        Takes only those elements that are in space
        :param space:
        :type space:
        :return:
        :rtype:
        """

        new_states = self.states.find(space)
        return type(self)(
            self.hams,
            self.states.take_subspace(new_states),
            self.coupled_states.take_subspace(new_states),
            self.total_basis,
            self.energy_corrs[new_states],
            [w[new_states, :] for w in self.wfn_corrections],
            # not sure what to do with all this...
            degenerate_states=self.degenerate_states,
            degenerate_transformation=self.degenerate_transf,
            degenerate_energies=self.degenerate_energies
        )



    def operator_representation(self, operator_expansion, order=None, subspace=None):
        """
        Generates the representation of the operator in the basis of stored states

        :param operator_expansion: the expansion of the operator
        :type operator_expansion: Iterable[float] | Iterable[np.ndarray]
        :param order: the order of correction to go up to
        :type order: Iterable[float] | Iterable[np.ndarray]
        :param subspace: the subspace of terms in which the operator expansion is defined
        :type subspace: None | BasisStateSpace
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

        # we stopped supporting indexing based on the total set of inds...
        if subspace is None:
            wfn_corrs = self.wfn_corrections[:order]
        else:
            # need to make the subspace good with the subspace in which the corrections are defined...
            subspace_sel = self.coupled_states.find(subspace)
            wfn_corrs = []
            for k in range(order):
                wfn_corrs.append(self.wfn_corrections[k][:, subspace_sel])

        # generalizes the dot product so that we can use 0 as a special value...
        def dot(a, b):
            if isinstance(a, (int, np.integer, float, np.floating)) and a == 0:
                return 0

            if isinstance(a, np.ndarray):
                return np.dot(a, b)
            else:
                return a.dot(b)

        # does the dirty work of acutally applying the rep...
        reps = [np.zeros(1)] * order
        for k in range(order):
            op = None
            # apply each thing up to requested order...
            for a in range(k+1):
                for b in range(k-a+1):
                    c = k - (a + b)
                    rop = operator_expansion[c]
                    if isinstance(rop, (int, float, np.integer, np.floating)): # constant reps...
                        if rop != 0: # cheap easy check
                            subrep = rop * dot(wfn_corrs[a], wfn_corrs[b].T)
                            if op is None:
                                op = subrep
                            else:
                                op += subrep
                    else:
                        subrep = dot(dot(wfn_corrs[a], rop), wfn_corrs[b].T)
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
                 coriolis_coupling=True,
                 parallelizer=None,
                 log=None,
                 checkpoint=None
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
        :param parallelizer: parallelism manager
        :type parallelizer: Parallelizer
        :param log: log file or logger to write to
        :type log: str | Logger
        :param checkpoint: checkpoint file or checkpointer to store intermediate results
        :type checkpoint: str | Checkpointer
        """

        if isinstance(log, Logger):
            self.logger = log
        elif log is True:
            self.logger = Logger()
        elif log is False:
            self.logger = NullLogger()
        else:
            self.logger = Logger(log)

        if parallelizer is None:
            parallelizer = "VPT"
        self.parallelizer = parallelizer

        if checkpoint is None:
            checkpoint = NullCheckpointer(None)
        elif isinstance(checkpoint, str):
            checkpoint = Checkpointer.from_file(checkpoint)
        self.checkpointer = checkpoint

        if molecule is None:
            raise PerturbationTheoryException("{} requires a Molecule to do its dirty-work")
        # molecule = molecule.get_embedded_molecule()
        self.molecule = molecule
        if modes is None:
            modes = molecule.normal_modes
            # this basically presupposes we've got a held fe...might need to think
            # abotu the input format we want for this
            try:
                phases = molecule.get_fchk_normal_mode_rephasing()
            except NotImplementedError:
                pass
            else:
                modes = modes.rescale(phases)
        mode_n = modes.basis.matrix.shape[1] if mode_selection is None else len(mode_selection)
        self.mode_n = mode_n
        if n_quanta is None:
            n_quanta = 10 # dunno yet how I want to handle this since it should really be defined by the order of state requested...
        self.n_quanta = np.full((mode_n,), n_quanta) if isinstance(n_quanta, (int, np.int)) else tuple(n_quanta)
        self.modes = modes
        self.V_terms = PotentialTerms(self.molecule, modes=modes, mode_selection=mode_selection,
                                      logger=self.logger, parallelizer=self.parallelizer, checkpointer=self.checkpointer)
        self.G_terms = KineticTerms(self.molecule, modes=modes, mode_selection=mode_selection,
                                      logger=self.logger, parallelizer=self.parallelizer, checkpointer=self.checkpointer)
        if coriolis_coupling and (self.molecule.internal_coordinates is None):
            self.coriolis_terms = CoriolisTerm(self.molecule, modes=modes,
                                      logger=self.logger, parallelizer=self.parallelizer, checkpointer=self.checkpointer)
        else:
            self.coriolis_terms = None
        self.watson_term = PotentialLikeTerm(self.molecule, modes=modes,
                                      logger=self.logger, parallelizer=self.parallelizer, checkpointer=self.checkpointer)

        self._h0 = self._h1 = self._h2 = None

        self.basis = HarmonicOscillatorProductBasis(self.n_quanta)

        # from ..BasisReps import SimpleProductBasis, HarmonicOscillatorBasis
        # self.basis = SimpleProductBasis(HarmonicOscillatorBasis, self.n_quanta)

    @classmethod
    def from_fchk(cls,
                  file,
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
            if isinstance(self.basis, HarmonicOscillatorProductBasis):
                iphase = 1
            else:
                iphase = -1
            self._h0 = (
                    (iphase * 1 / 2) * self.basis.representation('p', 'p',
                                                                 coeffs=self.G_terms[0],
                                                                 logger=self.logger,
                                                                 parallelizer=self.parallelizer
                                                                 )
                    + 1 / 2 * self.basis.representation('x', 'x',
                                                        coeffs=self.V_terms[0],
                                                        logger=self.logger,
                                                        parallelizer=self.parallelizer
                                                        )
            )

        return self._h0

    @property
    def H1(self):
        """
        Provides the representation for H1 in this basis
        """
        if self._h1 is None:
            if isinstance(self.basis, HarmonicOscillatorProductBasis):
                iphase = 1
            else:
                iphase = -1
            self._h1 = (
                    (iphase * 1 / 2) * self.basis.representation('p', 'x', 'p',
                                                                 coeffs=self.G_terms[1],
                                                                 axes=[[0, 1, 2], [1, 0, 2]],
                                                                 logger=self.logger,
                                                                 parallelizer=self.parallelizer
                                                                 )
                    + 1 / 6 * self.basis.representation('x', 'x', 'x',
                                                        coeffs=self.V_terms[1],
                                                        logger=self.logger,
                                                        parallelizer=self.parallelizer
                                                        )
            )
        return self._h1

    @property
    def H2(self):
        """
        Provides the representation for H2 in this basis
        """
        if self._h2 is None:
            if isinstance(self.basis, HarmonicOscillatorProductBasis):
                iphase = 1
            else:
                iphase = -1
            self._h2 = (
                    (iphase * 1 / 4) * self.basis.representation('p', 'x', 'x', 'p',
                                                                 coeffs=self.G_terms[2],
                                                                 axes=[[0, 1, 2, 3], [2, 0, 1, 3]],
                                                                 logger=self.logger,
                                                                 parallelizer=self.parallelizer
                                                                 )
                    + 1 / 24 * self.basis.representation('x', 'x', 'x', 'x',
                                                         coeffs=self.V_terms[2],
                                                         logger=self.logger,
                                                         parallelizer=self.parallelizer
                                                         )
            )
            if self.coriolis_terms is not None:
                total_cor = self.coriolis_terms[0] + self.coriolis_terms[1] + self.coriolis_terms[2]
                self._h2 += iphase * self.basis.representation('x', 'p', 'x', 'p',
                                                               coeffs=total_cor,
                                                               logger=self.logger,
                                                               parallelizer=self.parallelizer
                                                               )
            else:
                self._h2 += 0 * self.basis.representation(coeffs=0,
                                                          logger=self.logger,
                                                          parallelizer=self.parallelizer
                                                          )

            self._h2 += 1 / 8 * self.basis.representation(coeffs=self.watson_term[0],
                                                          logger=self.logger,
                                                          parallelizer=self.parallelizer
                                                          )

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
                + v3[s, s, s] * v3[s, t, t] / (2 * w[s])
                + v3[t, t, t] * v3[t, s, s] / (2 * w[t])
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

    @classmethod
    def _get_Nielsen_xmat(cls, freqs, v3, v4, zeta, Be):
        ndim = len(freqs)
        x_mat_linear = np.array([
            cls._Nielsen_xss(s, freqs, v3, v4, zeta, Be, ndim) if s == t else
            cls._Nielsen_xst(s, t, freqs, v3, v4, zeta, Be, ndim)
            for s in range(ndim) for t in range(s, ndim)
        ]).T
        x_mat = np.zeros((3, ndim, ndim))
        ri, ci = np.triu_indices(ndim)
        x_mat[:, ri, ci] = x_mat_linear
        x_mat[:, ci, ri] = x_mat_linear
        return x_mat

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

        x_mat = cls._get_Nielsen_xmat(freqs, v3, v4, zeta, Be)

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        raise Exception(x_mat * h2w)

        states = np.array(states) + 1/2 # n+1/2 for harmonic vibrations

        x_mat = np.sum(x_mat, axis=0)
        e_harm = np.tensordot(freqs, states, axes=[0, 1])
        outer_states = vec_outer(states, states)
        # raise Exception(states, outer_states)
        e_anharm = np.tensordot(x_mat, outer_states, axes=[[0, 1], [1, 2]])

        return e_harm, e_anharm

    def get_Nielsen_xmatrix(self):

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # TODO: figure out WTF the units on this have to be...

        freqs = self.modes.freqs
        v3 = self.V_terms[1]
        v4 = self.V_terms[2]

        # raise Exception(np.round( 6 * v3 * h2w))

        zeta, Be = self.coriolis_terms.get_zetas_and_momi()

        x = self._get_Nielsen_xmat(freqs, v3, v4, zeta, Be)

        return x

    def get_Nielsen_energies(self, states):
        """

        :param states:
        :type states:
        :return:
        :rtype:
        """

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

    def get_coupled_space(self, states, order):
        """
        Returns the set of states that couple the given states up to the given order at each level of perturbation (beyond zero order).
        We keep track of how each individual state in states is transformed, as we only need to compute elements within those
        blocks, allowing for relatively dramatic speed-ups.

        :param state: the states of interest
        :type state: BasisStateSpace
        :param order: the order of perturbation theory we're doing
        :type order: int
        :param freqs: the zero-order frequencies in each vibrational mode being coupled
        :type freqs: Iterable[float]
        :param freq_threshold: the threshold for the maximum frequency difference between states to be considered
        :type freq_threshold: None | float
        :return: the sets of coupled states
        :rtype: tuple[BasisMultiStateSpace]
        """

        # the states that can be coupled through H1
        transitions_h1 = self.basis.selection_rules("x", "x", "x")
        h1_space = states.apply_selection_rules(
            transitions_h1,
            iterations=(order - 1)
        )

        # from second order corrections
        transitions_h2 = self.basis.selection_rules("x", "x", "x", "x")
        h2_space = states.apply_selection_rules(
            transitions_h2,
            iterations=(order - 1)
        )

        return h1_space, h2_space

    @classmethod
    def _get_VPT_representations(
            cls,
            h_reps,
            states,
            coupled_states,
            logger=None,
            parallelizer=None,
            freq_threshold=None
    ):
        """
        Gets the sparse representations of h_reps inside the basis of coupled states.
        Doesn't really need to be a classmethod but ah well.

        :param h_reps:
        :type h_reps: Iterable[Representation]
        :param states:
        :type states: BasisStateSpace
        :param coupled_states:
        :type coupled_states: Iterable[BasisStateSpace | BasisMultiStateSpace]
        :return:
        :rtype:
        """

        if logger is None:
            logger = NullLogger()

        par = Parallelizer.lookup(parallelizer)
        with par: # we put an outermost block here to just make sure everything is clean
            if len(coupled_states) != len(h_reps) - 1:
                raise ValueError("coupled states must be specified for all perturbations (got {}, expected {})".format(
                    len(coupled_states),
                    len(h_reps) - 1
                ))

            # determine the total coupled space
            coupled_spaces = []
            input_state_space = states

            with logger.block(tag="generating total space"):
                start = time.time()
                space_list = [input_state_space] + list(coupled_states)
                total_state_space = BasisMultiStateSpace(np.array(space_list,  dtype=object))
                flat_total_space = total_state_space.to_single().take_unique()
                diag_inds = BraKetSpace(flat_total_space, flat_total_space)
                N = len(flat_total_space)


                end = time.time()
                logger.log_print(
                    [
                        "total coupled space dimensions: {d}",
                        "took {t:.3f}s"
                     ],
                    d=N,
                    t=end-start
                )

            H = [np.zeros(1)] * len(h_reps)
            with logger.block(tag="getting H0"):
                start = time.time()
                logger.log_print(["calculating diagonal elements"])
                diag = h_reps[0][diag_inds]
                logger.log_print(["constructing sparse representation"])
                H[0] = SparseArray.from_diag(diag)
                end = time.time()
                logger.log_print("took {t:.3f}s", t=end-start)

            # print(flat_total_space.indices)
            for i, h in enumerate(h_reps[1:]):
                # calculate matrix elements in the coupled subspace
                cs = total_state_space[i+1]

                with logger.block(tag="getting H" + str(i + 1)):
                    m_pairs = cs.get_representation_brakets(freq_threshold=freq_threshold)

                    start = time.time()
                    if len(m_pairs) > 0:
                        logger.log_print(["coupled space dimension {d}"], d=len(m_pairs))
                        sub = h[m_pairs]
                        SparseArray.clear_ravel_caches()
                    else:
                        logger.log_print('no states to couple!')
                        sub = 0

                    logger.log_print("constructing sparse representation...")

                    if isinstance(sub, (int, np.integer, np.floating, float)):
                        if sub == 0:
                            sub = SparseArray.empty((N, N), dtype=float)
                        else:
                            raise ValueError("Using a constant shift of {} will force Hamiltonians to be dense...".format(sub))
                            sub = np.full((N, N), sub)
                    else:
                        # figure out the appropriate inds for this data in the sparse representation
                        row_inds = flat_total_space.find(m_pairs.bras)
                        col_inds = flat_total_space.find(m_pairs.kets)
                        # raise Exception([row_inds, col_inds,
                        #                  np.unique(m_pairs.bras.excitations, axis=0), np.average(sub)])

                        # upper triangle of indices
                        up_tri = np.array([row_inds, col_inds]).T
                        # lower triangle is made by transposition
                        low_tri = np.array([col_inds, row_inds]).T
                        # but now we need to remove the duplicates, because many sparse matrix implementations
                        # will sum up any repeated elements
                        full_inds = np.concatenate([up_tri, low_tri])
                        full_dat = np.concatenate([sub, sub])

                        _, idx = np.unique(full_inds, axis=0, return_index=True)
                        sidx = np.sort(idx)
                        full_inds = full_inds[sidx]
                        full_dat = full_dat[sidx]
                        sub = SparseArray((full_dat, full_inds.T), shape=(N, N))
                    end = time.time()
                    logger.log_print("took {t:.3f}s", t=end - start)

                H[i+1] = sub #type: np.ndarray
                # raise Exception(np.min(H[1].toarray()), np.max(H[1].toarray()))

                # import McUtils.Plots as plt
                # plt.ArrayPlot(H[1].toarray()).show()
                # raise Exception("...")

            return H, total_state_space

    @classmethod
    def _get_degenerate_rotation(cls, degenerate_space, corrs, H, logger):
        """
        Pull groups of non-degenerate state, use to rotate H to the non-degenerate basis,
        diagonalize
        :type corrs: PerturbationTheoryCorrections
        """

        with logger.block(tag="applying degenerate PT"):
            with logger.block(tag="states"):
                logger.log_print(
                    str(
                        corrs.states.take_states(degenerate_space).excitations
                    ).splitlines()
                )

            subcorrs = corrs.take_subspace(degenerate_space)
            # with logger.block(tag="first order corrections...."):
            #     logger.log_print(
            #         np.array2string(
            #             subcorrs.wfn_corrections[1].toarray()
            #         ).replace("\n", "").replace("]", "]\n").splitlines()
            #     )
            H_nd = sum([x.toarray() for x in subcorrs.operator_representation(H)])
            # logger.log_print("non-degenerate energies {e}",
            #                  e=np.diag(H_nd)*UnitsData.convert("Hartrees", "Wavenumbers"))

            with logger.block(tag="non-degenerate Hamiltonian"):
                logger.log_print(
                    str(
                        np.round(H_nd*UnitsData.convert("Hartrees", "Wavenumbers"))
                    ).splitlines()
                )

            deg_engs, deg_transf = np.linalg.eigh(H_nd)

            for i in range(len(deg_transf)):
                max_ov = np.max(deg_transf[:, i] ** 2)
                ov_thresh = .5
                if max_ov < ov_thresh:  # there must be a single mode that has more than 50% of the initial state character?
                    logger.log_print(
                        "    state {i} is more than 50% mixed",
                        i=i
                    )
                #     raise PerturbationTheoryException("mode {} is has no contribution of greater than {}".format(
                #         i, ov_thresh
                #     ))

            # we pick the terms with the max contribution from each input state
            # and zero out the contributions so that two states can't map
            # to the same input state
            sort_transf = np.abs(deg_transf.copy())
            sorting = [-1] * len(deg_transf)
            for i in range(len(deg_transf)):
                o = np.argmax(sort_transf[i, :])
                sorting[i] = o
                sort_transf[:, o] = 0.#np.zeros(len(sort_transf))

            with logger.block(tag='contributions'):
                logger.log_print(
                        str(np.round(100 * (deg_transf**2)).astype(int)).splitlines()
                )

            logger.log_print('sorting: {s}', s=sorting)

            # sorting = np.argsort(sorting)

            #
            # print(sorting)
            # # if len(sorting) != len(np.unique(sorting)):
            # #     raise PerturbationTheoryException("After diagonalizing can't distinguish modes...")
            deg_engs = deg_engs[sorting,]

            logger.log_print("degenerate energies {e}",
                             e=np.round(deg_engs * UnitsData.convert("Hartrees", "Wavenumbers")))

            deg_transf = deg_transf[:, sorting]

        return deg_engs, deg_transf



    @classmethod
    def _apply_degenerate_PT(cls,
                             H,
                             corrs,
                             degenerate_states,
                             logger=None,
                             checkpointer=None
                             ):
        """
        Applies degenerate perturbation theory by building a representation
        for the degenerate terms in the Hamiltonian.
        This is then diagonalized, allowing the degenerate states to be expressed
        in the basis of non-degenerate states

        :param H:
        :type H: Iterable[SparseArray]
        :param corrs: the standard PerturbationTheory Corrections object that comes out of the application of non-deg PT
        :type corrs: PerturbationTheoryCorrections
        :param degenerate_states: population of degenerate states
        :type degenerate_states:
        :param logger:
        :type logger: Logger
        :return:
        :rtype:
        """
        # we pull the states and total states from the corrections object

        total_state_space = corrs.states # type: BasisStateSpace

        # set up space to store the degenerate energies and rotations coming from the
        # sets of diagonalizations
        energies = np.zeros(len(total_state_space))
        base_energies = corrs.energies # for when we're not rotating

        # this will be built from a series of block-diagonal matrices
        # so we store the relevant values and indices to compose the SparseArray
        rotation_vals = []
        rotation_row_inds = []
        rotation_col_inds = []

        for group in degenerate_states:
            # we apply the degenerate PT on a group-by-group basis
            # by transforming the H reps into the non-degenerate basis
            deg_inds = total_state_space.find(group)
            if len(deg_inds) > 1:
                deg_engs, deg_rot = cls._get_degenerate_rotation(group, corrs, H, logger)
                energies[deg_inds] = deg_engs
                rotation_vals.append(deg_rot.flatten())
                deg_rows, deg_cols = np.array([p for p in itertools.product(deg_inds, deg_inds)]).T
                rotation_row_inds.append(deg_rows)
                rotation_col_inds.append(deg_cols)
            else:
                # we'll be a little bit inefficient for now and speed up later
                energies[deg_inds[0]] = base_energies[deg_inds[0]]
                rotation_vals.append([1.])
                rotation_row_inds.append(deg_inds)
                rotation_col_inds.append(deg_inds)

        rotation_vals = np.concatenate(rotation_vals)
        rotation_row_inds = np.concatenate(rotation_row_inds)
        rotation_col_inds = np.concatenate(rotation_col_inds)

        rotations = SparseArray(
            (
                rotation_vals,
                (
                    rotation_row_inds,
                    rotation_col_inds
                )
            ),
            shape=(len(energies), len(energies))
        )

        return energies, rotations

    @classmethod
    def _get_degenerate_transformation(cls, perturbations, deg_inds):
        """
        Diagonalizes the perturbations in the degenerate subspace to
        get a cleaner basis in which to do the perturbation theory.
        This comes from Sakurai.

        :param perturbations:
        :type perturbations:
        :param deg_inds:
        :type deg_inds:
        :return:
        :rtype: tuple[Iterable[SparseArray], SparseArray, SparseArray]
        """

        inds = np.array(list(itertools.product(deg_inds, deg_inds))).T
        subham = sum(
            np.reshape(H[inds], (len(deg_inds), len(deg_inds)))
            for H in perturbations[1:]
        )

        import McUtils.Plots as plt
        plt.ArrayPlot(subham).show()

        new_eng, deg_transf = np.linalg.eigh(subham)

        # we sort now to get the "best" mapping back onto the OG states

        sort_transf = np.abs(deg_transf.copy())
        sorting = [-1] * len(deg_transf)
        for i in range(len(deg_transf)):
            o = np.argmax(sort_transf[i, :])
            sorting[i] = o
            sort_transf[:, o] = 0.  # np.zeros(len(sort_transf))

        new_eng = new_eng[sorting]
        deg_transf = deg_transf.T[sorting]

        # with logger.block(tag='contributions'):
        #     logger.log_print(
        #         str(np.round(100 * (deg_transf ** 2)).astype(int)).splitlines()
        #     )

        # logger.log_print('sorting: {s}', s=sorting)

        return new_eng, deg_transf

    @classmethod
    def _apply_VPT_equations(cls,
                             H, order,
                             state_index,
                             degenerate_space_indices,
                             degenerate_energy,
                             degenerate_corrs,
                             total_state_space,
                             non_zero_cutoff=1.0e-14
                             ):
        """
        Does the dirty work of doing the VPT iterative equations.

        :return:
        :rtype:
        """

        n = state_index
        # this is relatively cheap to get and totally general
        e_vec_full = H[0].diag if isinstance(H[0], SparseArray) else np.diag(H[0])#zero_order_energies
        deg_inds = degenerate_space_indices

        energies = np.zeros((order+1,), dtype=float)
        overlaps = np.zeros((order+1,), dtype=float)
        corrs = np.zeros((order+1, len(total_state_space)), dtype=float)

        # find the state index in the coupled subspace
        n_ind = total_state_space.find(n)
        # generate the perturbation operator
        E0 = e_vec_full[n_ind]
        e_vec = e_vec_full - E0
        e_vec[deg_inds] = 1
        zero_checks = np.where(np.abs(e_vec) < non_zero_cutoff)[0]
        if len(zero_checks) > 0:
            bad_vec = np.concatenate([[E0], e_vec_full[zero_checks]])
            raise ValueError(
                "degeneracies encountered: state {} and {} other states are degenerate (average energy: {} stddev: {})".format(
                    n_ind,
                    len(zero_checks),
                    np.average(bad_vec),
                    np.std(bad_vec)
                ))
        pi = 1 / e_vec
        pi[deg_inds] = 0
        pi = SparseArray.from_diag(pi)

        energies[0] = E0
        if degenerate_energy is not None:
            energies[1] = degenerate_energy
        overlaps[0] = 1
        if degenerate_corrs is not None:
            corrs[0, deg_inds] = degenerate_corrs
        else:
            corrs[0, n_ind] = 1

        # generalizes the dot product so that we can use 0 as a special value...
        def dot(a, b):
            if isinstance(a, (int, np.integer, float, np.floating)) and a == 0:
                return 0

            if isinstance(a, np.ndarray):
                return np.dot(a, b)
            else:
                return a.dot(b)

        for k in range(1, order + 1):  # to actually go up to k
            #         En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
            if k >1 or degenerate_corrs is None:
                Ek = (
                        (H[k][n_ind, n_ind] if not isinstance(H[k], (int, np.integer, float, np.floating)) else 0.)
                        + sum(
                    dot(
                        H[k - i][n_ind] if not isinstance(H[k - i], (int, np.integer, float, np.floating)) else 0.,
                        corrs[i]
                    )
                    - energies[k - i] * overlaps[i]
                    for i in range(1, k)
                ))
                energies[k] = Ek
                #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
            ok = -1 / 2 * sum(dot(corrs[i], corrs[k - i]) for i in range(1, k))
            overlaps[k] = ok
            #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=0...k-1) + <n^(0)|n^(k)> |n^(0)>
            corrs[k] = sum(
                dot(pi, energies[k - i] * corrs[i] - dot(H[k - i], corrs[i]))
                for i in range(0, k)
            )
            corrs[k][n_ind] = ok  # pi (the perturbation operator) ensures it's zero before this

        return energies, overlaps, corrs

    @classmethod
    def _apply_nondegenerate_VPT(cls,
                                 perturbations,
                                 states,
                                 order,
                                 total_state_space,
                                 degenerate_states=None,
                                 degeneracy_mode=None,
                                 logger=None,
                                 checkpointer=None,
                                 non_zero_cutoff=1.0e-14
                                 ):
        """
        :param perturbations: the zero-order Hamiltonian and perturbations to be applied
        :type perturbations: Iterable[SparseArray]
        :param states: the set of states for which the perturbation theory will be applied
        :type states: BasisStateSpace
        :param order: the order of perturbation theory to apply
        :type order: int
        :param total_state_space:
        :type total_state_space: BasisMultiStateSpace
        :param degenerate_states: the groups of degenerate states
        :type degenerate_states: None | Iterable[BasisStateSpace]
        :param logger: logger object to print out debug info
        :type logger: Logger
        :param checkpointer: checkpointer to store checkpoint data
        :type checkpointer: Checkpointer
        :return:
        :rtype:
        """
        # We use the iterative equations
        #            En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
        #     <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
        #           |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=1...k-1) + <n^(0)|n^(k)> |n^(0)>
        #  where Pi_n is the perturbation operator [1/(E_m-E_n) for m!=n]

        flat_total_space = total_state_space.to_single().take_unique()
        N = len(flat_total_space)

        checkpointer['indices'] = total_state_space
        checkpointer['representations'] = perturbations

        all_energies = np.zeros((len(states), order + 1))
        all_overlaps = np.zeros((len(states), order + 1))
        all_corrs = np.zeros((len(states), order + 1, N))
        # all_wfns = [] #np.zeros((len(states), order + 1, total_states))

        H0 = perturbations[0]
        e_vec_full = np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag
        if isinstance(e_vec_full, SparseArray):
            e_vec_full = e_vec_full.toarray()
            # raise Exception(e_vec_full)

        with logger.block(tag="applying non-degenerate PT"):
            logger.log_print(
                [
                    "order: {o}",
                    "states: {n}"
                    ],
                o=order,
                n=len(states.indices)
            )
            start = time.time()

            # loop over the degenerate sets
            for deg_group in degenerate_states:
                # we use this to build a pertubation operator that removes
                # then entire set of degenerate states
                deg_inds = flat_total_space.find(deg_group)
                if hasattr(deg_group, 'indices'):
                    deg_group = deg_group.indices

                # if len(deg_inds) > 1:
                #     H, transf, transfinv, denseinv = cls._get_degenerate_transformation(perturbations, deg_inds)
                # else:
                #     H = perturbations

                if degeneracy_mode == "PT" and len(deg_inds) > 1:
                    deg_engs, deg_corrs = cls._get_degenerate_transformation(perturbations, deg_inds)
                    with logger.block(tag="handling degeneracies"):
                        logger.log_print(
                            [
                                "degenerate space: {s}",
                                "first-order energies: {e}"
                            ],
                            s=deg_inds,
                            e=deg_engs*UnitsData.convert("Hartrees", "Wavenumbers")
                        )
                        with logger.block(tag='zero-order states'):
                            logger.log_print(
                                np.array_str(deg_corrs).splitlines()
                            )

                else:
                    deg_engs = deg_corrs = [None]*len(deg_group)

                for n, de, dc in zip(deg_group, deg_engs, deg_corrs):
                    energies, overlaps, corrs = cls._apply_VPT_equations(perturbations,
                                                                         order, n,
                                                                         deg_inds,
                                                                         de, dc,
                                                                         # e_vec_full,
                                                                         flat_total_space,
                                                                         non_zero_cutoff=non_zero_cutoff
                                                                         )

                    res_index = states.find(n)
                    all_energies[res_index] = energies
                    all_corrs[res_index] = corrs
                    all_overlaps[res_index] = overlaps

            end = time.time()
            logger.log_print(
                "took {t}s",
                t=round(end - start, 3)
            )

        # now we recompute reduced state spaces for use in results processing
        # and we also convert the correction vectors to sparse representations
        tci = flat_total_space.indices
        N = len(tci)
        nstates = len(all_corrs)

        corr_inds = [[] for i in range(nstates)]
        corr_mats = [None] * (order + 1)
        # si = state_inds

        # checkpointer['corrections'] = all_corrs

        for o in range(order + 1):
            non_zeros = []
            for i, corr in enumerate(all_corrs):
                # we find the non-zero elements within the o level of correction for the ith state
                nonzi = np.where(np.abs(corr[o]) > non_zero_cutoff)[0]
                # print(nonzi)
                # then we pull these out
                vals = corr[o][nonzi,]
                # and we add the values and indices to the list
                non_zeros.append(
                    (
                        vals,
                        np.column_stack([
                            np.full(len(nonzi), i),
                            nonzi
                        ])
                    )
                )

                # and then we add the appropriate basis indices to the list of basis data
                wat = tci[nonzi,]
                corr_inds[i].append(wat)

            # now we build the full mat rep for this level of correction
            vals = np.concatenate([x[0] for x in non_zeros])
            inds = np.concatenate([x[1] for x in non_zeros], axis=0).T
            # print(inds, N)
            corr_mats[o] = SparseArray(
                (
                    vals,
                    inds
                ),
                shape=(nstates, N)
            )

        # now we build state reps from corr_inds
        for i, dat in enumerate(corr_inds):
            cat = np.concatenate(dat)
            _, upos = np.unique(cat, return_index=True)
            full_dat = cat[np.sort(upos)]
            corr_inds[i] = flat_total_space.take_states(full_dat) # BasisStateSpace(states.basis, full_dat, mode="indices")

        cs_states = SelectionRuleStateSpace(states, corr_inds, None)
        total_states = total_state_space
        corrs = PerturbationTheoryCorrections.from_dicts(
            {
                "states": states,
                "coupled_states": cs_states,
                "total_states": total_states,
                "degenerate_states": degenerate_states
            },
            {
                "energies": all_energies,
                "wavefunctions": corr_mats,
                "degenerate_transformation": None,
                "degenerate_energies": None
            },
            perturbations
        )

        # subcorrs = corrs.take_subspace(
        #     BasisStateSpace(
        #         states.basis,
        #         [[0, 0, 2], [0, 2, 0], [0, 1, 1]]
        #     )
        # )
        #
        # H_reps = subcorrs.operator_representation(H)
        # import json, os
        # with open(os.path.expanduser("~/Desktop/H_embed_subspace_test_deg.json"), "w+") as corr_dump:
        #     json.dump([x.toarray().tolist() for x in H_reps], corr_dump)
        #
        # with open(os.path.expanduser("~/Desktop/H_dump.json"), "w+") as corr_dump:
        #     json.dump([x.toarray().tolist() for x in H], corr_dump)

        return corrs

    @classmethod
    def _martin_test(cls, h_reps, states, threshold, total_coupled_space):
        """
        Applies the Martin Test to a set of states and perturbations to determine which resonances need to be
        treated variationally. Everything is done within the set of indices for the representations.

        :param h_reps: The representation matrices of the perturbations we're applying.
        :type h_reps: Iterable[np.ndarray | SparseArray]
        :param states: The indices of the states to which we're going apply to the Martin test.
        :type states: np.ndarray
        :param threshold: The threshold for what should be treated variationally (in the same energy units as the Hamiltonians)
        :type threshold: float
        :return: Pairs of coupled states
        :rtype: tuple[BasisStateSpace, BasisStateSpace]
        """

        raise NotImplementedError("This is fucked up :weep:; need to do full non-degenerate calc per pair of states")

        H0 = h_reps[0]
        H1 = h_reps[1]
        energies = np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag

        # the 'states' should already be indices within the space over which we do the H1 calculation
        # basically whichever states we need to treat as degenerate for
        state_energies = energies[states]
        diffs = state_energies[:, np.newaxis] - energies[np.newaxis, :]
        for n, s in enumerate(states):
            diffs[n, s] = 1

        deg_states = []
        for s in states:
            # pull the blocks out of H1 that correspond to each the `states` we fed in...
            H1_block = H1[s, :]
            if isinstance(H1_block, SparseArray):
                nzvals = H1_block.block_vals
                nzinds, _ = H1_block.block_inds
                H1_block = nzvals
                diffs = energies[s] - energies[nzinds] # do I need an abs ?
            else:
                # compute the energy differences
                diffs = energies[s] - energies # do I need an abs ?
                nzinds = np.arange(len(energies))

            s_pos = np.where(nzinds == s)[0]
            H1_block[s_pos] = 0
            diffs[s_pos] = 1

            anh_eff = (np.abs(H1_block) ** 4) / (diffs ** 3)
            big = np.where(np.abs(anh_eff) > threshold)[0]
            if len(big) > 0:
                deg_states.extend((s, nzinds[d]) for d in big)

        if len(deg_states) == 0:
            return None
        else:
            new_degs = np.array(deg_states).T

            # raise Exception(new_degs)

            # we now have indices inside the space of coupled states...
            # so now we need to broadcast these back into their indices in the overall basis of states
            tc_inds = total_coupled_space.indices
            basis = total_coupled_space.basis

            degs = (
                BasisStateSpace(basis, tc_inds[new_degs[0]], mode='indices'),
                BasisStateSpace(basis, tc_inds[new_degs[1]], mode='indices')
            )

            return degs

    @classmethod
    def _group_states_by_energy_cutoff(cls, H, states, total_state_space, cutoff):
        """
        :type H: Iterable[SparseArray]
        :type states: BasisStateSpace
        :type total_state_space: BasisMultiStateSpace
        :type cutoff: float
        :rtype: Iterable[BasisStateSpace]
        """
        # we look for states with energies within a range...
        # so initially we pull the sets of energies
        energies = (H[0].diag if isinstance(H[0], SparseArray) else np.diag(H[0]))
        degenerate_groups = []
        # then we look through the input states
        for n in total_state_space.find(states):
            # we only want to apply this once per degenerate group
            # NOTE: this is a path to subtlety, since
            #   if state a is within 50 cm^-1 of state b, and state b is within of c,
            #   you might argue a and c are degenerate
            #   we are wagering that states are distinct _enough_ such that this is not
            #   an issue, but if it is a different strategy will be required
            if all(n not in d for d in degenerate_groups):
                e_diffs = np.abs(energies - energies[n])
                inds = np.where(e_diffs < cutoff)[0]
                degenerate_groups.append(set(inds))
        # raise Exception(degenerate_groups)
        degenerate_groups = [total_state_space.take_subspace(np.array(list(d), dtype=int)) for d in degenerate_groups]
        return degenerate_groups

    @classmethod
    def _group_states_by_nt_spec(cls, H, states, total_state_space, q_vec):
        """
        :type H: Iterable[SparseArray]
        :type states: BasisStateSpace
        :type total_state_space: BasisMultiStateSpace
        :type cutoff: Iterable[int]
        :rtype: Iterable[BasisStateSpace]
        """
        # we build the total N_t to compare against once...
        tot_n_t = np.dot(total_state_space.excitations, q_vec)
        degenerate_groups = {}
        # then we look through the input states
        for vec in states.excitations:
            # base n_t
            n_t = np.dot(q_vec, vec)
            if n_t not in degenerate_groups:
                degenerate_groups[n_t] = np.where(tot_n_t == n_t)[0]
        degenerate_groups = [total_state_space.take_subspace(np.array(d)) for d in degenerate_groups.values()]
        return degenerate_groups

    @classmethod
    def _get_degenerate_state_sets(cls,
                                   H,
                                   states,
                                   degenerate_states,
                                   total_state_space,
                                   logger=None,
                                   martin_test=False
                                   ):
        """

        :param H:
        :type H:
        :param states: the states we're planning on getting PT results for
        :type states: BasisStateSpace
        :param degenerate_states: some kind of spec for defining how we're getting degeneracies
        :type degenerate_states:
        :param total_state_space: total state space
        :type total_state_space: BasisMultiStateSpace
        :param logger:
        :type logger: Logger
        :return:
        :rtype:
        """

        total_state_space = total_state_space.to_single().take_unique()

        if degenerate_states is not None:

            with logger.block(tag="getting degeneracies"):
                if isinstance(degenerate_states, dict):
                    if 'MartinTest' in degenerate_states:
                        martin_test = degenerate_states['MartinTest']
                    if 'states' in degenerate_states:
                        degenerate_states = degenerate_states['states']
                    elif 'NT' in degenerate_states:
                        logger.log_print(
                            "NT vector: {s}",
                            s=degenerate_states
                        )
                        degenerate_states = cls._group_states_by_nt_spec(H, states, total_state_space,
                                                                         degenerate_states['NT'])
                    elif 'energy_cutoff' in degenerate_states:
                        logger.log_print(
                            "energy cutoff: {s}",
                            s=degenerate_states
                        )
                        degenerate_states = cls._group_states_by_energy_cutoff(H, states, total_state_space,
                                                                               degenerate_states['energy_cutoff'])
                    elif 'generator' in degenerate_states:
                        logger.log_print(
                            "callable: {s}",
                            s=degenerate_states
                        )

                        # we assume we have some kind of callable
                        try:
                            degenerate_states = degenerate_states['generator'](H, states)
                        except (TypeError, ValueError):
                            pass

                        if not isinstance(degenerate_states[0], (BasisStateSpace, BasisMultiStateSpace)):
                            raise NotImplementedError("can't deal with non-BasisStateSpace specs for degeneracies")

                    else:
                        raise NotImplementedError("unsure what to do with degeneracy spec {}".format(degenerate_states))

                else:

                    def _is_degenerate_NT_spec(spec):
                        test1 = isinstance(spec, np.ndarray) and spec.dtype == np.dtype(int)
                        if test1:
                            return test1
                        else:
                            try:
                                it = iter(spec)
                            except TypeError:
                                return False
                            else:
                                return all(isinstance(x, int) for x in it)

                    def _is_degenerate_state_spec(spec):
                        test1 = all(isinstance(x, (BasisStateSpace, BasisMultiStateSpace)) for x in spec)
                        if test1:
                            return test1
                        else:
                            spec = np.asarray(spec)
                            return (
                                    spec.dtype == int
                                    and (spec.ndim == 2 or spec.ndim == 3)
                            )

                    # we dispatch on the types of degeneracy specs we support
                    if isinstance(degenerate_states, (int, np.integer, np.floating, float)):
                        logger.log_print(
                            "energy cutoff: {s}",
                            s=degenerate_states
                        )
                        degenerate_states = cls._group_states_by_energy_cutoff(H, states, total_state_space, degenerate_states)
                    elif _is_degenerate_NT_spec(degenerate_states):
                        logger.log_print(
                            "N_T vector: {s}",
                            s=degenerate_states
                        )
                        degenerate_states = cls._group_states_by_nt_spec(H, states, total_state_space, degenerate_states)
                    elif _is_degenerate_state_spec(degenerate_states):
                        if not isinstance(degenerate_states[0], (BasisStateSpace, BasisMultiStateSpace)):
                            degenerate_states = [BasisStateSpace(states.basis, x) for x in degenerate_states]
                    elif degenerate_states is not None:
                        logger.log_print(
                            "callable: {s}",
                            s=degenerate_states
                        )

                        # we assume we have some kind of callable
                        try:
                            degenerate_states = degenerate_states(H, states)
                        except (TypeError, ValueError):
                            pass

                        if not isinstance(degenerate_states[0], (BasisStateSpace, BasisMultiStateSpace)):
                            raise NotImplementedError("can't deal with non-BasisStateSpace specs for degeneracies")

                logger.log_print(
                    "degenerate state sets found {n}",
                    n=len([x for x in degenerate_states if len(x) > 1])
                )

                if martin_test:
                    # need to apply Martin test to every pair of states to figure out if they are truly
                    # going to be significantly affected by the near degeneracy
                    raise NotImplementedError("Don't have Martin test applying cleanly yet")
                    logger.log_print(
                        "applying Martin test with threshold {}",
                        thresh
                    )
                    degenerate_states = cls._martin_test(
                        H,
                        state_inds,  # state indices in the coupled_states
                        thresh,
                        total_state_space
                    )
                else:
                    logger.log_print(
                        "skipping Martin test"
                    )

        # build groups of degenerate states for use later
        if degenerate_states is None:
            groups = [[x] for x in states.indices] # we're gonna loop through this later so why not destructure now...
        else:
            groups = [[]] * len(degenerate_states)
            deg_sets = [set(d.indices) for d in degenerate_states]
            for x in states.indices:
                for i, d in enumerate(deg_sets):
                    if x in d:
                        if len(groups[i]) == 0:
                            groups[i] = []
                        groups[i].append(x)
                        break
                else:
                    groups.append([x])

        return groups

    @classmethod
    def _apply_VPT(cls,
                   H,
                   states,
                   coupled_states,
                   order,
                   total_state_space,
                   degenerate_states=None,
                   degeneracy_mode="PT",
                   logger=None,
                   checkpointer=None
                   ):
        """
        Takes the set of perturbations, `H`, state spaces, and
        applies both non-degenerate and (if requested) degenerate
        perturbation theory

        :param H:
        :type H:
        :param states:
        :type states:
        :param coupled_states: ...not entirely sure what this is anymore?
        :type coupled_states:
        :param order: order of PT to apply
        :type order: int
        :param total_state_space: total state space of interest
        :type total_state_space: BasisMultiStateSpace
        :param degenerate_states: sets of states to be treated through degenerate PT
        :type degenerate_states: None | Iterable[BasisStateSpace]
        :param logger:
        :type logger:
        :return:
        :rtype:
        """

        # state_inds = total_state_space.find(total_state_space[0].indices)

        degenerate_states = cls._get_degenerate_state_sets(
            H, states, degenerate_states, total_state_space,
            logger=logger
        )

        corrs = cls._apply_nondegenerate_VPT(
            H,
            states,
            order,
            total_state_space,
            degenerate_states=degenerate_states,
            degeneracy_mode=degeneracy_mode,
            logger=logger,
            checkpointer=checkpointer
        )

        if (
                degeneracy_mode == "basis" and
                degenerate_states is not None
                and any(len(x) > 1 for x in degenerate_states)
        ):

            # if logger is not None:
            #     logger.log_print(
            #         "handling degeneracies...",
            #     )
            deg_engs, deg_transf = cls._apply_degenerate_PT(
                H,
                corrs,
                degenerate_states,
                logger=logger,
                checkpointer=checkpointer
            )
            corrs.degenerate_energies = deg_engs
            corrs.degenerate_transf = deg_transf

        return corrs

    @classmethod
    def _get_VPT_corrections(
            cls,
            h_reps,
            states,
            coupled_states,
            # total_states,
            order,
            degenerate_states=None,
            logger=None,
            checkpointer=None
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

        H, total_state_space = cls._get_VPT_representations(h_reps, states, coupled_states, logger)

        corrs = cls._apply_VPT(
                   H,
                   states,
                   coupled_states,
                   order,
                   total_state_space,
                   degenerate_states=degenerate_states,
                   logger=logger,
                   checkpointer=checkpointer
                   )

        return corrs

    def _prep_coupled_states(self, states, coupled_states, order):
        """
        Preps coupled states as input for `get_wavefunctions` and `get_representations`

        :param states:
        :type states:
        :param coupled_states:
        :type coupled_states:
        :param order:
        :type order:
        :return:
        :rtype:
        """
        if coupled_states is None or isinstance(coupled_states, (int, np.integer, float, np.floating)):
            # pull the states that we really want to couple
            coupled_states = self.get_coupled_space(states,
                                                    order
                                                    # freqs=self.modes.freqs,
                                                    # freq_threshold=coupled_states
                                                    )

        elif isinstance(coupled_states, (BasisStateSpace, BasisMultiStateSpace)):
            # same spec for both perturbations
            coupled_states = [coupled_states, coupled_states]
        elif isinstance(coupled_states[0], (BasisStateSpace, BasisMultiStateSpace)):
            pass

        elif isinstance(coupled_states[0], (int, np.integer)): # got a single `BasisStateSpace` as indices
            coupled_states = BasisStateSpace(self.basis, coupled_states, mode="indices")
            coupled_states = [coupled_states, coupled_states]

        elif isinstance(coupled_states[0][0], (int, np.integer)): # got a single `BasisStateSpace` as excitations
            coupled_states = BasisStateSpace(self.basis, coupled_states, mode="excitations")
            coupled_states = [coupled_states, coupled_states]

        else:
            raise ValueError("Not sure what to do with coupled space spec {}".format(
                coupled_states
            ))

        return coupled_states

    def _prep_degeneracies_spec(self, degeneracies):
        if (
                degeneracies is not None
                and not isinstance(degeneracies, (int, float, np.integer, np.floating))
        ):
            if isinstance(degeneracies[0], (int, np.integer)):
                degs = BasisStateSpace(self.basis, degeneracies, mode="indices")
                degeneracies = (degs, degs)
            elif isinstance(degeneracies[0][0], (int, np.integer)):
                degs = BasisStateSpace(self.basis, degeneracies, mode="excitations")
                degeneracies = (degs, degs)
            else:
                degeneracies = (
                                   BasisStateSpace(self.basis, degeneracies[0]),
                                   BasisStateSpace(self.basis, degeneracies[1])
                )

        return degeneracies

    def get_representations(self,
                            states,
                            coupled_states=None,
                            degeneracies=None,
                            order=2
                            ):
        """
        Returns the representations of the perturbations in the basis of coupled states

        :param coupled_states:
        :type coupled_states:
        :return:
        :rtype:
        """

        with self.logger.block(tag='Computing PT corrections:'):
            with self.logger.block(tag='getting coupled states'):
                start = time.time()

                states, coupled_states, degeneracies = self.get_input_state_spaces(states, coupled_states, degeneracies)

                end = time.time()
                self.logger.log_print("took {}s...", round(end - start, 3))

            H, tot_space = self._get_VPT_representations(
                self.perturbations,
                states,
                coupled_states,
                logger=self.logger
            )

        return H

    def get_input_state_spaces(self,
                               states,
                               coupled_states=None,
                               degeneracies=None,
                               order=2
                               ):
        """
        Converts the input state specs into proper `BasisStateSpace` specs that
        will directly feed into the code

        :param states:
        :type states:
        :param coupled_states:
        :type coupled_states:
        :param degeneracies:
        :type degeneracies:
        :return:
        :rtype:
        """

        # raise Exception(states)
        # need to rewrite this to work better with BasisStateSpace
        if not isinstance(states, BasisStateSpace):
            states = BasisStateSpace(self.basis, states)

        coupled_states = self._prep_coupled_states(states, coupled_states, order)

        # degeneracies = self._prep_degeneracies_spec(degeneracies)

        return states, coupled_states, degeneracies

    def get_wavefunctions(self,
                          states,
                          coupled_states=None,
                          degeneracies=None,
                          order=2
                          ):
        """
        Gets a set of `PerturbationTheoryWavefunctions` from the perturbations defined by the Hamiltonian

        :param states: the states to get the index for, given either as indices or as a numbers of quanta
        :type states: BasisStateSpace | Iterable[int] | Iterable[Iterable[int]]
        :param coupled_states: the list of states to explicitly allow to couple in
        :type coupled_states: BasisStateSpace | Iterable[int] | Iterable[Iterable[int]]
        :param degeneracies: the pairs of states to be treated via degenerate perturbation theory
        :type degeneracies: (Iterable[int], Iterable[int])  | (Iterable[Iterable[int]], Iterable[Iterable[int]])
        :return: generated wave functions
        :rtype: PerturbationTheoryWavefunctions
        """

        with self.checkpointer:

            from .Wavefunctions import PerturbationTheoryWavefunctions

            with self.logger.block(tag='Computing PT corrections:'):
                with self.logger.block(tag='getting coupled states'):
                    start = time.time()
                    # raise Exception("???", states)
                    states, coupled_states, degeneracies = self.get_input_state_spaces(states, coupled_states, degeneracies)

                    end = time.time()
                    self.logger.log_print("took {t}s...", t=round(end - start, 3))

                h_reps = self.perturbations
                if self.logger is not None:
                    bs = []
                    for b in coupled_states:
                        bs.append(b.unique_len)
                    self.logger.log_print(
                        [
                            "perturbations: {pert_num}",
                            "order: {ord}",
                            "states: {state_num}",
                            "basis sizes {basis_size}"
                        ],
                        pert_num=len(h_reps) - 1,
                        ord=order,
                        state_num=len(states),
                        basis_size=bs
                    )

                corrs = self._get_VPT_corrections(
                    h_reps,
                    states,
                    coupled_states,
                    order,
                    degenerate_states=degeneracies,
                    logger=self.logger,
                    checkpointer=self.checkpointer
                    )

        return PerturbationTheoryWavefunctions(self.molecule, self.basis, corrs, logger=self.logger)

    @classmethod
    def _invert_action_expansion_tensors(cls,
                                         states,
                                         energies,
                                         order
                                         ):
        """
        Obtains expansions of the relevant tensors in terms of their classical actions.
        Only applicable to a Harmonic PT approach, really, but quite useful.

        :param states: states up to `n` quanta of excitation, where n=order of expansion
        :type states: BasisStateSpace
        :param energies: energies from PT for the states
        :type energies: np.ndarray
        :param order: the order of perturbation theory we applied
        :type order: int
        :return:
        :rtype: list[np.array | float]
        """

        nmodes = states.ndim
        exc = states.excitations

        c_mat = np.zeros((len(states), len(states)), dtype=float)  # to invert

        col = 0
        blocks = []  # to more easily recompose the tensors later
        nterms = 1 + order // 2 # second order should be [2, 1], 4th order should be [3, 2, 1], 6th should be [4, 3, 2, 1]
        for k in range(nterms, 0, -1):
            ninds = []
            # generate the index tensors to loop over
            inds = np.indices((nmodes,) * k)
            inds = inds.transpose(tuple(range(1, k + 1)) + (0,))
            inds = np.reshape(inds, (-1, k))
            for perm in inds:
                if (np.sort(perm) != perm).any():
                    continue # only want the _unique_ permutations
                # generate the action coefficients
                coeffs = np.prod(exc[:, perm] + 1 / 2, axis=1)
                c_mat[:, col] = coeffs
                col += 1
                ninds.append(perm)
            blocks.append(ninds)
        # finally we add in the coefficient from k=0
        c_mat[:, col] = 1

        # get the solutions to the linear equation
        tensor_terms = np.linalg.solve(c_mat, energies)

        # reconstruct the tensors
        tens = [np.zeros(1)] * (nterms + 1)
        where_am_i = 0
        for i, b in enumerate(blocks):
            s = where_am_i
            nb = len(b)
            vec = tensor_terms[s:s+nb]
            k = nterms - i
            term = np.zeros((nmodes,) * k)
            bi = tuple(np.transpose(np.array(b)))
            term[bi] = vec
            where_am_i += nb
            tens[k] = term

        tens[0] = tensor_terms[-1]

        return tens

    def get_action_expansion(self, order=2):
        """
        Gets the expansion of the energies in terms of Miller's "classical actions" by
        doing just enough PT to invert the matrix

        :param order:
        :type order:
        :return:
        :rtype:
        """

        ndim = len(self.n_quanta)
        nterms = 1 + order // 2

        # # clumsy buy w/e it works for now
        # def get_states(n_quanta, n_modes, max_quanta=None):
        #     import itertools as ip
        #
        #     if max_quanta is None:
        #         max_quanta = n_quanta
        #     return tuple(sorted(
        #         [p for p in ip.product(*(range(n_quanta + 1) for i in range(n_modes))) if
        #          all(x <= max_quanta for x in p) and sum(p) <= n_quanta],
        #         key=lambda p: (
        #                 sum(p)
        #                 + sum(1 for v in p if v != 0) * n_quanta ** (-1)
        #                 + sum(v * n_quanta ** (-i - 2) for i, v in enumerate(p))
        #         )
        #     ))

        states = BasisStateSpace.from_quanta(HarmonicOscillatorProductBasis(ndim), range(nterms)).excitations

        wfns = self.get_wavefunctions(states)

        expansion = self._invert_action_expansion_tensors(wfns.corrs.states, wfns.energies, order)

        return expansion, wfns


    def get_breakdown(self,
                      states,
                      coupled_states=None,
                      degeneracies=None,
                      order=2
                      ):

            from collections import OrderedDict
            from .Wavefunctions import PerturbationTheoryWavefunctions

            if self.logger is not None:
                self.logger.log_print("Computing PT breakdown:", padding="")
                start = time.time()
                self.logger.log_print("getting coupled states...")

            states, coupled_states, degeneracies = self.get_input_state_spaces(states, coupled_states, degeneracies)

            if self.logger is not None:
                end = time.time()
                self.logger.log_print("took {}s...", round(end - start, 3))

            h_reps = self.perturbations
            if self.logger is not None:
                bs = []
                for b in coupled_states:
                    bs.append(len(b))
                self.logger.log_print(
                    [
                        "perturbations: {pert_num}",
                        "order: {ord}",
                        "states: {state_num}",
                        "basis sizes {basis_size}"
                    ],
                    pert_num=len(h_reps) - 1,
                    ord=order,
                    state_num=len(states),
                    basis_size=bs
                )

            H, total_state_space = self._get_VPT_representations(h_reps, states, coupled_states, self.logger)

            specs = OrderedDict((
                ("Harmonic",   (True, False, False)),
                ("Cubic",      (True, True,  False)),
                ("Quartic",    (True, False, True)),
                ("Full",       (True, True,  True))
            ))

            for k in specs:
                this_h = [H[i] if len(H) > i and s else 0 for i, s in enumerate(specs[k])]
                if self.logger is not None:
                    self.logger.log_print(
                        [
                            "getting breakdown for {} terms...",
                            "(non-zero terms {})"
                            ],
                        k,
                        len(this_h) - this_h.count(0)
                    )
                corrs = self._apply_VPT(
                    this_h,
                    states,
                    coupled_states,
                    order,
                    total_state_space,
                    degenerate_states=degeneracies,
                    logger=self.logger
                )

                wfns = PerturbationTheoryWavefunctions(self.molecule, self.basis, corrs, logger=self.logger)
                specs[k] = wfns

            return specs