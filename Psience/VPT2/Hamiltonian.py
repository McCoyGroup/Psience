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
from .Solver import PerturbationTheorySolver, PerturbationTheoryCorrections

__all__ = [
    'PerturbationTheoryHamiltonian',
    'PerturbationTheoryCorrections'
]

__reload_hook__ = [ "..BasisReps" , '.Terms', ".Solver", "..Molecools", ]

class PerturbationTheoryHamiltonian:
    """
    Represents the main Hamiltonian used in the perturbation theory calculation.
    Uses a harmonic oscillator basis for representing H0, H1, and H2 (and only goes up to H2 for now).
    Will before too long be split into a PerturbationTheoryHandler and a PerturbationTheoryHamiltonian
    where the Hamiltonian just implements the split of the perturbation and the Handler manages the equations.
    """

    def __init__(self,
                 molecule=None,
                 n_quanta=None,
                 modes=None,
                 mode_selection=None,
                 potential_derivatives=None,
                 coriolis_coupling=True,
                 parallelizer=None,
                 log=None,
                 checkpoint=None,
                 operator_chunk_size=None
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
                if phases is not None:
                    modes = modes.rescale(phases)
        if mode_selection is not None:
            mode_selection = tuple(mode_selection)
        mode_n = modes.basis.matrix.shape[1] if mode_selection is None else len(mode_selection)
        self.mode_n = mode_n
        if n_quanta is None:
            n_quanta = 10 # dunno yet how I want to handle this since it should really be defined by the order of state requested...
        self.n_quanta = np.full((mode_n,), n_quanta) if isinstance(n_quanta, (int, np.int)) else tuple(n_quanta)
        self.modes = modes
        self.V_terms = PotentialTerms(self.molecule, modes=modes, mode_selection=mode_selection,
                                      potential_derivatives=potential_derivatives,
                                      logger=self.logger, parallelizer=self.parallelizer, checkpointer=self.checkpointer)
        self.G_terms = KineticTerms(self.molecule, modes=modes, mode_selection=mode_selection,
                                    logger=self.logger, parallelizer=self.parallelizer, checkpointer=self.checkpointer)
        if coriolis_coupling and (self.molecule.internal_coordinates is None):
            self.coriolis_terms = CoriolisTerm(self.molecule, modes=modes, mode_selection=mode_selection,
                                      logger=self.logger, parallelizer=self.parallelizer, checkpointer=self.checkpointer)
        else:
            self.coriolis_terms = None
        self.watson_term = PotentialLikeTerm(self.molecule, modes=modes, mode_selection=mode_selection,
                                      logger=self.logger, parallelizer=self.parallelizer, checkpointer=self.checkpointer)

        self._h0 = self._h1 = self._h2 = None

        self.basis = HarmonicOscillatorProductBasis(self.n_quanta)

        self.operator_settings = {
            'chunk_size': operator_chunk_size,
            'logger': self.logger,
            'parallelizer': self.parallelizer
        }

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

        molecule = Molecule.from_file(file,
                                      zmatrix=internals,
                                      mode='fchk'
                                      )
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
                                                                 **self.operator_settings
                                                                 )
                    + 1 / 2 * self.basis.representation('x', 'x',
                                                        coeffs=self.V_terms[0],
                                                        **self.operator_settings
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
                                                                 **self.operator_settings
                                                                 )
                    + 1 / 6 * self.basis.representation('x', 'x', 'x',
                                                        coeffs=self.V_terms[1],
                                                        **self.operator_settings
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
                                                                 **self.operator_settings
                                                                 )
                    + 1 / 24 * self.basis.representation('x', 'x', 'x', 'x',
                                                         coeffs=self.V_terms[2],
                                                         **self.operator_settings
                                                         )
            )
            if self.coriolis_terms is not None:
                total_cor = self.coriolis_terms[0] + self.coriolis_terms[1] + self.coriolis_terms[2]
                self._h2 += iphase * self.basis.representation('x', 'p', 'x', 'p',
                                                               coeffs=total_cor,
                                                               **self.operator_settings
                                                               )
            else:
                self._h2 += 0 * self.basis.representation(coeffs=0,
                                                          **self.operator_settings
                                                          )

            self._h2 += 1 / 8 * self.basis.representation(coeffs=self.watson_term[0],
                                                          **self.operator_settings
                                                          )

        return self._h2

    @property
    def perturbations(self):
        return (self.H0, self.H1, self.H2)

    #region Nielsen energies

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

    #endregion Nielsen energies

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

        nits = order - 1
        if nits >= 0:
            # the states that can be coupled through H1
            self.logger.log_print('getting states coupled through H^(1)')
            transitions_h1 = self.basis.selection_rules("x", "x", "x")
            h1_space = states.apply_selection_rules(
                transitions_h1,
                iterations=(order - 1)
            )

            # from second order corrections
            self.logger.log_print('getting states coupled through H^(2)')
            transitions_h2 = self.basis.selection_rules("x", "x", "x", "x")
            h2_space = states.apply_selection_rules(
                transitions_h2,
                iterations=(order - 1)
            )
        else:
            h1_space = states.take_states([])
            h2_space = states.take_states([])

        return h1_space, h2_space

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
            coupled_states = self.get_coupled_space(states, order
                                                    # freqs=self.modes.freqs,
                                                    # freq_threshold=coupled_states
                                                    )

        elif isinstance(coupled_states, (BasisStateSpace, BasisMultiStateSpace)):
            # same spec for both perturbations
            coupled_states = [coupled_states, coupled_states]
        elif isinstance(coupled_states[0], (BasisStateSpace, BasisMultiStateSpace)):
            pass

        elif isinstance(coupled_states[0], (int, np.integer)):  # got a single `BasisStateSpace` as indices
            coupled_states = BasisStateSpace(self.basis, coupled_states, mode="indices")
            coupled_states = [coupled_states, coupled_states]

        elif isinstance(coupled_states[0][0], (int, np.integer)):  # got a single `BasisStateSpace` as excitations
            coupled_states = BasisStateSpace(self.basis, coupled_states, mode="excitations")
            coupled_states = [coupled_states, coupled_states]

        else:
            raise ValueError("Not sure what to do with coupled space spec {}".format(
                coupled_states
            ))

        return coupled_states

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

                self.logger.log_print(
                    [
                        "states: {state_num}",
                        "order: {ord}",
                    ],
                    ord=order,
                    state_num=len(states),
                )

                h_reps = self.perturbations
                self.logger.log_print(
                        "perturbations: {pert_num}",
                    pert_num=len(h_reps) - 1,
                )

                with self.logger.block(tag='getting coupled states'):
                    start = time.time()
                    states, coupled_states, degeneracies = self.get_input_state_spaces(states, coupled_states, degeneracies,
                                                                                       order=order)
                    end = time.time()
                    self.logger.log_print("took {t}s...", t=round(end - start, 3))

                bs = []
                for b in coupled_states:
                    bs.append(b.unique_len)
                self.logger.log_print(
                    [
                        "basis sizes {basis_size}"
                    ],
                    basis_size=bs
                )

                solver = PerturbationTheorySolver(h_reps, states, coupled_states, order,
                                                 degeneracy_spec=degeneracies,
                                                 logger=self.logger,
                                                 checkpointer=self.checkpointer,
                                                 allow_sakurai_degs=True,
                                                 allow_post_PT_calc=True
                                                 )
                corrs = solver.apply_VPT()

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

        raise NotImplementedError("changed up form of Solver and need to include these changes here too")

        from collections import OrderedDict
        from .Wavefunctions import PerturbationTheoryWavefunctions

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