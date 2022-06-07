"""
Provides support for build perturbation theory Hamiltonians
"""

import numpy as np, itertools, time

from McUtils.Numputils import SparseArray, vec_outer
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer, ParameterManager
from McUtils.Parallelizers import Parallelizer
from McUtils.Combinatorics import CompleteSymmetricGroupSpace

from ..Molecools import Molecule
from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, BraKetSpace, HarmonicOscillatorProductBasis

from .Common import PerturbationTheoryException
from .Terms import PotentialTerms, KineticTerms, CoriolisTerm, PotentialLikeTerm
from .Solver import PerturbationTheorySolver, PerturbationTheoryCorrections
from .Wavefunctions import PerturbationTheoryWavefunctions

__all__ = [
    'PerturbationTheoryHamiltonian',
    'PerturbationTheoryCorrections'
]

__reload_hook__ = [ "..BasisReps" , "..Molecools", '.Terms', ".Solver", ".Wavefunctions" ]

class PerturbationTheoryHamiltonian:
    """
    Represents the main Hamiltonian used in the perturbation theory calculation.
    Uses a harmonic oscillator basis for representing H0, H1, H2 etc.
    The PT process is split into a PerturbationTheorySolver and a PerturbationTheoryHamiltonian
    where the Hamiltonian just implements the split of the perturbation and the Solver manages the equations.
    """

    def __init__(self,
                 molecule=None,
                 n_quanta=None,
                 modes=None,
                 mode_selection=None,
                 potential_derivatives=None,
                 include_coriolis_coupling=True,
                 include_pseudopotential=True,
                 potential_terms=None,
                 allow_higher_potential_terms=False,
                 kinetic_terms=None,
                 coriolis_terms=None,
                 pseudopotential_terms=None,
                 selection_rules=None,
                 operator_chunk_size=None,
                 logger=None,
                 checkpoint=None,
                 results=None,
                 parallelizer=None,
                 **expansion_options
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
        :param include_coriolis_coupling: whether to add coriolis coupling if not in internals
        :type include_coriolis_coupling: bool
        :param parallelizer: parallelism manager
        :type parallelizer: Parallelizer
        :param logger: log file or logger to write to
        :type logger: str | Logger
        :param checkpoint: checkpoint file or checkpointer to store intermediate results
        :type checkpoint: str | Checkpointer
        """

        if isinstance(logger, Logger):
            self.logger = logger
        elif logger is True:
            self.logger = Logger()
        elif logger is False:
            self.logger = NullLogger()
        else:
            self.logger = Logger(logger)

        if parallelizer is None:
            parallelizer = "VPT"
        self.parallelizer = parallelizer

        self.results = Checkpointer.build_canonical(results) if results is not None else results
        self.checkpointer = Checkpointer.build_canonical(checkpoint)

        if molecule is None:
            raise PerturbationTheoryException("{} requires a Molecule to do its dirty-work")
        # molecule = molecule.get_embedded_molecule()
        self.molecule = molecule
        if modes is None:
            modes = molecule.normal_modes.modes
        if mode_selection is not None:
            mode_selection = tuple(mode_selection)
        mode_n = modes.basis.matrix.shape[1] if mode_selection is None else len(mode_selection)
        self.mode_n = mode_n
        if n_quanta is None:
            # This is a basically a historical option. We keep it but there's really no reason.
            n_quanta = 15 # dunno yet how I want to handle this since it should really be defined by the order of state requested...
        self.n_quanta = np.full((mode_n,), n_quanta) if isinstance(n_quanta, (int, np.int)) else tuple(n_quanta)
        self.modes = modes
        self.mode_selection = mode_selection

        expansion_options['logger'] = self.logger
        expansion_options['checkpointer'] = self.results if self.results is not None else self.checkpointer
        expansion_options['parallelizer'] = self.parallelizer
        expansion_params = ParameterManager(expansion_options)

        self._input_potential = potential_terms
        self.V_terms = PotentialTerms(self.molecule,
                                      modes=modes, mode_selection=mode_selection,
                                      potential_derivatives=potential_derivatives,
                                      allow_higher_potential_terms=allow_higher_potential_terms,
                                      **expansion_params.filter(PotentialTerms)
                                      )

        self._input_kinetic = kinetic_terms
        self.G_terms = KineticTerms(self.molecule, modes=modes, mode_selection=mode_selection,
                                    **expansion_params.filter(KineticTerms)
                                    )

        self._input_coriolis = coriolis_terms
        if (
                include_coriolis_coupling and
                (self.molecule.internal_coordinates is None or any(
                    k in expansion_options and expansion_options[k] for k in [
                        'use_cartesian_kinetic_energy',
                        'backpropagate_internals'
                        ]
                ))
        ):
            self.coriolis_terms = CoriolisTerm(self.molecule, modes=modes, mode_selection=mode_selection,
                                               **expansion_params.filter(CoriolisTerm)
                                               )
        else:
            # raise Exception(self.molecule.internal_coordinates)
            self.coriolis_terms = None

        self._input_pseudopotential = pseudopotential_terms
        if include_pseudopotential:
            self.pseudopotential_term = PotentialLikeTerm(self.molecule, modes=modes, mode_selection=mode_selection,
                                                          **expansion_params.filter(PotentialLikeTerm)
                                                          )
        else:
            self.pseudopotential_term = None

        self._expansions = []
        # self._h0 = self._h1 = self._h2 = None
        self._selection_rules = selection_rules

        self.basis = HarmonicOscillatorProductBasis(self.n_quanta)

        self.expansion_options = expansion_options

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
        o = 0
        if len(self._expansions) < o+1:
            self._expansions += [None] * (o+1 - len(self._expansions))

        if self._expansions[0] is None:
            if isinstance(self.basis, HarmonicOscillatorProductBasis):
                iphase = 1
            else:
                iphase = -1

            T0 = self._input_kinetic[0] if self._input_kinetic is not None and len(self._input_kinetic) > 0 else None
            if T0 is None:
                T0 = self.G_terms[0]
            V0 = self._input_potential[0] if self._input_potential is not None and len(self._input_potential) > 0 else None
            if V0 is None:
                V0 = self.V_terms[0]
            self._expansions[0] = (
                    (iphase * 1 / 2) * self.basis.representation('p', 'p',
                                                                 coeffs=T0,
                                                                 name='T(0)',
                                                                 **self.operator_settings
                                                                 )
                    + 1 / 2 * self.basis.representation('x', 'x',
                                                        coeffs=V0,
                                                        name='V(0)',
                                                        **self.operator_settings
                                                        )
            )
            self._expansions[0].name = "H(0)"

        return self._expansions[0]

    @property
    def H1(self):
        """
        Provides the representation for H1 in this basis
        """

        o = 1
        if len(self._expansions) < o + 1:
            self._expansions += [None] * (o + 1 - len(self._expansions))

        if self._expansions[1] is None:
            if isinstance(self.basis, HarmonicOscillatorProductBasis):
                iphase = 1
            else:
                iphase = -1

            T = self._input_kinetic[1] if self._input_kinetic is not None and len(self._input_kinetic) > 1 else None
            if T is None:
                T = self.G_terms[1]
            V = self._input_potential[1] if self._input_potential is not None and len(self._input_potential) > 1 else None
            if V is None:
                V = self.V_terms[1]
            self._expansions[1] = (
                    (iphase * 1 / 2) * self.basis.representation('p', 'x', 'p',
                                                                 coeffs=T,
                                                                 name='T(1)',
                                                                 axes=[[0, 1, 2], [1, 0, 2]],
                                                                 **self.operator_settings
                                                                 )
                    + 1 / 6 * self.basis.representation('x', 'x', 'x',
                                                        coeffs=V,
                                                        name='V(1)',
                                                        **self.operator_settings
                                                        )
            )


            if self._selection_rules is not None and len(self._selection_rules) > 0:
                self._expansions[1].selection_rules = self._selection_rules[0]

            self._expansions[1].name = "H(1)"

        return self._expansions[1]

    @property
    def H2(self):
        """
        Provides the representation for H2 in this basis
        """

        o = 2
        if len(self._expansions) < o + 1:
            self._expansions += [None] * (o + 1 - len(self._expansions))

        if self._expansions[2] is None:
            if isinstance(self.basis, HarmonicOscillatorProductBasis):
                iphase = 1
            else:
                iphase = -1

            T = self._input_kinetic[2] if self._input_kinetic is not None and len(self._input_kinetic) > 2 else None
            if T is None:
                T = self.G_terms[2]
            V = self._input_potential[2] if self._input_potential is not None and len(self._input_potential) > 2 else None
            if V is None:
                V = self.V_terms[2]
            self._expansions[2] = (
                    (iphase * 1 / 4) * self.basis.representation('p', 'x', 'x', 'p',
                                                                 coeffs=T,
                                                                 name='T(2)',
                                                                 axes=[[0, 1, 2, 3], [2, 0, 1, 3]],
                                                                 **self.operator_settings
                                                                 )
                    + 1 / 24 * self.basis.representation('x', 'x', 'x', 'x',
                                                         coeffs=V,
                                                         name='V(2)',
                                                         **self.operator_settings
                                                         )
            )
            if self.coriolis_terms is not None or self._input_coriolis is not None:
                Z = self._input_coriolis[0] if self._input_coriolis is not None and len(self._input_coriolis) > 0 else None
                if Z is None:
                    Z = np.sum(self.coriolis_terms[0], axis=0)
                self._expansions[2] += iphase * self.basis.representation('x', 'p', 'x', 'p',
                                                               coeffs=Z,
                                                               name='Coriolis(0)',
                                                               **self.operator_settings
                                                               )
            # else:
            #     self._h2 += 0 * self.basis.representation(coeffs=0,
            #                                               **self.operator_settings
            #                                               )

            if self.pseudopotential_term is not None or self._input_pseudopotential is not None:
                U = self._input_pseudopotential[0] if self._input_pseudopotential is not None and len(self._input_pseudopotential) > 0 else None
                if U is None:
                    U = self.pseudopotential_term[0]
                self._expansions[2] += 1 / 8 * self.basis.representation(coeffs=U,
                                                              name="V'(0)",
                                                              **self.operator_settings
                                                              )

            if self._selection_rules is not None and len(self._selection_rules) > 1:
                self._expansions[2].selection_rules = self._selection_rules[1]

            self._expansions[2].name = "H(2)"

        return self._expansions[2]

    def _get_H(self, o,
               include_potential=True,
               include_gmatrix=True,
               include_coriolis=True,
               include_pseudopotential=True
               ):
        """
        Provides the representation for H(i) in this basis
        """

        if len(self._expansions) < o + 1:
            self._expansions += [None] * (o + 1 - len(self._expansions))

        if self._expansions[o] is None:
            if isinstance(self.basis, HarmonicOscillatorProductBasis):
                iphase = 1
            else:
                iphase = -1

            if include_gmatrix:
                T = self._input_kinetic[o] if self._input_kinetic is not None and len(self._input_kinetic) > o else None
                if T is None:
                    T = self.G_terms[o]
            else:
                T = None

            if include_potential:
                V = self._input_potential[o] if self._input_potential is not None and len(
                    self._input_potential
                ) > o else None
                if V is None:
                    V = self.V_terms[o]
            else:
                V = None

            self.logger.log_print("V({o}) = {V}", o=o, V=V, log_level=self.logger.LogLevel.Never)
            self.logger.log_print("T({o}) = {T}", o=o, T=T, log_level=self.logger.LogLevel.Never)

            if T is not None:
                p_expansion = ['p'] + ['x'] * o + ['p']
                T = (iphase / (2*np.math.factorial(o))) * self.basis.representation(*p_expansion,
                                                                                        coeffs=T,
                                                                                        name='T({})'.format(o),
                                                                                        axes=[
                                                                                            list(range(o+2)),
                                                                                            [o] + list(range(o)) + [o+1]
                                                                                        ],
                                                                                        **self.operator_settings
                                                                                        )
            else:
                T = None

            if V is not None:
                v_expansion = ['x'] * (o + 2)
                V = 1 / np.math.factorial(o + 2) * self.basis.representation(*v_expansion,
                                                                               coeffs=V,
                                                                               axes=[
                                                                                    list(range(o+2)),
                                                                                    [o] + list(range(o)) + [o+1]
                                                                                ],
                                                                               name='V({})'.format(o),
                                                                               **self.operator_settings
                                                                               )
            else:
                V = None

            if T is None and V is None:
                self._expansions[o] = 0
            elif T is None:
                self._expansions[o] = V
            elif V is None:
                self._expansions[o] = T
            else:
                self._expansions[o] = T + V
            if o > 1:
                oz = o - 2
                if include_coriolis and (self.coriolis_terms is not None or self._input_coriolis is not None):
                    Z = self._input_coriolis[oz] if self._input_coriolis is not None and len(self._input_coriolis) > oz else None
                    if Z is None:
                        Z = self.coriolis_terms[oz]
                        for ax in range(2):
                            Z = np.sum(Z, axis=0)

                        # raise Exception(Z)

                    z_exp = ['x', 'p'] + ['x' for _ in range(oz)] + ['x', 'p']
                    Z = iphase * self.basis.representation(*z_exp,
                                                                              coeffs=Z,
                                                                              name='Coriolis({})'.format(oz),
                                                                              **self.operator_settings
                                                                              )
                    if isinstance(self._expansions[o], int):
                        self._expansions[o] = Z
                    else:
                       self._expansions[o] += Z

                    self.logger.log_print("Z({o}) = {Z}", o=oz, Z=Z, log_level=self.logger.LogLevel.Never)

                if include_pseudopotential and (self.pseudopotential_term is not None or self._input_pseudopotential is not None):
                    U = self._input_pseudopotential[oz] if self._input_pseudopotential is not None and len(self._input_pseudopotential) > oz else None
                    if U is None:
                        U = self.pseudopotential_term[oz]
                    u_exp = ['x' for _ in range(oz)]
                    U = 1 / 8 * self.basis.representation(*u_exp,
                                                             coeffs=U,
                                                             name="V'({})".format(oz),
                                                             **self.operator_settings
                                                             )
                    if isinstance(self._expansions[o], int):
                        self._expansions[o] = U
                    else:
                        self._expansions[o] += U

                    self.logger.log_print("U({o}) = {U}", o=oz, U=U, log_level=self.logger.LogLevel.Never)

            if self._selection_rules is not None and len(self._selection_rules) > o-1:
                self._expansions[o].selection_rules = self._selection_rules[o-1]

            self._expansions[o].name = "H({})".format(o)

        return self._expansions[o]

    def get_perturbations(self, expansion_orders):
        """
        Gets the `Representation` objects for the perturbations up through second order

        :param order:
        :type order:
        :return:
        :rtype:
        """
        # we get a benefit from going high first
        if isinstance(expansion_orders, int):
            expansion_orders = self._get_expansion_orders(None, expansion_orders)
        order = max(expansion_orders.values())
        perts = []
        for i in range(order, -1, -1):
            perts.append(
                self._get_H(i,
                            include_potential=expansion_orders['potential'] >= i,
                            include_gmatrix=expansion_orders['gmatrix'] >= i,
                            include_coriolis=expansion_orders['coriolis'] >= i,
                            include_pseudopotential=expansion_orders['pseudopotential'] >= i
                            )
            )
        return tuple(reversed(perts))

    # @property
    # def perturbations(self):
    #     """
    #     Returns `Representation` objects for the different perturbations through second order
    #
    #     :return:
    #     :rtype:
    #     """
    #     return self.get_perturbations(2)

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
                                  Be[a] * (zeta[a, s, t] ** 2) * (w[t] / w[s] + w[s] / w[t])
                          ) for a in range(3))

        return [xst_3, xst_4, xst_cor]

    @classmethod
    def _get_Nielsen_xmat(cls, freqs, v3, v4, zeta, Be):

        # raise Exception(UnitsData.convert("Hartrees", "Wavenumbers") * Be, np.round(zeta, 5))

        ndim = len(freqs)
        if v3 is None or isinstance(v3, (int, np.integer, float, np.floating)) and v3==0:
            v3 = np.zeros((ndim, ndim, ndim))
        if v4 is None or isinstance(v4, (int, np.integer, float, np.floating)) and v4==0:
            v4 = np.zeros((ndim, ndim, ndim, ndim))
        if zeta is None or isinstance(zeta, (int, np.integer, float, np.floating)) and zeta==0:
            zeta = np.zeros((ndim, ndim, ndim))
        if Be is None or isinstance(Be, (int, np.integer, float, np.floating)) and Be==0:
            Be = np.zeros((3,))

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
    def _get_Nielsen_energies_from_x(cls, states, freqs, x_mat, return_split=False):
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

        # h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        # raise Exception(np.sum(x_mat, axis=0) * h2w)
        # raise Exception(x_mat * h2w)

        states = np.array(states) + 1/2 # n+1/2 for harmonic vibrations

        e_harm = np.tensordot(freqs, states, axes=[0, 1])

        outer_states = vec_outer(states, states)

        if x_mat.ndim > 2:
            weights = np.full(x_mat[0].shape, 1/2)
            np.fill_diagonal(weights, 1)
            x_mat = [x * weights for x in x_mat]

            if return_split:
                e_anharm = np.array([np.tensordot(x, outer_states, axes=[[0, 1], [1, 2]]) for x in x_mat])
            else:
                x_mat = np.sum(x_mat, axis=0)
                e_anharm = np.tensordot(x_mat, outer_states, axes=[[0, 1], [1, 2]])
        else:
            weights = np.full(x_mat.shape, 1 / 2)
            np.fill_diagonal(weights, 1)
            x_mat = x_mat * weights
            e_anharm = np.tensordot(x_mat, outer_states, axes=[[0, 1], [1, 2]])

        return e_harm, e_anharm

    def get_Nielsen_xmatrix(self):
        """
        Provides Nielsen's X-Matrix when working in Cartesian coordinates

        :return:
        :rtype:
        """

        # freqs = self.modes.freqs
        # if self.mode_selection is not None:
        #     freqs = freqs[self.mode_selection,]
        freqs = np.diag(self.V_terms[0])
        v3 = self.V_terms[1]
        v4 = self.V_terms[2]

        # raise Exception(np.round( 6 * v3 * h2w))

        if self.coriolis_terms is not None:
            zeta, Be = self.coriolis_terms.get_zetas_and_momi()
        else:
            zeta = Be = None

        x = self._get_Nielsen_xmat(freqs, v3, v4, zeta, Be)

        return x

    def get_Nielsen_energies(self, states, x_mat=None, return_split=False):
        """

        :param states:
        :type states:
        :return:
        :rtype:
        """

        # h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # TODO: figure out WTF the units on this have to be...

        freqs = self.modes.freqs
        if self.mode_selection is not None:
            freqs = freqs[self.mode_selection,]
        if x_mat is None:
            x_mat = self.get_Nielsen_xmatrix()
        else:
            x_mat = np.asanyarray(x_mat)

        harm, anharm = self._get_Nielsen_energies_from_x(states, freqs, x_mat, return_split=return_split)

        # harm = harm / h2w
        # anharm = anharm

        if return_split:
            return harm, anharm, x_mat
        else:
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
        raise NotImplementedError("dead code path")

        nits = order - 1
        if nits >= 0:

            transitions_h1 = self.basis.selection_rules("x", "x", "x")
            transitions_h2 = self.basis.selection_rules("x", "x", "x", "x")

            # the states that can be coupled through H1
            self.logger.log_print('getting states coupled through H^(1)')
            h1_space = states.apply_selection_rules(
                transitions_h1,
                iterations=(order - 1)
            )

            # from second order corrections
            self.logger.log_print('getting states coupled through H^(2)')
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

        raise NotImplementedError("gotta fix this...")

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
                               order=2,
                               deg_extra_order=2
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
        raise NotImplementedError("dead code path")

        # raise Exception(states)
        # need to rewrite this to work better with BasisStateSpace
        if not isinstance(states, BasisStateSpace):
            states = BasisStateSpace(self.basis, states)

        no_states_to_start = coupled_states is None

        coupled_states = self._prep_coupled_states(states, coupled_states, order)

        # raise Exception(coupled_states)

        space_list = [states] + list(coupled_states)
        total_space = BasisMultiStateSpace(np.array(space_list, dtype=object))
        flat_total_space = total_space.to_single().take_unique()

        degeneracies = DegenerateMultiStateSpace.from_spec(self, states, flat_total_space, degeneracies)

        if no_states_to_start and degeneracies is not None:
            for g in degeneracies:
                if len(g) > 1: # these states can now couple to many more states through indirect couplings
                    extra_states = self.get_coupled_space(g, order+2) # we
                    g_inds = tuple(states.find(g))
                    for c, e in zip(coupled_states, extra_states):
                        for i, s in zip(g_inds, e):
                            c[i] = s

        space_list = [states] + list(coupled_states)
        total_space = BasisMultiStateSpace(np.array(space_list, dtype=object))
        flat_total_space = total_space.to_single().take_unique()

        return states, coupled_states, total_space, flat_total_space, degeneracies

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
            coupled_states = self.get_coupled_space(states, order,
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

    def _get_expansion_orders(self, expansion_order, order):
        if expansion_order is None:
            expansion_order = order
        if isinstance(expansion_order, int):
            expansion_order = {
                "kinetic":expansion_order,
                "potential":expansion_order
            }
        if "default" in expansion_order:
            default = expansion_order["default"]
        else:
            default = expansion_order if isinstance(expansion_order, int) else order
        if "kinetic" in expansion_order:
            ke_default = expansion_order["kinetic"]
        else:
            ke_default = default
        for k in ["gmatrix", "pseudopotential", "coriolis"]:
            if k not in expansion_order:
                expansion_order[k] = ke_default
        if "potential" not in expansion_order:
            expansion_order["potential"] = default
        return expansion_order

    def get_solver(self,
                   states,
                   degeneracies=None,
                   allow_post_PT_calc=True,
                   ignore_odd_order_energies=True,
                   use_full_basis=True,
                   order=2,
                   expansion_order=None,
                   memory_constrained=None,
                   target_property_rules=None,
                   **opts
                ):

        expansion_order = self._get_expansion_orders(expansion_order, order)
        h_reps = self.get_perturbations(expansion_order)

        if memory_constrained is None:
            memory_constrained = states.ndim > 20 if memory_constrained is None else memory_constrained

        if use_full_basis and states.full_basis is None:
            states = BasisStateSpace(self.basis, states.indices,
                                     mode=BasisStateSpace.StateSpaceSpec.Indices,
                                     full_basis=CompleteSymmetricGroupSpace(states.ndim, memory_constrained)
                                     )

        # if memory_constrained is True:
        #     memory_constrained = False

        if 'degenerate_states' not in opts:
            opts['degenerate_states'] = degeneracies

        return PerturbationTheorySolver(h_reps, states,
                                        order=order,
                                        logger=self.logger,
                                        checkpointer=self.checkpointer,
                                        parallelizer=self.parallelizer,
                                        allow_post_PT_calc=allow_post_PT_calc,
                                        memory_constrained=memory_constrained,
                                        ignore_odd_order_energies=ignore_odd_order_energies,
                                        target_property_rules=target_property_rules,
                                        **opts
                                        )

    def get_wavefunctions(self,
                          states,
                          degeneracies=None,
                          allow_post_PT_calc=True,
                          ignore_odd_order_energies=True,
                          use_full_basis=True,
                          order=2,
                          expansion_order=None,
                          memory_constrained=None,
                          target_property_rules=None,
                          results=None,
                          degenerate_transformation_layout=None,
                          **opts
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

        # print([self.V_terms[i] for i in range(5)])

        with self.checkpointer:

            with self.logger.block(tag='Computing PT corrections:', printoptions={'linewidth':int(1e8)}):

                if expansion_order is None:
                    expansion_order = order
                self.logger.log_print(
                    [
                        "states: {state_num}",
                        "order: {ord}",
                        "expansion order: {eord}"
                    ],
                    ord=order,
                    state_num=len(states),
                    eord=order if expansion_order is None else expansion_order
                )

                if not isinstance(states, BasisStateSpace):
                    states = BasisStateSpace(self.basis, states)

                if memory_constrained is None:
                    memory_constrained = states.ndim > 20 if memory_constrained is None else memory_constrained

                if results is None:
                    results = self.results
                if results is None:
                    results = NullCheckpointer(None)
                elif isinstance(results, str):
                    results = Checkpointer.from_file(results)

                solver = self.get_solver(
                    states,
                    degeneracies=degeneracies,
                    allow_post_PT_calc=allow_post_PT_calc,
                    ignore_odd_order_energies=ignore_odd_order_energies,
                    use_full_basis=use_full_basis,
                    memory_constrained=memory_constrained,
                    order=order,
                    expansion_order=expansion_order,
                    target_property_rules=target_property_rules,
                    results=results,
                    **opts
                )

                corrs = solver.apply_VPT()

                if memory_constrained:
                    solver.representations = None

        expansion_options = self.expansion_options.copy()
        expansion_options['expansion_order'] = expansion_order

        operator_settings = self.operator_settings.copy()
        operator_settings['memory_constrained'] = memory_constrained

        return PerturbationTheoryWavefunctions(self.molecule, self.basis, corrs,
                                               modes=self.modes,
                                               mode_selection=self.mode_selection,
                                               logger=self.logger,
                                               checkpoint=self.checkpointer,
                                               operator_settings=operator_settings,
                                               expansion_options=expansion_options,
                                               degenerate_transformation_layout=degenerate_transformation_layout
                                               )

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

        #TODO: add a check that makes sure that the number of states is sufficient to fully invert the tensor
        #       i.e. make sure that there are as many states as there are upper-triangle indices
        #       then probably want to check that the states are properly independent, too...
        col = 0
        blocks = []  # to more easily recompose the tensors later
        nterms = 1 + order // 2 # second order should be [2, 1], 4th order should be [3, 2, 1], 6th should be [4, 3, 2, 1]
        for k in range(nterms, 0, -1):
            ninds = []
            # generate the index tensors to loop over
            # we later pull only the upper triangle indices
            inds = np.indices((nmodes,) * k)
            inds = inds.transpose(tuple(range(1, k + 1)) + (0,))
            inds = np.reshape(inds, (-1, k))
            for perm in inds:
                if (np.sort(perm) != perm).any():
                    continue # only want the upper triangle
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
    def get_action_expansion(self,
                             coupled_states=None,
                             degeneracies=None,
                             allow_sakurai_degs=False,
                             allow_post_PT_calc=True,
                             modify_degenerate_perturbations=False,
                             intermediate_normalization=False,
                             ignore_odd_order_energies=True,
                             zero_element_warning=True,
                             state_space_iterations=None,
                             order=2
                             ):
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

        states = BasisStateSpace.from_quanta(HarmonicOscillatorProductBasis(ndim), range(nterms + 1))
        # raise Exception(states)

        wfns = self.get_wavefunctions(states,
                                      coupled_states=coupled_states,
                                      degeneracies=degeneracies,
                                      allow_sakurai_degs=allow_sakurai_degs,
                                      allow_post_PT_calc=allow_post_PT_calc,
                                      modify_degenerate_perturbations=modify_degenerate_perturbations,
                                      intermediate_normalization=intermediate_normalization,
                                      ignore_odd_order_energies=ignore_odd_order_energies,
                                      zero_element_warning=zero_element_warning,
                                      state_space_iterations=state_space_iterations
                                      )

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

            wfns = PerturbationTheoryWavefunctions(self.molecule, self.basis, corrs, modes=self.modes, mode_selection=self.mode_selection, logger=self.logger)
            specs[k] = wfns

        return specs