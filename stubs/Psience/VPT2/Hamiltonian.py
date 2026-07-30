"""
Provides support for build perturbation theory Hamiltonians
"""
import numpy as np, itertools, time, math
import McUtils.Numputils as nput
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer, ParameterManager
from McUtils.Parallelizers import Parallelizer
from McUtils.Combinatorics import CompleteSymmetricGroupSpace
from ..Modes import MixtureModes
from ..Molecools import Molecule
from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, BraKetSpace, HarmonicOscillatorProductBasis
from .Common import PerturbationTheoryException
from .Terms import PotentialTerms, KineticTerms, CoriolisTerm, PotentialLikeTerm, DipoleTerms, OperatorTerms
from .Solver import PerturbationTheorySolver, PerturbationTheoryCorrections
from .Wavefunctions import PerturbationTheoryWavefunctions
__all__ = ['PerturbationTheoryHamiltonian', 'PerturbationTheoryCorrections']
__reload_hook__ = ['..BasisReps', '..Molecools', '.Terms', '.Solver', '.Wavefunctions']

class PerturbationTheoryHamiltonian:
    """
    Represents the main Hamiltonian used in the perturbation theory calculation.
    Uses a harmonic oscillator basis for representing H0, H1, H2 etc.
    The PT process is split into a PerturbationTheorySolver and a PerturbationTheoryHamiltonian
    where the Hamiltonian just implements the split of the perturbation and the Solver manages the equations.
    """

    def __init__(self, molecule=None, n_quanta=None, modes=None, mode_selection=None, mode_transformation=None, rephase_modes=None, local_mode_couplings=False, local_mode_coupling_order=None, full_surface_mode_selection=None, potential_derivatives=None, include_potential=True, include_gmatrix=True, include_coriolis_coupling=True, include_pseudopotential=True, include_only_mode_couplings=None, potential_terms=None, allow_higher_potential_terms=False, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, include_dipole=True, dipole_terms=None, selection_rules=None, operator_chunk_size=None, operator_coefficient_threshold=None, matrix_element_threshold=None, logger=None, checkpoint=None, results=None, parallelizer=None, **expansion_options):
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
        ...

    @classmethod
    def from_fchk(cls, file, internals=None, mode_selection=None, **kw):
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
        ...

    class TermGetter:

        def __init__(self, base_terms, input_terms, mode_selection=None):
            """
            **LLM Docstring**

            Wrap a pair of term sources (a set of "base" expansion terms and/or a set of explicitly supplied "input" terms) so that indexing into this object transparently prefers the input terms where available, falling back to the base terms otherwise.

            :param base_terms: the fallback expansion terms (e.g. computed potential/kinetic terms), indexable by order
            :type base_terms: object | None
            :param input_terms: explicitly supplied terms to prefer, indexable by order
            :type input_terms: object | None
            :param mode_selection: a subset of mode indices to restrict `input_terms` to when they're used
            :type mode_selection: Iterable[int] | None
            :return: None
            :rtype: None
            """
            ...

        def __getitem__(self, o):
            """
            **LLM Docstring**

            Fetch the term at order `o`: prefers the (mode-subset-restricted) input term if one is available at that order, otherwise falls back to the (adjusted) base term.

            :param o: the perturbative order to fetch
            :type o: int
            :return: the term at that order, or `None` if neither source has one
            :rtype: object | None
            """
            ...

        def take_mode_subset(self, V, sel):
            """
            **LLM Docstring**

            Restrict a term tensor to a subset of modes along every axis, if a mode selection is configured.

            :param V: the term tensor to restrict
            :type V: np.ndarray | object
            :param sel: the mode-index subset to restrict to; if `None`, `V` is returned unchanged
            :type sel: Iterable[int] | None
            :return: the restricted tensor (or `V` unchanged if `sel`/`V` isn't a plain array)
            :rtype: np.ndarray | object
            """
            ...

        def adjust_base_term(self, V):
            """
            **LLM Docstring**

            Hook for post-processing a base term before it's returned; on the base `TermGetter` this is a no-op, returning `V` unchanged. Subclasses (e.g. `CoriolisTermGetter`) override this to apply term-specific adjustments.

            :param V: the base term to adjust
            :type V: object
            :return: `V`, unchanged
            :rtype: object
            """
            ...

    class CoriolisTermGetter(TermGetter):

        def adjust_base_term(self, Z):
            """
            **LLM Docstring**

            Collapse the Coriolis term's 3 rotational-axis components down to a single tensor by summing the `xx`, `yy`, and `zz` diagonal blocks.

            :param Z: the Coriolis term, with its leading two axes indexing the rotational axes
            :type Z: np.ndarray
            :return: the summed (axis-collapsed) tensor
            :rtype: np.ndarray
            """
            ...

    @property
    def dipole_terms(self):
        """
        **LLM Docstring**

        The (lazily constructed and cached) `DipoleTerms` object used to expand the dipole surface, or `None` if dipole terms weren't requested (`include_dipole=False` at construction).

        :return: the dipole-terms expansion object, or `None`
        :rtype: DipoleTerms | None
        """
        ...

    @classmethod
    def prep_local_couplings(cls, local_mode_couplings):
        """
        **LLM Docstring**

        Normalize the `local_mode_couplings` constructor argument into a `[v0, g0]` pair of local-mode potential/kinetic coupling matrices: `False`/falsy stays `False` (no local coupling), `True` becomes `[None, None]` (defer to defaults), a bare coupling matrix is split evenly between potential and kinetic contributions, and an explicit `(v0, g0)` pair is validated and passed through.

        :param local_mode_couplings: the raw local-mode-coupling specification
        :type local_mode_couplings: bool | np.ndarray | tuple
        :return: `False`, or the resolved `[v0, g0]` coupling-matrix pair
        :rtype: bool | list[np.ndarray | None]
        :raises ValueError: if a 2-element specification doesn't resolve to a valid 2D coupling matrix
        """
        ...

    def _get_H(self, o, include_potential=True, include_gmatrix=True, include_coriolis=True, include_pseudopotential=True, include_modes=None, local_mode_couplings=None, return_reps=True):
        """
        Provides the representation for H(i) in this basis
        """
        ...

    def prep_operator_terms(self, coeffs, order):
        """
        **LLM Docstring**

        Build the perturbative expansion terms for an arbitrary operator given as a `[constant, deriv1, deriv2, ...]` list of coefficients, padding any missing lower-order derivative tensors with zeros (inferred from the dimensionality of the first non-numeric entry) before expanding the whole thing through `OperatorTerms`.

        :param coeffs: the operator's raw coefficient list: a constant term followed by successive derivative tensors (which may start at a higher order, with lower orders implicitly zero)
        :type coeffs: list
        :param order: the highest derivative order to expand to
        :type order: int
        :return: `[const] + expansion_terms`, the constant term followed by the expanded (mode-basis) derivative terms
        :rtype: list
        :raises ValueError: if every entry in `coeffs[1:]` is purely numeric (so the intended tensor dimensionality can't be inferred)
        """
        ...

    def get_perturbations(self, expansion_orders, return_reps=True, order=None):
        """
        Gets the `Representation` objects for the perturbations up through second order

        :param order:
        :type order:
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def _Nielsen_xss(s, w, v3, v4, zeta, Be, ndim, return_components=False):
        """
        **LLM Docstring**

        Compute the same-mode (`X_ss`) anharmonicity constant contribution from Harald Nielsen's second-order perturbation-theory formulas, combining the quartic force-constant term with the (commented-out, currently unused) cubic and Coriolis contributions.

        :param s: the mode index (negative indices wrap via `ndim`)
        :type s: int
        :param w: the harmonic frequencies
        :type w: np.ndarray
        :param v3: the cubic force-constant tensor
        :type v3: np.ndarray
        :param v4: the quartic force-constant tensor
        :type v4: np.ndarray
        :param zeta: the Coriolis zeta-constant tensor
        :type zeta: np.ndarray
        :param Be: the rotational constants
        :type Be: np.ndarray
        :param ndim: the number of vibrational modes
        :type ndim: int
        :param return_components: whether to return the individual cubic/quartic/Coriolis contributions instead of their (partial) sum
        :type return_components: bool
        :return: the `X_ss` contribution, or its `[cubic, quartic, coriolis]` components if `return_components` is set
        :rtype: float | list
        """
        ...

    @staticmethod
    def _Nielsen_xst(s, t, w, v3, v4, zeta, Be, ndim, return_components=False):
        """
        **LLM Docstring**

        Compute the cross-mode (`X_st`) anharmonicity constant contribution from Harald Nielsen's second-order perturbation-theory formulas, combining cubic, quartic, and Coriolis force-constant terms.

        :param s: the first mode index (negative indices wrap via `ndim`)
        :type s: int
        :param t: the second mode index (negative indices wrap via `ndim`)
        :type t: int
        :param w: the harmonic frequencies
        :type w: np.ndarray
        :param v3: the cubic force-constant tensor
        :type v3: np.ndarray
        :param v4: the quartic force-constant tensor
        :type v4: np.ndarray
        :param zeta: the Coriolis zeta-constant tensor
        :type zeta: np.ndarray
        :param Be: the rotational constants
        :type Be: np.ndarray
        :param ndim: the number of vibrational modes
        :type ndim: int
        :param return_components: whether to return the individual cubic/quartic/Coriolis contributions instead of their sum
        :type return_components: bool
        :return: `[cubic, quartic, coriolis]` -- the `X_st` contribution split into its component terms (summed if `return_components` is `False`, listed per-term if `True`)
        :rtype: list
        """
        ...

    @classmethod
    def _get_Nielsen_xmat(cls, freqs, v3, v4, zeta, Be):
        """
        **LLM Docstring**

        Build the full symmetric `X` anharmonicity-constant matrix (used in Nielsen's second-order vibrational energy formula) from the harmonic frequencies and cubic/quartic force constants, Coriolis zeta constants, and rotational constants, computing each same-mode entry via `_Nielsen_xss` and each cross-mode entry via `_Nielsen_xst`.

        :param freqs: the harmonic frequencies
        :type freqs: np.ndarray
        :param v3: the cubic force-constant tensor; treated as all-zero if `None`/`0`
        :type v3: np.ndarray | None
        :param v4: the quartic force-constant tensor; treated as all-zero if `None`/`0`
        :type v4: np.ndarray | None
        :param zeta: the Coriolis zeta-constant tensor; treated as all-zero if `None`/`0`
        :type zeta: np.ndarray | None
        :param Be: the rotational constants; treated as all-zero if `None`/`0`
        :type Be: np.ndarray | None
        :return: the `(3, ndim, ndim)` X-matrix (one symmetric matrix per rotational axis... actually the leading axis holds the cubic/quartic/Coriolis split before being collapsed into a single symmetric `ndim x ndim` matrix per mode pair)
        :rtype: np.ndarray
        """
        ...

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
        ...

    def get_Nielsen_xmatrix(self, freqs=None, v3=None, v4=None, zeta_Be=None):
        """
        Provides Nielsen's X-Matrix when working in Cartesian coordinates

        :return:
        :rtype:
        """
        ...

    def get_Nielsen_energies(self, states, x_mat=None, freqs=None, v3=None, v4=None, zeta_Be=None, return_split=False, return_X=False):
        """

        :param states:
        :type states:
        :return:
        :rtype:
        """
        ...

    def get_2nd_order_freqs(self, states, *, freqs=None, V_terms=None, G_terms=None):
        """

        :param states:
        :type states:
        :return:
        :rtype:
        """
        ...

    def _get_expansion_orders(self, expansion_order, order):
        """
        **LLM Docstring**

        Normalize the `expansion_order` specification into a full per-term-type dict (`'potential'`, `'gmatrix'`, `'pseudopotential'`, `'coriolis'`), filling in unspecified entries from a `'default'` key (or the overall PT `order` if no default/kinetic entry is given) and from a shared `'kinetic'` fallback for the three kinetic-energy-related term types.

        :param expansion_order: the raw expansion-order specification: `None` (use `order` for everything), a bare integer (used for both kinetic and potential), or a partial dict of per-term-type orders
        :type expansion_order: int | dict | None
        :param order: the overall perturbation-theory order, used as the ultimate fallback
        :type order: int
        :return: the fully populated per-term-type expansion-order dict
        :rtype: dict
        """
        ...

    def get_solver(self, states, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, **opts):
        """
        **LLM Docstring**

        Build a `PerturbationTheorySolver` for the given target states: resolves the per-term expansion orders, builds the corresponding Hamiltonian perturbation representations (via `get_perturbations`), coerces `states` into a `BasisStateSpace` (optionally attaching a complete symmetric-group full basis), and constructs the solver with this Hamiltonian's logger/checkpointer/parallelizer and local-mode-coupling settings.

        :param states: the target states to solve for
        :type states: BasisStateSpace | Iterable[int] | Iterable[Iterable[int]]
        :param degeneracies: the degenerate-state specification, forwarded as `opts['degenerate_states']` if not already present in `opts`
        :type degeneracies: object | None
        :param allow_post_PT_calc: whether to allow post-perturbation-theory (e.g. degenerate) energy corrections
        :type allow_post_PT_calc: bool
        :param ignore_odd_order_energies: whether to skip odd-order energy corrections (which should vanish by symmetry)
        :type ignore_odd_order_energies: bool
        :param use_full_basis: whether to attach a complete symmetric-group full basis to `states` if one isn't already present
        :type use_full_basis: bool
        :param order: the perturbation-theory order to solve to
        :type order: int
        :param expansion_order: the per-term expansion orders; resolved via `_get_expansion_orders` if not already a full dict
        :type expansion_order: int | dict | None
        :param memory_constrained: whether to use a memory-constrained solving strategy; defaults to `True` if the state space has more than 20 dimensions
        :type memory_constrained: bool | None
        :param target_property_rules: property-specific selection rules to restrict the solve to
        :type target_property_rules: object | None
        :param opts: extra options forwarded to the `PerturbationTheorySolver` constructor
        :type opts: dict
        :return: the constructed solver
        :rtype: PerturbationTheorySolver
        """
        ...

    def get_wavefunctions(self, states, initial_states=None, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, results=None, degenerate_transformation_layout=None, return_solver=False, **opts):
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
        ...

    @classmethod
    def _invert_action_expansion_tensors(cls, states, energies, order):
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
        ...

    def get_action_expansion(self, coupled_states=None, degeneracies=None, allow_sakurai_degs=False, allow_post_PT_calc=True, modify_degenerate_perturbations=False, intermediate_normalization=False, ignore_odd_order_energies=True, zero_element_warning=True, state_space_iterations=None, order=2):
        """
        Gets the expansion of the energies in terms of Miller's "classical actions" by
        doing just enough PT to invert the matrix

        :param order:
        :type order:
        :return:
        :rtype:
        """
        ...

    def get_breakdown(self, states, coupled_states=None, degeneracies=None, order=2):
        """
        **LLM Docstring**

        Intended to compute a term-by-term breakdown of the VPT energies (harmonic-only, +cubic, +quartic, full) for a set of states, but currently disabled -- immediately raises, noting the surrounding solver machinery has changed and this method hasn't been updated to match, so the remaining implementation below is unreachable legacy code.

        :param states: the target states
        :type states: object
        :param coupled_states: the coupled-state space to use
        :type coupled_states: object | None
        :param degeneracies: the degenerate-state specification
        :type degeneracies: object | None
        :param order: the perturbation-theory order
        :type order: int
        :return: never returns
        :rtype: dict
        :raises NotImplementedError: always, noting the solver's form has changed and this method needs corresponding updates
        """
        ...