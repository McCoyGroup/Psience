import numpy as np, itertools, collections, dataclasses
from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer
from ..Spectra import DiscreteSpectrum
from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, HarmonicOscillatorProductBasis, StateMaker
from .Common import PerturbationTheoryException, _safe_dot
__all__ = ['PerturbationTheoryCorrections', 'AnalyticPerturbationTheoryCorrections', 'BasicAPTCorrections']

class PerturbationTheoryCorrections:
    """
    Represents a set of corrections from perturbation theory.
    Can be used to correct other operators in the basis of the original calculation.

    """

    def __init__(self, states, coupled_states, total_basis, energy_corrs, wfn_corrections, all_energy_corrections=None, degenerate_states=None, degenerate_transformation=None, degenerate_energies=None, degenerate_hamiltonians=None, nondeg_hamiltonian_precision=3, logger=None):
        """
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
        ...

    @classmethod
    def from_dicts(cls, states, corrections, **opts):
        """
        :param states: a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
        :type states: dict
        :param corrections: the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT
        :type corrections: dict
        """
        ...

    @property
    def degenerate(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def energies(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def order(self):
        """

        :return:
        :rtype:
        """
        ...

    def take_subspace(self, space):
        """
        Takes only those elements that are in space
        :param space:
        :type space:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def create_coupling_matrix(cls, corrs, states: BasisStateSpace, flat_total_space: BasisStateSpace, nstates, order, filters=None, non_zero_cutoff=1e-14, logger=None):
        """

        :param corrs:
        :type corrs:
        :param states:
        :type states:
        :param flat_total_space:
        :type flat_total_space:
        :param nstates:
        :type nstates:
        :param order:
        :type order:
        :param filters:
        :type filters:
        :param non_zero_cutoff:
        :type non_zero_cutoff:
        :return:
        :rtype:
        """
        ...

    def prune(self, threshold=0.1, in_place=False):
        """
        Returns corrections with couplings less than the given cutoff set to zero

        :param threshold:
        :type threshold:
        :return:
        :rtype:
        """
        ...

    def _take_subham(self, rep, inds):
        """
        Builds a subsampled version of a representation Hamiltonian
        to allow equations to be efficiently solved in subspaces.

        :param rep: representation matrix from which to pull the subspace
        :type rep: SparseArray
        :param inds: indices for the subspace
        :type inds: np.ndarray
        :return:
        :rtype:
        """
        ...

    def get_transformed_Hamiltonians(self, hams, deg_group=None):
        """

        :param corrs:
        :type corrs:
        :param deg_group:
        :type deg_group:
        :return:
        :rtype:
        """
        ...

    def get_degenerate_rotation(self, deg_group, hams, label=None, zero_point_energy=None, local_coupling_hamiltonian=None, local_coupling_order=None):
        """

        :param deg_group:
        :type deg_group:
        :param corrs:
        :type corrs:
        :return:
        :rtype:
        """
        ...

    def get_degenerate_transformation(self, group, hams, gaussian_resonance_handling=False, label=None, zero_point_energy=None, local_coupling_hamiltonian=None, local_coupling_order=None):
        """
        **LLM Docstring**

        Compute the degenerate-perturbation-theory rotation for a single group of resonant/degenerate states: finds the group's positions within this object's overall state space (warning about any states in the group that aren't present), and, if the group has more than one matched state (and isn't skipped due to Gaussian-style high-order resonance handling), diagonalizes the corresponding non-degenerate Hamiltonian block via `get_degenerate_rotation`.

        :param group: the group of mutually resonant/degenerate states to build the transformation for
        :type group: BasisStateSpace
        :param hams: the Hamiltonian correction matrices to build the block from
        :type hams: list[np.ndarray]
        :param gaussian_resonance_handling: whether to skip building a rotation for groups whose states have more than 2 quanta of excitation (mimicking Gaussian's resonance-handling convention)
        :type gaussian_resonance_handling: bool
        :param label: a label used for logging
        :type label: str | None
        :param zero_point_energy: the zero-point energy, used when building/logging the non-degenerate Hamiltonian block
        :type zero_point_energy: float | None
        :param local_coupling_hamiltonian: an explicit local coupling Hamiltonian to use instead of building one from `hams`
        :type local_coupling_hamiltonian: np.ndarray | None
        :param local_coupling_order: the perturbative order to build the local coupling Hamiltonian to, if not given explicitly
        :type local_coupling_order: int | None
        :return: `(deg_inds, H_nd, deg_rot, deg_engs)` -- the group's indices within this state space, the non-degenerate Hamiltonian block, the diagonalizing rotation, and the resulting (sorted) degenerate energies; `H_nd`/`deg_rot`/`deg_engs` are all `None` if the group doesn't need (or isn't eligible for) degenerate treatment
        :rtype: tuple
        """
        ...

    @staticmethod
    def default_state_filter(state, couplings, energy_cutoff=None, energies=None, basis=None, target_modes=None):
        """
        Excludes modes that differ in only one position, prioritizing states with fewer numbers of quanta
        (potentially add restrictions to high frequency modes...?)

        :param input_state:
        :type input_state:
        :param couplings:
        :type couplings:
        :return:
        :rtype:
        """
        ...

    def find_strong_couplings(self, threshold=0.1, state_filter=None):
        """
        Finds positions in the expansion matrices where the couplings are too large

        :param threshold:
        :type threshold:
        :return:
        :rtype:
        """
        ...

    def format_strong_couplings_report(self, couplings=None, threshold=0.1, int_fmt='{:>3.0f}', padding='{:<8}', join=True, use_excitations=True):
        """
        **LLM Docstring**

        Format a human-readable report of the states found by `find_strong_couplings` (or an explicitly supplied `couplings` dict), listing each state alongside the other states it's strongly coupled to at each perturbative order.

        :param couplings: the strong-coupling data to format; computed via `find_strong_couplings` if not given
        :type couplings: dict | None
        :param threshold: the coupling-strength threshold forwarded to `find_strong_couplings` if `couplings` isn't given
        :type threshold: float
        :param int_fmt: the format string used for each quantum-number column
        :type int_fmt: str
        :param padding: the format string used for row labels/indentation
        :type padding: str
        :param join: whether to join the report lines into a single string, or return them as a list
        :type join: bool
        :param use_excitations: whether to display states as their excitation-quantum-number vectors (rather than raw basis indices)
        :type use_excitations: bool
        :return: the formatted report, as a single string or a list of lines
        :rtype: str | list[str]
        """
        ...

    def collapse_strong_couplings(self, sc: dict):
        """

        :param sc:
        :type sc:
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def _fmt_operator_rep(full_ops, operator_symbol, conversion, real_fmt='{:>.8e}', padding_fmt='{:>16}'):
        """
        **LLM Docstring**

        Format a set of `<a|A(c)|b>`-labeled operator sub-representations into aligned text columns, for use in logging an operator's matrix representation.

        :param full_ops: an iterable of `((a, b, c), subrep)` pairs, each giving a bra/ket/order label and the corresponding operator sub-matrix (dense, sparse, or a scalar zero placeholder)
        :type full_ops: Iterable[tuple]
        :param operator_symbol: the symbol used to label the operator (e.g. `"A"`) in each column header
        :type operator_symbol: str
        :param conversion: a unit-conversion factor to multiply each sub-representation by before formatting
        :type conversion: float | None
        :param real_fmt: the format string used for each numeric matrix element
        :type real_fmt: str
        :param padding_fmt: the format string used to pad each formatted element to a fixed width
        :type padding_fmt: str
        :return: the formatted lines, with the header/tag line first
        :rtype: list[str]
        :raises ValueError: if a sub-representation is a nonzero scalar (not a recognized representation), or if a zero scalar appears before any array-valued sub-representation has established the operator's dimension
        """
        ...

    def operator_representation(self, operator_expansion, order=None, subspace=None, contract=True, logger_symbol='A', logger_conversion=None):
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
        ...

    def get_overlap_matrices(self):
        """
        Returns the overlap matrices for the set of corrections
        at each order of correction

        :return:
        :rtype:
        """
        ...

    def savez(self, file):
        """
        **LLM Docstring**

        Intended to serialize the corrections to an `npz` file, but currently disabled -- immediately raises, noting the implementation is outdated; use `to_state`/`from_state` instead.

        :param file: the target file path
        :type file: str
        :return: never returns
        :rtype: None
        :raises NotImplementedError: always, with a message noting this method is outdated
        """
        ...

    @classmethod
    def loadz(cls, file):
        """
        **LLM Docstring**

        Intended to reconstruct corrections from an `npz` file previously written by `savez`, but currently disabled -- immediately raises, noting the implementation is outdated; use `to_state`/`from_state` instead.

        :param file: the source file path
        :type file: str
        :return: never returns
        :rtype: PerturbationTheoryCorrections
        :raises NotImplementedError: always, with a message noting this method is outdated
        """
        ...

    def to_state(self, serializer=None):
        """
        **LLM Docstring**

        Serialize this object's core data (states, coupled states, total basis, energy/wavefunction corrections, and any degenerate-perturbation-theory data) into a plain dict.

        :param serializer: accepted for interface consistency but not used in this method's body
        :type serializer: object | None
        :return: the serialized state dict
        :rtype: dict
        """
        ...

    @classmethod
    def from_state(cls, data, serializer=None):
        """
        **LLM Docstring**

        Reconstruct a `PerturbationTheoryCorrections` object from a previously serialized state dict, deserializing the state-space objects via the given `serializer` and delegating to `from_dicts`.

        :param data: the serialized state, as produced by `to_state`
        :type data: dict
        :param serializer: the serializer used to deserialize the state-space objects
        :type serializer: object
        :return: the reconstructed corrections object
        :rtype: PerturbationTheoryCorrections
        """
        ...
BasicAPTCorrections = collections.namedtuple('PTCorrections', ['initial_states', 'final_states', 'corrections'])

@dataclasses.dataclass
class AnalyticPerturbationTheoryCorrections:
    states: BasisStateSpace
    state_lists: 'list[tuple[np.ndarray, np.ndarray]]'
    _energies: np.ndarray = None
    _transition_moments: 'Iterable[np.ndarray]' = None
    _spectra: 'Iterable[DiscreteSpectrum]' = None
    _deperturbed_energies: np.ndarray = None
    _deperturbed_transition_moments: 'Iterable[np.ndarray]' = None
    _deperturbed_spectra: DiscreteSpectrum = None
    degenerate_states: 'Iterable[BasisStateSpace]' = None
    only_degenerate_terms: 'bool' = True
    _degenerate_hamiltonians: 'Iterable[np.ndarray]' = None
    _degenerate_coefficients: 'Iterable[np.ndarray]' = None
    _degenerate_state_list_transformations: 'Iterable[list[np.ndarray, np.ndarray]]' = None
    energy_corrections: BasicAPTCorrections = None
    transition_moment_corrections: 'Iterable[BasicAPTCorrections]' = None
    degenerate_hamiltonian_corrections: 'Iterable[BasicAPTCorrections]' = None
    operator_corrections: 'Iterable[BasicAPTCorrections]' = None
    _deperturbed_operator_values: 'Iterable[np.ndarray]' = None
    _operator_values: 'Iterable[np.ndarray]' = None
    operator_keys: 'Iterable[Any]' = None
    logger: 'Logger' = None
    _zpe_pos: int = None

    def get_zpe_pos(self) -> int:
        """
        **LLM Docstring**

        Find (and cache) the index of the zero-point-energy (ground) state within `self.states`, falling back to index `0` if the all-zeros excitation vector can't be found.

        :return: the ZPE state's index
        :rtype: int
        """
        ...

    @property
    def energies(self) -> np.ndarray:
        """
        **LLM Docstring**

        The final state energies: the sum of the (per-order) energy corrections if there are no degenerate states, otherwise the degenerate-perturbation-theory-corrected energies (computed via `get_degenerate_transformations`, which also populates the cached degenerate Hamiltonians/coefficients).

        :return: the state energies
        :rtype: np.ndarray
        """
        ...

    @property
    def deperturbed_energies(self) -> np.ndarray:
        """
        **LLM Docstring**

        The (cached) deperturbed state energies -- the sum of the per-order energy corrections, without any degenerate-perturbation-theory rotation applied.

        :return: the deperturbed energies
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def handle_degenerate_transformation(cls, degenerate_ham):
        """
        **LLM Docstring**

        Diagonalize a degenerate-block Hamiltonian and reorder its eigenvalues/eigenvectors so that each output state is matched (by maximum overlap) to a distinct input state, avoiding two output states mapping to the same input state.

        :param degenerate_ham: the (Hermitian) degenerate-block Hamiltonian matrix to diagonalize
        :type degenerate_ham: np.ndarray
        :return: `(deg_engs, deg_transf)` -- the reordered eigenvalues and eigenvector (transformation) matrix
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def get_degenerate_transformations(self, basis, energies):
        """
        **LLM Docstring**

        Apply degenerate perturbation theory block by block: for each group of degenerate states, build its Hamiltonian block (summing the per-order corrections, with the diagonal set to the current deperturbed energies), diagonalize it via `handle_degenerate_transformation`, and write the resulting energies back into the running energy array.

        :param basis: the state space the energies are indexed against
        :type basis: BasisStateSpace
        :param energies: the deperturbed energies to correct, indexed against `basis`
        :type energies: np.ndarray
        :return: `(energies, (hams, transf))` -- the corrected energies, and the list of per-block Hamiltonians and diagonalizing transformations
        :rtype: tuple[np.ndarray, tuple[list[np.ndarray], list[np.ndarray]]]
        """
        ...

    @property
    def degenerate_hamiltonians(self):
        """
        **LLM Docstring**

        The (cached) per-block degenerate Hamiltonians, computed as a side effect of evaluating `energies` if not already cached.

        :return: the list of degenerate-block Hamiltonians
        :rtype: list[np.ndarray]
        """
        ...

    @property
    def degenerate_coefficients(self):
        """
        **LLM Docstring**

        The (cached) per-block degenerate-perturbation-theory mixing coefficients, computed as a side effect of evaluating `energies` if not already cached.

        :return: the list of degenerate-block mixing-coefficient matrices
        :rtype: list[np.ndarray]
        """
        ...

    def get_freqs(self):
        """
        **LLM Docstring**

        Compute the vibrational transition frequencies (final energies) relative to the zero-point energy.

        :return: the frequencies
        :rtype: np.ndarray
        """
        ...

    def get_deperturbed_freqs(self):
        """
        **LLM Docstring**

        Compute the deperturbed vibrational transition frequencies relative to the deperturbed zero-point energy, if degenerate states are present; otherwise falls back to `get_freqs`.

        :return: the deperturbed frequencies
        :rtype: np.ndarray
        """
        ...

    @property
    def degenerate_transformation_pairs(self):
        """
        **LLM Docstring**

        The (cached) per-block `(row_tf, col_tf)` transformation pairs mapping each `state_lists` block's raw initial/final states onto their degenerate-perturbation-theory-mixed counterparts, computed lazily via `_get_degenerate_tfs_mats`.

        :return: the list of per-block `[row_tf, col_tf]` transformation pairs
        :rtype: list[list[np.ndarray]]
        """
        ...

    def _get_degenerate_tfs_mats(self, logger=None):
        """
        **LLM Docstring**

        Build, for each `(initial_states, final_states)` block in `state_lists`, the pair of transformation matrices that map the raw (unmixed) initial/final states onto their degenerate-perturbation-theory-mixed linear combinations -- identity everywhere except in rows/columns corresponding to states that participate in a degenerate block, which are replaced by the corresponding column of that block's mixing coefficients -- validating that each resulting transformation is (approximately) unitary.

        :param logger: logger to report the degenerate blocks/transformations to; defaults to `self.logger`
        :type logger: Logger | str | None
        :return: the list of per-block `[row_tf, col_tf]` transformation-matrix pairs
        :rtype: list[list[np.ndarray]]
        :raises ValueError: if a computed row or column transformation isn't (approximately) unitary
        """
        ...

    def _apply_degs_to_corrs(self, corrs, logger=None):
        """
        **LLM Docstring**

        Apply the per-block degenerate-perturbation-theory transformation pairs (from `degenerate_transformation_pairs`) to a set of per-block correction matrices (or, for multi-axis quantities like transition moments, a list of such per-block matrices per axis), sandwiching each block between the row and column transformations.

        :param corrs: the per-block correction data to transform, either a flat list of per-block matrices or a list of such lists (one per axis/operator component)
        :type corrs: list
        :param logger: logger, accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | str | None
        :return: the transformed correction data, in the same nested structure as `corrs`
        :rtype: list
        """
        ...

    @property
    def transition_moments(self):
        """
        **LLM Docstring**

        The (cached) final transition-dipole moments: the deperturbed transition moments if there are no degenerate states, otherwise those moments rotated by the degenerate-perturbation-theory transformation (via `_apply_degs_to_corrs`).

        :return: the per-axis, per-block transition moments
        :rtype: list
        """
        ...

    @property
    def harmonic_transition_moments(self):
        """
        **LLM Docstring**

        The purely harmonic (zeroth-order) transition-dipole moments, extracted from the first entry of each block's transition-moment corrections.

        :return: the per-axis, per-block harmonic transition moments
        :rtype: list
        """
        ...

    @property
    def deperturbed_transition_moments(self):
        """
        **LLM Docstring**

        The (cached) deperturbed transition-dipole moments, summing each block's transition-moment corrections over all perturbative orders.

        :return: the per-axis, per-block deperturbed transition moments
        :rtype: list
        """
        ...

    def get_spectra(self, energies, transition_moments):
        """
        **LLM Docstring**

        Build a per-block list of discrete IR spectra from a set of state energies and transition moments, one spectrum per initial state in each `state_lists` block, with transition frequencies computed relative to that initial state's energy.

        :param energies: the state energies to compute transition frequencies from
        :type energies: np.ndarray
        :param transition_moments: the per-axis, per-block transition-moment data (as returned by `transition_moments`/`deperturbed_transition_moments`/`harmonic_transition_moments`)
        :type transition_moments: list
        :return: the per-block lists of per-initial-state discrete spectra
        :rtype: list[list[DiscreteSpectrum]]
        """
        ...

    @property
    def harmonic_spectra(self):
        """
        **LLM Docstring**

        The purely harmonic (zeroth-order) IR spectra, built from the zeroth-order energy corrections and harmonic transition moments via `get_spectra`.

        :return: the per-block lists of per-initial-state harmonic spectra
        :rtype: list[list[DiscreteSpectrum]]
        """
        ...

    @property
    def deperturbed_spectra(self):
        """
        **LLM Docstring**

        The (cached) deperturbed IR spectra, built from the deperturbed energies and transition moments via `get_spectra`.

        :return: the per-block lists of per-initial-state deperturbed spectra
        :rtype: list[list[DiscreteSpectrum]]
        """
        ...

    @property
    def spectra(self):
        """
        **LLM Docstring**

        The (cached) final IR spectra: the deperturbed spectra if there are no degenerate states, otherwise the spectra built from the degenerate-perturbation-theory-corrected energies and transition moments.

        :return: the per-block lists of per-initial-state spectra
        :rtype: list[list[DiscreteSpectrum]]
        """
        ...

    @property
    def deperturbed_operator_values(self):
        """
        **LLM Docstring**

        The (cached) deperturbed values of any extra tracked operators, summing each operator's per-block corrections over all perturbative orders.

        :return: the per-operator, per-block deperturbed operator values
        :rtype: list
        """
        ...

    @property
    def operator_values(self):
        """
        **LLM Docstring**

        The (cached) final values of any extra tracked operators: the deperturbed operator values if there are no degenerate states, otherwise those values rotated by the degenerate-perturbation-theory transformation (via `_apply_degs_to_corrs`).

        :return: the per-operator, per-block operator values
        :rtype: list
        """
        ...