"""
Provides classes to support wave functions coming from VPT calculations
"""
import numpy as np, itertools as ip, time, enum, math
from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Scaffolding import ParameterManager, Checkpointer, NullCheckpointer
from McUtils.Data import UnitsData
from ..Wavefun import Wavefunctions
from ..BasisReps import SimpleProductBasis, BraKetSpace
from .Terms import DipoleTerms, OperatorTerms
from .Solver import PerturbationTheoryCorrections
__all__ = ['PerturbationTheoryWavefunctions']
__reload_hook__ = ['..Wavefun']

class PerturbationTheoryWavefunctions(Wavefunctions):
    """
    These things are fed the first and second order corrections
    """

    def __init__(self, mol, basis, corrections, initial_states=None, modes=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, logger=None, checkpoint=None, results=None, operator_settings=None, expansion_options=None, degenerate_transformation_layout=None):
        """
        :param mol: the molecule the wavefunction is for
        :type mol: Molecule
        :param basis: the basis the expansion is being done in
        :type basis: SimpleProductBasis
        :param corrections: the corrections to the terms
        :type corrections: PerturbationTheoryCorrections
        """
        ...

    def get_dimension(self):
        ...

    @property
    def energies(self):
        ...

    @energies.setter
    def energies(self, engs):
        ...

    def to_state(self, serializer=None):
        ...

    @classmethod
    def from_state(cls, data, serializer=None):
        ...

    @property
    def degenerate_transformation(self):
        ...

    @property
    def initial_state_indices(self):
        ...

    def energies_to_order(self, order=None):
        ...

    @property
    def deperturbed_energies(self):
        ...

    def deperturbed_frequencies(self, order=None):
        ...

    @property
    def order(self):
        ...

    def expectation(self, operator, other=None):
        ...

    @property
    def zero_order_energies(self):
        ...

    def get_M0(self, mu_0):
        ...

    def get_M1(self, mu_1):
        ...

    def get_M2(self, mu_2):
        ...

    def get_M3(self, mu_3):
        ...

    def get_Mi(self, i, mu, base_sym='M'):
        ...

    def _get_rep_states(self, m, k, order, ket_spaces, bra_spaces):
        ...

    @staticmethod
    def _is_zero(x):
        ...

    def _mu_representations(self, mu, M, space, bra_spaces, ket_spaces, order, partitioning, rep_inds, allow_higher_dipole_terms=False):
        ...

    def _transition_moments(self, mu_x, mu_y, mu_z, correction_terms=None, lower_states=None, excited_states=None, degenerate_transformation=None, partitioning=None, order=None):
        """
        Calculates the x, y, and z components of the transition moment between
        the wavefunctions stored.
        By default, only calculates moments involving the ground states.

        :param mu_x: dipole x components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_x: Iterable[np.ndarray]
        :param mu_y: dipole y components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_y: Iterable[np.ndarray]
        :param mu_z: dipole z components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_z: Iterable[np.ndarray]
        :param partitioning: whether to put linear and constant in mu_0 (standard) or linear in mu_1 (intuitive)
        :type partitioning: None | str
        :return:
        :rtype:
        """
        ...

    @classmethod
    def _compute_tmom_to_order(cls, transition_moment_components, order):
        ...

    def _operator_representation(self, ders, M, space, bra_spaces, ket_spaces, order, rep_inds, base_sym='O'):
        ...

    def _operator_corrections(self, operator, correction_terms=None, lower_states=None, excited_states=None, degenerate_transformation=None, order=None):
        """
        Calculates operator corrections between the wavefunctions stored.

        :param operator: operator components (1st, 2nd, 3rd derivatives in normal modes)
        :type operator: Iterable[np.ndarray]
        :return:
        :rtype:
        """
        ...

    @classmethod
    def _compute_corr_to_order(cls, components, order):
        ...

    class TermHolder(tuple):
        """symbolic wrapper on Tuple so we can know that we've canonicalized some term"""

    @property
    def dipole_terms(self):
        ...

    @dipole_terms.setter
    def dipole_terms(self, v):
        ...

    class DipolePartitioningMethod(enum.Enum):
        Standard = 'standard'
        Intuitive = 'intuitive'

    @property
    def dipole_partitioning(self):
        ...

    @dipole_partitioning.setter
    def dipole_partitioning(self, p):
        ...

    def _load_tm_dat(self):
        ...

    @property
    def transition_moments(self):
        """
        Computes the transition moments between wavefunctions stored in the object

        :return:
        :rtype:
        """
        ...

    @property
    def deperturbed_transition_moments(self):
        """
        Computes the transition moments between wavefunctions stored in the object

        :return:
        :rtype:
        """
        ...

    @property
    def transition_moment_corrections(self):
        """
        Computes the transition moment corrections between wavefunctions stored in the object

        :return:
        :rtype:
        """
        ...

    @property
    def deperturbed_transition_moment_corrections(self):
        """
        Computes the transition moment corrections between wavefunctions stored in the object

        :return:
        :rtype:
        """
        ...

    @property
    def oscillator_strengths(self):
        """
        Computes the oscillator strengths for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        ...

    @property
    def deperturbed_oscillator_strengths(self):
        """
        Computes the oscillator strengths for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        ...

    def oscillator_strengths_to_order(self, order):
        """

        :param tms:
        :type tms:
        :return:
        :rtype:
        """
        ...

    def deperturbed_oscillator_strengths_to_order(self, order):
        """

        :param tms:
        :type tms:
        :return:
        :rtype:
        """
        ...

    def _oscillator_strengths(self, tms):
        ...

    @property
    def intensities(self):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        ...

    @property
    def deperturbed_intensities(self):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        ...

    def intensities_to_order(self, order, return_freqs=False):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        ...

    def deperturbed_intensities_to_order(self, order, return_freqs=False):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        ...

    def _get_int(self, eng, oscs):
        ...

    def _dep_intensities(self, oscs, energy_order=None):
        ...

    def _intensities(self, oscs, energy_order=None):
        ...

    @property
    def zero_order_intensities(self):
        """
        Computes the harmonic intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        ...

    def prep_operator_terms(self, coeffs):
        ...

    def operator_correction_data(self, operator_coeffs, order=None):
        ...

    def generate_intensity_breakdown(self, include_wavefunctions=True):
        """
        Generates a breakdown of the terms that contribute to the intensity
        Returns in a format that can be directly exported to JSON if desired.

        :return:
        :rtype:
        """
        ...

    @classmethod
    def write_CSV_breakdown(cls, file, intensity_breakdown, padding=None):
        """
        Writes an intensity breakdown to a CSV by annoyingly flattening all the arrays

        :param file:
        :type file:
        :param intensity_breakdown:
        :type intensity_breakdown:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def _format_energies_table(cls, gs, zpe, states, freqs, real_fmt='{:>12.5f}'):
        ...

    def _format_group_energies_table(self, engs, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'):
        ...

    def format_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'):
        ...

    def format_deperturbed_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'):
        ...

    def format_property_matrices(self, states, prop_corrs, real_fmt='{:>.8e}', padding_fmt='{:>16}'):
        ...

    def format_dipole_contribs_tables(self):
        ...

    def format_deperturbed_dipole_contribs_tables(self):
        ...

    def format_energy_corrections_table(self, real_fmt='{:>12.5f}'):
        ...

    def _format_intensities_table(self, states, freqs, ints, harm_freqs, harm_ints, real_fmt='{:>12.5f}'):
        ...

    def _format_group_intensities_table(self, engs, ints, states=None, real_fmt='{:>12.5f}'):
        ...

    def format_intensities_table(self, real_fmt='{:>12.5f}'):
        ...

    def format_deperturbed_intensities_table(self, real_fmt='{:>12.5f}'):
        ...

    def _get_spectrum(self, engs, intensities):
        ...

    def get_spectrum(self):
        ...

    def get_deperturbed_spectrum(self):
        ...

    def _format_operator_table(self, states, val, zo_val, real_fmt='{:>12.5f}'):
        ...

    def _format_group_operator_table(self, operator_data, states=None, real_fmt='{:>12.5f}'):
        ...

    def format_operator_table(self, operators, real_fmt='{:>12.5f}'):
        ...