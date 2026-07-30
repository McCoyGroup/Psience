import os
import uuid
import numpy as np
import weakref
import tempfile as tf
import subprocess
import io
from McUtils.Jupyter import JHTML
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
import McUtils.Iterators as itut
from McUtils.ExternalPrograms import ExternalProgramRunner
__all__ = ['ezXML', 'ezFCFInterface']

class ezXML:

    class ezTag(JHTML.XML.DeclarativeElement):
        all_remapped_attrs = weakref.WeakKeyDictionary()

        def __init__(self, *elems, flag=None, **attrs):
            ...

        @classmethod
        def sanitize_value(cls, val):
            ...

        def tostring(self, renormalize_padding=True):
            ...

    class job_parameters(ezTag):
        known_opts = {'temperature', 'spectrum_intensity_threshold'}

    class do_not_excite_subspace(ezTag):
        known_opts = {'size', 'normal_modes'}

    class initial_state(ezTag):
        known_opts = {}

    class target_state(ezTag):
        known_opts = {}

    class the_only_initial_state(ezTag):
        known_opts = {}

    class energy_thresholds(ezTag):
        known_opts = {'units'}

    class print_franck_condon_matrices(ezTag):
        known_opts = {'flag'}

    class print_normal_modes(ezTag):
        known_opts = {'flag'}

    class parallel_approximation(ezTag):
        ...

    class max_vibr_to_store(ezTag):
        known_opts = {'target_el_state'}

    class single_excitation(ezTag):
        known_opts = {'ini', 'targ'}

    class dushinsky_rotations(ezTag):
        known_opts = {'target_state', 'max_vibr_excitations_in_initial_el_state', 'max_vibr_excitations_in_target_el_state'}

    class geometry(ezTag):
        known_opts = {'number_of_atoms', 'linear', 'units', 'text'}

        def __init__(self, *, number_of_atoms, text, linear='false', units='angstr', **opts):
            ...

    class manual_atoms_reordering(ezTag):
        known_opts = {'new_order'}

    class manual_coordinates_transformation(ezTag):
        known_opts = {'shift_along_x', 'shift_along_y', 'shift_along_z', 'rotate_around_x', 'rotate_around_y', 'rotate_around_z'}

    class normal_modes(ezTag):
        known_opts = {'if_mass_weighted', 'text', 'atoms'}

        def __init__(self, *, atoms, if_mass_weighted, text, **opts):
            ...

    class frequencies(ezTag):
        known_opts = {'text'}

    class excitation_energy(ezTag):
        known_opts = {'units'}

    class manual_normal_modes_reordering(ezTag):
        known_opts = {'new_order'}

    class input(ezTag):
        known_opts = {'job'}

class ezFCFInterface:

    def __init__(self, atoms, gs_nms, es_nms, excitations, masses=None, ground_states=None, include_rotation=True, rotation_order='gs', rotation_blocks=None, rotation_center=None, logger=None, mode_reordering=None, rotation_method='duschinsky', mass_weight=False, dimensionless=False, always_run_parallel=True, print_all=True, embed=False):
        ...

    @classmethod
    def _check_listable(cls, excitations):
        ...

    def format(self, job_type='harmonic_pes', temperature=0, spectrum_intensity_threshold=1e-08):
        ...

    def prep_masses(self, atoms, masses, units=None):
        ...

    def format_masses_file(self, atom_map):
        ...

    class _write_dir:

        def __init__(self, dir=None, dir_prefix=None, dir_suffix=None, delete=True):
            ...

        def __enter__(self):
            ...

        def __exit__(self, exc_type, exc_val, exc_tb):
            ...

    class ezFCFRunner(ExternalProgramRunner):
        default_opts = dict(prefix='ezFCF-', suffix='.xml')

        def __init__(self, job, binary, **opts):
            ...

        def prep_dir(self, dir):
            ...
        text_file_extensions = ExternalProgramRunner.text_file_extensions + ['.spectrum_parallel', '.spectrum_duschinsky']

    def run(self, ezFCF_binary, dir=None, dir_prefix=None, dir_suffix=None, mode='w', prefix='ezFCF-', suffix='.xml', delete=True, raise_errors=True, **job_opts):
        ...
    _exc_aliases = {'max_quanta': 'max_vibr_excitations_in_target_el_state'}
    _gsstate_aliases = {'max_quanta': 'max_vibr_excitations_in_initial_el_state'}

    @classmethod
    def canonicalize_excitation_options(cls, nms, threshold=None, fixed_modes=None, ground_states=None, **opts):
        ...

    @classmethod
    def prep_excitations(cls, nms, ground_states, excitations):
        """
        Dispatcher to get appropriate state spaces
        :param excitations:
        :param check:
        :return:
        """
        ...

    @classmethod
    def prep_state_str(cls, nms, state):
        ...

    def prep_parallel(self, *elems, rotation_order='gs', max_vibr_excitations_in_initial_el_state=0, max_vibr_excitations_in_target_el_state=0, **etc):
        ...

    def prep_duschinsky(self, *elems, rotation_order='gs', max_vibr_excitations_in_initial_el_state=0, max_vibr_excitations_in_target_el_state=0, **etc):
        ...

    @classmethod
    def prep_initial_state(cls, atoms, nms, autoembed=False, excitation_energy_units='cm-1'):
        ...

    @classmethod
    def prep_target_state(cls, atoms, nms, autoembed=False, mode_reordering=None, excitation_energy_units='cm-1'):
        """

        :param atoms:
        :param nms:
        :param excitation_energy_units:
        :return:
        """
        ...

    @classmethod
    def prep_state(cls, tag, atoms, nms, *elems):
        ...

    @classmethod
    def prep_geometry(cls, atoms, nms, linear=False, units='BohrRadius'):
        ...

    @classmethod
    def format_modes_block(cls, modes):
        ...

    @classmethod
    def parse_modes_block(cls, modes_str, num_atoms=None):
        ...

    @classmethod
    def prep_normal_modes(cls, atoms, nms):
        ...

    @classmethod
    def format_freqs_block(cls, freqs):
        ...

    @classmethod
    def parse_freqs_block(cls, freqs):
        ...

    @classmethod
    def prep_frequencies(cls, nms):
        ...

    @classmethod
    def parse_state(self, state_xml):
        ...

    @classmethod
    def parse_fc_model(cls, input_xml, **opts):
        ...