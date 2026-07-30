"""
Provides analyzer class to handle common VPT analyses
"""
import enum, weakref, functools, numpy as np, itertools as ip, io
import itertools
import os.path
import tempfile
from McUtils.Data import UnitsData
import McUtils.Numputils as nput
import McUtils.Devutils as dev
from McUtils.Scaffolding import Checkpointer, Logger, LogParser
from McUtils.Parsers import StringLineByLineReader
from ..Spectra import DiscreteSpectrum
from .Wavefunctions import PerturbationTheoryWavefunctions
from .Runner import VPTRunner
__all__ = ['VPTResultsLoader', 'VPTResultsSource', 'VPTAnalyzer']
__reload_hook__ = ['..Spectra', '.Runner', '.Wavefunctions']

class VPTResultsSource(enum.Enum):
    """
    Enum of sources to load PT results from
    """
    Wavefunctions = 'wavefunctions'
    Checkpoint = 'checkpoint'
    LogFile = 'log_file'

def property_dispatcher(basefn):
    """
    **LLM Docstring**

    Class decorator factory that turns a method into a source-dependent "virtual property": the decorated method's body is never actually run (it typically just raises to document the property's meaning), and calling the resulting bound method instead dispatches to whichever implementation was registered (via the returned `.register` decorator) for the current `self.res_type`.

    :param basefn: the base method being decorated, used only for its name/docstring and as the fallback if no implementation is registered for the current source
    :type basefn: callable
    :return: the dispatching wrapper function, with `.register` and `.registry` attached
    :rtype: callable
    """
    ...

class VPTAnalyzerLogParser(LogParser):

    def __init__(self, log_file, **opts):
        """
        **LLM Docstring**

        Set up a log-file parser for a VPT run's text log output, accepting either a file path or an in-memory `io.StringIO` (which is first written out to a temporary file, since the underlying `LogParser` needs a real file).

        :param log_file: the path to the log file, or an in-memory string buffer containing its contents
        :type log_file: str | io.StringIO
        :param opts: extra options forwarded to the base `LogParser.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        :raises ValueError: if `log_file` is neither a string path nor a `StringIO` buffer
        """
        ...

    @property
    def tree(self):
        """
        **LLM Docstring**

        The (cached) parsed block-tree structure of the log file, collapsed down to just the "Computing PT corrections:" subtree (or otherwise condensed) if the raw parse produced multiple top-level blocks.

        :return: the parsed (and condensed) log-file tree
        :rtype: object
        """
        ...

    class EnergiesBlockParser(StringLineByLineReader):

        def __init__(self, spec_str, **opts):
            """
            **LLM Docstring**

            Set up a line-by-line reader for an energies table block, tracking the label column width as lines are parsed.

            :param spec_str: the raw text of the energies block to parse
            :type spec_str: str
            :param opts: extra options forwarded to the base `StringLineByLineReader.__init__`
            :type opts: dict
            :return: None
            :rtype: None
            """
            ...

        def check_tag(self, line: str, depth: int=0, active_tag=None, label: str=None, history: list[str]=None):
            """
            **LLM Docstring**

            Skip blank lines, the table header row (starting with `"State"`), and deeply-indented continuation lines, treating everything else as a data row.

            :param line: the current line being classified
            :type line: str
            :param depth: the current nesting depth, unused
            :type depth: int
            :param active_tag: the currently active block tag, unused
            :type active_tag: object
            :param label: the current block label, unused
            :type label: str | None
            :param history: the tag history, unused
            :type history: list[str] | None
            :return: `LineReaderTags.SKIP` for lines that should be ignored, or `None` to treat the line as ordinary block content
            :rtype: object | None
            """
            ...

        def handle_block_line(self, label, line, depth=0, history: list[str]=None):
            """
            **LLM Docstring**

            Split a data row of the energies table into the state label and its numeric columns, inferring the label column width from the first parsed row.

            :param label: the block's label, unused directly
            :type label: str
            :param line: the raw data-row line to parse
            :type line: str
            :param depth: the current nesting depth, unused
            :type depth: int
            :param history: the tag history, unused
            :type history: list[str] | None
            :return: `(state_label, value_tokens)` -- the fixed-width state label substring and the whitespace-split remainder
            :rtype: tuple[str, list[str]]
            """
            ...

    @classmethod
    def reformat_eng_block(cls, sb):
        """
        **LLM Docstring**

        Convert a parsed energies block (a list of `(state_label, value_tokens)` rows) into parallel arrays of state labels, harmonic energies, and anharmonic energies, handling the alternate row format that includes a placeholder `'-'` degeneracy-group marker column.

        :param sb: the parsed block data, as returned by `EnergiesBlockParser`
        :type sb: list
        :return: `(states, harm, anh)` -- the state labels and their harmonic/anharmonic energy arrays
        :rtype: tuple[list[str], np.ndarray, np.ndarray]
        """
        ...

    def parse_energies_blocks(self, spec_str):
        """
        **LLM Docstring**

        Parse a raw energies-table block of log text into `(states, harmonic, anharmonic)` arrays, via `EnergiesBlockParser` and `reformat_eng_block`.

        :param spec_str: the raw text of the energies block
        :type spec_str: str
        :return: `(states, harm, anh)`
        :rtype: tuple[list[str], np.ndarray, np.ndarray]
        """
        ...

    @property
    def harmonic_energies(self):
        """
        **LLM Docstring**

        The `(states, harmonic_energies)` pair parsed from the log's "States Energies" (or, for a degenerate run, "Degenerate Energies") table block, cached after the first access.

        :return: `(states, harmonic_energies)`
        :rtype: tuple[list[str], np.ndarray]
        """
        ...

    @property
    def energies(self):
        """
        **LLM Docstring**

        The `(states, anharmonic_energies)` pair parsed from the log's "States Energies"/"Degenerate Energies" table block, cached after the first access.

        :return: `(states, anharmonic_energies)`
        :rtype: tuple[list[str], np.ndarray]
        """
        ...

    @property
    def zero_order_energies(self):
        """
        **LLM Docstring**

        The `(states, harmonic_energies)` pair parsed from the log's "States Energies"/"Degenerate Energies" table block; equivalent to `harmonic_energies`.

        :return: `(states, harmonic_energies)`
        :rtype: tuple[list[str], np.ndarray]
        """
        ...

    @property
    def deperturbed_energies(self):
        """
        **LLM Docstring**

        The `(states, deperturbed_energies)` pair parsed from the log's "Deperturbed Energies" table block (falling back to "States Energies" if no degenerate treatment was run), cached after the first access.

        :return: `(states, deperturbed_energies)`
        :rtype: tuple[list[str], np.ndarray]
        """
        ...

    @property
    def spectra(self):
        """
        **LLM Docstring**

        The (cached) IR spectrum data parsed from the log's "IR Data" block, via `parse_spectrum_blocks`.

        :return: the parsed per-initial-state spectrum data
        :rtype: list[dict]
        """
        ...

    @property
    def deperturbed_spectra(self):
        """
        **LLM Docstring**

        The (cached) deperturbed IR spectrum data parsed from the log's "Deperturbed IR Data" block, via `parse_spectrum_blocks`.

        :return: the parsed per-initial-state deperturbed spectrum data
        :rtype: list[dict]
        """
        ...

    class SpectrumBlockParser(StringLineByLineReader):

        def __init__(self, spec_str, **opts):
            """
            **LLM Docstring**

            Set up a line-by-line reader for an IR-spectrum table block, whose sub-blocks are delimited by `" Initial State:"` header lines.

            :param spec_str: the raw text of the spectrum block to parse
            :type spec_str: str
            :param opts: extra options forwarded to the base `StringLineByLineReader.__init__`
            :type opts: dict
            :return: None
            :rtype: None
            """
            ...

        def check_tag(self, line: str, depth: int=0, active_tag=None, label: str=None, history: list[str]=None):
            """
            **LLM Docstring**

            Skip blank lines, header rows, and deeply-indented continuation lines, and recognize `" Initial State:"` lines as the start of a new per-initial-state sub-block.

            :param line: the current line being classified
            :type line: str
            :param depth: the current nesting depth, unused
            :type depth: int
            :param active_tag: the currently active block tag, unused
            :type active_tag: object
            :param label: the current block label, unused
            :type label: str | None
            :param history: the tag history, unused
            :type history: list[str] | None
            :return: `LineReaderTags.SKIP`, a `(BLOCK_START, label, None)` triple for a new initial-state block, or `None` for ordinary content
            :rtype: object | tuple | None
            """
            ...

        def handle_block_line(self, label, line, depth=0, history: list[str]=None):
            """
            **LLM Docstring**

            Split a data row of a spectrum sub-block into the final-state label and its numeric columns, using the label's own width to determine the split point.

            :param label: the enclosing sub-block's initial-state label, used to determine the column width
            :type label: str
            :param line: the raw data-row line to parse
            :type line: str
            :param depth: the current nesting depth, unused
            :type depth: int
            :param history: the tag history, unused
            :type history: list[str] | None
            :return: `(state_label, value_tokens)`
            :rtype: tuple[str, list[str]]
            """
            ...

    @classmethod
    def reformat_spec_block(cls, sb):
        """
        **LLM Docstring**

        Convert a parsed spectrum block (per-initial-state groups of `(state_label, value_tokens)` rows) into a list of per-initial-state dicts holding the final-state labels and their harmonic/anharmonic `(frequency, intensity)` value pairs.

        :param sb: the parsed block data, as returned by `SpectrumBlockParser`
        :type sb: dict
        :return: a single dict (if there was only one initial state) or list of dicts, each with `'states'`, `'harmonic'`, and `'anharmonic'` entries
        :rtype: dict | list[dict]
        """
        ...

    def parse_spectrum_blocks(self, spec_str):
        """
        **LLM Docstring**

        Parse a raw multi-initial-state spectrum-table block of log text into a list of reformatted per-block spectrum dicts, via `SpectrumBlockParser` and `reformat_spec_block`.

        :param spec_str: the raw text of the spectrum block
        :type spec_str: str
        :return: the list of reformatted spectrum blocks
        :rtype: list
        """
        ...

    class TransitionMomentBlockParser(StringLineByLineReader):

        def __init__(self, spec_str, **opts):
            """
            **LLM Docstring**

            Set up a line-by-line reader for a transition-moment-correction table block, whose sub-blocks are delimited by `" Initial State:"` header lines.

            :param spec_str: the raw text of the transition-moment block to parse
            :type spec_str: str
            :param opts: extra options forwarded to the base `StringLineByLineReader.__init__`
            :type opts: dict
            :return: None
            :rtype: None
            """
            ...

        def check_tag(self, line: str, depth: int=0, active_tag=None, label: str=None, history: list[str]=None):
            """
            **LLM Docstring**

            Skip blank lines, header rows, and deeply-indented continuation lines, and recognize `" Initial State:"` lines as the start of a new per-initial-state sub-block.

            :param line: the current line being classified
            :type line: str
            :param depth: the current nesting depth, unused
            :type depth: int
            :param active_tag: the currently active block tag, unused
            :type active_tag: object
            :param label: the current block label, unused
            :type label: str | None
            :param history: the tag history, unused
            :type history: list[str] | None
            :return: `LineReaderTags.SKIP`, a `(BLOCK_START, label, None)` triple for a new initial-state block, or `None` for ordinary content
            :rtype: object | tuple | None
            """
            ...

        def handle_block_line(self, label, line, depth=0, history: list[str]=None):
            """
            **LLM Docstring**

            Split a data row of a transition-moment sub-block into the final-state label and its numeric correction columns, using the label's own width to determine the split point.

            :param label: the enclosing sub-block's initial-state label, used to determine the column width
            :type label: str
            :param line: the raw data-row line to parse
            :type line: str
            :param depth: the current nesting depth, unused
            :type depth: int
            :param history: the tag history, unused
            :type history: list[str] | None
            :return: `(state_label, value_tokens)`
            :rtype: tuple[str, list[str]]
            """
            ...
    _indexer = None

    @classmethod
    def load_term_counts(cls, nterms):
        """
        **LLM Docstring**

        Compute how many perturbative-order transition-moment correction terms exist below a given total column count, using a `SymmetricGroupGenerator`'s cumulative term totals as the source of per-order term counts.

        :param nterms: the total number of numeric columns present in a transition-moment data row
        :type nterms: int
        :return: the list of cumulative term counts below `nterms`, one entry per perturbative order
        :rtype: list[int]
        """
        ...

    @classmethod
    def reformat_tm_block(cls, sb):
        """
        **LLM Docstring**

        Convert a parsed transition-moment-correction block (per-initial-state groups of `(state_label, value_tokens)` rows) into a list of per-initial-state dicts holding the final-state labels and their per-order correction arrays (split from the flat value-token list using `load_term_counts`).

        :param sb: the parsed block data, as returned by `TransitionMomentBlockParser`
        :type sb: dict
        :return: a single dict (if there was only one initial state) or list of dicts, each with `'states'` and `'corrections'` (a list of per-order arrays) entries
        :rtype: dict | list[dict]
        """
        ...

    def parse_tm_blocks(self, spec_str):
        """
        **LLM Docstring**

        Parse a raw multi-initial-state transition-moment-table block of log text into a list of reformatted per-block correction dicts, via `TransitionMomentBlockParser` and `reformat_tm_block`.

        :param spec_str: the raw text of the transition-moment block
        :type spec_str: str
        :return: the list of reformatted transition-moment blocks
        :rtype: list
        """
        ...

    @property
    def transition_moment_corrections(self):
        """
        **LLM Docstring**

        The (cached) transition-dipole-moment corrections parsed from the log's per-axis "X/Y/Z Dipole Contributions" blocks, combined into a single list of per-initial-state dicts each holding the `[x, y, z]` correction data.

        :return: the combined per-initial-state transition-moment correction data
        :rtype: list[dict]
        """
        ...

    @property
    def deperturbed_transition_moment_corrections(self):
        """
        **LLM Docstring**

        The (cached) deperturbed transition-dipole-moment corrections parsed from the log's per-axis "X/Y/Z Deperturbed Dipole Contributions" blocks, combined into a single list of per-initial-state dicts each holding the `[x, y, z]` correction data.

        :return: the combined per-initial-state deperturbed transition-moment correction data
        :rtype: list[dict]
        """
        ...

class VPTResultsLoader:
    """
    Provides tools for loading results into canonical
    sources from a simulation, both from checkpoint files and from
    `PerturbationTheoryWavefunctions` objects and potentially more
    """

    def __init__(self, res, res_type=None):
        """
        :param res:
        :type res:
        :param res_type:
        :type res_type:
        """
        ...
    checkpoint_file_types = {'.npz', '.json', '.hdf5'}

    @classmethod
    def resolve_file_res_type(cls, res):
        """
        **LLM Docstring**

        Infer whether a given file path refers to a checkpoint file or a text log file, based on its file extension.

        :param res: the file path to classify
        :type res: str
        :return: the inferred result source
        :rtype: VPTResultsSource
        """
        ...

    def get_res_type(self, res):
        """
        :param res:
        :type res:
        :return:
        :rtype:
        """
        ...

    @property_dispatcher
    def potential_terms(self):
        """
        Returns the expansion of the potential

        :return:
        :rtype:
        """
        ...

    @potential_terms.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `potential_terms`: reads the stored potential expansion directly out of the checkpoint dict.

        :return: the potential-energy expansion terms
        :rtype: list
        """
        ...

    @potential_terms.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `potential_terms`: pulls the potential expansion off the live `PerturbationTheoryWavefunctions`' Hamiltonian.

        :return: the potential-energy expansion terms
        :rtype: list
        """
        ...

    @potential_terms.register('log_file')
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `potential_terms`: unsupported, since the raw potential expansion isn't written to the text log.

        :return: never returns
        :rtype: None
        :raises ValueError: always, noting potential terms aren't available from a log file
        """
        ...

    @property_dispatcher
    def kinetic_terms(self):
        """
        Returns the expansion of the kinetic energy

        :return:
        :rtype:
        """
        ...

    @kinetic_terms.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `kinetic_terms`: assembles the G-matrix, Coriolis, and pseudopotential expansions from the checkpoint dict, tolerating missing Coriolis/pseudopotential entries.

        :return: a dict with `'G_matrix'`, `'coriolis'`, and `'pseudopotential'` keys
        :rtype: dict
        """
        ...

    @kinetic_terms.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `kinetic_terms`: pulls the G-matrix and Coriolis expansions off the live `PerturbationTheoryWavefunctions`' Hamiltonian. Note the `'pseudopotential'` entry is populated from `coriolis_terms` again rather than a distinct pseudopotential-terms attribute, which looks like a copy-paste bug.

        :return: a dict with `'G_matrix'`, `'coriolis'`, and `'pseudopotential'` keys
        :rtype: dict
        """
        ...

    @kinetic_terms.register('log_file')
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `kinetic_terms`: unsupported, since the raw kinetic expansion isn't written to the text log.

        :return: never returns
        :rtype: None
        :raises ValueError: always, noting kinetic terms aren't available from a log file
        """
        ...

    @property_dispatcher
    def dipole_terms(self):
        """
        Returns the expansion of the dipole moment

        :return:
        :rtype:
        """
        ...

    @dipole_terms.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `dipole_terms`: reads the stored dipole expansion directly out of the checkpoint dict.

        :return: the dipole-moment expansion terms
        :rtype: list
        """
        ...

    @dipole_terms.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `dipole_terms`: pulls the dipole expansion off the live `PerturbationTheoryWavefunctions`.

        :return: the dipole-moment expansion terms
        :rtype: list
        """
        ...

    @dipole_terms.register('log_file')
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `dipole_terms`: unsupported, since the raw dipole expansion isn't written to the text log.

        :return: never returns
        :rtype: None
        :raises ValueError: always, noting dipole terms aren't available from a log file
        """
        ...

    @property_dispatcher
    def basis(self):
        """
        Returns the basis for the calculation

        :return:
        :rtype:
        """
        ...

    @basis.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `basis`: reads the total (coupled) basis state space out of the checkpoint's stored corrections dict.

        :return: the total basis state space
        :rtype: object
        """
        ...

    @basis.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `basis`: pulls the total basis off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the total basis state space
        :rtype: object
        """
        ...

    @basis.register('log_file')
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `basis`: unsupported, since the total basis isn't written to the text log.

        :return: never returns
        :rtype: None
        :raises ValueError: always, noting the total basis isn't available from a log file
        """
        ...

    @property_dispatcher
    def target_states(self):
        """
        Returns the target states for the calculation

        :return:
        :rtype:
        """
        ...

    @target_states.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `target_states`: reads the target state list out of the checkpoint's stored corrections dict.

        :return: the target states
        :rtype: object
        """
        ...

    @target_states.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `target_states`: pulls the target states' excitation vectors off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the target states' excitation vectors
        :rtype: np.ndarray
        """
        ...

    @property_dispatcher
    def spectrum(self):
        """
        Returns the IR spectrum calculated from perturbation theory

        :return:
        :rtype:
        """
        ...

    @spectrum.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `spectrum`: builds a `DiscreteSpectrum` from the frequency/intensity pair stored in the checkpoint, converting frequencies from Hartrees to wavenumbers.

        :return: the IR spectrum
        :rtype: DiscreteSpectrum
        """
        ...

    @spectrum.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `spectrum`: builds a `DiscreteSpectrum` directly from the live `PerturbationTheoryWavefunctions`' frequencies/intensities.

        :return: the IR spectrum
        :rtype: DiscreteSpectrum
        """
        ...

    @spectrum.register('log_file')
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `spectrum`: builds a `DiscreteSpectrum` from the anharmonic frequency/intensity columns of the parsed log's first spectrum block.

        :return: the IR spectrum
        :rtype: DiscreteSpectrum
        """
        ...

    @property_dispatcher
    def zero_order_spectrum(self):
        """
        Returns the zero-order IR spectrum calculated from perturbation theory

        :return:
        :rtype:
        """
        ...

    @zero_order_spectrum.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `zero_order_spectrum`: builds a `DiscreteSpectrum` from the zero-order frequency/intensity pair stored in the checkpoint, converting frequencies from Hartrees to wavenumbers.

        :return: the zero-order (harmonic) IR spectrum
        :rtype: DiscreteSpectrum
        """
        ...

    @zero_order_spectrum.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `zero_order_spectrum`: builds a `DiscreteSpectrum` from the live `PerturbationTheoryWavefunctions`' zeroth-order deperturbed frequencies/intensities.

        :return: the zero-order (harmonic) IR spectrum
        :rtype: DiscreteSpectrum
        """
        ...

    @zero_order_spectrum.register('log_file')
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `zero_order_spectrum`: builds a `DiscreteSpectrum` from the harmonic frequency/intensity columns of the parsed log's first spectrum block.

        :return: the zero-order (harmonic) IR spectrum
        :rtype: DiscreteSpectrum
        """
        ...

    @property_dispatcher
    def energy_corrections(self):
        """
        Returns the corrections to the energies

        :return:
        :rtype:
        """
        ...

    @energy_corrections.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `energy_corrections`: reads the per-order energy corrections out of the checkpoint's stored corrections dict.

        :return: the per-order energy corrections
        :rtype: np.ndarray
        """
        ...

    @energy_corrections.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `energy_corrections`: pulls the per-order energy corrections off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the per-order energy corrections
        :rtype: np.ndarray
        """
        ...

    def energies(self):
        """

        :return:
        :rtype:
        """
        ...

    def zero_order_energies(self):
        """

        :return:
        :rtype:
        """
        ...

    def deperturbed_energies(self):
        """

        :return:
        :rtype:
        """
        ...

    @property_dispatcher
    def wavefunction_corrections(self):
        """
        Returns the corrections to the wavefunctions

        :return:
        :rtype:
        """
        ...

    @wavefunction_corrections.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `wavefunction_corrections`: reads the per-order wavefunction corrections out of the checkpoint's stored corrections dict.

        :return: the per-order wavefunction corrections
        :rtype: list
        """
        ...

    @wavefunction_corrections.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `wavefunction_corrections`: pulls the per-order wavefunction corrections off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the per-order wavefunction corrections
        :rtype: list
        """
        ...

    @property_dispatcher
    def transition_moment_corrections(self):
        """
        Returns the corrections to the wavefunctions

        :return:
        :rtype:
        """
        ...

    @transition_moment_corrections.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `transition_moment_corrections`: reads the transition-moment corrections out of the checkpoint dict.

        :return: the transition-moment corrections
        :rtype: object
        """
        ...

    @transition_moment_corrections.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `transition_moment_corrections`: pulls the transition-moment corrections off the live `PerturbationTheoryWavefunctions`.

        :return: the transition-moment corrections
        :rtype: object
        """
        ...

    @transition_moment_corrections.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Duplicate wavefunction-backed implementation of `transition_moment_corrections`, registered under the same `'wavefunctions'` key as the implementation just above (and with identical behavior) -- this second registration simply overwrites the first with an equivalent function, so it looks like an accidental copy-paste duplicate rather than intentional.

        :return: the transition-moment corrections
        :rtype: object
        """
        ...

    def transition_moments(self):
        """
        **LLM Docstring**

        Sum the transition-dipole-moment corrections (over all perturbative orders) into a single `(nstates, 3)` array of final transition moments.

        :return: the summed transition-dipole moments
        :rtype: np.ndarray
        """
        ...

    @property_dispatcher
    def deperturbed_transition_moment_corrections(self):
        """
        Returns the corrections to the wavefunctions

        :return:
        :rtype:
        """
        ...

    @deperturbed_transition_moment_corrections.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `deperturbed_transition_moment_corrections`: reads the non-degenerate (deperturbed) transition-moment corrections out of the checkpoint dict.

        :return: the deperturbed transition-moment corrections
        :rtype: object
        """
        ...

    @deperturbed_transition_moment_corrections.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `deperturbed_transition_moment_corrections`: pulls the deperturbed transition-moment corrections off the live `PerturbationTheoryWavefunctions`.

        :return: the deperturbed transition-moment corrections
        :rtype: object
        """
        ...

    def deperturbed_transition_moments(self):
        """
        :return:
        :rtype:
        """
        ...

    @property_dispatcher
    def degenerate_states(self):
        """
        Returns the deperturbed states used to make the degenerate transform

        :return:
        :rtype:
        """
        ...

    @degenerate_states.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `degenerate_states`: reads the degenerate state groupings out of the checkpoint's stored degenerate-data dict.

        :return: the degenerate state groupings
        :rtype: list
        """
        ...

    @degenerate_states.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `degenerate_states`: pulls the degenerate state groupings off the live `PerturbationTheoryWavefunctions`' corrections, filtering out any singleton (non-degenerate) groups.

        :return: the (filtered) degenerate state groupings, or `None` if none are present
        :rtype: list | None
        """
        ...

    @property_dispatcher
    def deperturbed_hamiltonians(self):
        """
        Returns the deperturbed Hamiltonians used to make the degenerate transform

        :return:
        :rtype:
        """
        ...

    @deperturbed_hamiltonians.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `deperturbed_hamiltonians`: reads the degenerate-block Hamiltonians out of the checkpoint's stored degenerate-data dict.

        :return: the degenerate-block Hamiltonians
        :rtype: list
        """
        ...

    @deperturbed_hamiltonians.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `deperturbed_hamiltonians`: pulls the degenerate-block Hamiltonians off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the degenerate-block Hamiltonians
        :rtype: list
        """
        ...

    @property_dispatcher
    def degenerate_energies(self):
        """
        Returns the deperturbed states used to make the degenerate transform

        :return:
        :rtype:
        """
        ...

    @degenerate_energies.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `degenerate_energies`: reads the degenerate-corrected energies out of the checkpoint's stored degenerate-data dict.

        :return: the degenerate-corrected energies
        :rtype: np.ndarray
        """
        ...

    @degenerate_energies.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `degenerate_energies`: pulls the degenerate-corrected energies off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the degenerate-corrected energies
        :rtype: np.ndarray
        """
        ...

    @degenerate_energies.register('log_file')
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `degenerate_energies`: falls back to the parsed log's own `energies` (which already reflect any degenerate treatment applied during the run).

        :return: the degenerate-corrected energies
        :rtype: np.ndarray
        """
        ...

    @property_dispatcher
    def degenerate_rotations(self):
        """
        Returns the deperturbed states used to make the degenerate transform

        :return:
        :rtype:
        """
        ...

    @degenerate_rotations.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `degenerate_rotations`: reads the degenerate-perturbation-theory mixing/rotation matrices out of the checkpoint dict.

        :return: the degenerate mixing/rotation matrices
        :rtype: list
        """
        ...

    @degenerate_rotations.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `degenerate_rotations`: pulls the degenerate-perturbation-theory mixing/rotation matrices off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the degenerate mixing/rotation matrices
        :rtype: list
        """
        ...

    @property_dispatcher
    def log_file(self):
        """
        Returns the log_file for the run

        :return:
        :rtype:
        """
        ...
    log_file_extension = '.txt'

    @log_file.register('checkpoint')
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `log_file`: derives the expected log-file path from the checkpoint's own file path (swapping its extension for `log_file_extension`), returning `None` if no such file actually exists on disk.

        :return: the log-file path, or `None` if it doesn't exist
        :rtype: str | None
        """
        ...

    @log_file.register('wavefunctions')
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `log_file`: reads the log-file path off the live `PerturbationTheoryWavefunctions`' logger.

        :return: the log-file path
        :rtype: str | None
        """
        ...

    @log_file.register('log_file')
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `log_file`: the data itself already is the log file/parser, so it's returned unchanged.

        :return: `self.data`
        :rtype: VPTAnalyzerLogParser
        """
        ...

    @property
    def log_parser(self):
        """
        **LLM Docstring**

        A `VPTAnalyzerLogParser` for this loader's data: the data itself, if it already is one, otherwise a freshly built parser wrapping the resolved log file (or `None` if no log file is available).

        :return: the log parser, or `None` if unavailable
        :rtype: VPTAnalyzerLogParser | None
        """
        ...

def loaded_prop(fn):
    """
    **LLM Docstring**

    Decorator for `VPTAnalyzer` convenience properties that mirror a `VPTResultsLoader` method of the same name: copies that method's docstring onto the decorated function and wraps it as a `property`.

    :param fn: the function to wrap (typically calling the matching `VPTResultsLoader` method on `self.loader`)
    :type fn: callable
    :return: the wrapped property
    :rtype: property
    """
    ...

class VPTAnalyzer:
    """
    Provides analysis tools on VPT results
    """

    def __init__(self, res):
        """

        :param res:
        :type res:
        """
        ...

    @classmethod
    def run_VPT(cls, *args, logger=None, **kwargs):
        """
        Runs a VPT calculation through `VPTRunner.run_simple` and
        stores the output wave functions to use

        :param args:
        :type args:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def potential_terms(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def kinetic_terms(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def dipole_terms(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def basis(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def target_states(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def spectrum(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def energy_corrections(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def energies(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def frequencies(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def zero_order_spectrum(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def deperturbed_spectrum(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def deperturbed_frequencies(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def wavefunction_corrections(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def transition_moment_corrections(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def transition_moments(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def deperturbed_transition_moment_corrections(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def deperturbed_transition_moments(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def deperturbed_hamiltonians(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def zero_order_energies(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def deperturbed_energies(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def degenerate_states(self):
        """

        :return:
        :rtype:
        """
        ...

    @loaded_prop
    def degenerate_energies(self):
        """

        :return:
        :rtype:
        """
        ...

    def shift_and_transform_hamiltonian(self, hams, shifts):
        """

        :param hams:
        :type hams:
        :param shifts:
        :type shifts:
        :return:
        :rtype:
        """
        ...

    def get_shifted_transformed_transition_moments(self, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose'):
        """

        :param deg_states:
        :type deg_states:
        :param target_states:
        :type target_states:
        :param hams:
        :type hams:
        :param shifts:
        :type shifts:
        :param tmoms:
        :type tmoms:
        :param handling_mode:
        :type handling_mode:
        :return:
        :rtype:
        """
        ...

    def get_shifted_transformed_spectrum(self, zpe, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose'):
        """

        :param zpe:
        :type zpe:
        :param deg_states:
        :type deg_states:
        :param target_states:
        :type target_states:
        :param hams:
        :type hams:
        :param shifts:
        :type shifts:
        :param tmoms:
        :type tmoms:
        :param handling_mode:
        :type handling_mode:
        :return:
        :rtype:
        """
        ...

    def shifted_transformed_spectrum(self, deg_states, hams, shifts, return_transformation=False, handling_mode='transpose'):
        """

        :param deg_states:
        :type deg_states:
        :param hams:
        :type hams:
        :param shifts:
        :type shifts:
        :param return_transformation:
        :type return_transformation:
        :param handling_mode:
        :type handling_mode:
        :return:
        :rtype:
        """
        ...

    def transition_data(self, states, keys=['frequency', 'transition_moment', 'intensity'], data='deperturbed'):
        """

        :param states:
        :type states:
        :param keys:
        :type keys:
        :param data:
        :type data:
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def _aggregate_tmoms(tmom_corrs, inds, terms):
        """
        **LLM Docstring**

        Sum a selected subset of perturbative transition-moment-correction terms (identified by `(i, j, k)`-order-triple membership in `terms`) across all orders and all three Cartesian axes, for a given set of state indices.

        :param tmom_corrs: the per-axis, per-order transition-moment correction data
        :type tmom_corrs: list
        :param inds: the state indices to extract contributions for
        :type inds: np.ndarray
        :param terms: the set of `(i, j, k)` order-triples (mechanical/electrical/mixed-derivative order combination) to include in the sum
        :type terms: set[tuple] | list[tuple]
        :return: the summed `(len(inds), 3)` transition-moment contributions
        :rtype: np.ndarray
        """
        ...

    def transition_moment_term_sums(self, states, terms=None, rotation=None, data='deperturbed'):
        """

        :param states:
        :type states:
        :param terms:
        :type terms:
        :param rotation:
        :type rotation:
        :param data:
        :type data:
        :return:
        :rtype:
        """
        ...

    def transition_moment_term_sums_first_order(self, states, rotation=None, data='deperturbed'):
        """

        :param states:
        :type states:
        :param rotation:
        :type rotation:
        :param data:
        :type data:
        :return:
        :rtype:
        """
        ...

    def intensity_breakdown(self, states, terms=None, data='deperturbed'):
        """

        :param states:
        :type states:
        :param terms:
        :type terms:
        :param data:
        :type data:
        :return:
        :rtype:
        """
        ...

    def degenerate_coupling_element(self, state1, state2):
        """

        :param state1:
        :type state1:
        :param state2:
        :type state2:
        :return:
        :rtype:
        """
        ...

    def format_deperturbed_hamiltonian(self, which):
        """

        :param which:
        :type which:
        :return:
        :rtype:
        """
        ...

    @property
    def log_parser(self):
        """
        **LLM Docstring**

        The underlying `VPTAnalyzerLogParser` for this analyzer's loaded results, via `self.loader.log_parser`.

        :return: the log parser, or `None` if unavailable
        :rtype: VPTAnalyzerLogParser | None
        """
        ...

    def print_output_tables(self, print_energy_corrections=False, print_energies=False, print_transition_moments=False, print_intensities=True, **kwargs):
        """
        **LLM Docstring**

        Print the standard VPT output tables (energy corrections, energies, transition moments, intensities) for the loaded results: if the underlying data is a live `PerturbationTheoryWavefunctions` with a logger, delegates to `VPTRunner.print_output_tables` (temporarily disabling the logger so this call's own prints aren't duplicated into the log); otherwise, if the data came from a log file, simply reprints the "IR Data" block(s) from the parsed log.

        :param print_energy_corrections: whether to print the energy-correction breakdown table
        :type print_energy_corrections: bool
        :param print_energies: whether to print the energies table
        :type print_energies: bool
        :param print_transition_moments: whether to print the transition-moment breakdown table
        :type print_transition_moments: bool
        :param print_intensities: whether to print the IR intensity table
        :type print_intensities: bool
        :param kwargs: extra options forwarded to `VPTRunner.print_output_tables`
        :type kwargs: dict
        :return: None
        :rtype: None
        """
        ...