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

__all__ = [
    "VPTResultsLoader",
    "VPTResultsSource",
    "VPTAnalyzer"
]

__reload_hook__ = ["..Spectra", ".Runner", ".Wavefunctions"]

class VPTResultsSource(enum.Enum):
    """
    Enum of sources to load PT results from
    """
    Wavefunctions = "wavefunctions"
    Checkpoint = "checkpoint"
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
    registry = {}#weakref.WeakValueDictionary()
    def register(key, registry=registry):
        """
        **LLM Docstring**

        Build a decorator that registers an implementation function under a given result-source key, for later dispatch by the enclosing `property_dispatcher`.

        :param key: the result-source key to register under
        :type key: str | VPTResultsSource
        :param registry: the shared registry dict to register into, captured from the enclosing `property_dispatcher` call
        :type registry: dict
        :return: the `bind` decorator used to actually register a function
        :rtype: callable
        """
        if not isinstance(key, VPTResultsSource):
            key = VPTResultsSource(key)
        def bind(fn, key=key, registry=registry):
            """
            **LLM Docstring**

            Register `fn` in the registry under the enclosing `key`, copying the base method's name/docstring onto it so it presents consistently regardless of which implementation ends up being called.

            :param fn: the source-specific implementation function to register
            :type fn: callable
            :param key: the result-source key to register under, from the enclosing scope
            :type key: VPTResultsSource
            :param registry: the shared registry dict, from the enclosing scope
            :type registry: dict
            :return: `fn`, unchanged (returned so it remains usable/reassignable as a plain method)
            :rtype: callable
            """
            fn.__name__ = basefn.__name__
            fn.__doc__ = basefn.__doc__
            registry[key] = fn
            return fn
        return bind

    @functools.wraps(basefn)
    def dispatch(self, *args, **kwargs):
        """
        **LLM Docstring**

        The actual dispatching call: looks up and invokes whichever implementation was registered for `self.res_type`.

        :param self: the `VPTResultsLoader` instance the property is being accessed on
        :type self: VPTResultsLoader
        :param args: positional arguments forwarded to the resolved implementation
        :type args: tuple
        :param kwargs: keyword arguments forwarded to the resolved implementation
        :type kwargs: dict
        :return: the resolved implementation's return value
        :rtype: object
        :raises KeyError: if no implementation is registered for `self.res_type`
        """
        return registry[self.res_type](self, *args, **kwargs)
    dispatch.registry = registry
    dispatch.register = register

    return dispatch

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
        if isinstance(log_file, str):
            lf = log_file
        elif isinstance(log_file, io.StringIO):
            with tempfile.NamedTemporaryFile(delete=False) as tf:
                lf = tf.name
            with open(lf, 'w+') as tf:
                log_file.seek(0)
                tf.write(log_file.read())
        else:
            lf = None

        if lf is None:
            raise ValueError("log file {} can't be parsed".format(log_file))

        super().__init__(lf, **opts)
        self._tree = None
        self._energies = None
        self._deperturbed_energies = None
        self._spectra = None
        self._deperturbed_spectra = None
        self._tms = None
        self._deperturbed_tms = None

    @property
    def tree(self):
        """
        **LLM Docstring**

        The (cached) parsed block-tree structure of the log file, collapsed down to just the "Computing PT corrections:" subtree (or otherwise condensed) if the raw parse produced multiple top-level blocks.

        :return: the parsed (and condensed) log-file tree
        :rtype: object
        """
        if self._tree is None:
            with self:
                self._tree = self.to_tree(depth=-1)
                if len(self._tree) > 1:
                    if (
                                    self._tree.keys() is not None
                                    and list(self._tree.keys())[0] != "Computing PT corrections:"
                        ):
                        #TODO: decide how a TreeWrapper will handle this kind of slice operation
                        self._tree = type(self._tree)(self._tree[-1])
                    else:
                        self._tree = self._tree.condense_subtrees()
        return self._tree

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
            self.label_width = None
            super().__init__(spec_str, max_nesting_depth=0, **opts)
        def check_tag(self, line:str, depth: int = 0, active_tag=None, label: str = None, history: list[str] = None):
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
            ls = line.lstrip()
            if len(ls) == 0 or ls.startswith("State") or line.startswith(" "*5):
                return self.LineReaderTags.SKIP
        def handle_block_line(self, label, line, depth=0, history:list[str]=None):
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
            if self.label_width is None:
                self.label_width = len(line.split("  ", 1)[0])
            n = self.label_width
            line = line.strip()
            state = line[:n]
            return state, line[n:].split()

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
        states = []
        sb = sb[0]
        harm = np.empty((len(sb),), dtype=float)
        anh = np.empty((len(sb),), dtype=float)
        for i,(state, sublist) in enumerate(sb):
            states.append(state)
            if sublist[0] == '-':
                harm[i] = float(sublist[2])
                anh[i] = float(sublist[3])
            else:
                harm[i] = float(sublist[0])
                anh[i] = float(sublist[1])
        return states, harm, anh

    def parse_energies_blocks(self, spec_str):
        """
        **LLM Docstring**

        Parse a raw energies-table block of log text into `(states, harmonic, anharmonic)` arrays, via `EnergiesBlockParser` and `reformat_eng_block`.

        :param spec_str: the raw text of the energies block
        :type spec_str: str
        :return: `(states, harm, anh)`
        :rtype: tuple[list[str], np.ndarray, np.ndarray]
        """
        with self.EnergiesBlockParser(spec_str) as parser:
            base_specs = list(parser)

        return self.reformat_eng_block(base_specs)
    @property
    def harmonic_energies(self):
        """
        **LLM Docstring**

        The `(states, harmonic_energies)` pair parsed from the log's "States Energies" (or, for a degenerate run, "Degenerate Energies") table block, cached after the first access.

        :return: `(states, harmonic_energies)`
        :rtype: tuple[list[str], np.ndarray]
        """
        if self._energies is None:
            try:
                eng_str = self.tree["States Energies"]
            except IndexError:
                eng_str = self.tree[0, "Degenerate Energies"]
            self._energies = self.parse_energies_blocks(eng_str)
        return self._energies[0], self._energies[1]

    @property
    def energies(self):
        """
        **LLM Docstring**

        The `(states, anharmonic_energies)` pair parsed from the log's "States Energies"/"Degenerate Energies" table block, cached after the first access.

        :return: `(states, anharmonic_energies)`
        :rtype: tuple[list[str], np.ndarray]
        """
        if self._energies is None:
            try:
                eng_str = self.tree["States Energies"]
            except IndexError:
                eng_str = self.tree["Degenerate Energies"]
            self._energies = self.parse_energies_blocks(eng_str)
        return self._energies[0], self._energies[2]

    @property
    def zero_order_energies(self):
        """
        **LLM Docstring**

        The `(states, harmonic_energies)` pair parsed from the log's "States Energies"/"Degenerate Energies" table block; equivalent to `harmonic_energies`.

        :return: `(states, harmonic_energies)`
        :rtype: tuple[list[str], np.ndarray]
        """
        if self._energies is None:
            try:
                eng_str = self.tree["States Energies"]
            except IndexError:
                eng_str = self.tree["Degenerate Energies"]
            self._energies = self.parse_energies_blocks(eng_str)
        return self._energies[0], self._energies[1]

    @property
    def deperturbed_energies(self):
        """
        **LLM Docstring**

        The `(states, deperturbed_energies)` pair parsed from the log's "Deperturbed Energies" table block (falling back to "States Energies" if no degenerate treatment was run), cached after the first access.

        :return: `(states, deperturbed_energies)`
        :rtype: tuple[list[str], np.ndarray]
        """
        if self._deperturbed_energies is None:
            try:
                eng_str = self.tree["Deperturbed Energies"]
            except IndexError:
                eng_str = self.tree["States Energies"]
            self._deperturbed_energies = self.parse_energies_blocks(eng_str)
        return self._deperturbed_energies[0], self._deperturbed_energies[2]

    @property
    def spectra(self):
        """
        **LLM Docstring**

        The (cached) IR spectrum data parsed from the log's "IR Data" block, via `parse_spectrum_blocks`.

        :return: the parsed per-initial-state spectrum data
        :rtype: list[dict]
        """
        if self._spectra is None:
            self._spectra = self.parse_spectrum_blocks(self.tree["IR Data"])
        return self._spectra

    @property
    def deperturbed_spectra(self):
        """
        **LLM Docstring**

        The (cached) deperturbed IR spectrum data parsed from the log's "Deperturbed IR Data" block, via `parse_spectrum_blocks`.

        :return: the parsed per-initial-state deperturbed spectrum data
        :rtype: list[dict]
        """
        if self._deperturbed_spectra is None:
            self._deperturbed_spectra = self.parse_spectrum_blocks(self.tree["Deperturbed IR Data"])
        return self._deperturbed_spectra

    # @property
    # def harmonic_spectra(self):
    #     if self._spectra is None:
    #
    #         self._harmonic_spectra, self._spectra = self.parse_spectrum_blocks(
    #             self.tree[0, "IR Data"]
    #         )
    #     return self._spectra

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
            super().__init__(spec_str, max_nesting_depth=0, **opts)
        def check_tag(self, line:str, depth: int = 0, active_tag=None, label: str = None, history: list[str] = None):
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
            block_tag = ' Initial State:'
            if len(line) == 0 or line.startswith("State") or line.startswith(" "*5):
                return self.LineReaderTags.SKIP
            elif line.startswith(block_tag):
                return self.LineReaderTags.BLOCK_START, line[len(block_tag):].strip(), None
        def handle_block_line(self, label, line, depth=0, history:list[str]=None):
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
            n = len(label)
            line = line.strip()
            state = line[:n]
            return state, line[n:].split()

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
        res = []
        for init, sublist in sb.items():
            harm = np.empty((len(sublist), 2), dtype=float)
            anh = np.empty((len(sublist), 2), dtype=float)
            states = []
            for i,(state, vals) in enumerate(sublist):
                states.append(state)
                harm[i][0] = float(vals[0])
                harm[i][1] = float(vals[1])
                anh[i][0] = float(vals[2])
                anh[i][1] = float(vals[3])

            res.append({
                'states':states,
                'harmonic':harm,
                'anharmonic':anh
            })
        if len(res) == 1:
            res = res[0]
        return res
    def parse_spectrum_blocks(self, spec_str):
        """
        **LLM Docstring**

        Parse a raw multi-initial-state spectrum-table block of log text into a list of reformatted per-block spectrum dicts, via `SpectrumBlockParser` and `reformat_spec_block`.

        :param spec_str: the raw text of the spectrum block
        :type spec_str: str
        :return: the list of reformatted spectrum blocks
        :rtype: list
        """
        with self.SpectrumBlockParser(spec_str) as parser:
            base_specs = list(parser)

        return [self.reformat_spec_block(sb) for sb in base_specs]

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
            super().__init__(spec_str, max_nesting_depth=0, **opts)
        def check_tag(self, line:str, depth: int = 0, active_tag=None, label: str = None, history: list[str] = None):
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
            block_tag = ' Initial State:'
            if len(line) == 0 or line.startswith("State") or line.startswith(" "*5):
                return self.LineReaderTags.SKIP
            elif line.startswith(block_tag):
                return self.LineReaderTags.BLOCK_START, line[len(block_tag):].strip(), None
        def handle_block_line(self, label, line, depth=0, history:list[str]=None):
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
            n = len(label)
            line = line.strip()
            state = line[:n]
            return state, line[n:].split()

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
        import McUtils.Combinatorics as comb

        if cls._indexer is None:
            cls._indexer = comb.SymmetricGroupGenerator(3)

        cls._indexer.load_to_size(nterms)
        return [c for c in cls._indexer._cumtotals if c < nterms]

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
        res = []
        for init, sublist in sb.items():
            counts = cls.load_term_counts(len(sublist[0][1]))
            blocks = [
                np.empty((len(sublist), n), dtype=float)
                for n in counts
            ]
            states = []
            for i,(state, vals) in enumerate(sublist):
                states.append(state)
                s = 0
                for n,p in enumerate(counts):
                    blocks[n][i, :] = [float(v) for v in vals[s:s+p]]
                    s += p
            res.append({
                'states':states,
                'corrections':blocks
            })
        if len(res) == 1:
            res = res[0]
        return res
    def parse_tm_blocks(self, spec_str):
        """
        **LLM Docstring**

        Parse a raw multi-initial-state transition-moment-table block of log text into a list of reformatted per-block correction dicts, via `TransitionMomentBlockParser` and `reformat_tm_block`.

        :param spec_str: the raw text of the transition-moment block
        :type spec_str: str
        :return: the list of reformatted transition-moment blocks
        :rtype: list
        """
        with self.TransitionMomentBlockParser(spec_str) as parser:
            base_specs = list(parser)

        return [self.reformat_tm_block(sb) for sb in base_specs]

    @property
    def transition_moment_corrections(self):
        """
        **LLM Docstring**

        The (cached) transition-dipole-moment corrections parsed from the log's per-axis "X/Y/Z Dipole Contributions" blocks, combined into a single list of per-initial-state dicts each holding the `[x, y, z]` correction data.

        :return: the combined per-initial-state transition-moment correction data
        :rtype: list[dict]
        """
        if self._tms is None:
            blocks = [
                self.parse_tm_blocks(self.tree["X Dipole Contributions"]),
                self.parse_tm_blocks(self.tree["Y Dipole Contributions"]),
                self.parse_tm_blocks(self.tree["Z Dipole Contributions"]),
            ]
            self._tms = [
                {
                    'states':x['states'],
                    'corrections':[
                        x, y, z
                    ]
                }
                for x,y,z in zip(*blocks)
            ]
        return self._tms

    @property
    def deperturbed_transition_moment_corrections(self):
        """
        **LLM Docstring**

        The (cached) deperturbed transition-dipole-moment corrections parsed from the log's per-axis "X/Y/Z Deperturbed Dipole Contributions" blocks, combined into a single list of per-initial-state dicts each holding the `[x, y, z]` correction data.

        :return: the combined per-initial-state deperturbed transition-moment correction data
        :rtype: list[dict]
        """
        if self._deperturbed_tms is None:
            blocks = [
                self.parse_tm_blocks(self.tree["X Deperturbed Dipole Contributions"]),
                self.parse_tm_blocks(self.tree["Y Deperturbed Dipole Contributions"]),
                self.parse_tm_blocks(self.tree["Z Deperturbed Dipole Contributions"]),
            ]
            self._deperturbed_tms = [
                {
                    'states': x['states'],
                    'corrections': [
                        x, y, z
                    ]
                }
                for x, y, z in zip(*blocks)
            ]
        return self._deperturbed_tms


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
        if isinstance(res_type, str):
            res_type = VPTResultsSource(res_type)
        if isinstance(res, str):
            if res_type is None:
                res_type = self.resolve_file_res_type(res)
            if res_type == VPTResultsSource.Checkpoint:
                res = Checkpointer.from_file(res)
            else:
                res = VPTAnalyzerLogParser(res)
        self.data = res
        self.res_type = self.get_res_type(res) if res_type is None else res_type

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
        _, ext = os.path.splitext(res)
        if ext in cls.checkpoint_file_types:
            return VPTResultsSource.Checkpoint
        else:
            return VPTResultsSource.LogFile

    def get_res_type(self, res):
        """
        :param res:
        :type res:
        :return:
        :rtype:
        """
        if isinstance(res, PerturbationTheoryWavefunctions):
            return VPTResultsSource.Wavefunctions
        elif isinstance(res, Checkpointer):
            return VPTResultsSource.Checkpoint
        else:
            raise ValueError("do not know how to load PT results from {}".format(
                res
            ))

    @property_dispatcher
    def potential_terms(self):
        """
        Returns the expansion of the potential

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @potential_terms.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `potential_terms`: reads the stored potential expansion directly out of the checkpoint dict.

        :return: the potential-energy expansion terms
        :rtype: list
        """
        return self.data["potential_terms"]
    @potential_terms.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `potential_terms`: pulls the potential expansion off the live `PerturbationTheoryWavefunctions`' Hamiltonian.

        :return: the potential-energy expansion terms
        :rtype: list
        """
        return self.data.corrs.hamiltonian.V_terms
    @potential_terms.register("log_file")
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `potential_terms`: unsupported, since the raw potential expansion isn't written to the text log.

        :return: never returns
        :rtype: None
        :raises ValueError: always, noting potential terms aren't available from a log file
        """
        raise ValueError("potential terms not in log file")

    @property_dispatcher
    def kinetic_terms(self):
        """
        Returns the expansion of the kinetic energy

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @kinetic_terms.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `kinetic_terms`: assembles the G-matrix, Coriolis, and pseudopotential expansions from the checkpoint dict, tolerating missing Coriolis/pseudopotential entries.

        :return: a dict with `'G_matrix'`, `'coriolis'`, and `'pseudopotential'` keys
        :rtype: dict
        """
        try:
            corr = self.data["coriolis_terms"]
        except KeyError:
            corr = None
        try:
            u = self.data["pseudopotential_terms"]
        except KeyError:
            u = None
        return {
            "G_matrix": self.data["gmatrix_terms"],
            "coriolis": corr,
            "pseudopotential": u
        }
    @kinetic_terms.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `kinetic_terms`: pulls the G-matrix and Coriolis expansions off the live `PerturbationTheoryWavefunctions`' Hamiltonian. Note the `'pseudopotential'` entry is populated from `coriolis_terms` again rather than a distinct pseudopotential-terms attribute, which looks like a copy-paste bug.

        :return: a dict with `'G_matrix'`, `'coriolis'`, and `'pseudopotential'` keys
        :rtype: dict
        """
        return {
            "G_matrix": self.data.corrs.hamiltonian.G_terms,
            "coriolis": self.data.corrs.hamiltonian.coriolis_terms,
            "pseudopotential": self.data.corrs.hamiltonian.coriolis_terms
        }
    @kinetic_terms.register("log_file")
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `kinetic_terms`: unsupported, since the raw kinetic expansion isn't written to the text log.

        :return: never returns
        :rtype: None
        :raises ValueError: always, noting kinetic terms aren't available from a log file
        """
        raise ValueError("kinetic terms not in log file")

    @property_dispatcher
    def dipole_terms(self):
        """
        Returns the expansion of the dipole moment

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @dipole_terms.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `dipole_terms`: reads the stored dipole expansion directly out of the checkpoint dict.

        :return: the dipole-moment expansion terms
        :rtype: list
        """
        return self.data["dipole_terms"]
    @dipole_terms.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `dipole_terms`: pulls the dipole expansion off the live `PerturbationTheoryWavefunctions`.

        :return: the dipole-moment expansion terms
        :rtype: list
        """
        return self.data.dipole_terms
    @dipole_terms.register("log_file")
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `dipole_terms`: unsupported, since the raw dipole expansion isn't written to the text log.

        :return: never returns
        :rtype: None
        :raises ValueError: always, noting dipole terms aren't available from a log file
        """
        raise ValueError("dipole terms not in log file")

    @property_dispatcher
    def basis(self):
        """
        Returns the basis for the calculation

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @basis.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `basis`: reads the total (coupled) basis state space out of the checkpoint's stored corrections dict.

        :return: the total basis state space
        :rtype: object
        """
        return self.data["corrections"]["total_states"]
    @basis.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `basis`: pulls the total basis off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the total basis state space
        :rtype: object
        """
        data = self.data #type:PerturbationTheoryWavefunctions
        return data.corrs.total_basis
    @basis.register("log_file")
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `basis`: unsupported, since the total basis isn't written to the text log.

        :return: never returns
        :rtype: None
        :raises ValueError: always, noting the total basis isn't available from a log file
        """
        raise ValueError("total basis not in log file")

    @property_dispatcher
    def target_states(self):
        """
        Returns the target states for the calculation

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @target_states.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `target_states`: reads the target state list out of the checkpoint's stored corrections dict.

        :return: the target states
        :rtype: object
        """
        return self.data["corrections"]["states"]
    @target_states.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `target_states`: pulls the target states' excitation vectors off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the target states' excitation vectors
        :rtype: np.ndarray
        """
        data = self.data  # type:PerturbationTheoryWavefunctions
        return data.corrs.states.excitations

    @property_dispatcher
    def spectrum(self):
        """
        Returns the IR spectrum calculated from perturbation theory

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @spectrum.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `spectrum`: builds a `DiscreteSpectrum` from the frequency/intensity pair stored in the checkpoint, converting frequencies from Hartrees to wavenumbers.

        :return: the IR spectrum
        :rtype: DiscreteSpectrum
        """
        freq, ints = self.data["spectrum"]
        return DiscreteSpectrum(freq * UnitsData.convert("Hartrees", "Wavenumbers"), ints)
    @spectrum.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `spectrum`: builds a `DiscreteSpectrum` directly from the live `PerturbationTheoryWavefunctions`' frequencies/intensities.

        :return: the IR spectrum
        :rtype: DiscreteSpectrum
        """
        return DiscreteSpectrum(self.data.frequencies() * UnitsData.convert("Hartrees", "Wavenumbers"), self.data.intensities[1:])
    @spectrum.register("log_file")
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `spectrum`: builds a `DiscreteSpectrum` from the anharmonic frequency/intensity columns of the parsed log's first spectrum block.

        :return: the IR spectrum
        :rtype: DiscreteSpectrum
        """
        freq, ints = self.data.spectra[0]["anharmonic"].T
        return DiscreteSpectrum(freq, ints)

    @property_dispatcher
    def zero_order_spectrum(self):
        """
        Returns the zero-order IR spectrum calculated from perturbation theory

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @zero_order_spectrum.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `zero_order_spectrum`: builds a `DiscreteSpectrum` from the zero-order frequency/intensity pair stored in the checkpoint, converting frequencies from Hartrees to wavenumbers.

        :return: the zero-order (harmonic) IR spectrum
        :rtype: DiscreteSpectrum
        """
        freq, ints = self.data["zero_order_spectrum"]
        return DiscreteSpectrum(freq * UnitsData.convert("Hartrees", "Wavenumbers"), ints)
    @zero_order_spectrum.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `zero_order_spectrum`: builds a `DiscreteSpectrum` from the live `PerturbationTheoryWavefunctions`' zeroth-order deperturbed frequencies/intensities.

        :return: the zero-order (harmonic) IR spectrum
        :rtype: DiscreteSpectrum
        """
        return DiscreteSpectrum(self.data.deperturbed_frequencies(order=0) * UnitsData.convert("Hartrees", "Wavenumbers"),
                                self.data.deperturbed_intensities_to_order(order=0)[1:]
                                )
    @zero_order_spectrum.register("log_file")
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `zero_order_spectrum`: builds a `DiscreteSpectrum` from the harmonic frequency/intensity columns of the parsed log's first spectrum block.

        :return: the zero-order (harmonic) IR spectrum
        :rtype: DiscreteSpectrum
        """
        freq, ints = self.data.spectra[0]["harmonic"].T
        return DiscreteSpectrum(freq, ints)

    @property_dispatcher
    def energy_corrections(self):
        """
        Returns the corrections to the energies

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @energy_corrections.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `energy_corrections`: reads the per-order energy corrections out of the checkpoint's stored corrections dict.

        :return: the per-order energy corrections
        :rtype: np.ndarray
        """
        return self.data["corrections"]["energies"]
    @energy_corrections.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `energy_corrections`: pulls the per-order energy corrections off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the per-order energy corrections
        :rtype: np.ndarray
        """
        data = self.data #type: PerturbationTheoryWavefunctions
        return data.corrs.energy_corrs

    def energies(self):
        """

        :return:
        :rtype:
        """
        if hasattr(self.data, 'energies'):
            return self.data.energies
        else:
            return np.sum(self.energy_corrections(), axis=1)
    def zero_order_energies(self):
        """

        :return:
        :rtype:
        """
        if hasattr(self.data, 'zero_order_energies'):
            return self.data.zero_order_energies
        else:
            return self.energy_corrections()[:, 0]
    def deperturbed_energies(self):
        """

        :return:
        :rtype:
        """
        if hasattr(self.data, 'deperturbed_energies'):
            return self.data.deperturbed_energies
        else:
            return np.sum(self.energy_corrections(), axis=1)

    @property_dispatcher
    def wavefunction_corrections(self):
        """
        Returns the corrections to the wavefunctions

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @wavefunction_corrections.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `wavefunction_corrections`: reads the per-order wavefunction corrections out of the checkpoint's stored corrections dict.

        :return: the per-order wavefunction corrections
        :rtype: list
        """
        return self.data["corrections"]["wavefunctions"]
    @wavefunction_corrections.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `wavefunction_corrections`: pulls the per-order wavefunction corrections off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the per-order wavefunction corrections
        :rtype: list
        """
        return self.data.corrs.wfn_corrections

    @property_dispatcher
    def transition_moment_corrections(self):
        """
        Returns the corrections to the wavefunctions

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @transition_moment_corrections.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `transition_moment_corrections`: reads the transition-moment corrections out of the checkpoint dict.

        :return: the transition-moment corrections
        :rtype: object
        """
        return self.data["transition_moments"]
    @transition_moment_corrections.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `transition_moment_corrections`: pulls the transition-moment corrections off the live `PerturbationTheoryWavefunctions`.

        :return: the transition-moment corrections
        :rtype: object
        """
        return self.data.transition_moment_corrections
    @transition_moment_corrections.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Duplicate wavefunction-backed implementation of `transition_moment_corrections`, registered under the same `'wavefunctions'` key as the implementation just above (and with identical behavior) -- this second registration simply overwrites the first with an equivalent function, so it looks like an accidental copy-paste duplicate rather than intentional.

        :return: the transition-moment corrections
        :rtype: object
        """
        return self.data.transition_moment_corrections
    def transition_moments(self):
        """
        **LLM Docstring**

        Sum the transition-dipole-moment corrections (over all perturbative orders) into a single `(nstates, 3)` array of final transition moments.

        :return: the summed transition-dipole moments
        :rtype: np.ndarray
        """
        transition_moment_components = self.transition_moment_corrections()
        order = len(transition_moment_components[0])
        tmom = np.array([
            sum(
                sum(ip.chain(transition_moment_components[i][j]))
                if not isinstance(transition_moment_components[i][j], (float, int, np.floating, np.integer))
                else transition_moment_components[i][j]
                for i in range(order)  # correction order
            ) for j in range(3)  # xyz
        ])

        return tmom

    @property_dispatcher
    def deperturbed_transition_moment_corrections(self):
        """
        Returns the corrections to the wavefunctions

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @deperturbed_transition_moment_corrections.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `deperturbed_transition_moment_corrections`: reads the non-degenerate (deperturbed) transition-moment corrections out of the checkpoint dict.

        :return: the deperturbed transition-moment corrections
        :rtype: object
        """
        return self.data["nondegenerate_transition_moments"]
    @deperturbed_transition_moment_corrections.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `deperturbed_transition_moment_corrections`: pulls the deperturbed transition-moment corrections off the live `PerturbationTheoryWavefunctions`.

        :return: the deperturbed transition-moment corrections
        :rtype: object
        """
        return self.data.deperturbed_transition_moment_corrections
    def deperturbed_transition_moments(self):
        """
        :return:
        :rtype:
        """
        transition_moment_components = self.deperturbed_transition_moment_corrections()
        order = len(transition_moment_components[0])
        tmom = np.array([
            sum(
                sum(ip.chain(transition_moment_components[i][j]))
                if not isinstance(transition_moment_components[i][j], (float, int, np.floating, np.integer))
                else transition_moment_components[i][j]
                for i in range(order)  # correction order
            ) for j in range(3)  # xyz
        ])

        return tmom

    @property_dispatcher
    def degenerate_states(self):
        """
        Returns the deperturbed states used to make the degenerate transform

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @degenerate_states.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `degenerate_states`: reads the degenerate state groupings out of the checkpoint's stored degenerate-data dict.

        :return: the degenerate state groupings
        :rtype: list
        """
        return self.data["degenerate_data"]["states"]
    @degenerate_states.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `degenerate_states`: pulls the degenerate state groupings off the live `PerturbationTheoryWavefunctions`' corrections, filtering out any singleton (non-degenerate) groups.

        :return: the (filtered) degenerate state groupings, or `None` if none are present
        :rtype: list | None
        """
        data = self.data  # type:PerturbationTheoryWavefunctions
        states = data.corrs.degenerate_states
        if states is not None:
            states = [x for x in data.corrs.degenerate_states if len(x) > 1]
        return states

    @property_dispatcher
    def deperturbed_hamiltonians(self):
        """
        Returns the deperturbed Hamiltonians used to make the degenerate transform

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @deperturbed_hamiltonians.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `deperturbed_hamiltonians`: reads the degenerate-block Hamiltonians out of the checkpoint's stored degenerate-data dict.

        :return: the degenerate-block Hamiltonians
        :rtype: list
        """
        return self.data["degenerate_data"]["hamiltonians"]
    @deperturbed_hamiltonians.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `deperturbed_hamiltonians`: pulls the degenerate-block Hamiltonians off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the degenerate-block Hamiltonians
        :rtype: list
        """
        data = self.data  # type:PerturbationTheoryWavefunctions
        return data.corrs.degenerate_hamiltonians

    @property_dispatcher
    def degenerate_energies(self):
        """
        Returns the deperturbed states used to make the degenerate transform

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @degenerate_energies.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `degenerate_energies`: reads the degenerate-corrected energies out of the checkpoint's stored degenerate-data dict.

        :return: the degenerate-corrected energies
        :rtype: np.ndarray
        """
        return self.data["degenerate_data"]["energies"]
    @degenerate_energies.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `degenerate_energies`: pulls the degenerate-corrected energies off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the degenerate-corrected energies
        :rtype: np.ndarray
        """
        data = self.data  # type:PerturbationTheoryWavefunctions
        return data.corrs.degenerate_energies
    @degenerate_energies.register("log_file")
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `degenerate_energies`: falls back to the parsed log's own `energies` (which already reflect any degenerate treatment applied during the run).

        :return: the degenerate-corrected energies
        :rtype: np.ndarray
        """
        self.data: VPTAnalyzerLogParser
        return self.data.energies

    @property_dispatcher
    def degenerate_rotations(self):
        """
        Returns the deperturbed states used to make the degenerate transform

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @degenerate_rotations.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `degenerate_rotations`: reads the degenerate-perturbation-theory mixing/rotation matrices out of the checkpoint dict.

        :return: the degenerate mixing/rotation matrices
        :rtype: list
        """
        return self.data["degenerate_states"]["rotations"]
    @degenerate_rotations.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `degenerate_rotations`: pulls the degenerate-perturbation-theory mixing/rotation matrices off the live `PerturbationTheoryWavefunctions`' corrections.

        :return: the degenerate mixing/rotation matrices
        :rtype: list
        """
        data = self.data  # type:PerturbationTheoryWavefunctions
        return data.corrs.degenerate_transf

    @property_dispatcher
    def log_file(self):
        """
        Returns the log_file for the run

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    log_file_extension = '.txt'
    @log_file.register("checkpoint")
    def _(self):
        """
        **LLM Docstring**

        Checkpoint-backed implementation of `log_file`: derives the expected log-file path from the checkpoint's own file path (swapping its extension for `log_file_extension`), returning `None` if no such file actually exists on disk.

        :return: the log-file path, or `None` if it doesn't exist
        :rtype: str | None
        """
        tf = self.data.checkpoint_file
        tf, _ = os.path.splitext(tf)
        tf = tf + self.log_file_extension #could be more sophisticated in the future maybe...
        if not os.path.exists(tf):
            tf = None
        return tf
    @log_file.register("wavefunctions")
    def _(self):
        """
        **LLM Docstring**

        Wavefunction-backed implementation of `log_file`: reads the log-file path off the live `PerturbationTheoryWavefunctions`' logger.

        :return: the log-file path
        :rtype: str | None
        """
        return self.data.logger.log_file
    @log_file.register("log_file")
    def _(self):
        """
        **LLM Docstring**

        Log-file-backed implementation of `log_file`: the data itself already is the log file/parser, so it's returned unchanged.

        :return: `self.data`
        :rtype: VPTAnalyzerLogParser
        """
        return self.data

    @property
    def log_parser(self):
        """
        **LLM Docstring**

        A `VPTAnalyzerLogParser` for this loader's data: the data itself, if it already is one, otherwise a freshly built parser wrapping the resolved log file (or `None` if no log file is available).

        :return: the log parser, or `None` if unavailable
        :rtype: VPTAnalyzerLogParser | None
        """
        if isinstance(self.data, VPTAnalyzerLogParser):
            return self.data
        else:
            lf = self.log_file()
            if lf is None: return None
            return VPTAnalyzerLogParser(lf)

def loaded_prop(fn):
    """
    **LLM Docstring**

    Decorator for `VPTAnalyzer` convenience properties that mirror a `VPTResultsLoader` method of the same name: copies that method's docstring onto the decorated function and wraps it as a `property`.

    :param fn: the function to wrap (typically calling the matching `VPTResultsLoader` method on `self.loader`)
    :type fn: callable
    :return: the wrapped property
    :rtype: property
    """
    name = fn.__name__
    fn.__doc__ = getattr(VPTResultsLoader, name).__doc__
    return property(fn)
class VPTAnalyzer:
    """
    Provides analysis tools on VPT results
    """

    def __init__(self, res):
        """

        :param res:
        :type res:
        """
        if not isinstance(res, VPTResultsLoader):
            res = VPTResultsLoader(res)
        self.loader = res

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

        if logger is None:
            logger = Logger(io.StringIO())

        wfns = VPTRunner.run_simple(
            *args,
            logger=logger,
            **kwargs
        )
        return cls(wfns)

    @loaded_prop
    def potential_terms(self):
        """

        :return:
        :rtype:
        """
        return self.loader.potential_terms()
    @loaded_prop
    def kinetic_terms(self):
        """

        :return:
        :rtype:
        """
        return self.loader.kinetic_terms()
    @loaded_prop
    def dipole_terms(self):
        """

        :return:
        :rtype:
        """
        return self.loader.dipole_terms()

    @loaded_prop
    def basis(self):
        """

        :return:
        :rtype:
        """
        return self.loader.basis()
    @loaded_prop
    def target_states(self):
        """

        :return:
        :rtype:
        """
        return self.loader.target_states()

    @loaded_prop
    def spectrum(self):
        """

        :return:
        :rtype:
        """
        return self.loader.spectrum()
    @loaded_prop
    def energy_corrections(self):
        """

        :return:
        :rtype:
        """
        return self.loader.energy_corrections()
    @loaded_prop
    def energies(self):
        """

        :return:
        :rtype:
        """
        engs = self.degenerate_energies
        if engs is None:
            engs = self.loader.energies()
        return engs
    @property
    def frequencies(self):
        """

        :return:
        :rtype:
        """
        return self.energies[1:] - self.energies[0]

    @loaded_prop
    def zero_order_spectrum(self):
        """

        :return:
        :rtype:
        """
        return self.loader.zero_order_spectrum()
    @property
    def deperturbed_spectrum(self):
        """

        :return:
        :rtype:
        """
        freqs = self.deperturbed_frequencies
        # if np.max(freqs) < 1:
        #     freqs = freqs * UnitsData.convert("Hartrees", "Wavenumbers")
        tmom = self.deperturbed_transition_moments[:, 0]
        osc = np.linalg.norm(tmom, axis=0) ** 2
        units = UnitsData.convert("OscillatorStrength", "KilometersPerMole")
        ints = units * freqs * osc[1:]
        return DiscreteSpectrum(freqs * UnitsData.convert("Hartrees", "Wavenumbers"), ints)
    @property
    def deperturbed_frequencies(self):
        """

        :return:
        :rtype:
        """
        return self.deperturbed_energies[1:] - self.deperturbed_energies[0]
    @loaded_prop
    def wavefunction_corrections(self):
        """

        :return:
        :rtype:
        """
        return self.loader.transition_moment_corrections()
    @loaded_prop
    def transition_moment_corrections(self):
        """

        :return:
        :rtype:
        """
        return self.loader.transition_moment_corrections()
    @loaded_prop
    def transition_moments(self):
        """

        :return:
        :rtype:
        """
        return self.loader.transition_moments()
    @loaded_prop
    def deperturbed_transition_moment_corrections(self):
        """

        :return:
        :rtype:
        """
        return self.loader.deperturbed_transition_moment_corrections()
    @loaded_prop
    def deperturbed_transition_moments(self):
        """

        :return:
        :rtype:
        """
        return self.loader.deperturbed_transition_moments()
    @loaded_prop
    def deperturbed_hamiltonians(self):
        """

        :return:
        :rtype:
        """
        return self.loader.deperturbed_hamiltonians()
    @property
    def zero_order_energies(self):
        """

        :return:
        :rtype:
        """
        return self.loader.zero_order_energies()
    @property
    def deperturbed_energies(self):
        """

        :return:
        :rtype:
        """
        return self.loader.deperturbed_energies()
    @loaded_prop
    def degenerate_states(self):
        """

        :return:
        :rtype:
        """
        return self.loader.degenerate_states()
    @loaded_prop
    def degenerate_energies(self):
        """

        :return:
        :rtype:
        """
        return self.loader.degenerate_energies()

    def shift_and_transform_hamiltonian(self, hams, shifts):
        """

        :param hams:
        :type hams:
        :param shifts:
        :type shifts:
        :return:
        :rtype:
        """

        shifts = np.asanyarray(shifts)
        if shifts.ndim == 1:
            shifts = np.diag(shifts)

        ham = sum(hams) + shifts
        deg_engs, deg_transf = np.linalg.eigh(ham)

        # ov_thresh = .5
        # for i in range(len(deg_transf)):
        #     max_ov = np.max(deg_transf[:, i] ** 2)
        #     if max_ov < ov_thresh:  # there must be a single mode that has more than 50% of the initial state character?
        #         logger.log_print(
        #             "    state {i} is more than 50% mixed",
        #             i=i
        #         )
        #     #     raise PerturbationTheoryException("mode {} is has no contribution of greater than {}".format(
        #     #         i, ov_thresh
        #     #     ))

        # print(
        #     *str(
        #         np.round(ham * UnitsData.convert("Hartrees", "Wavenumbers")).astype(int)
        #     ).splitlines(),
        #     sep="\n"
        # )

        # we pick the terms with the max contribution from each input state
        # and zero out the contributions so that two states can't map
        # to the same input state
        sort_transf = np.abs(deg_transf.copy())
        sorting = [-1] * len(deg_transf)
        for i in range(len(deg_transf)):
            o = np.argmax(sort_transf[i, :])
            sorting[i] = o
            sort_transf[:, o] = 0.  # np.zeros(len(sort_transf))
        # print(deg_transf)

        deg_engs = deg_engs[sorting,]
        deg_transf = deg_transf[:, sorting]

        return deg_engs, deg_transf

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

        deg_e, deg_t = self.shift_and_transform_hamiltonian(hams, shifts)
        inds, _ = nput.find(target_states, deg_states)
        subtms = [tm[0][inds] for tm in tmoms]

        tf = deg_t
        if handling_mode == 'transpose':
            tf = deg_t.T

        return deg_e, deg_t, [np.dot(tf, tm) for tm in subtms]
        # full_tm =

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


        eng, deg_transf, deg_tmom = self.get_shifted_transformed_transition_moments(
            deg_states, target_states, hams, shifts, tmoms,
            handling_mode=handling_mode
            )
        osc = np.linalg.norm(deg_tmom, axis=0) ** 2
        units = UnitsData.convert("OscillatorStrength", "KilometersPerMole")
        freqs = (eng - zpe)
        ints = units * freqs * osc#[1:]

        return DiscreteSpectrum(freqs * UnitsData.convert("Hartrees", "Wavenumbers"), ints), deg_transf

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

        spec, tf = self.get_shifted_transformed_spectrum(
            self.energies[0], deg_states, self.target_states, hams, shifts, self.deperturbed_transition_moments,
            handling_mode=handling_mode
            )
        if not return_transformation:
            return spec
        else:
            return spec, tf

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

        states = np.asanyarray(states)

        single = states.ndim == 1
        if single:
            states = np.expand_dims(states, 0)
        inds, _ = nput.find(self.target_states, states)

        freq = self.deperturbed_frequencies if data == 'deperturbed' else self.frequencies

        res = {'index':inds}
        if 'frequency' in keys:
            res['frequency'] = freq[inds-1] * UnitsData.convert("Hartrees", "Wavenumbers")
        if 'max_contribs' in keys:
            raise NotImplementedError("whoops")
            res['max_contribs'] = ...
        if 'transition_moment' in keys or 'intensity' in keys:
            if data == 'deperturbed':
                tm_base = self.deperturbed_transition_moments
            else:
                tm_base = self.transition_moments
            res['transition_moment'] = np.array([tm[0][inds] for tm in tm_base]).T
            res['intensity'] = freq[inds-1] * np.linalg.norm(res['transition_moment'], axis=1)**2 * UnitsData.convert("OscillatorStrength", "KilometersPerMole")

        new_res = []
        for i in range(len(states)):
            new_res.append({k:res[k][i] for k in res.keys()})

        if single:
            new_res = new_res[0]

        return new_res

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

        corrs = []
        order = len(tmom_corrs[0])
        for a in range(3):
            cur = 0
            for o in range(order):
                n = 0
                for i, j, k in ip.product(range(o + 1), range(o + 1), range(o + 1)):
                    if i + j + k == o:
                        if (i, k, j) in terms:
                            cur += tmom_corrs[o][a][n][0][inds]

                        n += 1
            corrs.append(cur)

        res = np.array(corrs).T
        return res

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
        if terms is None:
            terms = {
                'harmonic': [(0, 0, 0)],
                'electrical': [(0, 1, 0), (0, 2, 0)],
                'mechanical': [(1, 0, 0), (0, 0, 1), (1, 0, 1), (2, 0, 0), (0, 0, 2)],
                'mixed': [(1, 1, 0), (0, 1, 1)],
            }

        states = np.asanyarray(states)

        single = states.ndim == 1
        if single:
            states = np.expand_dims(states, 0)
        inds, _ = nput.find(self.target_states, states)

        tmom_corrs = self.deperturbed_transition_moment_corrections if data == 'deperturbed' else self.transition_moment_corrections

        res = {}
        for k,t in terms.items():
            res[k] = self._aggregate_tmoms(tmom_corrs, inds, t)

        new_res = []
        for i in range(len(states)):
            if rotation is None:
                new_res.append({k:res[k][i] for k in res.keys()})
            else:
                new_res.append({k:rotation@res[k][i] for k in res.keys()})

        if single:
            new_res = new_res[0]

        return new_res

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
        return self.transition_moment_term_sums(
            states,
            terms = {
                    'harmonic': [(0, 0, 0)],
                    'electrical': [(0, 1, 0)],
                    'mechanical': [(1, 0, 0), (0, 0, 1)]
                },
            rotation=rotation,
            data=data
        )

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

        woop = self.transition_moment_term_sums(states, terms=terms, data=data)
        single = isinstance(woop, dict)
        if single:
            woop = [woop]

        breakdowns = [{k:np.linalg.norm(v)**2 for k,v in w.items()} for w in woop]

        if single:
            breakdowns = breakdowns[0]

        return breakdowns

    def degenerate_coupling_element(self, state1, state2):
        """

        :param state1:
        :type state1:
        :param state2:
        :type state2:
        :return:
        :rtype:
        """

        for g,h in zip(self.degenerate_states, self.deperturbed_hamiltonians):
            try:
                pos = nput.find(g, np.array([state1]))[0]
            except IndexError:
                continue
            else:
                pos2 = nput.find(g, np.array([state2]))[0]
                return sum(h)[pos[0], pos2[0]] * UnitsData.convert("Hartrees", "Wavenumbers")
        else:
            raise IndexError("not in any degenerate group")


    def format_deperturbed_hamiltonian(self, which):
        """

        :param which:
        :type which:
        :return:
        :rtype:
        """
        with np.printoptions(linewidth=1e9):
            return str(sum(self.deperturbed_hamiltonians[which]) * UnitsData.convert("Hartrees", "Wavenumbers"))

    @property
    def log_parser(self):
        """
        **LLM Docstring**

        The underlying `VPTAnalyzerLogParser` for this analyzer's loaded results, via `self.loader.log_parser`.

        :return: the log parser, or `None` if unavailable
        :rtype: VPTAnalyzerLogParser | None
        """
        return self.loader.log_parser

    def print_output_tables(self,
                            print_energy_corrections=False,
                            print_energies=False,
                            print_transition_moments=False,
                            print_intensities=True,
                            **kwargs):
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
        no_print = False
        try:
            logger = self.loader.data.logger
        except AttributeError:
            no_print = True
            with self.log_parser:
                for block in self.log_parser.get_blocks():
                    if not isinstance(block, str) and block.tag == 'IR Data':
                        print(block)
        else:
            self.loader.data.logger = None
            VPTRunner.print_output_tables(
                self.loader.data,
                print_energy_corrections=print_energy_corrections,
                print_energies=print_energies,
                print_transition_moments=print_transition_moments,
                print_intensities=print_intensities,
                **kwargs
            )
        finally:
            if not no_print:
                self.loader.data.logger = logger