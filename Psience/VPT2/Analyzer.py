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
    registry = {}#weakref.WeakValueDictionary()
    def register(key, registry=registry):
        if not isinstance(key, VPTResultsSource):
            key = VPTResultsSource(key)
        def bind(fn, key=key, registry=registry):
            fn.__name__ = basefn.__name__
            fn.__doc__ = basefn.__doc__
            registry[key] = fn
            return fn
        return bind

    @functools.wraps(basefn)
    def dispatch(self, *args, **kwargs):
        return registry[self.res_type](self, *args, **kwargs)
    dispatch.registry = registry
    dispatch.register = register

    return dispatch

class VPTAnalyzerLogParser(LogParser):
    def __init__(self, log_file, **opts):
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
            self.label_width = None
            super().__init__(spec_str, max_nesting_depth=0, **opts)
        def check_tag(self, line:str, depth: int = 0, active_tag=None, label: str = None, history: list[str] = None):
            ls = line.lstrip()
            if len(ls) == 0 or ls.startswith("State") or line.startswith(" "*5):
                return self.LineReaderTags.SKIP
        def handle_block_line(self, label, line, depth=0, history:list[str]=None):
            if self.label_width is None:
                self.label_width = len(line.split("  ", 1)[0])
            n = self.label_width
            line = line.strip()
            state = line[:n]
            return state, line[n:].split()

    @classmethod
    def reformat_eng_block(cls, sb):
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
        with self.EnergiesBlockParser(spec_str) as parser:
            base_specs = list(parser)

        return self.reformat_eng_block(base_specs)
    @property
    def harmonic_energies(self):
        if self._energies is None:
            try:
                eng_str = self.tree["State Energies"]
            except IndexError:
                eng_str = self.tree[0, "Degenerate Energies"]
            self._energies = self.parse_energies_blocks(eng_str)
        return self._energies[0], self._energies[1]

    @property
    def energies(self):
        if self._energies is None:
            try:
                eng_str = self.tree["State Energies"]
            except IndexError:
                eng_str = self.tree["Degenerate Energies"]
            self._energies = self.parse_energies_blocks(eng_str)
        return self._energies[0], self._energies[2]

    @property
    def deperturbed_energies(self):
        if self._deperturbed_energies is None:
            self._deperturbed_energies = self.parse_energies_blocks(self.tree["Deperturbed Energies"])
        return self._deperturbed_energies[0], self._deperturbed_energies[2]

    @property
    def spectra(self):
        if self._spectra is None:
            self._spectra = self.parse_spectrum_blocks(self.tree["IR Data"])
        return self._spectra

    @property
    def deperturbed_spectra(self):
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
            super().__init__(spec_str, max_nesting_depth=0, **opts)
        def check_tag(self, line:str, depth: int = 0, active_tag=None, label: str = None, history: list[str] = None):
            block_tag = ' Initial State:'
            if len(line) == 0 or line.startswith("State") or line.startswith(" "*5):
                return self.LineReaderTags.SKIP
            elif line.startswith(block_tag):
                return self.LineReaderTags.BLOCK_START, line[len(block_tag):].strip(), None
        def handle_block_line(self, label, line, depth=0, history:list[str]=None):
            n = len(label)
            line = line.strip()
            state = line[:n]
            return state, line[n:].split()

    @classmethod
    def reformat_spec_block(cls, sb):
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
        with self.SpectrumBlockParser(spec_str) as parser:
            base_specs = list(parser)

        return [self.reformat_spec_block(sb) for sb in base_specs]

    class TransitionMomentBlockParser(StringLineByLineReader):
        def __init__(self, spec_str, **opts):
            super().__init__(spec_str, max_nesting_depth=0, **opts)
        def check_tag(self, line:str, depth: int = 0, active_tag=None, label: str = None, history: list[str] = None):
            block_tag = ' Initial State:'
            if len(line) == 0 or line.startswith("State") or line.startswith(" "*5):
                return self.LineReaderTags.SKIP
            elif line.startswith(block_tag):
                return self.LineReaderTags.BLOCK_START, line[len(block_tag):].strip(), None
        def handle_block_line(self, label, line, depth=0, history:list[str]=None):
            n = len(label)
            line = line.strip()
            state = line[:n]
            return state, line[n:].split()

    _indexer = None
    @classmethod
    def load_term_counts(cls, nterms):
        import McUtils.Combinatorics as comb

        if cls._indexer is None:
            cls._indexer = comb.SymmetricGroupGenerator(3)

        cls._indexer.load_to_size(nterms)
        return [c for c in cls._indexer._cumtotals if c < nterms]

    @classmethod
    def reformat_tm_block(cls, sb):
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
        with self.TransitionMomentBlockParser(spec_str) as parser:
            base_specs = list(parser)

        return [self.reformat_tm_block(sb) for sb in base_specs]

    @property
    def transition_moment_corrections(self):
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
        return self.data["potential_terms"]
    @potential_terms.register("wavefunctions")
    def _(self):
        return self.data.corrs.hamiltonian.V_terms
    @potential_terms.register("log_file")
    def _(self):
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
        return {
            "G_matrix": self.data.corrs.hamiltonian.G_terms,
            "coriolis": self.data.corrs.hamiltonian.coriolis_terms,
            "pseudopotential": self.data.corrs.hamiltonian.coriolis_terms
        }
    @kinetic_terms.register("log_file")
    def _(self):
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
        return self.data["dipole_terms"]
    @dipole_terms.register("wavefunctions")
    def _(self):
        return self.data.dipole_terms
    @dipole_terms.register("log_file")
    def _(self):
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
        return self.data["corrections"]["total_states"]
    @basis.register("wavefunctions")
    def _(self):
        data = self.data #type:PerturbationTheoryWavefunctions
        return data.corrs.total_basis
    @basis.register("log_file")
    def _(self):
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
        return self.data["corrections"]["states"]
    @target_states.register("wavefunctions")
    def _(self):
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
        freq, ints = self.data["spectrum"]
        return DiscreteSpectrum(freq * UnitsData.convert("Hartrees", "Wavenumbers"), ints)
    @spectrum.register("wavefunctions")
    def _(self):
        return DiscreteSpectrum(self.data.frequencies() * UnitsData.convert("Hartrees", "Wavenumbers"), self.data.intensities[1:])
    @spectrum.register("log_file")
    def _(self):
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
        freq, ints = self.data["zero_order_spectrum"]
        return DiscreteSpectrum(freq * UnitsData.convert("Hartrees", "Wavenumbers"), ints)
    @zero_order_spectrum.register("wavefunctions")
    def _(self):
        return DiscreteSpectrum(self.data.deperturbed_frequencies(order=0) * UnitsData.convert("Hartrees", "Wavenumbers"),
                                self.data.deperturbed_intensities_to_order(order=0)[1:]
                                )
    @zero_order_spectrum.register("log_file")
    def _(self):
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
        return self.data["corrections"]["energies"]
    @energy_corrections.register("wavefunctions")
    def _(self):
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
        return self.data["corrections"]["wavefunctions"]
    @wavefunction_corrections.register("wavefunctions")
    def _(self):
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
        return self.data["transition_moments"]
    @transition_moment_corrections.register("wavefunctions")
    def _(self):
        return self.data.transition_moment_corrections
    @transition_moment_corrections.register("wavefunctions")
    def _(self):
        return self.data.transition_moment_corrections
    def transition_moments(self):
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
        return self.data["nondegenerate_transition_moments"]
    @deperturbed_transition_moment_corrections.register("wavefunctions")
    def _(self):
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
        return self.data["degenerate_data"]["states"]
    @degenerate_states.register("wavefunctions")
    def _(self):
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
        return self.data["degenerate_data"]["hamiltonians"]
    @deperturbed_hamiltonians.register("wavefunctions")
    def _(self):
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
        return self.data["degenerate_data"]["energies"]
    @degenerate_energies.register("wavefunctions")
    def _(self):
        data = self.data  # type:PerturbationTheoryWavefunctions
        return data.corrs.degenerate_energies
    @degenerate_energies.register("log_file")
    def _(self):
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
        return self.data["degenerate_states"]["rotations"]
    @degenerate_rotations.register("wavefunctions")
    def _(self):
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
        tf = self.data.checkpoint_file
        tf, _ = os.path.splitext(tf)
        tf = tf + self.log_file_extension #could be more sophisticated in the future maybe...
        if not os.path.exists(tf):
            tf = None
        return tf
    @log_file.register("wavefunctions")
    def _(self):
        return self.data.logger.log_file
    @log_file.register("log_file")
    def _(self):
        return self.data

    @property
    def log_parser(self):
        if isinstance(self.data, VPTAnalyzerLogParser):
            return self.data
        else:
            lf = self.log_file()
            if lf is None: return None
            return VPTAnalyzerLogParser(lf)

def loaded_prop(fn):
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
    def deperturbed_energies(self):
        """

        :return:
        :rtype:
        """
        return self.loader.energies()
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
        return self.loader.log_parser

    def print_output_tables(self,
                            print_energy_corrections=False,
                            print_energies=False,
                            print_transition_moments=False,
                            print_intensities=True,
                            **kwargs):
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