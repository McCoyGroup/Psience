"""
Provides analyzer class to handle common VPT analyses
"""

import enum, weakref, functools
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Checkpointer
from ..Spectra import DiscreteSpectrum
from .Wavefunctions import PerturbationTheoryWavefunctions

__all__ = [
    "VPTResultsLoader",
    "VPTResultsSource",
    "VPTAnalyzer"
]

class VPTResultsSource(enum.Enum):
    """
    Enum of sources to load PT results from
    """
    Wavefunctions = "wavefunctions"
    Checkpoint = "checkpoint"

def property_dispatcher(basefn):
    registry = weakref.WeakKeyDictionary()
    def register(key):
        if not isinstance(key, VPTResultsSource):
            key = VPTResultsSource(key)
        def bind(fn, key=key):
            fn.__name__ = basefn.__name__
            fn.__doc__ = basefn.__doc__
            registry[key] = fn
            return fn
        return bind

    @functools.wraps(basefn)
    def dispatch(self, *args, **kwargs):
        return registry[self.res_type](self, *args, **kwargs)
    dispatch.register = register

    return dispatch

class VPTResultsLoader:
    """
    Provides tools for loading results into canonical
    sources from a simulation, both from checkpoint files and from
    `PerturbationTheoryWavefunctions` objects and potentially more
    """

    def __init__(self, res, res_type=None):
        if isinstance(res, str):
            res = Checkpointer.from_file(res)
        self.data = res
        self.res_type = self.get_res_type(res) if res_type is None else res_type

    def get_res_type(self, res):
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
        return {
            "G_matrix": self.data["gmatrix_terms"],
            "coriolis": self.data["coriolis_terms"],
            "pseudopotential": self.data["pseudopotential_terms"]
        }
    @kinetic_terms.register("wavefunctions")
    def _(self):
        return {
            "G_matrix": self.data.corrs.hamiltonian.G_terms,
            "coriolis": self.data.corrs.hamiltonian.coriolis_terms,
            "pseudopotential": self.data.corrs.hamiltonian.coriolis_terms
        }

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
        return DiscreteSpectrum(*self.data["dipole_terms"])
    @dipole_terms.register("wavefunctions")
    def _(self):
        return self.data.dipole_terms

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
        return DiscreteSpectrum(*self.data["spectrum"])
    @spectrum.register("wavefunctions")
    def _(self):
        return DiscreteSpectrum(self.data.frequencies() * UnitsData.convert("Hartrees", "Wavenumbers"), self.data.intensities)

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
        return self.data.corrs.energies

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

def loaded_prop(fn):
    name = fn.__name__
    fn.__doc__ = getattr(VPTResultsLoader, name).__doc__
    return property(fn)
class VPTAnalyzer:
    """
    Provides analysis tools on VPT results
    """

    def __init__(self, res):
        if not isinstance(res, VPTResultsLoader):
            res = VPTResultsLoader(res)
        self.loader = res

    @loaded_prop
    def potential_terms(self):
        return self.loader.potential_terms()
    @loaded_prop
    def kinetic_terms(self):
        return self.loader.kinetic_terms()
    @loaded_prop
    def dipole_terms(self):
        return self.loader.dipole_terms()

    @loaded_prop
    def spectrum(self):
        return self.loader.spectrum()
    @loaded_prop
    def energy_corrections(self):
        return self.loader.energy_corrections()
    @loaded_prop
    def wavefunction_corrections(self):
        return self.loader.transition_moment_corrections()
    @loaded_prop
    def transition_moment_corrections(self):
        return self.loader.transition_moment_corrections()
