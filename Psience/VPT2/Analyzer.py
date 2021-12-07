"""
Provides analyzer class to handle common VPT analyses
"""

import enum, weakref, functools, numpy as np, itertools as ip
from McUtils.Data import UnitsData
import McUtils.Numputils as nput
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
        raise NotImplementedError("missing")

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
        raise NotImplementedError("missing")

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
        return DiscreteSpectrum(self.data.frequencies(), self.data.intensities)

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

    def energies(self):
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
        Returns the deperturbed Hamiltonians used to make the degenerate transform

        :return:
        :rtype:
        """
        raise ValueError("no dispatch")
    @degenerate_states.register("checkpoint")
    def _(self):
        return self.data["degenerate_states"]
    @degenerate_states.register("wavefunctions")
    def _(self):
        raise NotImplementedError("missing")

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
        return self.data["nondegenerate_hamiltonians"]
    @transition_moment_corrections.register("wavefunctions")
    def _(self):
        raise NotImplementedError("missing")

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
    def basis(self):
        return self.loader.basis()
    @loaded_prop
    def target_states(self):
        return self.loader.target_states()

    @loaded_prop
    def spectrum(self):
        return self.loader.spectrum()
    @loaded_prop
    def energy_corrections(self):
        return self.loader.energy_corrections()
    @loaded_prop
    def energies(self):
        return self.loader.energies()
    @property
    def frequencies(self):
        return self.energies - self.energies[0]
    @loaded_prop
    def wavefunction_corrections(self):
        return self.loader.transition_moment_corrections()
    @loaded_prop
    def transition_moment_corrections(self):
        return self.loader.transition_moment_corrections()
    @loaded_prop
    def transition_moments(self):
        return self.loader.transition_moments()
    @loaded_prop
    def deperturbed_transition_moment_corrections(self):
        return self.loader.deperturbed_transition_moment_corrections()
    @loaded_prop
    def deperturbed_transition_moments(self):
        return self.loader.deperturbed_transition_moments()
    @loaded_prop
    def deperturbed_hamiltonians(self):
        return self.loader.deperturbed_hamiltonians()
    @loaded_prop
    def degenerate_states(self):
        return self.loader.degenerate_states()

    def shift_and_transform_hamiltonian(self, hams, shifts):

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

        # we pick the terms with the max contribution from each input state
        # and zero out the contributions so that two states can't map
        # to the same input state
        sort_transf = np.abs(deg_transf.copy())
        sorting = [-1] * len(deg_transf)
        for i in range(len(deg_transf)):
            o = np.argmax(sort_transf[i, :])
            sorting[i] = o
            sort_transf[:, o] = 0.  # np.zeros(len(sort_transf))

        deg_engs = deg_engs[sorting,]
        deg_transf = deg_transf[:, sorting]

        return deg_engs, deg_transf

    def get_shifted_transformed_transition_moments(self, deg_states, target_states, hams, shifts, tmoms):

        deg_e, deg_t = self.shift_and_transform_hamiltonian(hams, shifts)
        inds, _ = nput.find(target_states, deg_states)
        subtms = [tm[0][inds] for tm in tmoms]

        return deg_e, deg_t, [np.dot(deg_t.T, tm) for tm in subtms]
        # full_tm =

    def get_shifted_transformed_spectrum(self, zpe, deg_states, target_states, hams, shifts, tmoms):

        eng, deg_transf, deg_tmom = self.get_shifted_transformed_transition_moments(deg_states, target_states, hams, shifts, tmoms)
        osc = np.linalg.norm(deg_tmom, axis=0) ** 2
        units = UnitsData.convert("OscillatorStrength", "KilometersPerMole")
        freqs = (eng - zpe)
        ints = units * freqs * osc

        return DiscreteSpectrum(freqs * UnitsData.convert("Hartrees", "Wavenumbers"), ints), deg_transf

    def shifted_transformed_spectrum(self, deg_states, hams, shifts, return_transformation=False):

        spec, tf = self.get_shifted_transformed_spectrum(self.energies[0], deg_states, self.target_states, hams, shifts, self.deperturbed_transition_moments)
        if not return_transformation:
            return spec
        else:
            return spec, tf

    def transition_data(self, states, keys=['frequency', 'transition_moment', 'intensity'], data='deperturbed'):

        states = np.asanyarray(states)

        single = states.ndim == 1
        if single:
            states = np.expand_dims(states, 0)
        inds, _ = nput.find(self.target_states, states)

        res = {'index':inds}
        if 'frequency' in keys:
            res['frequency'] = self.frequencies[inds] * UnitsData.convert("Hartrees", "Wavenumbers")
        if 'max_contribs' in keys:
            raise NotImplementedError("whoops")
            res['max_contribs'] = ...
        if 'transition_moment' in keys or 'intensity' in keys:
            if data == 'deperturbed':
                tm_base = self.deperturbed_transition_moments
            else:
                tm_base = self.transition_moments
            res['transition_moment'] = np.array([tm[0][inds] for tm in tm_base]).T
            res['intensity'] = self.frequencies[inds] * np.linalg.norm(res['transition_moment'], axis=1)**2 * UnitsData.convert("OscillatorStrength", "KilometersPerMole")

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

    def transition_moment_term_sums(self, states, terms=None, data='deperturbed'):
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
            new_res.append({k:res[k][i] for k in res.keys()})

        if single:
            new_res = new_res[0]

        return new_res

    def intensity_breakdown(self, states, terms=None, data='deperturbed'):

        woop = self.transition_moment_term_sums(states, terms=terms, data=data)
        single = isinstance(woop, dict)
        if single:
            woop = [woop]

        breakdowns = [{k:np.linalg.norm(v)**2 for k,v in w.items()} for w in woop]

        if single:
            breakdowns = breakdowns[0]

        return breakdowns

    def degenerate_coupling_element(self, state1, state2):

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
        with np.printoptions(linewidth=1e9):
            return str(sum(self.deperturbed_hamiltonians[which]) * UnitsData.convert("Hartrees", "Wavenumbers"))
