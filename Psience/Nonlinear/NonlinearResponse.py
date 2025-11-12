import collections
import itertools
import numpy as np
import scipy.fft
import enum

from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from .. import BasisReps as breps

from McUtils.Data import UnitsData
import McUtils.Numputils as nput
import McUtils.Combinatorics as comb

__all__ = [
    "liouville_pathways",
    "nonlinear_response_generators",
    "experimental_response_generator"
]

def nested_commutator_expansion(k, side='left'):
    if side == 'left':
        comm = 0
        for i in range(1, k):
            comm = [comm, i]
    else:
        comm = k - 1
        for i in range(k - 2, -1, -1):
            comm = [i, comm]
    return comm

def closed_system_propagator(interaction_generator, time_delays):
    iteractions = [interaction_generator(t) for t in time_delays]
    interaction_n = len(iteractions)
    comms = nput.commutator_terms(nested_commutator_expansion(interaction_n))
    #TODO: allow for pathway filtering on the commutator terms
    return nput.commutator_evaluate(
        comms, iteractions,
        direct=False,
        normalized=True
    )

def liouville_pathways(k):
    final_states = [[i, k-i] for i in range(k+1)]
    subperms = [ # just lattice paths, but people don't call them that
        comb.UniquePermutations([0]*l + [1]*r).permutations()
        for l,r in final_states
    ]
    return np.concatenate(subperms)

def four_wave_averaging_function(polarization_vectors):
    rows, cols = np.triu_indices(4, k=1)
    polarization_cosines = nput.vec_dots(polarization_vectors[rows,], polarization_vectors[cols])

    # can easily be generalized
    term_orderings = [
        # in the basis (ab, ac, ad, bc, bd, cd)
        [0, 5, 1, 4, 2, 3],
        [1, 4, 0, 5, 2, 3],
        [2, 3, 1, 4, 0, 5]
    ]

    invariant_terms = {
        (i, j):
            (4 * polarization_cosines[i] * polarization_cosines[j]
             - polarization_cosines[k] * polarization_cosines[l]
             - polarization_cosines[m] * polarization_cosines[n]
             )

        for i, j, k, l, m, n
        in term_orderings
    }

    def average(unit_dipoles):
        unit_dipoles = np.asanyarray(unit_dipoles)
        dipole_cosines = nput.vec_dots(unit_dipoles[rows,], unit_dipoles[cols,])
        terms = [
            dipole_cosines[i]*dipole_cosines[j] * v
            for (i,j),v in invariant_terms.items()
        ]
        return (1 / 30) * sum(terms)

    return average

def interpret_polarization(pol):
    eye = np.eye(3)
    if isinstance(pol, str):
        _ = []
        for p in pol.lower():
            if p == 'x':
                _.append(eye[0])
            elif p == 'y':
                _.append(eye[1])
            elif p == 'z':
                _.append(eye[2])
        pol = _
    return np.asanyarray(pol)

TransitionData = collections.namedtuple("TransitionData",
                                        ['states', 'frequencies', 'transition_moments', 'couplings'])
def prep_nonlinear_transition_data(transition_dict: dict,
                                   states=None,
                                   couplings=None,
                                   frequencies=None,
                                   transition_moments=None,
                                   band_coherences=None,
                                   realign_transition_moments=False
                                   ):
    if states is None:
        if couplings is not None:
            raise ValueError(f"can't supply ordered `couplings` matrix without `states`")
        if frequencies is not None:
            raise ValueError(f"can't supply ordered `frequencies` matrix without `states`")
        if transition_moments is not None:
            raise ValueError(f"can't supply ordered `transition_moments` matrix without `states`")

        states = []
        for (si, sj) in transition_dict.keys():
            states.append(si)
            states.append(sj)

        if nput.is_int(states[0]):
            states = np.sort(states)
        else:
            # order the states lexicographically and by number of quanta
            states = BasisStateSpace(
                HarmonicOscillatorProductBasis(len(states[0])),
                states
            ).as_sorted()

    if not hasattr(states, 'as_excitations'): # check that we were given a proper basis object
        int_states = nput.is_int(states[0])
        if not int_states:
            states = BasisStateSpace(
                HarmonicOscillatorProductBasis(len(states[0])),
                states
            ).as_sorted()
    else:
        int_states = False

    if hasattr(states, 'as_excitations'):
        state_vector_tuples = [tuple(s) for s in states.excitations.tolist()]
        state_index_map = {s:i for i,s in enumerate(state_vector_tuples)}
    else:
        state_vector_tuples = []
        state_index_map = {s: i for i, s in enumerate(states)}

    nstates = len(states)
    needs_couplings = couplings is None
    if needs_couplings:
        couplings = np.zeros((nstates, nstates), dtype=float)
    needs_frequencies = frequencies is None
    if needs_frequencies:
        frequencies = np.zeros((nstates, nstates), dtype=float)
    needs_tms = transition_moments is None
    if needs_tms:
        transition_moments = np.zeros((nstates, nstates, 3), dtype=float)

    for (si, sj), tdata in transition_dict.items():
        # canonicalize index and state vector
        if int_states:
            i = state_index_map[si]
            j = state_index_map[sj]
        else:
            if nput.is_int(si):
                i = si
                si = state_vector_tuples[si]
            else:
                i = state_index_map[si]

            if nput.is_int(sj):
                j = sj
                sj = state_vector_tuples[sj]
            else:
                j = state_index_map[sj]

        if needs_couplings:
            coupling = tdata.get('coupling')
            if coupling is None and band_coherences is not None:
                if callable(band_coherences):
                    coupling = band_coherences(si, sj)
                elif int_states:
                    coupling = 0.0
                else:
                    nquanta_1 = sum(si)
                    nquanta_2 = sum(sj)
                    coupling = band_coherences.get((nquanta_1, nquanta_2), 0)
            couplings[i, j] = coupling
            if abs(couplings[j, i]) < 1e-12:
                couplings[j, i] = couplings[i, j]
        if needs_tms:
            tm = tdata.get('transition_moment', np.zeros(3))
            if not int_states and realign_transition_moments:
                # try to align 1-> 2 transition moments with the 0 -> 1
                s_diff = np.array(sj) - np.array(si)
                if np.all(s_diff > 0):
                    s0 = tuple([0] * len(s_diff))
                    sk = tuple(s_diff.astype(int))
                    x = state_index_map.get(s0)
                    y = state_index_map.get(sk)
                    if x is not None and y is not None:
                        tm_base = transition_moments[x, y]
                        if np.linalg.norm(tm_base) < 1e-8: # might not have been initialized
                            tm_base = transition_dict.get((s0, sk), {}).get('transition_moment', np.zeros(3))
                        cos_sign = np.sign(np.dot(tm, tm_base))
                        if cos_sign == 0: cos_sign = 1
                        tm = tm * cos_sign
            transition_moments[i, j] = tm
            if np.linalg.norm(transition_moments[j, i]) < 1e-12:
                transition_moments[j, i] = transition_moments[i, j]
        if needs_frequencies:
            frequencies[i, j] = tdata.get('frequency', 0)
            if abs(frequencies[j, i]) < 1e-12:
                frequencies[j, i] = -frequencies[i, j]

    if needs_couplings:
        for si in state_vector_tuples:
            if int_states:
                i = state_index_map[si]
            else:
                i = state_index_map[si]

            if abs(couplings[i, i]) < 1e-8 and band_coherences is not None:
                if callable(band_coherences):
                    coupling = band_coherences(si, si)
                elif int_states:
                    coupling = 0.0
                else:
                    nquanta = sum(si)
                    coupling = band_coherences.get((nquanta, nquanta),
                                                   band_coherences.get(nquanta, 0))
                couplings[i, i] = coupling

    return TransitionData(states, frequencies, transition_moments, couplings)

def get_interaction_basis(initial_states:BasisStateSpace, *, selection_rules, **filter_opts):
    def _apply_rules(space, rules):
        if hasattr(space, 'representative_space'):
            space = space.to_single(include_representative=False).take_unique()
        return space.apply_selection_rules(rules, **filter_opts)

    if nput.is_int(selection_rules[0][0][0]): # one path supplied
        bases = [initial_states]
        space = initial_states
        for rules in selection_rules:
            new = _apply_rules(bases[-1], rules)
            bases.append(new)
            space = space.union(new.to_single())
        total_space = space.to_single().take_unique().as_sorted()
    else:
        raise ValueError('disabled code path')
        bases = []
        spaces = []
        for rule_list in selection_rules:
            basis = [initial_states]
            space = initial_states
            for rules in rule_list:
                new = _apply_rules(basis[-1], rules)
                basis.append(new)
                space = space.union(new.to_single())
            bases.append(basis)
            spaces.append(space.to_single().take_unique())
        total_space = spaces[0]
        for s in spaces:
            total_space = total_space.union(s)
        total_space = total_space.take_unique().as_sorted()
    return total_space, bases

def get_basis_indices(total_basis, sel_space):
    if hasattr(sel_space, 'representative_space'):
        rep_inds = total_basis.find(sel_space.representative_space)
        trans_inds = [total_basis.find(s) for s in sel_space.spaces]
        return {
            i:s
            for i,s in zip(rep_inds, trans_inds)
        }
    else:
        return total_basis.find(sel_space)

def _enumerate_path_indices(paths):
    npath, ninter = paths.shape
    path_blocks = np.zeros((npath, ninter+1, 2), dtype=int)
    x = np.arange(npath)
    for i,p in enumerate(paths.T):
        path_blocks[:, i+1] = path_blocks[:, i]
        path_blocks[x, i+1, p] += 1
    return path_blocks

def prep_liouville_spaces(initial_states, paths, num_interactions=None, selection_rules=((1,), (-1,)), **filter_opts):
    smol = not isinstance(paths, dict)
    if smol:
        paths = {None:paths}
    paths = {k:np.asanyarray(p) for k,p in paths.items()}
    smol_rules = not isinstance(selection_rules, dict)
    if smol_rules:
        selection_rules = {k:selection_rules for k in paths.keys()}
    all_bases = {}
    path_map = {}
    if num_interactions is None:
        num_interactions = len(next(iter(paths.values()))[0]) + 1
    for k,rules in selection_rules.items():
        # a small efficiency, we only enumerate states for unique transition paths
        # and then use indexing to duplicate them
        kpaths = paths[k]
        if nput.is_int(rules[0][0]): # same rules for everything
            num_ints = len(kpaths[0])
            rules = [
                [rules] * num_ints
                for _ in kpaths
            ]
        elif nput.is_int(rules[0][0][0]):
            rules = [
                rules
                for _ in kpaths
            ]

        rule_groups = []
        for p,rg in zip(kpaths, rules):
            lr_rules = [], []
            for i,r in zip(p,rg):
                lr_rules[i].append(tuple(tuple(x) for x in sorted(r)))
            # final interaction can happen on right or left
            if len(rg) < num_interactions:
                lr_rules[0].append(((-1,), (1,)))
                lr_rules[1].append(((-1,), (1,)))
            elif len(p) < num_interactions:
                term = tuple(tuple(x) for x in sorted(rg[-1]))
                lr_rules[0].append(term)
                lr_rules[1].append(term)
            rule_groups.append(tuple(tuple(r) for r in lr_rules))


        flat_groups = [x for rg in rule_groups for x in rg]
        for rules_tuple in flat_groups:
            if rules_tuple not in all_bases:
                if len(rules_tuple) > 0:
                    all_bases[rules_tuple] = get_interaction_basis(initial_states, selection_rules=rules_tuple, **filter_opts)
                else:
                    all_bases[()] = (initial_states, [initial_states])
        for p,g in zip(kpaths, rule_groups):
            path_map[(k, tuple(p))] = g

    total_space = None
    for k,(space,_) in all_bases.items():
        if total_space is None:
            total_space = space
        else:
            total_space = total_space.union(space)
    total_space = total_space.take_unique().as_sorted()

    all_inds = {}
    for k, (_, interaction_basis) in all_bases.items():
        if hasattr(interaction_basis[0], 'as_excitations'): # one path
            all_inds[k] = [get_basis_indices(total_space, b) for b in interaction_basis]
        else:
            all_inds[k] = [
                [get_basis_indices(total_space, b) for b in subbasis]
                for subbasis in interaction_basis
            ]

    path_paths = {}
    for k,paths in paths.items():
        path_sampling = _enumerate_path_indices(paths)

        subinds = []
        for path, old in zip(path_sampling, paths):
            sign = (-1)**path[-1, -1]
            groups = [all_inds[x] for x in path_map[(k, tuple(old))]]
            inds = [[g[i] for g, i in zip(groups, p)] for p in path]
            if len(inds) < num_interactions + 1:
                inds.append([groups[0][path[-1][0]+1], groups[1][path[-1][1]+1]])

            subinds.append([sign, old, inds])
        path_paths[k] = subinds
    if smol:
        path_paths = path_paths[None]
    return total_space, path_paths

def enumerate_state_paths(liou_path, index_set, num_interactions=None):
    queue = collections.deque([])
    left, right = index_set[0]
    for l,r in zip(left,right):
        queue.append([[(l, r),], [], 1])
    termini = [0, 0]
    if num_interactions is None:
        m = len(liou_path)
    else:
        m = num_interactions - 1
    while queue:
        cur_path, transitions, n = queue.pop()
        done = n == m + 1
        if not done:
            i = liou_path[n-1]
            termini[i] += 1 # should reuse this from prior calcs but so cheap
            nxt = index_set[n][i]

            rem = cur_path[-1][(i+1)%2]
            last = cur_path[-1][i]
            for new in nxt[last]:
                new_path = cur_path + (
                    [(new, rem)]
                        if i == 0 else
                    [(rem, new)]
                )
                new_trans = transitions + [(last, new)]
                queue.append([new_path, new_trans, n+1])
        else:
            left, right = cur_path[-1]
            left_nxt, right_nxt = index_set[n]
            if len(liou_path) < n:
                use_left = termini[0] > termini[1]
            else:
                use_left = liou_path[n-1] == 0
            if use_left:
                for l in left_nxt[left]:
                    new_path = cur_path + [(l, right)]
                    new_trans = transitions + [(left, l)]
                    yield new_path, new_trans
            else:
                for r in right_nxt[right]:
                    new_path = cur_path + [(left, r)]
                    new_trans = transitions + [(right, r)]
                    yield new_path, new_trans

def expand_transition_data(td:TransitionData, full_space:BasisStateSpace):
    cur_pos = full_space.find(td.states, missing_val=-1)
    good_pos = np.where(cur_pos > -1)[0]
    if len(good_pos) == 0:
        raise ValueError("no transitions are in the full basis accessed")
    take_sel = np.ix_(good_pos, good_pos)
    cur_pos = cur_pos[good_pos,]
    set_sel = np.ix_(cur_pos, cur_pos)

    nstates = len(full_space)
    frequencies = np.zeros((nstates, nstates))
    couplings = np.zeros((nstates, nstates))
    transition_moments = np.zeros((nstates, nstates, 3))
    frequencies[set_sel] = td.frequencies[take_sel]
    couplings[set_sel] = td.couplings[take_sel]
    transition_moments[set_sel] = td.transition_moments[take_sel]

    return TransitionData(full_space, frequencies, transition_moments, couplings)

def _identify_response_tensor_paths(path_trees, total_space=None,
                                    response_tensor_elements=None,
                                    num_interactions=None):
    if response_tensor_elements is None:
        response_tensor_elements = [(i, i) for i in range(len(total_space))]
    elif not nput.is_int(response_tensor_elements[0]):
        response_inds_left = total_space.find(
            [l for l,r in response_tensor_elements]
        )
        response_inds_right = total_space.find(
            [r for l,r in response_tensor_elements]
        )
        response_tensor_elements = list(zip(response_inds_left, response_inds_right))

    response_tensor_elements = {
        e:[]
        for e in response_tensor_elements
    }

    for sign, p, inds in path_trees:
        for spath, trans in enumerate_state_paths(p, inds, num_interactions=num_interactions):
            if spath[-1] in response_tensor_elements:
                response_tensor_elements[spath[-1]].append((trans, sign))

    return response_tensor_elements

class NonlinearReponseApplicationDomain(enum.Enum):
    Time = "time"
    Frequency = "frequency"

def nonlinear_path_response_function(frequencies, dipole_magnitudes, dipole_directions, couplings, averaging, path,
                                    zero_dipole_cutoff=1e-8, application_domain="time"):
    mag = np.prod([
        dipole_magnitudes[i, j]
        for i,j in path
    ])
    if abs(zero_dipole_cutoff) < 1e-8:
        return None
    avg = averaging(
        dipole_directions[tuple(p[0] for p in path), tuple(p[1] for p in path)]
    )
    prefactor = mag * avg
    if abs(prefactor) < 1e-8:
        return None
    freqs = [
        frequencies[i, j]
        for i,j in path[:-1]
    ]
    cups = [
        couplings[i, j]
        for i, j in path[:-1]
    ]
    npulse = len(freqs)

    application_domain = NonlinearReponseApplicationDomain(application_domain)
    if application_domain == NonlinearReponseApplicationDomain.Time:
        def response(*times):
            if len(times) < npulse:
                raise ValueError(f"expected {npulse} interaction times, got {len(times)}")
            return prefactor * np.prod([
                np.exp((f * 1j - c) * t)
                for f, c, t in zip(freqs, cups, times)
            ], axis=0)
    elif application_domain == NonlinearReponseApplicationDomain.Frequency:
        def response(*fs):
            if len(fs) < npulse:
                raise ValueError(f"expected {npulse} interaction frequencies, got {len(fs)}")
            # cum_freqs = np.cumsum(fs, axis=0)
            return prefactor * np.prod([
                (w - f - c*1j)/((w - f)**2 + c**2 + 1e-14) # just make it huge if we're at a pole
                for f, c, w in zip(freqs, cups, fs)
            ], axis=0)
    else:
        raise ValueError(f"unknown application domain {application_domain}")
    return response

def nonlinear_multipath_response_function(frequencies, dipole_magnitudes, dipole_directions, couplings, averaging,
                                          paths_and_signs, zero_dipole_cutoff=1e-8, application_domain="time"):
    signs = []
    response_functions = []
    for path, sign in paths_and_signs:
        response = nonlinear_path_response_function(
            frequencies, dipole_magnitudes, dipole_directions, couplings, averaging, path,
            zero_dipole_cutoff=zero_dipole_cutoff,
            application_domain=application_domain
        )
        if response is not None:
            response_functions.append(response)
            signs.append(sign)

    if len(response_functions) == 0: return None
    def response(*times):
        response_values = None
        for s,r in zip(signs, response_functions):
            contrib = s * r(*times)
            if response_values is None:
                response_values = contrib
            else:
                response_values += contrib
        return response_values

    return response

def nonlinear_multipath_response_function_generators(tensor_element_paths, application_domain="time", aggregate=True):
    def response_function_generator(frequencies, dipole_magnitudes, dipole_directions, couplings, averaging,
                                    zero_dipole_cutoff=1e-8):
        response_tensor_generators = {}
        for element, path_data in tensor_element_paths.items():
            generator = nonlinear_multipath_response_function(
                frequencies, dipole_magnitudes, dipole_directions, couplings, averaging,
                path_data, zero_dipole_cutoff=zero_dipole_cutoff,
                application_domain=application_domain
            )
            if generator is not None:
                response_tensor_generators[element] = generator
        if aggregate:
            def total_response(*times):
                return sum(
                    f(*times)
                    for f in response_tensor_generators.values()
                )
            return total_response
        else:
            return response_tensor_generators
    return response_function_generator

class NonlinearResponseFunction:
    def __init__(self, response_generator, central_frequency=0,
                 frequency_unit='Hartrees',
                 time_unit='PicoSeconds',
                 application_domain: NonlinearReponseApplicationDomain = "time"):
        self.caller = response_generator
        self.center = central_frequency
        self._fu = frequency_unit
        self._tu = time_unit
        self._dom = NonlinearReponseApplicationDomain(application_domain)
        self._conv = None

    def call_times(self, t1, t2, t3):
        if self._dom == NonlinearReponseApplicationDomain.Time:
            t1_steps, t1_dt = t1
            t3_steps, t3_dt = t3
            t1_times = np.arange(t1_steps) * t1_dt
            t3_times = np.arange(t3_steps) * t3_dt
            t1, t3 = np.meshgrid(t1_times, t3_times)
            t2 = np.full(t1.shape, t2)
            if isinstance(self.caller, dict):
                signals = {k:f(t1, t2, t3) for k,f in self.caller.items()}
            else:
                signals = self.caller(t1, t2, t3)
            return signals
        else:
            raise NotImplementedError("calling with times with a frequency domain response function not currently supprted")

    @classmethod
    def _freq_conv(cls, freq_unit, time_unit):
        return UnitsData.convert(freq_unit, "Hertz") * UnitsData.convert(time_unit, 'Seconds')
    @property
    def conv(self):
        if self._conv is None:
            self._conv = self._freq_conv(self._fu, self._tu)
        return self._conv

    max_samples = 1024
    @classmethod
    def _determine_fft_samples(cls,
                               min_freq, max_freq,
                               num_periods,
                               sampling_fidelity,
                               time_step=None,
                               num_samples=None):
        if time_step is None:
            if num_samples is None:
                time_step = 1 / (sampling_fidelity*max_freq)
                num_samples = int(np.ceil((num_periods/min_freq) / time_step))
                if num_samples > cls.max_samples:
                    num_samples, time_step = cls._determine_fft_samples(
                        min_freq, max_freq,
                        num_periods,
                        sampling_fidelity,
                        time_step=None,
                        num_samples=cls.max_samples
                    )
            else:
                time_step = (num_periods/min_freq) / num_samples
        elif num_samples is None:
            num_samples = min([
                int(np.ceil((num_periods/min_freq) / time_step)),
                cls.max_samples
                ])
        return num_samples, time_step

    max_min_freq_wn = 1
    def prep_sampling_ranges(self,
                             w1_freqs, t2, w3_freqs,
                             num_periods=10,
                             sampling_fidelity=5,
                             time_step=None,
                             num_samples=None
                             ):
        w1_freqs = w1_freqs - self.center
        w3_freqs = w3_freqs - self.center

        if time_step is None or nput.is_numeric(time_step):
            time_step = [time_step, time_step]
        t1, t3 = time_step
        if num_samples is None or nput.is_numeric(num_samples):
            num_samples = [num_samples, num_samples]
        n1, n3 = num_samples

        min_freq1 = max([
            np.min(np.abs(w1_freqs)),
            self.max_min_freq_wn*UnitsData.convert("Wavenumbers", self._fu)
        ])
        max_freq1 = np.max(np.abs(w1_freqs))
        n1, t1 = self._determine_fft_samples(min_freq1 * self.conv, max_freq1* self.conv,
                                             num_periods,
                                             sampling_fidelity,
                                             time_step=t1,
                                             num_samples=n1
                                             )

        min_freq3 = max([
            np.min(np.abs(w3_freqs)),
            self.max_min_freq_wn * UnitsData.convert("Wavenumbers", self._fu)
        ])
        max_freq3 = np.max(np.abs(w3_freqs))
        n3, t3 = self._determine_fft_samples(min_freq3 * self.conv, max_freq3 * self.conv,
                                             num_periods,
                                             sampling_fidelity,
                                             time_step=t3,
                                             num_samples=n3
                                             )

        return (w1_freqs, w3_freqs), ((n1, t1), t2, (n3, t3))

    def call_freqs(self, w1, t2, w3,
                   *,
                   apply_ifft=True,
                   shift_ifft=True,
                   return_freqs=True,
                   default_frequency_divisions=100,
                   **sampling_options):
        if self._dom == NonlinearReponseApplicationDomain.Time:
            w1_freq_min, w1_freq_max = w1
            w3_freq_min, w3_freq_max = w3
            (w1_freqs, w3_freqs), (t1, t2, t3) = self.prep_sampling_ranges(
                [w1_freq_min, w1_freq_max],
                t2,
                [w3_freq_min, w3_freq_max],
                **sampling_options
            )
            signal = self.call_times(t1, t2, t3)
            if apply_ifft:
                return self._apply_ifft(t1, t3, signal, shift=shift_ifft, return_freqs=return_freqs)
            else:
                return (w1_freqs, w3_freqs), (t1, t3), signal
        elif self._dom == NonlinearReponseApplicationDomain.Frequency:
            w2 = t2 / self._freq_conv(self._fu, self._tu)
            if len(w1) == 2:
                w1_min, w1_max = w1
                w1_steps = default_frequency_divisions
            else:
                w1_min, w1_max, w1_steps = w1
            if len(w3) == 2:
                w3_min, w3_max = w3
                w3_steps = default_frequency_divisions
            else:
                w3_min, w3_max, w3_steps = w3
            w1_freqs = np.linspace(w1_min-self.center, w1_max-self.center, w1_steps)
            w3_freqs = np.linspace(w3_min-self.center, w3_max-self.center, w3_steps)
            w1, w3 = np.meshgrid(w1_freqs, w3_freqs)
            w2 = np.full(w1.shape, w2)
            if isinstance(self.caller, dict):
                signals = {k:f(w1, w2, w3) for k,f in self.caller.items()}
            else:
                signals = self.caller(w1, w2, w3)
            if return_freqs:
                return (w1_freqs+self.center, w3_freqs+self.center), signals
            else:
                return signals
        else:
            raise ValueError(f"unknown application domain {self._dom}")

    def _get_fft_freqs(self, n, dt):
        freqs = scipy.fft.fftfreq(n, dt)
        freqs = scipy.fft.fftshift(freqs)
        conv = 1 / self._freq_conv(self._fu, self._tu) # maybe should do this directly for stability...
        return conv * freqs + self.center

    def _apply_ifft(self,
                    t1, t3, signal,
                    shift=True, return_freqs=True
                    ):
        if isinstance(signal, dict):
            keys = list(signal.keys())
            S_fs = {
                k: self._apply_ifft(
                    t1, t3, signal[k], shift=shift,
                    return_freqs=return_freqs if i == 0 else False
                )
                for i, k in enumerate(keys)
            }
            if return_freqs:
                freqs, S_f = S_fs[keys[0]]
                S_fs[keys[0]] = S_f
                return freqs, S_fs
            else:
                return S_fs
        else:
            t1_steps, t1_dt = t1
            t3_steps, t3_dt = t3
            S_f = scipy.fft.ifft2(signal, s=(t1_steps, t3_steps))
            if shift:
                S_f = scipy.fft.fftshift(S_f)

            if return_freqs:
                freqs1 = self._get_fft_freqs(t1_steps, t1_dt)
                freqs3 = self._get_fft_freqs(t3_steps, t3_dt)
                return (freqs1, freqs3), S_f
            else:
                return S_f

def nonlinear_response_generators(transition_data,
                                  paths,
                                  num_interactions=None,
                                  selection_rules=((1,),(-1,)),
                                  polarization=None,
                                  initial_states=None,
                                  state_filter_opts=None,
                                  response_tensor_elements=None,
                                  central_frequency=True,
                                  frequency_unit="Hartrees",
                                  time_unit="PicoSeconds",
                                  application_domain="time",
                                  response_function_class=None,
                                  **state_opts
                                  ):
    td = prep_nonlinear_transition_data(transition_data, **state_opts)
    if central_frequency is True:
        freqs = td.frequencies
        sel = np.abs(freqs) > 1e-8
        central_frequency = np.average(np.abs(freqs[sel]))

    if central_frequency is False or central_frequency is None:
        central_frequency = 0
    else:
        freqs = td.frequencies.copy()
        sel = np.abs(freqs) > 1e-8
        vals = freqs[sel]
        freqs[sel] = np.sign(vals) * (np.abs(vals) - central_frequency)
        td = td._replace(frequencies=freqs)

    states = td.states
    if not hasattr(states, 'as_excitations'): # path enumeration requires
        raise ValueError("state vectors are required")

    if polarization is not None:
        polarization = interpret_polarization(polarization) # TODO: ensure norm
        orientational_averaging = four_wave_averaging_function(polarization)
    else:
        orientational_averaging = lambda d: 1

    needs_paths = (
        # TODO: handle the most generic case too?
        nput.is_int(next(iter(paths.values()))[0][0])
            if isinstance(paths, dict) else
        nput.is_int(paths[0][0])
    )
    if needs_paths: # can also just supply bases directly
        if not hasattr(initial_states, 'as_excitations'):
            if initial_states is None:
                initial_states = 0
            state_mode = None
            smol_init = nput.is_int(initial_states)
            if smol_init:  # TODO: handle state basis element?
                initial_states = [initial_states]
                state_mode = 'indices'

            initial_states = BasisStateSpace(
                breps.HarmonicOscillatorProductBasis(states.ndim),
                initial_states,
                mode=state_mode
            )

    if state_filter_opts is None:
        state_filter_opts = {}
    total_space, subblocks = prep_liouville_spaces(initial_states, paths,
                                                   selection_rules=selection_rules,
                                                   num_interactions=num_interactions,
                                                   **state_filter_opts)
    _, frequencies, transition_moments, couplings = expand_transition_data(td, total_space)
    dipole_directions, dipole_magnitudes = nput.vec_normalize(transition_moments, return_norms=True)

    application_domain = NonlinearReponseApplicationDomain(application_domain)
    if isinstance(subblocks, dict):
        nlf = {}
        for k,sububs in subblocks.items():
            response_tensor_paths = _identify_response_tensor_paths(
                sububs,
                total_space=total_space,
                response_tensor_elements=response_tensor_elements,
                num_interactions=num_interactions
            )

            nlf[k] = nonlinear_multipath_response_function_generators(response_tensor_paths, application_domain=application_domain)(
                frequencies,  dipole_magnitudes, dipole_directions, couplings, orientational_averaging,
            )
    else:
        response_tensor_paths = _identify_response_tensor_paths(
            subblocks,
            total_space=total_space,
            response_tensor_elements=response_tensor_elements,
            num_interactions=num_interactions
        )

        nlf = nonlinear_multipath_response_function_generators(response_tensor_paths, application_domain=application_domain)(
            frequencies, dipole_magnitudes, dipole_directions, couplings, orientational_averaging,

        )

    if response_function_class is None:
        response_function_class = NonlinearResponseFunction

    return response_function_class(
        nlf,
        central_frequency=central_frequency,
        frequency_unit=frequency_unit,
        time_unit=time_unit,
        application_domain=application_domain
    )

class TwoDimensionalIRResponseFunction(NonlinearResponseFunction):
    def call_freqs(self, w1, t2, w3,
                   *,
                   apply_ifft=True,
                   return_freqs=True,
                   shift_ifft=True,
                   **sampling_opts):
        signal = super().call_freqs(w1, t2, w3,
                                    apply_ifft=apply_ifft,
                                    return_freqs=return_freqs,
                                    shift_ifft=False,
                                    **sampling_opts)
        if self._dom == NonlinearReponseApplicationDomain.Time:
            if not apply_ifft:
                return signal
            else:
                if return_freqs:
                    freqs, signal = signal
                else:
                    freqs = None
                # rearrange
                signals = []
                if 'rephasing' in signal:
                    signals.append(np.flip(signal['rephasing'], axis=0))
                if 'non-rephasing' in signal:
                    signals.append(signal['non-rephasing'])
                S = sum(signals)
                if shift_ifft:
                    S = scipy.fft.fftshift(S)

                if return_freqs:
                    return freqs, S
                else:
                    return S
        else:
            if return_freqs:
                freqs, signal = signal
            else:
                freqs = None
            signals = []
            if 'rephasing' in signal:
                signals.append(np.flip(signal['rephasing'], axis=0))
            if 'non-rephasing' in signal:
                signals.append(signal['non-rephasing'])
            S = sum(signals)
            if return_freqs:
                return freqs, S
            else:
                return S

    def get_spectrum(self, w1, t2, w3, **sampling_opts):
        from ..Spectra import TwoDimensionalSpectrum
        (w1, w3), signal = self.call_freqs(w1, t2, w3, **sampling_opts)
        return TwoDimensionalSpectrum(w1, w3, np.real(signal))

experiment_defaults = {
    "2dir": {
        'num_interactions': 4,
        'paths': {
            'rephasing': [
                [1, 0, 1, 0],
                [1, 1, 0, 0],
                [1, 0, 0, 0],
            ],
            'non-rephasing': [
                [0, 1, 0, 0],
                [0, 0, 0, 0],
                [0, 1, 1, 0],
            ]
        },
        'selection_rules' : {
            'rephasing': [
                # have to match the paths
                # one set of selection rules per interactions
                [[(1,)], [(1,)], [(-1,)], [(1,)]],
                [[(1,)], [(-1,)], [(1,)], [(-1,)]],
                [[(1,)], [(1,)], [(1,)], [(1,)]]
            ],
            'non-rephasing': [
                # see above
                [[(1,)], [(1,)], [(-1,)], [(-1,)]],
                [[(1,)], [(-1,)], [(1,)], [(-1,)]],
                [[(1,)], [(1,)], [(1,)], [(-1,)]]
            ]
        },
        'polarization': 'XXYY',
        'response_function_class': TwoDimensionalIRResponseFunction
    }
}
def experimental_response_generator(
        transition_data,
        experiment_type='2dir',
        included_signals=None,
        **opts
):
    if experiment_type not in experiment_defaults:
        known_experiments = list(experiment_defaults.keys())
        raise ValueError(f"unknown experiment type `{experiment_type}`, available types are `{known_experiments}`")
    opts = collections.ChainMap(opts, experiment_defaults[experiment_type])
    if included_signals is not None:
        opts = dict(opts)
        paths = opts['paths']
        if not isinstance(paths, dict):
            raise ValueError(f"can't include only signals `{included_signals}` if `paths` is not a dict")
        opts['paths'] = {k:paths[k] for k in included_signals}
        selection_rules = opts.get('selection_rules')
        if isinstance(selection_rules, dict):
            opts['selection_rules'] = {k:selection_rules[k] for k in included_signals}
    return nonlinear_response_generators(
        transition_data,
        **opts
    )
