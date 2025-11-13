import collections
import itertools
import numpy as np
import scipy.fft
import enum

from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from .. import BasisReps as breps

import McUtils.Devutils as dev
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
                                   realign_transition_moments=False,
                                   include_implied_frequencies=None
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
            ).take_unique().as_sorted()

    if not hasattr(states, 'as_excitations'): # check that we were given a proper basis object
        int_states = nput.is_int(states[0])
        if not int_states:
            states = BasisStateSpace(
                HarmonicOscillatorProductBasis(len(states[0])),
                states
            ).take_unique().as_sorted()
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

    if include_implied_frequencies is None:
        include_implied_frequencies = needs_frequencies
    if include_implied_frequencies:
        for _ in range(5): # max 5 iters for now
            new_freqs = set()
            connections = np.where(np.abs(frequencies) > 1e-10)
            connection_graph = dict(zip(*nput.group_by(connections[1], connections[0])[0]))
            for i in connection_graph.keys():
                for j in connection_graph[i]:
                    for k in connection_graph[j]:
                        if i != k and abs(frequencies[i, k]) < 1e-10:
                            # E_k - E_j - (E_i - E_j)
                            frequencies[i, k] = frequencies[j, k] - frequencies[j, i]
                            frequencies[k, i] = -frequencies[i, k]
                            new_freqs.add((i,k))
            if len(new_freqs) == 0:
                break

    return TransitionData(states, frequencies, transition_moments, couplings)

def get_interaction_basis(initial_states:BasisStateSpace, *, selection_rules, **filter_opts):
    def _apply_rules(space, rules):
        if hasattr(space, 'representative_space'):
            space = space.to_single(include_representative=False).take_unique()
        if len(space) == 0:
            return None
        else:
            return space.apply_selection_rules(rules, **filter_opts)

    if nput.is_int(selection_rules[0][0][0]): # one path supplied
        bases = [initial_states]
        space = initial_states
        for rules in selection_rules:
            new = _apply_rules(bases[-1], rules)
            if new is None: return None
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

default_selection_rules = ((1,), (-1,))
def prep_liouville_spaces(initial_states, paths, num_interactions=None,
                          phases=None,
                          selection_rules=None, **filter_opts):
    smol = not isinstance(paths, dict)
    if smol:
        paths = {None:paths}
    paths = {
        k: (
            liouville_pathways(p)
                if nput.is_int(p) else
            np.asanyarray(p)
        )
        for k, p in paths.items()
    }
    if selection_rules is None:
        if phases is not None:
            if not isinstance(phases, dict):
                phases = {k:phases for k in paths.keys()}
            selection_rules = {}
            for k,phase_set in phases.items():
                path_set = paths[k]
                sub_rules = []
                if phase_set[0] is None or nput.is_int(phase_set[0]):
                    phase_set = [phase_set] * len(path_set)
                for path, phase in zip(path_set, phase_set):
                    sel_rules = []
                    for pk, pi in zip(path, phase):
                        if pi is None:
                            sel_rules.append(((-1,), (1,)))
                        elif pi == 1:
                            if pk == 0:
                                sel_rules.append(((1,),))
                            else:
                                sel_rules.append(((-1,),))
                        elif pi == -1:
                            if pk == 0:
                                sel_rules.append(((-1,),))
                            else:
                                sel_rules.append(((1,),))
                        elif pi == 0:
                            if pk == 0:
                                sel_rules.append(())
                            else:
                                sel_rules.append(())
                        else:
                            raise ValueError(f"unknown phase descriptor {pi}")
                    sub_rules.append(tuple(sel_rules))
                selection_rules[k] = sub_rules
        else:
            selection_rules = default_selection_rules
    smol_rules = not isinstance(selection_rules, dict)
    if smol_rules:
        selection_rules = {k:selection_rules for k in paths.keys()}
    all_bases = {}
    path_map = {}
    if num_interactions is None:
        num_interactions = len(next(iter(paths.values()))[0])
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
                    basis = get_interaction_basis(initial_states, selection_rules=rules_tuple, **filter_opts)
                    if basis is not None:
                        all_bases[rules_tuple] = basis
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
            groups = [all_inds.get(x) for x in path_map[(k, tuple(old))]]
            if any(g is None for g in groups):
                continue
            inds = [[g[i] for g, i in zip(groups, p)] for p in path]
            if len(inds) < num_interactions + 1:
                inds.append([groups[0][path[-1][0]+1], groups[1][path[-1][1]+1]])

            subinds.append([sign, old, inds])
        path_paths[k] = subinds
    if smol:
        path_paths = path_paths[None]
    return total_space, path_paths

def enumerate_state_paths(liou_path, index_set, num_interactions=None, enforce_pure_states=True):
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
            if enforce_pure_states:
                for l in left_nxt[left]:
                    new_path = cur_path + [(l, right)]
                    new_trans = transitions + [(left, l)]
                    if l == right:
                        yield new_path, new_trans
                        break
                else:
                    for r in right_nxt[right]:
                        new_path = cur_path + [(left, r)]
                        new_trans = transitions + [(right, r)]
                        if left == r:
                            yield new_path, new_trans
                            break
            else:
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
                                    num_interactions=None,
                                    enforce_pure_states=True):
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
        for spath, trans in enumerate_state_paths(p, inds,
                                                  num_interactions=num_interactions,
                                                  enforce_pure_states=enforce_pure_states):
            if spath[-1] in response_tensor_elements:
                response_tensor_elements[spath[-1]].append((spath[1:-1], trans, sign))

    return response_tensor_elements

class NonlinearReponseApplicationDomain(enum.Enum):
    Time = "time"
    Frequency = "frequency"


default_zero_dipole_cutoff = 1e-12
def _get_transition_dict(
        states,
        dipole_magnitudes,
        allowed_bands=None,
        zero_dipole_cutoff=None):
    if zero_dipole_cutoff is None:
        zero_dipole_cutoff = default_zero_dipole_cutoff
    true_cutoff = np.sqrt(zero_dipole_cutoff)
    mask = dipole_magnitudes > true_cutoff
    if allowed_bands is not None:
        quanta = np.sum(states, axis=1)
        diffs = quanta[:, np.newaxis] - quanta[np.newaxis, :]
        alt_mask = np.full(mask.shape, False)
        for b in allowed_bands:
            alt_mask = alt_mask | (diffs == b)
        mask = mask & alt_mask
    dipole_couplings = np.where(mask)
    return dict(zip(*nput.group_by(dipole_couplings[1], dipole_couplings[0])[0]))
def _prep_default_path_response_data(
        response_function_data,
        frequencies, dipole_magnitudes, dipole_directions, couplings, orientational_averaging,
        zero_dipole_cutoff=None,
        keys=None
):
    if zero_dipole_cutoff is None:
        zero_dipole_cutoff = default_zero_dipole_cutoff
    rfd = {}
    if keys is not None:
        response_function_data = {k:response_function_data[k] for k in response_function_data.keys() & set(keys)}
    for k, tdata in response_function_data.items():
        new_terms = []
        for sign, path, transitions in tdata:
            mag = np.prod([dipole_magnitudes[i, j] for i, j in transitions])
            if mag < zero_dipole_cutoff: continue
            prefactor = sign * mag * orientational_averaging(
                dipole_directions[tuple(i for i, j in transitions), tuple(j for i, j in transitions)]
            )
            new_terms.append(
                [
                    prefactor,
                    [frequencies[i, j] for i, j in path],
                    [couplings[i, j] for i, j in path]
                ]
            )
        rfd[k] = new_terms

    return rfd

def _default_path_response_function(response_function_data, application_domain):
    smol = not isinstance(response_function_data, dict)
    if smol:
        response_function_data = {None:response_function_data}
    if application_domain == NonlinearReponseApplicationDomain.Time:
        def response(*times):
            response_tensors = {}
            for k,d in response_function_data.items():
                resp = 0
                for mag, freqs, cups in d:
                    resp += mag * np.prod([
                        np.exp((f * 1j - c) * t)
                        for f, c, t in zip(freqs, cups, times)
                    ], axis=0)
                response_tensors[k] = resp
            if smol:
                response_tensors = response_tensors[None]
            return response_tensors
        return response
    elif application_domain == NonlinearReponseApplicationDomain.Frequency:
        def response(*fs):
            #TODO: check number of pulses
            response_tensors = {}
            for k,d in response_function_data.items():
                resp = 0
                for mag, freqs, cups in d:
                    resp += mag * np.prod([
                        ((w + f)*1j + c) / ((w + f) ** 2 + c ** 2 + 1e-14)  # just make it huge if we're at a pole
                        for f, c, w in zip([freqs[0], freqs[-1]], [cups[0], cups[-1]], [fs[0], fs[-1]])
                    ], axis=0) * np.prod([
                        np.exp((f * 1j - c) * t)
                        for f, c, t in zip(freqs[1:-1], cups[1:-1], fs[1:-1])
                    ], axis=0)
                response_tensors[k] = resp
            if smol:
                response_tensors = response_tensors[None]
            return response_tensors
        return response
    else:
        raise ValueError(f"unknown application domain {application_domain}")

def _cleanup_tensor_reponse_paths(response_tensor_paths):
    res = {}
    for transition_paths in response_tensor_paths.values():
        for path, transitions, sign in transition_paths:
            res[tuple(path)] = [
                sign,
                path,
                transitions
            ]
    return sorted(res.values(), key=lambda x:x[1])

def _complete_repsonse_function_generator(total_space, liouville_paths,
                                          frequencies,
                                          dipole_magnitudes, dipole_directions,
                                          couplings,
                                          orientational_averaging,
                                          *,
                                          application_domain,
                                          response_tensor_elements,
                                          num_interactions,
                                          zero_dipole_cutoff=None
                                          ):
    if zero_dipole_cutoff is None:
        zero_dipole_cutoff = default_zero_dipole_cutoff
    if isinstance(liouville_paths, dict):
        subpaths = {}
        for k, sububs in liouville_paths.items():
            response_tensor_paths = _identify_response_tensor_paths(
                sububs,
                total_space=total_space,
                response_tensor_elements=response_tensor_elements,
                num_interactions=num_interactions
            )

            subpaths[k] = _cleanup_tensor_reponse_paths(response_tensor_paths)

        # print(len(subpaths['non-rephasing']), [s[1] for s in subpaths['non-rephasing']])
        response_function_data = _prep_default_path_response_data(
            subpaths,
            frequencies, dipole_magnitudes,
            dipole_directions,
            couplings,
            orientational_averaging,
            zero_dipole_cutoff=zero_dipole_cutoff
        )
    else:
        response_tensor_paths = _identify_response_tensor_paths(
            liouville_paths,
            total_space=total_space,
            response_tensor_elements=response_tensor_elements,
            num_interactions=num_interactions
        )

        subpaths = _cleanup_tensor_reponse_paths(response_tensor_paths)

        response_function_data = _prep_default_path_response_data(
            {None:subpaths},
            frequencies, dipole_magnitudes,
            dipole_directions,
            couplings,
            orientational_averaging,
            zero_dipole_cutoff=zero_dipole_cutoff
        )[None]

    return _default_path_response_function(
        response_function_data,
        application_domain
    )

def _simple_path_response(
        paths,
        frequencies,
        dipole_magnitudes, dipole_directions,
        couplings,
        orientational_averaging,
        *,
        application_domain,
        zero_dipole_cutoff=None,
        keys=None,
        **__
):
    if zero_dipole_cutoff is None:
        zero_dipole_cutoff = default_zero_dipole_cutoff
    response_function_data = _prep_default_path_response_data(
        paths,
        frequencies, dipole_magnitudes,
        dipole_directions,
        couplings,
        orientational_averaging,
        zero_dipole_cutoff=zero_dipole_cutoff,
        keys=keys
    )

    return _default_path_response_function(
        response_function_data,
        application_domain
    )

def _get_simple_2dir_paths(
        total_space,
        dipole_magnitudes,
        zero_dipole_cutoff=None
):
    # TODO: add example support for starting from a different initial state?
    state_vectors = total_space.excitations
    quanta = np.sum(state_vectors, axis=1)
    gs = np.where(quanta == 0)[0]
    if len(gs) == 0:
        raise ValueError("ground state is required for now")
    gs = gs[0]
    coupling_dict = _get_transition_dict(
        state_vectors,
        dipole_magnitudes,
        allowed_bands=[-1, 1],
        zero_dipole_cutoff=zero_dipole_cutoff
    )

    response_function_data = {
        "rephasing": [],
        "non-rephasing": []
    }
    for i in coupling_dict[gs]:
        # R1
        for j in coupling_dict[gs]:
            response_function_data["rephasing"].append(
                [  # R1
                    1,
                    [(gs, j), (i, j), (i, gs)],
                    [(gs, j), (gs, i), (j, gs), (i, gs)],
                ]
            )

            response_function_data["rephasing"].append(
                [  # R2
                    1,
                    [(gs, j), (gs, gs), (i, gs)],
                    [(gs, j), (j, gs), (gs, i), (i, gs)],
                ]
            )

            response_function_data["non-rephasing"].append(
                [  # R4
                    1,
                    [(j, gs), (j, i), (j, gs)],
                    [(gs, j), (gs, i), (i, gs), (j, gs)],
                ]
            )

            response_function_data["non-rephasing"].append(
                [  # R5
                    1,
                    [(j, gs), (gs, gs), (i, gs)],
                    [(gs, j), (j, gs), (gs, i), (i, gs)],
                ]
            )
            for k in np.setdiff1d(np.intersect1d(coupling_dict[i], coupling_dict[j]), [gs]):
                response_function_data["rephasing"].append(
                    [  # R3
                        -1,
                        [(gs, j), (i, j), (k, j)],
                        [(gs, j), (gs, i), (i, k), (j, k)],
                    ]
                )
                response_function_data["non-rephasing"].append(
                    [  # R6
                        -1,
                        [(j, gs), (j, i), (k, i)],
                        [(gs, j), (gs, i), (j, k), (i, k)],
                    ]
                )

    return response_function_data

def _simple_2d_ir_response(
        total_space,
        liouville_paths,
        frequencies,
        dipole_magnitudes, dipole_directions,
        couplings,
        orientational_averaging,
        *,
        application_domain,
        zero_dipole_cutoff=None,
        **__
):
    if zero_dipole_cutoff is None:
        zero_dipole_cutoff = default_zero_dipole_cutoff
    paths = _get_simple_2dir_paths(
        total_space,
        dipole_magnitudes,
        zero_dipole_cutoff=zero_dipole_cutoff
    )

    # print(len(paths['non-rephasing']), [s[1] for s in paths['non-rephasing']])
    return _simple_path_response(
        paths,
        frequencies,
        dipole_magnitudes, dipole_directions,
        couplings,
        orientational_averaging,
        application_domain=application_domain,
        zero_dipole_cutoff=zero_dipole_cutoff,
        keys=liouville_paths.keys()
    )

def _signed_transition_data_from_paths(raw_paths):
    paths = []
    for p in raw_paths:
        right_ops = np.diff([r[1] for r in p]) != 0
        sign = (-1)**np.sum(right_ops)
        transitions = [
            (p[i][1], p[i+1][1])
                if r else
            (p[i][0], p[i+1][0])
            for i,r in enumerate(right_ops)
        ]
        paths.append(
            [
                sign,
                p[1:],
                transitions
            ]
        )
    return paths

def _enumerate_test_paths(states, dipole_magnitudes, phase=None, ordering=None, zero_dipole_cutoff=None):
    if zero_dipole_cutoff is None:
        zero_dipole_cutoff = default_zero_dipole_cutoff
    # sample path function adapted from Eitan
    band = np.sum(states, axis=1)

    if phase is None:
        phase = (None, None, None)
    elif phase == 'rephasing':
        phase = ('-', '+', '+')
    elif phase == 'non-rephasing':
        phase = ('+', '-', '+')
    elif phase == 'anharm':
        phase = ('+', '+', '-')

    if ordering is None:
        ordering = (None, None, None)
    elif isinstance(ordering, str) and len(ordering) >= 2:
        if ordering[1] == '1':
            if len(ordering) == 2:
                ordering = ('L', 'R', 'R')
            else:
                ordering = ('R', 'L', 'L')
        elif ordering[1] == '2':
            if len(ordering) == 2:
                ordering = ('R', 'L', 'R')
            else:
                ordering = ('L', 'R', 'L')
        elif ordering[1] == '3':
            if len(ordering) == 2:
                ordering = ('R', 'R', 'L')
            else:
                ordering = ('L', 'L', 'R')
        elif ordering[1] == '4':
            if len(ordering) == 2:
                ordering = ('L', 'L', 'L')
            else:
                ordering = ('R', 'R', 'R')

    nstate = len(band)

    # get all possible next steps in the pathway
    def getAllNext(mn, j):
        m, n = mn
        ori = phase[j]
        pat = ordering[j]
        result = []

        if (pat is None) or (pat == 'L'):
            # for mm in range(nstate):
            for mm in coupling_dict.get(m, []):
                if (
                        (ori is None and abs(band[mm] - band[m]) == 1)
                        or (ori == '+' and band[mm] == band[m] + 1)
                        or (ori == '-' and band[mm] == band[m] - 1)
                ):
                    result.append((mm, n))

        if (pat is None) or (pat == 'R'):
            # for nn in range(nstate):
            for nn in coupling_dict.get(n, []):
                if (
                        (ori is None and abs(band[nn] - band[n]) == 1)
                        or (ori == '+' and band[nn] == band[n] - 1)
                        or (ori == '-' and band[nn] == band[n] + 1)
                ):
                    result.append((m, nn))
                del nn

        return result

    mask = dipole_magnitudes > zero_dipole_cutoff
    # add extra selection rule restrictions here if desired
    # mask = mask & np.abs(band[:, np.newaxis] - band[np.newaxis, :]) == 1
    dipole_couplings = np.where(mask)
    coupling_dict = dict(zip(*nput.group_by(dipole_couplings[1], dipole_couplings[0])[0]))

    # get all possible next steps in the pathway
    def next_path_steps(mn, ori, pat):
        m, n = mn
        if (pat is None) or (pat == 'L'):
            for mm in coupling_dict.get(m, []):
                if (
                        ori is None
                        # or (ori == '+' and band[mm] == band[m] + 1)
                        # or (ori == '-' and band[mm] == band[m] - 1)
                        or (ori == '+' and band[mm] > band[m])
                        or (ori == '-' and band[mm] < band[m])
                ):
                    yield (mm, n)

        if (pat is None) or (pat == 'R'):
            for nn in coupling_dict.get(n, []):
                if (
                        ori is None
                        # or (ori == '+' and band[nn] == band[n] - 1)
                        # or (ori == '-' and band[nn] == band[n] + 1)
                        or (ori == '+' and band[nn] < band[n])
                        or (ori == '-' and band[nn] > band[n])
                ):
                    yield (m, nn)

    # construct pathways
    queue = collections.deque()
    queue.append([1, 0, [(0, 0)]])
    nstep = len(ordering)
    while queue:
        dipole, cur, path = queue.pop()
        if cur == nstep:
            yield path
        else:
            prev = path[-1]
            # for new in next_path_steps(prev, phase[cur], ordering[cur]):
            for new in getAllNext(prev, cur):
                if new[0] == prev[0]: # operated on the right
                    transition = (prev[1], new[1])
                else:
                    transition = (prev[0], new[0])
                dipole_new = dipole*dipole_magnitudes[transition[0], transition[1]]
                if dipole_new < zero_dipole_cutoff:
                    continue
                else:
                    queue.append([dipole_new, cur + 1,path + [new]])

def _alternate_test_reponse(
        total_space,
        liouville_paths,
        frequencies,
        dipole_magnitudes, dipole_directions,
        couplings,
        orientational_averaging,
        *,
        application_domain,
        zero_dipole_cutoff=None,
        **__
):
    paths = {
        k:_signed_transition_data_from_paths(list(_enumerate_test_paths(total_space.excitations, dipole_magnitudes,  phase=k)))
        for k in liouville_paths.keys()
    }

    # print(len(paths['non-rephasing']), sorted([p[1] for p in paths['non-rephasing']]))

    return _simple_path_response(
        paths,
        frequencies,
        dipole_magnitudes, dipole_directions,
        couplings,
        orientational_averaging,
        application_domain=application_domain,
        zero_dipole_cutoff=zero_dipole_cutoff
    )

class NonlinearResponseFunction:
    def __init__(self, response_generator, driving_frequency=0,
                 frequency_unit='Hartrees',
                 time_unit='PicoSeconds',
                 application_domain: NonlinearReponseApplicationDomain = "time"):
        self.caller = response_generator
        self.center = driving_frequency
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
            t1, t3 = np.meshgrid(t1_times, t3_times, indexing='xy')
            t2 = np.full(t1.shape, t2)
            t1, t2, t3 = [t * self.conv for t in [t1, t2, t3]] # w[cm-1] * conv(cm-1, Hz) * (t[ps] * conv(ps, s)) -> unitless
            if isinstance(self.caller, dict):
                signals = {k:f(t1, t2, t3) for k,f in self.caller.items()}
            else:
                signals = self.caller(t1, t2, t3)
            return signals
        else:
            raise NotImplementedError("calling with times with a frequency domain response function not currently supprted")

    @classmethod
    def _freq_conv(cls, freq_unit, time_unit):
        return UnitsData.convert(freq_unit, "Hertz") * UnitsData.convert(time_unit, 'Seconds') * (2*np.pi)
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
                np.array([w1_freq_min, w1_freq_max]),
                t2,
                np.array([w3_freq_min, w3_freq_max]),
                **sampling_options
            )
            signal = self.call_times(t1, t2, t3)
            if apply_ifft:
                return self._apply_ifft(t1, t3, signal, shift=shift_ifft, return_freqs=return_freqs)
            else:
                return (w1_freqs, w3_freqs), (t1, t3), signal
        elif self._dom == NonlinearReponseApplicationDomain.Frequency:
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
            w1, w3 = np.meshgrid(w1_freqs, w3_freqs, indexing='xy')
            t2 = np.full(w1.shape, t2 * self.conv)
            if isinstance(self.caller, dict):
                signals = {k:f(w1, t2, w3) for k,f in self.caller.items()}
            else:
                signals = self.caller(w1, t2, w3)
            if return_freqs:
                return (w1_freqs+self.center, w3_freqs+self.center), signals
            else:
                return signals
        else:
            raise ValueError(f"unknown application domain {self._dom}")

    def _get_fft_freqs(self, n, dt):
        freqs = scipy.fft.fftfreq(n, dt)
        freqs = scipy.fft.fftshift(freqs)
        conv = (2*np.pi) / self.conv # maybe should do this directly for stability...
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
            signal = signal.copy()
            signal[:, 0] /= 2
            signal[0, :] /= 2
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
                                  phases=None,
                                  selection_rules=None,
                                  polarization=None,
                                  initial_states=None,
                                  state_filter_opts=None,
                                  response_tensor_elements=None,
                                  response_function_generator=None,
                                  driving_frequency=True,
                                  frequency_unit="Hartrees",
                                  time_unit="PicoSeconds",
                                  application_domain="time",
                                  response_function_class=None,
                                  **state_opts
                                  ):
    td = prep_nonlinear_transition_data(transition_data, **state_opts)
    if driving_frequency is True:
        bands = np.sum(td.states.excitations, axis=1)
        gs = np.where(bands == 0)[0]
        if len(gs) == 0:
            raise ValueError("ground state is required for driving frequency determination")
        fundamentals = np.where(bands == 1)[0]
        if len(gs) == 0:
            raise ValueError("at least one fundamental is required for driving frequency determination")
        subfreqs = td.frequencies[gs, fundamentals]
        driving_frequency = np.average(subfreqs[subfreqs > 1e-10])

    if driving_frequency is False or driving_frequency is None:
        driving_frequency = 0
    else:
        freqs = td.frequencies.copy()
        bands = np.sum(td.states.excitations, axis=1)
        band_diffs = bands[np.newaxis, :] - bands[:, np.newaxis]
        sel = np.abs(freqs) > 1e-10
        freqs[sel] = (freqs[sel] - band_diffs[sel] * driving_frequency)
        td = td._replace(frequencies=freqs)

    states = td.states
    if not hasattr(states, 'as_excitations'): # path enumeration requires
        raise ValueError("state vectors are required")

    if polarization is not None:
        polarization = interpret_polarization(polarization) # TODO: ensure norm
        orientational_averaging = four_wave_averaging_function(polarization)
    else:
        orientational_averaging = lambda d: 1

    if isinstance(paths, dict):
        p_0 = next(iter(paths.values()))
    else:
        p_0 = paths
    needs_paths = (
        # TODO: handle the most generic case too?
        nput.is_int(p_0) or nput.is_int([0][0])
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
                                                   phases=phases,
                                                   selection_rules=selection_rules,
                                                   num_interactions=num_interactions,
                                                   **state_filter_opts)

    _, frequencies, transition_moments, couplings = expand_transition_data(td, total_space)
    dipole_directions, dipole_magnitudes = nput.vec_normalize(transition_moments, return_norms=True)

    # print("!")
    # print(
    #     sorted(list(_enumerate_test_2d_paths(total_space.excitations, dipole_magnitudes, phase='nr')))
    # )

    application_domain = NonlinearReponseApplicationDomain(application_domain)
    if response_function_generator is None:
        response_function_generator = _complete_repsonse_function_generator
    elif dev.str_is(response_function_generator, "simple2dir"):
        response_function_generator = _simple_2d_ir_response
    elif dev.str_is(response_function_generator, "alt2d"):
        response_function_generator = _alternate_test_reponse
    nlf = response_function_generator(
        total_space,
        subblocks,
        frequencies,
        dipole_magnitudes,
        dipole_directions,
        couplings,
        orientational_averaging,
        application_domain=application_domain,
        response_tensor_elements=response_tensor_elements,
        num_interactions=num_interactions
    )

    if response_function_class is None:
        response_function_class = NonlinearResponseFunction

    return response_function_class(
        nlf,
        driving_frequency=driving_frequency,
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
                   default_frequency_divisions=100,
                   **sampling_opts):
        if self._dom == NonlinearReponseApplicationDomain.Time:
            signal = super().call_freqs(w1, t2, w3,
                                        apply_ifft=apply_ifft,
                                        return_freqs=return_freqs,
                                        shift_ifft=False,
                                        **sampling_opts)
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
                    signals.append(np.flip(signal['rephasing'], axis=1))
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
            w1_min, w1_max = w1[:2]
            if len(w1) == 3:
                w1_step = w1[2]
            else:
                w1_step = default_frequency_divisions
            w3_min, w3_max = w3[:2]
            if len(w3) == 3:
                w3_step = w3[2]
            else:
                w3_step = default_frequency_divisions
            dw = min([(w1_max - w1_min)/w1_step, (w3_max - w3_min)/w3_step])
            w = max([
                abs(w1_max - self.center),
                abs(w1_min - self.center),
                abs(w3_max - self.center),
                abs(w3_min - self.center),
            ])
            nstep = int(np.ceil(2*w/dw))
            w1 = w3 = (self.center - w, self.center + w, nstep)
            signal = super().call_freqs(w1, t2, w3,
                                        apply_ifft=apply_ifft,
                                        return_freqs=return_freqs,
                                        shift_ifft=False,
                                        **sampling_opts)
            if return_freqs:
                freqs, signal = signal
            else:
                freqs = None
            signals = []
            if 'rephasing' in signal:
                signals.append(np.flip(signal['rephasing'], axis=1))
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
            'rephasing': 3,
            # 'rephasing': [
            #     [1, 0, 1],
            #     [1, 1, 0],
            #     [1, 0, 0],
            # ],
            'non-rephasing': 3
            # 'non-rephasing': [
            #     [0, 1, 0],
            #     [0, 0, 0],
            #     [0, 1, 1]
            # ]
        },
        # 'selection_rules' : {
        #     'rephasing': [
        #         # have to match the paths
        #         # one set of selection rules per interactions
        #         [[(1,)], [(1,) ], [(-1,)], [(-1,) ]],
        #         [[(1,)], [(-1,)], [(1,) ], [(-1,)]],
        #         [[(1,)], [(1,) ], [(1,) ], [(-1,) ]]
        #     ],
        #     'non-rephasing': [
        #         # see above
        #         [[(1,)], [(1,) ], [(-1,)], [(-1,)]],
        #         [[(1,)], [(-1,)], [(1,) ], [(-1,)]],
        #         [[(1,)], [(1,) ], [(1,) ], [(-1,)]]
        #     ]
        # },
        'phases': {
            'rephasing':[-1, 1, 1],
            'non-rephasing':[1, -1, 1],
        },
        # 'selection_rules' : {
        #     'rephasing': [
        #         # have to match the paths
        #         # one set of selection rules per interactions
        #         [[(1,)], [(1,) ], [(-1,)]],# [(-1,) ]],
        #         [[(1,)], [(-1,)], [(1,) ]],# [(-1,)]],
        #         [[(1,)], [(1,) ], [(1,) ]],# [(-1,) ]]
        #     ],
        #     'non-rephasing': [
        #         # see above
        #         [[(1,)], [(-1,) ], [(-1,)]],# [(-1,)]],
        #         [[(1,)], [(1,) ], [(-1,)]],# [(-1,)]],
        #         [[(1,)], [(-1,)], [(1,) ]],# [(-1,)]]
        #     ]
        # },
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
