from __future__ import annotations

import collections
import itertools
import os
import json
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
    "experimental_response_generator",
    "prep_vpt_response_data",
    "prep_vpt_response_data_from_log"
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

def _vpt_response_data_to_records(transition_dict):
    """
    Converts a `transition_dict` (as returned by `prep_vpt_response_data`) into
    a JSON-serializable list of records, one per transition, each of the form
    `{"state": [state_i, state_j], "frequency": ..., "transition_moment": [...], ...}`
    -- i.e. the `(state_i, state_j)` dict key gets folded into the record itself
    under a `"state"` key rather than kept as a (non-JSON-safe) tuple key.

    :param transition_dict: mapping of `(state_i, state_j)` state-vector tuples to
        per-transition data (at least `frequency` and `transition_moment`)
    :type transition_dict: dict
    :return: a list of JSON-safe transition records
    :rtype: list[dict]
    """
    records = []
    for (si, sj), data in transition_dict.items():
        rec = dict(data)
        if 'transition_moment' in rec and rec['transition_moment'] is not None:
            rec['transition_moment'] = np.asarray(rec['transition_moment']).tolist()
        if 'frequency' in rec and rec['frequency'] is not None:
            rec['frequency'] = float(rec['frequency'])
        rec['state'] = [list(si), list(sj)]
        records.append(rec)
    return records

def _vpt_response_data_from_records(records):
    """
    Inverse of `_vpt_response_data_to_records`: reconstitutes a `transition_dict`
    (state-tuple-pair keys, `transition_moment` as an `np.ndarray`) from the
    JSON-safe record list that gets written to/read from `output_file`.

    :param records: the JSON-decoded list of transition records
    :type records: list[dict]
    :return: the reconstituted `transition_dict`
    :rtype: dict
    """
    transition_dict = {}
    for rec in records:
        rec = dict(rec)
        si, sj = rec.pop('state')
        si = tuple(int(x) for x in si)
        sj = tuple(int(x) for x in sj)
        if 'transition_moment' in rec and rec['transition_moment'] is not None:
            rec['transition_moment'] = np.array(rec['transition_moment'])
        transition_dict[(si, sj)] = rec
    return transition_dict

def _prep_vpt_target_states(freqs, target_states=None, max_freq=None, max_quanta=2):
    """Resolve explicit or generated target states for ``prep_vpt_response_data``."""
    freqs = np.asanyarray(freqs)
    ndim = len(freqs)
    ground_state = tuple([0] * ndim)

    if target_states is None or isinstance(target_states, dict):
        state_opts = {} if target_states is None else target_states.copy()
        if 'max_freq' not in state_opts:
            if max_freq is None:
                max_freq = max_quanta * np.max(np.abs(freqs))
            state_opts['max_freq'] = max_freq
        if 'max_quanta' not in state_opts:
            # ``states_under_freq_threshold`` uses an exclusive upper bound,
            # while this helper's long-standing ``max_quanta`` option is inclusive.
            state_opts['max_quanta'] = max_quanta + 1
        raw_states = BasisStateSpace.states_under_freq_threshold(freqs, **state_opts)
    elif isinstance(target_states, BasisStateSpace):
        raw_states = target_states.excitations
    else:
        raw_states = target_states

    if not isinstance(raw_states, np.ndarray):
        raw_states = list(raw_states)
    raw_states = np.asanyarray(raw_states)
    if raw_states.size == 0:
        raw_states = np.empty((0, ndim), dtype=int)
    if raw_states.ndim == 1 and raw_states.shape == (ndim,):
        raw_states = raw_states[np.newaxis, :]
    if raw_states.ndim != 2 or raw_states.shape[1] != ndim:
        raise ValueError(
            f"target states must be a two-dimensional array with {ndim} columns; "
            f"got shape {raw_states.shape}"
        )
    if not np.issubdtype(raw_states.dtype, np.number):
        raise TypeError("target states must contain numeric quantum numbers")
    if np.any(~np.isfinite(raw_states)) or np.any(raw_states < 0):
        raise ValueError("target states must contain finite, non-negative quantum numbers")
    if np.any(raw_states != np.floor(raw_states)):
        raise ValueError("target states must contain integer quantum numbers")

    # Normalize to hashable tuples and remove duplicates without disturbing the
    # caller/BasisStateSpace ordering used by the downstream VPT runners.
    state_list = []
    seen = set()
    for state in raw_states:
        state = tuple(int(x) for x in state)
        if state not in seen:
            state_list.append(state)
            seen.add(state)
    if ground_state not in seen:
        state_list.insert(0, ground_state)
    return state_list

def prep_vpt_response_data(system,
                            max_freq=None,
                            max_quanta=2,
                            initial_quanta=(0, 1),
                            return_wavefunctions=False,
                            output_file=None,
                            overwrite=False,
                            use_analytic=False,
                            target_states=None,
                            **vpt_opts
                            ):
    """
    Runs a VPT calculation over the states reachable from the ground state
    within `max_quanta` quanta of excitation and returns the resulting state
    energies/transition moments as a `transition_dict` in the format expected
    by `prep_nonlinear_transition_data`/`experimental_response_generator`.

    By default (`use_analytic=False`) this runs `VPTRunner.run_simple`, which
    computes a single dense "every initial state x every final state"
    transition-moment matrix (so e.g. a direct ground-state overtone
    transition ends up included alongside the fundamentals, if VPT gives it
    a nonzero moment). Passing `use_analytic=True` instead runs
    `AnalyticVPTRunner.run_simple`, which has a different calling convention:
    rather than one flat state list plus an `initial_states` seed, it wants
    an explicit list of `[initial_space, target_space]` block pairs, and only
    ever computes transition moments *within* each block -- there's no
    implicit dense any-initial-to-any-final matrix. To keep the same overall
    coverage as the classic branch as closely as that block-based API allows,
    one block is built per consecutive pair of quantum shells reachable from
    `initial_quanta` (e.g. ground -> one-quantum fundamentals, one-quantum ->
    two-quantum states, and so on up to `max_quanta`), mirroring the
    ground -> fundamentals -> overtones/combinations cascade a 2D-IR
    calculation actually needs. One consequence of this block structure:
    a *non-adjacent*-shell transition (like a direct ground -> two-quantum
    overtone) is only included by the classic branch, not the analytic one,
    since no block connects those two shells directly.

    The two branches are independent VPT implementations, so their results
    agree closely but not bit-for-bit: energies for shared states typically
    match to a small fraction of a wavenumber, and per-mode transition
    moments can come back with an overall sign flipped relative to the
    classic branch's (an arbitrary phase-convention difference between the
    two evaluators, verified against real water VPT data to have no effect
    on downstream intensities/spectra, which only ever depend on these
    moments through even, sign-invariant combinations).

    By default, the target state list comes from
    `BasisStateSpace.states_under_freq_threshold`, run over the system's
    harmonic normal-mode frequencies and capped at `max_quanta` total quanta
    of excitation. Passing explicit state vectors through `target_states`
    restricts the calculation to precisely those states. A `BasisStateSpace`
    may be passed directly, or a dictionary may be supplied as keyword options
    for `BasisStateSpace.states_under_freq_threshold`; dictionary values
    override the `max_freq`/`max_quanta` defaults. The ground state is always
    included. The classic branch seeds `initial_states` with every target
    state whose total quantum number is in `initial_quanta` (by default the
    ground state and every singly-excited fundamental); the analytic branch
    uses the same `initial_quanta` values to build quantum-shell blocks.

    If `output_file` is given and already exists, the cached `transition_dict`
    is loaded from it directly and returned *without running any VPT
    calculation* (unless `overwrite=True`, which always reruns and rewrites
    the file). If `output_file` is given and doesn't yet exist (or
    `overwrite=True`), the calculation is run as usual and the resulting
    `transition_dict` is saved to `output_file` as JSON before being returned.
    A loaded-from-cache result has no associated wavefunctions/corrections
    object, so `return_wavefunctions` yields `None` in its place in that case.

    :param system: a molecule/system spec (path, `Molecule`, or `VPTSystem`) to run VPT on
    :type system: str | list | Molecule | VPTSystem
    :param max_freq: the maximum total (harmonic) excitation energy to include when generating
        the target state list, in the same units as the system's normal-mode
        frequencies (Hartrees). Defaults to `max_quanta` times the largest
        normal-mode frequency, which is always enough to admit every state
        satisfying the `max_quanta` cutoff below.
    :type max_freq: float | None
    :param max_quanta: the largest total number of vibrational quanta (summed over all
        modes) a target state is allowed to carry (inclusive)
    :type max_quanta: int
    :param initial_quanta: the total quantum numbers (again summed over modes) that qualify a
        state to be used as an `initial_states` seed for the VPT run (classic
        branch) or as a quantum-shell boundary to build a block across
        (analytic branch) -- by default the ground state (0 quanta) and every
        fundamental (1 quantum)
    :type initial_quanta: int | Iterable[int]
    :param target_states: optional target-state specification. May be an explicit iterable/array
        of full-dimensional excitation vectors, a `BasisStateSpace`, or a dictionary of keyword
        options passed to `BasisStateSpace.states_under_freq_threshold` (for example
        `{'max_freq': ..., 'max_quanta': 3, 'fixed_modes': [...]}`). When omitted, the legacy
        `max_freq`/inclusive-`max_quanta` generation behavior is retained. Dictionary-provided
        `max_quanta` is passed directly to `BasisStateSpace` and therefore uses its exclusive
        upper-bound convention. Duplicate states are removed and the ground state is added.
    :type target_states: Iterable[Iterable[int]] | BasisStateSpace | dict | None
    :param return_wavefunctions: if `True`, also return the raw results object from the VPT
        run alongside the `transition_dict` -- a `VPTWavefunctions` for the
        classic branch, or an `AnalyticPerturbationTheoryCorrections` for the
        analytic branch (`use_analytic=True`) -- or `None` if the result was
        loaded from `output_file` instead of computed
    :type return_wavefunctions: bool
    :param output_file: optional path to cache the resulting `transition_dict` as JSON
        (a list of `{"state": [state_i, state_j], "frequency":..., "transition_moment":...}`
        records, since JSON object keys can't be tuples). When this file already exists,
        it's loaded and returned as-is instead of rerunning the VPT calculation, unless
        `overwrite=True`
    :type output_file: str | None
    :param overwrite: if `True`, always (re)run the calculation and overwrite `output_file`,
        even if it already exists
    :type overwrite: bool
    :param use_analytic: if `True`, run `AnalyticVPTRunner.run_simple` (symbolic/analytic VPT)
        instead of `VPTRunner.run_simple`, building one `[initial_space, target_space]`
        block per consecutive pair of quantum shells reachable from `initial_quanta`
        (see above for how this differs from the classic branch's dense matrix)
    :type use_analytic: bool
    :param vpt_opts: extra options forwarded to `VPTRunner.run_simple`/`VPTRunner.construct`
        (classic branch) or `AnalyticVPTRunner.run_simple`/`AnalyticVPTRunner.construct`
        (analytic branch, e.g. `full_surface_mode_selection`,
        `mixed_derivative_handling_mode`, `mixed_derivative_handle_zeros`,
        `mixed_derivative_warning_threshold`, `corrected_fundamental_frequencies`, `logger`)
    :type vpt_opts: dict
    :return: a `transition_dict` of the form `{(state_i, state_j): {'frequency':..., 'transition_moment':...}}`
        (in wavenumbers/a.u., respectively), suitable for `prep_nonlinear_transition_data`,
        or `(transition_dict, wfns)` if `return_wavefunctions` is set
    :rtype: dict | tuple[dict, 'VPTWavefunctions' | 'AnalyticPerturbationTheoryCorrections']
    """

    if output_file is not None and not overwrite and os.path.isfile(output_file):
        with open(output_file, 'r') as woof:
            transition_dict = _vpt_response_data_from_records(json.load(woof))
        if return_wavefunctions:
            return transition_dict, None
        return transition_dict

    from ..VPT2 import VPTRunner, VPTSystem

    vpt_system = system if isinstance(system, VPTSystem) else VPTSystem(system)
    freqs = vpt_system.mol.normal_modes.modes.freqs
    ndim = len(freqs)
    ground_state = tuple([0] * ndim)
    state_list = _prep_vpt_target_states(
        freqs,
        target_states=target_states,
        max_freq=max_freq,
        max_quanta=max_quanta
    )

    if nput.is_int(initial_quanta):
        initial_quanta = (initial_quanta,)
    initial_states = [s for s in state_list if sum(s) in initial_quanta]
    if ground_state not in initial_states:
        initial_states = [ground_state] + initial_states

    h2w = UnitsData.convert("Hartrees", "Wavenumbers")
    transition_dict = {}

    if use_analytic:
        from ..VPT2 import AnalyticVPTRunner

        # `AnalyticVPTRunner` wants an explicit list of `[initial_space,
        # target_space]` block pairs rather than a flat state list + an
        # `initial_states` seed, and it only computes transition moments
        # *within* each block -- build one block per consecutive quantum
        # shell boundary in `initial_quanta` (e.g. 0->1, 1->2, ...)
        by_quanta = collections.defaultdict(list)
        for s in state_list:
            by_quanta[sum(s)].append(list(s))

        blocks = []
        for k in sorted(set(initial_quanta)):
            if k not in by_quanta or (k + 1) not in by_quanta:
                continue
            blocks.append([by_quanta[k], by_quanta[k + 1]])

        if len(blocks) == 0:
            raise ValueError(
                f"no consecutive-quanta-shell blocks to run between "
                f"`initial_quanta={initial_quanta}` and `max_quanta={max_quanta}`"
            )

        corrs = AnalyticVPTRunner.run_simple(
            vpt_system,
            blocks,
            **vpt_opts
        )

        energies = corrs.energies * h2w
        # `corrs.transition_moments` is `[axis][block_idx] -> (n_init, n_final)`,
        # not the classic branch's single dense `(3, n_initial, n_total)` tensor
        tms = corrs.transition_moments

        for block_idx, (init_block, final_block) in enumerate(corrs.state_lists):
            init_block = [tuple(int(x) for x in s) for s in init_block]
            final_block = [tuple(int(x) for x in s) for s in final_block]
            init_inds = corrs.states.find(init_block)
            final_inds = corrs.states.find(final_block)
            for i, (si, ii) in enumerate(zip(init_block, init_inds)):
                for j, (sj, jj) in enumerate(zip(final_block, final_inds)):
                    if ii == jj:
                        continue
                    freq = energies[jj] - energies[ii]
                    if freq <= 0:
                        # keep only the "upward" direction of each pair; `prep_nonlinear_transition_data`
                        # infers the reverse (negative-frequency) transition automatically
                        continue
                    key = (si, sj)
                    if key in transition_dict:
                        continue
                    tm = np.array([tms[c][block_idx][i][j] for c in range(3)])
                    transition_dict[key] = {'frequency': float(freq), 'transition_moment': tm}

        wfns = corrs
    else:
        wfns = VPTRunner.run_simple(
            vpt_system,
            state_list,
            initial_states=initial_states,
            **vpt_opts
        )

        state_tuples = [tuple(int(x) for x in s) for s in wfns.corrs.states.excitations]
        energies = wfns.energies * h2w
        tms = wfns.transition_moments

        for n, init_idx in enumerate(wfns.initial_state_indices):
            si = state_tuples[init_idx]
            for j, sj in enumerate(state_tuples):
                if j == init_idx:
                    continue
                freq = energies[j] - energies[init_idx]
                if freq <= 0:
                    # keep only the "upward" direction of each pair; `prep_nonlinear_transition_data`
                    # infers the reverse (negative-frequency) transition automatically
                    continue
                key = (si, sj)
                if key in transition_dict:
                    continue
                tm = np.array([tms[k][n][j] for k in range(3)])
                transition_dict[key] = {'frequency': float(freq), 'transition_moment': tm}

    if output_file is not None:
        out_dir = os.path.dirname(output_file)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(output_file, 'w') as woof:
            json.dump(_vpt_response_data_to_records(transition_dict), woof)

    if return_wavefunctions:
        return transition_dict, wfns
    return transition_dict

def _parse_vpt_state_label(label):
    """
    Parse a `VPTAnalyzer`/`VPTAnalyzerLogParser` state label (e.g. `"1 0 0"`) into an
    excitation-quanta tuple (e.g. `(1, 0, 0)`), matching the tuple format used for
    `transition_dict` keys elsewhere in this module.

    Nothing in `Psience.VPT2.Analyzer` currently does this conversion -- `VPTAnalyzerLogParser`
    hands back state labels as the raw, whitespace-joined digit strings taken verbatim from the
    log table (e.g. via `.spectra`/`.transition_moment_corrections`), and every consumer is left
    to parse them itself. This is one of the concrete small gaps found while building
    `prep_vpt_response_data_from_log` below; a `state_label` <-> excitation-tuple helper like this
    one would be a reasonable thing to add directly to `VPTAnalyzer`/`VPTAnalyzerLogParser`.

    :param label: whitespace-separated per-mode quanta, e.g. `"1 0 0"`
    :type label: str
    :return: excitation-quanta tuple, e.g. `(1, 0, 0)`
    :rtype: tuple[int]
    """
    return tuple(int(x) for x in label.split())

def prep_vpt_response_data_from_log(log_file, max_freq=None, initial_quanta=(0, 1)):
    """
    Reconstructs a `transition_dict` (in the same format returned by `prep_vpt_response_data`)
    purely from a saved VPT2 *text log*, using the existing `Psience.VPT2.Analyzer.VPTAnalyzer`
    class to do the parsing, rather than from an in-memory `VPTWavefunctions`/`AnalyticPerturbationTheoryCorrections`
    result.

    This only works for logs produced by the **classic** `VPTRunner` (i.e. `logger=<path>` passed
    to `VPTRunner.run_simple`), and even then only after two real bugs in `VPTAnalyzerLogParser`
    were found and fixed while building this function (see `claude_drafts/vpt_analyzer_log_parsing_fixes.patch`):
    the `.tree` property never unwrapped the single outer `">>--- Starting Perturbation Theory Runner ---<<"`
    banner block that every such log is wrapped in, so no named table (`"IR Data"`, `"X Dipole Contributions"`,
    etc) was ever reachable; and `SpectrumBlockParser`/`TransitionMomentBlockParser.check_tag` only
    recognized a *leading-space* `" Initial State:"` header as the start of a new per-initial-state
    sub-block and didn't skip the dashed separator line between sub-blocks, which happened to work by
    accident for single-initial-state logs (e.g. the one existing reference fixture, `methanol_vpt_3.out`)
    but silently mis-parsed (or crashed on) *any* log with more than one initial state -- which is the
    normal case for this module, since `prep_vpt_response_data`'s classic branch always requests both
    ground- and one-quantum initial states.

    Two further limitations remain, and are NOT fixed here (see the patch notes above for the full
    writeup):

    - `AnalyticVPTRunner` logs cannot be parsed by `VPTAnalyzerLogParser` **at all** -- confirmed by
      actually generating one (checked in as `ci/tests/TestData/water_vpt_analytic.log`) and attempting to load
      it: `AnalyticVPTRunner.run_VPT` never emits the named log blocks (`"IR Data"`, `"X/Y/Z Dipole
      Contributions"`, etc) that `VPTAnalyzerLogParser` looks for by exact tag string; it instead logs
      a combined `"Transition Moments:"` table (via `format_transition_moment_table`) and leaves the
      energies/spectrum output untagged. Supporting this would require either teaching `AnalyticVPTRunner`
      to emit `VPTRunner`-compatible tagged blocks, or writing an entirely separate parser for its log
      format. This function raises a clear `ValueError` (rather than a bare `IndexError`) if pointed at
      such a log.
    - The *transition moments* this function recovers are only approximately correct for combination-band
      and overtone transitions (verified exact for all pure fundamentals, and within ~1e-3 for about
      70% of all transitions tested against the in-memory `prep_vpt_response_data(water_freq.fchk)`
      ground truth, with the rest off by up to ~30%). This traces to a third, separate bug: `VPTAnalyzerLogParser.reformat_tm_block`
      (via `load_term_counts`/`McUtils.Combinatorics.SymmetricGroupGenerator`) mis-slices the raw
      per-order dipole-correction columns of a `"X/Y/Z Dipole Contributions"` table row, silently
      dropping roughly half of the printed correction terms for a 10-column row -- so the *frequencies*
      this function returns are exact (verified against all 28 ground-truth transitions for water),
      but the transition moments should be treated as approximate unless/until that slicing bug is
      also fixed. That fix needs to trace through exactly how many dipole-derivative-order correction
      terms `VPTWavefunctions.format_dipole_contribs_tables` prints for a given expansion order, which
      is out of scope here.

    :param log_file: path to a text log produced by `VPTRunner.run_simple(..., logger=log_file)`
    :type log_file: str
    :param max_freq: if given, drop any reconstructed transition whose frequency (in cm^-1) exceeds this
    :type max_freq: float | None
    :param initial_quanta: which initial-state total-quanta values to keep transitions from (matches the
        same-named parameter of `prep_vpt_response_data`); the ground state and one-quantum blocks (`(0, 1)`)
        are always present when the log was generated the way `prep_vpt_response_data`'s classic branch
        generates them
    :type initial_quanta: tuple[int]
    :return: a `transition_dict` of the form `{(state_i, state_j): {'frequency':..., 'transition_moment':...}}`,
        in the same format as `prep_vpt_response_data`
    :rtype: dict
    """
    from ..VPT2 import VPTAnalyzer

    analyzer = VPTAnalyzer(log_file)
    parser = analyzer.log_parser

    try:
        spectra = parser.spectra
        tm_corrections = parser.transition_moment_corrections
    except (IndexError, KeyError) as e:
        raise ValueError(
            "could not parse a transition_dict from log file '{}': {} "
            "(note: VPTAnalyzerLogParser currently only supports logs from the classic `VPTRunner` -- "
            "`AnalyticVPTRunner` logs use a different, untagged table format and aren't supported)".format(
                log_file, e
            )
        ) from e

    if isinstance(spectra, dict):
        spectra = [spectra]
    if isinstance(tm_corrections, dict):
        tm_corrections = [tm_corrections]

    if len(spectra) != len(tm_corrections):
        raise ValueError(
            "parsed {} spectrum block(s) but {} transition-moment block(s) from '{}'; "
            "can't reliably pair these up".format(len(spectra), len(tm_corrections), log_file)
        )

    transition_dict = {}
    for spec_block, tm_block in zip(spectra, tm_corrections):
        fin_labels = spec_block['states']
        tm_axes = tm_block['corrections']  # [x, y, z] axis dicts, each from `reformat_tm_block`
        tm_labels = tm_axes[0]['states']

        # `VPTAnalyzerLogParser` doesn't currently record which initial state a parsed
        # transition-moment block belongs to; each such block's raw table does include exactly
        # one extra row beyond what's in the matching spectrum block, though -- the block's own
        # diagonal <initial|mu|initial> self-term (a real, nonzero permanent-dipole matrix element,
        # but not a "transition" so it's excluded from the "IR Data" spectrum table). Whichever
        # label appears in the TM block but not in the spectrum block's final-state list is
        # therefore this block's initial state.
        fin_label_set = set(fin_labels)
        init_candidates = [s for s in tm_labels if s not in fin_label_set]
        if len(init_candidates) != 1:
            raise ValueError(
                "couldn't uniquely infer the initial state of a parsed transition-moment block in "
                "'{}' (candidates: {}) -- VPTAnalyzer doesn't currently label these blocks directly, "
                "so this had to be inferred, and the inference failed here".format(log_file, init_candidates)
            )
        init_label = init_candidates[0]
        init_state = _parse_vpt_state_label(init_label)
        if sum(init_state) not in initial_quanta:
            continue

        rows = [k for k, s in enumerate(tm_labels) if s != init_label]

        for row, fin_label, (freq, _intensity) in zip(rows, fin_labels, spec_block['anharmonic']):
            if freq <= 0:
                # keep only the "upward" direction of each pair, matching `prep_vpt_response_data`;
                # the reverse (negative-frequency) transition is inferred downstream automatically
                continue
            if max_freq is not None and freq > max_freq:
                continue
            fin_state = _parse_vpt_state_label(fin_label)
            key = (init_state, fin_state)
            if key in transition_dict:
                continue
            tm = np.array([
                sum(arr[row].sum() for arr in axis_block['corrections'])
                for axis_block in tm_axes
            ])
            transition_dict[key] = {'frequency': float(freq), 'transition_moment': tm}

    return transition_dict

def get_interaction_basis(initial_states:BasisStateSpace, *, selection_rules, **filter_opts):
    def _apply_rules(space, rules, filter_opts):
        # `apply_selection_rules` returns a bare space when it was called without
        # `filter_space`, but returns a `(space, updated_filter)` tuple whenever
        # `filter_space` is supplied (see `SelectionRuleStateSpace.from_rules` and
        # the analogous unpacking in `AbstractStateSpace.get_representation_indices`).
        # We unpack that here and thread the (progressively narrowed) filter forward
        # into subsequent calls instead of reusing the caller's original filter_space.
        if hasattr(space, 'representative_space'):
            space = space.to_single(include_representative=False).take_unique()
        if len(space) == 0:
            return None
        else:
            new = space.apply_selection_rules(rules, **filter_opts)
            if not isinstance(new, breps.AbstractStateSpace):
                new, updated_filter = new
                if 'filter_space' in filter_opts:
                    filter_opts = dict(filter_opts, filter_space=updated_filter)
            return new, filter_opts

    if nput.is_int(selection_rules[0][0][0]): # one path supplied
        bases = [initial_states]
        space = initial_states
        cur_filter_opts = filter_opts
        for rules in selection_rules:
            res = _apply_rules(bases[-1], rules, cur_filter_opts)
            if res is None: return None
            new, cur_filter_opts = res
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
            cur_filter_opts = filter_opts
            for rules in rule_list:
                res = _apply_rules(basis[-1], rules, cur_filter_opts)
                if res is None: break
                new, cur_filter_opts = res
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
                                  ) -> NonlinearResponseFunction:
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
