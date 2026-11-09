"""Selection-rule screening of analytic WFCs used by a target property.

The screen contains no Hamiltonian matrix construction or numerical VPT solver.
Retained coefficients use the runner's existing analytic evaluator unchanged.
"""
import numpy as np
from McUtils import Numputils as nput
from ..BasisReps import BasisStateSpace, BasisStateSpaceFilter
from .Corrections import BasicAPTCorrections

__all__ = ['AnalyticPropertyDiscovery']


def uses_automatic_discovery(spec):
    from .DegeneracySpecs import StronglyCoupledDegeneracySpec
    if spec is True or isinstance(spec, StronglyCoupledDegeneracySpec):
        return True
    if isinstance(spec, (int, float, np.integer, np.floating)) and spec is not False:
        return True
    if isinstance(spec, str):
        return spec == 'auto'
    if isinstance(spec, dict):
        return 'wfc_threshold' in spec
    if isinstance(spec, (list, tuple)):
        return any(uses_automatic_discovery(s) for s in spec)
    return False


class AnalyticPropertyDiscovery:
    """Check property use before evaluating an analytic WFC magnitude.

    Property paths are described by the existing selection-rule planner, with
    cumulative perturbative order retained at every WFC dependency. The planner
    works with quanta manifolds, as does the existing classic property filter;
    it is a conservative structural test, not a test of numerical dipole zeros.
    Actual endpoint changes must also satisfy the analytic WFC selection rules.
    """
    def __init__(self, runner, packet, settings=None):
        self.runner = runner
        self.eval = runner.eval
        self.logger = runner.logger
        self.packet = packet
        self.settings = {} if settings is None else dict(settings)
        self.target = self.settings.get('target_property') or 'intensities'
        if self.target not in {'intensities', 'frequencies', 'wavefunctions'}:
            raise ValueError('unknown discovery target property: '+str(self.target))
        self.order = runner.order
        self._rules = {}
        self._uses = None

    def for_degeneracy_spec(self, spec):
        enabled = getattr(spec, 'filter_by_property', None)
        if enabled is None:
            return self
        return type(self)(self.runner, self.packet,
            dict(self.settings, strong_coupling_property_filter=enabled))

    def log_strong_couplings(self, *args, **kwargs):
        return self.runner.log_strong_couplings(*args, **kwargs)

    @staticmethod
    def _selection_rules(term):
        # Analytic evaluation also caches (change, shift) entries in .changes.
        # They are numerical cache keys, not additional selection rules.
        return tuple(change for change in term.changes
            if all(isinstance(c, (int, np.integer)) for c in change))

    def _prepare(self):
        if self._uses is not None:
            return
        solver = self.eval.solver
        hamiltonian = [self._selection_rules(h)
            for h in solver.hamiltonian_expansion[1:self.order+1]]
        if self.target == 'intensities':
            operators = solver.operator_expansion_terms(self.order)
            properties = [self._selection_rules(term) for term in operators]
        else:
            properties = None
        all_uses = {}
        for initial, target in self.packet.state_list_pairs:
            uses = BasisStateSpaceFilter.property_wavefunction_quanta(
                initial, target, hamiltonian, properties, order=self.order)
            for key, quanta in uses.items():
                all_uses.setdefault(key, set()).update(quanta)
        self._uses = all_uses

    def _wavefunction_rules(self, order):
        if order not in self._rules:
            self._rules[order] = self._selection_rules(
                self.eval.solver.wavefunction_correction(order))
        return self._rules[order]

    def endpoint_quanta(self, source, order):
        self._prepare()
        source_quanta = int(np.sum(source))
        reachable = {source_quanta+sum(rule) for rule in self._wavefunction_rules(order)
            if source_quanta+sum(rule) >= 0}
        if self.target == 'wavefunctions':
            return reachable
        return reachable.intersection(self._uses.get((source_quanta, order), ()))

    def allowed_endpoints(self, source, endpoints, order):
        """Return a cheap structural mask; never evaluate a coefficient here."""
        from .Analytic import PerturbationTheoryTerm
        endpoints = np.asarray(endpoints, dtype=int)
        allowed = self.endpoint_quanta(source, order)
        keep = np.isin(np.sum(endpoints, axis=-1), list(allowed))
        if not np.any(keep):
            return keep
        changes = self._wavefunction_rules(order)
        signatures = {tuple(PerturbationTheoryTerm.sorted_changes(rule)) for rule in changes
            if len(rule) <= endpoints.shape[-1]}
        positions = np.flatnonzero(keep)
        if not signatures:
            keep[positions] = False
            return keep
        diffs = endpoints[positions]-source
        width = max(map(len, signatures), default=0)
        maximum = max((sum(abs(c) for c in rule) for rule in signatures), default=0)
        # Exact integer selection-rule bounds, not a numerical magnitude cutoff.
        possible = ((np.count_nonzero(diffs, axis=-1) <= width)
                    & (np.sum(np.abs(diffs), axis=-1) <= maximum))
        keep[positions[~possible]] = False
        positions = positions[possible]
        diffs = diffs[possible]
        if not len(positions):
            return keep
        if width == 0:
            return keep
        permutations = PerturbationTheoryTerm.change_sort(diffs)
        diffs = nput.vector_take(diffs, permutations, shared=1)[:, :width]
        patterns = np.zeros((len(signatures), width), dtype=diffs.dtype)
        for row, rule in enumerate(signatures):
            patterns[row, :len(rule)] = rule
        if len(patterns):
            matches, _ = nput.find(patterns, diffs, missing_val=len(patterns))
            keep[positions] = matches < len(patterns)
        else:
            keep[positions] = False
        return keep

    _candidate_batch_size = 32

    def candidate_spaces(self, input_states, energy_window, orders, *,
                         block_size=None, max_bytes=64*1024**2, filter_at_generation=False):
        """Apply the native vectorized endpoint screen to bounded dense blocks.

        filter_at_generation opts into the existing BasisStateSpaceFilter
        interface to carry exact signed-rule progress during generation. The
        default keeps the separately profiled vectorized dense-batch path.

        The iterator bounds unfiltered working storage; only accepted rows
        enter the existing sparse results. Retained result storage still grows
        with accepted endpoints. max_bytes bounds the emitted dense block and
        tags; selection/sorting workspaces add bounded overhead proportional
        to that block, rather than to the complete unfiltered search.
        """
        freqs = np.asarray(self.eval.freqs)
        sorting = np.argsort(freqs)
        inverse = np.argsort(sorting)
        sources = input_states.excitations
        energies = np.dot(sources, freqs)
        result = [{} for _ in sources]
        for a in orders:
            quanta = [self.endpoint_quanta(n, a) for n in sources]
            active = [i for i in np.argsort(energies, kind='stable') if quanta[i]]
            rules = self._wavefunction_rules(a)
            max_change = max((abs(c) for rule in rules for c in rule), default=0)
            for start in range(0, len(active), self._candidate_batch_size):
                rows = active[start:start+self._candidate_batch_size]
                windows = [[max(float(energies[i])-energy_window, 0),
                            float(energies[i])+energy_window] for i in rows]
                max_state = np.max(sources[rows], axis=0)+max_change
                pieces = [[] for _ in rows]
                generation_filter = None
                if filter_at_generation:
                    filter_sources = BasisStateSpace(input_states.basis,sources[rows][:,sorting])
                    by_quanta = {}
                    for position,i in enumerate(rows):
                        by_quanta.setdefault((int(np.sum(sources[i])),tuple(sorted(quanta[i]))),[]).append(position)
                    prefilters = [(BasisStateSpace(input_states.basis,
                                      filter_sources.excitations[positions]),
                                   [rule for rule in rules if len(rule)<=len(freqs)
                                    and sum(rule) in {q-source_quanta for q in allowed}])
                                  for (source_quanta,allowed),positions in by_quanta.items()]
                    generation_filter = BasisStateSpaceFilter(filter_sources,prefilters,None)
                for tags, block in BasisStateSpace.states_in_windows(
                        freqs[sorting], windows, max_state=max_state[sorting],
                        min_quantas=[min(quanta[i]) for i in rows],
                        max_quantas=[max(quanta[i]) for i in rows],
                        block_size=block_size, max_bytes=max_bytes, sparse=filter_at_generation,
                        sparse_traversal=filter_at_generation, return_iterator=True,
                        **({'filter':generation_filter} if generation_filter is not None else {})):
                    for window in np.unique(tags):
                        i = rows[window]
                        if filter_at_generation:
                            pieces[int(window)].append(nput.SparseArray.from_data(
                                block.data[:,inverse]))
                            continue
                        candidates = block[tags == window][:,inverse]
                        keep = self.allowed_endpoints(sources[i], candidates, a)
                        if np.any(keep):
                            accepted = candidates[keep]
                            dtype = np.promote_types(BasisStateSpace.excitations_dtype,
                                nput.infer_int_dtype(int(np.max(accepted))))
                            pieces[int(window)].append(nput.SparseArray.from_data(
                                accepted.astype(dtype, copy=False), layout='csr'))
                for i, chunks in zip(rows, pieces):
                    sparse = (chunks[0].concatenate(*chunks[1:],axis=0) if len(chunks)>1
                              else chunks[0] if chunks else
                              nput.SparseArray.empty((0,len(freqs)),dtype=int))
                    chunks.clear()
                    result[i][a] = BasisStateSpace(input_states.basis, sparse)
        return result

    def get_test_wfn_corrs(self, input_states, energy_window, order=None,
                           degenerate_states=None, target_orders=None,
                           state_specific_degeneracies=True):
        from .DegeneracySpecs import DegeneracySpec
        if self.target != 'intensities' or not self.settings.get('strong_coupling_property_filter', True):
            return self.runner.get_test_wfn_corrs(input_states, energy_window,
                order=order, degenerate_states=degenerate_states,
                target_orders=target_orders,
                state_specific_degeneracies=state_specific_degeneracies)
        if order is None:
            order = self.order
        orders = tuple(range(1, order+1) if target_orders is None else target_orders)
        if not orders or min(orders)<1 or max(orders)>self.order:
            raise ValueError('discovery orders must be positive and within the property order')
        candidates = self.candidate_spaces(input_states, energy_window, orders)
        if degenerate_states is not None and len(degenerate_states):
            if state_specific_degeneracies and not isinstance(degenerate_states, dict):
                degenerate_states = DegeneracySpec.state_specific_pairs(
                    DegeneracySpec.merge_state_blocks(degenerate_states))
        else:
            degenerate_states = None
        # Evaluate only the retained endpoints at their requested order. Normal
        # analytic compilation, dispatch, projection and numerical defaults apply.
        evaluated = {}
        for a in orders:
            pairs = [[source, spaces[a].excitations]
                for source,spaces in zip(input_states.excitations, candidates)
                if a in spaces and len(spaces[a])]
            if not pairs:
                continue
            corrections = self.eval.get_wavefunction_corrections(pairs, order=order,
                target_orders=(a,), degenerate_states=degenerate_states)
            for source, finals, values in zip(corrections.initial_states.excitations,
                    corrections.final_states, corrections.corrections):
                evaluated.setdefault(tuple(source), {})[a] = (finals, values[a])
        initial = []
        finals = []
        values = []
        for source in input_states.excitations:
            rows = evaluated.get(tuple(source))
            if not rows:
                continue
            space = next(iter(rows.values()))[0]
            for other, _ in rows.values():
                space = space.union(other)
            block = np.zeros((order+1, len(space)))
            for a, (other, row) in rows.items():
                block[a, space.find(other)] = row
            initial.append(source)
            finals.append(space)
            values.append(block)
        if not initial:
            return None
        return BasicAPTCorrections(BasisStateSpace(input_states.basis, initial), finals, values)


def resolve_discovery_space(runner, states, degeneracy_specs):
    from .Runner import MultiVPTStateSpace
    if isinstance(states, MultiVPTStateSpace) or not uses_automatic_discovery(degeneracy_specs):
        return states
    system = runner._dummy_system(runner)
    packet = MultiVPTStateSpace(states, system, frequencies=runner.eval.freqs, degeneracy_specs=None)
    view = AnalyticPropertyDiscovery(runner, packet, getattr(runner, '_property_discovery_settings', None))
    return MultiVPTStateSpace(states, system, frequencies=runner.eval.freqs,
        degeneracy_specs=degeneracy_specs, evaluator=view)
