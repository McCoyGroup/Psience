"""Manual pyinstrument profiles for the analytic evaluator.

These are opt-in because their purpose is to retain a reproducible performance
workload, not to impose a wall-time threshold on CI.  Run with, for example::

    PSIENCE_RUN_ANALYTIC_PROFILES=1 \
    PSIENCE_ANALYTIC_PROFILE_STATES=12 \
    PSIENCE_ANALYTIC_PROFILE_PERMS=24 \
    python -m unittest ci.tests.AnalyticEvaluationProfilingTests -v

Set ``PSIENCE_ANALYTIC_PROFILE_OUTPUT`` to a directory to retain text and HTML
profiles.  The default only prints the text reports.  Use
``PSIENCE_ANALYTIC_PROFILE_COMPONENTS=main`` to isolate the ordinary coupling
term, or a comma-separated selection of ``main,left,right,both``.  A shortened
resonance list can be selected with ``PSIENCE_ANALYTIC_PROFILE_CHANGE_LIMIT``.
``PSIENCE_ANALYTIC_PROFILE_EVALUATION_MODES`` selects the comma-separated
evaluator variants; it defaults to comparing the legacy and bounded-PolyPath-
cache materialized implementations.  ``dag_legacy,dag`` compares the DAG
streamer's original single-level cache with the shared PolyPath cache.

``test_profile_tbhp_sized_evaluator`` is the molecule-sized workload.  The
available TBHP checkpoint lacks cubic force derivatives, so the test uses its
actual normal-mode frequencies with deterministic synthetic VPT coefficients.
Set
``PSIENCE_ANALYTIC_PROFILE_SECONDS`` to bound each mode (default: five minutes)
and ``PSIENCE_ANALYTIC_PROFILE_TBHP_STATES`` to control its excitation depth.
``PSIENCE_ANALYTIC_PROFILE_INTERVAL`` controls the pyinstrument sampling period
(default: 1 ms); use 10--20 ms for long profiles to bound profiler storage.
An interrupted profile is still written and reported as an intentional sample.
"""

import os
import pathlib
import signal
import time
import unittest

import numpy as np

try:
    from pyinstrument import Profiler
except ImportError:
    Profiler = None

try:
    import Psience.VPT2.Analytic as Analytic
    from Psience.Molecools import Molecule
except ModuleNotFoundError:
    import Psience.Psience.VPT2.Analytic as Analytic
    from Psience.Psience.Molecools import Molecule


RUN_PROFILES = os.environ.get('PSIENCE_RUN_ANALYTIC_PROFILES', '').lower() in {
    '1', 'true', 'yes', 'on'
}


class _ProfileTimeout(RuntimeError):
    pass


@unittest.skipUnless(RUN_PROFILES, 'manual pyinstrument profiling workload')
@unittest.skipIf(Profiler is None, 'pyinstrument is not installed')
class AnalyticMaterializedEvaluationProfilingTests(unittest.TestCase):

    _expression_cache = None

    target_change = (1, -1)
    degenerate_change_classes = (
        (1, -1, -1, -1),
        (-2, 1, 1, 1, -1, -1),
        (-2, 1, -1),
        (-2, 1, 1),
        (2, 1, 1, -1, -1, -1, -1),
        (2, 1, 1, -1, -1, -1),
        (-2, -2, 1, 1, 1),
        (2, -2, 1, 1, -1, -1),
        (1, -1, -1),
        (2, 2, -2, -1, -1),
        (2, 2, -1, -1, -1),
        (2, -1, -1),
        (2, -2, -2),
        (2, 1, -1, -1, -1),
        (1, 1, 1, -1, -1, -1),
        (-2, 1, 1, -1, -1),
        (-2, 1, 1, -1),
        (-2, 1),
        (-2, 1, 1, 1),
        (1, 1, -1, -1, -1),
        (2, -1, -1, -1, -1),
        (2, 2, -2),
        (2, 2, -1, -1, -1, -1),
        (1, 1, 1, -1),
        (2, -2, 1, -1),
        (1, 1, 1, 1, -1, -1, -1),
        (1, 1, 1, -1, -1, -1, -1),
        (-2, -2, 1, 1, 1, 1),
        (2, -2),
        (2, -1),
        (2, -2, 1, 1),
        (2, -2, -2, 1, 1),
        (2, 1, -1, -1),
        (1, 1, 1, -1, -1),
        (2, 2, -1, -1),
        (2, 2, -2, -2),
        (2, 1, -1),
        (1, 1, -1, -1),
        (-2, 1, 1, 1, -1),
        (-2, 1, 1, 1, 1),
        (1, -1),
        (-2, -2, 1, 1),
        (1, 1, -1, -1, -1, -1),
        (-2, 1, 1, 1, 1, -1, -1),
        (2, -1, -1, -1),
        (2, -2, -1, -1),
        (1, 1, -1),
        (2, 1, 1, -1, -1),
        (1, 1, 1, 1, -1, -1),
        (1, 1, 1, 1, -1, -1, -1, -1)
    )

    @classmethod
    def _profile_size(cls):
        nmodes = int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_MODES', 8))
        nmodes = max(nmodes, max(map(len, cls.degenerate_change_classes)))
        return (
            nmodes,
            int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_STATES', 4)),
            int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_PERMS', 8))
        )

    @classmethod
    def _profile_change_classes(cls):
        limit = int(os.environ.get(
            'PSIENCE_ANALYTIC_PROFILE_CHANGE_LIMIT',
            len(cls.degenerate_change_classes)
        ))
        return cls.degenerate_change_classes[:limit]

    @classmethod
    def _profile_components(cls):
        components = os.environ.get(
            'PSIENCE_ANALYTIC_PROFILE_COMPONENTS',
            'main,left,right,both'
        )
        return tuple(component.strip() for component in components.split(',') if component.strip())

    @classmethod
    def _profile_evaluation_modes(cls):
        modes = os.environ.get(
            'PSIENCE_ANALYTIC_PROFILE_EVALUATION_MODES',
            'materialized_legacy,materialized'
        )
        return tuple(mode.strip() for mode in modes.split(',') if mode.strip())

    @classmethod
    def _build_expressions(cls, force=False):
        if cls._expression_cache is not None and not force:
            return cls._expression_cache
        Analytic.AnalyticPerturbationTheorySolver.clear_caches()
        solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
            4, polynomial_representation='path'
        )
        components = cls._profile_components()
        expressions = {}
        if 'main' in components:
            expressions['main'] = solver.reexpressed_hamiltonian(2)(cls.target_change)
        if any(component in components for component in ('left', 'right', 'both')):
            degenerate = solver.reexpressed_hamiltonian_degenerate_correction(
                2, cls._profile_change_classes()
            )
            if 'left' in components:
                expressions['left'] = degenerate.left(cls.target_change)
            if 'right' in components:
                expressions['right'] = degenerate.right(cls.target_change)
            if 'both' in components:
                expressions['both'] = degenerate.both(cls.target_change)
        cls._expression_cache = expressions
        return cls._expression_cache

    @staticmethod
    def _materialize(expressions):
        materialized = {}
        for name, evaluator in expressions.items():
            expr = evaluator.expr
            if Analytic.nput.is_numeric(expr):
                pass
            elif isinstance(expr, Analytic.SqrtChangePoly):
                expr.poly_obj.to_eager()
            materialized[name] = evaluator
        return materialized

    @staticmethod
    def _leaf_coefficient_keys(expression):
        if Analytic.nput.is_numeric(expression):
            return set()
        if isinstance(expression, Analytic.SqrtChangePoly):
            expression = expression.poly_obj
        if not isinstance(expression, Analytic.PTTensorCoeffProductDAG):
            return {
                coefficient
                for product in expression.terms
                for coefficient in product
            }

        keys = set()
        seen = set()
        stack = [expression]
        while stack:
            node = stack.pop()
            if id(node) in seen:
                continue
            seen.add(id(node))
            kind, *args = node._node
            if kind == 'leaf':
                keys.update(
                    coefficient
                    for product in args[0].terms
                    for coefficient in product
                )
            else:
                stack.extend(
                    arg for arg in args
                    if isinstance(arg, Analytic.PTTensorCoeffProductDAG)
                )
        return keys

    @classmethod
    def _build_numerical_data(cls, expressions):
        nmodes, nstates, nperms = cls._profile_size()
        coefficient_keys = set().union(*(
            cls._leaf_coefficient_keys(evaluator.expr)
            for evaluator in expressions.values()
        ))
        rng = np.random.default_rng(81723)
        expansion = []
        if len(coefficient_keys) == 0:
            expansion = [[0]]
        else:
            max_order = max(key[0] for key in coefficient_keys)
            for order in range(max_order + 1):
                types = [key[1] for key in coefficient_keys if key[0] == order]
                row = []
                for coefficient_type in range(max(types, default=0) + 1):
                    ranks = [
                        len(key) - 2
                        for key in coefficient_keys
                        if key[:2] == (order, coefficient_type)
                    ]
                    row.append(
                        0 if len(ranks) == 0 else
                        rng.normal(scale=.05, size=(nmodes,) * max(ranks))
                    )
                expansion.append(row or [0])

        permutations = []
        seen = set()
        while len(permutations) < nperms:
            permutation = tuple(rng.permutation(nmodes))
            if permutation not in seen:
                seen.add(permutation)
                permutations.append(permutation)
        permutations = np.array(permutations, dtype=int)
        state_perms = [
            [rng.integers(4, 9, size=nmodes), permutations]
            for _ in range(nstates)
        ]
        frequencies = np.linspace(.73, 2.41, nmodes) ** 1.17
        degenerate_changes = [
            [np.arange(len(change), dtype=int), np.array(change, dtype=int)]
            for change in cls._profile_change_classes()
        ]
        return state_perms, expansion, frequencies, degenerate_changes

    @staticmethod
    def _profile_sample(label, function, max_seconds=None):
        sample_interval = float(os.environ.get(
            'PSIENCE_ANALYTIC_PROFILE_INTERVAL', .001
        ))
        profiler = Profiler(interval=sample_interval, async_mode='disabled')
        previous_handler = None
        timed_out = False
        result = None

        if max_seconds is not None:
            def timeout_handler(signum, frame):
                raise _ProfileTimeout(
                    '{} exceeded {:.1f}s profiling window'.format(label, max_seconds)
                )

            previous_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.setitimer(signal.ITIMER_REAL, max_seconds)

        start = time.perf_counter()
        profiler.start()
        try:
            result = function()
        except _ProfileTimeout:
            timed_out = True
        finally:
            profiler.stop()
            if max_seconds is not None:
                signal.setitimer(signal.ITIMER_REAL, 0)
                signal.signal(signal.SIGALRM, previous_handler)
        elapsed = time.perf_counter() - start
        report = profiler.output_text(
            unicode=True,
            color=False,
            show_all=os.environ.get('PSIENCE_ANALYTIC_PROFILE_SHOW_ALL') == '1',
            timeline=False
        )
        print('\n{}: {:.3f}s{}\n{}'.format(
            label,
            elapsed,
            ' (sample timed out)' if timed_out else '',
            report
        ))

        output_dir = os.environ.get('PSIENCE_ANALYTIC_PROFILE_OUTPUT')
        if output_dir:
            output_dir = pathlib.Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / '{}.txt'.format(label)).write_text(report)
            profiler.write_html(output_dir / '{}.html'.format(label))
        return result, report, timed_out

    @classmethod
    def _profile(cls, label, function):
        result, report, _ = cls._profile_sample(label, function)
        return result, report

    def test_profile_materialized_construction(self):
        expressions, report = self._profile(
            'analytic_materialized_construction',
            lambda: self._materialize(self._build_expressions(force=True))
        )
        self.assertEqual(set(expressions), set(self._profile_components()))
        self.assertIn('Analytic.py', report)

    def test_profile_materialized_degenerate_evaluation(self):
        expressions = self._materialize(self._build_expressions())
        state_perms, expansion, frequencies, degenerate_changes = (
            self._build_numerical_data(expressions)
        )
        modes = {
            'main': False,
            'left': [True, False],
            'right': [False, True],
            'both': [True, True]
        }

        def evaluate_all(evaluation_mode):
            values = {}
            for name, evaluator in expressions.items():
                values[name] = evaluator.evaluate(
                    state_perms,
                    expansion,
                    frequencies,
                    degenerate_changes=degenerate_changes,
                    only_degenerate_terms=modes[name],
                    evaluation_mode=evaluation_mode
                )
            return values

        reference = None
        for evaluation_mode in self._profile_evaluation_modes():
            # Do not let exact-result hits from an earlier mode conceal the
            # amount of work performed by the PolyPath cache under test.
            Analytic.PerturbationTheoryExpressionEvaluator._poly_cache = (
                Analytic.PerturbationTheoryExpressionEvaluator.get_cache()
            )
            values, report = self._profile(
                'analytic_{}_degenerate_evaluation'.format(evaluation_mode),
                lambda mode=evaluation_mode: evaluate_all(mode)
            )
            if evaluation_mode in {'materialized', 'materialized_cached'}:
                print(
                    'materialized path cache: {}'.format(
                        Analytic.PerturbationTheoryExpressionEvaluator
                        .get_last_materialized_evaluation_stats()
                    )
                )
            elif evaluation_mode in {'dag', 'dag_legacy'}:
                print(
                    'dag cache: {}'.format(
                        Analytic.PerturbationTheoryExpressionEvaluator
                        .get_last_dag_evaluation_stats()
                    )
                )
            self.assertEqual(set(values), set(expressions))
            self.assertIn('Analytic.py', report)
            if reference is None:
                reference = values
            else:
                for name in reference:
                    np.testing.assert_allclose(
                        values[name], reference[name], rtol=2e-12, atol=2e-12
                    )

    def test_profile_dag_cache(self):
        Analytic.AnalyticPerturbationTheorySolver.clear_caches()
        solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
            4, polynomial_representation='path'
        )
        evaluator = solver.energy_correction(2)([])
        coefficient_keys = self._leaf_coefficient_keys(evaluator.expr)
        nmodes = int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_MODES', 4))
        nstates = int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_STATES', 6))
        nperms = int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_PERMS', 12))
        chunk_size = int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_DAG_CHUNK', 32))
        rng = np.random.default_rng(61492)

        expansion = []
        for order in range(max(key[0] for key in coefficient_keys) + 1):
            coefficient_types = [
                key[1] for key in coefficient_keys if key[0] == order
            ]
            order_expansion = []
            for coefficient_type in range(max(coefficient_types, default=0) + 1):
                ranks = [
                    len(key) - 2 for key in coefficient_keys
                    if key[:2] == (order, coefficient_type)
                ]
                order_expansion.append(
                    0 if len(ranks) == 0 else
                    rng.normal(scale=.05, size=(nmodes,) * max(ranks))
                )
            expansion.append(order_expansion or [0])

        permutations = []
        seen = set()
        while len(permutations) < nperms:
            permutation = tuple(rng.permutation(nmodes))
            if permutation not in seen:
                seen.add(permutation)
                permutations.append(permutation)
        permutations = np.array(permutations, dtype=int)
        state_perms = [
            [rng.integers(0, 7, size=nmodes), permutations]
            for _ in range(nstates)
        ]
        frequencies = np.linspace(.7, 2.3, nmodes)

        reference = None
        for evaluation_mode in self._profile_evaluation_modes():
            if evaluation_mode not in {'dag', 'dag_legacy'}:
                continue
            Analytic.PerturbationTheoryExpressionEvaluator._poly_cache = (
                Analytic.PerturbationTheoryExpressionEvaluator.get_cache()
            )
            values, report = self._profile(
                'analytic_{}_cache'.format(evaluation_mode),
                lambda mode=evaluation_mode: evaluator.evaluate(
                    state_perms,
                    expansion,
                    frequencies,
                    evaluation_mode=mode,
                    dag_chunk_size=chunk_size
                )
            )
            print(
                'dag cache: {}'.format(
                    Analytic.PerturbationTheoryExpressionEvaluator
                    .get_last_dag_evaluation_stats()
                )
            )
            self.assertIn('Analytic.py', report)
            if reference is None:
                reference = values
            else:
                for actual, expected in zip(values, reference):
                    np.testing.assert_allclose(
                        actual, expected, rtol=2e-12, atol=2e-12
                    )

    def test_profile_tbhp_sized_evaluator(self):
        data_file = pathlib.Path(__file__).with_name('TestData') / 'tbhp_180.fchk'
        molecule = Molecule.from_file(str(data_file))
        frequencies = molecule.get_normal_modes().freqs
        nmodes = min(
            len(frequencies),
            int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_TBHP_MODES', len(frequencies)))
        )
        frequencies = frequencies[:nmodes]
        nstates = int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_TBHP_STATES', 4))
        nperms = int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_TBHP_PERMS', 8))
        chunk_size = int(os.environ.get('PSIENCE_ANALYTIC_PROFILE_DAG_CHUNK', 256))
        max_seconds = float(os.environ.get(
            'PSIENCE_ANALYTIC_PROFILE_SECONDS', 300
        ))
        Analytic.AnalyticPerturbationTheorySolver.clear_caches()
        solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
            4, polynomial_representation='path'
        )
        evaluator = solver.energy_correction(2)([])
        coefficient_keys = self._leaf_coefficient_keys(evaluator.expr)
        rng = np.random.default_rng(81173)
        expansion = []
        for order in range(max(key[0] for key in coefficient_keys) + 1):
            coefficient_types = [
                key[1] for key in coefficient_keys if key[0] == order
            ]
            order_expansion = []
            for coefficient_type in range(max(coefficient_types, default=0) + 1):
                ranks = [
                    len(key) - 2 for key in coefficient_keys
                    if key[:2] == (order, coefficient_type)
                ]
                order_expansion.append(
                    0 if len(ranks) == 0 else
                    rng.normal(scale=.01, size=(nmodes,) * max(ranks))
                )
            expansion.append(order_expansion or [0])

        permutations = np.array([
            rng.permutation(nmodes) for _ in range(nperms)
        ], dtype=int)
        states = np.zeros((nstates, nmodes), dtype=int)
        for state_index in range(1, nstates):
            states[state_index, (state_index - 1) % nmodes] = 1
        state_perms = [[state, permutations] for state in states]

        for evaluation_mode in self._profile_evaluation_modes():
            if evaluation_mode not in {'dag', 'dag_legacy'}:
                continue
            Analytic.PerturbationTheoryExpressionEvaluator._poly_cache = (
                Analytic.PerturbationTheoryExpressionEvaluator.get_cache()
            )
            result, report, timed_out = self._profile_sample(
                'analytic_tbhp_sized_{}'.format(evaluation_mode),
                lambda mode=evaluation_mode: evaluator.evaluate(
                    state_perms,
                    expansion,
                    frequencies,
                    evaluation_mode=mode,
                    dag_chunk_size=chunk_size
                ),
                max_seconds=max_seconds
            )
            self.assertIn('Analytic.py', report)
            self.assertTrue(timed_out or result is not None)


if __name__ == '__main__':
    unittest.main()
