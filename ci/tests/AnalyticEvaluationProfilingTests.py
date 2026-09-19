"""Manual pyinstrument profiles for the materialized analytic evaluator.

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
"""

import os
import pathlib
import time
import unittest

import numpy as np

try:
    from pyinstrument import Profiler
except ImportError:
    Profiler = None

try:
    import Psience.VPT2.Analytic as Analytic
except ModuleNotFoundError:
    import Psience.Psience.VPT2.Analytic as Analytic


RUN_PROFILES = os.environ.get('PSIENCE_RUN_ANALYTIC_PROFILES', '').lower() in {
    '1', 'true', 'yes', 'on'
}


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
    def _profile(label, function):
        profiler = Profiler(interval=.001, async_mode='disabled')
        start = time.perf_counter()
        profiler.start()
        try:
            result = function()
        finally:
            profiler.stop()
        elapsed = time.perf_counter() - start
        report = profiler.output_text(
            unicode=True,
            color=False,
            show_all=os.environ.get('PSIENCE_ANALYTIC_PROFILE_SHOW_ALL') == '1',
            timeline=False
        )
        print('\n{}: {:.3f}s\n{}'.format(label, elapsed, report))

        output_dir = os.environ.get('PSIENCE_ANALYTIC_PROFILE_OUTPUT')
        if output_dir:
            output_dir = pathlib.Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / '{}.txt'.format(label)).write_text(report)
            profiler.write_html(output_dir / '{}.html'.format(label))
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

        def evaluate_all():
            values = {}
            for name, evaluator in expressions.items():
                values[name] = evaluator.evaluate(
                    state_perms,
                    expansion,
                    frequencies,
                    degenerate_changes=degenerate_changes,
                    only_degenerate_terms=modes[name],
                    evaluation_mode='materialized'
                )
            return values

        values, report = self._profile(
            'analytic_materialized_degenerate_evaluation', evaluate_all
        )
        self.assertEqual(set(values), set(expressions))
        self.assertIn('Analytic.py', report)


if __name__ == '__main__':
    unittest.main()
