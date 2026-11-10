import types
import unittest

import numpy as np
import Psience.VPT2.Analytic as Analytic


class AnalyticDAGConstructionTests(unittest.TestCase):

    def test_long_path_term_sums_remain_shallow_and_evaluable(self):
        primitive = Analytic.OperatorExpansionTerm(
            [['x']], order=0, identities=[0], logger=False
        )
        term = Analytic.PerturbationTheoryTermSum(*([primitive] * 1500))
        poly = term.get_core_poly((1,))
        self._check_long_path_sum(primitive, poly)

    def test_long_path_product_sums_remain_shallow_and_evaluable(self):
        primitive = Analytic.OperatorExpansionTerm(
            [['x']], order=0, identities=[0], logger=False
        )
        poly = primitive.get_poly_terms((1,))
        product = types.SimpleNamespace(
            changes={(1,): [((), (), [], [], [], [], [])] * 1500},
            polynomial_representation='path', gen1=primitive, gen2=primitive,
            get_poly_product_terms=lambda *args: poly
        )
        result = Analytic.PerturbationTheoryTermProduct.get_core_poly(product, (1,))
        self._check_long_path_sum(primitive, result)

    def _check_long_path_sum(self, primitive, poly):
        # Inspect without recursively walking the very tree under test.
        depths = {}
        pending = [(poly.poly_obj, False)]
        while pending:
            node, ready = pending.pop()
            if id(node) in depths:
                continue
            children = [arg for arg in node._node[1:]
                        if isinstance(arg, Analytic.PTTensorCoeffProductDAG)]
            if ready:
                depths[id(node)] = 1 + max(
                    (depths[id(child)] for child in children), default=0
                )
            else:
                pending.append((node, True))
                pending.extend((child, False) for child in children)
        self.assertEqual(poly.poly_obj.operator_keys, frozenset({(0, 0)}))
        self.assertLess(depths[id(poly.poly_obj)], 30)
        # Pruning and streaming must remain safe too; key lookup alone is
        # insufficient to protect a linearly deep tensor DAG.
        pruned = poly.poly_obj.prune_operators([(9, 9)])
        self.assertEqual(pruned.operator_keys, poly.poly_obj.operator_keys)
        evaluator = Analytic.PerturbationTheoryExpressionEvaluator(
            primitive, poly, (1,)
        )
        result = evaluator.evaluate(
            [[np.array([0]), np.array([[0]])]],
            [[np.array([.1])]], np.array([.04]),
            evaluation_mode='dag', evaluation_backend='indexed',
            zero_cutoff=1e-30
        )[0][0]
        self.assertAlmostEqual(result, 1500 * .1 / np.sqrt(2), places=10)


if __name__ == '__main__':
    unittest.main()
