import types
import unittest

import numpy as np
from McUtils.Scaffolding import NullLogger
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace,BasisMultiStateSpace
from Psience.VPT2.Wavefunctions import PerturbationTheoryWavefunctions

class VPT2DipolePairingTests(unittest.TestCase):
    def test_only_required_pairs_for_single_and_multiple_initial_states(self):
        basis=HarmonicOscillatorProductBasis(3)
        def space(states):return BasisStateSpace(basis,states)
        def nested(items):
            a=np.empty(len(items),dtype=object)
            for i,x in enumerate(items):a[i]=x
            return BasisMultiStateSpace(a)
        ground=[space([[0,0,0]]),space([[1,0,0],[0,1,0],[0,0,1]]),space([[2,0,0],[1,1,0],[0,2,0]])]
        extra=[space([[0,0,2]]),space([[1,0,2],[0,1,2],[0,0,3]]),space([[2,0,2],[1,1,2],[0,2,2]])]
        target=[space([[1,0,0],[0,1,0]]),space([[0,0,0],[2,0,0],[0,2,0],[1,1,0],[0,0,2],[1,0,1]]),space([[1,0,0],[0,1,0],[0,0,1],[3,0,0],[0,3,0],[2,1,0]])]
        obj=types.SimpleNamespace(logger=NullLogger())
        def keys(pairs):return {tuple(sorted([int(a),int(b)])) for a,b in zip(pairs.bras.indices,pairs.kets.indices)}
        for bra_orders in [[ground],[ground,extra]]:
            for degree in [1,2,3]:
                k=degree-1
                op=basis.representation(*(['x']*degree))
                got=PerturbationTheoryWavefunctions._get_rep_states(obj,op,k,3,nested([nested(target)]),nested([nested(b) for b in bra_orders]))
                expected=set()
                for bras in bra_orders:
                    for i in range(3):
                        for j in range(3):
                            if i+j+k<=2:
                                p=bras[i].get_representation_brakets(other=target[j],selection_rules=op.selection_rules)
                                expected.update(keys(p))
                self.assertEqual(keys(got),expected)

if __name__=='__main__':unittest.main(verbosity=2)
