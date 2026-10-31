from pathlib import Path
import tempfile,unittest
import numpy as np
from McUtils.Scaffolding import NullLogger,NullCheckpointer
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace
from Psience.VPT2.Solver import PerturbationTheorySolver

class SolverTests(unittest.TestCase):
    def run_case(self, method, external=False, disk=None):
        basis=HarmonicOscillatorProductBasis(2)
        targets=BasisStateSpace(basis,[[0,0],[1,0],[0,1]])
        def rep(word,c):return basis.representation(*word,coeffs=c,contraction_method=method)
        h0=rep(('x','x'),np.diag([.04,.063])/2)+rep(('p','p'),np.diag([.04,.063])/2)
        c3=np.array([[[.0001,-.00003],[-.00003,.00002]],[[-.00003,.00002],[.00002,-.00004]]])
        h1=rep(('x','x','x'),c3/6)
        h2=rep(('x','x','x','x'),np.ones((2,)*4)*.00001/24)
        s=PerturbationTheorySolver([h0,h1,h2],targets,order=2,logger=NullLogger(),
            checkpointer=NullCheckpointer(),memory_constrained=True,
            state_space_storage='external' if external else 'memory',hamiltonian_storage=disk)
        s.identify_strong_couplings=lambda result: ({},None)
        result=s.apply_VPT()
        return s.flat_total_space.indices,result

    def test_second_order_solver_and_disk(self):
        indices,baseline=self.run_case('tensor')
        for external in [False,True]:
            with tempfile.TemporaryDirectory() as directory:
                ids,result=self.run_case('direct',external,directory)
                np.testing.assert_array_equal(indices,ids)
                np.testing.assert_allclose(baseline.energy_corrs,result.energy_corrs,atol=1e-14)
                for a,b in zip(baseline.wfn_corrections,result.wfn_corrections):
                    np.testing.assert_allclose(a.asarray(),b.asarray(),atol=3e-13)

if __name__=='__main__':unittest.main()
