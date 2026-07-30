"""Extracted from VPT2Tests.test_BenzeneInternals via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_BenzeneInternals"""

import itertools
import tempfile
try:
    from Peeves.TestUtils import *
    from Peeves import BlockProfiler
except:
    pass
from unittest import TestCase
from Psience.VPT2 import *
from Psience.Molecools import Molecule
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace
from McUtils.Data import UnitsData
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.Scaffolding import *
from McUtils.Profilers import Timer
from McUtils.Parallelizers import Parallelizer, SerialNonParallelizer, MultiprocessingParallelizer
from McUtils.Zachary import FiniteDifferenceDerivative
import sys, os, numpy as np, itertools as ip

class VPT2Tests(TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(100000000.0))

    @validationTest
    def test_BenzeneInternals(self):
        import McUtils.Coordinerds as coordops
        zmat = coordops.parse_zmatrix_string('\n1          \n2   1\t\n3  \t2\t\t1\t\n4  \t2\t\t3\t\t1\t\n5  \t3\t\t2\t\t1\t\n6  \t3\t\t5\t\t2\t\n7  \t5\t\t3\t\t2\t\n8  \t5\t\t7\t\t3\t\n9  \t7\t\t5\t\t3\t\n10 \t7\t\t9\t\t5\t\n11 \t9\t\t1\t\t7\t\n12 \t1\t\t2\t\t9', has_values=False, atoms_are_order=True)[1]
        state = VPTStateMaker(30)
        runner, _ = VPTRunner.construct(TestManager.test_data('benzene.fchk'), [state()] + [state(i) for i in range(1, 7)], logger=True, degeneracy_specs='auto', mixed_derivative_handling_mode='old')
        import McUtils.Formatters as mfmt
        v4, v3, v0 = (runner.hamiltonian.V_terms[2], runner.hamiltonian.V_terms[1], runner.hamiltonian.V_terms[0])
        print(mfmt.format_symmetric_tensor_elements(v0 * UnitsData.hartrees_to_wavenumbers, symmetries=[(0, 1)], cutoff=0.01))
        print(mfmt.format_symmetric_tensor_elements(v3 * UnitsData.hartrees_to_wavenumbers, symmetries=[(1, 2)], cutoff=0.01))
        print(mfmt.format_symmetric_tensor_elements(v4 * UnitsData.hartrees_to_wavenumbers, symmetries=[(0, 1), (2, 3)], cutoff=0.01))
