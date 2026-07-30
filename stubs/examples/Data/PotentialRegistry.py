"""Extracted from DataTests.test_PotentialRegistry via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DataTests.test_PotentialRegistry"""

from Peeves.TestUtils import *
from Psience.Data import *
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Plots import *
import McUtils.Numputils as nput
from unittest import TestCase
import sys, h5py, math, numpy as np

class DataTests(TestCase):
    maxDiff = None

    @debugTest
    def test_PotentialRegistry(self):
        from Psience.Molecools import Molecule
        h2co_mod = PotentialRegistryAPI().get_potential('H2COPot')

        def cart_pot(coords, order=None):
            return h2co_mod.Potential.get_pot(coords)

        def internal_pot(coords, order=None):
            coords = coords[..., (0, 1, 3, 2, 4, 5)]
            vals = h2co_mod.InternalsPotential.get_pot(coords, threading_mode='serial')
            return vals
        ochh_base = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'), energy_evaluator={'potential_function': cart_pot, 'permutation': [2, 3, 1, 0], 'distance_units': 'Angstroms', 'energy_units': 'Wavenumbers'})
        cart_eng = ochh_base.calculate_energy()
        ochh_base = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'), energy_evaluator={'potential_function': internal_pot, 'distance_units': 'Angstroms', 'energy_units': 'Wavenumbers', 'strip_embedding': True}, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 1, 0, -1], [3, 1, 0, 2]])
        opt_ochh = ochh_base.optimize(method='conjugate-gradient', stencil=3, prevent_oscillations=3, restart_interval=15)
        b1 = ochh_base.calculate_energy()
        b2 = opt_ochh.calculate_energy()
        print(b1, b2)
