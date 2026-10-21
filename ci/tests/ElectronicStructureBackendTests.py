"""
Tests for the `Psience.ElectronicStructure` backends and the `Psience.Molecools` evaluators that
wrap them.

The interface and the wrapper machinery are tested with an analytic toy backend; xTB is tested
through a stand-in for the `xtb` python bindings (the real ones only come from conda-forge);
PySCF, ASE, and AIMNet2 run for real when installed (AIMNet2 also needs its weights, from
`$MODEL_CACHE_DIR/aimnet` or a download).
"""

import importlib.util
import os
import sys
import types
import unittest

import numpy as np

from McUtils.Data import UnitsData
from Psience.ElectronicStructure import (
    ElectronicStructureBackend, XTBBackend, PySCFBackend, ASECalculatorBackend, MACEPolarBackend
)
from Psience.Molecools.Evaluator import (
    ElectronicStructureBackendEvaluator, EnergyEvaluator, ChargeEvaluator, DipoleEvaluator,
    PySCFEnergyEvaluator, PySCFChargeEvaluator, PySCFDipoleEvaluator, PySCFDipolePolarizabilityEvaluator,
    XTBEnergyEvaluator, XTBChargeEvaluator, XTBDipoleEvaluator, ASECalcEnergyEvaluator, MACEEnergyEvaluator,
    MACEPolarElectrostaticPotentialEvaluator
)

A2B = UnitsData.convert("Angstroms", "BohrRadius")
EV2H = UnitsData.convert("ElectronVolts", "Hartrees")
WATER = ['O', 'H', 'H']
WATER_ANG = np.array([[0.0, 0.0, 0.117], [0.0, 0.757, -0.468], [0.0, -0.757, -0.468]])


def has_module(name):
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


class ToyBackend(ElectronicStructureBackend):
    """
    E = k/2 |x - x0|^2 + 0.1 charge (eV, Angstrom); q_i = charge/n + a (z_i - <z>);
    mu = sum_i q_i x_i (e Angstrom)
    """
    distance_units = 'Angstroms'
    property_units = {'energy': 'ElectronVolts', 'charges': 'ElementaryCharge',
                      'dipole': ('ElementaryCharge', 'Angstroms'), 'polarizability': None}
    analytic_derivative_orders = {'energy': 2, 'charges': 1, 'dipole': 0, 'polarizability': 0}
    property_shapes = {'charges': lambda n: (n,), 'dipole': (3,), 'polarizability': (3, 3)}
    k = 2.
    a = .3

    def __init__(self, atoms, x0=None, charge=None, multiplicity=None, quiet=True):
        super().__init__(atoms, charge=charge, multiplicity=multiplicity, quiet=quiet)
        self.x0 = np.zeros((len(atoms), 3)) if x0 is None else np.asanyarray(x0)
        self.calls = []

    def charges(self, coords, charge):
        z = coords[..., 2]
        return charge / self.natoms + self.a * (z - np.mean(z, axis=-1, keepdims=True))

    def evaluate_properties(self, coords, orders, charge, multiplicity, **opts):
        self.calls.append((len(coords), dict(orders), charge, multiplicity))
        nstruct, n = coords.shape[:2]
        d = (coords - self.x0).reshape(nstruct, -1)
        out = {}
        for p, o in orders.items():
            if p == 'energy':
                out[p] = [self.k / 2 * np.sum(d ** 2, axis=-1) + .1 * charge, self.k * d,
                          np.broadcast_to(self.k * np.eye(3 * n), (nstruct, 3 * n, 3 * n))][:o + 1]
            elif p == 'charges':
                dq = np.zeros((n, 3, n))
                dq[:, 2, :] = self.a * (np.eye(n) - 1 / n)
                out[p] = [self.charges(coords, charge), np.broadcast_to(dq.reshape(3 * n, n), (nstruct, 3 * n, n))][:o + 1]
            elif p == 'dipole':
                out[p] = [np.einsum('si,sij->sj', self.charges(coords, charge), coords)]
            elif p == 'polarizability':
                out[p] = [np.broadcast_to(np.diag([1., 2., 3.]), (nstruct, 3, 3))]
        return out


class ToyEnergyEvaluator(ElectronicStructureBackendEvaluator, EnergyEvaluator):
    backend_type = ToyBackend
    backend_property = 'energy'


class ToyChargeEvaluator(ElectronicStructureBackendEvaluator, ChargeEvaluator):
    backend_type = ToyBackend
    backend_property = 'charges'


class ToyDipoleEvaluator(ElectronicStructureBackendEvaluator, DipoleEvaluator):
    backend_type = ToyBackend
    backend_property = 'dipole'


class BackendInterfaceTests(unittest.TestCase):

    def setUp(self):
        self.x0 = WATER_ANG * .9
        self.backend = ToyBackend(WATER, x0=self.x0)

    def test_Shapes(self):
        b = self.backend
        x = WATER_ANG
        res = b.evaluate(x, order=2, properties='energy')['energy']
        self.assertEqual([r.shape for r in res], [(), (9,), (9, 9)])
        batch = np.stack([x, x * 1.1, x * 1.2, x]).reshape(2, 2, 3, 3)
        res = b.evaluate(batch, order={'energy': 1, 'charges': 1, 'dipole': 0})
        self.assertEqual([r.shape for r in res['energy']], [(2, 2), (2, 2, 9)])
        self.assertEqual([r.shape for r in res['charges']], [(2, 2, 3), (2, 2, 9, 3)])
        self.assertEqual([r.shape for r in res['dipole']], [(2, 2, 3)])
        # flat (..., 3N) coordinates work too, and everything came from one call
        self.assertEqual(b.evaluate(batch.reshape(2, 2, 9), properties='charges')['charges'][0].shape, (2, 2, 3))
        self.assertEqual(b.calls[-2][0], 4)
        self.assertTrue(np.allclose(res['energy'][0][0, 1], self.backend.k / 2 * np.sum((x * 1.1 - self.x0) ** 2)))
        with self.assertRaises(ValueError):
            b.evaluate(x, order=1, properties='dipole')  # no analytic dipole derivatives
        with self.assertRaises(ValueError):
            b.evaluate(x, properties='spin')
        with self.assertRaises(ValueError):
            b.evaluate(np.zeros((4, 3)))

    def test_StateOverrides(self):
        b = self.backend
        q = b.evaluate_property(WATER_ANG, 'charges', charge=1, multiplicity=2)[0]
        self.assertTrue(np.isclose(q.sum(), 1))
        self.assertEqual(b.calls[-1][2:], (1, 2))
        self.assertTrue(np.isclose(b.evaluate_property(WATER_ANG, 'charges')[0].sum(), 0))
        self.assertEqual(b.calls[-1][2:], (0, None))
        with self.assertRaises(ValueError):
            b.check_electron_count(0, 2)  # 10 electrons can't be a doublet
        b.check_electron_count(1, 2)

    def test_FieldResponse(self):
        mu, alpha = self.backend.evaluate_field_response(WATER_ANG, 0, 2)
        self.assertEqual(mu[0].shape, (3,))
        self.assertTrue(np.allclose(alpha[0], np.diag([1., 2., 3.])))

    def test_Registry(self):
        backend, opts = ElectronicStructureBackend.resolve_backend('pyscf:b3lyp/def2-svp')
        self.assertIs(backend, PySCFBackend)
        self.assertEqual(opts, {'level_of_theory': 'b3lyp', 'basis_set': 'def2-svp'})
        self.assertIs(ElectronicStructureBackend.resolve_backend('mace-polar')[0], MACEPolarBackend)
        self.assertEqual(ElectronicStructureBackend.resolve_backend('xtb:GFN1-xTB')[1], {'method': 'GFN1-xTB'})
        with self.assertRaises(ValueError):
            ElectronicStructureBackend.resolve_backend('not-a-model')
        names, takes_settings = ASECalculatorBackend.get_init_options()
        self.assertTrue(takes_settings)
        self.assertTrue({'charge', 'multiplicity', 'calculator'} <= names)
        names, takes_settings = PySCFBackend.get_init_options()
        self.assertFalse(takes_settings)
        self.assertIn('level_of_theory', names)  # keyword-only arguments count


class BackendEvaluatorTests(unittest.TestCase):

    def setUp(self):
        self.x0 = WATER_ANG * .9

    def test_EnergyWrapper(self):
        ev = ToyEnergyEvaluator(WATER, x0=self.x0)
        self.assertEqual(ev.analytic_derivative_order, 2)
        self.assertEqual(ev.property_units, 'ElectronVolts')
        self.assertEqual(ev.distance_units, 'Angstroms')
        self.assertEqual(ev.atoms, WATER)
        E, G, H, T = ev.evaluate(WATER_ANG, order=3)
        k = ToyBackend.k
        d = (WATER_ANG - self.x0).reshape(-1)
        # unit handling stays in the evaluator: Hartree and Bohr
        self.assertTrue(np.isclose(E, k / 2 * d @ d * EV2H))
        self.assertTrue(np.allclose(G, k * d * EV2H / A2B))
        self.assertTrue(np.allclose(H, k * np.eye(9) * EV2H / A2B ** 2))
        # third derivatives come from finite differences of the analytic Hessians
        self.assertEqual(T.shape, (9, 9, 9))
        self.assertTrue(np.allclose(T, 0, atol=1e-8))
        # a lower analytic order falls back to finite differences of gradients
        H_fd = ev.evaluate(WATER_ANG, order=2, analytic_derivative_order=1)[2]
        self.assertTrue(np.allclose(H_fd, H, atol=1e-6))

    def test_OptionSplitting(self):
        ev = ToyEnergyEvaluator(WATER, x0=self.x0, charge=1, mesh_spacing=.002, embedding=None)
        self.assertEqual(ev.charge, 1)
        self.assertEqual(ev.defaults, {'mesh_spacing': .002})  # evaluation defaults stay on the evaluator
        shared = ToyBackend(WATER, x0=self.x0)
        a, b = ToyEnergyEvaluator(backend=shared), ToyChargeEvaluator(backend=shared)
        self.assertIs(a.backend, b.backend)
        self.assertEqual(b.analytic_derivative_order, 1)
        # positional backend arguments follow the backend's signature
        self.assertIs(ToyEnergyEvaluator(WATER, self.x0).backend.x0, self.x0)
        with self.assertRaises(ValueError):
            ToyEnergyEvaluator(WATER, self.x0, backend=shared)

    def test_ChargeAndDipoleWrappers(self):
        q = ToyChargeEvaluator(WATER, x0=self.x0)
        q0, q1, q2 = q.evaluate(WATER_ANG, order=2)
        self.assertEqual((q0.shape, q1.shape, q2.shape), ((3,), (9, 3), (9, 9, 3)))
        self.assertTrue(np.allclose(q1[2::3], ToyBackend.a * (np.eye(3) - 1 / 3) / A2B))
        self.assertTrue(np.allclose(q2, 0, atol=1e-6))
        mu = ToyDipoleEvaluator(WATER, x0=self.x0)
        m0, m1 = mu.evaluate(WATER_ANG, order=1)  # e Bohr, derivatives by finite differences
        self.assertTrue(np.allclose(m0, q0 @ WATER_ANG * A2B))
        self.assertEqual(m1.shape, (9, 3))


class XTBStandInTests(unittest.TestCase):
    """the xTB backend's plumbing, against a stand-in for `xtb.interface`/`xtb.utils`"""

    def setUp(self):
        self.calcs = []
        tests = self

        class Result:
            def __init__(self, calc):
                self.calc = calc
            def get_energy(self):
                return -float(np.sum(self.calc.coords ** 2)) - self.calc.charge
            def get_gradient(self):
                return -2 * self.calc.coords
            def get_charges(self):
                return np.full(len(self.calc.numbers), self.calc.charge / len(self.calc.numbers))
            def get_dipole(self):
                return self.get_charges() @ self.calc.coords

        class XTBException(Exception):
            pass

        class Calculator:
            fail_first = False
            def __init__(self, param, numbers, coords, charge=0, uhf=0):
                self.param, self.numbers, self.coords, self.charge, self.uhf = param, numbers, coords, charge, uhf
                self.max_iterations = None
                tests.calcs.append(self)
            def set_verbosity(self, v):
                self.verbosity = v
            def set_max_iterations(self, n):
                self.max_iterations = n
            def singlepoint(self):
                if Calculator.fail_first and len(tests.calcs) == 1:
                    raise XTBException("SCF did not converge")
                return Result(self)

        self.Calculator = Calculator
        xtb = types.ModuleType('xtb')
        interface = types.ModuleType('xtb.interface')
        interface.Calculator = Calculator
        utils = types.ModuleType('xtb.utils')
        utils.get_method = lambda name: name if name.startswith('GFN') else None
        utils.get_solvent = lambda name: None
        xtb.interface, xtb.utils = interface, utils
        self._saved = {k: sys.modules.get(k) for k in ('xtb', 'xtb.interface', 'xtb.utils')}
        sys.modules.update({'xtb': xtb, 'xtb.interface': interface, 'xtb.utils': utils})

    def tearDown(self):
        for k, v in self._saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v

    def test_SinglePointProperties(self):
        b = XTBBackend(WATER)
        x = WATER_ANG * A2B
        res = b.evaluate(np.stack([x, x * 1.1]), order={'energy': 1, 'charges': 0, 'dipole': 0})
        self.assertEqual(len(self.calcs), 2)  # one single point per structure for everything
        self.assertTrue(np.allclose(res['energy'][0], [-np.sum(x ** 2), -np.sum((1.1 * x) ** 2)]))
        self.assertTrue(np.allclose(res['energy'][1][1], -2 * 1.1 * x.reshape(-1)))
        self.assertEqual(res['charges'][0].shape, (2, 3))
        self.assertEqual(self.calcs[0].verbosity, 'muted')
        b.evaluate(x, properties='charges', charge=1, multiplicity=2)
        self.assertEqual((self.calcs[-1].charge, self.calcs[-1].uhf), (1, 1))
        with self.assertRaises(ValueError):
            b.evaluate(x, properties='charges', charge=0, multiplicity=2)
        with self.assertRaises(ValueError):
            XTBBackend(WATER, multiplicity=2)
        with self.assertRaises(ValueError):
            b.evaluate(x, order=2, properties='energy')

    def test_Retry(self):
        self.Calculator.fail_first = True
        b = XTBBackend(WATER, retry_iterations=500)
        b.evaluate(WATER_ANG * A2B)
        self.assertEqual(len(self.calcs), 2)
        self.assertEqual(self.calcs[-1].max_iterations, 500)

    def test_Wrappers(self):
        x = WATER_ANG * A2B
        E, G = XTBEnergyEvaluator(WATER, method='GFN1-xTB').evaluate(x, order=1)
        self.assertEqual(self.calcs[-1].param, 'GFN1-xTB')
        self.assertTrue(np.allclose(G, 2 * -x.reshape(-1) / A2B * A2B))  # Bohr in, Bohr out
        q = XTBChargeEvaluator(WATER, charge=1, multiplicity=2).evaluate(x)[0]
        self.assertTrue(np.allclose(q, 1 / 3))
        mu = XTBDipoleEvaluator(WATER).evaluate(x)[0]
        self.assertEqual(mu.shape, (3,))


@unittest.skipUnless(has_module('pyscf'), "needs pyscf")
class PySCFBackendTests(unittest.TestCase):

    def setUp(self):
        self.x = WATER_ANG * A2B

    def test_EnergyAndGradient(self):
        b = PySCFBackend(WATER, level_of_theory='hf', basis_set='sto-3g')
        E, G = b.evaluate_property(self.x, 'energy', order=1)
        h = 1e-4
        fd = []
        for i in range(9):
            dx = np.zeros(9)
            dx[i] = h
            fd.append((b.evaluate_property(self.x + dx.reshape(3, 3), 'energy')[0]
                       - b.evaluate_property(self.x - dx.reshape(3, 3), 'energy')[0]) / (2 * h))
        self.assertTrue(np.allclose(G, fd, atol=1e-6))

    def test_SCFProperties(self):
        from pyscf import gto, scf
        b = PySCFBackend(WATER, level_of_theory='hf', basis_set='sto-3g')
        res = b.evaluate(self.x, properties=('energy', 'charges', 'dipole'))
        mol = gto.M(atom=list(zip(WATER, self.x)), basis='sto-3g', unit='Bohr', verbose=0)
        mf = scf.RHF(mol).run()
        self.assertTrue(np.isclose(res['energy'][0], mf.e_tot))
        self.assertTrue(np.allclose(res['charges'][0], mf.mulliken_pop(verbose=0)[1]))
        self.assertTrue(np.isclose(res['charges'][0].sum(), 0))
        self.assertTrue(np.allclose(res['dipole'][0], mf.dip_moment(unit='AU', verbose=0)))
        cation = b.evaluate(self.x, properties='charges', charge=1, multiplicity=2)['charges'][0]
        self.assertTrue(np.isclose(cation.sum(), 1))
        meta = PySCFBackend(WATER, level_of_theory='hf', basis_set='sto-3g', charge_model='meta_lowdin')
        self.assertTrue(np.isclose(meta.evaluate_property(self.x, 'charges')[0].sum(), 0))
        mp2 = PySCFBackend(WATER, level_of_theory='mp2', basis_set='sto-3g')
        with self.assertRaises(NotImplementedError):
            mp2.evaluate_property(self.x, 'dipole')

    def test_Wrappers(self):
        lot = dict(level_of_theory='hf', basis_set='sto-3g')
        E, G = PySCFEnergyEvaluator(WATER, **lot).evaluate(self.x, order=1)
        q = PySCFChargeEvaluator(WATER, **lot).evaluate(self.x)[0]
        mu = PySCFDipoleEvaluator(WATER, **lot).evaluate(self.x)[0]
        b = PySCFBackend(WATER, **lot)
        self.assertTrue(np.allclose(q, b.evaluate_property(self.x, 'charges')[0]))
        self.assertTrue(np.allclose(mu, b.evaluate_property(self.x, 'dipole')[0]))  # already e Bohr
        # backend options and SCF settings are split from evaluator options as before
        ev = PySCFEnergyEvaluator(WATER, **lot, conv_tol=1e-11, mesh_spacing=.002)
        self.assertEqual(ev.defaults, {'conv_tol': 1e-11, 'mesh_spacing': .002})
        self.assertTrue(np.isclose(ev.evaluate(self.x)[0], E, atol=1e-8))
        self.assertEqual(PySCFEnergyEvaluator(WATER, **lot, analytic_hessian=True).analytic_derivative_order, 2)

    @unittest.skipUnless(has_module('pyscf.prop'), "needs pyscf-properties")
    def test_Polarizability(self):
        lot = dict(level_of_theory='hf', basis_set='6-31g')
        alpha = PySCFBackend(WATER, **lot).evaluate_property(self.x, 'polarizability')[0]
        self.assertTrue(np.allclose(alpha, alpha.T, atol=1e-6))
        self.assertTrue(np.all(np.linalg.eigvalsh(alpha) > 0))
        # finite-field check: alpha = d mu / d F
        from pyscf import gto, scf
        mol = gto.M(atom=list(zip(WATER, self.x)), basis='6-31g', unit='Bohr', verbose=0)
        def mu(F):
            mf = scf.RHF(mol)
            h1 = mf.get_hcore()
            mf.get_hcore = lambda *a: h1 + np.einsum('x,xij->ij', F, mol.intor('int1e_r'))
            mf.conv_tol = 1e-12
            mf.run()
            return mf.dip_moment(unit='AU', verbose=0)
        f = 1e-4
        # the electrons' energy in a field F is +F.r, so the induced dipole is +alpha F
        ff = np.array([(mu(f * e) - mu(-f * e)) / (2 * f) for e in np.eye(3)]).T
        self.assertTrue(np.allclose(alpha, ff, atol=1e-3))
        mu0, alpha_ev = PySCFDipolePolarizabilityEvaluator(WATER, **lot).evaluate(self.x, order=0)
        self.assertTrue(np.allclose(alpha_ev[0], alpha))


@unittest.skipUnless(has_module('ase'), "needs ase")
class ASEBackendTests(unittest.TestCase):

    def make_calculator(self):
        from ase.calculators.calculator import Calculator, all_changes

        class Harmonic(Calculator):
            implemented_properties = ['energy', 'forces', 'charges', 'dipole']
            ncalls = 0
            def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
                super().calculate(atoms, properties, system_changes)
                Harmonic.ncalls += 1
                x = self.atoms.positions
                q = np.full(len(x), self.atoms.info.get('charge', 0) / len(x))
                self.results = {
                    'energy': .5 * np.sum(x ** 2), 'forces': -x, 'charges': q, 'dipole': q @ x,
                    'info_spin': float(self.atoms.info.get('spin', -1)),
                }
            def get_hessian(self, atoms):
                return np.eye(3 * len(atoms))
        return Harmonic

    def test_Properties(self):
        Harmonic = self.make_calculator()
        b = ASECalculatorBackend(WATER, calculator=Harmonic())
        x = np.stack([WATER_ANG, 1.1 * WATER_ANG])
        Harmonic.ncalls = 0
        res = b.evaluate(x, order={'energy': 1, 'charges': 0, 'dipole': 0}, charge=1)
        self.assertEqual(Harmonic.ncalls, 2)  # one calculator pass per structure for everything
        self.assertTrue(np.allclose(res['energy'][0], .5 * np.sum(x ** 2, axis=(1, 2))))
        self.assertTrue(np.allclose(res['energy'][1], x.reshape(2, 9)))
        self.assertTrue(np.allclose(res['charges'][0], 1 / 3))
        H = b.evaluate_property(WATER_ANG, 'energy', order=2)[2]
        self.assertTrue(np.allclose(H, np.eye(9)))
        E0 = b.evaluate_property(x, 'energy')[0]  # batched order 0
        self.assertEqual(E0.shape, (2,))
        self.assertEqual(b.evaluate_property(WATER_ANG, 'info_spin', multiplicity=3)[0], 3)

    def test_Wrapper(self):
        Harmonic = self.make_calculator()
        ev = ASECalcEnergyEvaluator(WATER, calc=Harmonic())
        E, G, H = ev.evaluate(WATER_ANG, order=2)
        self.assertTrue(np.allclose(G, WATER_ANG.reshape(-1) * EV2H / A2B))
        self.assertTrue(np.allclose(H, np.eye(9) * EV2H / A2B ** 2))
        self.assertIs(ev.to_ase(), ev.eval)

    def test_MACEPolarBackend(self):
        from ase.calculators.calculator import Calculator, all_changes
        sigma = .8

        class _Model:
            atomic_multipoles_smearing_width = sigma

        class FakePolar(Calculator):
            implemented_properties = ['energy', 'forces', 'density_coefficients', 'spin_charge_density']
            models = [_Model()]
            def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
                super().calculate(atoms, properties, system_changes)
                n = len(self.atoms)
                c = self.atoms.info.get('charge', 0)
                dens = np.zeros((n, 4))
                dens[:, 0] = c / n
                self.results = {'energy': -1., 'forces': np.zeros((n, 3)), 'density_coefficients': dens,
                                'spin_charge_density': np.stack([dens / 2, dens / 2], axis=-2),
                                'field': np.array(self.atoms.info.get('external_field', [0, 0, 0]))}

        b = MACEPolarBackend(WATER, calculator=FakePolar(), external_field=[0, 0, .01])
        self.assertEqual(b.smearing, sigma)
        res = b.evaluate_polar_properties(np.stack([WATER_ANG] * 2), properties=['density_coefficients', 'field'], charge=1)
        self.assertEqual(res['density_coefficients'].shape, (2, 3, 4))
        self.assertTrue(np.allclose(res['density_coefficients'][..., 0], 1 / 3))
        self.assertTrue(np.allclose(res['field'], [0, 0, .01]))
        # the Molecools wrappers share the backend
        esp = MACEPolarElectrostaticPotentialEvaluator(backend=b, points=np.array([[5., 0, 0]]))
        self.assertEqual(esp.model_smearing, sigma)
        self.assertIs(esp.energy_evaluator.backend, b)
        self.assertEqual(esp.evaluate(WATER_ANG)[0].shape, (1,))
        E, G = MACEEnergyEvaluator(backend=b).evaluate(WATER_ANG, order=1)
        self.assertTrue(np.isclose(E, -EV2H))


def _aimnet_available():
    if not (has_module('aimnet') and has_module('torch')):
        return False
    cache = os.environ.get('MODEL_CACHE_DIR')
    return cache is not None and os.path.exists(os.path.join(cache, 'aimnet', 'aimnet2_wb97m_d3_0.pt'))


@unittest.skipUnless(_aimnet_available(), "needs aimnet, torch, and the AIMNet2 weights in $MODEL_CACHE_DIR/aimnet")
class AIMNet2BackendTests(unittest.TestCase):

    def test_EnergiesAndCharges(self):
        from Psience.ElectronicStructure import AIMNet2Backend
        from Psience.Molecools.Evaluator import AIMNet2EnergyEvaluator, AIMNet2ChargeEvaluator
        b = AIMNet2Backend(WATER, device='cpu')
        res = b.evaluate(WATER_ANG, order={'energy': 2, 'charges': 1})
        self.assertEqual([r.shape for r in res['energy']], [(), (9,), (9, 9)])
        self.assertEqual([r.shape for r in res['charges']], [(3,), (9, 3)])
        self.assertTrue(np.isclose(res['charges'][0].sum(), 0, atol=1e-5))
        # batched charge derivatives are taken structure by structure
        batch = b.evaluate(np.stack([WATER_ANG, 1.01 * WATER_ANG]), order=1, properties='charges')['charges']
        single = b.evaluate(1.01 * WATER_ANG, order=1, properties='charges')['charges']
        self.assertTrue(np.allclose(batch[1][1], single[1]))
        cation = b.evaluate_property(WATER_ANG, 'charges', charge=1)[0]
        self.assertTrue(np.isclose(cation.sum(), 1, atol=1e-5))
        ev = AIMNet2EnergyEvaluator(backend=b)
        E, G = ev.evaluate(WATER_ANG, order=1)
        self.assertTrue(np.allclose(G, res['energy'][1] * EV2H / A2B))
        q0, q1, q2 = AIMNet2ChargeEvaluator(backend=b).evaluate(WATER_ANG, order=2)
        self.assertEqual(q2.shape, (9, 9, 3))


if __name__ == '__main__':
    unittest.main()
