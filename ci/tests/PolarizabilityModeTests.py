"""Response physics and public molecular routing, using deterministic torch graphs.

No pretrained weights or network access are required. The fixture deliberately
evaluates learned energy before exporting a squeezed charge view, reproducing
the graph distinction that matters for the total-energy response.
"""

import os as _test_os
import unittest as _test_unittest
_RUN_EXTRAS = _test_os.environ.get("MCUTILS_TEST_EXTRAS") == "1"

import unittest
import numpy as np
from McUtils.Data import UnitsData
from Psience.ElectronicStructure import AIMNet2Backend, FluctuatingChargeResponse
from Psience.Molecools import Molecule
from Psience.Molecools.Evaluator import AIMNet2DipolePolarizabilityEvaluator

A2B=UnitsData.convert('Angstroms','BohrRadius')
EV2H=UnitsData.convert('ElectronVolts','Hartrees')
XYZ=np.array([[0.,0.,0.],[1.,.1,0.],[-.3,.9,0.]])
def _load_torch_fixtures():
    global torch, backend
    import torch
    _TorchModule=torch.nn.Module

    class LRCoulomb(_TorchModule):
        method='simple'
        def coul_simple(self,data):
            q=data['charges'].reshape(-1)[:3]
            x=data['coord'][:3]
            return sum(1.5*q[i]*q[j]/torch.linalg.vector_norm(x[i]-x[j]) for i in range(3) for j in range(i))

    class SRCoulomb(_TorchModule):
        key_out='energy'
        def forward(self,data):
            correction=-.5*LRCoulomb().coul_simple(data)
            data[self.key_out]=data.get(self.key_out,0)+correction
            return data

    class NeuralFixture(_TorchModule):
        def __init__(self):
            super().__init__()
            self.parameter=torch.nn.Parameter(torch.tensor(1.,dtype=torch.float64))
            self.mlps=torch.nn.ModuleList([torch.nn.Identity() for _ in range(3)])
            self.outputs=torch.nn.ModuleDict({'srcoulomb':SRCoulomb()})
            self.fail=False

        def _update_q(self,data,x,delta_q=True):
            r=data['coord'][:3]
            native=.1*(r[:,0]**2+.3*r[:,1])
            native=native-native.mean()+data['charge']/3
            data['charges']=torch.cat((native,native.new_zeros(1))).reshape(-1,1)
            return data

        def forward(self,data):
            for _ in range(2):data=self._update_q(data,None)
            if self.fail:raise RuntimeError('synthetic forward failure')
            q=data['charges'][:3,0]
            r=data['coord'][:3]
            data['energy']=(8*q*q+(1+.1*r[:,0]**2)*q**4).sum()
            data['charges']=data['charges'].squeeze(-1)
            return self.outputs.srcoulomb(data)

    class CalculatorFixture:
        def __init__(self):
            self.model=NeuralFixture()
            self.external_coulomb=LRCoulomb()
            self.metadata={'coulomb_mode':'sr_embedded'}

    class FixtureBackend(AIMNet2Backend):
        @classmethod
        def setup_model(cls,model,**opts):return model

        def _prep_eval_unlocked(self,coords,order,charge=None,multiplicity=None,**opts):
            charge,multiplicity=self.resolve_state(charge,multiplicity)
            x=torch.tensor(np.asarray(coords).reshape(-1,3),dtype=torch.float64,requires_grad=True)
            coord=torch.cat((x,x.new_zeros(1,3)))
            data={'coord':coord,'charge':charge}
            data=self.model.model(data)
            data['energy']=data['energy']+self.model.external_coulomb.coul_simple(data)
            return np.shape(coords)[:-2],coords,data

        def _native_expansion(self,coords,order,key,charge,multiplicity):
            rows=[]
            for c in coords:
                _,_,data=self.prep_eval(c,0,charge=charge,multiplicity=multiplicity)
                value=data[key] if key=='energy' else data[key][:3]
                expansion=[value.detach().numpy()]
                if order:
                    d=[self._response_autograd(v,data['coord']).reshape(-1,3)[:3].reshape(-1).detach().numpy()
                       for v in value.reshape(-1)]
                    expansion.append(np.moveaxis(np.asarray(d).reshape(tuple(value.shape)+(9,)),-1,0))
                rows.append(expansion)
            return [np.stack([r[k] for r in rows]) for k in range(order+1)]

        def evaluate_charges(self,coords,order,charge=None,multiplicity=None,**opts):
            return self._native_expansion(coords,order,'charges',charge,multiplicity)

        def evaluate_energy(self,coords,order,charge=None,multiplicity=None,**opts):
            return self._native_expansion(coords,order,'energy',charge,multiplicity)

    def backend(**opts):
        return FixtureBackend(['O','H','H'],model=CalculatorFixture(),multiplicity=1,**opts)


class ConservedChargeSolverTests(unittest.TestCase):
    def test_singular_full_hardness_with_positive_neutral_curvature(self):
        Z=np.linalg.svd(np.ones((1,3)),full_matrices=True)[2][1:].T
        J=Z@np.diag([2.,3.])@Z.T
        alpha,diag=FluctuatingChargeResponse.constrained_polarizability(XYZ,J,return_diagnostics=True)
        np.testing.assert_allclose(alpha[0],XYZ.T@Z@np.diag([.5,1/3])@Z.T@XYZ)
        self.assertTrue(diag['stable'])
        np.testing.assert_allclose(diag['charge_response'].sum(axis=0),0,atol=1e-14)

    def test_inverse_equivalence_translation_rotation_and_planarity(self):
        J=np.diag([2.,3.,4.])+.2
        a=FluctuatingChargeResponse.constrained_polarizability(XYZ,J)[0]
        b=FluctuatingChargeResponse.constrained_polarizability(XYZ,J,solver='inverse')[0]
        np.testing.assert_allclose(a,b,atol=1e-14)
        translated=FluctuatingChargeResponse.constrained_polarizability(XYZ+[3,-2,4],J)[0]
        np.testing.assert_allclose(translated,a,atol=1e-14)
        Q=np.linalg.qr(np.random.default_rng(7).normal(size=(3,3)))[0]
        rotated=FluctuatingChargeResponse.constrained_polarizability(XYZ@Q,J)[0]
        np.testing.assert_allclose(rotated,Q.T@a@Q,atol=1e-14)
        self.assertEqual(a[2,2],0)

    def test_stability_singularity_and_derivative(self):
        with self.assertRaises(ValueError):FluctuatingChargeResponse.constrained_polarizability(XYZ,-np.eye(3))
        with self.assertRaises(np.linalg.LinAlgError):FluctuatingChargeResponse.constrained_polarizability(XYZ,np.zeros((3,3)),allow_indefinite=True)
        J=np.diag([2.,3.,4.]);dJ=np.zeros((9,3,3));dJ[1]=np.diag([.1,.2,.3])
        a=FluctuatingChargeResponse.constrained_polarizability(XYZ,J,dJ)
        h=1e-5;step=np.eye(9)[1].reshape(3,3)*h
        finite=(FluctuatingChargeResponse.constrained_polarizability(XYZ+step,J+h*dJ[1])[0]-
                FluctuatingChargeResponse.constrained_polarizability(XYZ-step,J-h*dJ[1])[0])/(2*h)
        np.testing.assert_allclose(a[1][1],finite,atol=1e-10)

    def test_inverse_requires_full_rank_and_single_atom_has_zero_response(self):
        J=np.diag([1.,1.,0.])
        FluctuatingChargeResponse.constrained_polarizability(XYZ,J)
        with self.assertRaises(np.linalg.LinAlgError):
            FluctuatingChargeResponse.constrained_polarizability(XYZ,J,solver='inverse')
        a=FluctuatingChargeResponse.constrained_polarizability(XYZ[:1],np.zeros((1,1)))
        np.testing.assert_array_equal(a[0],np.zeros((3,3)))

@_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
class AIMNetResponseModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _load_torch_fixtures()

    def test_total_captures_learned_energy_and_native_zero_increment(self):
        b=backend(polarizability_mode='total_energy')
        total=b.evaluate_charge_hardness(XYZ)[0]
        net=b.evaluate_charge_hardness(XYZ,polarizability_mode='net_coulomb')[0]
        ext=b.evaluate_charge_hardness(XYZ,polarizability_mode='coulomb')[0]
        np.testing.assert_allclose(net,.5*ext,atol=1e-13)
        self.assertTrue(np.all(np.linalg.eigvalsh(total-net)>15))
        _,_,native=b.prep_eval(XYZ,0)
        original=float(native['energy'].detach())
        with b._model_lock:
            E,delta,data=b._charge_energy_graph(XYZ,'total_energy',0,1)
        self.assertEqual(float(E.detach()),original)
        np.testing.assert_array_equal(data['charges'].detach(),native['charges'].detach())
        # The exported view cannot see upstream learned energy.
        g=b._response_autograd(E,data['charges'],True)
        naive=torch.stack([b._response_autograd(v,data['charges'])[:3] for v in g[:3]]).detach().numpy()
        np.testing.assert_allclose(naive,net,atol=1e-13)
        v=np.array([1.,-1.,0.])/np.sqrt(2);h=1e-3
        values=[]
        for sign in [-1,1]:
            increment=torch.tensor(sign*h*v,dtype=torch.float64,requires_grad=True)
            with b._independent_charge_increment(increment):_,_,d=b.prep_eval(XYZ,0)
            values.append(float(d['energy'].detach()))
        self.assertAlmostEqual((sum(values)-2*original)/h**2,v@total@v,places=5)

    def test_hook_restored_after_forward_error(self):
        b=backend(polarizability_mode='total_energy');model=b.model.model
        original=model._update_q;model.fail=True
        with self.assertRaisesRegex(RuntimeError,'synthetic forward failure'):
            b.evaluate_charge_hardness(XYZ)
        self.assertNotIn('_update_q',model.__dict__)
        self.assertEqual(model._update_q,original)

    def test_unsupported_and_invalid_modes_fail_explicitly(self):
        with self.assertRaises(ValueError):backend(polarizability_mode='guess')
        with self.assertRaises(ValueError):backend(response_solver='guess')
        b=backend();b.model.model.mlps=torch.nn.ModuleList([torch.nn.Identity()])
        with self.assertRaises(NotImplementedError):b.evaluate_charge_hardness(XYZ,polarizability_mode='total_energy')

    def test_stability_opt_in_and_analytic_nuclear_derivative(self):
        b=backend()
        with self.assertRaises(ValueError):b.evaluate_property(XYZ,'polarizability')
        for mode in b.polarizability_modes:
            a,diag=b.evaluate_polarizability(XYZ,1,polarizability_mode=mode,allow_indefinite=True,return_diagnostics=True)
            self.assertEqual(bool(diag['stable']),mode=='total_energy')
            h=1e-5;step=np.eye(9)[0].reshape(3,3)*h
            plus=b.evaluate_polarizability(XYZ+step,polarizability_mode=mode,allow_indefinite=True)[0]
            minus=b.evaluate_polarizability(XYZ-step,polarizability_mode=mode,allow_indefinite=True)[0]
            np.testing.assert_allclose(a[1][0],(plus-minus)/(2*h),rtol=2e-7,atol=1e-8)
            np.testing.assert_allclose(a[0][2],0,atol=1e-14)

    def test_batch_shapes_and_state_override(self):
        b=backend(polarizability_mode='total_energy')
        x=np.stack([XYZ,XYZ*1.01,XYZ*.99,XYZ]).reshape(2,2,3,3)
        r=b.evaluate(x,order={'energy':0,'charges':1,'dipole':1,'polarizability':1},charge=1,multiplicity=2)
        self.assertEqual(r['polarizability'][0].shape,(2,2,3,3))
        self.assertEqual(r['polarizability'][1].shape,(2,2,9,3,3))
        np.testing.assert_allclose(r['charges'][0].sum(axis=-1),1)
        response,diag=b.evaluate_polarizability(x,return_diagnostics=True,charge=1,multiplicity=2)
        self.assertEqual(diag['constrained_gradient_norm'].shape,(2,2))
        np.testing.assert_allclose(response[0],r['polarizability'][0])

    def test_standard_molecular_units_defaults_and_call_override(self):
        b=backend(polarizability_mode='total_energy')
        ev=AIMNet2DipolePolarizabilityEvaluator(backend=b)
        raw=b.evaluate_field_response(XYZ,1)
        molecular=ev.evaluate(XYZ,order=1)
        for field in [0,1]:
            np.testing.assert_allclose(molecular[field][0],raw[field][0])
            np.testing.assert_allclose(molecular[field][1],raw[field][1]/A2B)
        other=ev.evaluate(XYZ,order=0,polarizability_mode='coulomb',allow_indefinite=True)
        expected=b.evaluate_polarizability(XYZ,polarizability_mode='coulomb',allow_indefinite=True)
        np.testing.assert_allclose(other[1][0],expected[0])
        np.testing.assert_allclose(other[0][0],raw[0][0])
        constructed=AIMNet2DipolePolarizabilityEvaluator(['O','H','H'],model=b.model,
            polarizability_mode='total_energy')
        self.assertEqual(constructed.backend.polarizability_mode,'total_energy')
        mol=Molecule(['O','H','H'],XYZ*A2B,charge=0,spin=1)
        mu,alpha=mol.calculate_dipole_polarizability('aimnet2',backend=b)
        np.testing.assert_allclose(alpha,raw[1][0]);np.testing.assert_allclose(mu,raw[0][0])
        function=mol.get_dipole_polarizability_function('aimnet2',backend=b)
        fm,fa=function(mol.coords)
        np.testing.assert_allclose(fm,mu);np.testing.assert_allclose(fa,alpha)

    def test_molecular_second_derivatives_use_standard_finite_differences(self):
        b=backend(polarizability_mode='total_energy')
        ev=AIMNet2DipolePolarizabilityEvaluator(backend=b)
        value=ev.evaluate(XYZ,order=2,mesh_spacing=1e-3,stencil=3)
        self.assertEqual(value[0][2].shape,(9,9,3))
        self.assertEqual(value[1][2].shape,(9,9,3,3))
        h=1e-4;step=np.eye(9)[0].reshape(3,3)*h
        plus=ev.evaluate(XYZ+step,order=1)[1][1]
        minus=ev.evaluate(XYZ-step,order=1)[1][1]
        np.testing.assert_allclose(value[1][2][0],(plus-minus)/(2*h*A2B),rtol=2e-4,atol=1e-5)

    def test_field_order_and_mode_override_survive_finite_differences(self):
        b=backend()
        self.assertEqual(b.evaluate_field_response(XYZ,field_order=0),[])
        self.assertEqual(len(b.evaluate_field_response(XYZ,field_order=1)),1)
        self.assertEqual(len(b.evaluate_electrostatic_potential_derivatives(XYZ,0,1)),2)
        ev=AIMNet2DipolePolarizabilityEvaluator(backend=b)
        # The constructor default is unstable Coulomb. Dropping this per-call
        # override in finite differences would fail or use the wrong curvature.
        value=ev.evaluate(XYZ,order=2,polarizability_mode='total_energy',mesh_spacing=1e-3,stencil=3)
        self.assertEqual(value[1][2].shape,(9,9,3,3))

if __name__=='__main__':unittest.main()
