"""VPT2 fixture comparisons for direct molecular Hamiltonian expansions."""
import os
import itertools
import math
from pathlib import Path
import unittest
import numpy as np
from McUtils import Numputils as nput
from McUtils.Scaffolding import NullLogger
from Psience.Molecools import Molecule
from Psience.Molecools.Hamiltonian import ZetaExpansion, ReciprocalInertiaExpansion, PseudopotentialExpansion
from Psience.Modes import NormalModes
from Psience.VPT2.Terms import CoriolisTerm, PotentialLikeTerm, PotentialTerms, KineticTerms


class MolecularHamiltonianTermTests(unittest.TestCase):
    results = []
    fixtures = ('HOH_freq.fchk', 'HOD_freq.fchk', 'DOD_freq.fchk', 'HOT_freq.fchk',
                'OCHH_freq.fchk', 'OCHD_freq.fchk', 'OCHT_freq.fchk', 'HOONO_freq.fchk')

    @classmethod
    def molecule(cls, name):
        directory = Path(os.environ.get('PSIENCE_TEST_DATA', Path(__file__).parent/'TestData'))
        return Molecule.from_file(str(directory/name))

    def assertExpansion(self, actual, expected, label, atol=2e-13, rtol=2e-10):
        self.assertEqual(len(actual), len(expected), label)
        errors = []
        for o,(a,b) in enumerate(zip(actual,expected)):
            with self.subTest(component=label, order=o):
                np.testing.assert_allclose(a,b,atol=atol,rtol=rtol)
            errors.append(float(np.max(np.abs(np.asarray(a)-np.asarray(b)))))
        self.results.append({'comparison':label,'max_errors':errors})

    def test_cartesian_zeta_coriolis_and_watson_against_vpt(self):
        for name in self.fixtures:
            with self.subTest(system=name):
                mol = self.molecule(name)
                mol.internals = None
                ham = mol.hamiltonian
                embedding = ham.get_embedding(dimensionless=True)
                reference = CoriolisTerm(mol, logger=NullLogger())
                zeta = ZetaExpansion(embedding).get_terms(2)
                self.assertExpansion(zeta,[reference.get_zetas(),0,0],name+' zeta')
                inverse = ReciprocalInertiaExpansion(embedding).get_terms(2)
                expected_inverse = [np.moveaxis(t,(0,1),(-2,-1))
                                    for t in reference.moment_of_inertia_derivs(2)]
                self.assertExpansion(inverse,expected_inverse,name+' reciprocal inertia')
                expected = [t[0,0]+t[1,1]+t[2,2] for t in reference.get_terms(2)]
                actual = ham.get_VPT_expansions({'coriolis':2,'watson':2})
                self.assertExpansion(actual['coriolis'],expected,name+' coriolis')
                watson = PotentialLikeTerm(mol,logger=NullLogger()).get_terms(2)
                self.assertExpansion(actual['watson'],watson,name+' watson')

    def test_mode_signs_permutations_and_reduced_spaces(self):
        for name in ('HOD_freq.fchk','OCHD_freq.fchk'):
            mol = self.molecule(name)
            mol.internals = None
            original = mol.hamiltonian.get_embedding().modes
            n = len(original.freqs)
            variants = (original[np.arange(n)[::-1]],original[np.arange(min(n,3))],
                        original.modify(original.matrix*((-1.)**np.arange(n))[None],
                            inverse=original.inverse*((-1.)**np.arange(n))[:,None]))
            for index,modes in enumerate(variants):
                ham = mol.get_hamiltonian(modes=modes)
                reference = CoriolisTerm(mol,modes=modes,logger=NullLogger())
                expected = [t[0,0]+t[1,1]+t[2,2] for t in reference.get_terms(2)]
                self.assertExpansion(ham.get_VPT_expansions({'coriolis':2})['coriolis'],
                                     expected,name+' variant '+str(index))

    def test_cartesian_frequency_scaling(self):
        mol = self.molecule('OCHD_freq.fchk')
        mol.internals = None
        ham = mol.hamiltonian
        physical = ham.get_embedding(dimensionless=False)
        scaled = ham.get_embedding(dimensionless=True)
        self.assertExpansion(ZetaExpansion(physical).get_terms(2),
                             ZetaExpansion(scaled).get_terms(2),'physical zeta scaling')
        physical_terms = ham.coriolis_expansion(0,embedding=physical)
        scaled_terms = ham.coriolis_expansion(0,embedding=scaled)
        roots = np.sqrt(np.abs(scaled.modes.freqs))
        ratio = roots[None,:]/roots[:,None]
        expected = [physical_terms[0]*ratio[:,:,None,None]*ratio[None,None,:,:]]
        self.assertExpansion(scaled_terms,expected,'q/p frequency scaling')

    def test_complete_cartesian_vpt2_assembly(self):
        for name in ('HOH_freq.fchk','OCHD_freq.fchk'):
            mol = self.molecule(name)
            mol.internals = None
            actual = mol.hamiltonian.get_VPT_expansions(
                {'potential':4,'kinetic':2,'coriolis':0,'watson':0})
            self.assertExpansion(actual['potential'],PotentialTerms(mol,logger=NullLogger()).get_terms(2),
                                 name+' potential assembly')
            self.assertEqual(len(mol.hamiltonian.potential_expansion(2)),2)
            self.assertExpansion(actual['kinetic'],KineticTerms(mol,logger=NullLogger()).get_terms(2),
                                 name+' kinetic assembly')
            ref = CoriolisTerm(mol,logger=NullLogger()).get_terms(0)[0]
            self.assertExpansion(actual['coriolis'],[ref[0,0]+ref[1,1]+ref[2,2]],name+' Coriolis assembly')
            self.assertExpansion(actual['watson'],PotentialLikeTerm(mol,logger=NullLogger()).get_terms(0),
                                 name+' Watson assembly')
            automatic = mol.hamiltonian.get_VPT_expansions(2)
            for key in actual:
                self.assertExpansion(automatic[key],actual[key],name+' default '+key)
            for order in (0,1):
                lower = mol.hamiltonian.get_VPT_expansions(order)
                self.assertEqual(set(lower),{'potential','kinetic'})
                self.assertEqual(len(lower['potential']),order+1)
                self.assertEqual(len(lower['kinetic']),order+1)
            with self.assertRaises(ValueError):
                mol.hamiltonian.get_VPT_expansions(-1)

    def test_watson_coordinate_reexpansion(self):
        mol = self.molecule('HOD_freq.fchk')
        mol.internals = None
        ham = mol.hamiltonian
        generator = ham.watson_expansion(dimensionless=True)
        original = generator.get_terms(2)
        forward = [np.array([[1.,.1,0.],[0.,.8,.2],[.1,0.,1.2]]),
                   np.ones((3,3,3))*.01]
        inverse = nput.inverse_transformation(forward,2)
        expected = [original[0]]+nput.tensor_reexpand(forward,original[1:],2)
        self.assertExpansion(generator.get_terms(2,transformation=(forward,inverse)),
                             expected,'Watson coordinate transform')

    def test_watson_against_independent_matrix_polynomial(self):
        # Watson's coefficient model is B(q) I_e B(q), with
        # B(q) = inverse(I_e + I'_e(q)/2). Differentiate its independent
        # polynomial coefficients, including factorials and mixed slots.
        mol = self.molecule('HOD_freq.fchk')
        mol.internals = None
        embedding = mol.hamiltonian.get_embedding(dimensionless=True)
        generator = ReciprocalInertiaExpansion(embedding)
        base = generator.get_terms(0)[0]
        inertia = np.linalg.inv(base)
        first = embedding.get_inertia_tensor_expansion(1)[1]
        s = {(1,0):first[0]/2,(0,1):first[1]/2}
        inverse = {(0,0):base}
        keys = sorted((a for a in itertools.product(range(5),repeat=2) if sum(a)<=4),key=sum)
        for key in keys[1:]:
            residual = sum(t@inverse[tuple(a-b for a,b in zip(key,index))]
                for index,t in s.items() if all(a>=b for a,b in zip(key,index)))
            inverse[key] = -base@residual
        polynomial = {}
        for a,x in inverse.items():
            for b,y in inverse.items():
                key = tuple(i+j for i,j in zip(a,b))
                if sum(key)<=4:
                    polynomial[key] = polynomial.get(key,0)+x@inertia@y
        expected = []
        for o in range(5):
            derivative = np.empty((2,)*o+(3,3))
            for slots in itertools.product(range(2),repeat=o):
                key = tuple(slots.count(i) for i in range(2))
                derivative[slots] = polynomial[key]*math.prod(math.factorial(a) for a in key)
            expected.append(derivative)
        actual = [t[(slice(0,2),)*o] for o,t in enumerate(generator.get_terms(4))]
        self.assertExpansion(actual,expected,'independent Watson polynomial',atol=1e-13)
        reference = CoriolisTerm(mol,logger=NullLogger()).moment_of_inertia_derivs(4)
        reference = [np.moveaxis(t,(0,1),(-2,-1))[(slice(0,2),)*o]
                     for o,t in enumerate(reference)]
        self.assertExpansion(reference,expected,'VPT Watson polynomial',atol=1e-13)

    def test_coriolis_gaussian_reference_magnitudes(self):
        mol = self.molecule('HOH_freq.fchk')
        mol.internals = None
        zeta = ZetaExpansion(mol.hamiltonian.get_embedding()).get_terms(0)[0]
        expected = np.zeros((3,3))
        # Gaussian values from VPTComponentTests.test_HOHCoriolisCouplings,
        # remapped from decreasing to increasing frequency. Magnitudes remove
        # the arbitrary signs of normal modes and principal-axis vectors.
        expected[1,2] = expected[2,1] = .00960228
        expected[0,2] = expected[2,0] = .999954
        np.testing.assert_allclose(np.sqrt(np.sum(zeta*zeta,axis=0)),expected,atol=1e-3)

    def test_raw_vprime_through_fourth_derivatives(self):
        rng = np.random.default_rng(772)
        keys = [a for a in itertools.product(range(7),repeat=2) if sum(a)<=6]
        def multiply(a,b,degree):
            result = {}
            for i,x in a.items():
                for j,y in b.items():
                    key = tuple(u+v for u,v in zip(i,j))
                    if sum(key)<=degree:
                        result[key] = result.get(key,0)+x*y
            return result
        def determinant(matrix,size):
            result = {}
            for p in itertools.permutations(range(size)):
                inversions = sum(p[i]>p[j] for i in range(size) for j in range(i+1,size))
                term = {(0,0):np.array((-1.)**inversions)}
                for i,j in enumerate(p):
                    term = multiply(term,{k:v[i,j] for k,v in matrix.items()},6)
                for k,v in term.items():
                    result[k] = result.get(k,0)+v
            return result
        def logarithm(a):
            base = a[(0,0)]
            delta = {k:v/base for k,v in a.items() if sum(k)>0}
            result,term = {(0,0):np.log(base)},{(0,0):np.array(1.)}
            for o in range(1,7):
                term = multiply(term,delta,6)
                for k,v in term.items():
                    result[k] = result.get(k,0)+(-1)**(o+1)*v/o
            return result
        def derivative(a,axis):
            result = {}
            for k,v in a.items():
                if k[axis]:
                    index = list(k)
                    index[axis] -= 1
                    result[tuple(index)] = k[axis]*v
            return result
        def expand(a,order):
            base = np.asarray(a[(0,0)])
            result = []
            for o in range(order+1):
                value = np.empty((2,)*o+base.shape)
                for slots in itertools.product(range(2),repeat=o):
                    index = tuple(slots.count(i) for i in range(2))
                    value[slots] = a.get(index,0)*math.prod(math.factorial(i) for i in index)
                result.append(value)
            return result
        def matrix(size):
            result = {}
            for k in keys:
                value = rng.normal(size=(size,size))*.01
                result[k] = (value+value.T)/2
            result[(0,0)] = np.diag(np.arange(size)+2.)
            return result
        g,inertia = matrix(2),matrix(3)
        log_g,log_i = logarithm(determinant(g,2)),logarithm(determinant(inertia,3))
        gamma = {k:log_i.get(k,0)-log_g.get(k,0) for k in log_i.keys()|log_g.keys()}
        lambdas = [derivative(gamma,i) for i in range(2)]
        expected = {}
        for i in range(2):
            for j in range(2):
                metric = {k:v[i,j] for k,v in g.items()}
                products = (multiply(metric,derivative(lambdas[j],i),4),
                    {k:.25*v for k,v in multiply(metric,multiply(lambdas[i],lambdas[j],4),4).items()},
                    multiply(derivative(metric,i),lambdas[j],4))
                for terms in products:
                    for k,v in terms.items():
                        expected[k] = expected.get(k,0)+v
        self.assertExpansion(PseudopotentialExpansion.get_U(expand(g,6),expand(inertia,6)),
                             expand(expected,4),'independent fourth-order Vprime',atol=2e-10)

    def test_internal_vprime_against_vpt_test_systems(self):
        for name in ('HOH_freq.fchk','HOD_freq.fchk','DOD_freq.fchk','HOT_freq.fchk',
                     'OCHH_freq.fchk','OCHD_freq.fchk'):
            for chart in ('zmatrix','primitive'):
                with self.subTest(system=name,chart=chart):
                    mol = self.molecule(name)
                    if mol.num_atoms == 3:
                        zmatrix = [[0,-1,-1,-1],[1,0,-1,-1],[2,0,1,-1]]
                        specs = [[0,1],[0,2],[1,0,2]]
                    else:
                        zmatrix = [[0,-1,-1,-1],[1,0,-1,-1],[2,0,1,-1],[3,0,1,2]]
                        specs = [[0,1],[0,2],[2,0,1],[0,3],[3,0,1],[3,0,1,2]]
                    mol.internals = zmatrix if chart == 'zmatrix' else {'specs':specs}
                    mol.internal_coordinates.system.name = 'AuditInternals'
                    n = len(specs)
                    rng = np.random.default_rng(931+n)
                    rotation = np.linalg.qr(rng.normal(size=(n,n)))[0]
                    matrix = rotation*np.linspace(.7,1.3,n)[:,None]
                    modes = NormalModes(mol.internal_coordinates.system,matrix,
                        inverse=np.linalg.inv(matrix),freqs=np.ones(n),
                        mass_weighted=True,frequency_scaled=True)
                    ham = mol.get_hamiltonian(modes=modes)
                    reference = PotentialLikeTerm(mol,modes=modes,undimensionalize=False,
                        expansion_handling_mode='analytical',gmatrix_tolerance=None,
                        check_input_gmatrix=False,logger=NullLogger()).get_terms(2)
                    self.assertExpansion(ham.pseudopotential_expansion(2),reference,
                                         name+' '+chart+' Vprime',atol=2e-10,rtol=2e-8)
                    self.assertExpansion(ham.get_VPT_expansions({'pseudopotential':2})['pseudopotential'],
                                         reference,name+' '+chart+' assembled Vprime',atol=2e-10,rtol=2e-8)


if __name__ == '__main__':
    unittest.main(verbosity=2)
