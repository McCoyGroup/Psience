"""
Tests for the experimental `ContractedDVRHarmonicRepresentation` (and the
`ContractedBasisHarmonicRepresentation` machinery it sits on top of) in
`Psience.BasisReps.ContractedRepresentations`.

The "contracted basis" used throughout is the eigenbasis of a 1D Morse DVR
(parameterized like an OH stretch), coupled to a small `HarmonicOscillatorProductBasis`
"bath". Couplings are potential-like only (`coupling_potential`), taken through
cubic order in the H.O. coordinates; `coupling_kinetic` is never supplied, so
there is no kinetic coupling.

Rather than trusting the class's own internals, each numerical check rebuilds
its expectation independently: DVR matrix elements are recomputed via the
public `DVRWavefunctions.expectation` API, and the corresponding bare
H.O.-basis operator matrices are built directly with `basis.representation(...)`
(the same lower-level, already-tested entry point `ContractedDVRHarmonicRepresentation`
itself calls). This isolates what's actually new here: averaging the coupling
tensors over the DVR wavefunctions and assembling/caching the resulting blocks.
"""

import math
from unittest import TestCase

import numpy as np

from Peeves.TestUtils import *

from McUtils.Data import UnitsData, PotentialData, AtomData

from Psience.DVR import CartesianDVR
from Psience.BasisReps import ContractedDVRHarmonicRepresentation


class ContractedDVRHarmonicRepresentationTests(TestCase):

    #region setup helpers

    def get_morse_wavefunctions(self, n_states=4, domain=(0.8, 6.0), divs=200,
                                 w=3869.47, wx=84.0, re=1.82534):
        """
        Solves a 1D Colbert-Miller DVR for an OH-stretch-like Morse potential
        (parameterized from a harmonic frequency/anharmonicity pair, following
        the same `w`/`wx` -> `De`/`alpha` conversion used elsewhere in the
        package, e.g. `Psience.DGB.Runners.setupMorseFunction`) and returns
        the lowest `n_states` DVR eigenstates together with the potential
        parameters used to build them.
        """
        w = w / UnitsData.hartrees_to_wavenumbers
        wx = wx / UnitsData.hartrees_to_wavenumbers
        De = (w ** 2) / (4 * wx)

        m1 = AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        m2 = AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        mu = 1 / (1 / m1 + 1 / m2)
        alpha = np.sqrt(2 * wx / (1 / m1 + 1 / m2))

        morse = PotentialData["MorsePotential"]["Function"]
        def pot(r):
            return morse(r, re=re, alpha=alpha, De=De)

        dvr = CartesianDVR(domain=domain, divs=divs, mass=mu)
        wfns = dvr.run(potential_function=pot).wavefunctions

        return wfns[:n_states], dict(re=re, De=De, alpha=alpha, mass=mu)

    @staticmethod
    def outer_power(vec, order):
        """(n_ho,) -> (n_ho,)*order via repeated outer products (order 0 -> a bare scalar 1.0)"""
        t = np.array(1.0)
        for _ in range(order):
            t = np.multiply.outer(t, vec)
        return t

    def get_coupling_potential(self, re, v, scales):
        """
        Builds a `coupling_potential(grid)` callable whose order-`o` term is
        `scales[o] * (grid - re)**o * v(x)v...v` (`o` copies of `v`, via
        `outer_power`) -- i.e. every order is proportional to the same rank-1
        H.O.-side tensor, scaled by a coordinate-dependent prefactor and an
        order-specific magnitude. `scales` gives the order-0...order-k
        magnitudes (a `None` entry drops that order entirely).
        """
        def coupling_potential(grid):
            dr = grid - re
            terms = []
            for order, scale in enumerate(scales):
                if scale is None:
                    terms.append(None)
                elif order == 0:
                    terms.append(scale * np.ones_like(grid))
                else:
                    terms.append(scale * np.multiply.outer(dr ** order, self.outer_power(v, order)))
            return terms
        return coupling_potential

    @staticmethod
    def dense_operator_matrix(basis, state_space, *terms, coeffs=None, axes=None):
        """
        Independently builds the dense (n_ho_states, n_ho_states) matrix for
        a bare H.O.-basis operator (no DVR coupling involved at all), using
        only the lower-level `basis.representation(...)` entry point that
        `ContractedDVRHarmonicRepresentation` itself is built on top of.
        """
        rep = basis.representation(*terms, coeffs=coeffs, axes=axes)
        coupled_space = rep.get_transformed_space(state_space, filter_space=state_space)[0]
        return rep.get_representation_matrix(coupled_space, state_space, zero_element_warning=False).asarray()

    def get_full_matrix_reference(self, wfns, basis, state_space, re, v, scales):
        """
        Manual, order-by-order reference for the *full* joint
        (DVR-basis x H.O.-basis) coupling matrix -- built by hand from
        `DVRWavefunctions.expectation` and `dense_operator_matrix`, entirely
        independent of `ContractedDVRHarmonicRepresentation`'s own block
        assembly -- for comparison against `build_total_matrix()`.
        """
        n_states = wfns.wavefunctions.shape[1]
        n_ho_states = len(state_space.indices)
        N = n_states * n_ho_states
        ref = np.zeros((N, N))

        for order, scale in enumerate(scales):
            if scale is None:
                continue
            if order == 0:
                M = scale * wfns.expectation(np.ones_like(wfns.grid))
                op = np.eye(n_ho_states)
            else:
                coord_fn = (lambda grid, o=order: (grid - re) ** o)
                M = scale * wfns.expectation(coord_fn)
                op = self.dense_operator_matrix(
                    basis, state_space, *(['x'] * order), coeffs=self.outer_power(v, order)
                ) / math.factorial(order)

            for n in range(n_states):
                for m in range(n_states):
                    ref[n * n_ho_states:(n + 1) * n_ho_states, m * n_ho_states:(m + 1) * n_ho_states] \
                        += M[n, m] * op

        return ref

    #endregion

    #region tests

    @debugTest
    def test_Construction(self):
        """Basic smoke test: the representation builds, and its sizes/attributes line up."""
        n_states = 3
        n_ho = 2
        max_quanta = 2
        wfns, params = self.get_morse_wavefunctions(n_states=n_states)
        v = np.array([1.0, 0.5])
        coupling_potential = self.get_coupling_potential(params['re'], v, [7e-4, 1.1e-3, 3.3e-4, 1.5e-4])

        rep = ContractedDVRHarmonicRepresentation(
            wfns,
            n_ho=n_ho,
            max_quanta=max_quanta,
            coupling_potential=coupling_potential
        )

        self.assertEqual(rep.n_states, n_states)
        self.assertIsNone(rep.basis_kinetic_coupling)

        n_ho_states = len(rep.state_space.indices)
        mat = rep.build_total_matrix()
        self.assertEqual(mat.shape, (n_states * n_ho_states, n_states * n_ho_states))

    @validationTest
    def test_NoKineticCoupling(self):
        """
        With no `coupling_kinetic` supplied, `basis_kinetic_coupling` stays
        `None` and every block's scalar shift comes entirely from the
        order-0 potential term (kinetic order-0, `['p','p']`, would
        otherwise always contribute a nonzero real piece -- see
        `block_representation`'s docstring/comments).
        """
        n_states = 3
        wfns, params = self.get_morse_wavefunctions(n_states=n_states)
        v = np.array([1.0, 0.5])
        scale0 = 5e-4
        coupling_potential = self.get_coupling_potential(params['re'], v, [scale0, 8e-4, None, None])

        rep = ContractedDVRHarmonicRepresentation(
            wfns, n_ho=2, max_quanta=2, coupling_potential=coupling_potential
        )

        self.assertIsNone(rep.coupling_kinetic)
        self.assertIsNone(rep.basis_kinetic_coupling)

        for n in range(n_states):
            for m in range(n + 1):
                _, shift = rep.block_representation(n, m)
                expected_shift = scale0 if n == m else 0.0
                self.assertAlmostEqual(shift, expected_shift, places=10)

    @validationTest
    def test_TotalMatrixSymmetric(self):
        """The assembled coupling matrix should be exactly (numerically) symmetric."""
        wfns, params = self.get_morse_wavefunctions(n_states=4)
        v = np.array([1.0, 0.5])
        coupling_potential = self.get_coupling_potential(params['re'], v, [7e-4, 1.1e-3, 3.3e-4, 1.5e-4])

        rep = ContractedDVRHarmonicRepresentation(
            wfns, n_ho=2, max_quanta=2, coupling_potential=coupling_potential
        )
        mat = rep.build_total_matrix().toarray()

        self.assertLess(np.max(np.abs(mat - mat.T)), 1e-10)

    @validationTest
    def test_ConstantCouplingIsScalarShift(self):
        """
        A pure order-0 (r-independent) potential coupling has no operator content
        at all: every block should reduce to `scale * delta_{nm} * I`, i.e. the
        whole joint-space matrix should just be `scale * I`.
        """
        n_states = 4
        wfns, params = self.get_morse_wavefunctions(n_states=n_states)
        scale = 0.0021
        coupling_potential = self.get_coupling_potential(
            params['re'], np.array([1.0, 0.5]), [scale, None, None, None]
        )

        rep = ContractedDVRHarmonicRepresentation(
            wfns, n_ho=2, max_quanta=2, coupling_potential=coupling_potential
        )
        mat = rep.build_total_matrix().toarray()
        n_ho_states = len(rep.state_space.indices)
        N = n_states * n_ho_states

        self.assertLess(np.max(np.abs(mat - scale * np.eye(N))), 1e-8)

    @validationTest
    def test_LinearCouplingMatchesSingleModeOperator(self):
        """A pure order-1 coupling should reduce to `<n|scale*(r-re)|m> * (v . x)`."""
        wfns, params = self.get_morse_wavefunctions(n_states=4)
        v = np.array([1.0, 0.5])
        scale = 1.1e-3
        coupling_potential = self.get_coupling_potential(params['re'], v, [None, scale, None, None])

        rep = ContractedDVRHarmonicRepresentation(
            wfns, n_ho=2, max_quanta=2, coupling_potential=coupling_potential
        )
        mat = rep.build_total_matrix().toarray()

        ref = self.get_full_matrix_reference(
            wfns, rep.basis, rep.state_space, params['re'], v, [None, scale, None, None]
        )

        self.assertLess(np.max(np.abs(mat - ref)), 1e-8)

    @validationTest
    def test_QuadraticCouplingMatchesTwoModeOperator(self):
        """A pure order-2 coupling should reduce to `<n|scale*(r-re)^2|m> * (v x v) . xx / 2!`."""
        wfns, params = self.get_morse_wavefunctions(n_states=4)
        v = np.array([1.0, 0.5])
        scale = 3.3e-4
        coupling_potential = self.get_coupling_potential(params['re'], v, [None, None, scale, None])

        rep = ContractedDVRHarmonicRepresentation(
            wfns, n_ho=2, max_quanta=2, coupling_potential=coupling_potential
        )
        mat = rep.build_total_matrix().toarray()

        ref = self.get_full_matrix_reference(
            wfns, rep.basis, rep.state_space, params['re'], v, [None, None, scale, None]
        )

        self.assertLess(np.max(np.abs(mat - ref)), 1e-8)

    @validationTest
    def test_CubicCouplingMatchesThreeModeOperator(self):
        """A pure order-3 coupling should reduce to `<n|scale*(r-re)^3|m> * (v x v x v) . xxx / 3!`."""
        wfns, params = self.get_morse_wavefunctions(n_states=4)
        v = np.array([1.0, 0.5])
        scale = 1.5e-4
        coupling_potential = self.get_coupling_potential(params['re'], v, [None, None, None, scale])

        rep = ContractedDVRHarmonicRepresentation(
            wfns, n_ho=2, max_quanta=2, coupling_potential=coupling_potential
        )
        mat = rep.build_total_matrix().toarray()

        ref = self.get_full_matrix_reference(
            wfns, rep.basis, rep.state_space, params['re'], v, [None, None, None, scale]
        )

        self.assertLess(np.max(np.abs(mat - ref)), 1e-8)

    @validationTest
    def test_FullCubicExpansionIsAdditive(self):
        """
        With all four orders (constant through cubic) populated at once, the
        resulting matrix should just be the sum of what each order gives on
        its own -- checked here against the same order-by-order manual
        reference used in the single-order tests above.
        """
        wfns, params = self.get_morse_wavefunctions(n_states=4)
        v = np.array([1.0, 0.5])
        scales = [7e-4, 1.1e-3, 3.3e-4, 1.5e-4]
        coupling_potential = self.get_coupling_potential(params['re'], v, scales)

        rep = ContractedDVRHarmonicRepresentation(
            wfns, n_ho=2, max_quanta=2, coupling_potential=coupling_potential
        )
        mat = rep.build_total_matrix().toarray()

        ref = self.get_full_matrix_reference(wfns, rep.basis, rep.state_space, params['re'], v, scales)

        self.assertLess(np.max(np.abs(mat - ref)), 1e-8)

    @validationTest
    def test_BlockCachingAndTransposeConsistency(self):
        """
        `get_block` should cache (repeated calls return the identical sparse
        matrix object), and a block computed directly for `(m, n)` should
        agree with the transpose of the independently-computed `(n, m)`
        block -- `build_total_matrix` only ever builds `n >= m` directly and
        fills `n < m` in by transposition, so this checks that shortcut is
        actually valid rather than just assumed.
        """
        wfns, params = self.get_morse_wavefunctions(n_states=4)
        v = np.array([1.0, 0.5])
        coupling_potential = self.get_coupling_potential(params['re'], v, [7e-4, 1.1e-3, 3.3e-4, 1.5e-4])

        rep = ContractedDVRHarmonicRepresentation(
            wfns, n_ho=2, max_quanta=2, coupling_potential=coupling_potential
        )

        b01 = rep.get_block(0, 1)
        b01_again = rep.get_block(0, 1)
        self.assertIs(b01, b01_again)

        b10 = rep.get_block(1, 0)
        self.assertLess(np.max(np.abs(b01.toarray().T - b10.toarray())), 1e-10)

    #endregion
