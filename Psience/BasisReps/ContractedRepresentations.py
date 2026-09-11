"""
Provides `ContractedBasisHarmonicRepresentation` and friends, which couple
an external, already-indexed ("contracted") basis -- e.g. the eigenstates
of a 1D DVR -- to a `HarmonicOscillatorProductBasis`, following the same
general term-construction pattern VPT2 uses to build its potential (V) and
kinetic (G-matrix) term representations (see `Psience.VPT2.Hamiltonian`):

  potential-like term of order o:  ['x'] * o                (o factors of Q)
  kinetic-like term of order o:    ['p'] + ['x']*o + ['p']   (p Q...Q p)
"""

import math
import numpy as np
import scipy.sparse as sp

from .HarmonicOscillator import HarmonicOscillatorProductBasis
from .StateSpaces import BasisStateSpace

__all__ = [
    "ContractedBasisHarmonicRepresentation",
    "DVREigenbasisCouplingExpansion",
    "ContractedDVRHarmonicRepresentation"
]


# ===========================================================================
# ContractedBasisHarmonicRepresentation
# ===========================================================================

class ContractedBasisHarmonicRepresentation:
    """
    Couples an external "contracted" basis (indices n = 0 ... n_states-1;
    e.g. the eigenstates of a 1D DVR, but the class doesn't care what they
    actually are) to a `HarmonicOscillatorProductBasis`.

    For every pair (n, m) of contracted-basis indices, `basis_potential_coupling`
    and `basis_kinetic_coupling` each supply a `coupling_expansion`: a list
    of arrays, ordered 0, 1, ..., k, playing the role of successive Taylor
    orders in the H.O. coordinates:

        potential order o -> coefficient tensor of shape (n_ho,)*o,
                              contracted against o factors of 'x'
        kinetic   order o -> coefficient tensor of shape (n_ho,)*(o+2),
                              contracted against ['p'] + ['x']*o + ['p']
                              (a VPT2-style p Q...Q p sandwich)

    `basis_potential_coupling(n, m)` generalizes a bare linear
    `coupling_vector` (that's the order-1-only special case);
    `basis_kinetic_coupling(n, m)` is the analogous kinetic-energy-like
    expansion, following the same p/x term pattern used in
    `Psience.VPT2.Hamiltonian`.
    """

    def __init__(self,
                 n_states,
                 basis=None,
                 n_ho=None,
                 max_quanta=2,
                 state_space=None,
                 basis_potential_coupling=None,
                 basis_kinetic_coupling=None,
                 operator_settings=None
                 ):
        """
        :param n_states: size of the external contracted basis (e.g. number
            of DVR eigenstates retained)
        :param basis: a `HarmonicOscillatorProductBasis`; built from `n_ho`
            if not supplied
        :param n_ho: number of H.O. modes (used only if `basis` is None)
        :param max_quanta: total-quanta cutoff used to build `state_space`
            if it isn't supplied directly
        :param state_space: a `BasisStateSpace` over `basis`; built via
            `BasisStateSpace.from_quanta(basis, range(max_quanta+1))` if
            not supplied
        :param basis_potential_coupling: callable (n, m) -> coupling_expansion
            (list of arrays, order 0..k) for potential-like ('x'*o) terms
        :param basis_kinetic_coupling: callable (n, m) -> coupling_expansion
            (list of arrays, order 0..k) for kinetic-like
            ('p' + 'x'*o + 'p') terms
        :param operator_settings: extra kwargs forwarded to
            `basis.representation(...)`
        """
        if basis is None:
            if n_ho is None:
                raise ValueError("need either `basis` or `n_ho`")
            basis = HarmonicOscillatorProductBasis(n_ho)
        self.basis = basis

        if state_space is None:
            state_space = BasisStateSpace.from_quanta(self.basis, list(range(max_quanta + 1)))
        self.state_space = state_space
        self.max_quanta = max_quanta

        self.n_states = n_states
        self.basis_potential_coupling = basis_potential_coupling
        self.basis_kinetic_coupling = basis_kinetic_coupling
        self.operator_settings = {} if operator_settings is None else operator_settings

        self._coupled_space = None
        self._block_cache = {}

    # -- building a single order's `Representation` -------------------------

    def _potential_term_representation(self, order, coeffs, name=None):
        """order o -> o factors of 'x', axes = direct correspondence"""
        terms = ["x"] * order
        axes = [list(range(order)), list(range(order))]
        norm = 1.0 / math.factorial(order)
        return norm * self.basis.representation(
            *terms, coeffs=coeffs, axes=axes, name=name, **self.operator_settings
        )

    def _kinetic_term_representation(self, order, coeffs, name=None):
        """order o -> ['p'] + ['x']*o + ['p'], VPT2-style p Q...Q p sandwich"""
        terms = ["p"] + ["x"] * order + ["p"]
        naxes = order + 2
        axes = [list(range(naxes)), list(range(naxes))]
        norm = 1.0 / (2 * math.factorial(order)) if order > 0 else 0.5
        return norm * self.basis.representation(
            *terms, coeffs=coeffs, axes=axes, name=name, **self.operator_settings
        )

    # -- combining a whole coupling_expansion into one Representation ------

    def _expansion_representation(self, coupling_expansion, kind, n, m):
        """
        Sums the `Representation`s for every nonzero order in
        `coupling_expansion`. Potential order 0 (a bare constant -- no
        operators at all) can't be built as a `Representation` (there's no
        0-operator term to hand `basis.representation`), so it's split off
        and returned separately as a scalar to be added as
        `scalar * identity` when the final block matrix is realized.

        :return: (rep_or_None, scalar_shift)
        """
        rep = None
        scalar_shift = 0.0
        for order, coeffs in enumerate(coupling_expansion):
            if coeffs is None:
                continue
            coeffs = np.asanyarray(coeffs, dtype=float)
            if not np.any(coeffs):
                continue

            if kind == "potential" and order == 0:
                scalar_shift += float(coeffs)
                continue

            name = "{}[{},{}]({})".format("V" if kind == "potential" else "T", n, m, order)
            if kind == "potential":
                term_rep = self._potential_term_representation(order, coeffs, name=name)
            else:
                term_rep = self._kinetic_term_representation(order, coeffs, name=name)
            rep = term_rep if rep is None else rep + term_rep

        return rep, scalar_shift

    def block_representation(self, n, m):
        """
        Combined potential + kinetic `Representation` for contracted-basis
        indices (n, m) (`None` if there's no operator-bearing piece), plus
        any scalar (order-0 potential) shift.

        :return: (rep_or_None, scalar_shift)
        """
        pieces = []
        scalar_shift = 0.0

        if self.basis_potential_coupling is not None:
            v_expansion = self.basis_potential_coupling(n, m)
            if v_expansion is not None:
                v_rep, v_shift = self._expansion_representation(v_expansion, "potential", n, m)
                if v_rep is not None:
                    pieces.append(v_rep)
                scalar_shift += v_shift

        if self.basis_kinetic_coupling is not None:
            t_expansion = self.basis_kinetic_coupling(n, m)
            if t_expansion is not None:
                t_rep, t_shift = self._expansion_representation(t_expansion, "kinetic", n, m)
                if t_rep is not None:
                    pieces.append(t_rep)
                scalar_shift += t_shift  # kinetic order 0 is ['p','p'], always a real rep, so this stays 0

        rep = None
        if len(pieces) > 0:
            rep = pieces[0]
            for p in pieces[1:]:
                rep = rep + p

        return rep, scalar_shift

    # -- realizing blocks and the total matrix ------------------------------

    def _matrix_from_representation(self, rep):
        if self._coupled_space is None:
            self._coupled_space = rep.get_transformed_space(
                self.state_space, filter_space=self.state_space
            )[0]
        mat = rep.get_representation_matrix(
            self._coupled_space, self.state_space, zero_element_warning=False
        )
        return mat.asarray()

    def get_block(self, n, m):
        """
        The sparse (n_ho_states x n_ho_states) matrix block coupling
        contracted-basis states n and m (cached).
        """
        key = (n, m)
        if key not in self._block_cache:
            n_ho_states = len(self.state_space.indices)
            rep, scalar_shift = self.block_representation(n, m)
            if rep is None:
                mat = np.zeros((n_ho_states, n_ho_states))
            else:
                mat = self._matrix_from_representation(rep)
            if scalar_shift != 0.0:
                mat = mat + scalar_shift * np.eye(n_ho_states)
            self._block_cache[key] = sp.csr_matrix(mat)
        return self._block_cache[key]

    def build_total_matrix(self):
        """
        Assembles the full block-sparse matrix over the joint
        (contracted-basis) (x) (harmonic-oscillator) space. Only blocks
        n >= m are ever built directly; n < m blocks are filled in by
        transposition.
        """
        n_states = self.n_states
        blocks = [[None] * n_states for _ in range(n_states)]
        for n in range(n_states):
            for m in range(n + 1):
                block = self.get_block(n, m)
                blocks[n][m] = block
                if n != m:
                    blocks[m][n] = block.T
        return sp.bmat(blocks, format='csr')


# ===========================================================================
# DVREigenbasisCouplingExpansion
# ===========================================================================

class DVREigenbasisCouplingExpansion:
    """
    Builds `basis_potential_coupling`/`basis_kinetic_coupling` callables
    (suitable for `ContractedBasisHarmonicRepresentation`) from a proper
    `DVRWavefunctions` object and a set of explicitly-ordered terms, each a
    scalar function of the DVR coordinate together with a fixed H.O.-side
    coupling tensor:

        order o contributes   c_o(r) * K_o
    to the coupling_expansion, where `c_o` is a function of the DVR
    coordinate alone and `K_o` is a constant tensor of the appropriate
    shape (a scalar for order 0, a length-n_ho vector for order 1, a
    matrix for order 2, etc). The (n, m) matrix element of c_o(r) is
    obtained through the wavefunctions' own matrix-element machinery,
    `DVRWavefunctions.expectation`, rather than by manually contracting a
    bare eigenvector array:

        <n| c_o(r) |m> = wavefunctions.expectation(c_o)

    (`expectation` accepts `c_o` as a callable -- evaluating it on
    `wavefunctions.grid` internally -- and returns the full (n_states,
    n_states) matrix of <n|c_o(r)|m> values at once). Each order's matrix
    elements are computed once and cached, so repeated (n, m) queries are
    cheap. Terms are supplied as explicit `(order, coord_function, tensor)`
    triples (rather than assuming list-position == order) so that, e.g., a
    purely order-1 (linear) coupling can be specified without also having to
    supply a dummy order-0 term.
    """

    def __init__(self,
                 wavefunctions,
                 potential_terms=None,
                 kinetic_terms=None
                 ):
        """
        :param wavefunctions: a `DVRWavefunctions` object; matrix elements
            are obtained via `wavefunctions.expectation(coord_function)`
        :type wavefunctions: DVRWavefunctions
        :param potential_terms: list of (order, coord_function, tensor) triples
            for potential-like ('x'*order) terms; `coord_function` maps the
            grid -> array(n_grid,) (or any grid-resolved tensor `expectation`
            accepts), `tensor` has shape (n_ho,)*order (a plain scalar for
            order 0)
        :param kinetic_terms: same idea, for kinetic-like ('p'+'x'*order+'p')
            terms; `tensor` has shape (n_ho,)*(order + 2)
        """
        self.wavefunctions = wavefunctions

        self.potential_terms = list(potential_terms or [])
        self.kinetic_terms = list(kinetic_terms or [])

        self._potential_matrix_elements = {}
        self._kinetic_matrix_elements = {}

    def _matrix_elements(self, coord_function):
        """
        <n| c(r) |m> for every (n, m) at once, via the wavefunctions' own
        expectation-value machinery.
        """
        return self.wavefunctions.expectation(coord_function)

    def _get_expansion(self, terms, cache, n, m):
        if len(terms) == 0:
            return None
        max_order = max(order for order, _, _ in terms)
        expansion = [None] * (max_order + 1)
        for order, coord_function, tensor in terms:
            if order not in cache:
                cache[order] = self._matrix_elements(coord_function)
            contribution = cache[order][n, m] * np.asanyarray(tensor, dtype=float)
            expansion[order] = (
                contribution if expansion[order] is None else expansion[order] + contribution
            )
        return expansion

    def basis_potential_coupling(self, n, m):
        return self._get_expansion(self.potential_terms, self._potential_matrix_elements, n, m)

    def basis_kinetic_coupling(self, n, m):
        return self._get_expansion(self.kinetic_terms, self._kinetic_matrix_elements, n, m)


# ===========================================================================
# ContractedDVRHarmonicRepresentation
# ===========================================================================

class ContractedDVRHarmonicRepresentation(ContractedBasisHarmonicRepresentation):
    """
    A `ContractedBasisHarmonicRepresentation` whose "contracted basis" is
    literally the eigenstates of a 1D DVR: it stores the `DVRWavefunctions`
    object directly (rather than a bare grid + eigenvector array) and a
    harmonic state space, and builds `basis_potential_coupling` /
    `basis_kinetic_coupling` automatically from a `coupling_potential`
    (and, optionally, `coupling_kinetic`) function.

    `coupling_potential(grid)` is called once, on the whole DVR grid, and
    must return a list of arrays -- one per expansion order
    o = 0, 1, ..., k -- of shape (n_grid,) + (n_ho,)*o. This is more
    general than `DVREigenbasisCouplingExpansion`: there's no assumption
    that the order-o term factors into a scalar coordinate function times
    a fixed H.O.-side tensor; it can be an arbitrary grid-point-resolved
    tensor. `coupling_kinetic(grid)` works the same way, but returns
    tensors of shape (n_grid,) + (n_ho,)*(o+2), for the p Q...Q p sandwich
    terms.

    Every element of the resulting coupling_expansion is obtained by
    literally averaging the grid-resolved tensor over the DVR
    wavefunctions:

        M_o[n, m] = sum_k  psi_n(r_k) * T_o(r_k) * psi_m(r_k)

    (`np.einsum('kn,k...,km->nm...', wfn_data, T_o, wfn_data)`), computed
    once per order and cached; `basis_potential_coupling(n, m)`/
    `basis_kinetic_coupling(n, m)` then just index into that cache.
    """

    def __init__(self,
                 wavefunctions,
                 state_space=None,
                 basis=None,
                 n_ho=None,
                 max_quanta=2,
                 coupling_potential=None,
                 coupling_kinetic=None,
                 operator_settings=None
                 ):
        """
        :param wavefunctions: a `DVRWavefunctions` object (has `.grid` and
            `.wavefunctions`, the (n_grid, n_states) eigenvector matrix)
        :param state_space: harmonic `BasisStateSpace`; built from
            `max_quanta` if not supplied
        :param basis: a `HarmonicOscillatorProductBasis`; built from
            `n_ho` (or inferred from the coupling functions) if not supplied
        :param n_ho: number of H.O. modes, if not inferable/supplied via `basis`
        :param max_quanta: total-quanta cutoff for the default `state_space`
        :param coupling_potential: callable grid -> list of arrays, order
            0..k, shapes (n_grid,) + (n_ho,)*o
        :param coupling_kinetic: callable grid -> list of arrays, order
            0..k, shapes (n_grid,) + (n_ho,)*(o+2)
        :param operator_settings: extra kwargs forwarded to
            `basis.representation(...)`
        """
        self.wavefunctions = wavefunctions
        self.grid = wavefunctions.grid
        n_states = wavefunctions.wavefunctions.shape[1]

        self.coupling_potential = coupling_potential
        self.coupling_kinetic = coupling_kinetic
        self._potential_matrix_elements = None
        self._kinetic_matrix_elements = None

        if basis is None and n_ho is None:
            n_ho = self._infer_n_ho()

        super().__init__(
            n_states=n_states,
            basis=basis,
            n_ho=n_ho,
            max_quanta=max_quanta,
            state_space=state_space,
            basis_potential_coupling=(
                self._basis_potential_coupling if coupling_potential is not None else None
            ),
            basis_kinetic_coupling=(
                self._basis_kinetic_coupling if coupling_kinetic is not None else None
            ),
            operator_settings=operator_settings
        )

    def _infer_n_ho(self):
        for coupling, extra in ((self.coupling_potential, 0), (self.coupling_kinetic, 2)):
            if coupling is None:
                continue
            for T in coupling(self.grid):
                if T is not None and np.ndim(T) > 1:
                    return np.asanyarray(T).shape[1]
        raise ValueError(
            "couldn't infer `n_ho` from the coupling functions; please supply `basis` or `n_ho` explicitly"
        )

    def _average_over_wavefunctions(self, tensor_stack):
        """
        (n_grid,) + extra_dims  ->  (n_states, n_states) + extra_dims,
        via the wavefunctions' own `expectation` matrix-element machinery
        (it accepts a precomputed grid-resolved array directly, so there's
        no need to re-derive the contraction by hand).
        """
        return self.wavefunctions.expectation(np.asanyarray(tensor_stack, dtype=float))

    def _get_potential_matrix_elements(self):
        if self._potential_matrix_elements is None:
            terms = self.coupling_potential(self.grid)
            self._potential_matrix_elements = [
                self._average_over_wavefunctions(T) if T is not None else None
                for T in terms
            ]
        return self._potential_matrix_elements

    def _get_kinetic_matrix_elements(self):
        if self._kinetic_matrix_elements is None:
            terms = self.coupling_kinetic(self.grid)
            self._kinetic_matrix_elements = [
                self._average_over_wavefunctions(T) if T is not None else None
                for T in terms
            ]
        return self._kinetic_matrix_elements

    def _basis_potential_coupling(self, n, m):
        return [mat[n, m] if mat is not None else None for mat in self._get_potential_matrix_elements()]

    def _basis_kinetic_coupling(self, n, m):
        return [mat[n, m] if mat is not None else None for mat in self._get_kinetic_matrix_elements()]
