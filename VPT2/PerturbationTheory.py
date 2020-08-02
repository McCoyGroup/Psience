import numpy as np

from McUtils.Data import UnitsData

from ..Wavefun import Wavefunctions, Wavefunction
from ..Molecools import Molecule

from .Common import PerturbationTheoryException
from .Representations import SubHamiltonian, ProductOperator
from .Terms import PotentialTerms, KineticTerms

__all__ = [
    'PerturbationTheoryWavefunctions',
    'PerturbationTheoryHamiltonian'
]

class PerturbationTheoryWavefunction(Wavefunction):
    @property
    def coeffs(self):
        return self.data

class PerturbationTheoryWavefunctions(Wavefunctions):
    def __init__(self, basis, energies, coeffs, hamiltonian):
        """Energies for the wavefunctions generated

        :param energies:
        :type energies:
        :param coeffs:
        :type coeffs:
        :param hamiltonian:
        :type hamiltonian:
        """
        self.hamiltonian = hamiltonian
        self.basis = basis
        super().__init__(
            energies=energies,
            wavefunctions=coeffs,
            wavefunction_class=PerturbationTheoryWavefunction
        )

class PerturbationTheoryHamiltonian:
    """Represents the main handler used in the perturbation theory calculation
    I probably want it to rely on a Molecule object...

    :param molecule: The molecule we're doing the perturbation theory on
    :type molecule: Molecule
    :param internals: The internal coordinate specification (z-matrix layout)
    :type internals: None | Iterable[Iterable[int]]
    :param numbers_of_quanta: The numbers of quanta of excitation to use for every mode
    :type numbers_of_quanta: np.ndarray | Iterable[int]
    """
    def __init__(self,
                 *ignore,
                 molecule=None,
                 n_quanta=3
                 ):
        if len(ignore) > 0:
            raise PerturbationTheoryException("{} is keyword-only (i.e. takes no positional arguments)".format(
                type(self).__name__
            ))

        if molecule is None:
            raise PerturbationTheoryException("{} requires a Molecule to do its dirty-work")
        self.molecule = molecule

        modes = molecule.normal_modes
        mode_n = modes.basis.matrix.shape[1]
        self.mode_n = mode_n
        self.n_quanta = np.full((mode_n,), n_quanta) if isinstance(n_quanta, (int, np.int)) else tuple(n_quanta)
        # we assume that our modes are already in mass-weighted coordinates, so now we need to make them dimensionless
        # to make life easier, we include this undimensionalization differently for the kinetic and potential energy
        # the kinetic energy needs to be weighted by a sqrt(omega) term while the PE needs a 1/sqrt(omega)
        self.modes = modes

        self.V_terms = PotentialTerms(self.molecule)
        self.G_terms = KineticTerms(self.molecule)

    @classmethod
    def from_fchk(cls, file, internals=None, n_quanta=3):
        """

        :param file:
        :type file: str
        :param internals: internal coordinate specification as a Z-matrix ordering
        :type internals: Iterable[Iterable[int]]
        :param n_quanta:
        :type n_quanta: int | Iterable[int]
        :return:
        :rtype:
        """

        molecule = Molecule.from_file(file, zmatrix=internals, mode='fchk')
        return cls(molecule=molecule, n_quanta=n_quanta)

    @property
    def H0(self):
        def compute_H0(inds,
                       G=self.G_terms[0],
                       V=self.V_terms[0],
                       pp=ProductOperator.pp(self.n_quanta),
                       QQ=ProductOperator.QQ(self.n_quanta),
                       H=self._compute_h0
                       ):
            return H(inds, G, V, pp, QQ)

        return SubHamiltonian(compute_H0, self.n_quanta)

    def _compute_h0(self, inds, G, F, pp, QQ):
        """

        :param inds: which elements of H0 to compute
        :param G: The G-matrix
        :type G: np.ndarray
        :param F: The Hessian
        :type F: np.ndarray
        :param pp: Matrix representation of pp
        :type pp:
        :param QQ: Matrix representation of QQ
        :type QQ:
        :return:
        :rtype:
        """

        # print(type(gmatrix_derivs))
        if not isinstance(G, int):
            # takes a 5-dimensional SparseTensor and turns it into a contracted 2D one
            subKE = pp[inds]
            if isinstance(subKE, np.ndarray):
                ke = np.tensordot(subKE.squeeze(), G, axes=[[0, 1], [0, 1]])
            else:
                ke = subKE.tensordot(G, axes=[[0, 1], [0, 1]]).squeeze()
        else:
            ke = 0

        if not isinstance(F, int):
            subPE = QQ[inds]
            if isinstance(subPE, np.ndarray):
                pe = np.tensordot(subPE.squeeze(), F, axes=[[0, 1], [0, 1]])
            else:
                pe = subPE.tensordot(F, axes=[[0, 1], [0, 1]]).squeeze()

        else:
            pe = 0

        return -1/2*ke + 1/2*pe

    @property
    def H1(self):
        def compute_H1(inds,
                       G=self.G_terms[1],
                       V=self.V_terms[1],
                       pQp=ProductOperator.pQp(self.n_quanta),
                       QQQ=ProductOperator.QQQ(self.n_quanta),
                       H=self._compute_h1
                       ):
            return H(inds, G, V, pQp, QQQ)
        return SubHamiltonian(compute_H1, self.n_quanta)

    def _compute_h1(self, inds, gmatrix_derivs, V_derivs, pQp, QQQ):
        """

        :param inds: which elements of H1 to compute
        :param gmatrix_derivs: The derivatives of the G-matrix with respect to Q
        :type gmatrix_derivs: np.ndarray
        :param V_derivs: The derivatives of V with respect to QQQ
        :type V_derivs: np.ndarray
        :param pQp: Matrix representation of pQp
        :type pQp:
        :param QQQ: Matrix representation of QQQ
        :type QQQ:
        :return:
        :rtype:
        """

        if not isinstance(gmatrix_derivs, int):
            subpQp = pQp[inds]
            if isinstance(subpQp, np.ndarray):
                subpQp = subpQp.squeeze()
                ke = np.tensordot(subpQp, gmatrix_derivs, axes=[[0, 1, 2], [0, 1, 2]])
            else:
                ke = subpQp.tensordot(gmatrix_derivs, axes=[[0, 1, 2], [0, 1, 2]]).squeeze()
        else:
            ke = 0

        if not isinstance(V_derivs, int):
            subQQQ = QQQ[inds]
            if isinstance(subQQQ, np.ndarray):
                subQQQ = subQQQ.squeeze()
                pe = np.tensordot(subQQQ, V_derivs, axes=[[0, 1, 2], [0, 1, 2]])
            else:
                pe = subQQQ.tensordot(V_derivs, axes=[[0, 1, 2], [0, 1, 2]]).squeeze()

        else:
            pe = 0

        # test = UnitsData.convert("Hartrees", "Wavenumbers")/6*subQQQ*V_derivs

        return -1/2*ke + 1/6 * pe

    @property
    def H2(self):
        def compute_H2(inds,
                       G=self.G_terms[2],
                       V=self.V_terms[2],
                       KE=ProductOperator.pQQp(self.n_quanta),
                       PE=ProductOperator.QQQQ(self.n_quanta),
                       H=self._compute_h2
                       ):
            return H(inds, G, V, KE, PE)

        return SubHamiltonian(compute_H2, self.n_quanta)

    def _compute_h2(self, inds, gmatrix_derivs, V_derivs, KE, PE):
        """

        :param inds: which elements of H1 to compute
        :param gmatrix_derivs: The derivatives of the G-matrix with respect to QQ
        :type gmatrix_derivs: np.ndarray
        :param V_derivs: The derivatives of V with respect to QQQQ
        :type V_derivs: np.ndarray
        :param KE: Matrix representation of pQQp
        :type KE:
        :param PE: Matrix representation of QQQQ
        :type PE:
        :return:
        :rtype:
        """

        # print(type(gmatrix_derivs))
        if not isinstance(gmatrix_derivs, int):
            keTens = KE[inds]
            if isinstance(keTens, np.ndarray):
                ke = np.tensordot(keTens.squeeze(), gmatrix_derivs, axes=[[0, 1, 2, 3], [0, 1, 2, 3]])
            else:
                ke = keTens.tensordot(gmatrix_derivs, axes=[[0, 1, 2, 3], [0, 1, 2, 3]]).squeeze()
        else:
            ke = 0

        # print(inds)

        if not isinstance(V_derivs, int):
            peTens = PE[inds]
            if isinstance(peTens, np.ndarray):
                pe = np.tensordot(peTens.squeeze(), V_derivs, axes=[[0, 1, 2, 3], [0, 1, 2, 3]])
            else:
                pe = peTens.tensordot(V_derivs, axes=[[0, 1, 2, 3], [0, 1, 2, 3]]).squeeze()
        else:
            pe = 0

        return -1/4*ke + 1/24*pe

    def get_state_indices(self, states):
        if isinstance(states, (int, np.integer)):
            states = np.arange(min([np.prod(self.n_quanta), states]))
        if not isinstance(states, slice):
            if not isinstance(states[0], (int, np.integer)):
                states = np.ravel_multi_index(np.array(states).T, self.n_quanta)
            if isinstance(states, tuple):  # numpy is weird
                states = np.array(states)
        return states

    def get_state_quantum_numbers(self, states):
        if isinstance(states, slice):
            states = np.arange(np.prod(self.n_quanta))[states]
        elif isinstance(states, int):
            states = np.arange(min([np.prod(self.n_quanta), states]))
        qns = tuple(np.array(np.unravel_index(states, self.n_quanta)).T)
        return qns

    def _get_corrections(self, states=15, coupled_states=None, coeff_threshold=None, energy_threshold=None):
        if states is None:
            states = np.prod(self.n_quanta)
        states = self.get_state_indices(states)
        if coupled_states is None:
            coupled_states = slice(None, None, None)
        if isinstance(coupled_states, slice):
            coupled_states = np.arange(np.prod(self.n_quanta))
        coupled_states = self.get_state_indices(coupled_states)

        H0 = self.H0
        H1 = self.H1
        H2 = self.H2

        energies = H0.diag
        state_E = energies[states]
        energies = energies[coupled_states]
        if isinstance(coupled_states, slice):
            H1_blocks = H1[states, coupled_states]
        else:
            ixes = np.ix_(states, coupled_states)
            H1_blocks = H1[ixes]

        e_blocks = state_E[:, np.newaxis] - np.broadcast_to(energies[np.newaxis], (len(states), len(energies)))

        for n,s in enumerate(states):
            w = np.where(coupled_states==s)[0]
            if len(w) > 0:
                e_blocks[n, w[0]] = 1 # gotta prevent blowups

        if energy_threshold is not None:
            if isinstance(energy_threshold, (int, float, np.integer, np.floating)):
                energy_threshold = (energy_threshold, 1)
            dropped = np.abs(e_blocks) < energy_threshold[0]
            e_blocks[dropped] = np.sign(e_blocks[dropped]) * energy_threshold[1]

        coeffs = H1_blocks / e_blocks

        if coeff_threshold is not None:
            if isinstance(coeff_threshold, (int, float, np.integer, np.floating)):
                coeff_threshold = (coeff_threshold, 0)
            dropped = np.abs(coeffs) > coeff_threshold[0]
            coeffs[dropped] = np.sign(coeffs[dropped]) * coeff_threshold[1]

        for n,s in enumerate(states):
            w = np.where(coupled_states==s)[0]
            if len(w) > 0:
                coeffs[n, w[0]] = 0 # we don't want to double count this

        e_co = np.expand_dims(coeffs, axis=1)
        e_H1 = np.expand_dims(H1_blocks, axis=2)
        e1 = np.matmul(e_co, e_H1).squeeze()
        e2s = H2[states, states]

        # print(self._fmt_corr2_matrix(states, coupled_states, H1_blocks, coeffs))

        return coeffs, (state_E , e1,  e2s), e_blocks, H1_blocks, states, coupled_states

    def _fmt_corr_matrix(self, states, coupled_states, H1_blocks, coeffs):
        """A useful debug function"""
        h1_corr = UnitsData.convert("Hartrees", "Wavenumbers") * H1_blocks * coeffs
        return "\n".join(
            [
                "Corr:",
                "    " + (" ".join(" {2}{1}{0}".format(*x) for x in self.get_state_quantum_numbers(coupled_states))),
                *("{2}{1}{0} ".format(*s) + (" ".join("{:<+4.0f}".format(x) for x in h)) for s, h in
                  zip(self.get_state_quantum_numbers(states), h1_corr))
                ]
        )

    def _fmt_corr2_matrix(self, states, coupled_states, H1_blocks, coeffs):
        """A useful debug function"""
        h1_corr = UnitsData.convert("Hartrees", "Wavenumbers") * H1_blocks
        return "\n".join(
            [
                "Corr:",
                "    " + (" ".join(" {2}{1}{0}".format(*x) for x in self.get_state_quantum_numbers(coupled_states))),
                *("{2}{1}{0} ".format(*s) + (" ".join("{:<+4.0f}".format(x) for x in h)) for s, h in
                  zip(self.get_state_quantum_numbers(states), h1_corr))
                ]
        )

    def get_corrections(self, states=15, coupled_states=None, coeff_threshold=None, energy_threshold=None):
        """

        :param states:
        :type states:
        :param coeff_threshold: a hack for ditching near degeneracies
        :type coeff_threshold: float | Iterable[float]
        :return:
        :rtype:
        """

        return self._get_corrections(states=states, coupled_states=coupled_states,
                                     coeff_threshold=coeff_threshold, energy_threshold=energy_threshold)[:2]

    def get_wavefunctions(self, states=15, coupled_states=None, coeff_threshold=None, energy_threshold=None):
            """Computes perturbation expansion of the wavefunctions and energies

            :param states: the states to target
            :type states: int | iterable[int] | None
            :return: coeffs
            :rtype: np.ndarray
            """

            coeffs, corrs, e_blocks, H1_blocks, states, coupled_states = self._get_corrections(
                states, coeff_threshold=coeff_threshold, coupled_states=coupled_states
            )
            energies = sum(corrs)

            # basis = np.array(np.unravel_index(np.arange(num_ens), self.n_quanta)).T

            return PerturbationTheoryWavefunctions(
                None, # need to get the basis back int, but using the coupled_states args
                # (UnitsData.convert("AtomicUnitOfEnergy", "Wavenumbers")*energies, basis),
                energies,
                coeffs,
                self
            )

    def martin_test(self, states=15, coupled_states=None):
        """Applies the Martin Test to all of the specified states and returns the resulting correlation matrix

        :param states:
        :type states:
        :return:
        :rtype:
        """
        states = self.get_state_indices(states)
        if coupled_states is None:
            coupled_states = states#slice(None, None, None)
        else:
            coupled_states = self.get_state_indices(states)
        if isinstance(coupled_states, slice):
            H1_blocks = self.H1[states, coupled_states]
        else:
            H1_blocks = self.H1[np.ix_(states, coupled_states)]

        energies = self.H0.diag
        state_energies = energies[states]
        coupled_energies = energies[coupled_states]
        diffs = state_energies[:, np.newaxis] - coupled_energies[np.newaxis, :]
        for n,s in enumerate(states):
            w = np.where(coupled_states==s)[0]
            if len(w) > 0:
                diffs[n, w[0]] = 1

        return (H1_blocks**4)/(diffs**3)



