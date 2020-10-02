"""
Provides support for build perturbation theory Hamiltonians
"""

import numpy as np, collections

from ..Molecools import Molecule
from ..BasisReps import HarmonicOscillatorBasis, SimpleProductBasis

from .Common import PerturbationTheoryException
from .Terms import PotentialTerms, KineticTerms
from .Wavefunctions import PerturbationTheoryWavefunction, PerturbationTheoryWavefunctions

__all__ = [
    'PerturbationTheoryHamiltonian'
]

PerturbationTheoryCorrections = collections.namedtuple(
    "PerturbationTheoryCorrections",
    [
        "state",
        "coupled_states",
        "total_states",
        "energies",
        "wavefunctions",
        "hamiltonians"
    ]
)

class PerturbationTheoryHamiltonian:
    """
    Represents the main Hamiltonian used in the perturbation theory calculation.
    Uses a harmonic oscillator basis for representing H0, H1, and H2 (and only goes up to H2 for now)

    :param molecule: The molecule we're doing the perturbation theory on
    :type molecule: Molecule
    :param n_quanta: The numbers of quanta of excitation to use for every mode
    :type n_quanta: int | np.ndarray | Iterable[int]
    """
    def __init__(self, molecule=None, n_quanta=3):

        if molecule is None:
            raise PerturbationTheoryException("{} requires a Molecule to do its dirty-work")
        # molecule = molecule.get_embedded_molecule()
        self.molecule = molecule

        modes = molecule.normal_modes
        mode_n = modes.basis.matrix.shape[1]
        self.mode_n = mode_n
        self.n_quanta = np.full((mode_n,), n_quanta) if isinstance(n_quanta, (int, np.int)) else tuple(n_quanta)
        self.modes = modes

        self.V_terms = PotentialTerms(self.molecule)
        self.G_terms = KineticTerms(self.molecule)
        self._h0 = self._h1 = self._h2 = None

        self.basis = SimpleProductBasis(HarmonicOscillatorBasis, self.n_quanta)

    @classmethod
    def from_fchk(cls, file, internals=None, n_quanta=3):
        """
        :param file: fchk file to load from
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
        """
        Provides the representation for H0 in this basis
        """
        if self._h0 is None:
            self._h0 = (
                -1/2 * self.basis.representation('p', 'p', coeffs=self.G_terms[0])
                + 1/2 * self.basis.representation('x', 'x', coeffs=self.V_terms[0])
            )

        return self._h0

    @property
    def H1(self):
        """
        Provides the representation for H1 in this basis
        """
        if self._h1 is None:
            self._h1 = (
                -1/2 * self.basis.representation('p', 'x', 'p', coeffs=self.G_terms[1], axes=[[0, 1, 2], [1, 0, 2]])
                + 1/6 * self.basis.representation('x', 'x', 'x', coeffs=self.V_terms[1])
            )
        return self._h1

    @property
    def H2(self):
        """
        Provides the representation for H2 in this basis
        """
        if self._h2 is None:
            self._h2 = (
                -1/4 * self.basis.representation('p', 'x', 'x', 'p',
                                                coeffs=self.G_terms[2], axes=[[0, 1, 2, 3], [2, 0, 1, 3]])
                + 1/24 * self.basis.representation('x', 'x', 'x', 'x', coeffs=self.V_terms[2])
            )

        return self._h2

    @property
    def perturbations(self):
        return (self.H0, self.H1, self.H2)

    def _get_state_VPT_corrections(
                self,
                h_reps,
                states,
                coupled_states,
                total_states,
                order
        ):
        """
        Applies perturbation theory to the constructed representations of H0, H1, etc.

        :param h_reps: series of perturbations
        :type h_reps: Iterable[np.ndarray]
        :param states: index of the states to get corrections for
        :type states: Iterable[int]
        :param coupled_states: indices of states to couple when getting corrections
        :type coupled_states: Iterable[int]
        :param total_states: the full number of state indices
        :type total_states: int
        :param order: the order of perturbation theory to apply
        :type order: int
        :return: the vibrational perturbation theory corrections for a single target state
        :rtype: PerturbationTheoryCorrections
        """

        # We use the iterative equations
        #          En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
        #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
        #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=1...k-1) + <n^(0)|n^(k)> |n^(0)>
        # where Pi_n is the perturbation operator [1/(E_m-E_n) for m!=n]

        m = np.array(coupled_states).astype(int) # for safety
        for n in states: # safety again
            if n not in m:
                raise ValueError("state {} must be coupled to itself, but is not in the coupled subspace {}".format(
                    n, m
                ))

        # get explicit matrix reps inside the coupled subspace
        N = len(m)
        # import McUtils.Plots as plt
        wat = h_reps[0][np.ix_(m, m)]
        # plt.ArrayPlot(wat).show()
        H = [h[np.ix_(m, m)].reshape(N, N) for h in h_reps] #type: Iterable[np.ndarray]
        # raise Exception("profiling")

        all_energies = np.zeros((len(states), order + 1))
        all_overlaps = np.zeros((len(states), order + 1))
        all_corrs = np.zeros((len(states), order + 1, N))
        all_wfns = np.zeros((len(states), order + 1, total_states))

        for n, energies, overlaps, corrs, wfns in zip(
                states, all_energies, all_overlaps, all_corrs, all_wfns
        ):
            # taking advantage of mutability of the arrays here...

            # find the state index in the coupled subspace
            n_ind = np.where(m==n)[0][0]
            # generate the perturbation operator
            e_vec = np.diag(H[0])
            E0 = e_vec[n_ind]
            e_vec = e_vec - E0
            e_vec[n_ind] = 1
            pi = 1/e_vec
            pi[n_ind] = 0
            pi = np.diag(pi)
            # raise Exception(pi.shape)

            energies[0] = E0
            overlaps[0] = 1
            corrs[0, n_ind] = 1

            for k in range(1, order+1): # to actually go up to k
                #         En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
                Ek = (
                             H[k][n_ind, n_ind]
                             + sum(np.dot(H[k-i][n_ind], corrs[i]) - energies[k-i]*overlaps[i] for i in range(1, k))
                )
                energies[k] = Ek
                #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
                ok = -1/2 * sum(np.dot(corrs[i], corrs[k-i]) for i in range(1, k))
                overlaps[k] = ok
                #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=0...k-1) + <n^(0)|n^(k)> |n^(0)>
                corrs[k] = sum(np.dot(pi, energies[k-i]*corrs[i] - np.dot(H[k-i], corrs[i])) for i in range(0, k))
                corrs[k][n_ind] = ok # pi (the perturbation operator) ensures it's zero before this

            # now we broadcast the corrections back up so that they're good for the _entire_ population of states...
            for wfn, cor in zip(wfns, corrs):
                if isinstance(coupled_states[0], (int, np.integer)):
                    cs = (coupled_states,)
                else:
                    cs = coupled_states
                wfn[cs] = cor

        return PerturbationTheoryCorrections(
            states,
            coupled_states,
            total_states,
            all_energies,
            all_wfns,
            H
        )

    def get_wavefunctions(self, states, coupled_states=None, order=2):
        """
        Gets a set of `PerturbationTheoryWavefunctions` from the perturbations defined by the Hamiltonian

        :param states: the states to get the index for, given either as indices or as a numbers of quanta
        :type states: Iterable[int] | Iterable[Iterable[int]]
        :param coupled_states: the list of states to explicitly allow to couple in
        :type coupled_states: Iterable[int] | Iterable[Iterable[int]]
        :return: generated wave functions
        :rtype: PerturbationTheoryWavefunctions
        """

        if not isinstance(states[0], (int, np.integer)):
            states = self.basis.ravel_state_inds(states)
        total_states = int(np.prod(self.basis.quanta))
        if coupled_states is None:
            coupled_states = np.arange(total_states) # expensive but...I guess what are you gonna do?
        elif not isinstance(coupled_states[0], int):
            coupled_states = self.basis.ravel_state_inds(coupled_states)

        corrs = self._get_state_VPT_corrections(
            self.perturbations,
            states,
            coupled_states,
            total_states,
            order
            )

        return PerturbationTheoryWavefunctions(self.molecule, self.basis, corrs)

    # def martin_test(self, states=15, coupled_states=None):
    #     """Applies the Martin Test to all of the specified states and returns the resulting correlation matrix
    #
    #     :param states:
    #     :type states:
    #     :return:
    #     :rtype:
    #     """
    #     states = self.get_state_indices(states)
    #     if coupled_states is None:
    #         coupled_states = states#slice(None, None, None)
    #     else:
    #         coupled_states = self.get_state_indices(states)
    #     if isinstance(coupled_states, slice):
    #         H1_blocks = self.H1[states, coupled_states]
    #     else:
    #         H1_blocks = self.H1[np.ix_(states, coupled_states)]
    #
    #     energies = self.H0.diag
    #     state_energies = energies[states]
    #     coupled_energies = energies[coupled_states]
    #     diffs = state_energies[:, np.newaxis] - coupled_energies[np.newaxis, :]
    #     for n,s in enumerate(states):
    #         w = np.where(coupled_states==s)[0]
    #         if len(w) > 0:
    #             diffs[n, w[0]] = 1
    #
    #     return (H1_blocks**4)/(diffs**3)