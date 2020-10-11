"""
Provides support for build perturbation theory Hamiltonians
"""

import numpy as np, itertools

from ..Molecools import Molecule
from ..BasisReps import HarmonicOscillatorBasis, SimpleProductBasis

from .Common import PerturbationTheoryException
from .Terms import PotentialTerms, KineticTerms, CoriolisTerm, PotentialLikeTerm

__all__ = [
    'PerturbationTheoryHamiltonian',
    'PerturbationTheoryCorrections'
]


class PerturbationTheoryCorrections:
    """
    Represents a set of corrections from perturbation theory.
    Can be used to correct other operators in the basis of the original calculation.

    """
    def __init__(self,
                 states,
                 corrections,
                 hamiltonians
                 ):
        """
        :param states: a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
        :type states: dict
        :param corrections: the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT
        :type corrections: dict
        :param hamiltonians: the set of Hamiltonian matrices used as an expansion
        :type hamiltonians: Iterable[np.ndarray]
        """
        self.states = states['states']
        self.coupled_states = states['coupled_states']
        self.total_basis = states['total_states']
        self.energy_corrs = corrections['energies']
        self.wfn_corrections = corrections['wavefunctions']
        if 'degenerate_states' in states:
            self.degenerate_states = states['degenerate_states']
        else:
            self.degenerate_states = None

        if 'degenerate_transformation' in corrections:
            self.degenerate_transf = corrections['degenerate_transformation']
        else:
            self.degenerate_transf = None

        if 'degenerate_energies' in corrections:
            self.degenerate_energies = corrections['degenerate_energies']
        else:
            self.degenerate_energies = None
        self.hams = hamiltonians

    @property
    def degenerate(self):
        return self.degenerate_transf is not None

    @property
    def energies(self):
        if self.degenerate:
            return self.degenerate_energies
        else:
            return np.sum(self.energy_corrs, axis=1)

    @property
    def order(self):
        return len(self.energy_corrs[0])

    def operator_representation(self, operator_expansion, order=None, subspace=None):
        """
        Generates the representation of the operator in the basis of stored states

        :param operator_expansion: the expansion of the operator
        :type operator_expansion: Iterable[float] | Iterable[np.ndarray]
        :param order: the order of correction to go up to
        :type order: Iterable[float] | Iterable[np.ndarray]
        :param subspace: the subspace of terms in which the operator expansion is defined
        :type subspace: Iterable[int]
        :return: the set of representation matrices for this operator
        :rtype: Iterable[np.ndarray]
        """

        mordor = self.order - 1
        if order is None:
            order = mordor
        if order > mordor:
            raise PerturbationTheoryException("{}: can't correct up to order {} when zero-order states were only corrected up to order {}".format(
                type(self).__name__,
                order,
                mordor
            ))
        order = order + 1 # so that we actually do get up to the request order after accounting for the zeros...
        if len(operator_expansion) < order:
            operator_expansion = list(operator_expansion) + [0]*(order - len(operator_expansion))

        wfn_corrs = [
            self.wfn_corrections[:, k, subspace] if subspace is not None else self.wfn_corrections[:, k, subspace]
            for k in range(order)
        ]

        reps = [np.zeros(1)] * order
        for k in range(order):
            op = None
            for a in range(k+1):
                for b in range(k-a+1):
                    c = k - (a + b)
                    rop = operator_expansion[c]
                    if isinstance(rop, (int, float, np.integer, np.floating)): # constant reps...
                        if rop != 0: # cheap easy check
                            subrep = rop * np.dot(wfn_corrs[a], wfn_corrs[b].T)
                            if op is None:
                                op = subrep
                            else:
                                op += subrep
                    else:
                        subrep = np.dot(np.dot(wfn_corrs[a], rop), wfn_corrs[b].T)
                        if op is None:
                            op = subrep
                        else:
                            op += subrep
            reps[k] = op

        return reps

    def savez(self, file):
        keys = dict(
            states=self.states,
            coupled_states=self.coupled_states,
            total_states=self.total_basis,
            energies=self.energy_corrs,
            wavefunctions=self.wfn_corrections,
            hamiltonians=self.hams
        )
        if self.degenerate_states is not None:
            keys['degenerate_states'] = self.degenerate_states
        if self.degenerate_transf is not None:
            keys['degenerate_transformation'] = self.degenerate_transf
        if self.degenerate_energies is not None:
            keys['degenerate_energies'] = self.degenerate_energies
        np.savez(file, **keys)

    @classmethod
    def loadz(cls, file):
        keys = np.load(file)
        return cls(
            {
                "states":keys['states'],
                 "coupled_states":keys['coupled_states'],
                 "total_states":keys['total_states'],
                 "degenerate_states":keys['degenerate_states'] if 'degenerate_states' in keys else None
             },
            {
                "energies":keys['energies'],
                "wavefunctions":keys['wavefunctions'],
                "degenerate_transformation": keys['degenerate_transformation'] if 'degenerate_transformation' in keys else None,
                "degenerate_energies": keys['degenerate_energies'] if 'degenerate_energies' in keys else None
            },
            keys['hamiltonians']
        )


# PerturbationTheoryCorrections = collections.namedtuple(
#     "PerturbationTheoryCorrections",
#     [
#         "energies",
#         "wavefunctions",
#         "states",
#         "corrections",
#         "hamiltonians"
#     ]
# )

class PerturbationTheoryHamiltonian:
    """
    Represents the main Hamiltonian used in the perturbation theory calculation.
    Uses a harmonic oscillator basis for representing H0, H1, and H2 (and only goes up to H2 for now)
    """
    def __init__(self,
                 molecule=None,
                 n_quanta=3,
                 modes=None,
                 mode_selection=None,
                 coriolis_coupling = True
                 ):
        """
        :param molecule: the molecule on which we're doing perturbation theory
        :type molecule:  Molecule
        :param n_quanta: the numbers of quanta to use when building representations
        :type n_quanta: Iterable[int]
        :param modes: the set of modes to use as the basis
        :type modes: None | MolecularNormalModes
        :param mode_selection: the subset of modes to use when doing expansions
        :type mode_selection: None | Iterable[int]
        :param coriolis_coupling: whether to add coriolis coupling if not in internals
        :type coriolis_coupling: bool
        """
        if molecule is None:
            raise PerturbationTheoryException("{} requires a Molecule to do its dirty-work")
        # molecule = molecule.get_embedded_molecule()
        self.molecule = molecule
        if modes is None:
            modes = molecule.normal_modes
        mode_n = modes.basis.matrix.shape[1] if mode_selection is None else len(mode_selection)
        self.mode_n = mode_n
        self.n_quanta = np.full((mode_n,), n_quanta) if isinstance(n_quanta, (int, np.int)) else tuple(n_quanta)
        self.modes = modes

        self.V_terms = PotentialTerms(self.molecule, modes=modes, mode_selection=mode_selection)
        self.G_terms = KineticTerms(self.molecule, modes=modes, mode_selection=mode_selection)
        if coriolis_coupling and (self.molecule.internal_coordinates is None):
            self.coriolis_terms = CoriolisTerm(self.molecule, modes=modes, mode_selection=mode_selection)
        else:
            self.coriolis_terms = None
        self.watson_term = PotentialLikeTerm(self.molecule, modes=modes, mode_selection=mode_selection)

        self._h0 = self._h1 = self._h2 = None

        self.basis = SimpleProductBasis(HarmonicOscillatorBasis, self.n_quanta)

    @classmethod
    def from_fchk(cls, file,
                  internals=None,
                  n_quanta=3,
                  mode_selection=None
                  ):
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
        return cls(molecule=molecule, n_quanta=n_quanta, mode_selection=mode_selection)

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
            if self.coriolis_terms is not None:
                total_cor = self.coriolis_terms[0]+self.coriolis_terms[1]+self.coriolis_terms[2]
                # import McUtils.Plots as plt
                # plt.ArrayPlot(
                #     total_cor.reshape(total_cor.shape[0] ** 2, total_cor.shape[0] ** 2)
                # ).show()
                self._h2 += -1 * self.basis.representation('x', 'p', 'x', 'p', coeffs=total_cor)
            else:
                self._h2 += 0 * self.basis.representation(coeffs=0)

            self._h2 += 1/8 * self.basis.representation(coeffs=self.watson_term[0])

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
                order,
                degenerate_states = None
        ):
        """
        Applies perturbation theory to the constructed representations of H0, H1, etc.

        :param h_reps: series of perturbations as indexable objects
        :type h_reps: Iterable[np.ndarray]
        :param states: index of the states to get corrections for
        :type states: Iterable[int]
        :param coupled_states: indices of states to couple when getting corrections
        :type coupled_states: Iterable[int]
        :param total_states: the full number of state indices
        :type total_states: int
        :param order: the order of perturbation theory to apply
        :type order: int
        :param degenerate_states: the pairs of degeneracies to account for
        :type degenerate_states: None | (Iterable[int], Iterable[int])
        :return: the vibrational perturbation theory corrections for a single target state
        :rtype: PerturbationTheoryCorrections
        """

        # We use the iterative equations
        #          En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
        #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
        #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=1...k-1) + <n^(0)|n^(k)> |n^(0)>
        # where Pi_n is the perturbation operator [1/(E_m-E_n) for m!=n]

        m = np.array(coupled_states).astype(int) # for safety
        state_inds = [-1]*len(states)
        for i,n in enumerate(states):
            si = np.where(m == n)[0]
            if len(si) == 0:
                raise ValueError("requested state {} must be coupled to itself, but is not in the coupled subspace {}".format(
                    n, m
                ))
            elif len(si) > 1:
                raise ValueError("requested state {} appears multiple times in the coupled subspace {}".format(
                    n, m
                ))
            state_inds[i] = si[0]

        # get explicit matrix reps inside the coupled subspace
        N = len(m)
        # import McUtils.Plots as plt
        # wat = h_reps[1][np.ix_(m, m)]
        # plt.ArrayPlot(self.G_terms[0])
        # plt.ArrayPlot(wat).show()
        # raise Exception("...wat")

        H = [np.zeros(1)] * len(h_reps)
        for i,h in enumerate(h_reps):
            sub = h[np.ix_(m, m)]
            if isinstance(sub, (int, np.integer, np.floating, float)):
                sub = np.full((N, N), sub)
            H[i] = sub #type: np.ndarray

        # import McUtils.Plots as plt
        # plt.ArrayPlot(H[2]).show()

        if degenerate_states is not None:
            # we check this twice because the Martin test can return None
            if isinstance(degenerate_states, (int, np.integer, np.floating, float)):
                thresh = degenerate_states
                degenerate_states = self._martin_test(
                    H,
                    state_inds, # state indices in the coupled_states
                    thresh
                )
                # print("Got {} from the Martin test wth threshold {}".format(
                #     degenerate_states,
                #     thresh
                # ))
            elif all(isinstance(x, (int, np.integer)) for x in degenerate_states):
                # means we just got the degenerate subspace to work with
                pairs = [
                    l for l in itertools.product(degenerate_states, degenerate_states) if l[0] < l[1]
                    ]
                degenerate_states = np.array(pairs).T
                print('Got degenerate states {}'.format(pairs))
            elif degenerate_states is not None:
                try:
                    degenerate_states = degenerate_states(H, coupled_states)
                except (TypeError, ValueError):
                    pass

        if degenerate_states is None:
            deg_vals = None
        else:
            deg_i, deg_j = degenerate_states
            deg_vals = np.array([h[deg_i, deg_j] for h in H])

            deg_i, deg_j = degenerate_states
            H_non_deg = [np.zeros(1)] * len(H)
            for i,h in enumerate(H):
                h = h.copy()
                h[deg_i, deg_j] = 0.
                h[deg_j, deg_i] = 0.
                H_non_deg[i] = h
            H = H_non_deg

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

            energies[0] = E0
            overlaps[0] = 1
            corrs[0, n_ind] = 1

            for k in range(1, order+1): # to actually go up to k
                #         En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
                Ek = (
                        H[k][n_ind, n_ind]
                        + sum(
                             np.dot(H[k-i][n_ind], corrs[i])
                            - energies[k-i]*overlaps[i]
                            for i in range(1, k)
                        )
                )
                energies[k] = Ek
                #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
                ok = -1/2 * sum(np.dot(corrs[i], corrs[k-i]) for i in range(1, k))
                overlaps[k] = ok
                #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=0...k-1) + <n^(0)|n^(k)> |n^(0)>
                corrs[k] = sum(
                    np.dot(pi, energies[k-i]*corrs[i] - np.dot(H[k-i], corrs[i]))
                    for i in range(0, k)
                )
                corrs[k][n_ind] = ok # pi (the perturbation operator) ensures it's zero before this

            # now we broadcast the corrections back up so that they're good for the _entire_ population of states...
            for wfn, cor in zip(wfns, corrs):
                if isinstance(coupled_states[0], (int, np.integer)):
                    cs = (coupled_states,)
                else:
                    cs = coupled_states
                wfn[cs] = cor

        corrs = PerturbationTheoryCorrections(
            {
                "states":states,
                 "coupled_states":coupled_states,
                 "total_states":total_states,
                 "degenerate_states":degenerate_states
             },
            {
                "energies":all_energies,
                "wavefunctions":all_wfns,
                "degenerate_transformation":None,

                "degenerate_energies": None
            },
            H
        )

        if degenerate_states is not None:
            # means we need to do the second level of corrections
            # we're going to need to do a full diagonalization, but we handle this
            # by pulling only the coupled elements down to a dense matrix and
            # diagonalizing there before returning the degenerate rotation as a sparse
            # matrix

            # first we figure out what the space of degenerate states is
            # has to be an NxN, so we figure out what the total space of potential degeneracies is
            degenerate_space = np.unique(np.concatenate(degenerate_states))
            deg_dim = len(degenerate_space)
            # then we have to work out how indices in the larger space map onto those in the smaller one
            remap = {i:k for k,i in enumerate(degenerate_space)}
            mapped_inds = (
                np.array([remap[i] for i in degenerate_states[0]]),
                np.array([remap[j] for j in degenerate_states[1]])
            )
            # at this point, we can fill H_deg in the basis of degenerate states
            H_deg = [np.zeros(1)] * len(deg_vals[1:])
            for i, v in enumerate(deg_vals[1:]):
                H_deg[i] = np.zeros((deg_dim, deg_dim))
                H_deg[i] [mapped_inds[0], mapped_inds[1]] = v
                H_deg[i] [mapped_inds[1], mapped_inds[0]] = v
            # now we need to transform from the basis of zero-order states to the basis of non-degenerate states
            # which we do by using our previously built PerturbationTheoryCorrections
            H_deg = corrs.operator_representation(H_deg, subspace=degenerate_space)
            # now that we've corrected those elements, we add on the diagonal terms,
            # add things up, and diagonalize
            H_to_diag = np.sum(H_deg, axis=0)
            H_to_diag[np.diag_indices_from(H_to_diag)] = corrs.energies
            deg_engs, deg_transf = np.linalg.eigh(H_to_diag)

            # finally we re-sort so that the new wavefunctions look maximally like the non-degenerate ones
            deg_transf = deg_transf.T # easier for the moment...
            state_set = set(np.arange(len(deg_transf)))
            sort_transf = deg_transf.copy()
            for i in range(len(deg_transf)):
                max_ov = np.max(deg_transf[:, i]**2)
                ov_thresh = .5
                # if max_ov < ov_thresh: # there must be a single mode that has more than 50% of the initial state character?
                #     raise PerturbationTheoryException("mode {} is has no contribution of greater than {}".format(
                #         i, ov_thresh
                #     ))
            sorting = [ -1 ] * len(deg_transf)
            for i in range(len(deg_transf)):
                o = np.argmax(abs(sort_transf[:, i]))
                sorting[i] = o
                sort_transf[o] = np.zeros(len(sort_transf))
            # sorting = [ np.argmax(abs(sort_transf[:, i])) for i in range(len(deg_transf)) ]
            # if len(sorting) != len(np.unique(sorting)):
            #     raise PerturbationTheoryException("After diagonalizing can't distinguish modes...")
            deg_engs = deg_engs[sorting,]
            deg_transf = deg_transf[sorting, :]
            # import McUtils.Plots as plt
            # plt.ArrayPlot(deg_transf).show()
            # raise Exception(deg_transf)
            # # raise Exception(engs, H)
            # wfns = np.dot(deg_transf, wfns)
            corrs.degenerate_energies = deg_engs
            corrs.degenerate_transf = deg_transf

        return corrs

    def get_representations(self, coupled_states):
        """
        Returns the representations of the perturbations in the basis of coupled states

        :param coupled_states:
        :type coupled_states:
        :return:
        :rtype:
        """
        total_states = int(np.prod(self.basis.quanta))
        if coupled_states is None:
            coupled_states = np.arange(total_states)  # expensive but...I guess what are you gonna do?
        elif not isinstance(coupled_states[0], int):
            coupled_states = self.basis.ravel_state_inds(coupled_states)

        m = coupled_states
        N = len(m)
        h_reps = self.perturbations
        H = [None] * len(h_reps)
        for i, h in enumerate(h_reps):
            sub = h[np.ix_(m, m)]
            if isinstance(sub, (int, np.integer, np.floating, float)):
                sub = np.full((N, N), sub)
            H[i] = sub  # type: np.ndarray
        return H

    def get_wavefunctions(self,
                          states,
                          coupled_states=None,
                          degeneracies=None,
                          order=2
                          ):
        """
        Gets a set of `PerturbationTheoryWavefunctions` from the perturbations defined by the Hamiltonian

        :param states: the states to get the index for, given either as indices or as a numbers of quanta
        :type states: Iterable[int] | Iterable[Iterable[int]]
        :param coupled_states: the list of states to explicitly allow to couple in
        :type coupled_states: Iterable[int] | Iterable[Iterable[int]]
        :param degeneracies: the pairs of states to be treated via degenerate perturbation theory
        :type degeneracies: (Iterable[int], Iterable[int])  | (Iterable[Iterable[int]], Iterable[Iterable[int]])
        :return: generated wave functions
        :rtype: PerturbationTheoryWavefunctions
        """

        from .Wavefunctions import PerturbationTheoryWavefunction, PerturbationTheoryWavefunctions

        if not isinstance(states[0], (int, np.integer)):
            states = self.basis.ravel_state_inds(states)
        total_states = int(np.prod(self.basis.quanta))
        if coupled_states is None:
            coupled_states = np.arange(total_states) # expensive but...I guess what are you gonna do?
        elif not isinstance(coupled_states[0], int):
            coupled_states = self.basis.ravel_state_inds(coupled_states)
        if (
                degeneracies is not None
                and not isinstance(degeneracies, (int, float, np.integer, np.floating))
                and not isinstance(degeneracies[0], (int, np.integer))
                and not isinstance(degeneracies[0][0], (int, np.integer))
        ):
            degeneracies = (
                self.basis.ravel_state_inds(degeneracies[0]),
                self.basis.ravel_state_inds(degeneracies[1])
            )

        corrs = self._get_state_VPT_corrections(
            self.perturbations,
            states,
            coupled_states,
            total_states,
            order,
            degenerate_states=degeneracies
            )

        return PerturbationTheoryWavefunctions(self.molecule, self.basis, corrs)

    def _martin_test(self, h_reps, states, threshold):
        """
        Applies the Martin Test to all of the specified states and returns the resulting correlation matrix

        :param states:
        :type states:
        :return:
        :rtype:
        """

        energies = np.diag(h_reps[0])
        state_energies = energies[states]
        diffs = np.array([e - energies for e in state_energies])
        for n, s in enumerate(states):
            diffs[n, s] = 1
        H1_blocks = h_reps[1][states,]
        corr_mat = (H1_blocks**4)/(diffs**3)

        deg_states = []
        for s, block in zip(states, corr_mat):
            big = np.where(np.abs(block) > threshold)[0]
            if len(big) > 0:
                deg_states.extend((s, d) for d in big)

        if len(deg_states) == 0:
            return None
        else:
            return np.array(deg_states).T