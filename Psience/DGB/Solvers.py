
import numpy as np, scipy as sp

from .Wavefunctions import DGBWavefunctions

class DGBEigensolver:

    @classmethod
    def get_orthogonal_transform(self, S, min_singular_value=None, subspace_size=None):

        sig, evecs = np.linalg.eigh(S)

        # we'll ignore anything too large?
        # and anything negative?
        # mostly to remove the troubles from rotations and translations...
        base_loc = np.where(np.logical_and(sig > 0, np.abs(sig) < 1e30))
        sig = sig[base_loc]

        # sort = np.argsort(np.abs(sig))
        # sig = sig[sort]
        # evecs = evecs[:, sort]

        if subspace_size is not None:
            good_loc = slice(max(0, len(sig) - subspace_size), len(sig))
        elif min_singular_value is not None:
            good_loc = np.where(sig > min_singular_value)[0]
        else:
            good_loc = np.arange(len(S))
        d = np.diag(1 / np.sqrt(sig[good_loc]))

        gl = base_loc[0][good_loc]
        L = evecs[:, gl]
        Q = L @ d @ L.T

        # sorting = np.concatenate([good_loc, bad_loc[0]])
        # Lsort = evecs[:, sorting]
        Qinv = L @ np.diag(np.sqrt(sig[good_loc])) @ L.T

        if L.shape[0] == L.shape[1]: # no contraction
            return Q, Qinv, None

        Qe, QL = np.linalg.eigh(Q)
        qsub = np.where(Qe > 1e-8)[0]
        Qq = QL[:, qsub]
        qrest = np.where(Qe <= 1e-8)[0]
        qsort = np.concatenate([qsub, qrest])
        Qqinv = QL[:, qsort].T  # transforms back to the Q basis

        return Q, Qinv, (Qq, Qqinv)

    @classmethod
    def classic_eigensolver(self, H, S,
                            hamiltonian,
                            min_singular_value=None,
                            subspace_size=None,
                            nodeless_ground_state=True
                            ):
        Q, Qinv, proj = self.get_orthogonal_transform(
            S,
            min_singular_value=min_singular_value,
            subspace_size=subspace_size
        )

        if proj is None:
            # print(Q.shape, H.shape, self.S.shape, self.centers.shape)
            eigs, evecs = sp.linalg.eigh(H, S)
            Qq = np.eye(len(Q))
        else:
            Qq, Qqinv = proj
            hamiltonian.logger.log_print('solving with subspace size {}'.format(Qq.shape[1]))
            H = Qq.T @ Q @ H @ Q.T @ Qq  # in our projected orthonormal basis
            eigs, evecs = np.linalg.eigh(H)
            evecs = np.concatenate(
                [
                    evecs,
                    np.zeros((Qqinv.shape[1] - Qq.shape[1], len(evecs)))
                ],
                axis=0
            )
            evecs = Q.T @ Qqinv.T @ evecs
        if nodeless_ground_state:
            # fast enough if we have not that many points...
            gswfn = DGBWavefunctions(eigs, evecs, hamiltonian=hamiltonian)[0]
            gs = gswfn.evaluate(hamiltonian.gaussians.coords.centers)
            abs_gs = np.abs(gs)
            signs = np.sign(gs[abs_gs > np.max(abs_gs) * 1e-1])
            diffs = np.abs(np.diff(signs))
            # print(gs[np.argsort(np.abs(gs))][-5:], Qq.shape[1] - 1)
            if np.sum(diffs) > 0:  # had a sign flip
                eigs, evecs = self.classic_eigensolver(
                    H, S, hamiltonian,
                    subspace_size=Qq.shape[1] - 1,
                    nodeless_ground_state=Qq.shape[1] > 1  # gotta bottom out some time...
                )
        return eigs, evecs

    @classmethod
    def get_eigensimilarity_subspace_size(cls, H, S):
        eigs, Qs = np.linalg.eigh(S)
        eigh, Qh = np.linalg.eigh(H)

        similarity_matrix = Qs.T @ Qh

        dets = np.array([
            abs(np.linalg.det(similarity_matrix[-n:, -n:]))
            for n in range(1, len(eigs) + 1)
        ])

        # import McUtils.Plots as plt
        # smat = plt.MatrixPlot(similarity_matrix)
        # sdets = plt.Plot(list(range(1, len(eigs) + 1)), dets)
        # splot_dets = plt.Plot(range(1, len(eigs)), np.diff(dets)).show()
        # raise Exception(...)

        det_chunks = np.split(np.arange(len(dets)), np.where(np.abs(np.diff(dets)) > .01)[0] + 1)
        good_runs = [c for c in det_chunks if len(c) > 5 and np.mean(dets[c]) > .95]

        if len(good_runs) == 0 or len(good_runs[-1]) == 0:
            raise ValueError("couldn't find stable eigenspace")

        return good_runs[-1][-1]

        # raise Exception(good_runs)
        # raise Exception([dets[c] for c in det_chunks])


    @classmethod
    def similarity_mapped_solver(self, H, S, hamiltonian):
        subspace_size = self.get_eigensimilarity_subspace_size(H, S)

        return self.classic_eigensolver(H, S, hamiltonian=hamiltonian, subspace_size=subspace_size)

    @classmethod
    def fix_heiberger(self, H, S, hamiltonian, eps=1e-5):
        # raise NotImplementedError("I think my Fix-Heiberger is broken?"
        #                           "And also potentially not even the right algorithm for this problem"
        #                           "Like I potentially _want_ a large epsilon since I want to get the best"
        #                           "possible energies, but Fix-Heiberger will call that unstable I think")

        d, Q = np.linalg.eigh(S)
        d = np.flip(d)
        Q = np.flip(Q, axis=1)

        N = len(d)
        cutoff = np.max(d) * eps
        g1 = np.where(d >= cutoff)
        if len(g1) == 0:
            raise ValueError("totally zero S matrix")
        n1 = len(g1[0])
        n2 = N - n1

        if n2 == 0:
            print("Falling back on classic method...")
            return self.classic_eigensolver(H, S, hamiltonian)

        # D = np.diag(d[g1])
        # F = np.diag(d[g2])
        # B0 = np.diag(np.concatenate([
        #     d[g1],
        #     np.zeros(n2)
        # ]))
        R = np.diag(np.concatenate([
            1 / np.sqrt(d[g1]),
            np.ones(n2)
        ]))
        # with np.printoptions(linewidth=1e8):
        #     raise Exception(R@B0@R)
        A1 = R.T @ Q.T @ H @ Q @ R

        A22 = A1[n1:, n1:]
        d2, Q22 = np.linalg.eigh(A22)
        d2 = np.flip(d2)
        Q22 = np.flip(Q22, axis=1)

        cut2 = np.max(d2) * eps
        g3 = np.where(d2 >= cut2)[0]
        n3 = len(g3)
        if n3 == 0:
            # print("Early exiting after first iteration")
            if n1 <= n2:
                raise ValueError("singular problem (case 2)")
            A12 = A1[:n1, n1:]
            Q12, R13 = np.linalg.qr(A12, mode='complete')
            if np.linalg.matrix_rank(A12) < n2:  # singular
                raise ValueError("singular problem (case 3)")

            # print("N3 = 0???")

            Q2 = np.eye(N)
            Q2[:n1, :n1] = Q12
            A2 = Q2.T @ A1 @ Q2.T

            A12 = A2[:n2, n2:n1]
            A13 = A2[n1:, :n2]
            A22 = A2[n2:n1, n2:n1]

            eigs, U2 = np.linalg.eigh(A22)
            U1 = np.zeros((n2, n1 - n2))
            U3 = -np.linalg.inv(A13) @ A12 @ U2

            U = np.concatenate(
                [
                    U1, U2, U3
                ],
                axis=0
            )

            evecs = Q @ R @ Q2 @ U

        else:

            g4 = np.where(d2 < cut2)[0]
            n4 = len(g4)
            Q2 = np.eye(N)
            Q2[n1:, n1:] = Q22
            A2 = Q2.T @ A1 @ Q2
            # B2 = Q2.T@B1@Q2

            if n4 == 0:
                # print("Early, well conditioned after first iteration")
                # print("N4 = 0??? {d2}".format(d2=d2))
                A11 = A2[:n1, :n1]
                A12 = A2[:n1, n1:]
                Dinv = np.diag(1 / d2)

                eigs, U1 = np.linalg.eigh(A11 - A12 @ Dinv @ A12.T)
                U2 = -Dinv @ A12.T @ U1

                U = np.concatenate(
                    [
                        U1, U2
                    ],
                    axis=0
                )

                evecs = Q @ R @ Q2 @ U
            else:  # second iteration of this partitioning...
                if n1 <= n4:
                    raise ValueError("singular problem (case 5)")

                A2[-n4:, -n4:] = np.zeros((n4, n4))  #
                # B2 = B1

                A13 = A2[:n1, n1 + n3:]
                Q33, R14 = np.linalg.qr(A13, mode='complete')
                if np.linalg.matrix_rank(A13) < n4:  # singular
                    raise ValueError("singular problem (case 6)")
                # A14 = R14[:n4]

                Q3 = np.eye(N)
                Q3[:n1, :n1] = Q33
                A3 = Q3.T @ A2 @ Q3

                # B3 = B1

                n5 = n1 - n4
                A11 = A3[:n4, :n4]
                A12 = A3[:n4, n4:n1]
                A13 = A3[:n4, n1:n1 + n3]
                A14 = A3[:n4, -n4:]

                A22 = A3[n4:n1, n4:n1]
                A23 = A3[n4:n1, n1:n1 + n3]
                # A24 = A3[n4:n5, -n4:]
                #
                # raise Exception(A24)

                Dinv = np.diag(1 / d2[g3])

                U1 = np.zeros((n4, n5))
                eigs, U2 = np.linalg.eigh(
                    A22 - A23 @ Dinv @ A23.T
                )
                U3 = -Dinv @ A23.T @ U2
                U4 = -np.linalg.inv(A14) @ (A12 @ U2 + A13 @ U3)

                U = np.concatenate(
                    [
                        U1, U2, U3, U4
                    ],
                    axis=0
                )

                evecs = Q @ R @ Q2 @ Q3 @ U

        return eigs, evecs