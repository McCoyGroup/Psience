
import numpy as np, scipy as sp, collections

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
            good_loc = np.arange(len(sig))
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
                            nodeless_ground_state=False
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
            Hq = Qq.T @ Q @ H @ Q.T @ Qq  # in our projected orthonormal basis
            eigs, evecs = np.linalg.eigh(Hq)
            evecs = np.concatenate(
                [
                    evecs,
                    np.zeros((Qqinv.shape[1] - Qq.shape[1], len(evecs)))
                ],
                axis=0
            )
            evecs = Q.T @ Qqinv.T @ evecs
        if nodeless_ground_state:
            from .Wavefunctions import DGBWavefunctions
            # fast enough if we have not that many points...
            gswfn = DGBWavefunctions(eigs, evecs, hamiltonian=hamiltonian)[0]
            gs = gswfn.evaluate(hamiltonian.gaussians.coords.centers)
            abs_gs = np.abs(gs)
            signs = np.sign(gs[abs_gs > np.max(abs_gs) * 1e-1])
            diffs = np.abs(np.diff(signs))
            # print(gs[np.argsort(np.abs(gs))][-5:], Qq.shape[1] - 1)
            if np.sum(diffs) > 0:  # had a sign flip
                if subspace_size is not None:
                    subspace_size = min(subspace_size, Qq.shape[1])
                else:
                    subspace_size = Qq.shape[1]
                eigs, evecs = self.classic_eigensolver(
                    H, S, hamiltonian,
                    subspace_size=subspace_size - 1,
                    nodeless_ground_state=subspace_size > 1  # gotta bottom out some time...
                )
        return eigs, evecs

    eigensimilarity_cutoff = None
    eigensimilarity_chunk_size = 3
    similar_determinant_cutoff = 0.05
    @classmethod
    def get_eigensimilarity_subspace_size(cls, H, S,
                                          similarity_cutoff=None,
                                          similarity_chunk_size=None,
                                          similar_det_cutoff=None
                                          ):

        if similarity_cutoff is None:
            similarity_cutoff = cls.eigensimilarity_cutoff
        if similarity_chunk_size is None:
            similarity_chunk_size = cls.eigensimilarity_chunk_size

        eigs, Qs = np.linalg.eigh(S)
        eigh, Qh = np.linalg.eigh(H)

        similarity_matrix = Qs.T @ Qh

        # compute moving average of abs of blocks of dets
        prev_dets = collections.deque(maxlen=similarity_chunk_size)
        cur_sum = 0
        strike_det_cutoff = 0.05
        strikes = 0
        max_strikes = 10
        blocks = np.zeros(len(eigs))
        for n in range(1, len(eigs) + 1):
            subdet = abs(np.linalg.det(similarity_matrix[-n:, -n:]))
            if len(prev_dets) == similarity_chunk_size:
                oldest = prev_dets[0]
                cur_sum -= oldest
            prev_dets.append(subdet)
            cur_sum += subdet
            avg = cur_sum / len(prev_dets)
            blocks[n-1] = avg
            if avg < strike_det_cutoff:
                strikes += 1
                if strikes > max_strikes:
                    break
            else:
                strikes = 0

        # dets = np.array([
        #     abs(np.linalg.det(similarity_matrix[-n:, -n:]))
        #     for n in range(1, len(eigs) + 1)
        # ])
        #
        # w = similarity_chunk_size
        # avgs = np.cumsum(dets[:-1])
        # avgs[w:] = avgs[w:] - avgs[:-w]
        # blocks = avgs[w - 1:] / w

        if similarity_cutoff is None:
            similarity_cutoff = np.floor(np.max(blocks)*100/5)*5 / 100 # next lowest .5
        block_pos = np.where(blocks > similarity_cutoff)
        if len(block_pos) == 0 or len(block_pos[0]) == 0:
            raise ValueError(
                "couldn't find subspace size where average of previous {} dets greater than {}".format(
                    similarity_chunk_size,
                    similarity_cutoff
                )
            )

        return int(block_pos[0][-1] + 1) #+ np.ceil(w/2))


        # det_chunks = np.split(np.arange(len(dets)), np.where(np.abs(np.diff(dets)) > similar_det_cutoff)[0] + 1)
        # good_runs = [
        #     c for c in det_chunks
        #         if len(c) > similarity_chunk_size
        #             and
        #         np.mean(dets[c]) > similarity_cutoff
        # ]
        #
        # if len(good_runs) == 0 or len(good_runs[-1]) == 0:
        # raise ValueError("couldn't find stable eigenspace")

        # return good_runs[-1][-1]

        # raise Exception(good_runs)
        # raise Exception([dets[c] for c in det_chunks])


    @classmethod
    def similarity_mapped_solver(self, H, S,
                                 hamiltonian,
                                 similarity_cutoff=None,
                                 similarity_chunk_size=None,
                                 similar_det_cutoff=None
                                 ):


        subspace_size = self.get_eigensimilarity_subspace_size(H, S,
                                               similarity_cutoff=similarity_cutoff,
                                               similarity_chunk_size=similarity_chunk_size,
                                               similar_det_cutoff=similar_det_cutoff
                                               )

        hamiltonian.logger.log_print(
            'diagonalizing in the space of {} S functions'.format(np.sum(subspace_size))
        )

        eigh, Qh = np.linalg.eigh(H)

        eigs, Qs = np.linalg.eigh(S)
        Qs = Qs[:, -subspace_size:]
        eigs = eigs[-subspace_size:]

        A = Qs.T @ Qh

        if np.any(eigs < -1e-14):
            raise ValueError("negative eigenvalues from the S-matrix encountered")
        D = np.diag(1 / np.sqrt(np.abs(eigs)))
        C = D @ A @ np.diag(eigh) @ A.T @ D

        engs, Q = np.linalg.eigh(C)
        # if shift is not None:
        #     engs = engs - shift
        #
        # nz_pos = np.abs(engs) > 1e-12
        # engs = engs[nz_pos]
        # Q = Q[:, nz_pos]

        return engs, Qs @ D @ Q


    default_shift = 2
    @classmethod
    def shift_similarity_solver(self, H, S,
                                 hamiltonian,
                                 similarity_cutoff=None,
                                 similarity_chunk_size=None,
                                 similar_det_cutoff=None,
                                 similarity_shift=None
                                 ):
        if similarity_shift is None:
            similarity_shift = self.default_shift

        eigs, Qs = np.linalg.eigh(S)
        match_pos = np.full(len(S), True)
        subspace_size = len(S)
        shift = 0
        for shift in range(similarity_shift+1):
            shift_H = H + S * shift
            subspace_size = min(
                subspace_size,
                self.get_eigensimilarity_subspace_size(
                    shift_H, S,
                    similarity_cutoff=similarity_cutoff,
                    similarity_chunk_size=similarity_chunk_size,
                    similar_det_cutoff=similar_det_cutoff
                )
            )
            eigh, Qh = np.linalg.eigh(shift_H)
            A = Qs.T @ Qh
            match_pos = np.logical_and(match_pos, np.abs(np.diag(A)) > .9)

        eigh, Qh = np.linalg.eigh(H + S * similarity_shift)
        match_pos[:-subspace_size] = False

        hamiltonian.logger.log_print(
            'diagonalizing in the space of {} common eigenfunctions'.format(np.sum(match_pos))
        )

        eigh = eigh[match_pos]
        Qh = Qh[:, match_pos]
        Qs = Qs[:, match_pos]
        eigs = eigs[match_pos]

        A = Qs.T @ Qh
        D = np.diag(1 / np.sqrt(eigs))
        C = D @ A @ np.diag(eigh) @ A.T @ D

        engs, Q = np.linalg.eigh(C)
        engs = engs - similarity_shift

        nz_pos = np.abs(engs) > 1e-12
        engs = engs[nz_pos]
        Q = Q[:, nz_pos]

        return engs, Qs @ D @ Q
    @classmethod
    def low_rank_solver(cls, H, S, hamiltonian,
                        low_rank_energy_cutoff=None,
                        low_rank_overlap_cutoff=None,
                        low_rank_shift=None
                        ):
        from McUtils.Data import UnitsData

        if low_rank_energy_cutoff is None:
            low_rank_energy_cutoff = 1 / UnitsData.hartrees_to_wavenumbers #np.finfo(H.dtype).eps
        if low_rank_overlap_cutoff is None:
            low_rank_overlap_cutoff = np.finfo(S.dtype).eps
        if low_rank_shift is None:
            low_rank_shift = .5

        eigs, Qs = np.linalg.eigh(S)
        eigh, Qh = np.linalg.eigh(H + S*low_rank_shift)
        # import McUtils.Plots as plt
        # plt.ArrayPlot(Qs.T @ Qh).show()
        pos_H = np.abs(eigh) > low_rank_energy_cutoff
        pos_S = eigs > low_rank_overlap_cutoff
        while np.min(eigh[pos_H]) < 0 and low_rank_shift < 1e6:
            low_rank_shift = 2 * (low_rank_shift - np.min(eigh[pos_H]))
            eigh, Qh = np.linalg.eigh(H + S * low_rank_shift)
            pos_H = np.abs(eigh) > low_rank_energy_cutoff

        hamiltonian.logger.log_print(
            'solving with {} H values and {} overlap functions'.format(
                np.sum(pos_H),
                np.sum(pos_S)
            )
        )

        Qh = Qh[:, pos_H]
        eigh = eigh[pos_H]
        Qs = Qs[:, pos_S]
        eigs = eigs[pos_S]
        if np.min(eigh) >= 0:
            A = Qs.T @ Qh
            B = np.diag(1/np.sqrt(eigs)) @ A @ np.diag(np.sqrt(eigh))
            Q, sqrt_engs, _ = np.linalg.svd(B)
            eng_sort = np.argsort(sqrt_engs)
            engs = sqrt_engs[eng_sort]**2 - low_rank_shift
            Q = Q[eng_sort].T
            Q = Qs @ np.diag(1 / np.sqrt(eigs)) @ Q
        else:
            A = Qs.T @ Qh
            C = np.diag(1/np.sqrt(eigs)) @ A @ np.diag(eigh) @ A.T @ np.diag(1/np.sqrt(eigs))
            engs, Q = np.linalg.eigh(C)

            raise Exception(engs)
            # eng_sort = np.argsort(sqrt_engs)
            # engs = sqrt_engs[eng_sort]**2 - low_rank_shift


        return engs, Q

    @classmethod
    def cholesky_solver(cls, H, S, hamiltonian):

        L = np.linalg.cholesky(S)
        D = np.diag(np.diag(L))
        # Dinv = np.diag(1/np.diag(L))
        # L = L @ Dinv

        L_inv = np.linalg.inv(L)
        new_H = L_inv @ H @ L_inv.T

        engs = np.linalg.eigvalsh(new_H)

        return engs, Q

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