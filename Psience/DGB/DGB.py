
import numpy as np, scipy as sp

from McUtils.Zachary import DensePolynomial
from McUtils.Scaffolding import Logger

from . Components import *
from .Wavefunctions import DGBWavefunctions

__all__ = [
    "DGB"
]

__reload_hook__ = ['.Wavefunctions', '..Molecools', '..MixtureModes']


class DGB:
    """

    """

    @classmethod
    def run(cls,
            centers,
            potential_function,
            masses=None,
            mass_weighted=False,
            atoms=None,
            projection_indices=None,
            modes=None,
            transformations=None,
            alphas=None,
            logger=True,
            clustering_radius=None,
            min_singular_value=None,
            num_svd_vectors=None,
            svd_contrib_cutoff=None,
            optimize_centers=None,
            quadrature_degree=None,
            expansion_degree=None,
            expansion_type=None,
            reference_structure=None,
            pairwise_potential_functions=None # we can integrate this cleanly, not sure about anything else
            ):
        opts = dict(
            masses=masses,
            mass_weighted=mass_weighted,
            atoms=atoms,
            alphas=alphas,
            projection_indices=projection_indices,
            modes=modes,
            transformations=transformations,
            logger=logger,
            clustering_radius=clustering_radius,
            min_singular_value=min_singular_value,
            num_svd_vectors=num_svd_vectors,
            svd_contrib_cutoff=svd_contrib_cutoff,
            optimize_centers=optimize_centers,
            quadrature_degree=quadrature_degree,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            reference_structure=reference_structure,
            pairwise_potential_functions=pairwise_potential_functions
        )
        opts = {k:v for k,v in opts.items() if v is not None}
        logger = Logger.lookup(logger)
        with logger.block(tag="Running distributed Gaussian basis calculation"):
            opts['logger'] = logger
            ham = cls.construct(centers, potential_function, **opts)
            return ham.get_wavefunctions()

    @classmethod
    def construct(cls,
                  centers,
                  potential_function,
                  gmat_function=None,
                  masses=None,
                  # mass_weighted=False,
                  atoms=None,
                  alphas=None,
                  # coord_shape=None,
                  projection_indices=None,
                  transformations=None,
                  internals=None,
                  modes=None,
                  logger=False,
                  optimize_centers=False,
                  clustering_radius=-1,
                  min_singular_value=1e-4,
                  num_svd_vectors=None,
                  svd_contrib_cutoff=1e-3,
                  quadrature_degree=4,
                  expansion_degree=None,
                  expansion_type='multicenter',
                  reference_structure=None,
                  poly_coeffs=None,
                  pairwise_potential_functions=None
                  ):
        logger = Logger.lookup(logger)

        coords = cls.construct_gaussians(
            centers,
            alphas,
            potential_function,
            gmat_function=gmat_function,
            masses=masses,
            atoms=atoms,
            internals=internals,
            modes=modes,
            transformations=transformations,
            # projection_indices=projection_indices,
            poly_coeffs=poly_coeffs,
            logger=logger
        )
        if optimize_centers:
            coords = coords.optimize(
                clustering_radius=clustering_radius,
                min_singular_value=min_singular_value,
                num_svd_vectors=num_svd_vectors,
                svd_contrib_cutoff=svd_contrib_cutoff
            )

        potential = cls.construct_potential(
            potential_function,
            quadrature_degree=quadrature_degree,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            # reference_structure=reference_structure,
            pairwise_potential_functions=pairwise_potential_functions,
            logger=logger
        )

        return cls(
            coords,
            potential,
            logger=logger
        )

    @classmethod
    def construct_gaussians(cls,
                            centers,
                            alphas,
                            potential_function,
                            gmat_function=None,
                            masses=None,
                            atoms=None,
                            internals=None,
                            modes=None,
                            transformations=None,
                            # projection_indices=projection_indices,
                            poly_coeffs=None,
                            logger=None
                            ):
        # here to be overridden
        return DGBGaussians.construct(
            centers,
            alphas,
            potential_function,
            gmat_function=gmat_function,
            masses=masses,
            atoms=atoms,
            internals=internals,
            modes=modes,
            transformations=transformations,
            # projection_indices=projection_indices,
            poly_coeffs=poly_coeffs,
            logger=logger
        )

        # self._S, self._T, self._V = None, None, None
        # self._scaling_S = None
        # self.logger = Logger.lookup(logger)

        # self.potential_function = potential_function
        # self.quadrature_degree = quadrature_degree
        # self.expansion_degree = expansion_degree
        # self.expansion_type = expansion_type
        # self.ref = reference_structure
        #
        # self.pairwise_potential_functions = pairwise_potential_functions
        # self.modes = modes
        # self.internals = internals
    @classmethod
    def construct_potential(cls,
                            potential_function,
                            quadrature_degree=None,
                            expansion_degree=None,
                            expansion_type=None,
                            # reference_structure=None,
                            pairwise_potential_functions=None,
                            logger=None
                            ):
        return DGBPotentialEnergyEvaluator(
            potential_function,
            quadrature_degree=quadrature_degree,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            # reference_structure=reference_structure,
            pairwise_functions=pairwise_potential_functions,
            logger=logger
        )

    def __init__(self,
                 gaussians:DGBGaussians,
                 potential:DGBPotentialEnergyEvaluator,
                 logger=None
                 ):
        self.gaussians = gaussians
        self.pot = potential
        self._V = None
        self.logger = Logger.lookup(logger)

    @property
    def S(self):
        return self.gaussians.S
    @property
    def T(self):
        return self.S * self.gaussians.T
    @property
    def V(self):
        if self._V is None:
            self._V = self.S * self.pot.evaluate_pe(self.gaussians.overlap_data)
        return self._V

    def get_kinetic_polynomials(self):
        raise NotImplementedError("deprecated")
        # TODO: implement _multipolynomial_ support so that we can handle
        #       this in a simpler form

        ndim = self.centers.shape[-1]
        if self.transformations is None and self.inds is not None:
            ndim = len(self.inds)

        core_polys = self.get_base_polynomials()
        core_derivs = [ # these will be indexed by columns
            [poly.deriv(i) if isinstance(poly, DensePolynomial) else 0 for i in range(ndim)]
            for poly in core_polys
        ]
        core_d2s = [
            [p.deriv(i) if isinstance(p, DensePolynomial) else 0 for i, p in enumerate(poly_list)]
            for poly_list in core_derivs
        ]

        if self.transformations is None:
            centers = self.centers
            alphas = self.alphas
            if self.inds is not None:
                ndim = len(self.inds)
                centers = centers[:, self.inds]
                alphas = alphas[:, self.inds]
            exp_polys = [ # these will be indexed by columns as well
                DensePolynomial.from_tensors([0, np.zeros(ndim), np.diag(2*a)], shift=-pt, rescale=False)
                for pt,a in zip(centers, alphas)
            ]
        else:
            sigs = self.get_inverse_covariances()
            centers = self.centers
            if self.inds is not None:
                Lt, L = self.transformations
                Lt = Lt[:, :, self.inds]
                L = L[:, self.inds, :]
                centers = Lt@( L @ centers[:, :, np.newaxis] )
                centers = centers.reshape(self.centers.shape)

            raise Exception(sigs[2], centers[2])

            exp_polys = [
                DensePolynomial.from_tensors([0, np.zeros(ndim), sig], shift=-pt)
                for pt,sig in zip(
                    centers,
                    sigs
                )
            ]
        def _d(p, i):
            if isinstance(p, DensePolynomial):
                d = p.deriv(i)
                if isinstance(d, DensePolynomial):
                    d = d.clip()
            else:
                d = p
            return d
        exp_derivs = [
            [-1/2*_d(poly, i) for i in range(ndim)]
            for poly in exp_polys
        ]
 #        raise Exception(exp_derivs[0][1].coeffs)
 #        """
 #        [[[-0.24875972 -1.62256994]
 #  [-1.91803054 -0.        ]]
 #
 # [[ 3.58730135 -0.        ]
 #  [-0.         -0.        ]]]
 #  """
        exp_d2s = [
            [_d(p, i) for i,p in enumerate(poly_list)]
            for poly_list in exp_derivs
        ]
        # raise Exception(
        #     (
        #             (exp_derivs[0][0] * exp_derivs[0][0])
        #             + exp_d2s[0][0].coeffs*exp_d2s[0][0].scaling
        #     ).coeffs,
        #     core_polys[0].coeffs
        # )

        raise Exception(
            exp_polys[2].coeffs
        )

        # and now we can construct the polynomials for each center and axis
        def _k(p, pd1, pd2, ed1, ed2):
            t = pd2 + 2*pd1 * ed1 + p * ((ed1 * ed1) + ed2)
            if isinstance(t, DensePolynomial):
                t = t.clip(threshold=1e-8)
            return t
        kinetic_polys = [
            [
                _k(p, pd1, pd2, ed1, ed2)
                for pd1, pd2, ed1, ed2 in zip(pd1s, pd2s, ed1s, ed2s)
            ]
            for p, pd1s, pd2s, ed1s, ed2s in zip(
                core_polys,
                core_derivs,
                core_d2s,
                exp_derivs,
                exp_d2s
            )
        ]

        # raise Exception(
        #     kinetic_polys[0][0].coeffs,
        #     kinetic_polys[0][1],
        #     kinetic_polys[0][2].coefficient_tensors,
        #     # kinetic_polys[0][0].coeffs,
        # )

        return core_polys, kinetic_polys

    @staticmethod
    def mass_weighted_eval(function, coords, masses, deriv_order=None, **opts):
        raise NotImplementedError("deprecated")
        # expects just x and y coordinates for the atoms
        mass_vec = np.asanyarray(masses)

        coords = coords / np.sqrt(mass_vec[np.newaxis])
        coords = coords.reshape(-1, 3, 2)
        derivs = function(coords, deriv_order=deriv_order, **opts)
        new_d = []
        for n, do in enumerate(derivs):
            if n > 0:
                mv = np.sqrt(mass_vec)[np.newaxis]
                weight = mv
                for j in range(n - 1):
                    mv = np.expand_dims(mv, 0)
                    weight = np.expand_dims(weight, -1) * mv
                # with np.printoptions(linewidth=1e8):
                #     print(weight[0])
                new_d.append(do / weight)
            else:
                new_d.append(do)
        derivs = new_d
        return derivs

    def evaluate_multiplicative_operator(self,
                                         function,
                                         handler=None,
                                         expansion_degree=None,
                                         expansion_type=None,
                                         quadrature_degree=None,
                                         pairwise_functions=None # integrate out pairwise contribution
                                         ):
        raise NotImplementedError("phasing this out")

        pot_mat = self.evaluate_multiplicative_operator_base(
            function,
            handler=handler,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            quadrature_degree=quadrature_degree,
            pairwise_functions=pairwise_functions
        )

        # with np.printoptions(linewidth=1e8):
        #     print(pot_mat)

        if self._scaling_S is None:
            _ = self.S
        S = self._scaling_S
        for _ in range(pot_mat.ndim - 2):
            S = np.expand_dims(S, -1)
        return S * pot_mat

    def get_V(self, potential_handler=None, expansion_degree=None, expansion_type=None, quadrature_degree=None):
        raise NotImplementedError("deprecated")
        self.logger.log_print("calculating potential matrix")
        pot = self.evaluate_multiplicative_operator(
            self.potential_function,
            handler=potential_handler,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            quadrature_degree=quadrature_degree
        )

        return pot

    def get_orthogonal_transform(self, min_singular_value=None, subspace_size=None):
        if min_singular_value is None:
            min_singular_value = self.min_singular_value

        S = self.S
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

    default_min_singular_value=1e-4
    def diagonalize(self, print_debug_info=False, subspace_size=None, min_singular_value=None,
                    eps=5e-4,
                    mode='classic',
                    nodeless_ground_state=True
                    ):

        if min_singular_value is None:
            min_singular_value = self.default_min_singular_value
        if min_singular_value is not None:
            self.logger.log_print("solving with min_singular_value={ms}", ms=min_singular_value)

        H = self.T + self.V
        if print_debug_info:
            print('Condition Number:', np.linalg.cond(self.S))
            with np.printoptions(linewidth=1e8, threshold=1e8):
                print("Potential Matrix:")
                print(self.V)

        if min_singular_value is None:
            min_singular_value = 1e-4
        # mode = 'fix-heiberger'
        if mode == "classic":
            Q, Qinv, proj = self.get_orthogonal_transform(min_singular_value=min_singular_value, subspace_size=subspace_size)
            if proj is None:
                # print(Q.shape, H.shape, self.S.shape, self.centers.shape)
                eigs, evecs = sp.linalg.eigh(H, self.S)
                Qq = np.eye(len(Q))
            else:
                Qq, Qqinv = proj
                H = Qq.T @ Q @ H @ Q.T @ Qq # in our projected orthonormal basis
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
                gswfn = DGBWavefunctions(eigs, evecs, hamiltonian=self)[0]
                gs = gswfn.evaluate(self.gaussians.coords.centers)
                abs_gs = np.abs(gs)
                signs = np.sign(gs[abs_gs > np.max(abs_gs) * 1e-1])
                diffs = np.abs(np.diff(signs))
                # print(gs[np.argsort(np.abs(gs))][-5:], Qq.shape[1] - 1)
                if np.sum(diffs) > 0: # had a sign flip
                    eigs, evecs = self.diagonalize(
                        mode='classic',
                        subspace_size=Qq.shape[1] - 1,
                        nodeless_ground_state=Qq.shape[1] > 1 # gotta bottom out some time...
                    )

        elif mode == 'stable': # Implementation of the Fix-Heiberger algorithm
            # raise NotImplementedError("I think my Fix-Heiberger is broken?"
            #                           "And also potentially not even the right algorithm for this problem"
            #                           "Like I potentially _want_ a large epsilon since I want to get the best"
            #                           "possible energies, but Fix-Heiberger will call that unstable I think")

            S = self.S
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
                return self.diagonalize(mode='classic')

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
            A1 = R.T@Q.T@H@Q@R

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
                U3 = -np.linalg.inv(A13)@A12@U2

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
                A2 = Q2.T@A1@Q2
                # B2 = Q2.T@B1@Q2

                if n4 == 0:
                    # print("Early, well conditioned after first iteration")
                    # print("N4 = 0??? {d2}".format(d2=d2))
                    A11 = A2[:n1, :n1]
                    A12 = A2[:n1, n1:]
                    Dinv = np.diag(1 / d2)

                    eigs, U1 = np.linalg.eigh(A11 - A12@Dinv@A12.T)
                    U2 = -Dinv@A12.T@U1

                    U = np.concatenate(
                        [
                            U1, U2
                        ],
                        axis=0
                    )

                    evecs = Q @ R @ Q2 @ U
                else: # second iteration of this partitioning...
                    if n1 <= n4:
                        raise ValueError("singular problem (case 5)")

                    A2[-n4:, -n4:] = np.zeros((n4, n4)) #
                    # B2 = B1

                    A13 = A2[:n1, n1+n3:]
                    Q33, R14 = np.linalg.qr(A13, mode='complete')
                    if np.linalg.matrix_rank(A13) < n4: # singular
                        raise ValueError("singular problem (case 6)")
                    # A14 = R14[:n4]

                    Q3 = np.eye(N)
                    Q3[:n1, :n1] = Q33
                    A3 = Q3.T@A2@Q3

                    # B3 = B1

                    n5 = n1 - n4
                    A11 = A3[:n4, :n4]
                    A12 = A3[:n4, n4:n1]
                    A13 = A3[:n4, n1:n1+n3]
                    A14 = A3[:n4, -n4:]

                    A22 = A3[n4:n1, n4:n1]
                    A23 = A3[n4:n1, n1:n1+n3]
                    # A24 = A3[n4:n5, -n4:]
                    #
                    # raise Exception(A24)

                    Dinv = np.diag(1 / d2[g3])

                    U1 = np.zeros((n4, n5))
                    eigs, U2 = np.linalg.eigh(
                        A22 - A23@Dinv@A23.T
                    )
                    U3 = -Dinv@A23.T@U2
                    U4 = -np.linalg.inv(A14)@(A12@U2 + A13@U3)

                    U = np.concatenate(
                        [
                            U1, U2, U3, U4
                        ],
                        axis=0
                    )

                    evecs = Q @ R @ Q2 @ Q3 @ U

        else:
            raise ValueError("unknown solver {}".format(mode))

                # raise Exception((eigs[1:] - eigs[0])*219475.6)


        return eigs, evecs

    def get_wavefunctions(self, print_debug_info=False, min_singular_value=None, subspace_size=None, nodeless_ground_state=None,
                          mode=None, stable_epsilon=None
                          ):
        # print("======="*25)
        ops = dict(
            print_debug_info=print_debug_info,
            min_singular_value=min_singular_value,
            nodeless_ground_state=nodeless_ground_state,
            subspace_size=subspace_size,
            mode=mode,
            eps=stable_epsilon
        )
        ops = {k:v for k,v in ops.items() if v is not None}
        eigs, evecs = self.diagonalize(**ops)
        return DGBWavefunctions(eigs, evecs, hamiltonian=self)