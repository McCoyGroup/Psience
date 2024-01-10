
import numpy as np, scipy as sp

from McUtils.Zachary import DensePolynomial
from McUtils.Scaffolding import Logger

from . Components import *
from .Wavefunctions import DGBWavefunctions
from .Solvers import DGBEigensolver

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
                  quadrature_degree=4,
                  expansion_degree=None,
                  expansion_type='multicenter',
                  reference_structure=None,
                  poly_coeffs=None,
                  pairwise_potential_functions=None
                  ):
        logger = Logger.lookup(logger)

        coords, potential_function = cls.construct_gaussians(
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
            coords = coords.optimize(optimize_centers)

        potential = cls.construct_potential(
            potential_function,
            coords,
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
                            coords,
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
            pairwise_functions=(
                None
                    if pairwise_potential_functions is None else
                coords.coords.pairwise_potential_evaluator(pairwise_potential_functions)
            ),
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

    default_min_singular_value=1e-4
    def diagonalize(self,
                    *,
                    print_debug_info=False,
                    subspace_size=None,
                    min_singular_value=None,
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
            eigs, evecs = DGBEigensolver.classic_eigensolver(
                H, self.S, self,
                min_singular_value=min_singular_value,
                subspace_size=subspace_size,
                nodeless_ground_state=nodeless_ground_state
            )

        elif mode == 'fix-heiberger': # Implementation of the Fix-Heiberger algorithm
            eigs, evecs = DGBEigensolver.fix_heiberger(H, self.S, self, eps=eps)

        elif mode == 'similarity':
            eigs, evecs = DGBEigensolver.similarity_mapped_solver(H, self.S, self)

        elif callable(mode):
            eigs, evecs = mode(H, self.S, )
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