
import numpy as np, scipy as sp

from McUtils.Data import UnitsData
from McUtils.Zachary import DensePolynomial
from McUtils.Scaffolding import Logger, NullLogger

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
    def run_simple(cls,
            centers,
            potential_function,
            *,
            masses=None,
            # mass_weighted=False,
            atoms=None,
            # projection_indices=None,
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
            # reference_structure=None,
            pairwise_potential_functions=None, # we can integrate this cleanly, not sure about anything else
            dipole_function=None
            ):
        opts = dict(
            masses=masses,
            # mass_weighted=mass_weighted,
            atoms=atoms,
            alphas=alphas,
            # projection_indices=projection_indices,
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
            # reference_structure=reference_structure,
            pairwise_potential_functions=pairwise_potential_functions,
            dipole_function=dipole_function
        )
        opts = {k:v for k,v in opts.items() if v is not None}
        ham = cls.construct(centers, potential_function, **opts)
        return ham.run()

    def run(self,
            quiet=False,
            calculate_spectrum=True,
            **wavefunction_options
            ):
        """
        The default case...

        :return:
        """

        old_logger = self.logger
        try:
            if quiet:
                logger = NullLogger()
            else:
                logger = old_logger
            with logger.block(tag="Running distributed Gaussian basis calculation"):
                wfns = self.get_wavefunctions(**wavefunction_options)
                if not quiet and isinstance(logger, NullLogger):
                    logger = Logger.lookup(True)
                logger.log_print('ZPE: {zpe}', zpe=wfns.energies[0] * UnitsData.hartrees_to_wavenumbers)
                freqs = wfns.frequencies()
                if len(freqs) > 10:
                    freqs = freqs[:10]
                logger.log_print('Frequencies: {freqs}', freqs=freqs * UnitsData.hartrees_to_wavenumbers)
                if calculate_spectrum and wfns.dipole_function is not None:
                    spec = wfns.get_spectrum()
                    ints = spec.intensities
                    if len(ints) > 10:
                        ints = ints[:10]
                    logger.log_print('Intensities: {ints}', ints=ints)
                else:
                    spec = None
        finally:
            self.logger = old_logger
        return wfns, spec

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
                  # projection_indices=None,
                  transformations=None,
                  internals=None,
                  modes=None,
                  coordinate_selection=None,
                  cartesians=None,
                  logger=False,
                  optimize_centers=False,
                  quadrature_degree=3,
                  expansion_degree=None,
                  expansion_type='multicenter',
                  # reference_structure=None,
                  poly_coeffs=None,
                  pairwise_potential_functions=None,
                  dipole_function=None
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
            coordinate_selection=coordinate_selection,
            cartesians=cartesians,
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
            logger=logger,
            wavefunction_options={'dipole_function':dipole_function}
        )

    @classmethod
    def construct_gaussians(cls,
                            centers,
                            alphas,
                            potential_spec,
                            gmat_function=None,
                            masses=None,
                            atoms=None,
                            internals=None,
                            coordinate_selection=None,
                            cartesians=None,
                            modes=None,
                            transformations=None,
                            # projection_indices=projection_indices,
                            poly_coeffs=None,
                            logger=None
                            ):
        # here to be overridden
        if (
                (isinstance(potential_spec, dict) and 'values' in potential_spec)
                or (not callable(potential_spec) and all(isinstance(p, np.ndarray) for p in potential_spec))
        ):
            potential_function = None
            potential_expansion = potential_spec
        else:
            potential_function = potential_spec
            potential_expansion = None
        return DGBGaussians.construct(
            centers,
            alphas,
            potential_function=potential_function,
            potential_expansion=potential_expansion,
            gmat_function=gmat_function,
            masses=masses,
            atoms=atoms,
            internals=internals,
            coordinate_selection=coordinate_selection,
            cartesians=cartesians,
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
                 logger=None,
                 wavefunction_options=None
                 ):
        self.gaussians = gaussians
        self.pot = potential
        self.wfn_opts = wavefunction_options
        self._V = None
        self._T = None
        self._S = None
        self.logger = Logger.lookup(logger)
    def as_cartesian_dgb(self):
        if isinstance(self.gaussians.coords, DGBCartesians):
            return self

        new_gauss = self.gaussians.as_cartesians()
        if isinstance(self.pot.potential_function, DGBCoords.DGBEmbeddedFunction):
            import copy
            new_pot = copy.copy(self.pot)
            og_fn = self.pot.potential_function.og_fn
            embed_fn = self.pot.potential_function.embed_fn
            if isinstance(embed_fn, DGBWatsonInterpolator) and isinstance(og_fn, DGBGenericInterpolator):
                og_fn = DGBCartesianWatsonInterpolator(
                    og_fn.centers,
                    og_fn.derivs,
                    embed_fn.modes,
                    **og_fn.opts
                )
            new_pot.potential_function = new_gauss.coords.embed_function(og_fn)
        else:
            new_pot = self.pot # this probably won't work...?

        return type(self)(
            new_gauss,
            new_pot,
            logger=self.logger,
            wavefunction_options=self.wfn_opts
        )

    @property
    def S(self):
        if self._S is None:
            self._S = self.gaussians.S
        return self._S
    @property
    def T(self):
        if self._T is None:
            self._T = self.S * self.gaussians.T
        return self._T
    @T.setter
    def T(self, T):
        self._T = T
    @property
    def V(self):
        if self._V is None:
            self.logger.log_print("evaluating potential energy...")
            self._V = self.S * self.pot.evaluate_pe(self.gaussians.overlap_data)
        return self._V
    @V.setter
    def V(self, V):
        self._V = V

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

    def evaluate_multiplicative_operator(self,
                                         function,
                                         expansion_degree=None,
                                         expansion_type=None,
                                         quadrature_degree=None
                                         # pairwise_functions=None # integrate out pairwise contribution
                                         ):

        pot_mat = self.pot.evaluate_op(
            function,
            self.gaussians.overlap_data,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            quadrature_degree=quadrature_degree
        )

        # with np.printoptions(linewidth=1e8):
        #     print(pot_mat)

        # if self._scaling_S is None:
        #     _ = self.S
        # S = self._scaling_S
        S = self.S
        for _ in range(pot_mat.ndim - 2):
            S = np.expand_dims(S, -1)
        return S * pot_mat

    default_solver_mode = 'cholesky'
    def diagonalize(self,
                    *,
                    mode=None,
                    similarity_cutoff=None,
                    similarity_chunk_size=None,
                    similar_det_cutoff=None,
                    subspace_size=None,
                    min_singular_value=None,
                    nodeless_ground_state=True,
                    low_rank_energy_cutoff=None,
                    low_rank_overlap_cutoff=None,
                    low_rank_shift=None,
                    stable_eigenvalue_epsilon=None
                    ):

        if mode is None:
            if any(x is not None for x in [
                subspace_size,
                min_singular_value
            ]):
                mode = 'classic'
            elif any(x is not None for x in [
                low_rank_energy_cutoff,
                low_rank_overlap_cutoff,
                low_rank_shift,
            ]):
                mode = 'low-rank'
            elif stable_eigenvalue_epsilon is not None:
                mode = 'fix-heiberger'
            elif any(x is not None for x in [
                similarity_cutoff,
                similarity_chunk_size,
                similar_det_cutoff
            ]):
                mode = 'similarity'
            else:
                mode = self.default_solver_mode


        H = self.T + self.V

        if mode == "classic":
            if min_singular_value is not None:
                self.logger.log_print("solving with min_singular_value={ms}", ms=min_singular_value)
            eigs, evecs = DGBEigensolver.classic_eigensolver(
                H, self.S, self,
                min_singular_value=min_singular_value,
                subspace_size=subspace_size,
                nodeless_ground_state=nodeless_ground_state
            )

        elif mode == 'fix-heiberger': # Implementation of the Fix-Heiberger algorithm
            eigs, evecs = DGBEigensolver.fix_heiberger(H, self.S, self, eps=stable_eigenvalue_epsilon)

        elif mode == 'similarity':
            eigs, evecs = DGBEigensolver.similarity_mapped_solver(H, self.S, self,
                                                                  similarity_cutoff=similarity_cutoff,
                                                                  similarity_chunk_size=similarity_chunk_size,
                                                                  similar_det_cutoff=similar_det_cutoff
                                                                  )

        elif mode == 'low-rank':
            eigs, evecs = DGBEigensolver.low_rank_solver(H, self.S, self,
                                                         low_rank_energy_cutoff=low_rank_energy_cutoff,
                                                         low_rank_overlap_cutoff=low_rank_overlap_cutoff,
                                                         low_rank_shift=low_rank_shift
                                                         )

        elif mode == 'cholesky':
            eigs, evecs = DGBEigensolver.cholesky_solver(H, self.S, self)

        elif callable(mode):
            eigs, evecs = mode(H, self.S, )
        else:
            raise ValueError("unknown solver {}".format(mode))

                # raise Exception((eigs[1:] - eigs[0])*219475.6)


        return eigs, evecs

    def get_wavefunctions(self,
                          mode=None,
                          similarity_cutoff=None,
                          similarity_chunk_size=None,
                          similar_det_cutoff=None,
                          subspace_size=None,
                          min_singular_value=None,
                          nodeless_ground_state=None,
                          stable_eigenvalue_epsilon=None,
                          **wfn_opts
                          ):
        # print("======="*25)
        ops = dict(
            mode=mode,
            similarity_cutoff=similarity_cutoff,
            similarity_chunk_size=similarity_chunk_size,
            similar_det_cutoff=similar_det_cutoff,
            subspace_size=subspace_size,
            min_singular_value=min_singular_value,
            nodeless_ground_state=nodeless_ground_state,
            stable_eigenvalue_epsilon=stable_eigenvalue_epsilon
        )
        ops = {k:v for k,v in ops.items() if v is not None}
        eigs, evecs = self.diagonalize(**ops)
        if self.wfn_opts is not None:
            wfn_opts = dict(wfn_opts, **self.wfn_opts)
        return DGBWavefunctions(eigs, evecs, hamiltonian=self, **wfn_opts)