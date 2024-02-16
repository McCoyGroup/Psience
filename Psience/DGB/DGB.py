
import numpy as np, scipy as sp

from McUtils.Data import UnitsData
from McUtils.Zachary import DensePolynomial
from McUtils.Scaffolding import Logger, NullLogger
from McUtils.Parallelizers import Parallelizer

from .Gaussians import *
from .Evaluators import *
from .Coordinates import *
from .Interpolation import *
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
            dipole_degree=0,
            num_print=20,
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
                if len(freqs) > num_print:
                    freqs = freqs[:num_print]
                logger.log_print('Frequencies: {freqs}', freqs=freqs * UnitsData.hartrees_to_wavenumbers)
                if calculate_spectrum and wfns.dipole_function is not None:
                    spec = wfns.get_spectrum(
                        expansion_degree=dipole_degree
                    )
                    ints = spec.intensities
                    if len(ints) > num_print:
                        ints = ints[:num_print]
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
                  parallelizer=None,
                  optimize_centers=False,
                  quadrature_degree=3,
                  expansion_degree=None,
                  expansion_type='multicenter',
                  # reference_structure=None,
                  momenta=None,
                  poly_coeffs=None,
                  pairwise_potential_functions=None,
                  dipole_function=None,
                  kinetic_options=None
                  ):
        logger = Logger.lookup(logger)
        parallelizer = Parallelizer.lookup(parallelizer)

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
            momenta=momenta,
            pairwise_potential_functions=pairwise_potential_functions,
            poly_coeffs=poly_coeffs,
            logger=logger,
            parallelizer=parallelizer,
            kinetic_options=kinetic_options
        )
        if optimize_centers:
            if not isinstance(optimize_centers, (list, tuple)):
                optimize_centers = [optimize_centers]
            for optimizer in optimize_centers:
                coords = coords.optimize(optimizer, potential_function=potential_function)

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
            parallelizer=parallelizer,
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
                            kinetic_options=None,
                            transformations=None,
                            # projection_indices=projection_indices,
                            momenta=None,
                            pairwise_potential_functions=None,
                            poly_coeffs=None,
                            logger=None,
                            parallelizer=None
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
            kinetic_options=kinetic_options,
            momenta=momenta,
            pairwise_potential_functions=pairwise_potential_functions,
            poly_coeffs=poly_coeffs,
            logger=logger,
            parallelizer=parallelizer
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
                            logger=None,
                            parallelizer=None
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
                 parallelizer=None,
                 wavefunction_options=None
                 ):
        self.gaussians = gaussians
        self.pot = potential
        self.wfn_opts = wavefunction_options
        self._V = None
        self._T = None
        self._S = None
        self.logger = Logger.lookup(logger)
        self.parallelizer = Parallelizer.lookup(parallelizer)
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

    def get_S(self):
        od = self.gaussians.overlap_data
        S = self.gaussians.S
        return S
    @property
    def S(self):
        if self._S is None:
            self._S = self.get_S()
        return self._S

    def get_T(self):
        od = self.gaussians.overlap_data
        T = self.gaussians.T
        S = self.gaussians.prefactor
        return S * T
    @property
    def T(self):
        if self._T is None:
            self._T = self.get_T()
        return self._T
    @T.setter
    def T(self, T):
        self._T = T

    def get_V(self):
        od = self.gaussians.overlap_data
        V = self.pot.evaluate_pe(od)
        S = self.gaussians.prefactor
        return S * V
        #     print('why '*25)
        #     print(S_diff[0, 0], V_diff[0, 0])
        #     print(S_sum[0, 0], V_sum[0, 0])
        #     print(S_diff[0, 0] * V_diff[0, 0],  S_sum[0, 0] * V_sum[0, 0])
        #     return (S_diff * V_diff + S_sum * V_sum)
        # else:
        # return self.S * V
    @property
    def V(self):
        if self._V is None:
            with self.logger.block(tag="Evaluating potential energy matrix"):
                self._V = self.get_V()
        return self._V
    @V.setter
    def V(self, V):
        self._V = V

    def evaluate_multiplicative_operator(self,
                                         function,
                                         embed=True,
                                         expansion_degree=None,
                                         expansion_type=None,
                                         quadrature_degree=None,
                                         pairwise_functions=None # integrate out pairwise contribution
                                         ):

        if embed:
            function = self.gaussians.coords.embed_function(function)
            if pairwise_functions is not None:
                pairwise_functions = self.gaussians.coords.pairwise_potential_evaluator(pairwise_functions)

        pot_mat = self.pot.evaluate_op(
            function,
            self.gaussians.overlap_data,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            quadrature_degree=quadrature_degree,
            pairwise_functions=pairwise_functions
        )

        S = self.S
        for _ in range(pot_mat.ndim - 2):
            S = np.expand_dims(S, -1)
        return S * pot_mat

    default_solver_mode = 'similarity'
    def diagonalize(self,
                    *,
                    mode=None,
                    similarity_cutoff=None,
                    similarity_chunk_size=None,
                    similar_det_cutoff=None,
                    similarity_shift=None,
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
                low_rank_shift
            ]):
                mode = 'low-rank'
            elif stable_eigenvalue_epsilon is not None:
                mode = 'fix-heiberger'
            elif similarity_shift is not None:
                mode = 'shift'
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
        elif mode == 'shift':
            eigs, evecs = DGBEigensolver.shift_similarity_solver(H, self.S, self,
                                                                 similarity_cutoff=similarity_cutoff,
                                                                 similarity_chunk_size=similarity_chunk_size,
                                                                 similar_det_cutoff=similar_det_cutoff,
                                                                 similarity_shift=similarity_shift
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

    def get_similarity_matrix(self):
        eigs, Qs = np.linalg.eigh(self.S)
        eigh, Qh = np.linalg.eigh(self.V + self.T)

        similarity_matrix = Qs.T @ Qh

        # import McUtils.Plots as plt
        # plt.MatrixPlot(similarity_matrix).show()
        return similarity_matrix

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