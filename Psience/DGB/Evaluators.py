
import abc, math, gc

import numpy as np, itertools, functools

from McUtils.Zachary import DensePolynomial, TensorDerivativeConverter
from McUtils.Scaffolding import Logger
from McUtils.Parallelizers import Parallelizer
import McUtils.Numputils as nput

__all__ = [
    "DGBEvaluator",
    "DGBKineticEnergyEvaluator",
    "DGBCartesianEvaluator",
    "DGBWatsonEvaluator",
    "DGBPotentialEnergyEvaluator",
    "DGBPairwisePotentialEvaluator",
    "DGBCartesianPairwiseEvaluator",
    "DGBWatsonPairwiseEvaluator"
]


class OverlapGaussianData:
    def __init__(self,
                 input_data,
                 product_data
                 ):
        self._input_data = input_data
        self._product_data = product_data
    @property
    def npts(self):
        return self.initial_centers.shape[0]
    @property
    def ndim(self):
        return self.initial_centers.shape[1]
    @property
    def centers(self):
        return self._product_data['centers']
    @property
    def alphas(self):
        return self._product_data['alphas']
    @property
    def rotations(self):
        return self._product_data['rotations']
    @property
    def inverse_rotations(self):
        return self._product_data['inverse_rotations']
    @property
    def precision_matrices(self):
        return self._product_data['covariances']
    @property
    def covariance_matrices(self):
        return self._product_data['inverse']
    @property
    def shift_matrices(self):
        return self._product_data['sum_inverse']
    @property
    def phases_diff(self):
        return self._product_data['phases_diff']
    @property
    def phases_sum(self):
        return self._product_data['phases_sum']
    @property
    def momenta_sum(self):
        return self._product_data['diagonal_momentum_sum']
    @property
    def momenta_diff(self):
        return self._product_data['diagonal_momentum_diff']
    @property
    def rho_sum(self):
        return self._product_data['rho_sum']
    @property
    def rho_diff(self):
        return self._product_data['rho_diff']
    @property
    def correlation_factors_sum(self):
        return self.cos_correlation_sum, self.sin_correlation_sum
    @property
    def correlation_factors_diff(self):
        return self.cos_correlation_diff, self.sin_correlation_diff
    @property
    def cos_correlation_diff(self):
        return self._product_data['cos_correlation_diff']
    @property
    def sin_correlation_diff(self):
        return self._product_data['sin_correlation_diff']
    @property
    def cos_correlation_sum(self):
        return self._product_data['cos_correlation_sum']
    @property
    def sin_correlation_sum(self):
        return self._product_data['sin_correlation_sum']
    @property
    def decay_factor_diff(self):
        return self._product_data['phase_diff_decay']
    @property
    def decay_factor_sum(self):
        return self._product_data['phase_sum_decay']
    @property
    def center_difference(self):
        return self._product_data['center_diffs']
    @property
    def delta_position(self):
        return self._product_data['center_delta']
    @property
    def delta_phase_sum(self):
        return self._product_data['phase_delta_sum']
    @property
    def delta_phase_diff(self):
        return self._product_data['phase_delta_diff']
    @property
    def indices(self):
        return self.row_indices, self.col_indices
    @property
    def row_indices(self):
        return self._product_data['row_inds']
    @property
    def col_indices(self):
        return self._product_data['col_inds']
    @property
    def initial_centers(self):
        return self._input_data['init_centers']
    @property
    def initial_alphas(self):
        return self._input_data['init_alphas']
    @property
    def initial_precision_matrices(self):
        return self._input_data['init_covariances']
    @property
    def initial_momenta(self):
        return self._input_data['initial_momenta']
    @property
    def initial_phases(self):
        return self._input_data['initial_phases']

    @classmethod
    def from_gaussian_parameters(cls,
                                 centers, alphas, transformations,
                                 momenta,
                                 *,
                                 chunk_size=None,  # int(1e6),
                                 rows_cols=None,
                                 logger=None,
                                 parallelizer=None
                                 ):
        return cls(*cls.get_overlap_data(
            centers, alphas, transformations,
            momenta,
            chunk_size=chunk_size,  # int(1e6),
            rows_cols=rows_cols,
            logger=logger,
            parallelizer=parallelizer
        ))

    @classmethod
    def get_overlap_data(clget_s,
                          centers, alphas, transformations,
                          aligned_momenta,
                          *,
                          chunk_size=None,  # int(1e6),
                          rows_cols=None,
                          logger=None,
                          parallelizer=None
                          ):
        parallelizer = Parallelizer.lookup(parallelizer)
        logger = Logger.lookup(logger)

        if rows_cols is None:
            rows, cols = np.triu_indices(len(alphas))
        else:
            rows, cols = rows_cols

        # a bit inefficient if all tfs are identity, but generality is good
        sigs = DGBEvaluator.get_inverse_covariances(alphas, transformations)
        covs = DGBEvaluator.get_covariances(alphas, transformations)
        momenta = DGBEvaluator.get_phase_vectors(aligned_momenta, transformations)
        new_sigs = sigs[rows] + sigs[cols]

        if chunk_size is not None:
            num_segments = new_sigs.shape[0] // chunk_size + 1
            chunks = np.array_split(new_sigs, num_segments)
            row_chunks = np.array_split(rows, num_segments)
            col_chunks = np.array_split(cols, num_segments)
        else:
            chunks = [new_sigs]
            row_chunks = [rows]
            col_chunks = [cols]

        logger.log_print("getting {nT} overlap Gaussians over {nC} chunks", nT=len(rows), nC=len(chunks))
        chunk_data = []
        for news, r, c in zip(chunks, row_chunks, col_chunks):
            new_alphas, new_rots = np.linalg.eigh(news)  # eigenvalues of inverse tensor...
            new_rots_inv = new_rots.transpose(0, 2, 1)

            # I _could_ construct the inverse from the alphas and rotations
            # but I think it makes more sense to use a potentially more stable
            # inverse here...also it is apparently no slower...
            new_inv = np.linalg.inv(news)
            # new_inv = new_rots @ nput.vec_tensordiag(1 / new_alphas) @ new_rots.transpose(0, 2, 1)

            new_centers = new_inv @ (
                    sigs[r] @ centers[r][:, :, np.newaxis]
                    + sigs[c] @ centers[c][:, :, np.newaxis]
            )
            new_centers = new_centers.reshape(centers[c].shape)
            new_alphas = new_alphas / 2
            sum_sigs = sigs[r] @ new_inv @ sigs[c]

            # TODO: turn this into a proper object...
            chunk_data.append({
                'centers': new_centers,
                'alphas': new_alphas,
                'covariances': new_sigs,
                'inverse': new_inv,
                'sum_inverse': sum_sigs,
                'rotations': new_rots,
                'inverse_rotations': new_rots_inv
            })


        product_data = {
            k: None if chunk_data[0][k] is None else np.concatenate([d[k] for d in chunk_data])
            for k in chunk_data[0].keys()
        }

        new_cov = product_data['inverse']
        new_si = product_data['sum_inverse']
        if momenta is not None:
            new_centers = product_data['centers']
            new_rots = product_data['rotations']
            new_alphas = product_data['alphas']

            # expressed in terms of the global frame
            phases_sum = momenta[rows] + momenta[cols]
            phases_diff = momenta[rows] - momenta[cols]
            rho_sum = np.reshape(new_cov @ phases_sum[:, :, np.newaxis], phases_sum.shape)
            rho_diff = np.reshape(new_cov @ phases_diff[:, :, np.newaxis], phases_sum.shape)
            # dot momentum sum or difference into displacement from origin
            corr_sum = np.reshape(new_centers[:, np.newaxis, :] @ phases_sum[:, :, np.newaxis], -1)
            corr_diff = np.reshape(new_centers[:, np.newaxis, :] @ phases_diff[:, :, np.newaxis], -1)
            # compute the J+ terms we use for all derivations
            scj = covs @ momenta[:, :, np.newaxis]
            delta_j_sum = np.reshape(new_si @ (scj[rows] - scj[cols]), phases_sum.shape)
            delta_j_diff = np.reshape(new_si @ (scj[rows] + scj[cols]), phases_diff.shape)

            # now express these in terms of the rotated frame

            mom_sum = DGBEvaluator.get_momentum_vectors(phases_sum, (new_rots, new_rots.transpose(0, 2, 1)))
            mom_diff = DGBEvaluator.get_momentum_vectors(phases_diff, (new_rots, new_rots.transpose(0, 2, 1)))

            decay_sum = np.sum(mom_sum ** 2 / (4 * new_alphas), axis=1)
            decay_diff = np.sum(mom_diff ** 2 / (4 * new_alphas), axis=1)

            sum_prefac = np.exp(-decay_sum)  # * np.cos(phase_sum)
            diff_prefac = np.exp(-decay_diff)  # * np.cos(phase_diff)
        else:
            mom_sum = None
            mom_diff = None
            phases_sum = None
            phases_diff = None
            rho_sum = None
            rho_diff = None
            corr_sum = None
            corr_diff = None
            sum_prefac = None
            diff_prefac = None

            delta_j_sum = None
            delta_j_diff = None

        row_gammas = new_cov @ sigs[rows]
        col_gammas = new_cov @ sigs[cols]
        disps = centers[rows] - centers[cols]
        delta = np.reshape(new_si @ disps[:, :, np.newaxis], disps.shape)
        product_data.update({
            'row_gammas':row_gammas,
            'col_gammas':col_gammas,
            'gamma_diff':row_gammas - col_gammas,
            'center_diffs':disps,
            'center_delta':delta,
            'diagonal_momentum_sum': mom_sum,
            'diagonal_momentum_diff': mom_diff,
            'phases_sum': phases_sum,
            'phases_diff': phases_diff,
            'rho_sum': rho_sum,
            'rho_diff': rho_diff,
            'phase_sum_corr': corr_sum,
            'phase_diff_corr': corr_diff,
            'phase_sum_decay': sum_prefac,
            'phase_diff_decay': diff_prefac,
            'phase_delta_sum':delta_j_sum,
            'phase_delta_diff':delta_j_diff,
            'row_inds': rows,
            'col_inds': cols,
        })
        if mom_sum is not None:
            product_data.update({
                'cos_correlation_sum': np.cos(corr_sum),
                'cos_correlation_diff': np.cos(corr_diff),
                'sin_correlation_sum': np.sin(corr_sum),
                'sin_correlation_diff': np.sin(corr_diff)
            })
        input_data = {
            'init_centers': centers,
            'init_alphas': alphas,
            'init_covariances': sigs,
            'initial_phases': momenta,
            'initial_momenta': aligned_momenta
        }

        return input_data, product_data

    def take_subselection(self, positions):
        # we need to find the positions in the old stuff where we
        rows, cols = self.indices

        mask_1, _, _ = nput.contained(rows, positions)
        mask_2, _, _ = nput.contained(cols, positions)
        take_pos = np.where(np.logical_and(mask_1 > 0, mask_2 > 0))[0]
        new_prod = {
            k: v[take_pos] if v is not None else v
            for k, v in self._product_data.items()
        }
        rows, cols = np.triu_indices(len(positions))
        new_prod['row_inds'] = rows
        new_prod['col_inds'] = cols

        new_input = {
            k: v[positions] if v is not None else v
            for k, v in self._input_data.items()
        }

        return type(self)(new_input, new_prod)
class DGBEvaluator:
    """
    An object that supports evaluating matrix elements in a distributed Gaussian basis.
    Provides support for integrating a function via quadrature or as an expansion in a polynomial tensors
    """

    @classmethod
    def get_inverse_covariances(cls, alphas, transformations):
        """
        Transforms the alphas into proper inverse covariance matrices.
        Chosen so that in the case that the transformations, Q, diagonalize S we can write
            QT S Q = A

        :return:
        :rtype:
        """

        n = alphas.shape[-1]
        npts = len(alphas)
        covs = np.zeros((npts, n, n))
        diag_inds = (slice(None, None, None),) + np.diag_indices(n)
        covs[diag_inds] = 2 * alphas

        if transformations is not None:
            tfs, inv = transformations
            covs = tfs @ covs @ tfs.transpose((0, 2, 1))

        return covs

    @classmethod
    def get_covariances(cls, alphas, transformations):
        """
        Transforms the alphas into proper inverse covariance matrices.
        Chosen so that in the case that the transformations, Q, diagonalize S we can write
            QT S Q = A

        :return:
        :rtype:
        """

        n = alphas.shape[-1]
        npts = len(alphas)
        covs = np.zeros((npts, n, n))
        diag_inds = (slice(None, None, None),) + np.diag_indices(n)
        covs[diag_inds] = 1/(2 * alphas)

        if transformations is not None:
            tfs, inv = transformations
            covs = tfs @ covs @ tfs.transpose((0, 2, 1))

        return covs

    @classmethod
    def get_momentum_vectors(cls, phases, transformations):
        """
        Transforms the momenta so that they're aligned along the Gaussian axes

        :return:
        :rtype:
        """

        if phases is not None and transformations is not None:
            tfs, inv = transformations
            momenta = np.reshape(
                phases[:, np.newaxis, :] @ inv.transpose((0, 2, 1)),
                phases.shape
            )
        else:
            momenta = phases

        return momenta

    @classmethod
    def get_phase_vectors(cls, momenta, transformations):
        """
        Transforms the momenta so that they're aligned along the Gaussian axes

        :return:
        :rtype:
        """

        if momenta is not None and transformations is not None:
            tfs, inv = transformations
            phases = np.reshape(
                momenta[:, np.newaxis, :] @ tfs.transpose((0, 2, 1)),
                momenta.shape
            )
        else:
            phases = momenta

        return phases

    @classmethod
    def get_overlap_gaussians(cls,
                              centers, alphas, transformations,
                              momenta,
                              *,
                              chunk_size=None,#int(1e6),
                              rows_cols=None,
                              logger=None,
                              parallelizer=None
                              ) -> 'OverlapGaussianData':
        return OverlapGaussianData.from_gaussian_parameters(
            centers, alphas, transformations,
            momenta,
            chunk_size=chunk_size,  # int(1e6),
            rows_cols=rows_cols,
            logger=logger,
            parallelizer=parallelizer
        )

    @classmethod
    def poch(cls, n, m): #pochammer/generalized gamma
        nums = np.arange(n-2*m+1 if m < n/3 else m + 1, n+1)
        dens = np.arange(1, m+1 if m < n/3 else n-2*m+1)
        if len(dens) < len(nums): # pad on left so we can have most stable eval
            dens = np.concatenate([
                np.ones(len(nums) - len(dens)),
                dens
            ])
        elif len(nums) < len(dens):
            nums = np.concatenate([
                np.ones(len(dens) - len(nums)),
                nums
            ])
        return np.prod(nums/dens)
    @classmethod
    def polyint_1D(cls, centers, alphas, n):
        if n == 0:
            return np.ones(centers.shape[:2])
        c = np.sqrt(alphas) * centers
        term = sum(
            (cls.poch(n, l) * (1/2**(2*l-n) if 2*l > n else 2**(n-2*l)))*(c)**(n-2*l)
            for l in range(0, int(np.floor(n/2)) + 1)
        )
        return term

    @classmethod
    def momentum_coeffient(cls, k, n):
        m = k//2
        s = (-1) ** k
        fac = (-1)**(n) * np.prod(
            2*np.arange(n+1, m+1, 2) - s
        )
        bin = np.math.comb(m, n)
        return bin * fac
    @classmethod
    def momentum_integral(cls, p, a, k):
        var = (p**2)/(2*a) # the rest of the sqrt(a) term is included elsewhere
        expansion = sum(
            cls.momentum_coeffient(k, n) * (var**n)
            for n in range(0, int(np.ceil(k/2)) + 1)
        ) / (2*a)**np.ceil(k/2)
        if k%2 == 1:
            expansion = p * expansion
        return expansion
    @classmethod
    def simple_poly_int(cls, n):
        return np.prod(np.arange(1, n, 2)) / 2**(n/2) # double factorial/gamma/whatever
    @classmethod
    def tensor_expansion_integrate(cls,
                                   npts, derivs, overlap_data:'OverlapGaussianData',
                                   expansion_type='multicenter',
                                   logger=None, reweight=True
                                   ):
        """
        provides an integral from a polynomial expansion with derivs as an expansion in tensors

        :param npts:
        :param derivs:
        :param centers:
        :param alphas:
        :param inds:
        :param rot_data:
        :param expansion_type:
        :param logger:
        :return:
        """

        centers = overlap_data.centers
        ndim = centers.shape[-1]
        rotations = overlap_data.rotations
        if expansion_type == 'taylor':
            centers = np.reshape(rotations @ centers[:, :, np.newaxis], centers.shape)
        alphas = overlap_data.alphas
        m_diff = overlap_data.momenta_diff # rotated along the Gaussian axes like the alphas
        m_sum = overlap_data.momenta_sum

        # if self.
        fdim = derivs[0].ndim - 1 # to simplify dipole functions
        fshape = derivs[0].shape[1:]

        new_derivs = []
        # rotations = rotations[:, :, :, np.newaxis] # to test shapes
        for n,d in enumerate(derivs):
            for _ in range(n):
                d = nput.vec_tensordot(
                    d, rotations,
                    axes=[fdim+1, 1],
                    shared=1
                )
            new_derivs.append(d)
        derivs = new_derivs

        if logger is not None:
            logger.log_print("adding up all derivative contributions...")

        row_inds = overlap_data.row_indices
        col_inds = overlap_data.col_indices

        pot = np.zeros((npts, npts) + fshape)
        caches = [{} for _ in range(ndim)]
        # mom_sum_pot = np.zeros((npts, npts) + fshape) # for momentum sum
        ms_caches = [{} for _ in range(ndim)]
        for nd,d in enumerate(derivs): # add up all independent integral contribs...
            # iterate over upper triangle coordinates (we'll add bottom contrib by symmetry)
            inds = itertools.combinations_with_replacement(range(ndim), r=nd) if nd > 0 else [()]
            for idx in inds:
                count_map = {k: v for k, v in zip(*np.unique(idx, return_counts=True))}
                if (
                        expansion_type != 'taylor' and m_sum is None
                        and any(n%2 !=0 for n in count_map.values())
                ):
                    continue # odd contribs vanish

                contrib = 1
                ms_contrib = 1
                for k in range(ndim): # do each dimension of integral independently
                    n = count_map.get(k, 0)
                    if n not in caches[k]:
                        if expansion_type == 'taylor':
                            raise NotImplementedError("Taylor series in rotated basis not implemented yet")

                        if m_sum is not None:
                            caches[k][n] = cls.momentum_integral(
                                m_diff[..., k],
                                alphas[..., k],
                                n
                            )
                            ms_caches[k][n] = cls.momentum_integral(
                                m_sum[..., k],
                                alphas[..., k],
                                n
                            )
                        else:
                            caches[k][n] = (
                               cls.simple_poly_int(n)
                                 if expansion_type != 'taylor' else
                               cls.polyint_1D(centers[..., k], alphas[..., k], n)
                            ) / alphas[..., k]**np.ceil(n/2)
                    base_contrib = caches[k][n]
                    for _ in range(fdim):
                        base_contrib = np.expand_dims(base_contrib, -1)
                    if isinstance(contrib, int) and fdim > 0:
                        base_contrib = np.broadcast_to(base_contrib,  alphas.shape[:-1] + fshape)
                    contrib *= base_contrib

                    if m_sum is not None:
                        base_contrib = ms_caches[k][n]
                        for _ in range(fdim):
                            base_contrib = np.expand_dims(base_contrib, -1)
                        if isinstance(ms_contrib, int) and fdim > 0:
                            base_contrib = np.broadcast_to(base_contrib, alphas.shape[:-1] + fshape)
                        ms_contrib *= base_contrib

                dcont = d[(slice(None, None, None),)*(fdim+1) + idx] if len(idx) > 0 else d
                if reweight:
                    # compute multinomial coefficient for weighting purposes
                    _, counts = np.unique(idx, return_counts=True)
                    multicoeff = 1
                    for x in counts: multicoeff *= math.factorial(x)
                    multicoeff = multicoeff / math.factorial(len(idx))
                    scaling = multicoeff / math.factorial(nd)
                    contrib *= scaling * dcont
                    if m_sum is not None:
                        ms_contrib *= scaling * dcont
                else:
                    contrib *= dcont
                    if m_sum is not None:
                        ms_contrib *= dcont

                if m_sum is not None:
                    sum_prefac = overlap_data.decay_factor_sum
                    diff_prefac = overlap_data.decay_factor_diff

                    contrib *= np.expand_dims(
                        overlap_data.cos_correlation_diff if nd%2 == 0 else -overlap_data.sin_correlation_diff,
                        [-x for x in range(1, fdim+1)]
                    )
                    ms_contrib *= np.expand_dims(
                        overlap_data.cos_correlation_sum if nd%2 == 0 else -overlap_data.sin_correlation_sum,
                        [-x for x in range(1, fdim+1)]
                    )

                    if ms_contrib.ndim > 1:
                        sum_prefac = sum_prefac[:, np.newaxis]
                        diff_prefac = diff_prefac[:, np.newaxis]

                    ms_contrib *= sum_prefac
                    contrib *= diff_prefac
                    # print("dd>", contrib[:5])

                    contrib = contrib + ms_contrib

                pot[row_inds, col_inds] += contrib

        pot[col_inds, row_inds] = pot[row_inds, col_inds]
        # mom_sum_pot[col_inds, row_inds] = mom_sum_pot[row_inds, col_inds]
        #
        # if m_sum is not None:
        #     return pot, mom_sum_pot
        # else:
        return pot

    @classmethod
    def quad_weight_eval(cls, function, d_chunk, w_chunk, ndim, centers, squa):
        d = d_chunk[np.newaxis] / squa[:, np.newaxis]
        c = centers[:, np.newaxis, :] + d
        fv = function(c.reshape(-1, ndim))
        if isinstance(fv, list):
            chunk_vals = []
            fshape = fv[0].shape[1:]
            w_chunk = np.expand_dims(w_chunk, [0] + [-x for x in range(1, len(fshape) + 1)])
            for mv in fv:
                mv = mv.reshape(c.shape[:2] + fshape)
                chunk_vals.append(np.sum(w_chunk * mv, axis=1))
            return chunk_vals
        else:
            fshape = fv.shape[1:]
            fv = fv.reshape(c.shape[:2] + fshape)
            w_chunk = np.expand_dims(w_chunk, [0] + [-x for x in range(1, len(fshape) + 1)])
            chunk_val = np.sum(w_chunk * fv, axis=1)
            return chunk_val
    @classmethod
    def quad_nd(cls,
                centers, alphas, function,
                flatten=False,
                degree=3, chunk_size=int(1e6),
                normalize=True
                ):
        """
        N-dimensional quadrature

        :param centers:
        :param alphas:
        :param function:
        :param degree:
        :return:
        """

        centers = np.asanyarray(centers)
        if centers.ndim == 1:
            centers = centers[np.newaxis]
            alphas = [alphas]
        alphas = np.asanyarray(alphas)

        # Quadrature point displacements and weights (thanks NumPy!)
        disps, weights = np.polynomial.hermite.hermgauss(degree)

        ndim = centers.shape[-1]
        indices = np.moveaxis(np.array(
                np.meshgrid(*([np.arange(0, degree, dtype=int)] * ndim))
            ), 0, -1).reshape(-1, ndim)
        disps = disps[indices]

        w = np.prod(weights[indices], axis=-1)
        n_disps = alphas.shape[0] * disps.shape[0] # how many structures total
        num_segments = n_disps // chunk_size + 1
        disp_chunks = np.array_split(disps, num_segments)
        w_chunks = np.array_split(w, num_segments)

        val = None
        squa = np.sqrt(alphas)
        for d_chunk, w_chunk in zip(disp_chunks, w_chunks):
            chunk_val = cls.quad_weight_eval(function, d_chunk, w_chunk, ndim, centers, squa)
            if val is None:
                val = chunk_val
            else:
                if isinstance(val, list):
                    for i in range(len(val)):
                        val[i] += chunk_val[i]
                else:
                    val += chunk_val

        if normalize:
            normalization = 1 / np.prod(squa, axis=-1)
            if isinstance(val, list):
                for i in range(len(val)):
                    val[i] *= normalization
            else:
                val *= normalization
        return val

    @classmethod
    def _wrap_rotated_function(cls, func, rotations, inverse, momenta):

        # momenta = momenta@rotations
        def eval(coords, deriv_order=None, inverse=inverse, rotations=rotations):
            # coords are expressed in the rotated space, we need to transform back out
            quad_disps = coords.shape[0] // inverse.shape[0]
            inverse = np.broadcast_to(
                inverse[:, np.newaxis, :, :],
                (inverse.shape[0], quad_disps) + inverse.shape[1:]
            ).reshape((coords.shape[0],) + inverse.shape[1:])
            unrot_coords = np.reshape(inverse @ coords[:, :, np.newaxis], coords.shape)
            fdat = func(unrot_coords, deriv_order=deriv_order)
            if deriv_order is not None:
                rotations = np.broadcast_to(
                    rotations[:, np.newaxis, :, :],
                    (rotations.shape[0], quad_disps) + rotations.shape[1:]
                ).reshape((coords.shape[0],) + rotations.shape[1:])

                f_vals = fdat[0]
                f_derivs = []
                for n, d in enumerate(fdat[1:]):
                    # now we want to rexpress the derivatives in the rotated space (I think)
                    for _ in range(n):
                        d = nput.vec_tensordot(
                            d, rotations,
                            axes=[1, 1],
                            shared=1
                        )
                    f_derivs.append(d)
                fdat = [f_vals] + f_derivs
            if momenta is not None:
                # new_dat = []
                fvals = fdat if isinstance(fdat, np.ndarray) else fdat[0]
                cos_cor = 0
                for p in momenta:
                    p = np.broadcast_to(
                        p[:, np.newaxis, :],
                        (p.shape[0], quad_disps) + p.shape[1:]
                    ).reshape((unrot_coords.shape[0],) + p.shape[1:])
                    correlation = np.reshape(
                        np.matmul(
                            unrot_coords[:, np.newaxis, :],
                            p[:, :, np.newaxis]
                        ),
                        -1
                    )
                    cos_cor += np.cos(correlation)
                new_dat = cos_cor * fvals
                if isinstance(fdat, np.ndarray):
                    fdat = new_dat
                else:
                    fdat[0] = new_dat
            return fdat

        return eval

    @classmethod
    def rotated_gaussian_quadrature(cls,
                                    function,
                                    alphas,
                                    centers,
                                    rotations,
                                    inverse,
                                    momenta,
                                    normalize=True,
                                    degree=2,
                                    ):
        f = cls._wrap_rotated_function(function, rotations, inverse, momenta)
        centers = np.reshape(rotations @ centers[:, :, np.newaxis], centers.shape)

        vals = cls.quad_nd(centers, alphas, f, degree=degree, normalize=normalize)

        normalization = 1 / (np.sqrt(np.pi)) ** centers.shape[-1]  # part of the prefactor...
        vals *= normalization
            # diff_vals, sum_vals = vals
            # diff_phase, sum_phase = momenta
            # diff_corr = np.reshape(
            #     np.matmul(centers[:, np.newaxis, :], diff_phase[:, :, np.newaxis]),
            #     -1
            # )
            # sum_corr = np.reshape(
            #     np.matmul(centers[:, np.newaxis, :], sum_phase[:, :, np.newaxis]),
            #     -1
            # )
            # print(sum_vals)
            # vals = [
            #     diff_vals * np.cos(diff_corr),
            #     sum_vals * np.cos(sum_corr)
            # ]


        return vals

    @classmethod
    def quad_integrate(cls, function, overlap_data:'OverlapGaussianData', degree=2, logger=None):
        """
        Integrate potential over all pairs of Gaussians at once

        :param degree:
        :type degree:
        :return:
        :rtype:
        """

        use_momenta = overlap_data.momenta_diff is not None
        vals = cls.rotated_gaussian_quadrature(
            function,
            overlap_data.alphas,
            overlap_data.centers,
            overlap_data.inverse_rotations,
            overlap_data.rotations,
            [overlap_data.phases_diff, overlap_data.phases_sum]
                if use_momenta else None,
            degree=degree,
            normalize=False
        )

        npts = overlap_data.npts
        rows, cols = np.triu_indices(npts)

        pots = np.zeros((npts, npts) + vals.shape[1:])
        pots[rows, cols] = vals
        pots[cols, rows] = vals

        return pots

    @classmethod
    def _rot_base_S_components(cls, overlap_data:'OverlapGaussianData'):
        row_inds = overlap_data.row_indices
        col_inds = overlap_data.col_indices
        alphas = overlap_data.initial_alphas
        rotas = overlap_data.alphas

        ndim = overlap_data.ndim
        # we prefactor out the 2**ndim

        det_rat = np.prod(
            np.sort(alphas[row_inds]/rotas, axis=1)
            * np.flip(np.sort(alphas[col_inds]/rotas, axis=1), axis=1),
            axis=1
        )

        thing1 = np.prod(alphas[row_inds]/rotas, axis=1)
        n = len(alphas)
        rat_mat = np.zeros((n, n))
        rat_mat[row_inds, col_inds] = thing1

        disps = overlap_data.center_difference
        C = disps[:, np.newaxis, :] @ overlap_data.shift_matrices @ disps[:, :, np.newaxis]
        C = C.reshape(disps.shape[0])

        return row_inds, col_inds, ndim, det_rat, C

    @classmethod
    def _rot_base_S(cls, overlap_data:'OverlapGaussianData', logger=None, return_prefactor=False):

        logger.log_print('evaluating {nT} overlaps', nT=len(overlap_data.row_indices))
        n = overlap_data.npts
        S = np.eye(n)
        row_inds, col_inds, ndim, det_rat, C = cls._rot_base_S_components(overlap_data)

        scaled_rats = ( 2**(ndim/2) ) * ( (det_rat) ** (1 / 4) )
        exp_terms = np.exp(-C / 2)
        full_terms = scaled_rats * exp_terms

        # include phase factors
        if overlap_data.initial_momenta is not None:
            j_vecs = overlap_data.initial_phases

            phases = np.reshape(overlap_data.initial_centers[:, np.newaxis, :] @ j_vecs[:, :, np.newaxis], -1)
            decay_norm = np.sum(overlap_data.initial_momenta**2 / (2*overlap_data.initial_alphas), axis=1)
            norm = np.sqrt((1 + np.cos(2*phases)*np.exp(-decay_norm)))
            norms = norm[row_inds] * norm[col_inds]

            full_terms = full_terms / norms

            prefac = np.eye(n)
            prefac[row_inds, col_inds] = full_terms
            prefac[col_inds, row_inds] = prefac[row_inds, col_inds]

            decay_contrib = (
                    overlap_data.decay_factor_diff * overlap_data.cos_correlation_diff
                    + overlap_data.decay_factor_sum * overlap_data.cos_correlation_sum
            )
            S[row_inds, col_inds] = decay_contrib
            S[col_inds, row_inds] = decay_contrib
            S = S * prefac

        else:
            S[row_inds, col_inds] = full_terms
            S[col_inds, row_inds] = full_terms
            prefac = S

        if return_prefactor:
            return prefac, S
        else:
            return S

    @classmethod
    def evaluate_overlap(cls, overlap_data:'OverlapGaussianData', logger=None, return_prefactor=False):
        return cls._rot_base_S(overlap_data, logger=logger, return_prefactor=return_prefactor)

class DGBKineticEnergyEvaluator(DGBEvaluator):
    """
    """

    @classmethod
    def _evaluate_polynomial_ke(cls, overlap_data:'OverlapGaussianData', terms, prefactors):
        # allows us to create a tensor of prefactors and a tensor of terms and induce
        # a polynomial by the action of the terms and then multiply that elementwise by
        # the prefactors to build the total polynomial expansion we can integrate
        # rows = overlap_data.row_indices
        cols = overlap_data.col_indices
        centers = overlap_data.centers
        # alphas = overlap_data.alphas
        init_sigs = overlap_data.initial_precision_matrices
        init_cents = overlap_data.initial_centers

        nz_pos = np.where(prefactors != 0)
        if len(nz_pos) == 0 or len(nz_pos[0]) == 0:
            # short circuit
            return np.zeros((init_cents.shape[0], init_cents.shape[0]))

        # ndim = centers.shape[-1]
        poly = cls._induced_polys(terms, init_sigs[cols], init_cents[cols], centers)
        pcs = poly.coeffs
        if not isinstance(pcs, np.ndarray):
            pcs.asarray()
        # force prefactors to be broadcastable
        scaled_coeffs = np.expand_dims(
            prefactors[nz_pos],
            list(range(1, 1 + poly.coordinate_dim))
        ) * pcs[nz_pos]
        poly = DensePolynomial(scaled_coeffs, stack_dim=1)
        # print(poly.coeffs[0])
        # poly = poly.shift(-centers[nz_pos[0]]).make_sparse_backed()
        # raise Exception(centers[0], alphas[0], init_sigs[0])
        tcoeffs = poly.coefficient_tensors
        (keys, subsel), sorting = nput.group_by(np.arange(len(nz_pos[0])), nz_pos[0])
        if len(keys) != len(centers):
            raise NotImplementedError("no contrib case not yet handled...")
        # we assume keys is sorted, so we just need to iterate over the subsels
        tensors = [
            [np.sum(x[sel], axis=0) for x in tcoeffs]
            for sel in subsel
        ]
        # nput.group_by(, nz_pos[0])
        # poly = poly.shift(-np.expand_dims(centers, list(range(1, len(terms)+1))))

        # # raise Exception(poly.coeffs.shape)
        # for pc, prefs in zip(poly.coeffs, prefactors):
        #     nz_pos = np.where(prefs != 0)
        #     if len(nz_pos) == 0:
        #         raise NotImplementedError("no contrib case not yet handled...")
        #     sub_poly = DensePolynomial(prefs[nz_pos] * pc[nz_pos], stack_dim=1)
        #     print(sub_poly.coefficient_tensors[0], sub_poly.coefficient_tensors[1])
        #     tc_sum = [np.sum(x, axis=0) for x in sub_poly.coefficient_tensors]
        #     tensors.append(tc_sum)

        # for ri, ci, pref, new_cent in zip(rows, cols, prefactors, centers): #TODO: vectorize most of this
        #     sig = init_sigs[ci]
        #     cent = init_cents[ci]
        #     polys = prefactors * cls._induced_polys(terms, ndim, sig, cent)
        #     expansion = None
        #     for p in polys.flat:
        #         p = p.shift(-new_cent)
        #         if expansion is None:
        #             expansion = p.coefficient_tensors
        #         else:
        #             expansion = [p1+p2 for p1, p2 in zip(expansion, p.coefficient_tensors)]
        #     tensors.append(expansion)

        tensors = [
            np.array([t[i] for t in tensors])
            for i in range(len(tensors[0]))
        ]
        # print(tensors[0][:3].flatten())
        # print(tensors[2][:3].flatten())

        return -1/2*cls.tensor_expansion_integrate(init_cents.shape[0], tensors, overlap_data, reweight=False)

    @classmethod
    def _induced_polys(cls, terms, sigs, initial_centers, final_centers):
        centering_shift = initial_centers - final_centers

        sig_polys = DensePolynomial.from_tensors(
            [np.zeros(sigs.shape[0]), np.zeros(sigs.shape[:2]), -sigs/2]
        ).shift(centering_shift)#.make_sparse_backed(1e-10)
        sig_grad = sig_polys.grad() # for applying p

        if 'q' in terms:
            monomial = DensePolynomial.from_tensors( # for applying q
                [np.zeros(sig_polys.coordinate_dim), np.eye(sig_polys.coordinate_dim)]
            )
            mcs = monomial.coeffs
            if not isinstance(mcs, np.ndarray):
                mcs = mcs.asarray()
            monomial = DensePolynomial(
                np.broadcast_to(mcs[np.newaxis], (sigs.shape[0],) + mcs.shape),
                stack_dim=2
            )
            monomial = monomial.shift(-final_centers[:, np.newaxis, :])#.make_sparse_backed(1e-10)

        else:
            monomial = None

        poly = DensePolynomial(
            np.ones((len(sigs),) + (1,) * sig_polys.coordinate_dim),
            stack_dim=1
        )#.make_sparse_backed(1e-10) # stack_dim = 1
        for t in terms:
            # _ = []
            # from Peeves.Timer import Timer
            # with Timer(tag=f"term {t}"):
            if t == 'q':

                mon_coeffs = monomial.coeffs
                if isinstance(mon_coeffs, np.ndarray):
                    mon_coeffs = np.expand_dims(mon_coeffs, list(range(1, poly.stack_dim)))
                else:
                    mon_coeffs = mon_coeffs.expand_dims(list(range(1, poly.stack_dim)))

                poly_coeffs= poly.coeffs
                if isinstance(poly_coeffs, np.ndarray):
                    poly_coeffs = np.expand_dims(poly_coeffs, poly.stack_dim)
                else:
                    poly_coeffs = poly_coeffs.expand_dims(poly.stack_dim)

                poly = DensePolynomial(
                    poly_coeffs,
                    stack_dim=poly.stack_dim + 1
                ) * DensePolynomial(
                    mon_coeffs,
                    stack_dim=poly.stack_dim + 1
                )
                # for i in range(ndim):
                #     monomial = DensePolynomial(np.zeros((2,) * ndim))
                #     if coeffs is None:
                #         new = monomial
                #     else: # coeffs is a NumPy array of objects
                #         new = np.apply_along_axis(
                #             lambda poly_slice:cls._apply_q(poly_slice, monomial),
                #             -1,
                #             coeffs
                #         )
                #     _.append(new)

                # print(poly.coeffs)
            elif t == 'p':
                # first term is just taking derivs of coeffs
                # basic idea is each time we hit a p we have
                # d/dq_i(K * e^(sigma*q*q)) = (d/dq_i(K) + K * ??? )e^(sigma*q*q)
                poly_1 = poly.grad()

                sigg_coeffs = sig_grad.coeffs
                if isinstance(sigg_coeffs, np.ndarray):
                    sigg_coeffs = np.expand_dims(sigg_coeffs, list(range(1, poly.stack_dim)))
                else:
                    sigg_coeffs = sigg_coeffs.expand_dims(list(range(1, poly.stack_dim)))

                poly_coeffs = poly.coeffs
                if isinstance(poly_coeffs, np.ndarray):
                    poly_coeffs = np.expand_dims(poly_coeffs, poly.stack_dim)
                else:
                    poly_coeffs = poly_coeffs.expand_dims(poly.stack_dim)

                pad_coeffs = DensePolynomial(
                    poly_coeffs,
                    stack_dim=poly.stack_dim + 1
                )
                pad_grad = DensePolynomial(
                    sigg_coeffs,
                    stack_dim=poly.stack_dim + 1
                )
                # raise Exception(pad_coeffs, pad_grad)
                poly_2 = pad_coeffs * pad_grad
                poly = poly_1 + poly_2
                # cfs = poly.coeffs
                    # raise Exception(cfs.shape)
            # for i in range(ndim):
            #     if coeffs is None:
            #         new = sig_polys.deriv(i)
            #     else: # coeffs is a NumPy array of objects
            #         new = np.apply_along_axis(
            #             lambda poly_slice:cls._apply_p(poly_slice, sig_polys, i), -1, coeffs
            #         )
            #     _.append(new)
            else:
                raise ValueError("don't understand term {t}".format(t=t))

            # cshape = () if coeffs is None else coeffs.shape
            # coeffs = np.full((ndim,) + cshape, None, dtype=object)
            # for i, p in enumerate(_):
            #     coeffs[i] = p
            # print(coeffs.shape, coeffs)

        return poly

    @abc.abstractmethod
    def evaluate_ke(self, overlap_data:'OverlapGaussianData', logger=None, **opts):
        ...

    @classmethod
    def evaluate_diagonal_rotated_momentum_contrib(self, overlap_data:'OverlapGaussianData', masses):
        # polynomial eval is too slow, need to fix in general but this is faster in the meantime
        si_cov = overlap_data.shift_matrices
        diag_inds = np.diag_indices_from(si_cov[0])
        delta_vec = overlap_data.delta_position
        shift_contrib = (delta_vec)**2
        diag_inds = (slice(None),) + diag_inds
        diag_contrib = si_cov[diag_inds]

        base_contrib = diag_contrib - shift_contrib

        j = overlap_data.initial_phases
        if j is not None:
            jp_diff = overlap_data.delta_phase_diff
            jp_sum = overlap_data.delta_phase_sum

            cos_sum = overlap_data.cos_correlation_sum
            cos_diff = overlap_data.cos_correlation_diff
            sin_sum = overlap_data.sin_correlation_sum
            sin_diff = overlap_data.sin_correlation_diff

            diff_contrib = 1 / 2 * np.dot(
                cos_diff[:, np.newaxis]*(base_contrib+jp_diff**2)
                + sin_diff[:, np.newaxis] * (2*jp_diff*delta_vec), # factored out the negative sign already
                1 / masses
            )
            sum_contrib = 1 / 2 * np.dot(
                cos_sum[:, np.newaxis] * (base_contrib+jp_sum**2)
                + sin_sum[:, np.newaxis] * (2*jp_sum*delta_vec), # factored out the negative sign already
                1 / masses
            )

            sum_prefac = overlap_data.decay_factor_sum
            diff_prefac = overlap_data.decay_factor_diff
            contrib = diff_prefac * diff_contrib + sum_prefac * sum_contrib

        else:
            contrib = 1 / 2 * np.dot(base_contrib, 1 / masses)

        n = overlap_data.npts
        rows, cols = overlap_data.indices
        ke = np.zeros((n, n))
        ke[rows, cols] = contrib
        ke[cols, rows] = contrib
        return ke

class DGBCartesianEvaluator(DGBKineticEnergyEvaluator):

    def __init__(self, masses):
        self.masses = masses

    def evaluate_ke(self, overlap_data:'OverlapGaussianData', logger=None):
        logger = Logger.lookup(logger)
        logger.log_print("evaluating Cartesian kinetic energy contribution")

        masses = self.masses
        return self.evaluate_diagonal_rotated_momentum_contrib(
            overlap_data,
            masses
        )

        # prefactors = np.broadcast_to(
        #     np.diag(1 / (masses))[np.newaxis],
        #     (overlap_data.centers.shape[0], len(masses), len(masses))
        # )
        # return self._evaluate_polynomial_ke(overlap_data, ['p', 'p'], prefactors)

class DGBWatsonEvaluator(DGBKineticEnergyEvaluator):
    """
    """
    def __init__(self, modes, coriolis_inertia_function):
        self.modes = modes
        self.ci_func = coriolis_inertia_function

    @classmethod
    def _embed_watson_modes(cls, watson_data, centers): # for putting cartesians in the watson frame?
        origin = watson_data['origin'].reshape(-1)
        tf = watson_data['matrix']
        new_cent_vecs = (centers-origin[np.newaxis]).reshape(centers.shape[0], -1, centers.shape[1])@tf.T[np.newaxis]
        return new_cent_vecs.reshape(-1, new_cent_vecs.shape[-1]) # basically a squeeze

    @staticmethod
    def annoying_coriolis_term(
            n, u, m, v,

            Xc, Dx, Sc, Sp,
            Gi, Gj, DG
    ):
        return (
                (Xc[:, m] * Xc[:, n] + Sc[:, n, m]) * (Dx[:, v] * Dx[:, u] - Sp[:, u, v])
                - (Gi[:, m, v]*Gj[:, n, u] + Gi[:, m, u]*Gj[:, n, u])
                - (Dx[:, v]*Xc[:, n]*DG[:, m, u] + Dx[:, v]*Xc[:, m]*DG[:, n, u])
        )
    @staticmethod
    def annoying_coriolis_momentum_term(
            n, u, m, v,

            Xc, r, Jp, Dx,
            Sc, Sp, DG
            # Gi, Gj, DG
    ):
        return (
            Jp[:, u]*Jp[:, v]*(Xc[:, m]*Xc[:, n] - r[:, m]*r[:, m] + Sc[:, n, m])
            -r[:, m]*r[:, n]*(Dx[:, v] * Dx[:, u] - Sp[:, u, v])
            -(Xc[:, n]*r[:, m] + Xc[:, m]*r[:, n])*(Dx[:, v]*Jp[:, u] + Dx[:, u]*Jp[:, v])
            +Jp[:, v]*r[:, n]*DG[:, m, u] + Jp[:, v]*r[:, m]*DG[:, n, u]
        )
    @staticmethod
    def annoying_imaginary_momentum_term(
            n, u, m, v,

            Xc, r, Jp, Dx,
            Sc, Sp, DG
            # Gi, Gj, DG
    ):
        return (
            -(Xc[:, n]*r[:, m] + Xc[:, m]*r[:, n])*(Dx[:, v] * Dx[:, u] - Sp[:, u, v] - Jp[:, u]*Jp[:, v])
            -(Xc[:, m]*Xc[:, n] - r[:, n]*r[:, m] + Sc[:, n, m])*(Dx[:, v]*Jp[:, u] + Dx[:, u]*Jp[:, v])
            +(Jp[:, v]*Xc[:, n]+Dx[:, v]*r[:, n])*DG[:, m, u]
            +(Jp[:, v]*Xc[:, m]+Dx[:, v]*r[:, m])*DG[:, n, u]
        )
    @classmethod
    def evaluate_coriolis_contrib(cls, coriolis_tensors, overlap_data:'OverlapGaussianData'):

        Sc = overlap_data.covariance_matrices
        Sp = overlap_data.shift_matrices
        SI = overlap_data.initial_precision_matrices
        Xc = overlap_data.centers

        rows, cols = overlap_data.indices
        SIi = SI[rows]
        SIj = SI[cols]

        Gi = Sc @ SIi
        Gj = Sc @ SIj
        DG = Gi - Gj

        Dx = overlap_data.center_difference
        Dx = np.reshape(Sp @ overlap_data.center_difference[:, :, np.newaxis], Dx.shape)

        ndim = Xc.shape[-1]

        term = lambda n, u, m, v: coriolis_tensors[:, n, u, m, v] * cls.annoying_coriolis_term(
            n, u, m, v,

            Xc, Dx, Sc, Sp,
            Gi, Gj, DG
        )

        contrib = np.zeros(Sc.shape[0])
        inds = itertools.product(*[range(ndim)]*4)
        for n,u,m,v in inds:
            if n != u and m != v:
                t = term(n, u, m, v)
                contrib += t

        if overlap_data.initial_momenta is not None:
            jp_sum = overlap_data.delta_phase_sum
            jp_diff = overlap_data.delta_phase_diff
            rho_sum = overlap_data.rho_sum
            rho_diff = overlap_data.rho_diff

            sum_real_term = lambda n, u, m, v: coriolis_tensors[:, n, u, m, v] * cls.annoying_coriolis_momentum_term(
                    n, u, m, v,

                    Xc, rho_sum, jp_sum, Dx,
                    Sc, Sp, DG
            )
            sum_imag_term = lambda n, u, m, v: coriolis_tensors[:, n, u, m, v] * cls.annoying_imaginary_momentum_term(
                n, u, m, v,

                Xc, rho_sum, jp_sum, Dx,
                Sc, Sp, DG
            )

            diff_real_term = lambda n, u, m, v: coriolis_tensors[:, n, u, m, v] * cls.annoying_coriolis_momentum_term(
                n, u, m, v,

                Xc, rho_diff, jp_diff, Dx,
                Sc, Sp, DG
            )
            diff_imag_term = lambda n, u, m, v: coriolis_tensors[:, n, u, m, v] * cls.annoying_imaginary_momentum_term(
                n, u, m, v,

                Xc, rho_diff, jp_diff, Dx,
                Sc, Sp, DG
            )

            sum_real_contrib = np.zeros(Sc.shape[0])
            sum_imag_contrib = np.zeros(Sc.shape[0])
            diff_real_contrib = np.zeros(Sc.shape[0])
            diff_imag_contrib = np.zeros(Sc.shape[0])
            inds = itertools.product(*[range(ndim)] * 4)
            for n, u, m, v in inds:
                if n != u and m != v:
                    real_cont = sum_real_term(n, u, m, v)
                    imag_cont = sum_imag_term(n, u, m, v)
                    # print("Sum I:", (n, u, m, v), 219475*imag_cont)
                    # print("Sum R:", (n, u, m, v), 219475*real_cont)
                    sum_real_contrib += real_cont
                    sum_imag_contrib += imag_cont
                    real_cont = diff_real_term(n, u, m, v)
                    imag_cont = diff_imag_term(n, u, m, v)
                    # print("Diff I:", (n, u, m, v), 219475*imag_cont)
                    # print("Diff R:", (n, u, m, v), 219475*real_cont)
                    diff_real_contrib += real_cont
                    diff_imag_contrib += imag_cont

            cos_sum, sin_sum = overlap_data.correlation_factors_sum
            cos_diff, sin_diff = overlap_data.correlation_factors_diff
            sum_prefac = overlap_data.decay_factor_sum
            diff_prefac = overlap_data.decay_factor_diff

            diff_contrib = (
                    cos_diff * (contrib + diff_real_contrib)
                    + sin_diff * diff_imag_contrib
            )
            sum_contrib = (
                    cos_sum * (contrib + sum_real_contrib)
                    + sin_sum * sum_imag_contrib
            )

            contrib = diff_prefac * diff_contrib + sum_prefac * sum_contrib

        npts = overlap_data.npts
        ke = np.zeros((npts, npts))
        ke[rows, cols] = -contrib
        ke[cols, rows] = -contrib

        return ke

    @classmethod
    def evaluate_watson_term(cls, B_e, overlap_data:'OverlapGaussianData'):
        B_e = -1 / 4 * B_e
        if overlap_data.initial_momenta is not None:
            phase_sum = overlap_data.cos_correlation_sum
            phase_diff = overlap_data.cos_correlation_diff

            sum_prefac = overlap_data.decay_factor_sum
            diff_prefac = overlap_data.decay_factor_diff

            diff_wat = np.cos(phase_diff)
            sum_wat = np.cos(phase_sum)
            watson_contrib = B_e * (diff_prefac * diff_wat + sum_prefac * sum_wat)
        else:
            watson_contrib = B_e

        n = overlap_data.npts
        watson = np.zeros((n, n))
        rows, cols = overlap_data.indices
        watson[rows, cols] = watson_contrib
        watson[cols, rows] = watson_contrib

        return watson

    def evaluate_ke(self, overlap_data: 'OverlapGaussianData',
                    logger=None,
                    include_diagonal_contribution=True,
                    include_coriolis_coupling=True,
                    include_watson_term=True
                    ):

        logger = Logger.lookup(logger)
        logger.log_print("evaluating Watson diagonal momentum contribution")
        if include_diagonal_contribution:
            base = self.evaluate_diagonal_rotated_momentum_contrib(
                overlap_data,
                np.ones(overlap_data.ndim)
            )
        else:
            base = 0

        if include_watson_term or include_coriolis_coupling:
            if include_coriolis_coupling:
                logger.log_print("evaluating coriolis contribution")
            B_e, coriolis = self.ci_func(overlap_data.centers)
            if include_coriolis_coupling:
                coriolis = self.evaluate_coriolis_contrib(coriolis, overlap_data)
            else:
                coriolis = 0
        else:
            coriolis = 0
            B_e = 0

        if include_watson_term:
            logger.log_print("evaluating Watson term contribution")
        if include_watson_term:
            watson = self.evaluate_watson_term(B_e, overlap_data)
        else:
            watson = 0

        ke = base + watson + coriolis

        return ke

class DGBPotentialEnergyEvaluator(DGBEvaluator):
    """
    An evaluator designed
    """

    def __init__(self,
                 potential_function,
                 integral_handler=None,
                 expansion_degree=None,
                 expansion_type=None,
                 quadrature_degree=None,
                 pairwise_functions=None,
                 logger=None
                 ):
        self.potential_function = potential_function
        self.handler = integral_handler
        self.expansion_degree = expansion_degree
        self.expansion_type = expansion_type
        self.quadrature_degree = quadrature_degree
        self.pairwise_handler = pairwise_functions

        self.logger = Logger.lookup(logger)

    def analytic_integrate(self):
        raise NotImplementedError("flooped up")
        centers = [np.array(np.meshgrid(x, x)).T for x in self.centers.T]
        alphas = np.array(np.meshgrid(self.alphas, self.alphas)).T
        # raise Exception(alphas.shape)
        return self.potential_function['analytic_integrals'](
            centers, # ...no
            alphas
        )

    @classmethod
    def expansion_integrate(cls,
                            function,
                            overlap_data:'OverlapGaussianData',
                            expansion_type,
                            expansion_degree=2,
                            pairwise_functions=None,
                            logger=None
                            ):

        if pairwise_functions is not None:
            logger.log_print("evaluating pairwise potential contributions...")
            pot_contribs, deriv_corrs = pairwise_functions.evaluate_pairwise_contrib(
                overlap_data,
                expansion_degree=None if expansion_degree < 0 else expansion_degree
            )
        else:
            pot_contribs = None
            deriv_corrs = None

        if not isinstance(logger, Logger):
            logger = Logger.lookup(logger)

        if expansion_degree > -1:

            deriv_order = expansion_degree # renamed for consistency
            if expansion_type == 'taylor':
                raise NotImplementedError('need to fix up taylor expansion')
                logger.log_print("expanding as a Taylor series about the minumum energy geometry...")
                assert self.ref is None #TODO: centers need a displacement
                zero = np.zeros((1, centers.shape[-1])) if self.ref is None else np.array([self.ref])
                derivs = function(zero, deriv_order=deriv_order)
                if isinstance(derivs, np.ndarray): # didn't get the full list so we do the less efficient route
                    derivs = (
                            [function(zero)] +
                            [function(zero, deriv_order=d) for d in range(1, deriv_order + 1)]
                    )
                derivs = [
                    np.broadcast_to(
                        d[np.newaxis],
                        alphas.shape + d.squeeze().shape
                    )
                    for d in derivs
                ]
            else:
                alphas = overlap_data.alphas
                centers = overlap_data.centers
                if overlap_data.momenta_sum is None:
                    deriv_order = deriv_order - (deriv_order % 2) # odd orders don't contribute so why evaluate the derivatives...
                # if self.mass_weighted:
                #     derivs = self.mass_weighted_eval(function, centers, self.masses, deriv_order=deriv_order)
                # else:
                derivs = function(centers, deriv_order=deriv_order)
                if isinstance(derivs, np.ndarray):  # didn't get the full list so we do the less efficient route'
                    derivs = [function(centers)] + [
                        function(centers, deriv_order=d) for d in range(1, deriv_order+1)
                    ]

            if deriv_corrs is not None:
                _ = []
                for i,(d1,d2) in enumerate(zip(derivs, deriv_corrs)):
                    _.append(d1 - d2)
                derivs = _

            npts = overlap_data.npts
            pot = cls.tensor_expansion_integrate(
                npts,
                derivs,
                overlap_data,
                reweight=True,
                expansion_type=expansion_type,
                logger=logger
            )

        else:
            npts = overlap_data.npts
            pot = np.zeros((npts, npts))

        # if prefactors is not None:
        #     row_inds, col_inds = np.triu_indices(npts)
        #     pref = np.zeros((npts, npts))
        #     pref[row_inds, col_inds] = prefactors
        #     pref[col_inds, row_inds] = prefactors
        #     pot *= pref

        if pot_contribs is not None: # implies no rotation
            row_inds, col_inds = np.triu_indices(npts)
            pot2 = np.zeros((npts, npts))
            pot2[row_inds, col_inds] = pot_contribs
            pot2[col_inds, row_inds] = pot_contribs
            pot += pot2

        return pot

    @classmethod
    def evaluate_multiplicative(
            cls,
            function,
            overlap_data:'OverlapGaussianData',
            integral_handler=None,
            expansion_degree=None,
            expansion_type=None,
            quadrature_degree=None,
            pairwise_functions=None,
            logger=None
    ):

        if integral_handler is None:
            if isinstance(function, dict):
                if 'analytic_integrals' in function:
                    integral_handler = 'analytic'
            elif expansion_degree is not None:
                integral_handler = 'expansion'
            else:
                integral_handler = 'quad'

        if not isinstance(logger, Logger):
            logger = Logger.lookup(logger)


        if integral_handler == 'quad':
            if pairwise_functions is not None:
                raise ValueError("pairwise functions can't go with direct quadrature")
            logger.log_print("evaluating integrals with {n}-order quadrature about {npt} points",
                             n=quadrature_degree,
                             npt=overlap_data.npts
                             )
            pot_mat = cls.quad_integrate(
                function,
                overlap_data,
                degree=quadrature_degree,
                logger=logger
            )
        elif integral_handler == 'expansion':
            logger.log_print("evaluating integrals with {n}-degree expansions about {npt} points",
                             n=expansion_degree,
                             npt=overlap_data.npts
                             )
            pot_mat = cls.expansion_integrate(function,
                                              overlap_data,
                                              expansion_degree=expansion_degree,
                                              expansion_type=expansion_type,
                                              pairwise_functions=pairwise_functions,
                                              logger=logger
                                              )
        elif integral_handler == 'analytic':
            logger.log_print("evaluating integrals analytically", n=expansion_degree)
            pot_mat = cls.analytic_integrate()
        else:
            raise ValueError("unknown operator evaluation scheme {}".format(integral_handler))

        return pot_mat

    def evaluate_pe(self, overlap_data:'OverlapGaussianData', logger=None):
        if logger is None:
            logger = self.logger
        logger.log_print("evaluating potential energy contribution")
        return self.evaluate_multiplicative(
            self.potential_function,
            overlap_data,
            integral_handler=self.handler,
            quadrature_degree=self.quadrature_degree,
            expansion_degree=self.expansion_degree,
            expansion_type=self.expansion_type,
            pairwise_functions=self.pairwise_handler,
            logger=logger
        )

    def evaluate_op(self,
                    operator,
                    overlap_data:'OverlapGaussianData',
                    integral_handler=None,
                    expansion_degree=None,
                    expansion_type=None,
                    quadrature_degree=None,
                    pairwise_functions=None,
                    logger=None
                    ):
        no_type_passed = expansion_degree is None and quadrature_degree is None
        return self.evaluate_multiplicative(
            operator,
            overlap_data,
            integral_handler=self.handler if integral_handler is None else integral_handler,
            quadrature_degree=self.quadrature_degree if no_type_passed else quadrature_degree,
            expansion_degree=self.expansion_degree if no_type_passed else expansion_degree,
            expansion_type=self.expansion_type if expansion_type is None else expansion_type,
            pairwise_functions=pairwise_functions,
            logger=self.logger if logger is None else logger
        )

class DGBPairwisePotentialEvaluator(DGBEvaluator, metaclass=abc.ABCMeta):
    def __init__(self, coords, pairwise_potential_functions, quadrature_degree=3,
                 use_with_interpolation='ignored'
                 ):
        self.coords = coords
        self.pairwise_potential_functions = pairwise_potential_functions
        self.quadrature_degree = quadrature_degree

    @classmethod
    def get_bond_length_deltas(cls, natoms, ndim, i, j, full=False):
        if not full:
            mat = np.zeros((ndim, natoms, ndim))
            for k in range(ndim):
                mat[k][i][k] = 1
                mat[k][j][k] = -1
            mat = mat.reshape((mat.shape[0], natoms * ndim))
        else:
            mat = np.sqrt(2) * np.eye(natoms*ndim).reshape((natoms, ndim, natoms, ndim))
            for k in range(ndim):
                mat[i][k][i][k] = 1
                mat[i][k][j][k] = -1
                mat[j][k][i][k] = 1
                mat[j][k][j][k] = 1
            mat = mat.reshape((natoms * ndim, natoms * ndim))
        return mat

    @abc.abstractmethod
    def get_coordinate_bond_length_projection(self, i, j) -> 'tuple[np.ndarray, np.ndarray]':
        ...
    def get_coordinate_change_transformation(self, coordinate_projection_data) -> np.ndarray:
        v, s, w = np.linalg.svd(coordinate_projection_data)  # w is the important transformation
        # we get the dimension for every specified transformation (doesn't always need to be the same)
        if not np.allclose(np.abs(s), s):
            s_sorts = np.argsort(-np.abs(s), axis=-1)
            s = np.abs(s[s_sorts])
            w = w[s_sorts]
        w = w.transpose(0, 2, 1)
        dims = np.sum(np.abs(s) > 1e-14, axis=1)
        udims, upos = np.unique(dims, return_index=True)
        if len(udims) == 1:
            return w[:, :, :udims[0]]
        else:
            raise NotImplementedError("bond length transforms with different dimensions not yet supported...")

    def get_bond_length_change_transformation(self, overlap_data:'OverlapGaussianData', i, j) -> np.ndarray:
        d0, base_proj = self.get_coordinate_bond_length_projection(i, j)
        proj_data = np.broadcast_to(
            base_proj[np.newaxis],
            (len(overlap_data.rotations),) + base_proj.shape
        ) #@ overlap_data.rotations.transpose(0, 2, 1) # not sure if we need a .T or not...
        return self.get_coordinate_change_transformation(proj_data)

    def wrap_distance_function(self, i, j, overlap_data:'OverlapGaussianData', transformations, pairwise_function):
        d0, proj = self.get_coordinate_bond_length_projection(i, j) # dim x ncoordx
        tf_proj = proj[np.newaxis] @ transformations # N x dim x subdim
        @functools.wraps(pairwise_function)
        def eval_pairwise(coords, *, deriv_order=None, d0=d0, tf_proj=tf_proj):
            # coords is expressed in terms of some sort of transformed version
            # of the normal modes (or Cartesians) that were supplied to the program
            # and we know exactly how these change the deltas that feed into the bond
            # length using tf_proj from above, this means we can get the displacements of the deltas
            # for free and ... add these onto the initial deltas???
            base_shape = coords.shape[:-1]
            coords = coords.reshape(tf_proj.shape[0], -1, coords.shape[-1]) # N x k x subdim
            proj_d = tf_proj[:, np.newaxis, :, :] @ coords[:, :, :, np.newaxis] # N x k x dim x 1
            d = d0[np.newaxis, np.newaxis, :] + proj_d.reshape(proj_d.shape[:3])

            if d.shape[1] == 1:
                d = d.reshape((d.shape[0], d.shape[2]))

            if deriv_order is None or deriv_order == 0:
                r = np.linalg.norm(d, axis=-1)  # we ensure no scaling of difference coords
                r_derivs = None
            else:
                r_derivs = nput.vec_norm_derivs(d, order=deriv_order)
                r, r_derivs = r_derivs[0], r_derivs[1:]

            fdat = pairwise_function(r, deriv_order=deriv_order)
            if r_derivs is not None and len(r_derivs) > 0:
                # first we get df/d_delta by transforming df/dr.dr/d_delta
                f_vals = fdat[0]
                r_derivs = [r[..., np.newaxis] for r in r_derivs]
                f_derivs = [fdat[i].reshape(f_vals.shape + (1,) * i) for i in range(1, deriv_order + 1)]
                f_derivs = list(TensorDerivativeConverter(r_derivs, f_derivs).convert())
                # why bother? -> we're just gonna remove this part of the transformation anyway
                # then we get df/dq by transforming with the linear transformation d_del/dq
                q_derivs = []
                for i,d in enumerate(f_derivs):
                    for j in range(d.ndim-1):
                        d = np.tensordot(d, proj, axes=[1, 0])
                        # d = nput.vec_tensordot(d, tf_proj, shared=1, axes=[1, 1])
                    q_derivs.append(d)
                fdat = [f_vals] + q_derivs

            if deriv_order is not None:
                # print(base_shape, [f.shape for f in fdat])
                fdat = [f.reshape(base_shape + f.shape[r.ndim:]) for f in fdat]
            else:
                fdat = fdat.reshape(base_shape + fdat.shape[r.ndim:])
            return fdat

        return eval_pairwise

    # def get_angle_change_transformation(self, coords, i, j, k):
    #     ...

    # def wrap_angle_function(self, coord, i, j, k, transformations, pairwise_function):
    #     ...

    def evaluate_pairwise_contrib(self,
                                  overlap_data:'OverlapGaussianData',
                                  quadrature_degree=None,
                                  expansion_degree=2
                                  ):

        if quadrature_degree is None:
            quadrature_degree = self.quadrature_degree

        centers = overlap_data.centers
        covs = overlap_data.covariance_matrices
        mom_sum = overlap_data.momenta_sum
        mom_diff = overlap_data.momenta_diff

        npts = centers.shape[0]
        fdim = centers.shape[1]

        potential_contrib = 0
        if expansion_degree is not None:
            derivative_contribs = [
                np.zeros((npts,) + (fdim,) * i)
                    if i > 0 else
                np.zeros((npts,))
                for i in range(expansion_degree + 1)
            ]
        else:
            derivative_contribs = [np.zeros((npts,))]

        for index_pair, pairwise_func in self.pairwise_potential_functions.items():
            if len(index_pair) == 2:
                i,j = index_pair
            else:
                raise NotImplementedError("currently only pairwise distances potentials are supported")

            tfsT = self.get_bond_length_change_transformation(overlap_data, i, j)
            tfs = tfsT.transpose((0, 2, 1)) # tfs is unitary so this is the inverse
            f = self.wrap_distance_function(i, j, overlap_data, tfsT, pairwise_func)
            sub_cov = tfs @ covs @ tfsT
            tf_cov = np.linalg.inv(sub_cov)

            tf_alphas, tf_vecs = np.linalg.eigh(tf_cov)
            tf_alphas = tf_alphas / 2 # baked into the covariances...
            ndim = tf_alphas.shape[1]

            tf_centers = tfs @ centers[:, :, np.newaxis]
            tf_centers = tf_centers.reshape(tf_centers.shape[:2])

            if mom_sum is not None:
                tf_momenta = [
                    np.reshape(tf_cov @ tfs @ covs @ mom_diff[:, :, np.newaxis], sub_cov.shape[:2]),
                    np.reshape(tf_cov @ tfs @ covs @ mom_sum[:, :, np.newaxis], sub_cov.shape[:2]),
                ]
            else:
                tf_momenta = None

            pairwise_contrib = self.rotated_gaussian_quadrature(
                f,
                tf_alphas,
                tf_centers,
                tf_vecs,
                tf_vecs.transpose(0, 2, 1),
                tf_momenta,
                degree=quadrature_degree,
                normalize=False
            )

            potential_contrib += pairwise_contrib

            if expansion_degree is not None:
                # we make use of the inverse of the transformation that takes us from
                # the space of total coordinates (normal modes or Cartesians) to the space
                # of coordinates (deltas) that lead to a change in the coordinate of interest
                # we'll note that our wrapped function is expected to implement the appropriate
                # conversion from derivatives with respect to r to derivatives with respect to
                # the deltas

                tf_derivs_contribs = f(tf_centers, deriv_order=expansion_degree)
                for i, d in enumerate(tf_derivs_contribs):
                    if i == 0:
                        derivative_contribs[0] += d
                    else:
                        # if i == 2:
                        #     print("pre", d[0])
                        # for j in range(i):
                        #     d = nput.vec_tensordot(d, tfsT, shared=1, axes=[1, 2])
                        # if i == 2:
                        #     print("post", d[0])
                        derivative_contribs[i] += d
            else:
                derivative_contribs[0] += f(tf_centers[:, :ndim])

        return potential_contrib, derivative_contribs

class DGBCartesianPairwiseEvaluator(DGBPairwisePotentialEvaluator):

    def __init__(self, coords:'DGBCartesians', pairwise_functions, **opts):
        super().__init__(coords, pairwise_functions, **opts)

    def get_coordinate_bond_length_projection(self, i, j):
        natoms, ndim = self.coords.cart_shape[1:]
        base_mat = self.get_bond_length_deltas(natoms, ndim, i, j)
        # inverse_projection = base_mat.T / 2
        d0 = np.zeros(ndim)
        return d0, base_mat

class DGBWatsonPairwiseEvaluator(DGBPairwisePotentialEvaluator):
    def __init__(self, coords:'DGBWatsonModes', pairwise_functions, **opts):
        super().__init__(coords, pairwise_functions, **opts)

    def get_coordinate_bond_length_projection(self, i, j):
        natoms = len(self.coords.masses)
        modes = self.coords.modes.matrix
        ndim = modes.shape[0] // natoms # this will almost always be 3
        base_mat = self.get_bond_length_deltas(natoms, ndim, i, j)
        tf_base = base_mat @ modes
        d0 = np.dot(base_mat, self.coords.modes.origin.reshape(-1))
        return d0, tf_base

