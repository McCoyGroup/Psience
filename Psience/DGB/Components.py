
import abc
import math

import numpy as np, scipy as sp, itertools, functools

from McUtils.Zachary import RBFDInterpolator, DensePolynomial, TensorDerivativeConverter
from McUtils.Coordinerds import CartesianCoordinates1D, CartesianCoordinates2D, CartesianCoordinates3D
from McUtils.Scaffolding import Logger
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput

from ..Molecools import StructuralProperties
from ..MixtureModes import NormalModes

__all__ = [
    "DGBEvaluator",
    "DGBKineticEnergyEvaluator",
    "DGBCartesianEvaluator",
    "DGBWatsonEvaluator",
    "DGBPotentialEnergyEvaluator",
    "DGBPairwisePotentialEvaluator",
    "DGBCoords",
    "DGBCartesians",
    "DGBInternals",
    "DGBWatsonModes",
    "DGBGaussians"
]

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
    def get_overlap_gaussians(cls, centers, alphas, transformations):
        # a bit inefficient if all tfs are identity, but generality is good
        rows, cols = np.triu_indices(len(alphas))

        sigs = cls.get_inverse_covariances(alphas, transformations)
        new_sigs = sigs[rows] + sigs[cols]

        new_alphas, new_rots = np.linalg.eigh(new_sigs)  # eigenvalues of inverse tensor...
        new_rots_inv = new_rots.transpose(0, 2, 1)

        # I _could_ construct the inverse from the alphas and rotations
        # but I think it makes more sense to use a potentially more stable
        # inverse here...
        new_inv = np.linalg.inv(new_sigs)

        new_centers = new_inv @ (
                sigs[rows] @ centers[rows][:, :, np.newaxis]
                + sigs[cols] @ centers[cols][:, :, np.newaxis]
        )
        new_centers = new_centers.reshape(centers[cols].shape)
        new_alphas = new_alphas / 2
        sum_sigs = sigs[rows] @ new_inv @ sigs[cols]

        # TODO: turn this into a proper object...
        overlap_data = {
            'centers': new_centers,
            'alphas': new_alphas,
            'covariances': new_sigs,
            'inverse':new_inv,
            'sum_inverse': sum_sigs,
            'row_inds': rows,
            'col_inds': cols,
            'rotations': new_rots,
            'inverse_rotations': new_rots_inv,
            'init_centers':centers,
            'init_alphas':alphas,
            'init_covariances':sigs
        }

        return overlap_data

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
    def simple_poly_int(cls, n):
        return np.prod(np.arange(1, n, 2)) / 2**(n/2) # double factorial/gamma/whatever
    @classmethod
    def tensor_expansion_integrate(cls,
                                   npts, derivs, overlap_data,
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

        centers = overlap_data['centers']
        ndim = centers.shape[-1]
        rotations = overlap_data['rotations']
        if expansion_type == 'taylor':
            centers = np.reshape(rotations @ centers[:, :, np.newaxis], centers.shape)
        alphas = overlap_data['alphas']

        # if self.
        fdim = derivs[0].ndim - 1 # to simplify dipole functions
        fshape = derivs[0].shape[1:]

        new_derivs = []
        # rotations = rotations[:, :, :, np.newaxis] # to test shapes
        for n,d in enumerate(derivs):
            # print("???", d.shape, rotations.shape)
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

        row_inds = overlap_data['row_inds']
        col_inds = overlap_data['col_inds']

        pot = np.zeros((npts, npts) + fshape)
        caches = [{} for _ in range(ndim)]
        for nd,d in enumerate(derivs): # add up all independent integral contribs...
            # iterate over upper triangle coordinates (we'll add bottom contrib by symmetry)
            inds = itertools.combinations_with_replacement(range(ndim), r=nd) if nd > 0 else [()]
            for idx in inds:
                count_map = {k: v for k, v in zip(*np.unique(idx, return_counts=True))}
                if expansion_type != 'taylor' and any(n%2 !=0 for n in count_map.values()):
                    continue # odd contribs vanish

                contrib = 1
                for k in range(ndim): # do each dimension of integral independently
                    n = count_map.get(k, 0)
                    if n not in caches[k]:
                        if expansion_type == 'taylor':
                            raise NotImplementedError("Taylor series in rotated basis not implemented yet")
                        caches[k][n] = (
                           cls.simple_poly_int(n)
                             if expansion_type != 'taylor' else
                           cls.polyint_1D(centers[..., k], alphas[..., k], n)
                        )
                    base_contrib = caches[k][n] / alphas[..., k]**(n/2)
                    for _ in range(fdim):
                        base_contrib = np.expand_dims(base_contrib, -1)
                    if isinstance(contrib, int) and fdim > 0:
                        base_contrib = np.broadcast_to(base_contrib,  alphas.shape[:-1] + fshape)
                    contrib *= base_contrib

                dcont = d[(slice(None, None, None),)*(fdim+1) + idx] if len(idx) > 0 else d

                if reweight:
                    # compute multinomial coefficient for weighting purposes
                    _, counts = np.unique(idx, return_counts=True)
                    multicoeff = 1
                    for x in counts: multicoeff*=math.factorial(x)
                    multicoeff = multicoeff / math.factorial(len(idx))
                    scaling = multicoeff / math.factorial(nd)
                    contrib *= dcont * scaling
                else:
                    contrib *= dcont

                pot[row_inds, col_inds] += contrib

        pot[col_inds, row_inds] = pot[row_inds, col_inds]

        return pot

    @staticmethod
    def quad_nd(centers, alphas, function, flatten=False, degree=3, chunk_size=int(1e6), normalize=True):
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
            d = d_chunk[np.newaxis] / squa[:, np.newaxis]
            c = centers[:, np.newaxis, :] + d
            fv = function(c.reshape(-1, ndim))
            fshape = fv.shape[1:]
            fv = fv.reshape(c.shape[:2] + fshape)
            w_chunk = np.expand_dims(w_chunk, [0] + [-x for x in range(1, len(fshape)+1)])
            chunk_val = np.sum(w_chunk * fv, axis=1)
            if val is None:
                val = chunk_val
            else:
                val += chunk_val

        if normalize:
            normalization = 1 / np.prod(np.sqrt(alphas), axis=-1)
            val = val * normalization
        return val

    @classmethod
    def _wrap_rotated_function(cls, func, rotations, inverse):

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
            return fdat

        return eval

    @classmethod
    def rotated_gaussian_quadrature(cls,
                                   function,
                                   alphas,
                                   centers,
                                   rotations,
                                   inverse,
                                   normalize=True,
                                   degree=2,
                                   ):
        f = cls._wrap_rotated_function(function, rotations, inverse)
        centers = np.reshape(rotations @ centers[:, :, np.newaxis], centers.shape)

        vals = cls.quad_nd(centers, alphas, f, degree=degree, normalize=normalize)

        normalization = 1 / (np.sqrt(np.pi)) ** centers.shape[-1]  # part of the prefactor...
        vals *= normalization

        return vals

    @classmethod
    def quad_integrate(cls, function, overlap_data, degree=2, logger=None):
        """
        Integrate potential over all pairs of Gaussians at once

        :param degree:
        :type degree:
        :return:
        :rtype:
        """

        vals = cls.rotated_gaussian_quadrature(
            function,
            overlap_data['alphas'],
            overlap_data['centers'],
            overlap_data['rotations'],
            overlap_data['inverse_rotations'],
            degree=degree,
            normalize=False
        )

        npts = overlap_data['init_centers'].shape[0]
        rows, cols = np.triu_indices(npts)

        pots = np.zeros((npts, npts) + vals.shape[1:])
        pots[rows, cols] = vals
        pots[cols, rows] = vals

        return pots

    @classmethod
    def _rot_base_S_components(cls, rot_data):
        row_inds = rot_data['row_inds']
        col_inds = rot_data['col_inds']
        alphas = rot_data['init_alphas']
        centers = rot_data['init_centers']
        rotas = rot_data["alphas"]

        ndim = centers.shape[-1]
        dets = np.prod(rotas, axis=-1)
        # we prefactor out the 2**ndim
        rdets = np.prod(alphas[row_inds], axis=-1)
        cdets = np.prod(alphas[col_inds], axis=-1)  # literally a product of passed in alphas

        L = rot_data['inverse_rotations']  #TODO: make this make more sense...
        Lt = rot_data['rotations']
        disps = centers[row_inds] - centers[col_inds]
        C = disps[:, np.newaxis, :] @ rot_data['sum_inverse'] @ disps[:, :, np.newaxis]
        C = C.reshape(disps.shape[0])

        return alphas, rotas, L, Lt, row_inds, col_inds, ndim, rdets, cdets, dets, C

    @classmethod
    def _rot_base_S(cls, overlap_data):

        S = np.eye(len(overlap_data['init_centers']))
        alphas, rotas, L, Lt, row_inds, col_inds, ndim, rdets, cdets, dets, C = \
            cls._rot_base_S_components(overlap_data)

        S[row_inds, col_inds] = S[col_inds, row_inds] = (
                2 ** (ndim / 2) * ((rdets * cdets) / (dets ** 2)) ** (1 / 4) * np.exp(-C / 2)
        )
        return S

    @classmethod
    def evaluate_overlap(cls, overlap_data):
        return cls._rot_base_S(overlap_data) # maybe I should make this return the overlap params instead?

class DGBKineticEnergyEvaluator(DGBEvaluator):
    """
    """

    @classmethod
    def _evaluate_polynomial_ke(cls, overlap_data, terms, prefactors):
        # allows us to create a tensor of prefactors and a tensor of terms and induce
        # a polynomial by the action of the terms and then multiply that elementwise by
        # the prefactors to build the total polynomial expansion we can integrate
        # rows = overlap_data['row_inds']
        cols = overlap_data['col_inds']
        centers = overlap_data['centers']
        # alphas = overlap_data['alphas']
        init_sigs = overlap_data['init_covariances']
        init_cents = overlap_data['init_centers']

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

    # @staticmethod
    # def _apply_p(poly_slice, sig_poly, i):
    #     new = np.full(len(poly_slice), None)
    #     for n, p in enumerate(poly_slice):
    #         new[n] = p.deriv(i) + p * sig_poly.deriv(i)
    #     return new

    # @staticmethod
    # def _apply_q(poly_slice, monomial):
    #     new = np.full(len(poly_slice), None)
    #     for n,p in enumerate(poly_slice):
    #         new[n] = monomial*p
    #     return new

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
    def evaluate_ke(self, overlap_data):
        ...

    @classmethod
    def evaluate_diagonal_rotated_momentum_contrib(self, overlap_data, masses):
        # polynomial eval is too slow, need to fix in general but this is faster in the meantime
        init_covs = overlap_data['init_covariances']
        invs = overlap_data['inverse']
        si_cov = overlap_data['sum_inverse']
        init_cents = overlap_data['init_centers']

        rows, cols = np.triu_indices(len(init_covs))
        diag_inds = np.diag_indices_from(si_cov[0])
        disps = init_cents[cols] - init_cents[rows]
        shift_mat = (init_covs[rows] @ invs @ init_covs[cols])
        shift_contrib = (shift_mat @ disps[:, :, np.newaxis])**2
        shift_contrib = shift_contrib.reshape(shift_contrib.shape[:2])
        diag_inds = (slice(None),) + diag_inds
        diag_contrib = si_cov[diag_inds]

        # aa = overlap_data['init_alphas']
        # raise Exception(
        #     shift_mat[1], si_cov[1],
        #     2*(aa[0]*aa[1]/(aa[0] + aa[1]))
        # )

        base_contrib = diag_contrib - shift_contrib

        full_contrib = 1/2 * np.dot(base_contrib, 1/masses)

        ke = np.zeros((len(init_covs), len(init_covs)))
        ke[rows, cols] = full_contrib
        ke[cols, rows] = full_contrib

        return ke

    @classmethod
    def evaluate_classic_momentum_contrib(cls, overlap_data, masses):
        # we assume covariances are 1

        centers = overlap_data['init_centers'].view(np.ndarray)
        alphas = overlap_data['init_alphas']

        aouter = alphas[:, np.newaxis] * alphas[np.newaxis, :]
        aplus = alphas[:, np.newaxis] + alphas[np.newaxis, :]
        arat = aouter / aplus
        disps = centers[:, np.newaxis, :] - centers[np.newaxis, :, :]

        C = arat * np.power(disps, 2)

        # S_dim = (np.sqrt(2) * np.power(aouter, 1 / 4) / B) * np.exp(-C)

        # Really this is -1/2 (2*arat) ((2*arat)(disps^2) - 1) / m
        # but we factor out the initial -1/2
        T_dim = arat * (1 - 2*C) / masses[np.newaxis, np.newaxis, :]

        # S = np.prod(S_dim, axis=-1)
        T = np.sum(T_dim, axis=-1)

        return T

class DGBCartesianEvaluator(DGBKineticEnergyEvaluator):

    def __init__(self, masses):
        self.masses = masses

    def evaluate_ke(self, overlap_data):
        masses = self.masses

        """
        [357.41451351 343.01569053 303.52678279 248.33184886 188.6011872  133.58986444]
        [357.41451351 345.70911849 313.10549681 266.12068585 212.94433354 161.03916799]
        """

        # return self.evaluate_classic_momentum_contrib(
        #     overlap_data,
        #     masses
        # )

        return self.evaluate_diagonal_rotated_momentum_contrib(
            overlap_data,
            masses
        )

        prefactors = np.broadcast_to(
            np.diag(1 / (masses))[np.newaxis],
            (overlap_data['centers'].shape[0], len(masses), len(masses))
        )
        return self._evaluate_polynomial_ke(overlap_data, ['p', 'p'], prefactors)

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

            Sc,
            #Xc, Sp,
            Dx,
            # Gi, Gj,
            # DS,

            ScXc, DxSp, GjGi,
            DSX
    ):

        return (
            ScXc[:, n, m]*DxSp[:, u, v]
            - (GjGi[:, n, v, m, u] + GjGi[:, n, v, m, u])
            - np.reshape(
                Sc[:, :, n][:, np.newaxis, :] @ (DSX[:, :, v, u] + DSX[:, :, u, v])[:, :, np.newaxis] , (-1)
                ) * Dx[:, m]
        )
    @classmethod
    def evaluate_coriolis_contrib(cls, coriolis_tensors, overlap_data):

        Sc = overlap_data['inverse']
        Sp = overlap_data['sum_inverse']
        SI = overlap_data['init_covariances']
        X0 = overlap_data['init_centers']
        Xc = overlap_data['centers']

        rows = overlap_data['row_inds']
        cols = overlap_data['col_inds']
        SIi = SI[rows]
        SIj = SI[cols]
        Xi = X0[rows]
        Xj = X0[cols]

        Gi = Sc @ SIi
        Gj = Sc @ SIj

        DS = SIi - SIj
        Dx = np.reshape(Sp @ (Xi - Xj)[:, :, np.newaxis], Xi.shape)

        ScXc = Sc + Xc[:, :, np.newaxis] * Xc[:, np.newaxis, :]
        DxSp = Dx[:, :, np.newaxis] * Dx[:, np.newaxis, :] - Sp
        GjGi = Gj[:, :, :, np. newaxis, np. newaxis] * Gi[:, np. newaxis, np. newaxis, :, :]
        DSX = DS[:, :, :, np.newaxis] * Xc[:, np.newaxis, np.newaxis, :]

        ndim = Xc.shape[-1]

        term = lambda n, u, m, v: coriolis_tensors[:, n, u, m, v] * cls.annoying_coriolis_term(
            n, u, m, v,

            Sc, Dx,
            ScXc, DxSp, GjGi,
            DSX
        )

        contrib = np.zeros(Sc.shape[0])
        inds = itertools.combinations_with_replacement(range(ndim), r=4)
        for n,u,m,v in inds:
            # upper triangle so we need to do all the combos
            contrib += (
                term(n, u, m, v)
                + term(n, u, v, m)
                + term(u, n, v, m)
                + term(u, n, m, v)
            )

        npts = X0.shape[0]
        ke = np.zeros((npts, npts))
        ke[rows, cols] = contrib
        ke[cols, rows] = contrib

        return ke


    def evaluate_ke(self, overlap_data):

        # #[0.00881508 0.00808985 0.00661593 0.00479724 0.00300391 0.00148472]
        # return self.evaluate_classic_momentum_contrib(
        #     overlap_data,
        #     np.ones(overlap_data['init_alphas'].shape[-1])
        # )

        base = self.evaluate_diagonal_rotated_momentum_contrib(
            overlap_data,
            np.ones(overlap_data['init_alphas'].shape[-1])
        )

        B_e, coriolis = self.ci_func(overlap_data['centers'])
        coriolis = -self.evaluate_coriolis_contrib(coriolis, overlap_data)

        watson = np.zeros(base.shape)
        rows, cols = np.triu_indices_from(base)
        watson[rows, cols] = -B_e
        watson[cols, rows] = -B_e

        ke = base + coriolis + watson

        return ke

        # coriolis = self.ci_func(overlap_data['centers'])
        # watson = self._evaluate_polynomial_ke(overlap_data, ['q','p', 'q', 'p'], coriolis)
        #
        # n = coriolis.shape[-1]
        # prefactors = np.broadcast_to(
        #     np.eye(n),
        #     (overlap_data['centers'].shape[0], n, n)
        # )
        # base = self._evaluate_polynomial_ke(overlap_data, ['p', 'p'], prefactors)
        # ke = base + watson
        #
        # return ke

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
                            overlap_data,
                            expansion_type,
                            expansion_degree=2,
                            pairwise_functions=None,
                            logger=None
                            ):

        if pairwise_functions is not None:
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
                alphas = overlap_data['alphas']
                centers = overlap_data['centers']
                deriv_order = deriv_order - (deriv_order % 2) # odd orders don't contribute so why evaluate the derivatives...
                logger.log_print("expanding about {N} points...", N=len(alphas))
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
                    # if i == 0 or self.inds is None:
                    _.append(d1 - d2)
                    # else:
                    #     d1 = d1.copy()
                    #     idx = (slice(None, None, None),)+ np.ix_(*[self.inds,] * i)
                    #     d1[idx] -= d2
                    #     _.append(d1)
                derivs = _
            #     raise Exception(derivs[2][:5])
            # else:
            #     raise Exception(derivs[2][:5])

            npts = len(overlap_data['init_centers'])
            pot = cls.tensor_expansion_integrate(
                npts,
                derivs,
                overlap_data,
                reweight=True,
                expansion_type=expansion_type,
                logger=logger
            )

        else:
            npts = len(overlap_data['init_centers'])
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
            overlap_data,
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
            logger.log_print("evauating integrals with {n}-order quadrature", n=quadrature_degree)
            pot_mat = cls.quad_integrate(
                function,
                overlap_data,
                degree=quadrature_degree,
                logger=logger
            )
        elif integral_handler == 'expansion':
            logger.log_print("evauating integrals with {n}-degree expansions", n=expansion_degree)
            pot_mat = cls.expansion_integrate(function,
                                              overlap_data,
                                              expansion_degree=expansion_degree,
                                              expansion_type=expansion_type,
                                              pairwise_functions=pairwise_functions,
                                              logger=logger
                                              )
        elif integral_handler == 'analytic':
            logger.log_print("evauating integrals analytically", n=expansion_degree)
            pot_mat = cls.analytic_integrate()
        else:
            raise ValueError("unknown operator evaluation scheme {}".format(integral_handler))

        return pot_mat

    def evaluate_pe(self, overlap_data):
        return self.evaluate_multiplicative(
            self.potential_function,
            overlap_data,
            integral_handler=self.handler,
            quadrature_degree=self.quadrature_degree,
            expansion_degree=self.expansion_degree,
            expansion_type=self.expansion_type,
            pairwise_functions=self.pairwise_handler,
            logger=self.logger
        )

    def evaluate_op(self,
                    operator,
                    overlap_data,
                    integral_handler=None,
                    expansion_degree=None,
                    expansion_type=None,
                    quadrature_degree=None,
                    # pairwise_functions=None,
                    logger=None
                    ):
        return self.evaluate_multiplicative(
            operator,
            overlap_data,
            integral_handler=self.handler if integral_handler is None else integral_handler,
            quadrature_degree=self.quadrature_degree if quadrature_degree is None else quadrature_degree,
            expansion_degree=self.expansion_degree if expansion_degree is None else expansion_degree,
            expansion_type=self.expansion_type if expansion_type is None else expansion_type,
            # pairwise_functions=pairwise_functions
            logger=self.logger if logger is None else logger
        )

class DGBPairwisePotentialEvaluator(DGBEvaluator, metaclass=abc.ABCMeta):
    def __init__(self, coords, pairwise_potential_functions, quadrature_degree=3):
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

    def get_bond_length_change_transformation(self, overlap_data, i, j) -> np.ndarray:
        d0, base_proj = self.get_coordinate_bond_length_projection(i, j)
        proj_data = base_proj[np.newaxis] @ overlap_data['rotations'] # not sure if we need a .T or not...
        return self.get_coordinate_change_transformation(proj_data)

    def wrap_distance_function(self, i, j, overlap_data, transformations, pairwise_function):
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
                                  overlap_data,
                                  quadrature_degree=None,
                                  expansion_degree=2
                                  ):

        if quadrature_degree is None:
            quadrature_degree = self.quadrature_degree

        centers = overlap_data['centers']
        alphas = overlap_data['alphas']
        covs = overlap_data['inverse']

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
            tf_cov = np.linalg.inv(tfs @ covs @ tfsT)

            tf_alphas, tf_vecs = np.linalg.eigh(tf_cov)
            tf_alphas = tf_alphas / 2 # baked into the covariences...
            ndim = tf_alphas.shape[1]

            tf_centers = tfs @ centers[:, :, np.newaxis]
            tf_centers = tf_centers.reshape(tf_centers.shape[:2])

            # # we have, in principle, scaling = sqrt(pi^(d - k) |A| / |A_r| )
            # #           but in this |A| = prod(alphas*2)
            # #                       |A_r| = prod(tf_alphas*2)
            # tf_rats = np.prod(alphas[:, :ndim] / tf_alphas, axis=1) * np.prod(alphas[:, ndim:], axis=1)
            # # scaling = np.sqrt( tf_rats * (2*np.pi)**(alphas.shape[1] - ndim) ) # accounting for prefactor

            pairwise_contrib = self.rotated_gaussian_quadrature(
                f,
                tf_alphas,
                tf_centers,
                tf_vecs,
                tf_vecs.transpose(0, 2, 1),
                degree=quadrature_degree,
                normalize=False
            )

            # pairwise_contrib = self.quad_nd(
            #     tf_centers,
            #     tf_alphas,
            #     f,
            #     degree=quadrature_degree,
            #     normalize=False
            # ) / np.sqrt(np.pi**ndim)
            #
            # print(pairwise_contrib[:3])

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
        natoms = self.coords.natoms
        modes = self.coords.modes.matrix
        # invs = self.coords.modes.inverse
        ndim = modes.shape[0] // natoms # this will almost always be 3
        base_mat = self.get_bond_length_deltas(natoms, ndim, i, j)
        tf_base = base_mat @ modes
        # print("="*50)
        # print(base_mat)
        # print([modes, self.coords.modes.origin])
        # print(tf_base)
        # print("-"*50)
        # inverse_base = base_mat.T / 2
        # tf_inv = invs @ inverse_base
        d0 = np.dot(base_mat, self.coords.modes.origin.reshape(-1))
        return d0, tf_base

class DGBCoords(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def centers(self) -> 'np.ndarray':
        ...
    @property
    def shape(self):
        return self.centers.shape

    @property
    @abc.abstractmethod
    def kinetic_energy_evaluator(self) -> 'DGBKineticEnergyEvaluator':
        ...

    @property
    @abc.abstractmethod
    def pairwise_potential_evaluator_type(self) -> 'type[DGBPairwisePotentialEvaluator]':
        ...
    def pairwise_potential_evaluator(self, potential_functions) -> 'DGBPairwisePotentialEvaluator':
        if not isinstance(potential_functions, dict):
            raise ValueError("don't know what to do with pairwise functions {}".format(potential_functions))
        if 'functions' in potential_functions:
            opts = potential_functions.copy()
            potential_functions = potential_functions['functions']
            del opts['functions']
        else:
            opts = {}
        return self.pairwise_potential_evaluator_type(self, potential_functions, **opts)
    @abc.abstractmethod
    def __getitem__(self, item) -> 'DGBCoords':
        ...
    def take_indices(self, subinds) -> 'DGBCoords':
        return self[:, subinds] # default
    def drop_indices(self, subinds) -> 'DGBCoords':
        remaining = np.setdiff1d(np.arange(self.centers.shape[-1]), subinds)
        return self.take_indices(remaining) # default
    @abc.abstractmethod
    def gmatrix(self, coords:np.ndarray) -> np.ndarray:
        ...

    @classmethod
    def embedded_mode_function(cls, func, modes, natoms=None):
        @functools.wraps(func)
        def embedded_function(mode_coords, deriv_order=None):
            origin = modes.origin
            carts = (mode_coords[:, np.newaxis, :] @ modes.matrix.T[np.newaxis]).reshape(
                mode_coords.shape[:1] + origin.shape
            )
            carts = carts + origin
            if natoms is not None:
                carts = carts.reshape((carts.shape[0], natoms, -1))

            vals = func(carts, deriv_order=deriv_order)
            if deriv_order is not None:
                fshape = vals[0].ndim
                _ = []
                for n, d in enumerate(vals):
                    for j in range(n):
                        d = np.tensordot(d, modes.matrix, axes=[fshape, 0])
                    _.append(d)
                vals = _
            return vals

        return embedded_function

    @classmethod
    def embedded_subcoordinate_function(cls, func, sel, ndim):
        @functools.wraps(func)
        def embedded_function(subcoords, deriv_order=None):
            full_shape = (subcoords.shape[0], ndim)
            if sel is None:
                coords = subcoords.reshape(full_shape)
            else:
                coords = np.zeros(full_shape)
                coords[:, sel] = subcoords
            vals = func(coords, deriv_order=deriv_order)
            if sel is not None and deriv_order is not None:
                fshape = vals[0].ndim
                _ = []
                for n, d in enumerate(vals):
                    for j in range(n):
                        d = np.take(d, sel, axis=j + fshape)
                    _.append(d)
                vals = _
            return vals

        return embedded_function

    @classmethod
    def embedded_cartesian_function(cls, func, atom_sel, xyz_sel, natoms, ndim):
        @functools.wraps(func)
        def embedded_function(subcart_coords, deriv_order=None):
            subcart_coords = np.reshape(
                subcart_coords,
                (
                    subcart_coords.shape[0],
                    len(atom_sel) if atom_sel is not None else natoms,
                    len(xyz_sel) if xyz_sel is not None else ndim
                )
            )
            full_shape = (subcart_coords.shape[0], natoms, ndim)
            if atom_sel is None and xyz_sel is None:
                carts = subcart_coords.reshape(full_shape)
            else:
                carts = np.zeros(full_shape)
                if atom_sel is None:
                    carts[..., xyz_sel] = subcart_coords
                elif xyz_sel is None:
                    carts[..., atom_sel, :] = subcart_coords
                else:
                    carts[..., np.ix_(atom_sel, xyz_sel)] = subcart_coords
            vals = func(carts, deriv_order=deriv_order)
            if (atom_sel is not None or xyz_sel is not None) and deriv_order is not None:
                full_sel = np.arange(natoms * ndim).reshape(natoms, ndim)
                if atom_sel is None:
                    full_sel = full_sel[:, xyz_sel]
                elif xyz_sel is None:
                    full_sel = full_sel[atom_sel, :]
                else:
                    full_sel = full_sel[np.ix_(atom_sel, xyz_sel)]
                flat_sel = full_sel.flatten()

                fshape = vals[0].ndim
                _ = []
                for n, d in enumerate(vals):
                    for j in range(n):
                        d = np.take(d, flat_sel, axis=j + fshape)
                    _.append(d)
                vals = _
            return vals

        return embedded_function

    class DGBEmbeddedFunction:
        def __init__(self, embedded_function, original_function, coords):
            self.og_fn = original_function
            self.base_coords = coords
            self.embed_fn = embedded_function
        def __call__(self, coords, deriv_order=None):
            return self.embed_fn(coords, deriv_order=deriv_order)
    @abc.abstractmethod
    def embed_function(self, fn) -> 'DGBEmbeddedFunction':
        ...

    def as_cartesians(self) -> 'tuple[DGBCartesians, tuple[np.ndarray, np.ndarray]]':
        raise NotImplementedError("{} can't be converted to Cartesians".format(
            type(self).__name__
        ))

class DGBCartesians(DGBCoords):
    def __init__(self, coords, masses, *, natoms=None, atom_sel=None, ndim=None, xyz_sel=None):
        coords = np.asanyarray(coords)
        masses = np.asanyarray(masses)
        if masses.ndim > 1:
            raise ValueError("expected a vector of masses...")
        coords = coords.reshape((len(coords), len(masses), -1))

        self.coords = coords
        self.masses = masses
        self.natoms = len(masses) if natoms is None else natoms# the original numebr of atoms
        self.ndim = self.coords.shape[-1] if ndim is None else ndim
        self.atom_sel = atom_sel
        self.xyz_sel = xyz_sel
    @property
    def centers(self):
        return self.coords.reshape((self.coords.shape[0], self.coords.shape[1]*self.coords.shape[2]))
    @property
    def cart_shape(self):
        return self.coords.shape
    @property
    def kinetic_energy_evaluator(self):
        return DGBCartesianEvaluator(
            np.broadcast_to(self.masses[:, np.newaxis], self.coords.shape[1:]).flatten()
        )
    @property
    def pairwise_potential_evaluator_type(self):
        return DGBCartesianPairwiseEvaluator
    @classmethod
    def resolve_masses(cls, coords, masses=None, atoms=None):
        if masses is None:
            if atoms is not None:
                atoms = [
                    AtomData[a, "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
                        if isinstance(a, str) else
                    a
                    for a in atoms
                ]
                masses = np.array(atoms)
            else:
                raise ValueError("either atoms or masses must not be `None`")

        if isinstance(masses, (int, np.integer, float, np.floating)):
            masses = [masses] * coords.shape[1]
        return masses
    @classmethod
    def from_cartesians(cls, centers, masses=None, atoms=None):
        centers = np.asanyarray(centers)
        masses = cls.resolve_masses(centers, masses=masses, atoms=atoms)
        return cls(centers, masses)

    def infer_shape_sel(self, selector):
        n, na, nd = self.coords.shape
        test = np.arange(na*nd).reshape((1, na, nd))
        inferrable = np.broadcast_to(test, self.coords.shape)[selector][0]

        test_atoms = test[0, :, 0]
        sub_atoms = inferrable[:, 0]
        # we map the atoms back to their positions in the original list
        atom_sel = np.searchsorted(test_atoms, sub_atoms)
        if len(atom_sel) == na and np.all(atom_sel == np.arange(na)):
            atom_sel = None

        test_xyz = test[0, 0] % nd
        sub_xyz = inferrable[0] % nd
        xyz_sel, _ = nput.find(test_xyz, sub_xyz)
        if len(xyz_sel) == nd and np.all(xyz_sel == np.arange(nd)):
            xyz_sel = None

        return atom_sel, xyz_sel
    @staticmethod
    def _merge_sel(new_sel, old_sel):
        if old_sel is not None:
            new_sel = np.asanyarray(old_sel)[new_sel,]
        return new_sel
    def __getitem__(self, item) -> 'DGBCartesians':
        c = self.coords[item]
        if not isinstance(c, np.ndarray) or c.ndim != 3:
            raise ValueError("bad slice {}".format(item))

        if c.ndim < 3:
            raise ValueError("zero-length slice?")

        if c.shape[1:] == self.coords.shape[1:]:
            return type(self)(
                c,
                self.masses,
                atom_sel=self.atom_sel,
                xyz_sel=self.xyz_sel,
                natoms=self.natoms,
                ndim=self.ndim
            )
        else:
            atom_sel, xyz_sel = self.infer_shape_sel(item)
            return type(self)(
                c,
                self.masses if atom_sel is None else self.masses[atom_sel,],
                atom_sel=self._merge_sel(atom_sel, self.atom_sel),
                xyz_sel=self._merge_sel(xyz_sel, self.xyz_sel),
                natoms=self.natoms,
                ndim=self.ndim
            )
    def take_indices(self, subinds):
        subinds = np.asanyarray(subinds)
        atom_inds = np.unique(subinds // self.coords.shape[2])
        xyz_inds = np.unique(subinds % self.coords.shape[2])
        return self[:, atom_inds, :][:, :, xyz_inds]
    def embed_function(self, function):
        """
        Embeds assuming we got a function in Cartesians _before_ any selections happened

        :param function:
        :return:
        """
        return self.DGBEmbeddedFunction(
            self.embedded_cartesian_function(
                function,
                self.atom_sel,
                self.xyz_sel,
                self.natoms,
                self.ndim
            ),
            function,
            self
        )

    def gmatrix(self, coords:np.ndarray) -> np.ndarray:
        mass_spec = np.broadcast_to(self.masses[:, np.newaxis], (len(self.masses), self.coords.shape[2])).flatten()
        mass_spec = np.diag(1 / mass_spec)
        return np.broadcast_to(mass_spec[np.newaxis], (coords.shape[0],) + mass_spec.shape)
class DGBInternals(DGBCoords):
    def __init__(self, coords, gmat_function=None, vprime_function=None):
        raise NotImplementedError('internals coming soon?')
class DGBWatsonModes(DGBCoords):
    def __init__(self, coords, modes,
                 *,
                 coriolis_inertia_function=None,
                 natoms=None,
                 subselection=None
                 ):
        self.coords = coords
        self.natoms = natoms
        self.modes = modes
        self.subsel = subselection
        self.ci_func = coriolis_inertia_function

    @property
    def centers(self):
        return self.coords
    @property
    def kinetic_energy_evaluator(self):
        return DGBWatsonEvaluator(self.modes, self.ci_func)
    @property
    def pairwise_potential_evaluator_type(self):
        return DGBWatsonPairwiseEvaluator
    @staticmethod
    def zeta_momi(watson_coords, modes, masses):
        origin = modes.origin
        carts = (watson_coords[:, np.newaxis, :] @ modes.matrix.T[np.newaxis]).reshape(watson_coords.shape[:1] + origin.shape)
        carts = carts + origin[np.newaxis]
        carts = carts.reshape((carts.shape[0], len(masses), -1))
        if carts.shape[-1] < 3:
            # pad with zeros
            carts = np.concatenate(
                [
                    carts,
                    np.zeros(carts.shape[:2] + (3 - carts.shape[-1],))
                ],
                axis=-1
            )

        zeta, (B_e, eigs) = StructuralProperties.get_prop_coriolis_constants(carts,
                                                                             modes.inverse,
                                                                             masses
                                                                             )

        return zeta, B_e

    @classmethod
    def default_coriolis_intertia_function(cls, modes, masses):
        def coriolis_inertia_function(watson_coords, *, modes=modes, masses=masses, deriv_order=None):
            if deriv_order is not None:
                raise NotImplementedError("don't have support for Coriolis derivatives in DGB")
            zeta, mom_i = cls.zeta_momi(watson_coords, modes, masses)
            B_e = 1 / (2 * mom_i)  # * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass"))
            return np.sum(B_e, axis=-1), sum(
                B_e[:, a, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
                 * zeta[:, a, :, :, np. newaxis, np.newaxis]*zeta[:, a, np. newaxis, np.newaxis, :, :]
                for a in range(B_e.shape[1])
            )

        return coriolis_inertia_function

    @classmethod
    def from_cartesians(cls,
                        coords,
                        modes,
                        masses=None,
                        coriolis_inertia_function=None
                        ):
        flat_carts = (coords - modes.origin[np.newaxis]).reshape((len(coords), -1))
        coords = (flat_carts[:, np.newaxis, :] @ modes.inverse.T[np.newaxis]).reshape(
            flat_carts.shape[0],
            modes.inverse.shape[0]
        )
        if coriolis_inertia_function is None:
            coriolis_inertia_function = cls.default_coriolis_intertia_function(modes, masses)
        return cls(coords,
                   modes,
                   coriolis_inertia_function=coriolis_inertia_function,
                   natoms=len(masses)
                   )
    def as_cartesians(self, masses=None) -> 'tuple[DGBCartesians, tuple[np.ndarray, np.ndarray]]':
        mode_coords = self.coords
        modes = self.modes

        origin = modes.origin
        carts = (mode_coords[:, np.newaxis, :] @ modes.matrix.T[np.newaxis]).reshape(
            mode_coords.shape[:1] + origin.shape
        )
        carts = carts + origin

        natoms = self.natoms
        if natoms is None:
            raise ValueError("need `natoms` to be able to convert back to Cartesians")
        carts = carts.reshape((carts.shape[0], natoms, -1))

        if masses is None:
            masses = [-1] * natoms # so things break

        return DGBCartesians(carts, masses), (modes.matrix, modes.inverse)

    def __getitem__(self, item):
        c = self.coords[item]
        if not isinstance(c, np.ndarray) or c.ndim != self.coords.ndim:
            raise ValueError("bad slice {}".format(item))
        if c.shape[1] != self.coords.shape[1]:
            # means we took some subselection of modes
            test_sel = np.broadcast_to(
                np.arange((self.coords.shape[1]))[np.newaxis],
                self.coords.shape
            ) # something that can be indexed the same way
            subsel = test_sel[item][0]
            modes = self.modes[subsel]
            ci_func = self.embedded_subcoordinate_function(
                self.ci_func,
                subsel,
                self.coords.shape[1]
            )
        else:
            modes = self.modes
            subsel = None
            ci_func = self.ci_func
        if self.subsel is not None:
            subsel = self.subsel[subsel]
        return type(self)(
            c,
            modes,
            coriolis_inertia_function=ci_func,
            natoms=self.natoms,
            subselection=subsel
        )
    def gmatrix(self, coords:np.ndarray) -> np.ndarray:
        base_spec = np.eye(coords.shape[1])
        return np.broadcast_to(base_spec[np.newaxis], (coords.shape[0],) + base_spec.shape)

    def embed_function(self, fn):
        return self.DGBEmbeddedFunction(
            self.embedded_mode_function(  # a function that takes in normal mode coordinates
                fn,
                self.modes,
                natoms=self.natoms
            ),
            fn,
            self
        )

class DGBCovarianceTransformations:

    def __init__(self, tfs, invs=None, inverse_sum=None):
        self.tfs = tfs
        if invs is None and self.tfs is not None:
            dets = np.linalg.det(tfs)
            if not np.allclose(dets, np.ones(len(dets))):
                raise ValueError("inverse must be supplied for non-unitary transforms")
            else:
                invs = tfs.transpose(0, 2, 1)
        self.invs = invs
        self.inverse_sum = inverse_sum

    def transform(self, derivative_operators):
        raise NotImplementedError("proper transformation support still to come")

    def inverse(self):
        return type(self)(self.invs, self.tfs)

    def combine(self, other:'DGBCovarianceTransformations'):
        if self.tfs is None and other.tfs is None:
            tfs = None
            invs = None
            inverse_sum = None
        elif self.tfs is None or other.tfs is None:
            raise NotImplementedError("if one set of covariances is specified, both must be")
        else:
            inverse_sum = np.linalg.inv(self.invs + other.invs)
            tfs = self.tfs + other.tfs
            invs = other.invs - other.invs @ inverse_sum @ other.invs

        return type(self)(tfs, invs, inverse_sum=inverse_sum)

class DGBGaussians:
    """
    A class to set up the actual N-dimensional Gaussians used in a DGB
    """

    def __init__(self, coords, alphas, transformations=None, *, poly_coeffs=None, logger=None):
        self._S = None
        self._T = None

        if not isinstance(coords, DGBCoords):
            raise ValueError("need DGBCoords passed in")

        self.coords = coords
        if isinstance(alphas, (int, float, np.integer, np.floating)):
            alphas = [alphas] * self.coords.shape[0]
        self.alphas = np.asanyarray(alphas)

        if self.alphas.ndim == 1:
            self.alphas = np.broadcast_to(self.alphas[:, np.newaxis], self.coords.shape)
        self._transforms = self.canonicalize_transforms(self.alphas, transformations)

        self._overlap_data = None

        if poly_coeffs is not None:
            poly_coeffs = self.canonicalize_poly_coeffs(poly_coeffs, self.alphas)
        self._poly_coeffs = poly_coeffs

        self.logger = Logger.lookup(logger)

    @property
    def overlap_data(self):
        if self._overlap_data is None:
            self._overlap_data = DGBEvaluator.get_overlap_gaussians(
                self.coords.centers,
                self.alphas,
                self.transformations
            )
        return self._overlap_data
    def get_S(self):
        if self._poly_coeffs is not None:
            raise NotImplementedError("need to reintroduce polynomial support")
        return DGBEvaluator.evaluate_overlap(self.overlap_data)
    def get_T(self):
        if self._poly_coeffs is not None:
            raise NotImplementedError("need to reintroduce polynomial support")
        self.logger.log_print("evaluating kinetic energy...")
        return self.coords.kinetic_energy_evaluator.evaluate_ke(self.overlap_data)

    def get_base_polynomials(self):
        base_coeffs = self._poly_coeffs
        if self.inds is not None:
            comp = np.setdiff1d(np.arange(self.centers.shape[-1]), self.inds)
            new_pc = []
            for c in base_coeffs:
                if c is not None:
                    c = c.copy()
                    for idx in comp:
                        c = np.moveaxis(c, idx, 0)
                        c = c[:1].reshape(c.shape[1:])
                        # c = np.moveaxis(c, 0, idx)
                new_pc.append(c)
            base_coeffs = new_pc
        base_polys = [
            DensePolynomial(coeffs) if coeffs is not None else None
            for coeffs in base_coeffs
        ]
        centers = self.centers
        if self.transformations is not None:
            tfs, invs = self.transformations
            if self.inds is not None:
                tfs = tfs[:, :, self.inds]
                invs = invs[:, self.inds, :]
            # raise Exception(invs[0].shape)
            base_polys = [
                p.transform(inv) if p is not None else 1
                for p,inv in zip(base_polys, tfs)
            ]
        elif self.inds is not None:
            centers = centers[..., self.inds]
        base_polys = [
                p.shift(-c) if isinstance(p, DensePolynomial) else 1
                for p,c in zip(base_polys, centers)
            ]
        # raise Exception(base_polys[0].coeffs)
        return base_polys

    def optimize(self, optimizer_options, **opts):
        if optimizer_options is True:
            optimizer_options = 'gram-schmidt'
        if isinstance(optimizer_options, (int, float, np.integer, np.floating)):
            optimizer_options = {'method':'gram-schmidt', 'overlap_cutoff':optimizer_options}
        if isinstance(optimizer_options, str):
            optimizer_options = {'method':optimizer_options}

        # dispatch
        method = optimizer_options['method'].lower()
        opts = dict(optimizer_options, **opts)
        del opts['method']
        if method == 'svd':
            optimized_positions = self._optimize_svd(**opts)
        elif method == 'gram-schmidt':
            optimized_positions = self._optimize_gs(**opts)
        else:
            raise ValueError("don't know what to do with optimizer method '{}'".format(method))

        return self.take_gaussian_selection(optimized_positions)

    def take_gaussian_selection(self, full_good_pos):
        centers = self.coords[full_good_pos,]
        alphas = self.alphas[full_good_pos,]
        if self.transformations is not None:
            transformations = [t[full_good_pos,] for t in self.transformations]
        else:
            transformations = None

        if self._poly_coeffs is not None:
            _poly_coeffs = [t[full_good_pos,] for t in self._poly_coeffs]
        else:
            _poly_coeffs = None

        return type(self)(
            centers, alphas, transformations,
            poly_coeffs=_poly_coeffs, logger=self.logger
        )

    def _calc_gs_norm(self, new, S, old):
        # calculates a very particular form of norm where (for speed purposes)
        #   1. new is the new _unnormalized_ orthogonal vector
        #   2. old is the row of the overlap matrix pulled to orthogonalize
        #   3. S is the subblock of the overlap matrix corresponding to everything up to the newest element
        #
        # This means we get to calculate a basic norm corresponding to the previous block by using everything
        # up to the final element of new and then we add on a correction knowing that the missing piece of S
        # is encoded in old

        vec = new[:-1]
        rows, cols = np.triu_indices(len(vec), k=1)
        off_diag = 2 * np.dot(S[rows, cols], vec[rows] * vec[cols])
        diag = np.dot(vec, vec) # taking advantage of the fact that the diagonal of S is 1
        base_norm = diag + off_diag

        # correction from old itself, knowing that the final element of new is 1
        correction_norm = 1 + 2 * np.dot(old[:-1], vec)

        return base_norm + correction_norm
    def _calc_gs_next(self, vec, ortho):
        new = np.zeros(len(vec))
        new[-1] = 1
        for prev in ortho:
            n = len(prev)
            new[:n] -= np.dot(vec[:n], prev) * prev
        return new

    gs_optimization_overlap_cutoff=1e-3
    def _optimize_gs(self, *, S=None, overlap_cutoff=None, norm_truncation_cutoff=1e-6):
        if overlap_cutoff is None:
            overlap_cutoff = self.gs_optimization_overlap_cutoff
        if S is None:
            S=self.S
        mask = np.full(len(S), False, dtype=bool)
        mask[0] = True
        orthog_vecs = [np.array([1.])]
        norms = [] # debug
        for i in range(1, len(S)):
            mask[i] = True
            vec = S[i, :i+1][mask[:i+1]]
            new = self._calc_gs_next(vec, orthog_vecs)
            norm = self._calc_gs_norm(new, S[np.ix_(mask[:i], mask[:i])], vec)
            norms.append(norm)
            if norm > overlap_cutoff:
                if norm < norm_truncation_cutoff:
                    norm = norm_truncation_cutoff
                orthog_vecs.append(new / np.sqrt(norm))
            else:
                mask[i] = False

        # raise Exception(norms)

        return np.where(mask)[0]

    def _optimize_svd(self, min_singular_value=1e-4, num_svd_vectors=None, svd_contrib_cutoff=1e-3):

        if min_singular_value is None and num_svd_vectors is None:
            raise ValueError("either a minimum singular value or number of eigenvectors needs to be passes")
        sig, evecs = np.linalg.eigh(self.get_S())
        if num_svd_vectors:
            good_loc = slice(max(len(sig) - num_svd_vectors, 0), len(sig))
        else:
            # self.logger.log_print("most important center threshold: {t}", t=min_singular_value)
            good_loc = np.where(sig > min_singular_value)[0]
        # raise Exception(np.min(np.abs(U)))
        full_good_pos = np.unique(np.where(np.abs(evecs[:, good_loc]) > svd_contrib_cutoff)[0])

        return full_good_pos

    @classmethod
    def construct(cls,
                  coords,
                  alphas,
                  potential_function=None,
                  transformations=None,
                  masses=None,
                  atoms=None,
                  modes=None,
                  internals=None,
                  coordinate_selection=None,
                  cartesians=None,
                  gmat_function=None,
                  poly_coeffs=None,
                  logger=None
                  ):
        if not isinstance(coords, DGBCoords):
            if modes is not None:
                if isinstance(modes, str):
                    if modes == 'normal':
                        modes = cls.get_normal_modes(
                            coords,
                            potential_function,
                            masses=masses,
                            atoms=atoms,
                            internals=internals,
                            gmat_function=gmat_function
                        )

                    else:
                        raise ValueError("unknown mode spec {}".format(modes))

                if coords.ndim == 3: # Cartesian-like
                    # coords_old = coords

                    coords = DGBWatsonModes.from_cartesians(
                        coords,
                        modes,
                        masses=DGBCartesians.resolve_masses(
                            coords,
                            masses=masses,
                            atoms=atoms
                        )
                    )

                else:
                    if internals is not None:
                        raise NotImplementedError("internal modes not implemented")
                    else:
                        raise NotImplementedError("transformation not fully implemented")
                        # Is this right???
                        coords = DGBWatsonModes(
                            modes.matrix @ coords,
                            modes=modes
                        )

            elif internals is not None:
                raise NotImplementedError("internal coordinate mapping not supported yet...")
            else:
                coords = DGBCartesians.from_cartesians(coords, masses=masses, atoms=atoms)
                if cartesians is not None:
                    coords = coords[:, :, cartesians]

            if coordinate_selection is not None:
                coords = coords.take_indices(coordinate_selection)

            potential_function = coords.embed_function(potential_function)

        if transformations is not None:
            if isinstance(transformations, str):
                if transformations == 'reaction_path':
                    transformations = cls.get_reaction_path_transformations(
                        coords,
                        potential_function,
                        masses=masses,
                        atoms=atoms,
                        internals=internals,
                        gmat_function=gmat_function
                    )
                else:
                    raise ValueError("unknown transformation spec {}".format(modes))

        if isinstance(alphas, (str, dict)):
            if isinstance(alphas, str) and alphas == 'auto':
                if modes is not None:
                    alphas = 'virial'
                else:
                    alphas = 'masses'
            alphas = cls.dispatch_get_alphas(
                alphas,
                coords,
                masses=masses,
                gmat_function=coords.gmatrix if gmat_function is None else gmat_function,
                potential_function=potential_function,
                transformations=transformations
            )

        shp = coords.shape
        if isinstance(alphas, (int, float, np.integer, np.floating)):
            alphas = np.full(shp, alphas)
        else:
            alphas = np.asanyarray(alphas)
            alphas = np.broadcast_to(alphas, shp)

        return cls(coords, alphas, transformations,
                   poly_coeffs=poly_coeffs, logger=logger
                   ), potential_function
    @classmethod
    def get_normal_modes(cls,
                         coords,
                         potential_function,
                         masses=None,
                         atoms=None,
                         internals=None,
                         gmat_function=None,
                         reference_structure=None,
                         stationary_point_norm=1e-6
                         ):
        if internals is not None:
            raise ValueError("internal coordinate support still to come")
        if reference_structure is None:
            f_data = potential_function(coords, deriv_order=2)
            if isinstance(f_data, np.ndarray) and f_data.ndim == 3:
                hess = f_data
                grad = potential_function(coords, deriv_order=1)
                pot = potential_function(coords, deriv_order=1)
            else:
                pot, grad, hess = f_data
            sorting = np.argsort(pot)
            coords = coords[sorting,]
            grad = grad[sorting,]
            hess = hess[sorting,]
            grad_norms = np.linalg.norm(grad, axis=-1)
            stationary_structures = np.where(grad_norms < stationary_point_norm)
            if len(stationary_structures) == 0 or len(stationary_structures[0]) == 0:
                raise ValueError("no stationary state found?")
            stationary_structures = stationary_structures[0]
            reference_structure = coords[stationary_structures[0]]
            f_matrix = hess[stationary_structures[0]]
        else:
            f_matrix = potential_function(reference_structure, deriv_order=2)
            if not isinstance(f_matrix, np.ndarray):
                f_matrix = f_matrix[2]

        if gmat_function is not None:
            mass_spec = gmat_function(reference_structure)
        else:
            mass_spec = DGBCartesians.resolve_masses(coords, masses=masses, atoms=atoms)

        if coords.shape[-1] == 3:
            basis = CartesianCoordinates3D
        elif coords.shape[-1] == 2:
            basis = CartesianCoordinates2D
        elif coords.shape[-1] == 1:
            basis = CartesianCoordinates1D
        else:
            raise ValueError("higher than 3D Cartesians?")

        return NormalModes.from_fg(
            basis,
            f_matrix,
            mass_spec,
            origin=reference_structure
        )

    @classmethod
    def get_reaction_path_transformations(
            cls,
            coords,
            potential_function,
            gmat_function,
            stationary_point_norm=1e-6
    ):
        f_data = potential_function(coords, deriv_order=2)
        if isinstance(f_data, np.ndarray) and f_data.ndim == 3:
            hess = f_data
            grad = potential_function(coords, deriv_order=1)
            pot = potential_function(coords, deriv_order=1)
        else:
            pot, grad, hess = f_data

        gmats = gmat_function(coords)
        gvals, gvecs = np.linalg.eigh(gmats)
        if np.any((gvals <= 0).flatten()):
            raise ValueError("bad G-matrix?")
        g12_diags = np.zeros(gvecs.shape)
        diag_ings = (slice(None),) + np.diag_indices_from(gvecs[0])
        g12_diags[diag_ings] = np.sqrt(g12_diags)
        g12 = gvecs @ g12_diags @ gvecs.T

        grad = grad @ g12  # mass-weight
        hess = g12 @ hess @ g12


        grad_norms = np.linalg.norm(grad, axis=-1)
        # stationary_structures = np.where(grad_norms < stationary_point_norm)
        #
        # if len(stationary_structures) == 0:
        #     stationary_structures = np.array([], dtype=int)
        # else:
        #     stationary_structures = stationary_structures[0]
        # non_stationary = np.setdiff1d(np.arange(len(grad_norms)), stationary_structures)

        non_stationary = np.where(grad_norms >= stationary_point_norm)

        # in this case I'm _only_ projecting out the gradient for whatever that's worth...
        # which...probably has some rotation/translation component which
        # is why I don't get clean zero eigenvectors...
        rp_mode = grad[non_stationary,] / grad_norms[non_stationary,][:, np.newaxis]
        proj = nput.vec_outer(rp_mode, rp_mode, axes=[1, 1])
        hess[non_stationary,] -= proj@hess[non_stationary,]@proj

        freqs, g12_tfs = np.eigh(hess) # replace zero tf with mass-weighted rp_mode
        tfs = g12 @ g12_tfs

        g12_diags[diag_ings] = 1 / np.sqrt(g12_diags)
        g12_inv = gvecs @ g12_diags @ gvecs.T
        inv = g12_tfs.T @ g12_inv

        return tfs, inv

    @staticmethod
    def _filter_alpha_method_keys(method, opts, necessary_keys, optional_keys):
        new_opts = {}
        for k in necessary_keys:
            if k not in opts:
                raise ValueError(
                    "alpha inference method {} needs keys {}, missing {}".format(method, necessary_keys, k)
                )
            new_opts[k] = opts[k]
        for k in optional_keys:
            if k in opts: new_opts[k] = opts[k]
        return new_opts
    @classmethod
    def dispatch_get_alphas(self, alphas, centers, **extra_opts):
        if isinstance(alphas, str):
            alphas = {'method':alphas}
        opts = dict(extra_opts, **alphas)
        del opts['method']
        method = alphas['method']
        if method == 'virial':
            return self.get_virial_alphas(
                centers,
                **self._filter_alpha_method_keys(
                    'virial',
                    opts,
                    ['potential_function', 'gmat_function', 'transformations'],
                    ['scaling']
                )
            )
        elif method == 'masses':
            return self.get_mass_alphas(
                centers,
                **self._filter_alpha_method_keys(
                    'masses',
                    opts,
                    ['masses'],
                    ['scaling']
                )
            )
        # elif method == 'reaction_path':
        #     return self.get_virial_alphas(
        #         centers,
        #         **self._filter_alpha_method_keys(
        #             'virial',
        #             opts,
        #             ['potential_function', 'masses'],
        #             []
        #         )
        #     )
        elif method == 'min_dist':
            return self.get_min_distance_alphas(
                centers,
                **self._filter_alpha_method_keys(
                    'min_dist',
                    opts,
                    ['masses'],
                    ['scaling']
                )
            )
        else:
            raise ValueError("unknown method for getting alphas {}".format(alphas['method']))

    @classmethod
    def get_mass_alphas(cls, centers, *, masses, scaling=10, use_mean=False):
        h_mass =  AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        a = scaling * np.asanyarray(masses) / h_mass
        nats = len(masses)
        ndim = centers.shape[-1] // nats
        a = np.broadcast_to(a[:, np.newaxis], (nats, ndim)).flatten()
        a = np.broadcast_to(a[np.newaxis], centers.shape)
        return a

    @classmethod
    def get_min_distance_alphas(cls, masses, centers, scaling=1/4, use_mean=False):
        # if clustering_radius is None:
        #     clustering_radius = 1
        # np.sqrt(masses / min_dist)
        distances = np.linalg.norm(centers[:, np.newaxis, :] - centers[np.newaxis, :, :], axis=-1)
        # mean_dist = np.average(distances[distances > 1e-8], axis=None)
        distances[distances < 1e-8] = np.max(distances) # exclude zeros
        closest = np.min(distances, axis=1)
        # too hard to compute convex hull for now...so we treat the exterior
        # the same as the interior
        a = scaling * np.sqrt(masses[np.newaxis, :] / closest[:, np.newaxis])
        if use_mean:
            a = np.mean(a)
        # print("???", a[:5])
        return a

    # TODO: add method to rotate to point along bonds
    #       not sure exactly how to define the rest of the coordinate system
    #       but in principle it could be isotropic...
    @classmethod
    def get_virial_alphas(cls,
                          coords,
                          *,
                          potential_function,
                          gmat_function,
                          transformations,
                          scaling=1/2
                          ):
        # we're just computing local frequencies along the coordinates defined by
        # applying the linear transformations given in "transformations" to the supplied coords

        #TODO: avoid recomputing this for like the 15th time
        if isinstance(coords, DGBCoords):
            # raise Exception(coords, coords.centers.shape)
            coords = coords.centers
        f_data = potential_function(coords, deriv_order=2)
        if isinstance(f_data, np.ndarray) and f_data.ndim == 3:
            hess = f_data
            # grad = potential_function(coords, deriv_order=1)
            # pot = potential_function(coords, deriv_order=1)
        else:
            pot, grad, hess = f_data

        gmats = gmat_function(coords)

        # raise Exception(gmats.shape, hess.shape)

        #TODO: support the identity transformation
        # our transformations are expected to be _linear_
        if transformations is not None:
            hess = nput.vec_tensordot(
                nput.vec_tensordot(
                    hess,
                    transformations,
                    axes=[1, 2],
                    shared=1
                ),
                transformations,
                axes=[1, 2],
                shared=1
            )

            gmats = nput.vec_tensordot(
                nput.vec_tensordot(
                    gmats,
                    transformations,
                    axes=[1, 2],
                    shared=1
                ),
                transformations,
                axes=[1, 2],
                shared=1
            )

        # raise Exception(
        #     np.diagonal(gmats, axis1=1, axis2=2),
        #     np.diagonal(hess, axis1=1, axis2=2)
        # )

        alphas = scaling * np.sqrt(np.abs(np.diagonal(hess, axis1=1, axis2=2) / np.diagonal(gmats, axis1=1, axis2=2)))

        return alphas

    @staticmethod
    def _get_hermite_poly(coeff_dict, alphas):
        ndim = alphas.shape[-1]
        # raise Exception(
        #     np.flip(np.asarray(sp.special.hermite(coeff_dict.get(0, 0), monic=False))) *
        #         np.sqrt( (2 * alphas[0]) ** np.arange(coeff_dict.get(0, 0) + 1) )
        # )
        parts = [
            np.flip(np.asarray(sp.special.hermite(coeff_dict.get(k, 0), monic=False)))
            * np.sqrt(
                (2 * alphas[k] ) ** np.arange(coeff_dict.get(k, 0)+1)
                / ( 2**(coeff_dict.get(k, 0)) * np.math.factorial(coeff_dict.get(k, 0)) )
            )
            for k in range(ndim)
        ]
        coeffs = parts[0]
        for c in parts[1:]:
            coeffs = np.multiply.outer(coeffs, c)
        return coeffs
    @classmethod
    def canonicalize_poly_coeffs(cls, coeffs, alphas):
        if coeffs is None:
            return None

        poly_coeffs = [
            cls._get_hermite_poly(c, a) if isinstance(c, dict) else c
            for c, a in zip(coeffs, alphas)
        ]

        return poly_coeffs

    @property
    def transformations(self):
        return self._transforms
    @transformations.setter
    def transformations(self, tf):
        self._transforms = self.canonicalize_transforms(self.alphas, tf)
    @classmethod
    def canonicalize_transforms(self, alphas, tfs):
        npts = alphas.shape[0]
        ndim = alphas.shape[1]
        if tfs is None:
            return tfs
        if isinstance(tfs, np.ndarray):
            dets = np.linalg.det(tfs)
            if not np.allclose(dets, np.ones(len(dets))):
                raise ValueError("inverse must be supplied for non-unitary transforms")
            invs = tfs.transpose((0, 2, 1))
            tfs = (tfs, invs)
        if len(tfs) != 2:
            raise ValueError("need both transformations and inverses")

        tfs, invs = tfs
        if tfs.shape[0] != npts:
            if tfs.shape[0] == 1:
                tfs = np.broadcast_to(tfs, (npts,) + tfs.shape[1:])
            else:
                raise ValueError("wrong number of transformations supplied")
        if invs.shape[0] != npts:
            if invs.shape[0] == 1:
                invs = np.broadcast_to(invs, (npts,) + invs.shape[1:])
            else:
                raise ValueError("wrong number of inverses supplied")
        if tfs.shape[2] != ndim:
            raise ValueError("transformations have wrong dimension got {}[2], expected {}".format(
                tfs.shape,
                ndim
            ))
        if invs.shape[1] != ndim:
            raise ValueError("inverses have wrong dimension got {}[1], expected {}".format(
                invs.shape,
                ndim
            ))

        return (tfs, invs)

    @property
    def S(self):
        if self._S is None:
            self._S = self.get_S()
        return self._S
    @S.setter
    def S(self, smat):
        self._S = smat
    @property
    def T(self):
        if self._T is None:
            self._T = self.get_T()
        return self._T
    @T.setter
    def T(self, tmat):
        self._T = tmat

    def marginalize_out(self, indices):

        subcoords = self.coords.drop_indices(indices)

        ndim = self.coords.shape[-1]
        remaining = np.setdiff1d(np.arange(ndim), indices)

        if self.transformations is not None:
            full_covs = DGBEvaluator.get_inverse_covariances(
                1 / (4*self.alphas), # hacky but the inverse is expected to have 2a on the diagonal
                self.transformations
            )
            full_sel = np.arange(self.alphas.shape[0])

            dropcovs = full_covs[np.ix_(full_sel, indices, indices)]
            dropdets = np.linalg.det(dropcovs)
            zp = np.abs(dropdets) < 1e-3 #TODO: fix this in case someone is in the numerically unstable limit...
            dropdets[zp] = 1
            scaling = 1 / dropdets
            scaling[zp] = 1e3
            scaling *= np.power(np.pi, len(indices) / 4)

            subcovs = full_covs[np.ix_(full_sel, remaining, remaining)]
            inv_alphas, tfs = np.linalg.eigh(subcovs)
            tfs = (tfs, tfs.transpose(0, 2, 1))
            alphas = 1 / (2*inv_alphas) # we have 1/2a
        else:
            scaling = np.power(2 * np.pi, len(indices) / 4) / np.power(np.prod(self.alphas[:, indices], axis=-1), 1 / 4)
            alphas = self.alphas[:, remaining]
            tfs = None

        return scaling, type(self)(
            subcoords,
            alphas=alphas,
            transformations=tfs
        )

    def as_cartesians(self):
        if isinstance(self.coords, DGBCartesians):
            return self

        carts, (transformation, inverse) = self.coords.as_cartesians() # tfs is ncarts x nmodes
        if self.transformations is not None:
            tf, inv = self.transformations
            tfs = transformation[np.newaxis] @ tf
            invs = inv @ inverse[np.newaxis]
        else:
            tfs = np.broadcast_to(transformation[np.newaxis], (carts.shape[0],) + transformation.shape)
            invs = np.broadcast_to(inverse[np.newaxis], (carts.shape[0],) + inverse.shape)
        transformations = (tfs, invs)

        # num_extra = carts.shape[1] - self.alphas.shape[1]
        # alphas = np.concatenate(
        #     [self.alphas, np.zeros((carts.shape[0], num_extra))],
        #     axis=1
        # )
        alphas = self.alphas
        return type(self)(
            carts,
            alphas,
            transformations=transformations
        )

    def plot_centers(self,
                     figure=None,
                     xyz_sel=None,
                     **plot_styles
                     ):
        import McUtils.Plots as plt

        if isinstance(self.coords, DGBCartesians):
            coords = self.coords.coords
            if xyz_sel is None and coords.shape[-1] > 2:
                raise ValueError("need an `xyz_sel` to be plottable")
            coords = coords[:, :, xyz_sel]
            if coords.shape[-1] > 2:
                raise ValueError("can't plot centers for more than 2D data...?")
            for atom in range(coords.shape[1]):
                if coords.shape[-1] == 1:
                    figure = plt.ScatterPlot(
                        coords[:, atom, 0],
                        np.zeros(len(coords)),
                        figure=figure,
                        **plot_styles
                    )
                else:
                    figure = plt.ScatterPlot(
                        coords[:, atom, 0],
                        coords[:, atom, 1],
                        figure=figure,
                        **plot_styles
                    )

        else:
            coords = self.coords.centers
            if coords.shape[-1] > 2:
                raise ValueError("can't plot centers for more than 2D data...?")

            if coords.shape[-1] == 1:
                figure = plt.ScatterPlot(
                    coords[:, 0],
                    np.zeros(len(coords)),
                    figure=figure,
                    **plot_styles
                )
            else:
                figure = plt.ScatterPlot(
                    coords[:, 0],
                    coords[:, 1],
                    figure=figure,
                    **plot_styles
                )

        return figure