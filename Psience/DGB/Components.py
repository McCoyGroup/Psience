
import abc

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
        Transforms the alphas into proper inverse covariance matrices

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
        new_derivs = []
        # rotations = rotations[:, :, :, np.newaxis] # to test shapes
        for n,d in enumerate(derivs):
            # print("???", d.shape)
            for _ in range(n):
                d = nput.vec_tensordot(
                    d, rotations,
                    axes=[1, 1],
                    shared=1
                )
            new_derivs.append(d)
        derivs = new_derivs

        if logger is not None:
            logger.log_print("adding up all derivative contributions...")

        row_inds = overlap_data['row_inds']
        col_inds = overlap_data['col_inds']

        fdim = derivs[0].ndim - 1
        fshape = derivs[0].shape[1:]
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

                dcont = d[(slice(None, None, None),) + idx] if len(idx) > 0 else d
                facterms = np.unique([x for x in itertools.permutations(idx)], axis=0)
                nfac = len(facterms) # this is like a binomial coeff or something but my sick brain won't work right now...
                scaling = 2**(len(idx)/2) * np.prod([np.math.factorial(count_map.get(k, 0)) for k in range(ndim)])
                for _ in range(fdim):
                    scaling = np.expand_dims(scaling, -1)

                if reweight:
                    contrib *= nfac * dcont / scaling
                else:
                    contrib *= nfac * dcont

                pot[row_inds, col_inds] += contrib

        pot[col_inds, row_inds] = pot[row_inds, col_inds]

        return pot

    @staticmethod
    def quad_nd(centers, alphas, function, flatten=False, degree=3, normalize=True):
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

        w = np.prod(weights[indices], axis=-1)
        d = disps[indices][np.newaxis] / np.sqrt(alphas)[:, np.newaxis]
        c = centers[:, np.newaxis, :] + d
        fv = function(c.reshape(-1, ndim)).reshape(c.shape[:2])
        val = np.sum(w[np.newaxis, :] * fv, axis=1)

        if normalize:
            normalization = 1 / np.prod(np.sqrt(alphas), axis=-1)
            val = val * normalization
        return val

    @classmethod
    def _wrap_rotated_function(cls, func, overlap_data):
        rotations = overlap_data['rotations']
        inverse = overlap_data['inverse_rotations']
        # centers = np.reshape(rotations @ centers[:, :, np.newaxis], centers.shape)

        # if inds is not None:
        #     ndim = len(inds)
        #     rotations = rotations[:, :, inds]
        #     alphas = alphas[..., inds]
        #     centers = centers[..., inds]
        #
        # # if self.
        # new_derivs = []
        # # rotations = rotations[:, :, :, np.newaxis] # to test shapes
        # for n, d in enumerate(derivs):
        #     # print("???", d.shape)
        #     for _ in range(n):
        #         d = nput.vec_tensordot(
        #             d, rotations,
        #             axes=[1, 1],
        #             shared=1
        #         )
        #     new_derivs.append(d)
        # derivs = new_derivs

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
    def quad_integrate(cls, function, overlap_data, degree=2, logger=None):
        """
        Integrate potential over all pairs of Gaussians at once

        :param degree:
        :type degree:
        :return:
        :rtype:
        """

        f = cls._wrap_rotated_function(function, overlap_data)
        unrot_centers = overlap_data['centers']
        alphas = overlap_data['alphas']
        # expressed along ellipsoid axes
        centers = np.reshape(overlap_data['rotations'] @ unrot_centers[:, :, np.newaxis], unrot_centers.shape)

        npts = overlap_data['init_centers'].shape[0]
        pots = np.zeros((npts, npts))
        rows, cols = np.triu_indices(npts)

        vals = cls.quad_nd(centers, alphas, f, degree=degree, normalize=False)
        pots[rows, cols] = vals
        pots[cols, rows] = vals

        normalization = 1 / (np.sqrt(np.pi)) ** centers.shape[-1] # part of the prefactor...
        pots *= normalization

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

        # ndim = centers.shape[-1]
        poly = cls._induced_polys(terms, init_sigs[cols], init_cents[cols])
        pcs = poly.coeffs
        nz_pos = np.where(prefactors != 0)
        # force prefactors to be broadcastable
        scaled_coeffs = np.expand_dims(
            prefactors[nz_pos],
            list(range(1, 1 + poly.coordinate_dim))
        ) * pcs[nz_pos]
        poly = DensePolynomial(scaled_coeffs, stack_dim=1)
        # print(poly.coeffs[0])
        poly = poly.shift((-centers)[nz_pos[0]])
        # print(poly.coeffs[0])
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
    def _induced_polys(cls, terms, sigs, centers):
        sig_polys = DensePolynomial.from_tensors(
            [np.zeros(sigs.shape[0]), np.zeros(sigs.shape[:2]), -sigs/2]
        ).shift(centers)
        sig_grad = sig_polys.grad() # for applying p
        monomial = DensePolynomial.from_tensors( # for applying q
            [np.zeros(sig_polys.coordinate_dim), np.eye(sig_polys.coordinate_dim)]
        ) # stack_dim = 1
        poly = DensePolynomial(
            np.ones((len(sigs),) + (1,) * sig_polys.coordinate_dim),
            stack_dim=1
        ) # stack_dim = 1
        for t in terms:
            # _ = []
            if t == 'q':
                mon_coeffs = np.expand_dims(monomial.coeffs, list(range(poly.stack_dim)))
                poly_coeffs = np.expand_dims(poly.coeffs, poly.stack_dim)
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
            elif t == 'p':
                # first term is just taking derivs of coeffs
                # basic idea is each time we hit a p we have
                # d/dq_i(K * e^(sigma*q*q)) = (d/dq_i(K) + K * ??? )e^(sigma*q*q)
                poly_1 = poly.grad()
                sigg_coeffs = np.expand_dims(sig_grad.coeffs, list(range(1, poly.stack_dim)))
                poly_coeffs = np.expand_dims(poly.coeffs, poly.stack_dim)
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

class DGBCartesianEvaluator(DGBKineticEnergyEvaluator):

    def __init__(self, masses):
        self.masses = masses

    def evaluate_ke(self, overlap_data):
        masses = self.masses
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
    def _embed_watson_modes(cls, watson_data, centers):
        origin = watson_data['origin'].reshape(-1)
        tf = watson_data['matrix']
        new_cent_vecs = (centers-origin[np.newaxis]).reshape(centers.shape[0], -1, centers.shape[1])@tf.T[np.newaxis]
        return new_cent_vecs.reshape(-1, new_cent_vecs.shape[-1]) # basically a squeeze

    # def evaluate_ke(cls, overlap_data):
    #     raise NotImplementedError(...)
    #     prefactors = np.broadcast_to(
    #         np.diagonal(1 / masses)[np.newaxis],
    #         (overlap_data['centers'].shape[0], len(masses), len(masses))
    #     )
    #     return cls._evaluate_polynomial_ke(overlap_data, ['p', 'p'], prefactors) / 2

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
                 pairwise_function_options=None,
                 logger=None
                 ):
        self.potential_function = potential_function
        self.handler = integral_handler
        self.expansion_degree = expansion_degree
        self.expansion_type = expansion_type
        self.quadrature_degree = quadrature_degree
        self.pairwise_handler = None if pairwise_functions is None else DGBPairwisePotentialEvaluator(
            pairwise_functions,
            **({} if pairwise_function_options is None else pairwise_function_options)
        )

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
            raise NotImplementedError("not done delegating")
            pot_contribs, deriv_corrs = ...
        else:
            pot_contribs = None
            deriv_corrs = None

        if not isinstance(logger, Logger):
            logger = Logger.lookup(logger)

        deriv_order = expansion_degree # renamed for consistency
        if expansion_type == 'taylor':
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

        npts = len(overlap_data['init_centers'])
        pot = cls.tensor_expansion_integrate(
            npts,
            derivs,
            overlap_data,
            expansion_type=expansion_type,
            logger=logger
        )

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

class DGBPairwisePotentialEvaluator(DGBEvaluator):

    def __init__(self,
                 pairwise_potential_functions,
                 coord_shape,
                 quadrature_degree=4,
                 expansion_degree=4
                 ):
        self.pairwise_potential_functions = pairwise_potential_functions
        self.coord_shape = coord_shape
    @staticmethod
    def _wrap_pairwise(func):
        def eval(coords, deriv_order=None):
            if deriv_order is None:
                r = np.linalg.norm(coords, axis=-1) # we ensure no scaling of difference coords
                r_derivs = None
            else:
                r_derivs = nput.vec_norm_derivs(coords, order=deriv_order)
                r, r_derivs = r_derivs[0], r_derivs[1:]
            fdat = func(r, deriv_order=deriv_order)
            if r_derivs is not None:
                f_vals = fdat[0]
                r_derivs = [r[..., np.newaxis] for r in r_derivs]
                f_derivs = [fdat[i].reshape(f_vals.shape + (1,)*i) for i in range(1, deriv_order+1)]
                f_derivs = list(TensorDerivativeConverter(r_derivs, f_derivs).convert())
                fdat = [f_vals] + f_derivs
            return fdat

        return eval

    @classmethod
    def integrate_pairwise_potential_contrib(cls,
                                             functions,
                                             centers,
                                             alphas,
                                             coord_shape=None,
                                             quadrature_degree=4,
                                             expansion_degree=4
                                             ):

        if coord_shape is not None and len(coord_shape) == 1:
            raise ValueError("can't determine if coords are planar or not")

        if coord_shape is None:
            raise NotImplementedError("coord shape inference not implemented for pairwise potentials")
        else:
            ndim = coord_shape[1]

        npts = centers.shape[0]
        fdim = centers.shape[1]
        centers = centers.reshape(npts, -1, ndim)
        alphas = alphas.reshape(npts, -1, ndim)

        potential_contrib = 0
        if expansion_degree is not None:
            derivative_contribs = [
                np.zeros((npts,) + (fdim,)*i)
                    if i > 0 else
                np.zeros((npts,))
                for i in range(expansion_degree+1)
            ]
        else:
            derivative_contribs = [np.zeros((npts,))]
        for index_pair, pairwise_func in functions.items():
            f = cls._wrap_pairwise(pairwise_func)
            subcenters = centers[:, index_pair, :].reshape(npts, 2*ndim)
            subalphas = alphas[:, index_pair, :].reshape(npts, 2, ndim)

            # if ndim == 1:
            #     tf = np.array([
            #         [1, -1],
            #         [1,  1]
            #     ]) # / np.sqrt(2)
            #     tf = np.broadcast_to(tf[np.newaxis], (npts, 2, 2)).copy()
            #     tf[:, 1, 1] = subalphas[:, 0, 0] / subalphas[:, 1, 0]
            # elif ndim == 2:
            #     tf = np.array([
            #         [1, 0, -1,  0],
            #         [0, 1,  0, -1],
            #         [1, 0,  1,  0],
            #         [0, 1,  0,  1]
            #     ]) #/ np.sqrt(2)
            #     tf = np.broadcast_to(tf[np.newaxis], (npts, 4, 4)).copy()
            #     tf[:, 2, 2] = subalphas[:, 0, 0] / subalphas[:, 1, 0]
            #     tf[:, 3, 3] = subalphas[:, 0, 1] / subalphas[:, 1, 1]
            # elif ndim == 3:
            #     tf = np.array([
            #         [1, 0, 0, -1,  0,  0],
            #         [0, 1, 0,  0, -1,  0],
            #         [0, 0, 1,  0,  0, -1],
            #         [1, 0, 0,  1,  0,  0],
            #         [0, 1, 0,  0,  1,  0],
            #         [0, 0, 1,  0,  0,  1]
            #     ]) #/ np.sqrt(2)
            #     tf = np.broadcast_to(tf[np.newaxis], (npts, 6, 6)).copy()
            #     tf[:, 3, 3] = subalphas[:, 0, 0] / subalphas[:, 1, 0]
            #     tf[:, 4, 4] = subalphas[:, 0, 1] / subalphas[:, 1, 1]
            #     tf[:, 5, 5] = subalphas[:, 0, 2] / subalphas[:, 1, 2]
            # else: # implements the above transformations for arbitary dims

            tf = np.eye(2*ndim)
            off_diag_rows = np.arange(ndim)
            off_diag_cols = ndim + off_diag_rows
            tf[off_diag_rows, off_diag_cols] = -1
            tf[off_diag_cols, off_diag_rows] =  1
            tf = np.broadcast_to(tf[np.newaxis], (npts, 2*ndim, 2*ndim)).copy()
            for i in range(ndim):
                j = ndim + i
                tf[:, j, j] = subalphas[:, 0, i] / subalphas[:, 1, i]

            # tf_inv = np.eye(2*ndim)
            # tf_inv = np.broadcast_to(tf_inv[np.newaxis], (npts, 2 * ndim, 2 * ndim)).copy()
            # for i in range(ndim):
            #     j = ndim + i
            #     a = 1 + subalphas[:, 0, i] / subalphas[:, 1, i]
            #     tf_inv[:, i, i] = (a - 1) / a
            #     tf_inv[:, j, j] =  1 / a
            #     tf_inv[:, i, j] =  1 / a
            #     tf_inv[:, j, i] = -1 / a
            # tf_dets = np.linalg.det(tf)
            # print(tf[0])
            # print(tf_inv[0])
            # print(tf_dets)

            subcov = np.zeros((subalphas.shape[0], 2*ndim, 2*ndim))
            # fill diagonals across subcov
            diag_inds = (slice(None, None, None),) + np.diag_indices(2*ndim)
            subcov[diag_inds] = np.concatenate([subalphas[:, 0, :], subalphas[:, 1, :]], axis=1)
            tf_centers = (subcenters[:, np.newaxis, :] @ tf.transpose(0, 2, 1) ).reshape(-1, 2*ndim)
            tf_cov = tf@subcov@tf.transpose(0, 2, 1)
            # print(tf_cov[0])
            # print("...", tf_centers[:, 0])
            tf_alphas = tf_cov[diag_inds] / 4 # account for sqrt of 2 I omitted
            # print(tf_centers[:, :ndim])
            pairwise_contrib = cls.quad_nd(
                tf_centers[:, :ndim],
                tf_alphas[:, :ndim],
                f,
                degree=18,
                normalize=True
            ) * (
                np.prod(np.sqrt(np.pi / tf_alphas[:, ndim:]), axis=-1)
            ) / 2

            # now we scale by dgb prefactors
            prefactor_scaling = np.power(
                np.prod(subalphas.reshape(-1, 2*ndim), axis=-1),
                1/2
            ) * (1/np.pi)**ndim
            # print(prefactor_scaling)
            # print(
            #     np.prod(subalphas.reshape(-1, 2*ndim), axis=-1)
            # )

            potential_contrib += prefactor_scaling * pairwise_contrib
            # print(potential_contrib)
            # raise Exception(potential_contrib)

            # tf_inv = tf_inv * tf_dets[:, np.newaxis, np.newaxis]
            if expansion_degree is not None:
                # need to map index_pair (which correspond to 4 or 6 inds) to the
                # correct list of inds in full space
                idx_pair_full_inds = sum(
                    (tuple(range(ndim*x, ndim*(x+1))) for x in index_pair),
                    ()
                )
                tf_derivs_contribs = f(tf_centers[:, :ndim], deriv_order=expansion_degree)
                for i,d in enumerate(tf_derivs_contribs):
                    if i == 0:
                        derivative_contribs[0] += d
                    else:
                        subd_contrib = np.zeros((npts,) + (2*ndim,) * i)
                        subd_inds = (slice(None, None, None),) + (slice(None, ndim, None),) * i # contrib[:, :ndim, :ndim, etc]
                        subd_contrib[subd_inds] = d
                        # print(subd_contrib.shape, tf_inv.shape)
                        for j in range(i):
                            # print(subd_contrib.shape[j+1], tf_inv.shape[1])
                            subd_contrib = nput.vec_tensordot(subd_contrib, tf, shared=1, axes=[1, 1])
                        dc_inds = (slice(None, None, None),) + np.ix_(*[idx_pair_full_inds,]*i)
                        # print(dc_inds)
                        # print(derivative_contribs[i].shape, derivative_contribs[i][dc_inds].shape)
                        # print(subd_contrib.shape)
                        derivative_contribs[i][dc_inds] += subd_contrib
                # raise Exception(...)
            else:
                derivative_contribs[0] += f(tf_centers[:, :ndim])

        ## factor out prefactors from integrals
        # normalization = 1 / (np.sqrt(np.pi)) ** ndim
        # pairwise_contrib *= normalization
        # "1.19662817"
        # print("--->", pairwise_contrib[0], np.linalg.det(tf))
        # raise Exception(potential_contrib[:2])

        return potential_contrib, derivative_contribs

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
    @abc.abstractmethod
    def __getitem__(self, item) -> 'DGBCoords':
        ...
class DGBCartesians(DGBCoords):
    def __init__(self, coords, masses):
        coords = np.asanyarray(coords)
        masses = np.asanyarray(masses)
        if masses.ndim > 1:
            raise ValueError("expected a vector of masses...")
        coords = coords.reshape((len(coords), len(masses), -1))

        self.coords = coords
        self.masses = masses
    @property
    def centers(self):
        return self.coords.reshape((self.coords.shape[0], self.coords.shape[1]*self.coords.shape[2]))
    @property
    def kinetic_energy_evaluator(self):
        return DGBCartesianEvaluator(
            np.broadcast_to(self.masses, self.centers.shape[1:]).flatten()
        )
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

    def __getitem__(self, item):

        c = self.centers[item]
        if not isinstance(c, np.ndarray) or c.ndim != 3:
            raise ValueError("bad slice {}".format(item))

        if len(item) == 1:
            return type(self)(
                self.centers[item],
                self.masses
            )
        elif len(item) > 1:
            return type(self)(
                self.centers[item],
                self.masses[item[1]]
            )
        else:
            raise ValueError("zero-length slice?")
class DGBInternals(DGBCoords):
    def __init__(self, coords, gmat_function=None, vprime_function=None):
        raise NotImplementedError('internals coming soon?')
class DGBWatsonModes(DGBCoords):
    def __init__(self, coords, modes, coriolis_inertia_function=None):
        self.coords = coords
        self.modes = modes
        self.ci_func = coriolis_inertia_function

    @property
    def centers(self):
        return self.coords
    @property
    def kinetic_energy_evaluator(self):
        return DGBWatsonEvaluator(self.modes, self.ci_func)

    @classmethod
    def default_coriolis_intertia_function(cls, modes, masses):
        def coriolis_inertia_function(watson_coords):
            carts = watson_coords[:, np.newaxis, :] @ modes.inverse.transpose()[np.newaxis]
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

            B_e, eigs = StructuralProperties.get_prop_moments_of_inertia(carts, masses)
            J = modes.matrix.reshape()
            J = np.tensordot(J, eigs, axes=[2, 0])

            # coriolis terms are given by zeta = sum(JeJ^T, n)
            # which is also just zeta[x][r, s] = sum(J[r][n, y]*J[s][n, z] - J[r][n, z]*J[s][n, y], n)
            ce = -nput.levi_cevita3
            zeta = sum(
                np.tensordot(
                    np.tensordot(ce, J[n], axes=[0, 1]),
                    J[n],
                    axes=[1, 1]
                )
                for n in range(J.shape[0])
            )

            return zeta, B_e

        return coriolis_inertia_function

    @classmethod
    def from_cartesians(cls,
                        coords,
                        modes,
                        masses=None,
                        coriolis_inertia_function=None
                        ):
        flat_carts = (coords - modes.origin[np.newaxis]).reshape((len(coords), -1))
        coords = flat_carts[:, np.newaxis, :] @ modes.matrix[np.newaxis]
        if coriolis_inertia_function is None:
            coriolis_inertia_function = cls.default_coriolis_intertia_function(modes, masses)
        return cls(coords, modes, coriolis_inertia_function=coriolis_inertia_function)

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

    def __init__(self, coords, alphas, transformations, poly_coeffs=None, logger=None):
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
        # raise Exception(self._poly_coeffs[0].shape, base_coeffs[0].shape)
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

    def optimize(self, min_singular_value=None, num_svd_vectors=None, svd_contrib_cutoff=1e-3):
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

        centers = self.coords[full_good_pos]
        alphas = self.alphas[full_good_pos]
        if self.transformations is not None:
            transformations = [t[full_good_pos] for t in self.transformations]
        else:
            transformations = None

        if self._poly_coeffs is not None:
            _poly_coeffs = [t[full_good_pos] for t in self._poly_coeffs]
        else:
            _poly_coeffs = None

        return type(self)(
            centers,
            alphas,
            transformations,
            poly_coeffs=_poly_coeffs
        )

    def optimize_centers(self,
                         centers, alphas,
                         max_condition_number=1e16,
                         initial_custering=.005,
                         cluster_step_size=.005,
                         max_steps=50
                         ):
        raise NotImplementedError("ugh...")
        c, a = centers, alphas
        cr = initial_custering
        centers, alphas = self.initialize_gaussians(c, a, cr)
        _, S, T = self.get_ST(centers, alphas)
        n = 0
        while np.linalg.cond(S) > max_condition_number and n < max_steps:
            cr += cluster_step_size
            n += 1
            centers, alphas = self.initialize_gaussians(c, a, cr)
            S, T = self.get_ST(centers, alphas)

        return cr, centers, alphas, S, T

    def decluster_gaussians(self, centers, alphas, clustering_radius):
        raise NotImplementedError("need to ensure this works right")
        centers = np.asanyarray(centers)
        if centers.ndim == 1:
            centers = centers[:, np.newaxis]

        mask = None
        if clustering_radius is not None and clustering_radius >= 0:
            centers, _, _, mask = RBFDInterpolator.decluster_data(centers, np.empty(len(centers)), [], clustering_radius, return_mask=True)

        rots = None
        if alphas is None:
            alphas = self.get_alphas(self.masses, centers)#, clustering_radius)
        elif isinstance(alphas, (str, dict)):
            alphas, rots, opts = self.dispatch_get_alphas(alphas, centers)
            for k,v in opts.items():
                setattr(self, k, v)

        if isinstance(alphas, (int, float, np.integer, np.floating)):
            alphas = np.full(centers.shape, alphas)
        else:
            alphas = np.asanyarray(alphas)
            if mask is not None:
                alphas = alphas[mask]

        return centers, alphas, rots

    @classmethod
    def _embedded_function(cls, func, modes, natoms=None):
        @functools.wraps(func)
        def embedded_function(mode_coords, deriv_order=None):
            carts = mode_coords[:, np.newaxis, :] @ modes.inverse.transpose()[np.newaxis] + modes.origin[np.newaxis]
            if natoms is not None:
                carts = carts.reshape((carts.shape[0], natoms, -1))

            vals = func(carts, deriv_order=deriv_order)
            if deriv_order is not None:
                _ = []
                for n,d in enumerate(vals):
                    for j in range(n):
                        d = np.tensordot(d, modes.matrix, axes=[1, 1])
                    _.append(d)
                vals = _
            return vals

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

                        potential_function = cls._embedded_function( # a function that takes in normal mode coordinates
                            potential_function,
                            modes,
                            natoms=coords.shape[-1] if coords.ndim == 3 else None
                        )

                    else:
                        raise ValueError("unknown mode spec {}".format(modes))

                if coords.ndim == 3: # Cartesian-like
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
                    coords = DGBWatsonModes(
                        modes.matrix(coords),
                        modes=modes
                    )

            elif internals is not None:
                raise NotImplementedError("internal coordinate mapping not supported yet...")
            else:
                coords = DGBCartesians.from_cartesians(coords, masses=masses, atoms=atoms)

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
            alphas = cls.dispatch_get_alphas(
                alphas,
                coords,
                potential_function=potential_function,
                transformations=transformations
            )

        return cls(coords, alphas, transformations,
                   poly_coeffs=poly_coeffs, logger=logger
                   )
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
            if len(stationary_structures) == 0:
                raise ValueError("no stationary state found?")
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
    def dispatch_get_alphas(self, alphas, centers):
        if isinstance(alphas, str):
            alphas = {'method':alphas}
        opts = alphas.copy()
        del opts['method']
        method = alphas['method']
        if method == 'virial':
            return self.get_virial_alphas(
                centers,
                **self._filter_alpha_method_keys(
                    'virial',
                    opts,
                    ['potential_function'],
                    ['scaling']
                )
            )
        elif method == 'reaction_path':
            return self.get_virial_alphas(
                centers,
                **self._filter_alpha_method_keys(
                    'virial',
                    opts,
                    ['potential_function', 'masses']
                )
            )
        elif method == 'min_dist':
            return self.get_alphas(
                centers,
                **self._filter_alpha_method_keys(
                    'min_dist',
                    opts,
                    ['masses']
                )
            )
        else:
            raise ValueError("unknown method for getting alphas {}".format(alphas['method']))

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
                          potential_function,
                          gmat_function,
                          transformations,
                          scaling=1
                          ):
        # we're just computing local frequencies along the coordinates defined by
        # applying the linear transformations given in "transformations" to the supplied coords

        #TODO: avoid recomputing this for like the 15th time
        f_data = potential_function(coords, deriv_order=2)
        if isinstance(f_data, np.ndarray) and f_data.ndim == 3:
            hess = f_data
            # grad = potential_function(coords, deriv_order=1)
            # pot = potential_function(coords, deriv_order=1)
        else:
            pot, grad, hess = f_data

        gmats = gmat_function(coords)

        #TODO: support the identity transformation
        # our transformations are expected to be _linear_
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

        alphas = scaling * np.sqrt(np.abs(np.diagonal(gmats, axis1=1, axis2=2) * np.diagonal(hess, axis1=1, axis2=2)))

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
        if (
                len(tfs) != 2 or
                not all(
                    isinstance(x, np.ndarray) and x.shape == (npts, ndim, ndim)
                    for x in tfs
                )
        ):
            raise ValueError("transforms must have shape {s}".format(
                s=(npts, ndim, ndim)
            ))
        return tfs

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