import numpy as np, scipy as sp, itertools, functools
from McUtils.Zachary import RBFDInterpolator
from McUtils.Combinatorics import StirlingS1
from McUtils.Scaffolding import Logger

__all__ =  [
    "DGB"
]


class DGB:
    """

    """
    def __init__(self,
                 centers,
                 potential_function,
                 alphas=None,
                 logger=True,
                 clustering_radius=.005,
                 min_singular_value=1e-4,
                 num_svd_points=None,
                 optimize_centers=False,
                 quadrature_degree=4,
                 expansion_degree=None,
                 expansion_type='multicenter'
    ):
        self._S, self._T, self._V = None, None, None
        self.logger = Logger.lookup(logger)

        if optimize_centers:
            self.logger.log_print("optimizing DGB centers...")
            self.clustering_radius, self.centers, self.alphas, self._S, self._T = self.optimize_centers(
                centers, alphas,
                initial_custering=clustering_radius
            )
        else:
            self.logger.log_print("initializing Gaussians...")
            self.clustering_radius = clustering_radius
            self.centers, self.alphas = self.initialize_gaussians(
                centers, alphas,
                clustering_radius
            )

        if min_singular_value is not None or num_svd_points is not None:
            # Use SVD to prune out matrix rows that will be super ill conditioned
            U, sig, VT = np.linalg.svd(self.S)
            if num_svd_points:
                good_loc = np.arange(num_svd_points)
            else:
                self.logger.log_print("most important center threshold: {t}", t=min_singular_value)
                good_loc = np.where(sig > min_singular_value)[0]
            U = U[:, good_loc]
            VT = VT[good_loc, :]

            full_good_pos = np.unique(np.concatenate([
                np.where(np.abs(U) > 1e-14)[0],
                np.where(np.abs(VT) > 1e-14)[1]
            ]))
            self.centers = self.centers[full_good_pos]
            self.alphas = self.alphas[full_good_pos]

            self._S = self._S[full_good_pos, :][:, full_good_pos]
            self._T = self._T[full_good_pos, :][:, full_good_pos]


        self.potential_function = potential_function
        self.quadrature_degree = quadrature_degree
        self.expansion_degree = expansion_degree
        self.expansion_type = expansion_type

    def optimize_centers(self,
                         centers, alphas,
                         max_condition_number=1e16,
                         initial_custering=.005,
                         cluster_step_size=.005,
                         max_steps = 50
                         ):
        c, a = centers, alphas
        cr = initial_custering
        centers, alphas = self.initialize_gaussians(c, a, cr)
        S, T = self.get_ST(centers, alphas)
        n = 0
        while np.linalg.cond(S) > max_condition_number and n < max_steps:
            cr += cluster_step_size
            n += 1
            centers, alphas = self.initialize_gaussians(c, a, cr)
            S, T = self.get_ST(centers, alphas)

        return cr, centers, alphas, S, T

    def initialize_gaussians(self, centers, alphas, clustering_radius):
        centers = np.asanyarray(centers)
        if centers.ndim == 1:
            centers = centers[:, np.newaxis]
        if clustering_radius is not None and clustering_radius >= 0:
            centers, _, _ = RBFDInterpolator.decluster_data(centers, np.empty(len(centers)), [], clustering_radius)

        if alphas is None:
            alphas = self.get_alphas(centers, clustering_radius)
        if isinstance(alphas, (int, float, np.integer, np.floating)):
            alphas = np.full(len(centers), alphas)
        else:
            alphas = np.asanyarray(alphas)

        return centers, alphas

    @classmethod
    def get_alphas(cls, centers, clustering_radius=None):
        if clustering_radius is None:
            clustering_radius = 1
        distances = np.linalg.norm(centers[:, np.newaxis, :] - centers[np.newaxis, :, :], axis=-1)
        mean_dist = np.average(distances[distances > 1e-8], axis=None)
        distances[distances < 1e-8] = np.max(distances) # exclude zeros
        closest = np.min(distances, axis=1)
        # too hard to compute convex hull for now...so we treat the exterior
        # the same as the interior
        a = 1/15*(mean_dist/closest)**2
        print("???", a[:5])
        return a

    @property
    def S(self):
        if self._S is None:
            self._S, self._T = self.get_ST()
        return self._S
    @S.setter
    def S(self, smat):
        self._S = smat
    @property
    def T(self):
        if self._T is None:
            self._S, self._T = self.get_ST()
        return self._T
    @T.setter
    def T(self, tmat):
        self._T = tmat
    @property
    def V(self):
        if self._V is None:
            self._V = self.get_V()
        return self._V
    @V.setter
    def V(self, mat):
        self._V = mat

    def get_ST(self, centers=None, alphas=None):
        if centers is None:
            centers = self.centers
        if alphas is None:
            alphas = self.alphas

        aouter = alphas[:, np.newaxis] * alphas[np.newaxis, :]
        aplus = alphas[:, np.newaxis] + alphas[np.newaxis, :]
        arat = aouter / aplus

        disps = centers[:, np.newaxis, :] - centers[np.newaxis, :, :]

        # A = outer_tet / np.sqrt(np.pi)
        B = np.sqrt(aplus)
        C = arat[:, :, np.newaxis] * np.power(disps, 2)

        # Base components
        S_dim = (np.sqrt(2) * np.power(aouter, 1/4) / B)[:, :, np.newaxis] * np.exp(-C)
        T_dim = arat[:, :, np.newaxis] * (1 - 2*C)

        # Combine appropriately
        S = np.prod(S_dim, axis=-1)
        T = S * np.sum(T_dim, axis=-1)

        return S, T

    def get_overlap_gaussians(self):
        # find overlap gaussians
        new_alphas = self.alphas[:, np.newaxis] + self.alphas[np.newaxis, :]
        w_centers = self.alphas[:, np.newaxis]*self.centers
        # moving weighted average by alpha value
        return (w_centers[:, np.newaxis, :] + w_centers[np.newaxis, :, :])/new_alphas[:, :, np.newaxis], new_alphas

    def quad_integrate(self, degree=2):
        """
        Integrate potential over all pairs of Gaussians at once

        :param degree:
        :type degree:
        :return:
        :rtype:
        """

        # Quadrature point displacements and weights (thanks NumPy!)
        disps, weights = np.polynomial.hermite.hermgauss(degree)

        # I can do only the upper triangle in the future to speed things up
        centers, alphas = self.get_overlap_gaussians()
        pots = np.zeros(alphas.shape)
        ndim = centers.shape[-1]
        for disp_inds in itertools.product(*([range(degree)]*ndim)):
            disp_inds = np.array(disp_inds)
            w = np.prod(weights[disp_inds])
            c = centers + disps[disp_inds][np.newaxis, np.newaxis, :] / np.sqrt(alphas[:, :, np.newaxis])
            pots = pots + w * self.potential_function(c)

        normalization = 1 / (np.sqrt(np.pi)) ** self.centers.shape[-1]
        return pots * normalization

    @classmethod
    def morse_integral1d(cls, centers, alpha, de, a):
        # Centers: (n, n, 2)
        # Alphas: (n, n, 2)
        ...

    # @classmethod
    # def polyint_1D(cls, centers, alphas, order):
    #     ...


    @classmethod
    def polyint_1D(cls, centers, alphas, n):
        if n == 0:
            return np.ones(alphas.shape)
        c = centers*np.sqrt(alphas)
        prefac = cls.polyint_1D_prefac(centers, alphas, n)
        expr = cls.polyint_1D_poly_eval(n, c)
        # with np.printoptions(linewidth=1e8):
        #     print(prefac.shape)
        #     print(expr.shape)

        return prefac * expr

    # a complicated set of functions to handle integrating a monomial (x^n) between
    # two Gaussians in a Taylor series type of way...
    @classmethod
    def polyint_1D_prefac(cls, centers, alphas, n):
        return np.sqrt(np.pi) / np.sqrt(4 * alphas**(n+1))
    _stirs = None
    @classmethod
    def stirling(cls, n, m):
        if cls._stirs is None or len(cls._stirs) <= n:
            cls._stirs = StirlingS1(2*n)
        return cls._stirs[n, m]
    _poly_cache = {}
    @classmethod
    def polyint_1D_poly_coeff_gen(cls, n):
        # cls.stirling(4, 1),
        # with np.printoptions(linewidth=1e8):
        #     print(cls._stirs)
        # raise Exception( ... )
        if n not in cls._poly_cache:
            scaling = np.prod(2*(1+np.arange(n))) * (np.power(2, n-1) if n > 0 else 1/2 )
            terms = [(-1) ** i * cls.stirling(2*n, i) for i in range(1, 2*n+1)]
            cls._poly_cache[n] = (scaling, terms)
        return cls._poly_cache[n]
    @classmethod
    def polyint_1D_coeff(cls, i, n):
        scaling, terms = cls.polyint_1D_poly_coeff_gen(i)
        return np.dot(n**np.arange(1, 2*i+1), terms) / scaling
    @classmethod
    def polyint_1D_coeffs(cls, i):
        o = i % 2
        return [ # I could make this faster but I don't _think_ it'll cost me all that much
                cls.polyint_1D_coeff(
                    (i - (n-1-o) )//2,
                    n
                )
                for n in range(1+o, 1+i, 2)
            ] + [2]
    @classmethod
    def polyint_1D_poly_eval(cls, i, c):
        o = i%2
        exps = np.arange(o, i+1, 2)
        wat = sum(k*c**o for k,o in zip(cls.polyint_1D_coeffs(i), exps))
        return wat
    @classmethod
    def simple_poly_int(cls, alphas, n):
        from scipy.special import gamma
        return (1 + (-1)**n) / np.sqrt(4*alphas**(n+1)) * gamma((n+1)/2)

    def expansion_integrate(self, deriv_order=2, expansion_type=None):
        if expansion_type is None:
            expansion_type = self.expansion_type

        centers, alphas = self.get_overlap_gaussians()
        # this prefac allows us to reexpress things in terms of the new Gaussians
        # with appropriate weighting
        aprod = self.alphas[:, np.newaxis] * self.alphas[np.newaxis, :]
        # asum = self.alphas[:, np.newaxis] + self.alphas[np.newaxis, :]
        # cdiff = self.centers[:, np.newaxis, :] - self.centers[np.newaxis, :, :]
        ndim = centers.shape[-1]
        prefac = ( np.sqrt(2 / np.pi) * (aprod**(1/4)) ) ** (ndim if expansion_type != 'taylor' else 1) #* np.exp(-aprod/asum * np.sum(cdiff**2, axis=-1))

        if expansion_type == 'taylor':
            self.logger.log_print("expanding as a Taylor series about the minumum energy geometry...")
            zero = np.zeros((1, centers.shape[-1]))
            derivs = self.potential_function(zero, deriv_order=deriv_order)
            if isinstance(derivs, np.ndarray): # didn't get the full list so we do the less efficient route
                derivs = (
                        [self.potential_function(zero)] +
                        [self.potential_function(zero, deriv_order=d) for d in range(1, deriv_order + 1)]
                )
            derivs = [
                np.broadcast_to(
                    d[np.newaxis],
                    alphas.shape + d.shape[1:]
                )
                for d in derivs
            ]
        else:

            self.logger.log_print("expanding about {N} points...", N=len(np.triu_indices_from(alphas)[0]))
            # derivs = self.potential_function(centers, deriv_order=deriv_order)
            # if isinstance(derivs, np.ndarray):  # didn't get the full list so we do the less efficient route'
            #     derivs = [self.potential_function(centers)] + [self.potential_function(centers, deriv_order=d) for d in
            #                                                range(1, deriv_order + 1)]

            # use symmetry to cut down number of indices by 2
            row_inds, col_inds = np.triu_indices_from(alphas)
            ics = centers[row_inds, col_inds]
            derivs = self.potential_function(ics, deriv_order=deriv_order)
            if isinstance(derivs, np.ndarray):  # didn't get the full list so we do the less efficient route'
                derivs = [self.potential_function(ics)] + [
                    self.potential_function(ics, deriv_order=d) for d in range(1, deriv_order+1)
                ]
            new_derivs = [] # reverse symmetrization
            for d in derivs:
                full_mat = np.empty(
                    alphas.shape + d.shape[1:]
                )
                full_mat[row_inds, col_inds] = d
                full_mat[col_inds, row_inds] = d
                new_derivs.append(full_mat)
            derivs = new_derivs

        caches = [{} for _ in range(ndim)]
        pot = 0
        for d in derivs: # add up all independent integral contribs...
            # iterate over upper triangle coordinates (we'll add bottom contrib by symmetry)
            inds = itertools.combinations_with_replacement(range(ndim), r=d.ndim-2) if d.ndim > 2 else [()]
            for idx in inds:
                count_map = {k: v for k, v in zip(*np.unique(idx, return_counts=True))}
                contrib = 1
                for k in range(ndim): # do each dimension of integral independently
                    n = count_map.get(k, 0)
                    if n not in caches[k]:
                        caches[k][n] = (
                           self.simple_poly_int(alphas, n)
                            if expansion_type != 'taylor' else
                           self.polyint_1D(centers[..., k], alphas, n)
                        )

                    contrib *= caches[k][n]

                    # with np.printoptions(linewidth=1e8):
                    #     print(n, centers[0, 0], alphas[0, 0])
                    #     print(self.centers[0])
                    #     print(contrib[0, 0])

                dcont = d[(slice(None, None, None), slice(None, None, None)) + idx] if len(idx) > 0 else d
                facterms = np.unique([x for x in itertools.permutations(idx)], axis=0)
                nfac = len(facterms) # this is like a binomial coeff or something but my sick brain won't work right now...
                scaling = np.prod([np.math.factorial(count_map.get(k, 0)) for k in range(ndim)])

                contrib = contrib * nfac * dcont / scaling

                pot += contrib
        return pot * prefac


    def analytic_integrate(self):
        raise NotImplementedError("flooped up")
        centers = [np.array(np.meshgrid(x, x)).T for x in self.centers.T]
        alphas = np.array(np.meshgrid(self.alphas, self.alphas)).T
        # raise Exception(alphas.shape)
        return self.potential_function['analytic_integrals'](
            centers, # ...no
            alphas
        )

    def get_base_pot(self, potential_handler=None, expansion_degree=None, degree=None):
        if expansion_degree is None:
            expansion_degree = self.expansion_degree

        if potential_handler is None:
            if isinstance(self.potential_function, dict):
                if 'analytic_integrals' in self.potential_function:
                    potential_handler = 'analytic'
            elif expansion_degree is not None:
                potential_handler = 'expansion'
            else:
                potential_handler = 'quad'

        if potential_handler == 'quad':
            self.logger.log_print("evauating integrals with {n}-order quadrature", n=self.quadrature_degree)
            pot_mat = self.quad_integrate(degree=self.quadrature_degree if degree is None else degree)
        elif potential_handler == 'expansion':
            self.logger.log_print("evauating integrals with {n}-degree expansions", n=self.expansion_degree)
            pot_mat = self.expansion_integrate(deriv_order=expansion_degree)
        elif potential_handler == 'analytic':
            self.logger.log_print("evauating integrals analytically", n=self.expansion_degree)
            pot_mat = self.analytic_integrate()
        else:
            raise ValueError("woof")

        return pot_mat

    def get_V(self, potential_handler=None, expansion_degree=None, degree=None):

        self.logger.log_print("calculating potential matrix")
        pot_mat = self.get_base_pot(potential_handler=potential_handler, expansion_degree=expansion_degree, degree=degree)
        return self.S * pot_mat

    def get_wavefunctions(self, print_debug_info=False, min_singular_value=1e-4):
        H = self.T + self.V
        if print_debug_info:
            print('Condition Number:', np.linalg.cond(self.S))
            with np.printoptions(linewidth=1e8, threshold=1e8):
            #     # print(self.S)
            #     # print(self.T)
                print("Potential Matrix:")
                print(self.V)

        if min_singular_value is None:
            try:
                return sp.linalg.eigh(H, self.S)
            except np.linalg.LinAlgError:
                raise ValueError(
                    "Overlap matrix poorly conditioned ({}) usually means data is too clustered or a min singular value is required".format(np.linalg.cond(self.S))
                ) from None
        else:
            S = self.S

            # Use SVD to prune out matrix rows that will be super ill conditioned
            U, sig, VT = np.linalg.svd(S)
            good_loc = np.where(sig > min_singular_value)[0]
            d = np.diag(1/np.sqrt(sig[good_loc]))
            U = U[:, good_loc]
            VT = VT[good_loc, :]
            Q = U @ d @ VT
            Qinv = VT.T @ d @ U.T

            # shift_val = 1e5
            # shift = np.diag(np.full(len(H), shift_val)) # add a uniform shift to make it clear where the zeros are
            # H = Q @ (H + shift) @ Q

            H = Q @ H @ Q

            eigs, evals = np.linalg.eigh(H)

            # take non-zero eigvals
            zero_pos = np.where(np.abs(eigs) < 1e-8)[0] # drop the appropriate number of zeros
            nz = len(H) - len(good_loc)
            zero_pos = zero_pos[:nz]
            sel = np.setdiff1d(np.arange(len(H)), zero_pos)
            eigs = eigs[sel]
            evals = evals[:, sel]

            # now invert transformation to get back to OG basis
            # and remove uniform shift to eigenvalues
            eigs = eigs # - shift_val
            evals = Qinv @ evals

            return eigs, evals


