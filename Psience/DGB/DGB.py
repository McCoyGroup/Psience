import numpy as np, scipy as sp, itertools, functools
from McUtils.Zachary import RBFDInterpolator
from McUtils.Combinatorics import StirlingS1

__all__ =  [
    "DGB"
]


class DGB:
    """

    """
    def __init__(self,
                 centers,
                 potential_function,
                 alphas=1,
                 clustering_radius=.075,
                 quadrature_degree=4,
                 expansion_degree=None,
                 expansion_type='multicenter'
    ):
        centers = np.asanyarray(centers)
        if centers.ndim == 1:
            centers = centers[:, np.newaxis]
        if clustering_radius is not None and clustering_radius >= 0:
            centers, _, _ = RBFDInterpolator.decluster_data(centers, np.empty(len(centers)), [], clustering_radius)
        self.centers = centers

        if isinstance(alphas, (int, float, np.integer, np.floating)):
            alphas = np.full(len(centers), alphas)
        self.alphas = np.asanyarray(alphas)
        self.potential_function = potential_function
        self.quadrature_degree = quadrature_degree
        self.expansion_degree = expansion_degree
        self.expansion_type = expansion_type
        self._S, self._T, self._V = None, None, None

    @property
    def S(self):
        if self._S is None:
            self._S, self._T = self.get_ST()
        return self._S
    @property
    def T(self):
        if self._T is None:
            self._S, self._T = self.get_ST()
        return self._T
    @property
    def V(self):
        if self._V is None:
            self._V = self.get_V()
        return self._V

    def get_ST(self):

        aouter = self.alphas[:, np.newaxis] * self.alphas[np.newaxis, :]
        aplus = self.alphas[:, np.newaxis] + self.alphas[np.newaxis, :]
        arat = aouter / aplus

        disps = self.centers[:, np.newaxis, :] - self.centers[np.newaxis, :, :]

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
            return np.ones(alphas.shape[:-1])
        c = centers*np.sqrt(alphas)
        prefac = cls.polyint_1D_prefac(centers, alphas, n)
        expr = cls.polyint_1D_poly_eval(n, c)
        # with np.printoptions(linewidth=1e8):
        #     print(prefac)
        #     print(expr)

        return prefac * expr

    # a complicated set of functions to handle integrating a monomial (x^n) between
    # two Gaussians in a Taylor series type of way...
    @classmethod
    def polyint_1D_prefac(cls, centers, alphas, n):
        return np.sqrt( np.pi / (4 * alphas**(n+1)) )
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
        return np.dot(
            c[:, :, np.newaxis]**np.arange(o, i+1, 2)[np.newaxis, :],
            cls.polyint_1D_coeffs(i)
        )
    @classmethod
    def simple_poly_int(cls, alphas, n):
        """
        1/2 (1 + (-1)^n) \[Alpha]^(-(1/2) - n/2) Gamma[(1 + n)/2]
        """
        from scipy.special import gamma
        return 1/2 * (1 + (-1)**n) / np.sqrt(alphas**(n+1)) * gamma((n+1)/2)

    def expansion_integrate(self, deriv_order=2, expansion_type='multicenter'):
        centers, alphas = self.get_overlap_gaussians()
        # this prefac allows us to reexpress things in terms of the new Gaussians
        # with appropriate weighting
        aprod = self.alphas[:, np.newaxis] * self.alphas[np.newaxis, :]
        asum = self.alphas[:, np.newaxis] + self.alphas[np.newaxis, :]
        cdiff = self.centers[:, np.newaxis, :] - self.centers[np.newaxis, :, :]
        prefac = 4 * np.sqrt(2) / np.sqrt(np.pi) * (aprod**1/4) * np.exp(-aprod/asum * np.sum(cdiff**2, axis=-1))
        ndim = centers.shape[-1]
        if expansion_type == 'taylor':
            zero = np.zeros((1, centers.shape[-1]))
            derivs = (
                    [self.potential_function(zero)] +
                     [self.potential_function(zero, deriv_order=d) for d in range(1, deriv_order+1)]
            )
            derivs = [
                np.broadcast_to(
                    d[np.newaxis],
                    alphas.shape + d.shape[1:]
                )
                for d in derivs
            ]
        else:
            derivs = [self.potential_function(centers)] + [self.potential_function(centers, deriv_order=d) for d in range(1, deriv_order+1)]
        caches = [{} for _ in range(ndim)]
        pot = 0
        for d in derivs: # add up all independent integral contribs...
            # iterate over upper triangle coordinates (we'll add bottom contrib by symmetry)
            inds = itertools.combinations_with_replacement(range(ndim), r=d.ndim-2) if d.ndim > 2 else [()]
            for idx in inds:
                # print("--->", idx)
                count_map = {k: v for k, v in zip(*np.unique(idx, return_counts=True))}
                contrib = prefac
                for k in range(ndim): # do each dimension of integral independently
                    n = count_map.get(k, 0)
                    if n not in caches[k]:
                        caches[k][n] = (
                           self.simple_poly_int(alphas, n)
                            if expansion_type != 'taylor' else
                           self.polyint_1D(centers[..., k], alphas, n)
                        )

                    dcont = d[(slice(None, None, None), slice(None, None, None)) + idx] if len(idx) > 0 else d
                    # print(n, alphas[0, 0], centers[0, 0])
                    # with np.printoptions(linewidth=1e8):
                    #     print(caches[k][n])
                        # print(dcont/np.math.factorial(n))

                    contrib = contrib * dcont/np.math.factorial(n) * caches[k][n]
                # handle symmetry?
                pot += contrib
        return pot


    def analytic_integrate(self):
        raise NotImplementedError("flooped up")
        centers = [np.array(np.meshgrid(x, x)).T for x in self.centers.T]
        alphas = np.array(np.meshgrid(self.alphas, self.alphas)).T
        # raise Exception(alphas.shape)
        return self.potential_function['analytic_integrals'](
            centers, # ...no
            alphas
        )

    def get_V(self, potential_handler=None, expansion_degree=None, degree=None):
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
            pot_mat = self.quad_integrate(degree=self.quadrature_degree if degree is None else degree)
        elif potential_handler == 'expansion':
            pot_mat = self.expansion_integrate(deriv_order=expansion_degree)
        elif potential_handler == 'analytic':
            pot_mat = self.analytic_integrate()
        else:
            raise ValueError("woof")

        return self.S * pot_mat

    def get_wavefunctions(self):
        H = self.T + self.V
        print(np.linalg.cond(self.S))
        with np.printoptions(linewidth=1e8):
            # print(self.S)
            # print(self.T)
            print(self.V)

        try:
            return sp.linalg.eigh(H, self.S)
        except np.linalg.LinAlgError:
            raise ValueError(
                "Overlap matrix poorly conditioned ({}) usually means data is too clustered".format(np.linalg.cond(self.S))
            ) from None