import numpy as np, scipy as sp, itertools, functools
from McUtils.Zachary import RBFDInterpolator

__all__ =  [
    "DGB"
]


class DGB:
    """

    """
    def __init__(self, centers, potential_function,
                 alphas=1,
                 clustering_radius=None,
                 quadrature_degree=4
    ):
        centers = np.asanyarray(centers)
        if clustering_radius is not None and clustering_radius >= 0:
            centers, _, _ = RBFDInterpolator.decluster_data(centers, np.empty(len(centers)), [], clustering_radius)
        self.centers = centers

        if isinstance(alphas, (int, float, np.integer, np.floating)):
            alphas = np.full(len(centers), alphas)
        self.alphas = np.asanyarray(alphas)
        self.potential_function = potential_function
        self.quadrature_degree = quadrature_degree
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

    def analytic_integrate(self):
        centers = [np.array(np.meshgrid(x, x)).T for x in self.centers.T]
        alphas = np.array(np.meshgrid(self.alphas, self.alphas)).T
        # raise Exception(alphas.shape)
        return self.potential_function['analytic_integrals'](
            centers,
            alphas
        )

    def get_V(self, potential_handler=None, degree=None):

        if potential_handler is None:
            if isinstance(self.potential_function, dict):
                if 'analytic_integrals' in self.potential_function:
                    potential_handler = 'analytic'
            else:
                potential_handler = 'quad'

        if potential_handler == 'quad':
            pot_mat = self.quad_integrate(degree=self.quadrature_degree if degree is None else degree)
        elif potential_handler == 'analytic':
            pot_mat = self.analytic_integrate()
        else:
            raise ValueError("woof")

        return self.S * pot_mat

    def get_wavefunctions(self):
        H = self.T + self.V
        # print(np.linalg.cond(self.S))
        # with np.printoptions(linewidth=1e8):
        #     print(self.S)
        #     print(self.T)
        #     print(self.V)

        try:
            return sp.linalg.eigh(H, self.S)
        except np.linalg.LinAlgError:
            raise ValueError(
                "Overlap matrix poorly conditioned ({}) usually means data is too clustered".format(np.linalg.cond(self.S))
            ) from None