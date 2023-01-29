import numpy as np, scipy as sp, itertools, functools
from McUtils.Zachary import RBFDInterpolator
from McUtils.Combinatorics import StirlingS1
from McUtils.Scaffolding import Logger
from McUtils.Data import AtomData, UnitsData

__all__ =  [
    "DGB"
]

from .Wavefunctions import DGBWavefunctions

class DGB:
    """

    """

    @classmethod
    def run(cls,
            centers,
            potential_function,
            masses=None,
            atoms=None,
            projection_indices=None,
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
            reference_structure=None
            ):
        opts = dict(
            masses=masses,
            atoms=atoms,
            alphas=alphas,
            projection_indices=projection_indices,
            logger=logger,
            clustering_radius=clustering_radius,
            min_singular_value=min_singular_value,
            num_svd_vectors=num_svd_vectors,
            svd_contrib_cutoff=svd_contrib_cutoff,
            optimize_centers=optimize_centers,
            quadrature_degree=quadrature_degree,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            reference_structure=reference_structure
        )
        opts = {k:v for k,v in opts.items() if v is not None}
        logger = Logger.lookup(logger)
        with logger.block(tag="Running distributed Gaussian basis calculation"):
            opts['logger'] = logger
            ham = cls(centers, potential_function, **opts)
            return ham.get_wavefunctions()

    def __init__(self,
                 centers,
                 potential_function,
                 masses=None,
                 atoms=None,
                 alphas=None,
                 projection_indices=None,
                 transformations=None,
                 logger=False,
                 clustering_radius=.005,
                 min_singular_value=1e-4,
                 num_svd_vectors=None,
                 svd_contrib_cutoff=1e-3,
                 optimize_centers=False,
                 quadrature_degree=4,
                 expansion_degree=None,
                 expansion_type='multicenter',
                 reference_structure=None
    ):
        self._S, self._T, self._V = None, None, None
        self.logger = Logger.lookup(logger)

        centers = np.asanyarray(centers)
        if masses is None:
            if atoms is not None:
                atoms = [
                    [AtomData[a, "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")] * 3
                        if isinstance(a, str) else
                    a * 3
                    for a in atoms
                ]
                masses = np.array(atoms).flatten()
            else:
                masses = [1] * centers.shape[-1]
        if isinstance(masses, (int, np.integer, float, np.floating)):
            masses = [masses]
        self.masses = np.asanyarray(masses)
        self.inds = projection_indices
        # print("MAsSES:", masses)

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

        if self.alphas.ndim == 1:
            self.alphas = np.broadcast_to(self.alphas[:, np.newaxis], self.centers.shape)

        self.min_singular_value = min_singular_value
        if min_singular_value is not None or num_svd_vectors is not None:
            # Use SVD to prune out matrix rows that will be super ill conditioned
            U, sig, VT = np.linalg.svd(self.S)
            if num_svd_vectors:
                good_loc = np.arange(np.min([num_svd_vectors, len(sig)]))
            else:
                self.logger.log_print("most important center threshold: {t}", t=min_singular_value)
                good_loc = np.where(sig > min_singular_value)[0]
            U = U[:, good_loc]
            VT = VT[good_loc, :]
            # raise Exception(np.min(np.abs(U)))
            full_good_pos = np.unique(np.concatenate([
                np.where(np.abs(U) > svd_contrib_cutoff)[0],
                np.where(np.abs(VT) > svd_contrib_cutoff)[1]
            ]))
            self.centers = self.centers[full_good_pos]
            self.alphas = self.alphas[full_good_pos]

            self._S = None#self._S[full_good_pos, :][:, full_good_pos]
            self._T = None#self._T[full_good_pos, :][:, full_good_pos]

        self.logger.log_print("Number of centers: {N}", N=len(self.centers))

        self.potential_function = potential_function
        self.quadrature_degree = quadrature_degree
        self.expansion_degree = expansion_degree
        self.expansion_type = expansion_type
        self.ref = reference_structure

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

        mask = None
        if clustering_radius is not None and clustering_radius >= 0:
            centers, _, _, mask = RBFDInterpolator.decluster_data(centers, np.empty(len(centers)), [], clustering_radius, return_mask=True)

        if alphas is None:
            alphas = self.get_alphas(centers, clustering_radius)
        if isinstance(alphas, (int, float, np.integer, np.floating)):
            alphas = np.full(len(centers), alphas)
        else:
            alphas = np.asanyarray(alphas)
            if mask is not None:
                alphas = alphas[mask]

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
        # print("???", a[:5])
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
        C = arat * np.power(disps, 2)

        # Base components

        S_dim = (np.sqrt(2) * np.power(aouter, 1/4) / B) * np.exp(-C)
        T_dim = arat * (1 - 2*C) / self.masses[np.newaxis, np.newaxis, :]
        if self.inds is not None:
            S_dim = S_dim[:, :, self.inds]
            T_dim = T_dim[:, :, self.inds]

        # Combine appropriately
        S = np.prod(S_dim, axis=-1)
        T = S * np.sum(T_dim, axis=-1)

        return S, T

    def get_overlap_gaussians(self):
        # find overlap gaussians
        new_alphas = self.alphas[:, np.newaxis, :] + self.alphas[np.newaxis, :, :]
        w_centers = self.alphas*self.centers
        # moving weighted average by alpha value
        return (w_centers[:, np.newaxis, :] + w_centers[np.newaxis, :, :])/new_alphas, new_alphas

    def quad_integrate(self, function, degree=2):
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
        pots = np.zeros(alphas.shape[:2])
        ndim = centers.shape[-1]
        for disp_inds in itertools.product(*([range(degree)]*ndim)):
            disp_inds = np.array(disp_inds)
            w = np.prod(weights[disp_inds])
            c = centers + disps[disp_inds][np.newaxis, np.newaxis, :] / np.sqrt(alphas[:, :])
            pots = pots + w * function(c)

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
    def poch(cls, n, m):
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
    def polyint_1D(cls, centers, n):
        if n == 0:
            return np.ones(centers.shape[:2])
        term = sum(
            (cls.poch(n, l) * (1/2**(2*l-n) if 2*l > n else 2**(n-2*l)))*(centers)**(n-2*l)
            for l in range(0, int(np.floor(n/2)) + 1)
        )
        return term
    @classmethod
    def simple_poly_int(cls, n):
        return np.prod(np.arange(1, n, 2)) # double factorial/gamma/whatever

    def expansion_integrate(self, function, deriv_order=2, expansion_type=None):
        if expansion_type is None:
            expansion_type = self.expansion_type

        centers, alphas = self.get_overlap_gaussians()
        # this prefac allows us to reexpress things in terms of the new Gaussians
        # with appropriate weighting
        aprod = self.alphas[:, np.newaxis, :] * self.alphas[np.newaxis, :, :]
        if self.inds is not None:
            aprod = aprod[..., self.inds]
        # asum = self.alphas[:, np.newaxis] + self.alphas[np.newaxis, :]
        # cdiff = self.centers[:, np.newaxis, :] - self.centers[np.newaxis, :, :]
        ndim = centers.shape[-1]
        scaling_dim = ndim if self.inds is None else len(self.inds)
        # prefac = (
        #         ( np.sqrt(2 / np.pi) ) ** (scaling_dim if expansion_type != 'taylor' else 1)
        #         * np.product(aprod, axis=-1)**(1/4)
        # )

        if expansion_type == 'taylor':
            self.logger.log_print("expanding as a Taylor series about the minumum energy geometry...")
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
            row_inds, col_inds = np.triu_indices(alphas.shape[0])
            self.logger.log_print("expanding about {N} points...", N=len(row_inds))
            # derivs = self.potential_function(centers, deriv_order=deriv_order)
            # if isinstance(derivs, np.ndarray):  # didn't get the full list so we do the less efficient route'
            #     derivs = [self.potential_function(centers)] + [self.potential_function(centers, deriv_order=d) for d in
            #                                                range(1, deriv_order + 1)]

            # use symmetry to cut down number of indices by 2
            ics = centers[row_inds, col_inds]
            derivs = function(ics, deriv_order=deriv_order)
            if isinstance(derivs, np.ndarray):  # didn't get the full list so we do the less efficient route'
                derivs = [function(ics)] + [
                    function(ics, deriv_order=d) for d in range(1, deriv_order+1)
                ]
            new_derivs = [] # reverse symmetrization
            for d in derivs:
                full_mat = np.empty(
                    alphas.shape[:2] + d.shape[1:]
                )
                full_mat[row_inds, col_inds] = d
                full_mat[col_inds, row_inds] = d
                new_derivs.append(full_mat)
            derivs = new_derivs

        if self.inds is not None:
            ndim = len(self.inds)
            centers = centers[..., self.inds]
            alphas = alphas[..., self.inds]
            new_derivs = []
            for n,d in enumerate(derivs):
                if n > 0:
                    d_sel = (slice(None, None, None), slice(None, None, None),) + np.ix_(*[self.inds]*n)
                    new_derivs.append(d[d_sel])
                else:
                    new_derivs.append(d)
            derivs = new_derivs

        self.logger.log_print("adding up all derivative contributions...")
        caches = [{} for _ in range(ndim)]
        pot = 0
        for d in derivs: # add up all independent integral contribs...
            # iterate over upper triangle coordinates (we'll add bottom contrib by symmetry)
            inds = itertools.combinations_with_replacement(range(ndim), r=d.ndim-2) if d.ndim > 2 else [()]
            for idx in inds:
                count_map = {k: v for k, v in zip(*np.unique(idx, return_counts=True))}
                if expansion_type != 'taylor' and any(n%2 !=0 for n in count_map.values()):
                    continue # odd contribs vanish

                contrib = 1
                for k in range(ndim): # do each dimension of integral independently
                    n = count_map.get(k, 0)
                    if n not in caches[k]:
                        caches[k][n] = (
                           self.simple_poly_int(n)
                             if expansion_type != 'taylor' else
                           self.polyint_1D(centers[..., k], n)
                        )
                    contrib *= caches[k][n] / alphas[..., k]**(n/2)

                dcont = d[(slice(None, None, None), slice(None, None, None)) + idx] if len(idx) > 0 else d
                facterms = np.unique([x for x in itertools.permutations(idx)], axis=0)
                nfac = len(facterms) # this is like a binomial coeff or something but my sick brain won't work right now...
                scaling = 2**(len(idx)/2) * np.prod([np.math.factorial(count_map.get(k, 0)) for k in range(ndim)])

                contrib *= nfac * dcont / scaling

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
            pot_mat = self.quad_integrate(self.potential_function, degree=self.quadrature_degree if degree is None else degree)
        elif potential_handler == 'expansion':
            self.logger.log_print("evauating integrals with {n}-degree expansions", n=self.expansion_degree)
            pot_mat = self.expansion_integrate(self.potential_function, deriv_order=expansion_degree)
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

    def diagonalize(self, print_debug_info=False, min_singular_value=None):

        if min_singular_value is None:
            min_singular_value = self.min_singular_value
        self.logger.log_print("solving with min_singular_value={ms}", ms=min_singular_value)

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

    def get_wavefunctions(self, print_debug_info=False, min_singular_value=None):
        eigs, evals = self.diagonalize(print_debug_info=print_debug_info, min_singular_value=min_singular_value)
        return DGBWavefunctions(eigs, evals, hamiltonian=self)