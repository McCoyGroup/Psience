import numpy as np, scipy as sp, itertools, functools
from McUtils.Zachary import RBFDInterpolator
from McUtils.Scaffolding import Logger
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput

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
            reference_structure=None
            ):
        opts = dict(
            masses=masses,
            atoms=atoms,
            alphas=alphas,
            projection_indices=projection_indices,
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
                 optimize_centers=False,
                 clustering_radius=-1,
                 min_singular_value=1e-4,
                 num_svd_vectors=None,
                 svd_contrib_cutoff=1e-3,
                 quadrature_degree=4,
                 expansion_degree=None,
                 expansion_type='multicenter',
                 reference_structure=None
    ):
        self._S, self._T, self._V = None, None, None
        self.logger = Logger.lookup(logger)

        self.potential_function = potential_function
        self.quadrature_degree = quadrature_degree
        self.expansion_degree = expansion_degree
        self.expansion_type = expansion_type
        self.ref = reference_structure

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
            self.centers, self.alphas, tfs = self.initialize_gaussians(
                centers, alphas,
                clustering_radius
            )
            if transformations is None and tfs is not None:
                transformations = tfs

        if self.alphas.ndim == 1:
            self.alphas = np.broadcast_to(self.alphas[:, np.newaxis], self.centers.shape)
        self.transformations = transformations

        self.min_singular_value = min_singular_value
        if optimize_centers:
            if min_singular_value is not None or num_svd_vectors is not None:
                # Use SVD to prune out matrix rows that will be super ill conditioned
                sig, evecs = np.linalg.eigh(self.S)
                if num_svd_vectors:
                    good_loc = slice(max(len(sig)-num_svd_vectors, 0), len(sig))
                else:
                    self.logger.log_print("most important center threshold: {t}", t=min_singular_value)
                    good_loc = np.where(sig > min_singular_value)[0]
                # raise Exception(np.min(np.abs(U)))
                full_good_pos = np.unique(np.where(np.abs(evecs[:, good_loc]) > svd_contrib_cutoff)[0])
                self.centers = self.centers[full_good_pos]
                self.alphas = self.alphas[full_good_pos]
                if self.transformations is not None:
                    self.transformations = self.transformations[full_good_pos]

                self._S = None#self._S[full_good_pos, :][:, full_good_pos]
                self._T = None#self._T[full_good_pos, :][:, full_good_pos]

        self.logger.log_print("Number of centers: {N}", N=len(self.centers))

    def optimize_centers(self,
                         centers, alphas,
                         max_condition_number=1e16,
                         initial_custering=.005,
                         cluster_step_size=.005,
                         max_steps = 50
                         ):
        raise NotImplementedError("ugh...")
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

        rots = None
        if alphas is None:
            alphas = self.get_alphas(self.masses, centers)#, clustering_radius)
        elif isinstance(alphas, (str, dict)):
            alphas, rots = self.dispatch_get_alphas(alphas, centers)

        if isinstance(alphas, (int, float, np.integer, np.floating)):
            alphas = np.full(len(centers), alphas)
        else:
            alphas = np.asanyarray(alphas)
            if mask is not None:
                alphas = alphas[mask]

        return centers, alphas, rots

    def dispatch_get_alphas(self, alphas, centers):
        if isinstance(alphas, str):
            alphas = {'method':alphas}
        opts = alphas.copy()
        del opts['method']
        if alphas['method'] == 'virial':
            return self.get_virial_alphas(
                centers,
                self.masses,
                self.potential_function,
                **opts
            )
        elif alphas['method'] == 'min_dist':
            return self.get_alphas(
                self.masses,
                centers,
                **opts
            ), None
        else:
            raise ValueError("unknown method for getting alphas {}".format(alphas['method']))

    @classmethod
    def get_alphas(cls, masses, centers, scaling=1/4, use_mean=False):
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

    @classmethod
    def get_virial_alphas(cls, pts, masses, potential_function, allow_rotations=False, min_frequency=2e-3, scaling=2):
        """
        Provides a way to get alphas that satisfy the virial theorem locally
        for a quadratic expansion about each center

        :param pts:
        :type pts:
        :param potential_function:
        :type potential_function:
        :param masses:
        :type masses:
        :param allow_rotations:
        :type allow_rotations:
        :return:
        :rtype:
        """

        derivs = potential_function(pts, deriv_order=2)
        if len(derivs) != 3 or not (
                isinstance(derivs[0], np.ndarray) and derivs[0].shape == (len(pts),)
        ):
            hess = derivs
            grads = potential_function(pts, deriv_order=1)
        else:
            _, grads, hess = derivs

        if allow_rotations:
            npts = len(pts)
            ndim = pts.shape[-1]

            alphas = np.zeros((npts, ndim))
            rots = np.zeros((npts, ndim, ndim))

            rm = masses
            grads = grads / np.sqrt(rm[np.newaxis])  # mass-weight
            hess = hess / np.sqrt(rm[np.newaxis, :, np.newaxis] * rm[np.newaxis, np.newaxis, :])

            grad_norms = np.linalg.norm(grads, axis=1)
            non_stationary = grad_norms > 1e-6
            stationary = np.where(grad_norms <= 1e-6)
            if len(stationary) == 0:
                stationary = np.array([], dtype=int)
            else:
                stationary = stationary[0]

            # obviously kinda an adaptation...but the reduced masses are the smae
            # for both coords so it's just a scaling factor
            rp_mode = grads[non_stationary] / grad_norms[non_stationary][:, np.newaxis]
            num_rp = np.sum(non_stationary.astype(int))
            proj = np.broadcast_to(np.eye(ndim)[np.newaxis], (num_rp, ndim, ndim)) - nput.vec_outer(rp_mode, rp_mode)
            h2 = proj @ hess[non_stationary] @ proj
            freqs, modes = np.linalg.eigh(h2)
            modes[:, :, 1] = modes[:, :, 1] * np.linalg.det(modes)[:, np.newaxis]

            freq_cuts = np.abs(freqs) < 1e-8
            kill_pos = np.where(np.all(freq_cuts, axis=1))
            if len(kill_pos) > 0 and len(kill_pos[0]) > 0:
                sel = np.where(non_stationary)[0][kill_pos]
                stationary = np.unique(np.concatenate([stationary, sel]))

            zi = np.where(freq_cuts)
            for i, j in zip(*zi):
                m = modes[i, :, j][:, np.newaxis]
                f = m.T @ hess[i] @ m
                freqs[i, j] = f

            freqs = np.sqrt(np.abs(freqs))
            freqs[freqs < min_frequency] = min_frequency

            g = np.diag(1 / masses)
            g = modes.transpose(0, 2, 1) @ g[np.newaxis] @ modes
            rpms = 1 / np.diagonal(g, axis1=1, axis2=2)
            alphas[non_stationary] = scaling * rpms * freqs
            rots[non_stationary] = modes

            # rp_coords = rp_mode[:, np.newaxis, np.newaxis, :] @ pts[np.newaxis, :, :, np.newaxis]
            # rp_coords = rp_coords.reshape((len(non_stationary), len(pts)))
            # raise Exception(rp_coords)

            if len(stationary) > 0:
                freqs2, modes = np.linalg.eigh(hess[stationary])
                freqs = np.sqrt(np.abs(freqs2))
                g = np.diag(1 / masses)
                g = modes.transpose(0, 2, 1) @ g[np.newaxis] @ modes
                rpms = 1 / np.diagonal(g, axis1=1, axis2=2)  # note this is the same as masses...
                alphas[stationary] = scaling * rpms * freqs
                rots[stationary] = modes


        else:
            rots = None
            alphas = scaling * np.sqrt(np.abs(masses * np.diagonal(hess, axis1=1, axis2=2)))

        return alphas, rots

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

    def get_inverse_covariances(self, alphas=None, transformations=None, inds=None):
        """
        Transforms the alphas into proper inverse covariance matrices

        :return:
        :rtype:
        """

        if transformations is None:
            transformations = self.transformations
        if transformations is None:
            return None

        if alphas is None:
            alphas = self.alphas

        if inds is not None:

            alphas = alphas
            transformations = transformations[:, self.inds]
            n = alphas.shape[-1]
            npts = len(alphas)
            diag_covs = np.zeros((npts, n, n))
            diag_inds = (slice(None, None, None),) + np.diag_indices(n)
            diag_covs[diag_inds] = 1 / (2 * alphas)
            covs = transformations @ diag_covs @ transformations.transpose(0, 2, 1)
            covs[np.abs(covs) < 1e-12] = 0  # numerical garbage can be an issue...

            two_a_inv, transformations = np.linalg.eigh(covs)
            alphas = 1/(2*two_a_inv)

        n = alphas.shape[-1]
        npts = len(alphas)
        diag_covs = np.zeros((npts, n, n))
        diag_inds = (slice(None, None, None),) + np.diag_indices(n)
        diag_covs[diag_inds] = 2*alphas

        covs = transformations @ diag_covs @ transformations.transpose(0, 2, 1)
        covs[np.abs(covs) < 1e-12] = 0 # numerical garbage can be an issue...

        return covs

    def get_overlap_gaussians(self):
        rows, cols = np.triu_indices(len(self.alphas))
        if self.transformations is None:
            # find overlap gaussians

            centers = self.centers
            alphas = self.alphas
            if self.inds is not None:
                centers = centers[:, self.inds]
                alphas = alphas[:, self.inds]
            new_alphas = alphas[rows] + alphas[cols]
            w_centers = alphas*centers
            # moving weighted average by alpha value
            overlap_data = (w_centers[rows] + w_centers[cols])/new_alphas, new_alphas
        else:
            centers = self.centers
            if self.inds is not None:
                centers = centers[:, self.inds]

            sigs = self.get_inverse_covariances()
            new_sigs = sigs[rows] + sigs[cols]
            new_inv = np.linalg.inv(new_sigs)
            new_centers = new_inv@(
                sigs[rows] @ centers[rows][:, :, np.newaxis]
                + sigs[cols] @ centers[cols][:, :, np.newaxis]
            )
            new_centers = new_centers.reshape(self.centers[cols].shape)
            new_alphas, new_rots = np.linalg.eigh(new_sigs) # eigenvalues of inverse tensor...
            new_alphas = new_alphas/2
            sum_sigs = sigs[rows]@new_inv@sigs[cols]

            overlap_data = {
                'centers': new_centers,
                'alphas': new_alphas,
                'sigmas':new_sigs,
                'rotations': new_rots,
                'row_sigs': sigs[rows],
                'col_sigs': sigs[cols],
                'sum_inverse':sum_sigs
            }

        return overlap_data

    def get_ST(self, centers=None, alphas=None, transformations=None):
        if centers is None:
            centers = self.centers
        if alphas is None:
            alphas = self.alphas
        if transformations is None:
            transformations = self.transformations

        if transformations is None:

            if self.inds is not None:
                centers = centers[:, self.inds]
                alphas = alphas[:, self.inds]

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

            # if self.inds is not None:
            #     S_dim = S_dim[:, :, self.inds]
            #     T_dim = T_dim[:, :, self.inds]

            # Combine appropriately
            S = np.prod(S_dim, axis=-1)
            T = S * np.sum(T_dim, axis=-1)
        else:
            T = np.zeros((len(self.centers), len(self.centers)))
            S = np.zeros((len(self.centers), len(self.centers)))
            row_inds, col_inds = np.triu_indices(len(self.alphas))
            rot_data = self.get_overlap_gaussians()

            alphas = self.alphas
            centers = self.centers
            masses = self.masses
            if self.inds is not None:
                alphas = alphas[:, self.inds]
                centers = centers[:, self.inds]
                masses = masses[self.inds]


            ndim = centers.shape[-1]
            dets = np.prod(rot_data["alphas"], axis=-1)
            rows = rot_data['row_sigs']
            cols = rot_data['col_sigs']
            # we prefactor out the 2**ndim
            rdets = np.prod(alphas[row_inds], axis=-1)
            cdets = np.prod(alphas[col_inds], axis=-1) # literally a product of passed in alphas
            disps = centers[row_inds] - centers[col_inds]
            C = disps[:, np.newaxis, :]@rot_data['sum_inverse']@disps[:, :, np.newaxis]
            C = C.reshape(disps.shape[0])

            # raise Exception(
            #     (2**ndim),
            #     ((rdets*cdets)/(dets**2))
            # )

            S[row_inds, col_inds] = 2**(ndim/2) * ((rdets*cdets)/(dets**2))**(1/4) * np.exp(-C/2)

            L = np.transpose(rot_data['rotations'], (0, 2, 1))
            Lt = rot_data['rotations']
            zetas =L@(rot_data['centers'] - centers[col_inds])[:, :, np.newaxis] # seems weird but we get the symmetry back in the end
            zetas = zetas.reshape(rot_data['centers'].shape)

            minv = np.diag(1/masses)
            Sj = cols
            Amat = L@(Sj@minv@Sj)@Lt
            # msj = rot_data['rotations']@(Sj@minv)@np.transpose(rot_data['rotations'], (0, 1, 2))
            # raise Exception(
            #     np.diagonal(msj, axis1=1, axis2=2),# / self.masses[np.newaxis],
            #     - 1 / 2 * np.diagonal(Amat, axis1=1, axis2=2) / rot_data['alphas']
            # )

            T[row_inds, col_inds] = 1 / 2 * (
                    np.sum(
                        2 * alphas[col_inds] / masses[np.newaxis] -
                         1 / 2 * np.diagonal(Amat, axis1=1, axis2=2) / rot_data['alphas'],
                        axis=1
                    )
                    - 1/2 * sum(
                        # easier this way than the proper series summation...
                        Amat[..., k, kp] * (
                            4 * zetas[..., k] * zetas[..., kp]  # I discovered I was off by a factor of 4 in the 2D case...
                                if k != kp else
                            2 * zetas[..., k] ** 2
                        )
                        for k, kp in zip(*np.triu_indices(Amat.shape[-1]))
                    )  # could be a dot but this is fine
            )
            T[row_inds, col_inds] *= S[row_inds, col_inds]

            S[col_inds, row_inds] = S[row_inds, col_inds]
            T[col_inds, row_inds] = T[row_inds, col_inds]

        return S, T

    def quad_integrate(self, function, degree=2):
        """
        Integrate potential over all pairs of Gaussians at once

        :param degree:
        :type degree:
        :return:
        :rtype:
        """

        if self.transformations is not None:
            raise NotImplementedError("quadrature in rotated basis not implemented yet")

        # TODO: add in ability to do quadrature along ellipsoid axes

        # Quadrature point displacements and weights (thanks NumPy!)
        disps, weights = np.polynomial.hermite.hermgauss(degree)

        centers, alphas = self.get_overlap_gaussians() # Only upper triangle here
        npts = len(self.alphas)
        rows, cols = np.triu_indices(npts)
        pots = np.zeros((npts, npts))
        ndim = centers.shape[-1]
        for disp_inds in itertools.product(*([range(degree)]*ndim)):
            disp_inds = np.array(disp_inds)
            w = np.prod(weights[disp_inds])
            c = centers + disps[disp_inds][np.newaxis, :] / np.sqrt(alphas)
            pots[rows, cols] += w * function(c)
        pots[cols, rows] = pots[rows, cols]

        normalization = 1 / (np.sqrt(np.pi)) ** self.centers.shape[-1]
        return pots * normalization

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
        return np.prod(np.arange(1, n, 2)) # double factorial/gamma/whatever

    def expansion_integrate(self, function, deriv_order=2, expansion_type=None):
        if expansion_type is None:
            expansion_type = self.expansion_type

        if self.transformations is None:
            centers, alphas = self.get_overlap_gaussians()
            rotations = None
        else:
            rot_data = self.get_overlap_gaussians()
            centers = rot_data['centers']
            alphas = rot_data['alphas']
            rotations = rot_data['rotations'].transpose(0, 2, 1)

        ndim = centers.shape[-1]
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
            deriv_order = deriv_order - (deriv_order % 2) # odd orders don't contribute so why evaluate the derivatives...
            self.logger.log_print("expanding about {N} points...", N=len(alphas))
            derivs = function(centers, deriv_order=deriv_order)
            if isinstance(derivs, np.ndarray):  # didn't get the full list so we do the less efficient route'
                derivs = [function(centers)] + [
                    function(centers, deriv_order=d) for d in range(1, deriv_order+1)
                ]

        if rotations is not None:
            new_derivs = []
            # rotations = rotations[:, :, :, np.newaxis] # to test shapes
            for n,d in enumerate(derivs):
                for _ in range(n):
                    d = nput.vec_tensordot(
                        d, rotations,
                        axes=[1, 2],
                        shared=1
                    )
                new_derivs.append(d)
            derivs = new_derivs

        self.logger.log_print("adding up all derivative contributions...")
        row_inds, col_inds = np.triu_indices(len(self.centers))
        fdim = derivs[0].ndim - 1
        fshape = derivs[0].shape[1:]
        pot = np.zeros((len(self.centers), len(self.centers)) + fshape)
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
                           self.simple_poly_int(n)
                             if expansion_type != 'taylor' else
                           self.polyint_1D(centers[..., k], alphas[..., k], n)
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

                contrib *= nfac * dcont / scaling

                pot[row_inds, col_inds] += contrib

        pot[col_inds, row_inds] = pot[row_inds, col_inds]

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

    def evaluate_multiplicative_operator_base(self, function,
                                              handler=None,
                                              expansion_degree=None,
                                              expansion_type=None,
                                              quadrature_degree=None
                                              ):
        if expansion_degree is None:
            expansion_degree = self.expansion_degree

        if handler is None:
            if isinstance(function, dict):
                if 'analytic_integrals' in function:
                    handler = 'analytic'
            elif expansion_degree is not None:
                handler = 'expansion'
            else:
                handler = 'quad'

        if handler == 'quad':
            self.logger.log_print("evauating integrals with {n}-order quadrature", n=self.quadrature_degree)
            pot_mat = self.quad_integrate(function, degree=self.quadrature_degree if quadrature_degree is None else quadrature_degree)
        elif handler == 'expansion':
            self.logger.log_print("evauating integrals with {n}-degree expansions", n=self.expansion_degree)
            pot_mat = self.expansion_integrate(function, deriv_order=expansion_degree, expansion_type=expansion_type)
        elif handler == 'analytic':
            self.logger.log_print("evauating integrals analytically", n=self.expansion_degree)
            pot_mat = self.analytic_integrate()
        else:
            raise ValueError("unknown operator evaluation scheme {}".format(handler))

        return pot_mat

    def evaluate_multiplicative_operator(self, function,
                                         handler=None,
                                         expansion_degree=None,
                                         expansion_type=None,
                                         quadrature_degree=None
                                         ):

        pot_mat = self.evaluate_multiplicative_operator_base(
            function,
            handler=handler,
            expansion_degree=expansion_degree, expansion_type=expansion_type,
            quadrature_degree=quadrature_degree
        )
        S = self.S
        for _ in range(pot_mat.ndim - 2):
            S = np.expand_dims(S, -1)
        return S * pot_mat

    def get_V(self, potential_handler=None, expansion_degree=None, expansion_type=None, quadrature_degree=None):
        self.logger.log_print("calculating potential matrix")
        return self.evaluate_multiplicative_operator(
            self.potential_function,
            handler=potential_handler,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            quadrature_degree=quadrature_degree
        )

    def get_orthogonal_transform(self, min_singular_value=None, subspace_size=None):
        if min_singular_value is None:
            min_singular_value = self.min_singular_value

        S = self.S
        sig, evecs = np.linalg.eigh(S)

        if subspace_size is not None:
            good_loc = slice(max(0, len(S) - subspace_size), len(S))
        elif min_singular_value is not None:
            good_loc = np.where(sig > min_singular_value * np.max(sig))[0]
        else:
            good_loc = np.arange(len(S))
        d = np.diag(1 / np.sqrt(sig[good_loc]))
        L = evecs[:, good_loc]
        Q = L @ d @ L.T

        # sorting = np.concatenate([good_loc, bad_loc[0]])
        # Lsort = evecs[:, sorting]
        Qinv = L @ np.diag(np.sqrt(sig[good_loc])) @ L.T

        if L.shape[0] == L.shape[1]: # no contraction
            return Q, Qinv, None

        Qe, QL = np.linalg.eigh(Q)
        qsub = np.where(Qe > 1e-8)[0]
        Qq = QL[:, qsub]
        qrest = np.where(Qe <= 1e-8)[0]
        qsort = np.concatenate([qsub, qrest])
        Qqinv = QL[:, qsort].T  # transforms back to the Q basis

        return Q, Qinv, (Qq, Qqinv)

    def diagonalize(self, print_debug_info=False, subspace_size=None, min_singular_value=None,
                    eps=1e-12,
                    mode='classic',
                    nodeless_ground_state=True
                    ):

        if min_singular_value is None:
            min_singular_value = self.min_singular_value
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
            Q, Qinv, proj = self.get_orthogonal_transform(min_singular_value=min_singular_value, subspace_size=subspace_size)
            if proj is None:
                # print(Q.shape, H.shape, self.S.shape, self.centers.shape)
                eigs, evecs = sp.linalg.eigh(H, self.S)
                Qq = np.eye(len(Q))
            else:
                Qq, Qqinv = proj
                H = Qq.T @ Q @ H @ Q.T @ Qq # in our projected orthonormal basis
                eigs, evecs = np.linalg.eigh(H)
                evecs = np.concatenate(
                    [
                        evecs,
                        np.zeros((Qqinv.shape[1] - Qq.shape[1], len(evecs)))
                    ],
                    axis=0
                )
                evecs = Q.T @ Qqinv.T @ evecs
            if nodeless_ground_state:
                # fast enough if we have not that many points...
                gswfn = DGBWavefunctions(eigs, evecs, hamiltonian=self)[0]
                gs = gswfn.evaluate(self.centers)
                abs_gs = np.abs(gs)
                signs = np.sign(gs[abs_gs > np.max(abs_gs) * 1e-1])
                diffs = np.abs(np.diff(signs))
                # print(gs[np.argsort(np.abs(gs))][-5:], Qq.shape[1] - 1)
                if np.sum(diffs) > 0: # had a sign flip
                    eigs, evecs = self.diagonalize(
                        mode='classic',
                        subspace_size=Qq.shape[1] - 1,
                        nodeless_ground_state=True
                    )

        elif mode == 'stable': # Implementation of the Fix-Heiberger algorithm
            raise NotImplementedError("I think my Fix-Heiberger is broken?"
                                      "And also potentially not even the right algorithm for this problem"
                                      "Like I potentially _want_ a large epsilon since I want to get the best"
                                      "possible energies, but Fix-Heiberger will call that unstable I think")

            S = self.S
            d, Q = np.linalg.eigh(S)
            d = np.flip(d)
            Q = np.flip(Q, axis=1)

            N = len(d)
            cutoff = np.max(d) * eps
            g1 = np.where(d >= cutoff)
            if len(g1) == 0:
                raise ValueError("totally zero S matrix")
            n1 = len(g1[0])
            n2 = N - n1

            if n2 == 0:
                print("Falling back on classic method...")
                return self.diagonalize(mode='classic')

            # D = np.diag(d[g1])
            # F = np.diag(d[g2])
            A = Q.T@H@Q

            R = np.diag(np.concatenate([
                1 / np.sqrt(d[g1]),
                np.ones(n2)
            ]))
            A1 = R @ A @ R

            # B = Q.T@S@Q
            # B11 = R @ B @ R
            # B1 = np.diag(np.concatenate([
            #     np.ones(n1),
            #     np.zeros(n2) #d[g2]
            # ]))

            A22 = A1[n1:, n1:]
            d2, Q22 = np.linalg.eigh(A22)
            d2 = np.flip(d2)
            Q22 = np.flip(Q22, axis=1)

            cut2 = np.max(d2) * eps
            g3 = np.where(d2 >= cut2)[0]
            n3 = len(g3)
            if n3 == 0:
                # print("Early exiting after first iteration")
                if n1 <= n2:
                    raise ValueError("singular problem (case 2)")
                A12 = A1[:n1, n1:]
                Q12, R13 = np.linalg.qr(A12, mode='complete')
                if np.linalg.matrix_rank(A12) < n2:  # singular
                    raise ValueError("singular problem (case 3)")

                # print("N3 = 0???")

                Q2 = np.eye(N)
                Q2[:n1, :n1] = Q12
                A2 = Q2.T @ A1 @ Q2.T

                A12 = A2[:n2, n2:n1]
                A13 = A2[n1:, :n2]
                A22 = A2[n2:n1, n2:n1]

                eigs, U2 = np.linalg.eigh(A22)
                U1 = np.zeros((n2, n1 - n2))
                U3 = -np.linalg.inv(A13)@A12@U2

                U = np.concatenate(
                    [
                        U1, U2, U3
                    ],
                    axis=0
                )

                evecs = Q @ R @ Q2 @ U

            else:
                g4 = np.where(d2 < cut2)[0]
                n4 = len(g4)
                Q2 = np.eye(N)
                Q2[n1:, n1:] = Q22
                A2 = Q2.T@A1@Q2
                # B2 = Q2.T@B1@Q2

                if n4 == 0:
                    # print("Early, well conditioned after first iteration")
                    # print("N4 = 0??? {d2}".format(d2=d2))
                    A11 = A2[:n1, :n1]
                    A12 = A2[:n1, n1:]
                    Dinv = np.diag(1 / d2)

                    eigs, U1 = np.linalg.eigh(A11 - A12@Dinv@A12.T)
                    U2 = -Dinv@A12.T@U1

                    U = np.concatenate(
                        [
                            U1, U2
                        ],
                        axis=0
                    )

                    evecs = Q @ R @ Q2 @ U
                else: # second iteration of this partitioning...
                    if n1 <= n4:
                        raise ValueError("singular problem (case 5)")

                    A2[-n4:, -n4:] = np.zeros((n4, n4)) #
                    # B2 = B1

                    A13 = A2[:n1, n1+n3:]
                    Q33, R14 = np.linalg.qr(A13, mode='complete')
                    if np.linalg.matrix_rank(A13) < n4: # singular
                        raise ValueError("singular problem (case 6)")
                    # A14 = R14[:n4]

                    Q3 = np.eye(N)
                    Q3[:n1, :n1] = Q33
                    A3 = Q3.T@A2@Q3

                    # B3 = B1

                    n5 = n1 - n4
                    A11 = A3[:n4, :n4]
                    A12 = A3[:n4, n4:n1]
                    A13 = A3[:n4, n1:n1+n3]
                    A14 = A3[:n4, -n4:]

                    A22 = A3[n4:n1, n4:n1]
                    A23 = A3[n4:n1, n1:n1+n3]
                    # A24 = A3[n4:n5, -n4:]
                    #
                    # raise Exception(A24)

                    Dinv = np.diag(1 / d2[g3])

                    U1 = np.zeros((n4, n5))
                    eigs, U2 = np.linalg.eigh(
                        A22 - A23@Dinv@A23.T
                    )
                    U3 = -Dinv@A23.T@U2
                    U4 = -np.linalg.inv(A14)@(A12@U2 + A13@U3)

                    U = np.concatenate(
                        [
                            U1, U2, U3, U4
                        ],
                        axis=0
                    )

                    evecs = Q @ R @ Q2 @ Q3 @ U

        else:
            raise ValueError("unknown solver {}".format(mode))

                # raise Exception((eigs[1:] - eigs[0])*219475.6)


        return eigs, evecs

    def get_wavefunctions(self, print_debug_info=False, min_singular_value=None, subspace_size=None, nodeless_ground_state=True):
        # print("======="*25)
        eigs, evecs = self.diagonalize(print_debug_info=print_debug_info,
                                       min_singular_value=min_singular_value,
                                       nodeless_ground_state=nodeless_ground_state,
                                       subspace_size=subspace_size
                                       )
        return DGBWavefunctions(eigs, evecs, hamiltonian=self)