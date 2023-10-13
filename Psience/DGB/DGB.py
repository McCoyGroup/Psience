import numpy as np, scipy as sp, itertools, functools

from McUtils.Zachary import RBFDInterpolator, DensePolynomial, TensorDerivativeConverter
from McUtils.Scaffolding import Logger
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput

__all__ = [
    "DGB"
]

from .Wavefunctions import DGBWavefunctions
from ..Molecools import Molecule

class DGB:
    """

    """

    @classmethod
    def run(cls,
            centers,
            potential_function,
            masses=None,
            mass_weighted=False,
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
            reference_structure=None,
            pairwise_potential_functions=None # we can integrate this cleanly, not sure about anything else
            ):
        opts = dict(
            masses=masses,
            mass_weighted=mass_weighted,
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
            reference_structure=reference_structure,
            pairwise_potential_functions=pairwise_potential_functions
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
                 mass_weighted=False,
                 atoms=None,
                 alphas=None,
                 coord_shape=None,
                 projection_indices=None,
                 transformations=None,
                 internals=None,
                 modes=None,
                 logger=False,
                 optimize_centers=False,
                 clustering_radius=-1,
                 min_singular_value=1e-4,
                 num_svd_vectors=None,
                 svd_contrib_cutoff=1e-3,
                 quadrature_degree=4,
                 expansion_degree=None,
                 expansion_type='multicenter',
                 reference_structure=None,
                 poly_coeffs=None,
                 pairwise_potential_functions=None
    ):

        if internals is not None:
            raise NotImplementedError("internals coming in the future")

        self._S, self._T, self._V = None, None, None
        self._scaling_S = None
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
        if coord_shape is None:
            if centers.ndim > 2:
                coord_shape = centers.shape[1:]
            elif atoms is not None:
                coord_shape = (len(atoms), centers.shape[-1] // len(atoms))
            else:
                coord_shape = (centers.shape[-1],)
        self.coord_shape = coord_shape
        self.inds = projection_indices
        self.mass_weighted = mass_weighted
        # print("MAsSES:", masses)

        if centers.ndim > 2:
            centers = np.reshape(centers, (len(centers), -1))
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
        self._transforms = None
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
                    self.transformations = [t[full_good_pos] for t in self.transformations]

                self._S = None#self._S[full_good_pos, :][:, full_good_pos]
                self._scaling_S = None
                self._T = None#self._T[full_good_pos, :][:, full_good_pos]

        self.logger.log_print("Number of centers: {N}", N=len(self.centers))

        if poly_coeffs is not None:
            poly_coeffs = self.canonicalize_poly_coeffs(poly_coeffs, self.alphas)
        self._poly_coeffs = poly_coeffs

        self.pairwise_potential_functions = pairwise_potential_functions
        self.modes = modes
        self.internals = internals

    def canonicalize_transforms(self, tfs):
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
                    isinstance(x, np.ndarray) and x.shape == (
                    len(self.alphas), self.alphas.shape[-1], self.alphas.shape[-1])
                    for x in tfs
                )
        ):
            raise ValueError("transforms must have shape {s}".format(
                s=(len(self.alphas), self.alphas.shape[-1], self.alphas.shape[-1])
            ))
        return tfs
    @property
    def transformations(self):
        return self._transforms
    @transformations.setter
    def transformations(self, tf):
        self._transforms = self.canonicalize_transforms(tf)

    def _get_hermite_poly(self, coeff_dict, alphas):
        ndim = self.centers.shape[-1]
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

    def canonicalize_poly_coeffs(self, coeffs, alphas):
        if coeffs is None:
            return None

        poly_coeffs = [
            self._get_hermite_poly(c, a) if isinstance(c, dict) else c
            for c, a in zip(coeffs, alphas)
        ]

        return poly_coeffs

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
            ), None, {}
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

    # TODO: add method to rotate to point along bonds
    #       not sure exactly how to define the rest of the coordinate system
    #       but in principle it could be isotropic...

    @classmethod
    def get_virial_alphas(cls, pts, masses,
                          potential_function,
                          remove_translation_rotations=None,
                          translation_rotation_masses=None,
                          planar=None,
                          translation_rotation_frequency=1e-12, # tiny to keep well separated from others...
                          min_frequency=1e-12,
                          allow_rotations=False,
                          reference_modes=None,
                          construct_reference_modes=None,
                          base_rotation=None,
                          scaling=1
                          ):
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
            pots = potential_function(pts)
            grads = potential_function(pts, deriv_order=1)
        else:
            pots, grads, hess = derivs

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

            # in this case I'm _only_ projecting out the gradient for whatever that's worth...
            # which...probably has some rotation/translation component which
            # is why I don't get clean zero eigenvectors...
            rp_mode = grads[non_stationary] / grad_norms[non_stationary][:, np.newaxis]
            num_rp = np.sum(non_stationary.astype(int))
            proj = np.broadcast_to(np.eye(ndim)[np.newaxis], (num_rp, ndim, ndim)) - nput.vec_outer(rp_mode, rp_mode)

            transrot_modes = None
            tr_proj = None
            if remove_translation_rotations is None:
                remove_translation_rotations = True
            if remove_translation_rotations: # assumes pts are flattened molecular structure coords
                # unweight = translation_rotation_masses is not None
                # if not unweight:
                #     translation_rotation_masses = masses
                # trms = translation_rotation_masses
                trms = masses
                structs = None
                if len(trms) % 3 == 0:
                    mass_array = np.reshape(trms, (-1, 3))
                    structs = np.reshape(pts, (-1, mass_array.shape[0], 3))
                    if not np.abs(np.sum(np.diff(mass_array, axis=1))) < 1e-14:  # was actually planar
                        structs = None
                if planar is None:
                    planar = (structs is None and len(trms) >= 2 and trms[0] == trms[1] and len(trms) % 2 == 0)
                if planar: # planar, so we pad
                    mass_array = np.reshape(trms, (-1, 2))
                    structs = np.reshape(pts, (-1, mass_array.shape[0], 2))
                    # if unweight:
                    #     structs = structs / np.sqrt(mass_array)[np.newaxis]
                    structs = np.concatenate([structs, np.zeros((len(structs), mass_array.shape[0], 1))], axis=-1)
                else: #TODO: wat
                    raise ValueError("can't remove translations and rotation eigenvectors for inconsistent structures")

                # check that masses are consistent
                if not np.abs(np.sum(np.diff(mass_array, axis=1))) < 1e-14:  # I guess numerical garbage can happen
                    raise ValueError("can't remove translations with inconsistent mass array {}")

                mol = Molecule(["H"]*mass_array.shape[0], structs[non_stationary], masses=mass_array[:, 0])
                eig, tr = mol.translation_rotation_modes
                if planar:
                    # we now want to remove all eigenvectors involving the Z-coordinate (1 translation & 2 rotations)
                    tr = tr[:, (0, 1, 3, 4, 6, 7), :][:, :, (0, 1, 5)]
                transrot_modes = tr
                # tr_proj = np.broadcast_to(np.eye(ndim)[np.newaxis], (num_rp, ndim, ndim)) - nput.vec_tensordot(tr, tr, shared=1, axes=[2, 2])
                proj = proj - nput.vec_tensordot(tr, tr, shared=1, axes=[2, 2])

            h2 = proj @ hess[non_stationary] @ proj
            freqs, modes = np.linalg.eigh(h2)
            if remove_translation_rotations:
                # we need to forcibly make sure the
                # initial modes correspond to the correct axes
                # which we do by first forcing a sort

                mag_sort = np.argsort(np.abs(freqs), axis=1)
                idx = np.arange(len(freqs)).reshape(-1, 1)
                freqs = freqs[idx, mag_sort]
                modes = modes.transpose(0, 2, 1)[idx, mag_sort].transpose(0, 2, 1)

                ntrot = transrot_modes.shape[-1]
                modes[:, :, :ntrot] = transrot_modes
                modes[:, :, ntrot] = rp_mode
            modes[:, :, -1] = modes[:, :, -1] * np.linalg.det(modes)[:, np.newaxis] # fix inversions

            if remove_translation_rotations:
                # we know the first 3 or 6 frequencies are zeros
                # and after that we have our
                ntrot = transrot_modes.shape[-1]
                freqs[:, :ntrot] = translation_rotation_frequency # some huuuge frequency to make a super tight Gaussian?
                freqs[:, ntrot] = nput.vec_tensordot(
                    rp_mode,
                    nput.vec_tensordot(
                        rp_mode,
                        hess[non_stationary],
                        shared=1,
                        axes=[1, -1]
                    ),
                    shared=1,
                    axes=[1, -1]
                )
            else: # eh we need _something_
                ns_pos = np.where(non_stationary)[0]

                freq_cuts = np.abs(freqs) < 1e-12
                kill_pos = np.where(np.all(freq_cuts, axis=1))
                if len(kill_pos) > 0 and len(kill_pos[0]) > 0:
                    sel = ns_pos[kill_pos]
                    stationary = np.unique(np.concatenate([stationary, sel]))

                zi = np.where(freq_cuts)
                for i, j in zip(*zi):
                    m = modes[i, :, j][:, np.newaxis]
                    f = m.T @ hess[ns_pos[i]] @ m
                    freqs[i, j] = f

            freqs = np.sqrt(np.abs(freqs))

            freqs[freqs < min_frequency] = min_frequency

            # modes = modes / np.sqrt(masses)[np.newaxis, :, np.newaxis] # expressed in terms of mass-weighted coordinates, gotta remove that
            alphas[non_stationary] = scaling * freqs
            rots[non_stationary] = modes

            # rp_coords = rp_mode[:, np.newaxis, np.newaxis, :] @ pts[np.newaxis, :, :, np.newaxis]
            # rp_coords = rp_coords.reshape((len(non_stationary), len(pts)))
            # raise Exception(rp_coords)

            if len(stationary) > 0:
                freqs2, modes = np.linalg.eigh(hess[stationary]) # should automatically identify translations and rotations
                freqs = np.sqrt(np.abs(freqs2))
                freqs[freqs < min_frequency] = min_frequency
                # g = np.diag(1/np.ones(len(masses)))
                # g = modes.transpose(0, 2, 1) @ g[np.newaxis] @ modes
                # rpms = 1 / np.diagonal(g, axis1=1, axis2=2)  # note this is the same as masses...
                alphas[stationary] = scaling * freqs
                # modes[:, :, -1] = modes[:, :, -1] * np.linalg.det(modes)[:, np.newaxis] # fix inversions
                rots[stationary] = modes

            if (
                    reference_modes is None
                    and (construct_reference_modes is None or construct_reference_modes)
            ):
                # we choose a set of reference modes and reexpress all rotations in terms of these
                if len(stationary) > 0:
                    ref_pos = stationary[np.argmin(pots[stationary])]
                else:
                    ref_pos = np.argmin(pots) # no stationary points so...

                reference_modes = rots[ref_pos]

            rots_inv = rots.transpose(0, 2, 1)

            if reference_modes is not None:
                rots = rots@reference_modes[np.newaxis] # express rotations in this new basis
                rots_inv = reference_modes.T[np.newaxis]@rots_inv

            rots = (rots, rots_inv)

            opts = {'mass_weighted':True, 'modes':reference_modes}

        else:

            if remove_translation_rotations: # assumes pts are flattened molecular structure coords
                structs = None
                if len(masses) % 3 == 0:
                    mass_array = np.reshape(masses, (-1, 3))
                    structs = np.reshape(pts, (-1, mass_array.shape[0], 3))
                    if not np.abs(np.sum(np.diff(mass_array, axis=1))) < 1e-14:  # was actually planar
                        structs = None
                planar = (structs is None and len(masses) >= 2 and masses[0] == masses[1] and len(masses) % 2 == 0)
                if planar: # planar, so we pad
                    mass_array = np.reshape(masses, (-1, 2))
                    structs = np.reshape(pts, (-1, mass_array.shape[0], 2))
                    structs = np.concatenate([structs, np.zeros((len(structs), mass_array.shape[0], 1))], axis=-1)
                else:
                    raise ValueError("can't remove translations and rotation eigenvectors for inconsistent structures")

                # check that masses are consistent
                if not np.abs(np.sum(np.diff(mass_array, axis=1))) < 1e-14:  # I guess numerical garbage can happen
                    raise ValueError("can't remove translations with inconsistent mass array {}")

                mol = Molecule(["H"]*mass_array.shape[0], structs, masses=mass_array[:, 0])
                eig, tr = mol.translation_rotation_modes
                if planar:
                    # we now want to remove all eigenvectors involving the Z-coordinate (1 translation & 2 rotations)
                    tr = tr[:, (0, 1, 3, 4, 6, 7), :][:, :, (0, 1, 5)]

                npts = len(pts)
                ndim = pts.shape[-1]

                rm = masses
                hess = hess / np.sqrt(rm[np.newaxis, :, np.newaxis] * rm[np.newaxis, np.newaxis, :])
                tr_proj = np.broadcast_to(np.eye(ndim)[np.newaxis], (npts, ndim, ndim)) - nput.vec_tensordot(tr, tr, shared=1, axes=[2, 2])
                hess = tr_proj @ hess @ tr_proj
                hess *= np.sqrt(rm[np.newaxis, :, np.newaxis] * rm[np.newaxis, np.newaxis, :])

            if base_rotation is not None:
                rots = np.broadcast_to(np.array(base_rotation)[np.newaxis], hess.shape)
                hess = rots.transpose(0, 2, 1) @ hess @ rots
                # with np.printoptions(linewidth=1e8):
                #     raise Exception(str(hess[0]))
            else:
                rots = None

            alphas = scaling * np.sqrt(np.abs(masses * np.diagonal(hess, axis1=1, axis2=2)))

            opts = {}

        return alphas, rots, opts

    @property
    def S(self):
        if self._S is None:
            self._scaling_S, self._S, self._T = self.get_ST()
        return self._S
    @S.setter
    def S(self, smat):
        self._S = smat
    @property
    def T(self):
        if self._T is None:
            self._scaling_S, self._S, self._T = self.get_ST()
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
        else:
            transformations = self.canonicalize_transforms(transformations)
        if transformations is None:
            return None

        if alphas is None:
            alphas = self.alphas

        # if inds is not None:
        #
        #     alphas = alphas
        #     transformations = transformations[:, self.inds]
        #     n = alphas.shape[-1]
        #     npts = len(alphas)
        #     diag_covs = np.zeros((npts, n, n))
        #     diag_inds = (slice(None, None, None),) + np.diag_indices(n)
        #     diag_covs[diag_inds] = 1 / (2 * alphas)
        #     covs = transformations @ diag_covs @ transformations.transpose(0, 2, 1)
        #     covs[np.abs(covs) < 1e-12] = 0  # numerical garbage can be an issue...
        #
        #     two_a_inv, transformations = np.linalg.eigh(covs)
        #     alphas = 1/(2*two_a_inv)

        n = alphas.shape[-1]
        npts = len(alphas)
        diag_covs = np.zeros((npts, n, n))
        diag_inds = (slice(None, None, None),) + np.diag_indices(n)
        diag_covs[diag_inds] = 2*alphas

        if inds is None:
            inds = self.inds
        if inds is not None:
            tfs = transformations[0][:, :, inds]
            inv = transformations[1][:, inds, :]
            diag_covs = diag_covs[:, :, inds][:, inds, :]
            covs = tfs @ diag_covs @ tfs.transpose(0, 2, 1)
        else:
            tfs, inv = transformations
            covs = tfs @ diag_covs @ tfs.transpose((0, 2, 1))

        # if self.inds is not None:
        #     complement = np.setdiff1d(np.arange(self.centers.shape[-1]), self.inds)
        #     covs[np.ix_(complement, self.inds)] = 0
        #     covs[np.ix_(self.inds, complement)] = 0
        # covs[np.abs(covs) < 1e-14] = 0 # numerical garbage can be an issue...

        return covs

    def get_overlap_gaussians(self):
        rows, cols = np.triu_indices(len(self.alphas))

        centers = self.centers
        if self.mass_weighted:
            centers = centers * np.sqrt(self.masses)[np.newaxis, :]

        if self.transformations is None:
            # find overlap gaussians
            alphas = self.alphas
            if self.inds is not None:
                centers = centers[:, self.inds]
                alphas = alphas[:, self.inds]
            new_alphas = alphas[rows] + alphas[cols]
            w_centers = alphas*centers
            # moving weighted average by alpha value
            overlap_data = (w_centers[rows] + w_centers[cols])/new_alphas, new_alphas
        else:
            # if self.inds is not None:
            #     centers = centers[:, self.inds]

            sigs = self.get_inverse_covariances()
            new_sigs = sigs[rows] + sigs[cols]

            if np.allclose(new_sigs, new_sigs.transpose(0, 2, 1)):
                new_alphas, new_rots = np.linalg.eigh(new_sigs) # eigenvalues of inverse tensor...
                new_rots_inv = new_rots.transpose(0, 2, 1)
            else:
                new_alphas, new_rots = np.linalg.eig(new_sigs)  # eigenvalues of inverse tensor...

                comp_part = np.imag(new_alphas)
                if np.sum(comp_part) > 0:
                    raise ValueError(
                        "complex alphas obtained... ({})".format(np.sum(comp_part))
                    )

                new_alphas = np.real(new_alphas)
                new_rots = np.real(new_rots)
                inv_sort = np.argsort(np.abs(new_alphas), axis=-1)
                idx = np.arange(new_alphas.shape[0])[:, np.newaxis]
                new_alphas = new_alphas[idx, inv_sort]

                idx = np.arange(new_alphas.shape[0])[:, np.newaxis, np.newaxis]
                cdx = np.arange(new_alphas.shape[1])[np.newaxis, :, np.newaxis]
                inv_sort = inv_sort[:, np.newaxis, :]

                new_rots = new_rots[idx, cdx, inv_sort]
                new_rots_inv = np.linalg.inv(new_rots)

            # with np.printoptions(linewidth=1e8):
            #     raise Exception(str(
            #         new_alphas[0]
            #         # M @ new_rots[1] @ Minv
            #     ))

            # overlaps = self.transformations[1][rows] @ new_rots
            # sortings = np.argmax(overlaps ** 2, axis=2)  # should test against dupes....
            #
            # idx = np.arange(len(new_alphas)).reshape((-1, 1))
            # new_alphas = new_alphas[idx, sortings]
            # new_rots = new_rots.transpose(0, 2, 1)[idx, sortings].transpose(0, 2, 1)

            if self.inds is not None:
                dets = np.prod(new_alphas[:, self.inds], axis=1)
            else:
                dets = np.prod(new_alphas, axis=1)

            good_pos = np.where(np.abs(dets) > 1e-12)
            new_sigs = new_sigs[good_pos]
            rows = rows[good_pos]
            cols = cols[good_pos]
            new_alphas = new_alphas[good_pos]
            new_rots = new_rots[good_pos]
            new_rots_inv = new_rots_inv[good_pos]

            # if self.inds is not None:
            #     dets = dets[]
            # else:
            #     dets =
            #
            # # dets = np.linalg.det(new_sigs)
            # # # raise Exception(dets)
            # good_pos = np.where(np.abs(dets) > 1e-12)
            # new_sigs = new_sigs[good_pos]
            # rows = rows[good_pos]
            # cols = cols[good_pos]

            # with np.printoptions(linewidth=1e8):
            #     raise Exception(str(new_sigs[0]))
            #     print(np.round(
            #         self.transformations[0].T @
            #                 self.transformations[42],
            #         3))
            #     print(np.round(
            #         self.transformations[0].T @
            #             np.linalg.eigh(new_sigs[0])[1],
            #         3))
            #     raise Exception(
            #         np.linalg.eigh(new_sigs[0])[0],
            #             np.linalg.eigh(new_sigs[1])[0]
            #     )

            # with np.printoptions(linewidth=1e8):
            #     print(new_alphas[0])
            #     print(new_alphas[0] / self.alphas[0])
            #     print(new_alphas[0] / np.linalg.eigh(sigs[0])[0])
            #     print(new_rots[0])
            #     print(self.transformations[0])
            # raise Exception(...)

            if self.inds is not None:
                # build the inverse covariance in the basis of vibrational coordinates
                i = np.array(self.inds)
                d = np.zeros((new_alphas.shape[0], new_alphas.shape[1], new_alphas.shape[1]))
                diag_inds = (
                    np.arange(new_alphas.shape[0])[:, np.newaxis],
                    i[np.newaxis, :],
                    i[np.newaxis, :]
                )
                d[diag_inds] = 1/new_alphas[:, i]
                new_inv = new_rots @ d @ new_rots.transpose(0, 2, 1)
            else:
                # I _could_ construct the inverse from the alphas and rotations
                # but I think it makes more sense to use a potentially more stable
                # inverse here...
                new_inv = np.linalg.inv(new_sigs)

            new_centers = new_inv@(
                sigs[rows] @ centers[rows][:, :, np.newaxis]
                + sigs[cols] @ centers[cols][:, :, np.newaxis]
            )
            new_centers = new_centers.reshape(self.centers[cols].shape)
            new_alphas = new_alphas/2
            sum_sigs = sigs[rows]@new_inv@sigs[cols]

            overlap_data = {
                'row_inds': rows,
                'col_inds': cols,
                'centers': new_centers,
                'alphas': new_alphas,
                'sigmas': new_sigs,
                'rotations': new_rots,
                'inverse_rotations': new_rots_inv,
                'row_sigs': sigs[rows],
                'col_sigs': sigs[cols],
                'sum_inverse':sum_sigs
            }

        return overlap_data

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

    def get_kinetic_polynomials(self):
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

    @classmethod
    def _unrot_base_ST_components(cls, centers, alphas, masses, inds):
        if inds is not None:
            alphas = alphas[:, inds]
            centers = centers[:, inds]
            masses = masses[inds]

        aouter = alphas[:, np.newaxis] * alphas[np.newaxis, :]
        aplus = alphas[:, np.newaxis] + alphas[np.newaxis, :]
        arat = aouter / aplus

        disps = centers[:, np.newaxis, :] - centers[np.newaxis, :, :]

        # A = outer_tet / np.sqrt(np.pi)
        B = np.sqrt(aplus)
        C = arat * np.power(disps, 2)

        return arat, aouter, B, C, masses

    @classmethod
    def _unrot_base_S(cls, centers, alphas, masses, inds):

        arat, aouter, B, C, masses = cls._unrot_base_ST_components(centers, alphas, masses, inds)

        S_dim = (np.sqrt(2) * np.power(aouter, 1 / 4) / B) * np.exp(-C)
        S = np.prod(S_dim, axis=-1)

        return S

    @classmethod
    def _unrot_base_ST(cls, centers, alphas, masses, inds):

        arat, aouter, B, C, masses = cls._unrot_base_ST_components(centers, alphas, masses, inds)

        S_dim = (np.sqrt(2) * np.power(aouter, 1 / 4) / B) * np.exp(-C)
        T_dim = arat * (1 - 2 * C) / masses[np.newaxis, np.newaxis, :]

        S = np.prod(S_dim, axis=-1)
        T = S * np.sum(T_dim, axis=-1)

        return S, T

    @classmethod
    def _rot_base_S_components(cls, rot_data, centers, alphas, inds):
        row_inds = rot_data['row_inds']
        col_inds = rot_data['col_inds']

        # alphas = self.alphas
        # masses = self.masses
        rotas = rot_data["alphas"]
        if inds is not None:
            alphas = alphas[:, inds]
            rotas = rotas[:, inds]

        ndim = centers.shape[-1]
        dets = np.prod(rotas, axis=-1)
        # we prefactor out the 2**ndim
        rdets = np.prod(alphas[row_inds], axis=-1)
        cdets = np.prod(alphas[col_inds], axis=-1)  # literally a product of passed in alphas

        L = rot_data['inverse_rotations']  #
        Lt = rot_data['rotations']
        if inds is not None:
            ndim = len(inds)
            disps = L @ (centers[row_inds] - centers[col_inds])[:, :, np.newaxis]
            disps = np.reshape(disps, disps.shape[:2])[:, inds]
            si = L @ rot_data['sum_inverse'] @ Lt
            si = si[:, inds, :][:, :, inds]
            C = disps[:, np.newaxis, :] @ si @ disps[:, :, np.newaxis]
        else:
            disps = centers[row_inds] - centers[col_inds]
            C = disps[:, np.newaxis, :] @ rot_data['sum_inverse'] @ disps[:, :, np.newaxis]
        C = C.reshape(disps.shape[0])

        return alphas, rotas, L, Lt, row_inds, col_inds, ndim, rdets, cdets, dets, C

    @classmethod
    def _rot_base_S(cls, rot_data, centers, alphas, inds):

        S = np.eye(len(centers))
        alphas, rotas, L, Lt, row_inds, col_inds, ndim, rdets, cdets, dets, C = \
            cls._rot_base_S_components(rot_data, centers, alphas, inds)

        S[row_inds, col_inds] = S[col_inds, row_inds] = (
                2 ** (ndim / 2) * ((rdets * cdets) / (dets ** 2)) ** (1 / 4) * np.exp(-C / 2)
        )
        return S

    @classmethod
    def _rot_base_ST(cls, rot_data, centers, alphas, masses, inds):

        S = np.eye(len(centers))
        alphas, rotas, L, Lt, row_inds, col_inds, ndim, rdets, cdets, dets, C = \
            cls._rot_base_S_components(rot_data, centers, alphas, inds)

        S[row_inds, col_inds] = S[col_inds, row_inds] = (
                2 ** (ndim / 2) * ((rdets * cdets) / (dets ** 2)) ** (1 / 4) * np.exp(-C / 2)
        )

        T = np.zeros((len(centers), len(centers)))

        rows = rot_data['row_sigs']
        cols = rot_data['col_sigs']

        idx = col_inds
        Sj = cols

        zetas = L @ (rot_data['centers'] - centers[idx])[:, :, np.newaxis]
        zetas = zetas.reshape(rot_data['centers'].shape)
        if inds is not None:
            zetas = zetas[:, inds]

        minv = np.diag(1 / masses)
        Amat = L @ (Sj @ minv @ Sj) @ Lt
        M = L @ np.broadcast_to(minv[np.newaxis], L.shape) @ Lt
        # msj = rot_data['rotations']@(Sj@minv)@np.transpose(rot_data['rotations'], (0, 1, 2))
        # raise Exception(
        #     np.diagonal(msj, axis1=1, axis2=2),# / self.masses[np.newaxis],
        #     - 1 / 2 * np.diagonal(Amat, axis1=1, axis2=2) / rot_data['alphas']
        # )
        submasses = np.diagonal(M, axis1=1, axis2=2)
        sminv = np.diagonal(Amat, axis1=1, axis2=2)
        if inds is not None:
            submasses = submasses[:, inds]
            sminv = sminv[:, inds]
            Amat = Amat[:, inds, :][:, :, inds]
        # raise Exception(
        #     submasses,
        #     sminv,
        #     np.sum(
        #         2 * alphas[col_inds] * submasses -
        #         1 / 2 * sminv / rotas,
        #         axis=1
        #     )
        # )
        T[row_inds, col_inds] = 1 / 2 * (
                np.sum(
                    2 * alphas[idx] * submasses - 1 / 2 * sminv / rotas,
                    axis=1
                ) - 1 / 2 * sum(
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
        T[col_inds, row_inds] = T[row_inds, col_inds]

        return S, T

    @classmethod
    def _rot_poly_ST(cls, rot_data, base_polys, centers, alphas, masses, inds):
        raise ValueError("need to account for Coriolis effects...")

        # S = np.eye(len(centers))
        alphas, rotas, L, Lt, row_inds, col_inds, ndim, rdets, cdets, dets, C = \
            cls._rot_base_S_components(rot_data, centers, alphas, inds)

        # S[row_inds, col_inds] = S[col_inds, row_inds] = (
        #         2 ** (ndim / 2) * ((rdets * cdets) / (dets ** 2)) ** (1 / 4) * np.exp(-C / 2)
        # )

        T = np.zeros((len(centers), len(centers)))

        rows = rot_data['row_sigs']
        cols = rot_data['col_sigs']

        idx = col_inds
        Sj = cols

        tf_centers = L @ centers[idx][:, :, np.newaxis]
        tf_centers = tf_centers[:, :, inds]

        tf_sigs = L @ (Sj @ np.transpose(L, (0, 2, 1)))
        tf_sigs = tf[:, inds, :][:, :, inds]

        exp_polys = [
            DensePolynomial.from_tensors([0, np.zeros(ndim), sig]).shift(c)
            # , shift=-pt)
            for sig, c in zip(tf_sigs, tf_centers)
            # for pt, sig in zip(
            #     centers,
            #     sigs
            # )
        ]




        raise Exception('need to apply indices before floop')

        exp_d1 = [
            d.der
        ]

        zetas = L @ (rot_data['centers'] - centers[idx])[:, :, np.newaxis]
        zetas = zetas.reshape(rot_data['centers'].shape)
        if inds is not None:
            zetas = zetas[:, inds]

        minv = np.diag(1 / masses)
        Amat = L @ (Sj @ minv @ Sj) @ Lt
        M = L @ np.broadcast_to(minv[np.newaxis], L.shape) @ Lt
        # msj = rot_data['rotations']@(Sj@minv)@np.transpose(rot_data['rotations'], (0, 1, 2))
        # raise Exception(
        #     np.diagonal(msj, axis1=1, axis2=2),# / self.masses[np.newaxis],
        #     - 1 / 2 * np.diagonal(Amat, axis1=1, axis2=2) / rot_data['alphas']
        # )
        submasses = np.diagonal(M, axis1=1, axis2=2)
        sminv = np.diagonal(Amat, axis1=1, axis2=2)
        if inds is not None:
            submasses = submasses[:, inds]
            sminv = sminv[:, inds]
            Amat = Amat[:, inds, :][:, :, inds]
        # raise Exception(
        #     submasses,
        #     sminv,
        #     np.sum(
        #         2 * alphas[col_inds] * submasses -
        #         1 / 2 * sminv / rotas,
        #         axis=1
        #     )
        # )
        T[row_inds, col_inds] = 1 / 2 * (
                np.sum(
                    2 * alphas[idx] * submasses - 1 / 2 * sminv / rotas,
                    axis=1
                ) - 1 / 2 * sum(
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
        T[col_inds, row_inds] = T[row_inds, col_inds]

        return S, T

    def get_ST(self, centers=None, alphas=None, transformations=None):
        if centers is None:
            centers = self.centers
        if alphas is None:
            alphas = self.alphas
        if transformations is None:
            transformations = self.transformations

        masses = self.masses
        if self.mass_weighted:
            centers = centers * np.sqrt(masses)[np.newaxis]
            masses = np.ones(len(masses))

        if self._poly_coeffs is not None:

            rot_data = self.get_overlap_gaussians()
            base_polys, ke_polys = self.get_kinetic_polynomials() # polys for each dimension + other stuff

            if not isinstance(rot_data, dict):
                ov_centers, ov_alphas = rot_data
                rot_data = None
                row_inds, col_inds = np.triu_indices(len(self.centers))
            else:
                ov_centers = rot_data['centers']
                ov_alphas = rot_data['alphas']
                row_inds = rot_data['row_inds']
                col_inds = rot_data['col_inds']

            if rot_data is None:
                S = self._unrot_base_S(centers, alphas, masses, self.inds)
            else:
                S = self._rot_base_S(rot_data, centers, alphas, self.inds)

            ndim = self.centers.shape[-1]
            if rot_data is None and self.inds is not None:
                ndim = len(self.inds)

            # raise Exception(S)

            overlap_polys = [
                base_polys[r] * base_polys[c]
                for r, c in zip(row_inds, col_inds)
            ]
            overlap_polys = [
                p.shift(c) if isinstance(p, DensePolynomial) else p
                for p, c in zip(overlap_polys, ov_centers)
            ]
            # print(overlap_polys[0].coeffs)

            subtensors = [
                p.clip(threshold=1e-8).coefficient_tensors
                if isinstance(p, DensePolynomial) else None
                for p in overlap_polys
            ]

            # raise Exception(subtensors[0])
            i = 1
            # raise Exception([
            #         tt[i]
            #             if isinstance(tt, list) else
            #         np.zeros((nd,)*i)
            #             if i > 0 else
            #         1
            #         for tt in subtensors
            #     ])
            order = max(tuple(len(tt) for tt in subtensors if tt is not None) + (1,))
            derivs = [
                np.array([
                    tt[i]
                        if isinstance(tt, list) and i < len(tt) else
                    np.zeros((ndim,)*i)
                        if i > 0 else
                    1
                        if tt is None else
                    tt
                    for tt in subtensors
                ])
                for i in range(order)
            ]
            # raise Exception(derivs[2][0])
            #
            S_extra = self.tensor_expansion_integrate(
                    len(self.centers),
                    derivs,
                    ov_centers,
                    ov_alphas,
                    inds=None if rot_data is None else self.inds,
                    rot_data=rot_data,
                    expansion_type='multicenter',
                    logger=None,#self.logger
                    reweight=False
                )

            S_base = S
            S = S * S_extra

            if self.transformations is None and self.inds is None:
                T = np.zeros((len(self.centers), len(self.centers)))
                idxs = range(ndim) if rot_data is not None else self.inds
                for coord in idxs:

                    kinetic_polys = [
                        base_polys[r] * ke_polys[c][coord]
                        for r, c in zip(row_inds, col_inds)
                    ]
                    kinetic_polys = [
                        p.shift(c).clip(threshold=1e-4) if isinstance(p, DensePolynomial) else p
                        for p, c in zip(kinetic_polys, ov_centers)
                    ]
                    subtensors = [
                        p.coefficient_tensors if isinstance(p, DensePolynomial) else p
                        for p in kinetic_polys
                    ]

                    print("?", coord, subtensors[2])

                    order = max([
                                    len(tt) for tt in subtensors
                                    if tt is not None and not isinstance(tt, (int, float, np.integer, np.floating))
                                ] + [0])
                    derivs = [
                        np.array([
                            tt[i]
                                if isinstance(tt, list) and i < len(tt) else
                            np.zeros((ndim,) * i)
                                if i > 0 else
                            1
                                if tt is None else
                            tt
                            for tt in subtensors
                        ]
                        )
                        for i in range(order)
                    ]
                    if len(derivs) > 0:
                        contrib = self.tensor_expansion_integrate(
                            len(self.centers),
                            derivs,
                            ov_centers,
                            ov_alphas,
                            inds=None if rot_data is None else self.inds,
                            rot_data=rot_data,
                            expansion_type='multicenter',
                            logger=None,#self.logger
                            reweight=False
                        ) / self.masses[coord]

                        T = T + contrib

                    # raise Exception(T*S_base)

                T *= -1/2 * S_base
            else:
                T = self._rot_poly_ST(rot_data, base_polys, centers, alphas, masses, self.inds)


                raise Exception(...)

        else:

            if transformations is None:
                S, T = self._unrot_base_ST(centers, alphas, masses, self.inds)
            else:
                rot_data = self.get_overlap_gaussians()
                S, T = self._rot_base_ST(rot_data, centers, alphas, masses, self.inds)

            S_base = S

        return S_base, S, T

    @staticmethod
    def quad_nd(centers, alphas, function, degree=3, normalize=True):
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
        # disps, weights = np.polynomial.hermite.hermgauss(degree)
        centers, alphas = self.get_overlap_gaussians() # Only upper triangle here

        npts = len(self.alphas)
        pots = np.zeros((npts, npts))
        rows, cols = np.triu_indices(npts)
        #
        vals = self.quad_nd(centers, alphas, function, degree=degree, normalize=False)
        pots[rows, cols] = vals
        pots[cols, rows] = vals

        # pots = np.zeros((npts, npts))
        #
        # ndim = centers.shape[-1]
        # for disp_inds in itertools.product(*([range(degree)] * ndim)):
        #     disp_inds = np.array(disp_inds)
        #     w = np.prod(weights[disp_inds])
        #     c = centers + disps[disp_inds][np.newaxis, :] / np.sqrt(alphas)
        #     print(w, function(c))
        #     pots[rows, cols] += w * function(c)
        # pots[cols, rows] = pots[rows, cols]
        # raise Exception(pots[0, 0]/pots2[0, 0])

        normalization = 1 / (np.sqrt(np.pi)) ** self.centers.shape[-1]
        pots *= normalization

        return pots

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

    @staticmethod
    def mass_weighted_eval(function, coords, masses, deriv_order=None, **opts):
        # expects just x and y coordinates for the atoms
        mass_vec = np.asanyarray(masses)

        coords = coords / np.sqrt(mass_vec[np.newaxis])
        coords = coords.reshape(-1, 3, 2)
        derivs = function(coords, deriv_order=deriv_order, **opts)
        new_d = []
        for n, do in enumerate(derivs):
            if n > 0:
                mv = np.sqrt(mass_vec)[np.newaxis]
                weight = mv
                for j in range(n - 1):
                    mv = np.expand_dims(mv, 0)
                    weight = np.expand_dims(weight, -1) * mv
                # with np.printoptions(linewidth=1e8):
                #     print(weight[0])
                new_d.append(do / weight)
            else:
                new_d.append(do)
        derivs = new_d
        return derivs

    @classmethod
    def tensor_expansion_integrate(cls,
                                   npts, derivs, centers, alphas,
                                   inds=None, rot_data=None, expansion_type='multicenter',
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

        ndim = centers.shape[-1]
        if rot_data is not None:
            rotations = rot_data['rotations']
            centers = np.reshape(rotations @ centers[:, :, np.newaxis], centers.shape)

            if inds is not None:
                ndim = len(inds)
                rotations = rotations[:, :, inds]
                alphas = alphas[..., inds]
                centers = centers[..., inds]

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

        elif inds is not None:
            ndim = len(inds)
            if alphas.shape[-1] > ndim:
                alphas = alphas[..., inds]
            if centers.shape[-1] > ndim:
                centers = centers[..., inds]

            new_derivs = []
            for n, d in enumerate(derivs):
                for k in range(n):
                    d = np.take(d, inds, axis=k+1)
                new_derivs.append(d)
            derivs = new_derivs

        if logger is not None:
            logger.log_print("adding up all derivative contributions...")

        if rot_data is None:
            row_inds, col_inds = np.triu_indices(npts)
        else:
            row_inds = rot_data['row_inds']
            col_inds = rot_data['col_inds']

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

    def expansion_integrate(self, function, deriv_order=2,
                            expansion_type=None,
                            pairwise_functions=None,
                            pairwise_quadrature_degree=None
                            ):
        if expansion_type is None:
            expansion_type = self.expansion_type

        if self.transformations is None:
            centers, alphas = self.get_overlap_gaussians()
            rot_data = None
        else:
            if pairwise_functions is not None:
                raise ValueError("currently can't have pairwise functions & rotations")
            rot_data = self.get_overlap_gaussians()
            centers = rot_data['centers']
            alphas = rot_data['alphas']

        if pairwise_functions is not None:

            pot_contribs, deriv_corrs = self.integrate_pairwise_potential_contrib(
                pairwise_functions,
                centers,
                alphas,
                coord_shape=self.coord_shape,
                quadrature_degree=pairwise_quadrature_degree,
                expansion_degree=deriv_order
            )

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
            if self.mass_weighted:
                derivs = self.mass_weighted_eval(function, centers, self.masses, deriv_order=deriv_order)
            else:
                derivs = function(centers, deriv_order=deriv_order)
                if isinstance(derivs, np.ndarray):  # didn't get the full list so we do the less efficient route'
                    derivs = [function(centers)] + [
                        function(centers, deriv_order=d) for d in range(1, deriv_order+1)
                    ]


        pot = self.tensor_expansion_integrate(
            len(self.centers),
            derivs,
            centers,
            alphas,
            inds=self.inds,
            rot_data=rot_data,
            expansion_type=expansion_type,
            logger=self.logger
        )

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

    def integrate_pairwise_potential_contrib(self,
                                             functions,
                                             centers,
                                             alphas,
                                             coord_shape=None,
                                             quadrature_degree=4,
                                             expansion_degree=4
                                             ):

        # if self.transformations is None:
        #     centers, alphas = self.get_overlap_gaussians()
        #     rot_data = None
        # else:
        #     rot_data = self.get_overlap_gaussians()
        #     centers = rot_data['centers']
        #     alphas = rot_data['alphas']
        if coord_shape is not None and len(coord_shape) == 1:
            raise ValueError("can't determine if coords are planar or not")

        planar = coord_shape is not None and coord_shape[-1] == 2

        npts = centers.shape[0]
        ndim = 2 if planar else 3
        sdim = ndim // 2
        centers = centers.reshape(npts, -1, ndim)
        alphas = alphas.reshape(npts, -1, ndim)

        potential_contrib = 0
        if expansion_degree is not None:
            derivative_contribs = [
                np.zeros((npts,) + (centers.shape[-1],)*i)
                    if i > 0 else
                0
                for i in range(expansion_degree+1)
            ]
        else:
            derivative_contribs = [0]
        for index_pair, f in functions.items():
            subcenters = centers[:, index_pair, :].reshape(npts, 2, ndim)
            subalphas = alphas[:, index_pair, :].reshape(npts, 2, ndim)

            if planar:
                tf = np.array([
                    [1, 0, -1, 0],  # 0,  0],
                    [0, 1,  0, -1],  # 0,  0],
                    # [0,  0, 0,  0, 1, -1],
                    [1, 0, 1, 0],  # 0,  0],
                    [0, 1, 0, 1],  # 0,  0],
                    # [0,  0, 0,  0, 1,  1]
                ]) #/ np.sqrt(2)
                tf = np.broadcast_to(tf[np.newaxis], (npts, 4, 4))
                tf[:, 2, 2] = subalphas[:, 0, 0] / subalphas[:, 1, 0]
                tf[:, 3, 3] = subalphas[:, 0, 1] / subalphas[:, 1, 1]
            else:
                tf = np.array([
                    [1, 0, 0, -1, 0, 0],
                    [0, 1, 0,  0, -1, 0],
                    [0, 0, 1,  0, 0, -1],
                    [1, 0, 0,  1, 0, 0],
                    [0, 1, 0,  0, 1, 0],
                    [0, 0, 1,  0, 0, 1]
                ]) #/ np.sqrt(2)
                tf = np.broadcast_to(tf[np.newaxis], (npts, 6, 6))
                tf[:, 3, 3] = subalphas[:, 0, 0] / subalphas[:, 1, 0]
                tf[:, 4, 4] = subalphas[:, 0, 1] / subalphas[:, 1, 1]
                tf[:, 5, 5] = subalphas[:, 0, 2] / subalphas[:, 1, 2]

            subcenters = centers[:, index_pair, :].reshape(npts, 2 * ndim)
            subalphas = alphas[:, index_pair, :].reshape(npts, 2 * ndim)

            subcov = np.zeros((subalphas.shape[0], ndim, ndim))
            # fill diagonals across subcov
            diag_inds = (slice(None, None, None),) + np.diag_indices(ndim)
            subcov[diag_inds] = 2 * subalphas
            tf_centers = subcenters[:, np.newaxis, :]@tf
            tf_alphas = (tf@subcov@tf.T)[diag_inds] / 2
            # a2 = np.diag(np.dot(np.dot(tf, np.diag(2 * gauss_alpha)), tf.T)) / 2
            pairwise_contrib = self.quad_nd(
                tf_centers[:, :sdim],
                tf_alphas[:, :sdim],
                f,
                degree=quadrature_degree
            ) * (
                    np.sqrt(np.pi) ** sdim /
                        np.prod(np.sqrt(tf_alphas[:, sdim:]), axis=-1)
            ) #TODO: do I need to multiply by determinant????
            potential_contrib += pairwise_contrib
            if expansion_degree is not None:
                # need to map index_pair (which correspond to 4 or 6 inds) to the
                # correct list of inds in full space
                idx_pair_full_inds = sum(
                    (tuple(range(sdim*x, sdim*(x+1))) for x in index_pair),
                    ()
                )
                deriv_contrib = f(tf_centers[:, :sdim], deriv_order=expansion_degree)
                for i,d in deriv_contrib:
                    if i == 0:
                        derivative_contribs[0] += d
                    else:
                        subd_contrib = np.zeros((npts,) + (ndim,) * i)
                        subd_inds = (slice(None, None, None),) + (slice(None, sdim, None),)*i
                        subd_contrib[subd_inds] = d
                        for j in range(i):
                            subd_contrib = np.tensordot(subd_contrib, tf, axes=[j+1, 0])
                        dc_inds = (slice(None, None, None),) + (idx_pair_full_inds,)*i
                        derivative_contribs[i][dc_inds] += subd_contrib
            else:
                derivative_contribs[0] += f(tf_centers[:, :sdim])

        return potential_contrib, derivative_contribs

    class PairwisePotentialFunction:
        """
        A wrapper for a base function with derivatives that allows it to be evaluate
        """

        def __init__(self, pot_fun):
            self.pot = pot_fun
        def eval(self, delta_values, deriv_order=None):
            """
            We assume the delta_values can be used to directly construct the distances by taking their norms squared

            :param delta_values:
            :return:
            """

            if deriv_order > 2:
                raise NotImplementedError("need to use FD on the second derivs...")
            r_derivs = nput.vec_norm_derivs(delta_values, order=deriv_order)
            r_vals = r_derivs[0]
            r_derivs = r_derivs[1:]
            fun_derivs = self.pot(r_vals, deriv_order=deriv_order)
            # need to do tensor derivative conversion on fun_derivs
            fun_vals = fun_derivs[0]
            fun_derivs = fun_derivs[1:]

            fun_derivs = TensorDerivativeConverter(r_derivs, fun_derivs).convert()

            return fun_vals, fun_derivs

    class PairwiseMorsePotential:
        def __init__(self, re, alpha, de):
            self.pot = ... #TODO: load from data source



    # def wrap_pairwise_distance_function(self, ...):
    #     """
    #     Want to wrap have a way to wrap a function (with derivatives) that only
    #     depends on the norm of a transformed set of those coordinates (which will be the distance between two atoms)
    #     and have it return both the function value and, if requested, the derivatives evaluated at
    #     those points
    #
    #     :return:
    #     """
    #     ...

    def evaluate_multiplicative_operator_base(self,
                                              function,
                                              handler=None,
                                              expansion_degree=None,
                                              expansion_type=None,
                                              quadrature_degree=None,
                                              pairwise_functions=None
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
            if pairwise_functions is not None:
                raise ValueError("pairwise functions can't go with direct quadrature")
            self.logger.log_print("evauating integrals with {n}-order quadrature", n=self.quadrature_degree)
            pot_mat = self.quad_integrate(function, degree=self.quadrature_degree if quadrature_degree is None else quadrature_degree)
        elif handler == 'expansion':
            self.logger.log_print("evauating integrals with {n}-degree expansions", n=self.expansion_degree)
            pot_mat = self.expansion_integrate(function,
                                               deriv_order=expansion_degree,
                                               expansion_type=expansion_type,
                                               pairwise_functions=pairwise_functions
                                               )
        elif handler == 'analytic':
            self.logger.log_print("evauating integrals analytically", n=self.expansion_degree)
            pot_mat = self.analytic_integrate()
        else:
            raise ValueError("unknown operator evaluation scheme {}".format(handler))

        return pot_mat

    def evaluate_multiplicative_operator(self,
                                         function,
                                         handler=None,
                                         expansion_degree=None,
                                         expansion_type=None,
                                         quadrature_degree=None,
                                         pairwise_functions=None # integrate out pairwise contribution
                                         ):

        pot_mat = self.evaluate_multiplicative_operator_base(
            function,
            handler=handler,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            quadrature_degree=quadrature_degree,
            pairwise_functions=pairwise_functions
        )
        if self._scaling_S is None:
            _ = self.S
        S = self._scaling_S
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

        # we'll ignore anything too large?
        # and anything negative?
        # mostly to remove the troubles from rotations and translations...
        base_loc = np.where(np.logical_and(sig > 0, np.abs(sig) < 1e30))
        sig = sig[base_loc]

        # sort = np.argsort(np.abs(sig))
        # sig = sig[sort]
        # evecs = evecs[:, sort]

        if subspace_size is not None:
            good_loc = slice(max(0, len(sig) - subspace_size), len(sig))
        elif min_singular_value is not None:
            good_loc = np.where(sig > min_singular_value)[0]
        else:
            good_loc = np.arange(len(S))
        d = np.diag(1 / np.sqrt(sig[good_loc]))

        gl = base_loc[0][good_loc]
        L = evecs[:, gl]
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
                    eps=5e-4,
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
                        nodeless_ground_state=Qq.shape[1] > 1 # gotta bottom out some time...
                    )

        elif mode == 'stable': # Implementation of the Fix-Heiberger algorithm
            # raise NotImplementedError("I think my Fix-Heiberger is broken?"
            #                           "And also potentially not even the right algorithm for this problem"
            #                           "Like I potentially _want_ a large epsilon since I want to get the best"
            #                           "possible energies, but Fix-Heiberger will call that unstable I think")

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
            # B0 = np.diag(np.concatenate([
            #     d[g1],
            #     np.zeros(n2)
            # ]))
            R = np.diag(np.concatenate([
                1 / np.sqrt(d[g1]),
                np.ones(n2)
            ]))
            # with np.printoptions(linewidth=1e8):
            #     raise Exception(R@B0@R)
            A1 = R.T@Q.T@H@Q@R

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

    def get_wavefunctions(self, print_debug_info=False, min_singular_value=None, subspace_size=None, nodeless_ground_state=None,
                          mode=None, stable_epsilon=None
                          ):
        # print("======="*25)
        ops = dict(
            print_debug_info=print_debug_info,
            min_singular_value=min_singular_value,
            nodeless_ground_state=nodeless_ground_state,
            subspace_size=subspace_size,
            mode=mode,
            eps=stable_epsilon
        )
        ops = {k:v for k,v in ops.items() if v is not None}
        eigs, evecs = self.diagonalize(**ops)
        return DGBWavefunctions(eigs, evecs, hamiltonian=self)