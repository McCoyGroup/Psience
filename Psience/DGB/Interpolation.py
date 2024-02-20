
import abc

import numpy as np, scipy as sp

from McUtils.Zachary import InverseDistanceWeightedInterpolator, TensorDerivativeConverter
import McUtils.Numputils as nput

__all__ = [
    "DGBInterpolator",
    "DGBGenericInterpolator",
    "DGBWatsonInterpolator",
    "DGBCartesianWatsonInterpolator",
]

class DGBInterpolator(metaclass=abc.ABCMeta):
    default_neighborhood_size=10
    default_declustering_alpha=0
    default_declustering_overlap=.1
    def __init__(self, centers, potential_derivatives,
                 declustering_alpha=None,
                 declustering_overlap=None,
                 neighborhood_size=None,
                 logger=None,
                 pairwise_potential_functions=None,
                 **opts
                 ):
        self.centers = centers
        self.derivs = potential_derivatives
        if neighborhood_size is None:
            neighborhood_size = self.default_neighborhood_size
        neighborhood_size = min(neighborhood_size, len(centers))
        if declustering_alpha is None:
            declustering_alpha = self.default_declustering_alpha
        if declustering_overlap is None:
            declustering_overlap = self.default_declustering_overlap
        self.declustering_opts = dict(
            declustering_overlap=declustering_overlap,
            declustering_alpha=declustering_alpha
        )
        self.ppfs = pairwise_potential_functions
        self.logger = logger
        self.opts = dict(opts,
                         pairwise_potential_functions=pairwise_potential_functions,
                         neighborhood_size=neighborhood_size
                         )
    @abc.abstractmethod
    def __call__(self, points, deriv_order=None, **kwargs) -> 'np.ndarray | list[np.ndarray]':
        ...
class DGBGenericInterpolator(DGBInterpolator):

    def __init__(self, centers, potential_derivatives, **opts):
        super().__init__(centers, potential_derivatives, **opts)
        interp_opts = self.opts.copy()
        del interp_opts['pairwise_potential_functions']
        self.interp = InverseDistanceWeightedInterpolator(
            self.centers,
            *self.derivs,
            **interp_opts
        )
    def __call__(self, points, deriv_order=None, **kwargs):
        just_vals = deriv_order is None
        if deriv_order is None: deriv_order = 0
        hess_interp = self.interp(points.reshape(points.shape[0], -1), deriv_order=deriv_order, **kwargs)
        if just_vals:
            return hess_interp[0]
        else:
            return hess_interp

class DGBWatsonInterpolatorSingleRef(DGBInterpolator):
    def __init__(self, centers, potential_derivatives, modes, **opts):
        super().__init__(centers, potential_derivatives, **opts)
        self.modes = modes
        self.adjusted_derivs = self.take_remainder_potential(centers, potential_derivatives, modes)
        self.interp = InverseDistanceWeightedInterpolator(
                self.centers,
                *self.adjusted_derivs,
                **self.opts
            )
    def take_remainder_potential(self, centers, potential_derivatives, modes):
        freqs = modes.freqs
        harmonic_contribs = np.dot(centers**2, freqs**2) / 2
        harmonic_deriv = np.diag(freqs**2)
        remainder_vals = potential_derivatives[0] - harmonic_contribs
        remainder_hess = potential_derivatives[2] - harmonic_deriv[np.newaxis, :, :]
        return [remainder_vals, potential_derivatives[1], remainder_hess] + potential_derivatives[3:]
    def __call__(self, points, deriv_order=None, **kwargs):
        just_vals = deriv_order is None
        if deriv_order is None: deriv_order = 0
        interped_vals = self.interp(points.reshape(points.shape[0], -1), deriv_order=deriv_order, **kwargs)
        if isinstance(interped_vals, np.ndarray):
            interped_vals = [interped_vals]
        freqs = self.modes.freqs
        harmonic_contribs = np.dot(points**2, freqs**2) / 2
        harmonic_deriv = np.diag(freqs ** 2)
        base_vals = harmonic_contribs + interped_vals[0]
        if just_vals:
            return base_vals
        interped_vals[0] = base_vals
        if deriv_order > 1:
            interped_vals[2] = harmonic_deriv[np.newaxis, :, :] + interped_vals[2]
            # interped_vals[2] = np.broadcast_to(harmonic_deriv[np.newaxis, :, :], interped_vals[2].shape)
        return interped_vals

class WatsonPairwisePotential:
    """
    Partial mimic of the PairwisePotentialEvaluator, but just taking
    the bits needed to actually convert f(q) -> f(delta) -> f(r) and get the
    derivatives cleanly
    """
    def __init__(self, ppfs:'dict[(int, int), callable]', modes):
        if 'functions' in ppfs: ppfs = ppfs['functions']
        self.ppfs = ppfs
        self.modes = modes

    @classmethod
    def get_bond_length_deltas(cls, natoms, ndim, i, j, full=False):
        if not full:
            mat = np.zeros((ndim, natoms, ndim))
            for k in range(ndim):
                mat[k][i][k] = 1
                mat[k][j][k] = -1
            mat = mat.reshape((mat.shape[0], natoms * ndim))
        else:
            mat = np.sqrt(2) * np.eye(natoms * ndim).reshape((natoms, ndim, natoms, ndim))
            for k in range(ndim):
                mat[i][k][i][k] = 1
                mat[i][k][j][k] = -1
                mat[j][k][i][k] = 1
                mat[j][k][j][k] = 1
            mat = mat.reshape((natoms * ndim, natoms * ndim))
        return mat

    def get_coordinate_bond_length_projection(self, i, j, ndim=3):# TODO: relax this
        modes = self.modes.matrix
        natoms = modes.shape[0] // ndim
        base_mat = self.get_bond_length_deltas(natoms, ndim, i, j)
        tf_base = base_mat @ modes
        d0 = np.dot(base_mat, self.modes.origin.reshape(-1))
        return d0, tf_base

    def get_bond_lengths(self, coords, i, j, deriv_order=0):
        d0, proj = self.get_coordinate_bond_length_projection(i, j)
        delta_contribs = proj[np.newaxis] @ coords[:, :, np.newaxis]
        delta_contribs = np.reshape(delta_contribs, delta_contribs.shape[:2])

        d = d0[np.newaxis, :] + delta_contribs

        if deriv_order is None or deriv_order == 0:
            r = np.linalg.norm(d, axis=-1)  # we ensure no scaling of difference coords
            r_derivs = None
        else:
            r_derivs = nput.vec_norm_derivs(d, order=deriv_order)
            r, r_derivs = r_derivs[0], r_derivs[1:]

        if r_derivs is not None and len(r_derivs) > 0:
            q_derivs = []
            for i, d in enumerate(r_derivs):
                for j in range(d.ndim - 1):
                    d = np.tensordot(d, proj, axes=[1, 0])
                    # d = nput.vec_tensordot(d, tf_proj, shared=1, axes=[1, 1])
                q_derivs.append(d)

            return [r] + q_derivs
        else:
            return r

    def eval_ppf(self, f, coords, i, j, deriv_order=None):
        rdat = self.get_bond_lengths(coords, i, j, deriv_order=deriv_order)
        if deriv_order is not None:
            r = rdat[0]
        else:
            r = rdat

        fdat = f(r, deriv_order=deriv_order)

        if deriv_order is not None:
            f_vals = fdat[0]
            r_derivs = rdat[1:]
            r_derivs = [r[..., np.newaxis] for r in r_derivs]
            f_derivs = [fdat[i].reshape(f_vals.shape + (1,) * i) for i in range(1, deriv_order + 1)]
            f_derivs = list(TensorDerivativeConverter(r_derivs, f_derivs).convert())
            fdat = [f_vals] + f_derivs

        return fdat

    def __call__(self, normal_mode_coords, deriv_order=None):
        vals = None
        for inds, f in self.ppfs.items():
            if len(inds) > 2:
                raise NotImplementedError("only stretches currently supported")
            fdat = self.eval_ppf(f, normal_mode_coords, deriv_order=deriv_order, *inds)
            if vals is None:
                vals = fdat
            else:
                if deriv_order is not None:
                    vals = [v + f for v,f in zip(vals, fdat)]
                else:
                    vals += fdat
        return vals

class DGBWatsonTaylorInterpolator(DGBInterpolator):
    def __init__(self, centers, potential_derivatives, modes,
                 power=4,
                 include_harmonic_basis=False, # provably does nothing for now...
                 harmonic_distance_cutoff=None, # just to make plots look better...
                 pairwise_potential_functions=None,
                 **opts
                 ):

        super().__init__(centers, potential_derivatives, power=power,
                         pairwise_potential_functions=pairwise_potential_functions,
                         **opts)

        self.modes = modes
        if (
                pairwise_potential_functions is not None and
                pairwise_potential_functions.get('use_with_interpolation', True)
        ):
            pairwise_potential_functions = WatsonPairwisePotential(pairwise_potential_functions, modes)
        else:
            pairwise_potential_functions = None
        self.ppfs = pairwise_potential_functions

        alpha = self.declustering_opts['declustering_alpha']
        if alpha > 0: # I really don't think we need this...
            from .Evaluators import DGBEvaluator
            from .Gaussians import DGBGaussians
            overlap_data = DGBEvaluator.get_overlap_gaussians(
                    self.centers,
                    np.full(self.centers.shape, alpha),
                    None
                )
            S = DGBEvaluator.evaluate_overlap(overlap_data)
            inds = DGBGaussians._optimize_gs_block(S, self.declustering_opts['declustering_overlap'])
        else:
            inds = slice(None)

        if pairwise_potential_functions is not None:
            self.adjusted_derivs = self.take_ppf_remainder(centers, potential_derivatives)
        elif include_harmonic_basis:
            self.adjusted_derivs = self.take_remainder_potential(centers, potential_derivatives, modes)
        else:
            self.adjusted_derivs = self.derivs

        # if include_harmonic_basis:
        #     self.adjusted_derivs = self.take_remainder_potential(centers, potential_derivatives, modes)
        # else:
        #     self.adjusted_derivs = self.derivs

        self.include_harmonic_basis = include_harmonic_basis
        self.harmonic_distance_cutoff = harmonic_distance_cutoff
        self.declustering_inds = inds
        self.decluster_centers = self.centers[inds]
        self.decluster_derivs = [d[inds] for d in self.adjusted_derivs]
        self.kd = sp.spatial.KDTree(self.decluster_centers)
        self.power = power
    def take_remainder_potential(self, centers, potential_derivatives, modes):
        freqs = modes.freqs
        diag = freqs ** 2
        harmonic_contribs = np.dot(centers**2, diag) / 2
        harmonic_grad = diag[np.newaxis, :] * centers
        harmonic_deriv = np.diag(diag)
        remainder_vals = potential_derivatives[0] - harmonic_contribs
        remainder_grad = potential_derivatives[1] - harmonic_grad
        remainder_hess = potential_derivatives[2] - harmonic_deriv[np.newaxis, :, :]
        return [remainder_vals, remainder_grad, remainder_hess] + potential_derivatives[3:]

    def take_ppf_remainder(self, centers, potential_derivatives):
        deriv_order = len(potential_derivatives)-1
        ppf_data = self.ppfs(centers, deriv_order=deriv_order)
        if deriv_order > 0:
            vals = [pd - ppd for pd, ppd in zip(potential_derivatives, ppf_data)]
        else:
            vals = potential_derivatives - ppf_data
        return vals
    def taylor_interp(self, points, dists, neighbors, derivs, power=None, deriv_order=None):
        if power is None:
            power = self.power
        disps = points[:, np.newaxis, :] - neighbors

        idw_weights = InverseDistanceWeightedInterpolator.get_idw_weights(points, dists,
                                                                          disps=disps,
                                                                          zero_tol=1e-6,
                                                                          deriv_order=None,
                                                                          power=power
                                                                          )

        # if deriv_order is not None:
        #     idw_derivs = idw_weights
        #     idw_weights = idw_weights[0]
        # else:
        #     idw_derivs = None

        vals = derivs[0]
        for n,d in enumerate(derivs[1:]):
            for _ in range(n+1):
                d = nput.vec_tensordot(d, disps, shared=2, axes=[2, 2]) # one axis vanishes every time
            vals += d / np.math.factorial(n+1)

        vals = nput.vec_tensordot(
            vals,
            idw_weights,
            shared=1,
            axes=[1, 1]
        )
        if deriv_order is not None:
            # og_derivs = derivs

            # use taylor series to get derivatives
            new_derivs = []
            for i in range(1, deriv_order+1):
                # taylor series
                derv = derivs[i]
                for n,d in enumerate(derivs[i+1:]):
                    for _ in range(n+1):
                        d = nput.vec_tensordot(d, disps, shared=2, axes=[2, 2]) # one axis vanishes every time
                    derv += d / np.math.factorial(n+1)
                # interpolate
                new_derivs.append(
                    nput.vec_tensordot(idw_weights, derv, shared=1, axes=[1, 1])
                )
            # base_derivs = [np.moveaxis(b, 1, -1) for b in derivs[:deriv_order+1]]
            # derivs = [nput.vec_tensordot(idw_weights, d, shared=1, axes=[1, -1]) for d in base_derivs[1:]]
            # if deriv_order > 1:
            #     # we can correct the first deriv from the Hessian
            #     hess_contrib = np.moveaxis(
            #         nput.vec_tensordot(og_derivs[2], disps, shared=2, axes=[2, 2]),
            #         1, -1
            #     )
            #     derivs[0] =
            #     # derivs[0] = derivs[0] - nput.vec_tensordot(idw_weights, hess_contrib, shared=1, axes=[1, -1])
            # if deriv_order > 2:
            #     raise NotImplementedError(...)
            vals = [vals] + new_derivs
            # vals = [vals] + nput.tensordot_deriv(idw_derivs, base_derivs, deriv_order, shared=1, axes=[1, 1])

        return vals

    def __call__(self, points, deriv_order=None, *, neighborhood_size=None, power=None, **kwargs):
        if neighborhood_size is None: neighborhood_size = self.opts['neighborhood_size']
        if power is None: power = self.power

        dists, inds = self.kd.query(points, k=neighborhood_size)

        interp_data = self.taylor_interp(
            points,
            dists,
            self.decluster_centers[inds],
            [d[inds] for d in self.decluster_derivs],
            power=power,
            deriv_order=deriv_order
        )

        if self.ppfs is not None:
            ppd = self.ppfs(points, deriv_order=deriv_order)
            if deriv_order is None or deriv_order == 0:
                interp_data += ppd
            else:
                interp_data = [id + pd for id,pd in zip(interp_data, ppd)]

        if self.include_harmonic_basis:
            freqs = self.modes.freqs
            diag = freqs ** 2
            harmonic_contribs = np.dot(points ** 2, diag) / 2
            if deriv_order is not None:
                interp_data[0] = interp_data[0] + harmonic_contribs
                if deriv_order > 0:
                    harmonic_grad = diag[np.newaxis, :] * points
                    interp_data[1] = interp_data[1] + harmonic_grad
                if deriv_order > 1:
                    harmonic_hess = np.diag(diag)[np.newaxis, :, :]
                    interp_data[2] = interp_data[2] + harmonic_hess
            else:
                interp_data = interp_data + harmonic_contribs

        if self.harmonic_distance_cutoff is not None:
            min_dists = np.min(dists, axis=1)
            bad_spots = np.where(min_dists > self.harmonic_distance_cutoff)
            if len(bad_spots) > 0 and len(bad_spots[0]) > 0:
                bad_spots = bad_spots[0]
                bad_points = points[bad_spots]

                diag = self.modes.freqs ** 2
                harmonic_contribs = np.dot(bad_points ** 2, diag) / 2

                if deriv_order is not None:
                    interp_data[0][bad_spots] = harmonic_contribs
                    if deriv_order > 0:
                        harmonic_grad = diag[np.newaxis, :] * bad_points
                        interp_data[1][bad_spots] = harmonic_grad
                    if deriv_order > 1:
                        harmonic_hess = np.diag(diag)[np.newaxis, :, :]
                        interp_data[2][bad_spots] = np.broadcast_to(
                            harmonic_hess,
                            (len(bad_spots),) + harmonic_hess.shape[1:]
                        )
                else:
                    interp_data[bad_spots] = harmonic_contribs

        return interp_data

class DGBWatsonLeastSquaresInterpolator(DGBInterpolator):
    ...

DGBWatsonInterpolator = DGBWatsonTaylorInterpolator

class DGBCartesianWatsonInterpolator(DGBWatsonInterpolator):
    def __init__(self, centers, potential_derivatives, modes, **opts):
        from .Coordinates import DGBWatsonModes

        super().__init__(
            DGBWatsonModes.embed_coords(centers, modes),
            DGBWatsonModes.embed_derivs(potential_derivatives, modes),
            modes,
            **opts
        )
    def __call__(self, points, deriv_order=None, **kwargs):
        from .Coordinates import DGBWatsonModes

        if deriv_order is not None and deriv_order > 0:
            raise NotImplementedError("currently just for potential plotting...")

        emb_coords = DGBWatsonModes.embed_coords(points, self.modes)
        return super().__call__(
            emb_coords,
            deriv_order=deriv_order
        )