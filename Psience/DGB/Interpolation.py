
import abc

import numpy as np, scipy as sp

from McUtils.Zachary import InverseDistanceWeightedInterpolator
import McUtils.Numputils as nput

__all__ = [
    "DGBInterpolator",
    "DGBGenericInterpolator",
    "DGBWatsonInterpolator",
    "DGBCartesianWatsonInterpolator",
]

class DGBInterpolator(metaclass=abc.ABCMeta):
    default_neighborhood_size=15
    default_declustering_alpha=0
    default_declustering_overlap=.1
    def __init__(self, centers, potential_derivatives,
                 declustering_alpha=None,
                 declustering_overlap=None,
                 neighborhood_size=None, **opts):
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
        self.opts = dict(opts,
                         neighborhood_size=neighborhood_size
                         )
    @abc.abstractmethod
    def __call__(self, points, deriv_order=None, **kwargs) -> 'np.ndarray | list[np.ndarray]':
        ...
class DGBGenericInterpolator(DGBInterpolator):

    def __init__(self, centers, potential_derivatives, **opts):
        super().__init__(centers, potential_derivatives, **opts)
        self.interp = InverseDistanceWeightedInterpolator(
            self.centers,
            *self.derivs,
            **self.opts
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

class DGBWatsonTaylorInterpolator(DGBInterpolator):
    def __init__(self, centers, potential_derivatives, modes,
                 power=2,
                 include_harmonic_basis=False, # provably does nothing for now...
                 harmonic_distance_cutoff=None, # just to make plots look better...
                 **opts
                 ):
        from .Evaluators import DGBEvaluator
        from .Gaussians import DGBGaussians

        super().__init__(centers, potential_derivatives, power=power, **opts)
        self.modes = modes
        alpha = self.declustering_opts['declustering_alpha']
        if alpha > 0:
            overlap_data = DGBEvaluator.get_overlap_gaussians(
                    self.centers,
                    np.full(self.centers.shape, alpha),
                    None
                )
            S = DGBEvaluator.evaluate_overlap(overlap_data)
            inds = DGBGaussians._optimize_gs_block(S, self.declustering_opts['declustering_overlap'])
        else:
            inds = slice(None)
        if include_harmonic_basis:
            self.adjusted_derivs = self.take_remainder_potential(centers, potential_derivatives, modes)
        else:
            self.adjusted_derivs = self.derivs
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
    def taylor_interp(self, points, dists, neighbors, derivs, power=2, deriv_order=None):
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
            base_derivs = [np.moveaxis(b, 1, -1) for b in derivs[:deriv_order+1]]
            derivs = [nput.vec_tensordot(idw_weights, d, shared=1, axes=[1, -1]) for d in base_derivs[1:]]
            vals = [vals] + derivs
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
        return super().__call__(
            DGBWatsonModes.embed_coords(points, self.modes),
            deriv_order=deriv_order
        )