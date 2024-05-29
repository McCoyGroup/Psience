import numpy as np

from McUtils.Coordinerds import CartesianCoordinates1D, CartesianCoordinates2D, CartesianCoordinates3D
from McUtils.Scaffolding import Logger
from McUtils.Parallelizers import Parallelizer
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput

from ..Modes import NormalModes

from .Coordinates import *
from .Evaluators import *
from .Solvers import *
from .Interpolation import *

__all__ = [
    "DGBGaussians"
]

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

    def __init__(self, coords, alphas,
                 transformations=None, *,
                 momenta=None,
                 poly_coeffs=None,
                 kinetic_options=None,
                 logger=None,
                 parallelizer=None
                 ):
        self._pref = None
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

        if isinstance(momenta, (int, float, np.integer, np.floating)):
            momenta = [momenta] * self.coords.shape[0]
        if momenta is not None:
            momenta = np.asanyarray(momenta)
            if momenta.ndim == 1:
                momenta = np.broadcast_to(momenta[:, np.newaxis], self.alphas.shape)
        self.momenta = momenta

        self._overlap_data = None

        if poly_coeffs is not None:
            poly_coeffs = self.canonicalize_poly_coeffs(poly_coeffs, self.alphas)
        self._poly_coeffs = poly_coeffs

        if kinetic_options is None:
            kinetic_options = {}
        self.kinetic_options = kinetic_options

        self.parallelizer = Parallelizer.lookup(parallelizer)
        self.logger = Logger.lookup(logger)

    @property
    def overlap_data(self):
        if self._overlap_data is None:
            self._overlap_data = DGBEvaluator.get_overlap_gaussians(
                self.coords.centers,
                self.alphas,
                self.transformations,
                self.momenta,
                logger=self.logger,
                parallelizer=self.parallelizer
            )
        return self._overlap_data
    def get_S(self, return_prefactor=False):
        if self._poly_coeffs is not None:
            raise NotImplementedError("need to reintroduce polynomial support")
        with self.logger.block(tag="Evaluating S matrix"):
            return DGBEvaluator.evaluate_overlap(self.overlap_data, logger=self.logger, return_prefactor=return_prefactor)
    def get_T(self):
        if self._poly_coeffs is not None:
            raise NotImplementedError("need to reintroduce polynomial support")
        with self.logger.block(tag="Evaluating kinetic energy matrix"):
            return self.coords.kinetic_energy_evaluator.evaluate_ke(self.overlap_data, logger=self.logger, **self.kinetic_options)

    def optimize(self, optimizer_options, potential_function=None, logger=None, **opts):
        if logger is None:
            logger = self.logger
        if optimizer_options is True:
            optimizer_options = 'gram-schmidt'
        if isinstance(optimizer_options, (int, float, np.integer, np.floating)):
            optimizer_options = {'method':'gram-schmidt', 'norm_cutoff':optimizer_options}
        if isinstance(optimizer_options, str):
            optimizer_options = {'method':optimizer_options}

        # dispatch
        method = optimizer_options['method'].lower()
        opts = dict(optimizer_options, **opts)
        del opts['method']
        if method == 'svd':
            optimized_positions = self._optimize_svd(logger=logger, **opts)
        elif method == 'gram-schmidt':
            optimized_positions = self._optimize_gs(logger=logger, **opts)
        elif method == 'energy-cutoff':
            optimized_positions = self._prune_energy(potential_function, logger=logger,  **opts)
        elif method == 'decluster':
            optimized_positions = self._prune_dists(logger=logger, **opts)
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

        if self.momenta is not None:
            momenta = self.momenta[full_good_pos,]
        else:
            momenta = None

        new = type(self)(
            centers, alphas, transformations,
            poly_coeffs=_poly_coeffs, logger=self.logger,
            momenta=momenta,
            kinetic_options=self.kinetic_options
        )
        
        if self._overlap_data is not None:
            new._overlap_data = self.overlap_data.take_subselection(full_good_pos)
        if self._pref is not None:
            base = self._pref
            idx = np.ix_(full_good_pos, full_good_pos)
            if isinstance(base, tuple):
                new._pref = (base[0][idx], base[1][idx])
            else:
                new._pref = base[idx]
        if self._S is not None:
            base = self._S
            idx = np.ix_(full_good_pos, full_good_pos)
            if isinstance(base, tuple):
                new._S = (base[0][idx], base[1][idx])
            else:
                new._S = base[idx]
        if self._T is not None:
            base = self._T
            idx = np.ix_(full_good_pos, full_good_pos)
            if isinstance(base, tuple):
                new._T = (base[0][idx], base[1][idx])
            else:
                new._T = base[idx]

        return new

    @staticmethod
    def _calc_gs_norm(new, S, old):
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
    @staticmethod
    def _calc_gs_next(vec, ortho):
        new = np.zeros(len(vec))
        new[-1] = 1
        for prev in ortho:
            n = len(prev)
            new[:n] -= np.dot(vec[:n], prev) * prev
        return new

    _fast_outer = None
    @classmethod
    def _gs_off_diag_norms(cls, S, rows, cols, vec): # this is the slow step...
        if cls._fast_outer is None:
            import time
            # from McUtils.Misc import jit
            # @jit(nopython=True, cache=True, warn=True)
            def fast_mul(vec:np.ndarray, rows:np.ndarray, cols:np.ndarray):
                # wtf = time.time()
                mul = (vec[:, rows] * vec[:, cols])
                # elp = time.time() - wtf
                # print("--->", vec.shape, np.average(np.sum(np.abs(vec) > 1e-3, axis=1)), len(rows), f"{elp}s")
                return mul
            cls._fast_outer = fast_mul


        return 2 * np.matmul(
            S[rows, cols][np.newaxis, np.newaxis, :],
            cls._fast_outer(vec, rows, cols)[:, :, np.newaxis]
        ).reshape(-1)
    @staticmethod
    def _gs_diag_norms(vec):
        return np.matmul(
            vec[:, np.newaxis, :],
            vec[:, :, np.newaxis]
        ).reshape(-1)  # taking advantage of the fact that the diagonal of S is 1
    @staticmethod
    def _gs_extra_norm(old, vec):
        # correction from old itself, being a column of the overlap matrix, knowing that the final element of new is 1
        return 1 + 2 * np.matmul(
            old[:, np.newaxis, :], vec[:, :, np.newaxis]
        ).reshape(-1)
    # @staticmethod
    # def _gs_direct_norm(S, vec):
    #     print(S.shape, vec.shape)
    #     ugh = S[:len(vec), :len(vec)]
    #     Sv = np.matmul(ugh, vec).reshape(vec.shape)
    #     return np.matmul(Sv[:, np.newaxis, :], vec[:, :, np.newaxis]).reshape(-1)
    @classmethod
    def _gs_direct_norm(cls, S, vec):
        c1 = np.dot(vec, S)# I just assume dot is faster
        return np.matmul(
            c1[:, np.newaxis, :],
            vec[:, :, np.newaxis]
        ).flatten()
    @classmethod
    def _calc_gs_norm_block(cls, new, S, old):
        # generalizes the row-by-row form to calculate the norms for every vector at once so
        # we can figure out which ones to keep

        # This provides the contribution from (k-1)x(k-1)
        # of older terms
        vec = new[:, :-1]
        rows, cols = np.triu_indices(vec.shape[-1], k=1)

        # off_diag = cls._gs_off_diag_norms(S, rows, cols, vec)
        # diag = cls._gs_diag_norms(vec)
        # base_norm = diag + off_diag
        subs = np.eye(vec.shape[-1])
        subs[rows, cols] = S[rows, cols]
        subs[cols, rows] = S[rows, cols]
        base_norm = cls._gs_direct_norm(subs, vec)
        correction_norm = cls._gs_extra_norm(old, vec)

        return base_norm + correction_norm
    @staticmethod
    def _calc_gs_next_block(vecs, ortho):
        new = np.zeros(vecs.shape)
        new[:, -1] = 1
        for prev in ortho:
            n = len(ortho)
            new[:, :n] -= np.dot(vecs[:, :n], prev[:n])[:, np.newaxis] * prev[np.newaxis, :n]
        return new
    @staticmethod
    def _calc_gs_next_block_single(cur, overlaps, prev):
        cur = cur.copy()
        cur[:, -1] = 1
        n = len(prev)
        changes = np.dot(overlaps[:, :n], prev[:n])[:, np.newaxis] * prev[np.newaxis, :n]
        cur[:, :n] -= changes
        return changes, cur
    @classmethod
    def _optimize_gs_iter(cls, S, overlap_cutoff, norm_truncation_cutoff, logger=None):
        orthog_vecs = [np.array([1.])]
        norms = []  # debug
        mask = np.full(len(S), False, dtype=bool)
        mask[0] = True
        logger.log_print("pruning with unpivoted Gram-Schmidt norm > {n}", n=norm_truncation_cutoff)
        for i in range(1, len(S)):
            mask[i] = True
            vec = S[i, :i + 1][mask[:i + 1]]
            new = cls._calc_gs_next(vec, orthog_vecs)
            norm = cls._calc_gs_norm(new, S[np.ix_(mask[:i], mask[:i])], vec)
            norms.append(norm)
            if norm > overlap_cutoff:
                if norm < norm_truncation_cutoff:
                    norm = norm_truncation_cutoff
                orthog_vecs.append(new / np.sqrt(norm))
            else:
                mask[i] = False

        return np.where(mask)[0]

    @staticmethod
    def _pivot_vector(vec, i, pivot):
        new_pivot = vec[i + pivot]
        vec[i + 1:i + 1 + pivot] = vec[i:i + pivot]
        vec[i] = new_pivot
        return vec
    @classmethod
    def _pivot_matrix(cls, mat, i, pivot):
        pivot_row = cls._pivot_vector(mat[i + pivot, :].copy(), i, pivot)
        pivot_col = cls._pivot_vector(mat[:, i + pivot].copy(), i, pivot)
        pivot_pos = i + pivot

        # Unnecessarily complicated, maybe more efficient...
        mat[i + 1:pivot_pos + 1, :] = mat[i:pivot_pos, :]
        mat[:, i + 1:pivot_pos + 1] = mat[:, i:pivot_pos]
        mat[i, :] = pivot_row
        mat[:, i] = pivot_col
        mat[i, i+1:pivot_pos+1] = pivot_col[i+1:pivot_pos+1]
        mat[i+1:pivot_pos+1, i] = pivot_row[i+1:pivot_pos+1]

        return mat

    @classmethod
    def _optimize_gs_block(cls, S, overlap_cutoff, max_overlap_component,
                           # recalculate_cutoff_factor=0,
                           # chunk_size=5
                           logger=None
                           ):

        logger.log_print("pruning with pivoted Gram-Schmidt norm > {oc}", oc = overlap_cutoff)

        pivots = np.arange(len(S))
        # S_og = S
        S = S.copy()  # we'll use the upper triangle for norms and the lower triangle for storing vecs
        cols, rows = np.triu_indices_from(S, k=1)
        S[rows, cols] = 0 # remove lower triangle to redefine later
        if max_overlap_component < 1: # a fast pre-filter
            bad_pos = np.where(S[cols, rows] > max_overlap_component)
            if len(bad_pos[0]) > 0:
                # we contract S to reflect all of the things we've dropped
                kept = np.setdiff1d(np.arange(len(S)), rows[bad_pos[0]])
                pivots = pivots[kept]
                full_spec = kept
                S = S[np.ix_(full_spec, full_spec)]
                logger.log_print("dropped {n} Gaussians with overlaps greater than {mc}",
                                 n=len(np.unique(rows[bad_pos[0]])),
                                 mc=max_overlap_component)
        # all_norms = None
        for i in range(1, len(S)):

            all_vecs = S[i:, :i + 1]
            all_ovs = S[:i + 1, i:].T
            # all_news = cls._calc_gs_next_block(all_vecs, [S[:j+1] for j in range(i)])  # subtract off newest change
            changes, all_news = cls._calc_gs_next_block_single(all_vecs, all_ovs, S[i-1, :i])  # subtract off newest change
            all_norms = cls._calc_gs_norm_block(all_news, S, S[:i, i:].T)

            # rem = len(pivots) - i
            good_pos = np.where(all_norms > overlap_cutoff)
            # print("-->", len(good_pos[0]))
            if len(good_pos) == 0 or len(good_pos[0]) == 0:
                # we know everything else can be eliminated
                pivots = pivots[:i]
                break
            else:
                # we contract S to reflect all of the things we've dropped
                kept = good_pos[0]
                pivots = np.concatenate([pivots[:i], pivots[i + kept]])
                if len(kept) == 1:
                    break

                full_spec = np.concatenate([np.arange(i), i + kept])
                S = S[np.ix_(full_spec, full_spec)]
                # S_og = S_og[np.ix_(full_spec, full_spec)]
                all_news = all_news[kept,]
                all_norms = all_norms[kept,]

            S[i:, :i] = all_news[:, :i]  # update without normalizing for now
            pivot_pos = np.argmax(all_norms)  # find the position with the largest remaining norm
            new_vec = all_news[pivot_pos] / np.sqrt(all_norms[pivot_pos])  # normalize just this vector
            # print(all_norms[pivot_pos])

            # Rearrange the pivots now that we have the new ordering
            pivots = cls._pivot_vector(pivots, i, pivot_pos)
            # all_norms = cls._pivot_vector(all_norms, i, pivot_pos)

            # rearrange S to reflect this new ordering
            S = cls._pivot_matrix(S, i, pivot_pos)
            # S_og = cls._pivot_matrix(S_og, i, pivot_pos)

            # now insert the normalized vector
            S[i, :i+1] = new_vec
            # print(S[:i+1, :i+1])
            # print(S_og[:i+1, :i+1])
            # if i > 1:
            #     raise Exception(...)
        return np.sort(pivots)

    # def _S_generator(self,
    #                  chunk_start, chunk_size
    #                  ):
    #
    #     DGBEvaluator.evaluate_overlap(self.overlap_data, logger=self.logger)
    #     self.get_S()
    gs_optimization_overlap_cutoff=1e-3
    def _optimize_gs(self, *, S=None,
                     norm_cutoff=None,
                     norm_truncation_cutoff=0,
                     max_overlap_cutoff=1,
                     allow_pivoting=True,
                     chunk_size=None,
                     logger=None
                     ):

        with logger.block(tag='optimizing using Gram-Schmidt'):
            if norm_cutoff is None:
                norm_cutoff = self.gs_optimization_overlap_cutoff
            if S is None:
                if chunk_size is None:
                    S=self.S
                    if isinstance(S, tuple): # phases introduced
                        S = 1/2*(S[0] + S[1])
                else:
                    raise NotImplementedError('getting chunking working is hard...')
                    S=self.S

            if allow_pivoting:
                pivots = self._optimize_gs_block(S, norm_cutoff, max_overlap_cutoff, logger=logger)
            else:
                pivots = self._optimize_gs_iter(S, norm_cutoff, norm_truncation_cutoff, logger=logger)
            logger.log_print("pruned {nD} Gaussians",  nD=len(S)-len(pivots))

        return pivots

    def _optimize_svd(self, *, S=None,
                      min_value=1e-12,
                      num_vectors=None,
                      contrib_cutoff=1e-3,
                      logger=None
                      ):
        if S is None:
            S = self.S

        with logger.block(tag='optimizing by taking subspaces of the overlap matrix'):
            sig, evecs = np.linalg.eigh(S)
            # for i in range(5):
            #     print("wtf", i, len(np.unique(np.where(np.abs(evecs[:, -(i+1):]) > contrib_cutoff)[0])), sig[-(i+1)],
            #           np.min(evecs[:, -(i+1)]), np.max(evecs[:, -(i+1)]), np.std(evecs[:, -(i+1)]))
            if num_vectors:
                good_loc = slice(max(len(sig) - num_vectors, 0), len(sig))
            else:
                # self.logger.log_print("most important center threshold: {t}", t=min_singular_value)
                good_loc = np.where(sig > min_value)[0]
            full_good_pos = np.unique(np.where(np.abs(evecs[:, good_loc]) > contrib_cutoff)[0])
            logger.log_print("pruned {nD} Gaussians",  nD=len(S)-len(full_good_pos))

        return full_good_pos

    default_energy_cutoff = 1600 / UnitsData.hartrees_to_wavenumbers
    def _prune_energy(self, pot, *, cutoff=None, probabilities=None, logger=None, potential_values=None):
        if potential_values is None:
            potential_values = pot(self.coords.centers)
        if probabilities is not None:
            with logger.block(tag="pruning Gaussians with energies by distribution {d} cm-1",
                              d=[
                                  [np.round(c * UnitsData.hartrees_to_wavenumbers, 3), 1 if p is None else p]
                                  for c, p in probabilities
                              ]):
                pots = np.sort(potential_values)
                cur_pivot = 0
                pivots = []
                for next_e, next_prob in probabilities:
                    pivot_pos = np.searchsorted(pots[cur_pivot:], next_e)+1 # position past where the energy would occur
                    if next_prob is None:
                        new_pivs = np.arange(cur_pivot, pivot_pos)
                    else:
                        if next_prob < 1:
                            probs = np.random.uniform(0, 1, size=pivot_pos)
                            sel = np.where(probs < next_prob)
                        else: # take this many points from that region
                            if pivot_pos >= next_prob:
                                sel = [
                                    np.sort(np.random.choice(np.arange(pivot_pos-1), next_prob, replace=False))
                                ]
                            else:
                                sel = [np.arange(pivot_pos-1)]
                        if len(sel) > 0 and len(sel[0]) > 0:
                            new_pivs = cur_pivot + sel[0]
                    cur_pivot = pivot_pos
                    pivots.append(new_pivs)
                full_good_pos = np.concatenate(pivots)
                logger.log_print("pruned {nD} Gaussians", nD=len(potential_values) - len(full_good_pos))
        else:
            if cutoff is not None:
                cutoff = self.default_energy_cutoff
            with logger.block(tag="pruning Gaussians with energies larger than {ec} cm-1",
                              ec=np.round(cutoff * UnitsData.hartrees_to_wavenumbers, 3)
                              ):
                full_good_pos = np.where(potential_values < cutoff)[0]
                logger.log_print("pruned {nD} Gaussians", nD=len(potential_values) - len(full_good_pos))
        return full_good_pos

    def _prune_dists(self, *,  cluster_radius, logger=None):
        with logger.block(tag=logger.log_print('declustering data with a radius of {r}', r=cluster_radius)):
            points = self.coords.centers
            pivots = np.arange(len(points))
            dec_pts = points.reshape(points.shape[0], -1)
            for i in range(len(points)):
                cur_pos = pivots[i]
                dists = np.linalg.norm(
                    dec_pts[cur_pos][np.newaxis, :] - dec_pts[pivots[i + 1:], :],
                    axis=1
                )
                good_pos = np.where(dists > cluster_radius)
                if len(good_pos) == 0 or len(good_pos[0]) == 0:
                    break
                pivots = np.concatenate([pivots[:i + 1], pivots[i + 1:][good_pos]])
            logger.log_print("pruned {nD} Gaussians", nD=len(self.coords.centers) - len(pivots))
        return pivots

    @classmethod
    def construct(cls,
                  coords,
                  alphas,
                  *,
                  potential_expansion=None,
                  potential_function=None,
                  transformations=None,
                  masses=None,
                  atoms=None,
                  modes=None,
                  kinetic_options=None,
                  internals=None,
                  coordinate_selection=None,
                  cartesians=None,
                  gmat_function=None,
                  momenta=None,
                  poly_coeffs=None,
                  logger=None,
                  pairwise_potential_functions=None,
                  parallelizer=None
                  ):
        if potential_expansion is not None:
            if potential_function is not None:
                raise ValueError("can't have both `potential_function` and `potential_expansion`")
            if isinstance(potential_expansion, dict):
                expansion_opts = potential_expansion.copy()
                del expansion_opts['values']
                potential_expansion = potential_expansion['values']
                interp_coords = expansion_opts.get('centers', coords)
                if 'centers' in expansion_opts:
                    del expansion_opts['centers']
            else:
                interp_coords = coords
                expansion_opts = {}
            potential_function = DGBGenericInterpolator(
                interp_coords,
                potential_expansion,
                pairwise_potential_functions=pairwise_potential_functions,
                **expansion_opts
            )

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

                    if momenta is not None: # defined as d/dx, even though classically given by m dx/dt
                        momenta = [np.zeros(momenta.shape[0]), momenta.reshape((momenta.shape[0], -1))]
                        momenta = DGBWatsonModes.embed_derivs(momenta, modes)[1]

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

                gmat_function = coords.gmatrix

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
            if gmat_function is None:
                gmat_function = coords.gmatrix
            if isinstance(transformations, dict):
                tf_opts = dict({'gmat_function':gmat_function}, **transformations)
                transformations = tf_opts['method']
                del tf_opts['method']
            else:
                tf_opts = {'gmat_function':gmat_function}
            if isinstance(transformations, str):
                if transformations == 'rpath':
                    transformations = cls.get_reaction_path_transformations(
                        coords.coords,
                        potential_function,
                        # masses=masses,
                        # atoms=atoms,
                        # internals=internals,
                        **tf_opts#gmat_function=gmat_function
                    )
                elif transformations == 'diag':
                    transformations = cls.get_hessian_diagonalizing_transformations(
                        coords.coords,
                        potential_function,
                        # masses=masses,
                        # atoms=atoms,
                        # internals=internals,
                        **tf_opts
                    )
                else:
                    raise ValueError("unknown transformation spec {}".format(modes))
            transformations = cls.canonicalize_transforms(coords.centers, transformations)

            if momenta is not None:
                momenta = DGBEvaluator.get_momentum_vectors(momenta, transformations)

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
                   poly_coeffs=poly_coeffs, logger=logger,
                   momenta=momenta,
                   kinetic_options=kinetic_options
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
                         stationary_point_norm=1e-2
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
            origin=reference_structure,
            dimensionless=False
        )

    @classmethod
    def get_reaction_path_transformations(
            cls,
            coords,
            potential_function,
            gmat_function,
            stationary_point_norm=1e-4,
            sort_alphas=True
    ):
        f_data = potential_function(coords, deriv_order=2)
        if isinstance(f_data, np.ndarray) and f_data.ndim == 3:
            hess = f_data
            grad = potential_function(coords, deriv_order=1)
            pot = potential_function(coords, deriv_order=1)
        else:
            pot, grad, hess = f_data

        hess_og = hess
        gmats = gmat_function(coords)
        # g12 = sp.linalg.fractional_matrix_power(gmats, 1 / 2)
        # gi12 = sp.linalg.fractional_matrix_power(gmats, -1 / 2)

        gvals, gvecs = np.linalg.eigh(gmats)
        if np.any((gvals <= 0).flatten()):
            raise ValueError("bad G-matrix?")
        g12_diags = np.zeros(gvecs.shape)
        diag_ings = (slice(None),) + np.diag_indices_from(gvecs[0])
        g12_diags[diag_ings] = np.sqrt(gvals)
        g12 = gvecs @ g12_diags @ gvecs.transpose(0, 2, 1)

        grad = np.reshape(grad[:, np.newaxis, :] @ g12, grad.shape) # mass-weight
        hess = g12 @ hess @ g12

        grad_norms = np.linalg.norm(grad, axis=-1)
        non_stationary = np.where(grad_norms >= stationary_point_norm)

        # in this case I'm _only_ projecting out the gradient for whatever that's worth...
        # which...probably has some rotation/translation component which
        # is why I don't get clean zero eigenvectors...
        rp_mode = -grad[non_stationary] / grad_norms[non_stationary][:, np.newaxis]
        proj = np.eye(rp_mode.shape[-1])[np.newaxis] - nput.vec_outer(rp_mode, rp_mode, axes=[1, 1])
        hess[non_stationary] = proj@hess[non_stationary]@proj

        freqs, g12_tfs = np.linalg.eigh(hess) # replace zero tf with mass-weighted rp_mode
        # raise Exception(np.sqrt(np.abs(freqs[7:15])) * 219475)

        # g12_diags[diag_ings] = 1 / np.sqrt(gvals)
        # g12_inv = gvecs @ g12_diags @ gvecs.transpose(0, 2, 1)

        tfs = g12 @ g12_tfs
        tfs[non_stationary, :, 0] = rp_mode

        # now we calculate the unweighted hessian diagonals to sort by
        if sort_alphas:
            tf_hess = tfs.transpose(0, 2, 1) @ hess_og @ tfs
            diags = np.argsort(np.abs(np.diagonal(tf_hess, axis1=1, axis2=2)), axis=1)

            tfs = tfs[
                np.arange(tfs.shape[0])[:, np.newaxis, np.newaxis],
                np.arange(tfs.shape[1])[np.newaxis, :, np.newaxis],
                diags[:, np.newaxis, :]
            ]

        inv = np.linalg.inv(tfs) #g12_tfs.transpose(0, 2, 1) @ g12_inv

        # return inv.transpose(0, 2, 1), tfs.transpose(0, 2, 1)
        return tfs, inv

    @classmethod
    def get_hessian_diagonalizing_transformations(
            cls,
            coords,
            potential_function,
            gmat_function
    ):
        f_data = potential_function(coords, deriv_order=2)
        if isinstance(f_data, np.ndarray) and f_data.ndim == 3:
            hess = f_data
            grad = potential_function(coords, deriv_order=1)
            pot = potential_function(coords, deriv_order=1)
        else:
            pot, grad, hess = f_data

        hess_og = hess
        gmats = gmat_function(coords)
        # g12 = sp.linalg.fractional_matrix_power(gmats, 1 / 2)
        # gi12 = sp.linalg.fractional_matrix_power(gmats, -1 / 2)

        gvals, gvecs = np.linalg.eigh(gmats)
        if np.any((gvals <= 0).flatten()):
            raise ValueError("bad G-matrix?")
        g12_diags = np.zeros(gvecs.shape)
        diag_ings = (slice(None),) + np.diag_indices_from(gvecs[0])
        g12_diags[diag_ings] = np.sqrt(gvals)
        g12 = gvecs @ g12_diags @ gvecs.transpose(0, 2, 1)

        hess = g12 @ hess @ g12

        _, g12_tfs = np.linalg.eigh(hess)

        tfs = g12 @ g12_tfs
        g12_diags[diag_ings] = 1 / np.sqrt(gvals)
        g12_inv = gvecs @ g12_diags @ gvecs.transpose(0, 2, 1)
        inv = g12_tfs.transpose(0, 2, 1) @ g12_inv

        # return tfs, inv
        return inv.transpose(0, 2, 1), tfs.transpose(0, 2, 1)

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
        if transformations is not None:
            inv, transformations = transformations
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

        alphas = scaling * np.sqrt(np.abs(np.diagonal(hess, axis1=1, axis2=2) / np.diagonal(gmats, axis1=1, axis2=2)))
        # print(hess)
        # print(
        #     np.sqrt(np.abs(np.diagonal(hess, axis1=1, axis2=2) / np.diagonal(gmats, axis1=1, axis2=2))) * 219475
        # )
        # print(
        #     np.sum(np.diagonal(hess, axis1=1, axis2=2) * 219475, axis=1)
        # )
        #[ 9426.1752331  10317.82470038  9108.73436512  9714.64698828  8579.37664985  9764.6991486   9764.69922505  9197.73091615  9197.73093045 11040.69284146  8823.36978916  9911.15518915  7911.13613538 10086.93927288 10086.93926775  9071.32381428  9071.32391714 11374.03054902  8618.17129632  9967.37235458  7484.76311907 10195.248758   10195.24883265  8970.84073996  8970.84093599 11209.80670989  8524.28409883  9866.77241186  7325.03365781  9995.46064476  9995.46077198  8839.0596967   8839.05989418 10599.3940994   8556.00072571  9632.81612374  7438.73297697  9579.709646    9579.709779    8683.59676335  8683.59680886]
        #[ 9426.1752331  10317.97352922  9342.83704949 10024.66092923  8579.25827834  9748.71327209  9748.71335084  9210.28053516  9210.28054868 11040.48757709  9044.73381151 10280.23127104  7910.82531061 10003.62991937 10003.62993513  9078.13129445  9078.1314092  11373.51783935  8846.43367466 10397.23539095  7484.34647915 10082.52862084 10082.52870943  8972.7922302   8972.79246371 11210.08657008  8784.25964707 10318.11304479  7324.5703684   9957.98059176  9957.98072257  8823.24924889  8823.24949849 10602.16863697  8879.27410661  9971.67663237  7438.26430876  9726.22755942  9726.22773469  8638.84650038  8638.84655498]

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
    def canonicalize_transforms(self, coords, tfs):
        npts = coords.shape[0]
        ndim = coords.shape[1]
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
    def prefactor(self):
        if self._pref is None:
            self._pref, self._S = self.get_S(return_prefactor=True)
        return self._pref
    @property
    def S(self):
        if self._S is None:
            self._pref, self._S = self.get_S(return_prefactor=True)
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

    bad_alpha_limit = 1e-15
    bad_scaling_limit = 1e-3
    def marginalize_out(self, indices, *, bad_alpha_limit=None, bad_scaling_limit=None):
        if bad_scaling_limit is None:
            bad_scaling_limit = self.bad_scaling_limit
        if bad_alpha_limit is None:
            bad_alpha_limit = self.bad_alpha_limit

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
            zp = np.abs(dropdets) < bad_scaling_limit #TODO: fix this in case someone is in the numerically unstable limit...
            dropdets[zp] = 1
            scaling = 1 / dropdets
            scaling[zp] = bad_scaling_limit
            scaling *= np.power(np.pi, len(indices) / 4)

            subcovs = full_covs[np.ix_(full_sel, remaining, remaining)]
            inv_alphas, tfs = np.linalg.eigh(subcovs)
            inv_alphas[inv_alphas < bad_alpha_limit] = bad_alpha_limit

            alphas = 1 / (2*inv_alphas) # we have 1/2a
            tfs = (tfs, tfs.transpose(0, 2, 1))

            momenta = self.momenta
            if momenta is not None:
                # we use the fact that the momenta are along the same axes as the old alphas to get
                scaled_moms = self.alphas * self.momenta
                # now we have to deal with a S_A @ T @ p type term by noting that we can turn T into its
                # [indices x len(alphas)] subblock
                cov_mom = self.transformations[0][:, remaining, :] @ scaled_moms[:, :, np.newaxis]
                inv_cov = DGBEvaluator.get_inverse_covariances(alphas, tfs)
                momenta = np.reshape(inv_cov @ cov_mom, alphas.shape)

        else:
            scaling = np.power(2 * np.pi, len(indices) / 4) / np.power(np.prod(self.alphas[:, indices], axis=-1), 1 / 4)
            alphas = self.alphas[:, remaining]
            tfs = None
            momenta = self.momenta
            if momenta is not None:
                momenta = momenta[:, remaining]

        return scaling, type(self)(
            subcoords,
            alphas=alphas,
            transformations=tfs,
            momenta = momenta
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
        momenta = self.momenta
        return type(self)(
            carts,
            alphas,
            transformations=transformations,
            momenta=momenta
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