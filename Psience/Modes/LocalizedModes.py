"""
Provides a general class of localized modes using a potentially non-linear
transformation of normal modes
"""

import numpy as np, scipy
import McUtils.Numputils as nput

from .MixtureModes import *
from .NormalModes import *

__all__ = ['LocalizedModes']
class LocalizedModes(MixtureModes):

    @staticmethod
    def _mass_mat(mass_spec):
        return np.diag(np.repeat(mass_spec[:, np.newaxis], 3, axis=1).flatten())
    @staticmethod
    def _frac_mat_power(mat, power):
        vals, vecs = np.linalg.eigh(mat)
        new_vals = vals**power
        return vecs @ np.diag(new_vals) @ vecs.T
    @staticmethod
    def _order_matrix(mat):
        # finds the permutation that makes the diagonal of the matrix the most ordered
        abs_mat = np.abs(mat)
        sort_mat = np.argsort(abs_mat, axis=1)
        # check if we have a clean permutation
        max_inds = sort_mat[:, -1]
        if np.all(np.sort(max_inds) == np.arange(len(max_inds))):
            perm = np.argsort(max_inds)
        else:
            # kinda ill-defined, but in the case of collisions we choose
            # the sorting to prefer the maximum value
            sort_maxima = np.abs(nput.vector_take(abs_mat, sort_mat[:, -1:], shared=1)).flatten()
            perm = np.lexsort([
                -sort_maxima,
                max_inds
            ])
        return perm, mat[perm,]
        # perm_inds = []
        # sort_maxima = np.flip(np.argsort(nput.vector_take(abs_mat, sort_mat[:, -1:], shared=1)))

    @classmethod
    def localize_by_masses(cls, modes:NormalModes, positions,
                           scaling=None,
                           min_scaling=1, max_scaling=10, steps=20,
                           similarity_cutoff=.90):
        base_tf = modes.make_mass_weighted().matrix
        # base_inv = modes.make_mass_weighted().inverse

        mw_F = base_tf @ np.diag(modes.freqs**2) @ base_tf.T
        g12i = cls._mass_mat(np.sqrt(modes.masses))
        base_F = g12i @ mw_F @ g12i
        # base_G = cls._mass_mat(1/modes.masses)
        if scaling is not None:

            scaled_masses = modes.masses.copy()
            scaled_masses[positions] *= scaling
            # local_G = cls._mass_mat(1/scaled_masses)
            new_modes = NormalModes.get_normal_modes(base_F, scaled_masses)
        else:
            # we find the minimum scaling factor that gives high similarity
            # with the heavy-atom modes

            localized_modes = []
            scaling_factors = np.linspace(min_scaling, max_scaling, steps)
            for s in scaling_factors:
                scaled_masses = modes.masses.copy()
                scaled_masses[positions] *= s
                # local_G = cls._mass_mat(1/scaled_masses)
                localized_modes.append(NormalModes.get_normal_modes(base_F, scaled_masses))

            # g12 = cls._frac_mat_power(modes.matrix.T @ modes.matrix, -1/2)
            # base_mode = modes.matrix.T
            local_g12s = [
                cls._frac_mat_power(l.modes.T @ l.modes, -1/2)
                for l in localized_modes
            ]
            scaled_modes = [l.modes @ lg for l, lg in zip(localized_modes, local_g12s)]

            row, col = np.triu_indices(len(localized_modes), k=1)
            correlation_mats = [
                # [
                cls._order_matrix(m2.T @ m)
                # ]
                for i, m in enumerate(scaled_modes)
                for m2 in scaled_modes[i+1:]
            ] # flat upper triangle

            correlations = np.array([
                np.count_nonzero(np.abs(np.diag(m)) > similarity_cutoff)
                for p,m in correlation_mats
            ])


            good_pos = np.where(correlations > modes.matrix.shape[1] - .1)
            left_inds = row[good_pos[0]]
            right_inds = col[good_pos[0]]

            pos, counts = np.unique(np.concatenate([left_inds, right_inds]), return_counts=True)
            max_count_pos = np.argsort(counts)


            new_modes = localized_modes[max_count_pos[-1]]

        return new_modes

