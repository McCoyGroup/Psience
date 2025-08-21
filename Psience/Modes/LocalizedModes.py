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

    def __init__(self,
                 normal_modes: MixtureModes,
                 transformation,
                 inverse=None,
                 origin=None,
                 masses=None,
                 freqs=None,
                 mass_weighted=None,
                 frequency_scaled=None,
                 **etc
                 ):
        if len(transformation) == 2 and nput.is_array_like(transformation[0]):
            tf, inv = transformation
            tf = np.asanyarray(tf)
            if tf.ndim == 2:
                transformation = tf
                if inverse is None:
                    inverse = np.asanyarray(inv)
        else:
            transformation = np.asanyarray(transformation)
            if inverse is None:
                inverse = transformation.T
        mat = normal_modes.modes_by_coords @ transformation
        self.localizing_transformation = (transformation, inverse)
        self.base_modes = normal_modes
        inverse = inverse @ normal_modes.coords_by_modes
        if origin is None:
            origin = normal_modes.origin
        if masses is None:
            masses = normal_modes.masses
        super().__init__(
            normal_modes.basis,
            mat,
            inverse=inverse,
            origin=origin,
            masses=masses,
            freqs=normal_modes.freqs,
            mass_weighted=normal_modes.mass_weighted,
            frequency_scaled=frequency_scaled,
            **etc
        )
        self._freqs = freqs

    def __getitem__(self, item):
        """
        Takes a slice of the modes
        :param item:
        :type item:
        :return:
        :rtype:
        """

        if isinstance(item, (int, np.integer)):
            item = (item,)
        elif isinstance(item, slice):
            ...
        elif not isinstance(item[0], (int, np.integer)):
            item = tuple(item[0])

        # sub_modes = self.matrix[:, item]
        # inv = self._inv
        # if inv is not None:
        #     inv = inv[item, :]
        # freq = self.freqs[item,]
        tf, inv = self.localizing_transformation
        hmm = self.modify(
            transformation=(tf[:, item], inv[item, :])
        )

        return hmm

    @property
    def freqs(self):
        if self._freqs is None:
            self._freqs = self.compute_freqs()
        return self._freqs
    @freqs.setter
    def freqs(self, freqs):
        ...
        # if self._in_init:
        #     ...
        # else:
        #     self._freqs = freqs
    @property
    def mass_weighted(self):
        return self.base_modes.mass_weighted
    @mass_weighted.setter
    def mass_weighted(self, new):
        if (
                new and not self.base_modes.mass_weighted
                or not new and self.base_modes.mass_weighted
        ):
            raise ValueError("can't set `mass_weighted` directly")
    # @property
    # def frequency_scaled(self):
    #     return self.base_modes.frequency_scaled
    # @frequency_scaled.setter
    # def frequency_scaled(self, new):
    #     if (
    #             new and not self.base_modes.frequency_scaled
    #             or not new and self.base_modes.frequency_scaled
    #     ):
    #         raise ValueError("can't set `frequency_scaled` directly")
    @property
    def g_matrix(self):
        return self.base_modes.g_matrix
    #     if self.base_modes.g_matrix is not None:
    #         tf, inv = self.localizing_transformation
    #         base_gm = self.base_modes.g_matrix
    #         if base_gm is None:
    #             base_gm = self.base_modes.compute_gmatrix()
    #         pinv = self.modes_by_coords @ inv.T
    #         return self.coords_by_modes @ inv.T @ base_gm @ inv
    #     else:
    #         return None
    @g_matrix.setter
    def g_matrix(self, g):
        ...

    def modify(self,
               base_modes=None,
               *,
               transformation=None,
               freqs=None,
               origin=None,
               masses=None,
               inverse=None,
               name=None,
               mass_weighted=None,
               frequency_scaled=None,
               g_matrix=None
               ):
        if (
                transformation is not None
                and inverse is None
                and len(transformation) == 2
                and nput.is_array_like(transformation[0])
        ):
            tf, inv = transformation
            tf = np.asanyarray(tf)
            if tf.ndim == 2:
                transformation = tf
                inverse = inv
        return type(self)(
            self.base_modes if base_modes is None else base_modes,
            self.localizing_transformation[0] if transformation is None else transformation,
            inverse=self.localizing_transformation[1] if inverse is None else inverse,
            origin=origin,
            masses=masses,
            freqs=freqs,
            g_matrix=g_matrix,
            name=self.name if name is None else name,
            frequency_scaled=frequency_scaled
        )

    def _frequency_scaling(self, freqs=None):
        # L = self.matrix.shape.T
        if freqs is None:
            freqs = self.local_freqs
        conv = np.sqrt(freqs)
        return freqs, conv

    def make_mass_weighted(self, **kwargs):
        return self.modify(self.base_modes.make_mass_weighted(**kwargs))
    def remove_mass_weighting(self, **kwargs):
        return self.modify(self.base_modes.remove_mass_weighting(**kwargs))
    def make_frequency_scaled(self, freqs=None, **kwargs):
        if self.frequency_scaled: return self
        freqs, conv = self._frequency_scaling(freqs=freqs)
        tf, inv = self.localizing_transformation

        return self.modify(
            transformation=(tf@np.diag(conv), np.diag(1/conv)@inv),
            frequency_scaled=True
        )
    def remove_frequency_scaling(self, freqs=None, **kwargs):
        if not self.frequency_scaled: return self
        freqs, conv = self._frequency_scaling(freqs=freqs)
        tf, inv = self.localizing_transformation

        return self.modify(
            transformation=(tf@np.diag(1/conv), np.diag(conv)@inv),
            frequency_scaled=False
        )

    # @property
    # def local_hessian(self):
    #     tf, inv = self.localizing_transformation
    #     f = inv @ np.diag(self.freqs ** 2) @ inv.T
    #     g = self.g_matrix
    #     a = np.diag(np.power(np.diag(g) / np.diag(f), 1 / 4))
    #     return a @ f @ a

    def compute_hessian(self, system='modes'):
        if system == 'modes':
            tf, inv = self.localizing_transformation
            return inv @ np.diag(np.sign(self.base_modes.freqs) * (self.base_modes.freqs ** 2)) @ inv.T
        elif system == 'coords':
            return self.base_modes.compute_hessian('coords')
        else:
            raise ValueError(f'unknown system for normal modes "{system}", valid are "modes", "coords"')

    # def localize(self,
    #              method=None,
    #              **opts
    #              ):
    #     return NormalModes.localize(self, method=method, **opts)

    def apply_transformation(self, transformation, inverse=None, **opts):
        if len(transformation) == 2 and nput.is_array_like(transformation[0]):
            tf, inv = transformation
            tf = np.asanyarray(tf)
            if tf.ndim == 2:
                transformation = tf
                if inverse is None:
                    inverse = np.asanyarray(inv)
        else:
            transformation = np.asanyarray(transformation)
            if inverse is None:
                inverse = transformation.T
        tf_1, inv_1 = self.localizing_transformation
        transformation = tf_1 @ transformation
        inverse = inverse @ inv_1

        return type(self)(
            self.base_modes,
            transformation,
            inverse=inverse,
            **opts
        )

    def get_complement(self, concatenate=False):
        rem = nput.find_basis(
            nput.orthogonal_projection_matrix(self.localizing_transformation[0])
        )
        submodes = self.base_modes.apply_transformation(rem)
        if not concatenate:
            loc_freqs, loc_modes = scipy.linalg.eigh(submodes.compute_hessian(), submodes.compute_gmatrix(), type=3)
            return submodes.apply_transformation(loc_modes)
        else:
            g = submodes.compute_gmatrix()
            g12 = nput.fractional_power(g, 1/2)
            gi12 = nput.fractional_power(g, -1/2)
            f = submodes.compute_hessian()
            _, mw_modes = np.linalg.eigh(g12 @ f @ g12)
            modes, inv_loc = g12 @ mw_modes, mw_modes.T @ gi12
            # loc_freqs, loc_modes = scipy.linalg.eigh(submodes.compute_hessian(), submodes.compute_gmatrix(), type=3)
            tf_rem = rem @ modes
            tf, inv = self.localizing_transformation
            inv_rem = inv_loc @ rem.T

            tf_full = np.concatenate([tf, tf_rem], axis=1)
            inv_full = np.concatenate([inv, inv_rem], axis=0)

            return self.base_modes.apply_transformation((tf_full, inv_full))
        # return self.base_modes.localize(
        #     projections=[self.make_mass_weighted().modes_by_coords],
        #     orthogonal_projection=True
        # )