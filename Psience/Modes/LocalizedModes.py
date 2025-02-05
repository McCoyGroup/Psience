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
                 normal_modes: NormalModes,
                 transformation,
                 inverse=None,
                 origin=None,
                 masses=None,
                 freqs=None,
                 mass_weighted=None,
                 frequency_scaled=None,
                 **etc
                 ):
        mat = normal_modes.matrix @ transformation
        if inverse is None:
            inverse = transformation.T
        self.localizing_transformation = (transformation, inverse)
        self.base_modes = normal_modes
        inverse = inverse @ normal_modes.inverse
        if origin is None:
            origin = normal_modes.origin
        if masses is None:
            masses = normal_modes.masses
        if freqs is None:
            freqs = normal_modes.freqs  # this preserves information, but could be confusing...
        super().__init__(
            normal_modes.basis,
            mat,
            inverse=inverse,
            origin=origin,
            masses=masses,
            freqs=freqs,
            mass_weighted=normal_modes.mass_weighted,
            frequency_scaled=normal_modes.frequency_scaled,
            **etc
        )

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
    @property
    def frequency_scaled(self):
        return self.base_modes.frequency_scaled
    @frequency_scaled.setter
    def frequency_scaled(self, new):
        if (
                new and not self.base_modes.frequency_scaled
                or not new and self.base_modes.frequency_scaled
        ):
            raise ValueError("can't set `frequency_scaled` directly")
    @property
    def g_matrix(self):
        return self.matrix.T @ self.base_modes.g_matrix @ self.matrix
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
        return type(self)(
            self.base_modes if base_modes is None else base_modes,
            self.localizing_transformation[0] if transformation is None else transformation,
            inverse=self.localizing_transformation[1] if inverse is None else inverse,
            origin=origin,
            masses=masses,
            freqs=freqs,
            g_matrix=g_matrix,
            name=self.name if name is None else name
        )

    def make_mass_weighted(self, **kwargs):
        return self.modify(self.base_modes.make_mass_weighted(**kwargs))
    def remove_mass_weighting(self, **kwargs):
        return self.modify(self.base_modes.remove_mass_weighting(**kwargs))
    def make_frequency_scaled(self, **kwargs):
        return self.modify(self.base_modes.make_frequency_scaled(**kwargs))
    def remove_frequency_scaling(self, **kwargs):
        return self.modify(self.base_modes.remove_frequency_scaling(**kwargs))

    @property
    def local_freqs(self):
        return np.diag(self.local_hessian)

    @property
    def local_hessian(self):
        tf, inv = self.localizing_transformation
        f = inv @ np.diag(self.freqs ** 2) @ inv.T
        g = self.g_matrix
        a = np.diag(np.power(np.diag(g) / np.diag(f), 1/4))
        return a @ f @ a