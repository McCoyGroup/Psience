"""
Provides support for handling modes that arise from
"""

import numpy as np
from McUtils.Coordinerds import CoordinateSystem

# from .MoleculeInterface import AbstractMolecule
# from .Transformations import MolecularTransformation

__all__ = [
    "MixtureModes"
]
class MixtureModes(CoordinateSystem):
    """
    A `McUtils.Coordinerds.CoordinateSystem` object that expresses coordinates as
    a rotation on some base set of coordinates with some associated frequencies.
    """
    name="MixtureModes"
    def __init__(self,
                 basis,
                 coeffs,
                 freqs=None,
                 origin=None, inverse=None,
                 name=None
                 ):
        coeffs = np.asanyarray(coeffs)

        super().__init__(
            matrix=coeffs,
            inverse=inverse,
            name=self.name if name is None else name,
            basis=basis,
            dimension=(len(freqs),),
            origin=origin
        )
        self.freqs = freqs

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
        elif not isinstance(item[0], (int, np.integer)):
            item = tuple(item[0])

        sub_modes = self.matrix[:, item]
        inv = self._inv
        if inv is not None:
            inv = inv[item, :]
        freq = self.freqs[item,]
        return type(self)(
            self.basis,
            sub_modes,
            name=self.name,
            freqs=freq,
            origin=self._origin,
            inverse=inv
        )

    def rotate(self, rot, in_place=False):
        raise NotImplementedError("too confusing...")

    def transform(self, tf, inv=None):
        #TODO: handle Cartesian variant where tf gets broadcasted

        if inv is None:
            if tf.shape[0] != tf.shape[1]:
                raise ValueError("if given only a transformation, need corresponding inverse")
            else:
                inv = np.linalg.inv(tf)
        base_inv = self.inverse
        base_mat = self.matrix
        new_mat = tf@base_mat@inv
        new_inv = inv@base_inv@tf

        orig = self.origin
        if orig is not None:
            orig = tf@orig

        return type(self)(
            self.basis,
            new_mat,
            freqs=self.freqs,
            origin=orig,
            inverse=new_inv
        )

    @property
    def cartesian_modes(self):
        return self.origin.ndim == 2

    def embed_coords(self, carts):
        flat_carts = (carts - self.origin[np.newaxis]).reshape((len(carts), -1))
        return (flat_carts[:, np.newaxis, :] @ self.inverse.T[np.newaxis]).reshape(
            flat_carts.shape[0],
            self.matrix.shape[1]
        )
    def unembed_coords(self, mode_coords):
        origin = self.origin
        carts = (mode_coords[:, np.newaxis, :] @ self.matrix.T[np.newaxis]).reshape(
            mode_coords.shape[:1] + origin.shape
        )
        carts = carts + origin[np.newaxis]
        return carts
    def embed_derivs(self, derivs):
        fshape = derivs[0].ndim
        _ = []
        for n, d in enumerate(derivs):
            for j in range(n):
                d = np.tensordot(d, self.matrix, axes=[fshape, 0])
            _.append(d)
        return _
    def unembed_derivs(self, derivs):
        fshape = derivs[0].ndim
        _ = []
        for n, d in enumerate(derivs):
            for j in range(n):
                d = np.tensordot(d, self.inverse, axes=[fshape, 0])
            _.append(d)
        return _