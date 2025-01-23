"""
Provides support for handling modes that arise from
"""

import numpy as np
import McUtils.Numputils as nput
from McUtils.Coordinerds import CoordinateSystem, CartesianCoordinateSystem3D, InternalCoordinateSystem

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
                 origin=None,
                 masses=None,
                 inverse=None,
                 name=None,
                 ):
        if (
                isinstance(coeffs, np.ndarray) or
                (not nput.is_numeric(coeffs[0]) and nput.is_numeric(coeffs[0][0]))
        ):
            coeffs = [coeffs]
        coeffs = [np.asanyarray(c) for c in coeffs]
        full_coeffs = coeffs
        coeffs = coeffs[0]

        super().__init__(
            matrix=coeffs,
            inverse=inverse,
            name=self.name if name is None else name,
            basis=basis,
            dimension=(len(freqs),),
            origin=origin
        )
        self.freqs = freqs
        self.masses = masses
        self._extended_coeffs = full_coeffs
        self._inverse_coeffs = None

    @classmethod
    def prep_modes(cls, modes):
        if isinstance(modes, cls):
            return modes

        matrix = None
        inverse = None
        basis = None
        if isinstance(modes, dict):
            opts = modes.copy()
            for k in ["matrix", "modes"]:
                if k in opts:
                    matrix = opts[k]
                    del opts[k]
        else:
            matrix = np.asanyarray(modes)
            opts = {}

        if 'inverse' in opts:
            inverse = opts['inverse']
            del opts['inverse']
        if 'basis' in opts:
            basis = opts['basis']
            del opts['basis']

        if matrix is None and inverse is None:
            raise ValueError(f"can't prep {cls.__name__} without matrix or inverse")

        if basis is None:
            if matrix.shape[0] % 3 == 0:
                basis = CartesianCoordinateSystem3D
            else:
                basis = InternalCoordinateSystem(dimension=(None, matrix.shape[0]))

        cls(
            basis,
            matrix,
            inverse=inverse,
            **opts
        )

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
            masses=self.masses,
            inverse=inv
        )

    def rotate(self, rot, in_place=False):
        raise NotImplementedError("too confusing...")

    def transform(self, tf, inv=None, origin=None):
        raise NotImplementedError("ambiguous...")
        #TODO: handle Cartesian variant where tf gets broadcasted

        if inv is None:
            if tf.shape[0] != tf.shape[1]:
                raise ValueError("if given only a transformation, need corresponding inverse")
            else:
                inv = np.linalg.inv(tf)
        base_inv = self.inverse
        base_mat = self.matrix
        raise Exception(base_inv @ base_mat, base_mat.shape, base_inv.shape)
        new_mat = tf@base_mat
        new_inv = base_inv@inv

        if origin is None:
            origin = self.origin
            if origin is not None:
                origin = tf@origin.flatten()

        return type(self)(
            self.basis,
            new_mat,
            freqs=self.freqs,
            masses=self.masses,
            origin=origin,
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

    @property
    def total_transformation(self):
        return self._extended_coeffs
    @property
    def inverse_transformation(self):
        if self._inverse_coeffs is None:
            self._inverse_coeffs = nput.inverse_transformation(self.total_transformation,
                                                               len(self.total_transformation),
                                                               reverse_expansion=[self.inverse])
        return self._inverse_coeffs
    def embed_derivs(self, derivs):
        return nput.tensor_reexpand(self.total_transformation, derivs)
    def unembed_derivs(self, derivs):
        return nput.tensor_reexpand(self.inverse_transformation, derivs)