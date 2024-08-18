"""
Handles coordinate manipulation in a unified manner
"""
import numpy as np

from McUtils.Coordinerds import CoordinateSystem, CoordinateSystemConverter, CoordinateSystemConverters, ZMatrixCoordinates
import McUtils.Numputils as nput

class TotalCoordinateConverter:
    """
    A converter object that can apply normal mode transformations, non-linear transformations,
    and generic conversions (and does them in that order)
    """
    def __init__(self,
                 converter, inverse, reference_coords,
                 coordinate_expansion, inverse_expansion,
                 modes_matrix, modes_inverse
                 ):
        self.converter = converter
        self.inverse = inverse
        if converter is not None and inverse is None:
            raise ValueError("conversion function requires inverse")
        if inverse is not None and converter is None:
            raise ValueError("inverse conversion function requires forward function")
        self.ref = reference_coords
        self.coordinate_expansion, self.inverse_expansion = self.prep_expansions(coordinate_expansion, inverse_expansion)
        self.modes_matrix, self.modes_inverse = self.prep_modes(modes_matrix, modes_inverse)

    @classmethod
    def prep_modes(cls, modes_matrix, modes_inverse):
        if modes_matrix is not None and modes_inverse is None:
            modes_inverse = np.linalg.inv(modes_matrix)
        if modes_inverse is not None and modes_matrix is None:
            modes_matrix = np.linalg.inv(modes_inverse)

        return modes_matrix, modes_inverse

    @classmethod
    def prep_expansions(cls, coordinate_expansion, inverse_expansion):
        if coordinate_expansion is not None and inverse_expansion is None:
            inverse_expansion = nput.inverse_transformation(coordinate_expansion, len(coordinate_expansion))
        if inverse_expansion is not None and coordinate_expansion is None:
            coordinate_expansion = nput.inverse_transformation(inverse_expansion, len(inverse_expansion))

        return coordinate_expansion, inverse_expansion

    def convert(self, internals):
        internals = np.asanyarray(internals)
        smol = internals.ndim == 2
        if smol:
            internals = internals[np.newaxis]

        excess_shape = internals.shape[:-2]
        nc = internals.shape[-2]
        internals = internals.reshape(-1, internals.shape[-2]*internals.shape[-1])
        embedding_coords = [
            x for x in [0, 1, 2, 4, 5, 8]
            if x < internals.shape[-1]
        ]
        rem_coords = np.setdiff1d(np.arange(internals.shape[-1]), embedding_coords)
        embedding = internals[:, embedding_coords]
        internals = internals[:, rem_coords]
        if self.ref is not None:
            internals = internals - self.ref[np.newaxis]
        if self.modes_matrix is not None:
            internals = internals @ self.modes_matrix.T
        if self.coordinate_expansion is not None:
            coords = self.coordinate_expansion[0]
            for expansion_tensor in self.coordinate_expansion[1:]:
                contrib = expansion_tensor[np.newaxis]
                for _ in expansion_tensor.ndim:
                    contrib = nput.vec_tensordot(internals, expansion_tensor,
                                                 shared=1,
                                                 axes=[1, 1]
                                                 )
                coords += contrib
            internals = coords
        if self.converter is not None:
            internals = self.converter(internals)

        totals = np.empty((internals.shape[0], nc, 3), dtype=internals.dtype)
        totals[:, 0] = embedding[:, :3]
        if nc > 1:
            totals[:, 1, :1] = internals[:, 0]
            totals[:, 1, 1:] = embedding[:, 3:5]
        if nc > 2:
            totals[:, 2, :2] = internals[:, 1:3]
            totals[:, 2, 2] = embedding[:, 5]
        if nc > 3:
            totals[:, 3] = internals[:, 3:].reshape((internals.shape[0], nc-2, 3))

        if smol:
            totals = totals[0]

        return totals

    def invert(self, totals):
        totals = np.asanyarray(totals)
        smol = totals.ndim == 1
        if smol:
            totals = totals[np.newaxis]

        excess_shape = totals.shape[:-1]
        nc = totals.shape[-2]
        totals = totals.reshape(-1, totals.shape[-2] * totals.shape[-1])
        embedding_coords = [
            x for x in [0, 1, 2, 4, 5, 8]
            if x < totals.shape[-1]
        ]
        rem_coords = np.setdiff1d(np.arange(totals.shape[-1]), embedding_coords)
        embedding = totals[:, embedding_coords]
        totals = totals[:, rem_coords]

        if self.inverse is not None:
            totals = self.inverse(totals)

        if self.inverse_expansion is not None:
            coords = self.inverse_expansion[0]
            for expansion_tensor in self.inverse_expansion:
                contrib = expansion_tensor[np.newaxis]
                for _ in expansion_tensor.ndim:
                    contrib = nput.vec_tensordot(totals, expansion_tensor,
                                                 shared=1,
                                                 axes=[1, 1]
                                                 )
                coords += contrib
            totals = coords

        if self.modes_inverse is not None:
            totals = totals @ self.modes_inverse.T

        if self.ref is not None:
            totals = totals + self.ref[np.newaxis]

        internals = np.empty((totals.shape[0], nc, 3), dtype=totals.dtype)
        internals[:, 0] = embedding[:, :3]
        if nc > 1:
            internals[:, 1, :1] = totals[:, 0]
            internals[:, 1, 1:] = embedding[:, 3:5]
        if nc > 2:
            internals[:, 2, :2] = totals[:, 1:3]
            internals[:, 2, 2] = embedding[:, 5]
        if nc > 3:
            internals[:, 3] = totals[:, 3:].reshape((totals.shape[0], nc - 2, 3))

        internals = internals.reshape(excess_shape + (internals.shape[-1],))
        if smol:
            internals = internals[0]

        return internals

class TotalCoordinateSystem(CoordinateSystem):
    """
    A generic coordinate layer to unify coordinate handling,
    things written in terms of the basic `Coordinerds.CoordinateSet` should be
    converted over to this format bit-by-bit
    """

    # Design-wise, this is intended to solve the very common case
    # of having a set of coordinates that is defined as a function applied to a
    # set of base internal coordinates (or some linear combination of those base internal coordinates)
    # while still maintaining a reference (potentially) to some base set of Cartesian coordinates


    def __init__(self,
                 internals_spec,
                 *,
                 reference_coords=None,
                 converter=None,
                 inverse=None,
                 modes_matrix=None,
                 modes_inverse=None,
                 coordinate_expansion=None,
                 inverse_expansion=None,
                 **embedding_options
                 ):
        self.ints, self.embedding = self.canonicalize_internals(internals_spec, embedding_options)
        self.converter = TotalCoordinateConverter(
            converter, inverse, reference_coords,
            coordinate_expansion, inverse_expansion,
            modes_matrix, modes_inverse
        )
        super().__init__(
            "Total"+self.ints.name,
            dimension=(None, 3)
        )

    def canonicalize_internals(self, internals_spec, embedding):
        if isinstance(internals_spec, CoordinateSystem):
            embedding = dict(embedding, **internals_spec.converter_options)
        else:
            embedding['ordering'] = internals_spec
            internals_spec = ZMatrixCoordinates

        return internals_spec, embedding

class TotalCoordinatesToInternalCoordinateConverter(CoordinateSystemConverter):
    def __init__(self, totals:TotalCoordinateSystem, internals, **opts):
        self._types = (totals, internals)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types
    def convert(self, coords, **kw):
        return self._types[0].converter.invert(coords)

    def convert_many(self, coords, **kwargs):
        return self._types[0].converter.invert(coords)

class InternalCoordinateToTotalCoordinatesConverter(CoordinateSystemConverter):
    def __init__(self, totals: TotalCoordinateSystem, internals, **opts):
        self._types = (internals, totals)
        super().__init__(**opts)

    @property
    def types(self):
        return self._types

    def convert(self, coords, **kw):
        return self._types[0].converter.convert(coords)

    def convert_many(self, coords, **kwargs):
        return self._types[0].converter.convert(coords)
