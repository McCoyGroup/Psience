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
        """
        **LLM Docstring**

        Set up a converter that maps between a "total" coordinate representation (embedding coordinates plus a reduced set of internal coordinates) and a plain reduced internal-coordinate representation, chaining together an optional reference-shift, a linear mode transformation, a polynomial coordinate expansion, and a generic (possibly nonlinear) converter, applied in that order on the way in and reversed on the way out.

        :param converter: a generic forward-conversion function to apply to the coordinates, after the mode transformation and expansion
        :type converter: callable | None
        :param inverse: the inverse of `converter`; required if `converter` is given
        :type inverse: callable | None
        :param reference_coords: a reference geometry to shift the internal coordinates relative to
        :type reference_coords: np.ndarray | None
        :param coordinate_expansion: a Taylor-series coordinate expansion to apply after the mode transformation
        :type coordinate_expansion: list[np.ndarray] | None
        :param inverse_expansion: the inverse of `coordinate_expansion`; derived from it if not given
        :type inverse_expansion: list[np.ndarray] | None
        :param modes_matrix: a linear mode transformation matrix to apply to the (reference-shifted) internal coordinates
        :type modes_matrix: np.ndarray | None
        :param modes_inverse: the inverse of `modes_matrix`; derived from it if not given
        :type modes_inverse: np.ndarray | None
        :return: None
        :rtype: None
        :raises ValueError: if only one of `converter`/`inverse` is given
        """
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
        """
        **LLM Docstring**

        Fill in whichever of `modes_matrix`/`modes_inverse` wasn't supplied by inverting the other.

        :param modes_matrix: the mode transformation matrix, or `None` to derive it from `modes_inverse`
        :type modes_matrix: np.ndarray | None
        :param modes_inverse: the inverse mode transformation matrix, or `None` to derive it from `modes_matrix`
        :type modes_inverse: np.ndarray | None
        :return: `(modes_matrix, modes_inverse)`, both resolved (or both `None` if neither was given)
        :rtype: tuple[np.ndarray | None, np.ndarray | None]
        """
        if modes_matrix is not None and modes_inverse is None:
            modes_inverse = np.linalg.inv(modes_matrix)
        if modes_inverse is not None and modes_matrix is None:
            modes_matrix = np.linalg.inv(modes_inverse)

        return modes_matrix, modes_inverse

    @classmethod
    def prep_expansions(cls, coordinate_expansion, inverse_expansion):
        """
        **LLM Docstring**

        Fill in whichever of `coordinate_expansion`/`inverse_expansion` wasn't supplied by inverting the other via `nput.inverse_transformation`.

        :param coordinate_expansion: the forward Taylor-series coordinate expansion, or `None` to derive it from `inverse_expansion`
        :type coordinate_expansion: list[np.ndarray] | None
        :param inverse_expansion: the inverse Taylor-series expansion, or `None` to derive it from `coordinate_expansion`
        :type inverse_expansion: list[np.ndarray] | None
        :return: `(coordinate_expansion, inverse_expansion)`, both resolved (or both `None` if neither was given)
        :rtype: tuple[list[np.ndarray] | None, list[np.ndarray] | None]
        """
        if coordinate_expansion is not None and inverse_expansion is None:
            inverse_expansion = nput.inverse_transformation(coordinate_expansion, len(coordinate_expansion))
        if inverse_expansion is not None and coordinate_expansion is None:
            coordinate_expansion = nput.inverse_transformation(inverse_expansion, len(inverse_expansion))

        return coordinate_expansion, inverse_expansion

    def convert(self, internals):
        """
        **LLM Docstring**

        Convert reduced internal coordinates into the "total" coordinate representation: splits off the fixed embedding coordinates, shifts/transforms/expands/converts the remaining internal coordinates through the configured chain (reference shift, mode transformation, coordinate expansion, generic converter), and reassembles the result together with the untouched embedding coordinates into the Z-matrix-style `(natoms, 3)` total-coordinate array.

        :param internals: the reduced internal coordinates to convert
        :type internals: np.ndarray
        :return: the total-coordinate representation
        :rtype: np.ndarray
        """
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
        """
        **LLM Docstring**

        Convert the "total" coordinate representation back into reduced internal coordinates: splits off the fixed embedding coordinates, applies the inverse of the generic converter/coordinate expansion/mode transformation/reference shift (in reverse order to `convert`) to the remaining coordinates, and reassembles the result together with the untouched embedding coordinates.

        :param totals: the total-coordinate representation to invert
        :type totals: np.ndarray
        :return: the reduced internal coordinates
        :rtype: np.ndarray
        """
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
        """
        **LLM Docstring**

        Build a `TotalCoordinateSystem` wrapping a base internal-coordinate specification together with the `TotalCoordinateConverter` chain (reference shift, mode transformation, coordinate expansion, generic converter) used to map to/from a reduced internal representation.

        :param internals_spec: the base internal-coordinate specification (a `CoordinateSystem`, or a Z-matrix ordering)
        :type internals_spec: CoordinateSystem | np.ndarray
        :param reference_coords: a reference geometry to shift the internal coordinates relative to
        :type reference_coords: np.ndarray | None
        :param converter: a generic forward-conversion function to layer on top of the base internal coordinates
        :type converter: callable | None
        :param inverse: the inverse of `converter`
        :type inverse: callable | None
        :param modes_matrix: a linear mode transformation matrix
        :type modes_matrix: np.ndarray | None
        :param modes_inverse: the inverse of `modes_matrix`
        :type modes_inverse: np.ndarray | None
        :param coordinate_expansion: a Taylor-series coordinate expansion
        :type coordinate_expansion: list[np.ndarray] | None
        :param inverse_expansion: the inverse of `coordinate_expansion`
        :type inverse_expansion: list[np.ndarray] | None
        :param embedding_options: extra embedding options merged into the base internal-coordinate system's converter options
        :type embedding_options: dict
        :return: None
        :rtype: None
        """
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
        """
        **LLM Docstring**

        Normalize the `internals_spec` constructor argument into a `(coordinate_system, embedding_options)` pair: if already a `CoordinateSystem`, merges in its own converter options; otherwise treats `internals_spec` as a raw Z-matrix ordering and wraps it as `ZMatrixCoordinates`.

        :param internals_spec: the base internal-coordinate specification
        :type internals_spec: CoordinateSystem | np.ndarray
        :param embedding: the embedding options dict to merge into/extend
        :type embedding: dict
        :return: `(internals_spec, embedding)` -- the resolved coordinate-system class/instance and the merged embedding options
        :rtype: tuple[CoordinateSystem, dict]
        """
        if isinstance(internals_spec, CoordinateSystem):
            embedding = dict(embedding, **internals_spec.converter_options)
        else:
            embedding['ordering'] = internals_spec
            internals_spec = ZMatrixCoordinates

        return internals_spec, embedding

class TotalCoordinatesToInternalCoordinateConverter(CoordinateSystemConverter):
    def __init__(self, totals:TotalCoordinateSystem, internals, **opts):
        """
        **LLM Docstring**

        Store the `(totals, internals)` coordinate-system type pair this converter handles.

        :param totals: the source `TotalCoordinateSystem`
        :type totals: TotalCoordinateSystem
        :param internals: the target internal-coordinate system
        :type internals: CoordinateSystem
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (totals, internals)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(totals, internals)`
        :rtype: tuple
        """
        return self._types
    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Convert a single frame of total coordinates back to reduced internal coordinates, via the source system's own `TotalCoordinateConverter.invert`.

        :param coords: the total coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, unused
        :type kw: dict
        :return: the reduced internal coordinates
        :rtype: np.ndarray
        """
        return self._types[0].converter.invert(coords)

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Convert a batch of total coordinates back to reduced internal coordinates, via the source system's own `TotalCoordinateConverter.invert`.

        :param coords: the batch of total coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, unused
        :type kwargs: dict
        :return: the reduced internal coordinates
        :rtype: np.ndarray
        """
        return self._types[0].converter.invert(coords)

class InternalCoordinateToTotalCoordinatesConverter(CoordinateSystemConverter):
    def __init__(self, totals: TotalCoordinateSystem, internals, **opts):
        """
        **LLM Docstring**

        Store the `(internals, totals)` coordinate-system type pair this converter handles.

        :param totals: the target `TotalCoordinateSystem`
        :type totals: TotalCoordinateSystem
        :param internals: the source internal-coordinate system
        :type internals: CoordinateSystem
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (internals, totals)
        super().__init__(**opts)

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(internals, totals)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Convert a single frame of reduced internal coordinates to the total-coordinate representation, via the source system's own `TotalCoordinateConverter.convert`.

        :param coords: the internal coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, unused
        :type kw: dict
        :return: the total-coordinate representation
        :rtype: np.ndarray
        """
        return self._types[0].converter.convert(coords)

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Convert a batch of reduced internal coordinates to the total-coordinate representation, via the source system's own `TotalCoordinateConverter.convert`.

        :param coords: the batch of internal coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, unused
        :type kwargs: dict
        :return: the total-coordinate representation
        :rtype: np.ndarray
        """
        return self._types[0].converter.convert(coords)
