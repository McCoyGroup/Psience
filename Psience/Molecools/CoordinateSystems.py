"""
Defines useful extended internal coordinate frames
"""


import numpy as np
import McUtils.Devutils as dev
import McUtils.Numputils as nput
from McUtils.Parallelizers import Parallelizer
import McUtils.Coordinerds
from McUtils.Coordinerds import (
    ZMatrixCoordinateSystem, CartesianCoordinateSystem, CoordinateSystemConverter,
    ZMatrixCoordinates, CartesianCoordinates3D, CoordinateSet, CompositeCoordinateSystem,
    GenericInternalCoordinateSystem, GenericInternalCoordinates,
    CartesianToGICSystemConverter, GICSystemToCartesianConverter,
    PrimitiveCoordinatePicker, RedundantCoordinateGenerator
)

from ..Modes import NormalModes
from .Properties import StructuralProperties
# from .MoleculeInterface import AbstractMolecule

__all__ = [
    "MolecularEmbedding",
    "ModeEmbedding",
    "MolecularZMatrixCoordinateSystem",
    "MolecularCartesianCoordinateSystem"
]

__reload_hook__ = [".MoleculeInterface"]

class MolecularEmbedding:

    def __init__(self,
                 masses,
                 coords,
                 internals,
                 internal_fd_opts=None,
                 cartesian_fd_opts=None
                 ):
        """
        **LLM Docstring**

        Set up a molecule's coordinate embedding: wraps the Cartesian coordinates in a `MolecularCartesianCoordinateSystem`, canonicalizes the internal-coordinate specification (or stores it directly if already a `CoordinateSet`), and initializes the Jacobian cache and finite-difference option overrides.

        :param masses: the atomic masses
        :type masses: np.ndarray
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param internals: the internal-coordinate specification (Z-matrix, generic-internal specs, a callable conversion, or an already-built `CoordinateSet`)
        :type internals: object
        :param internal_fd_opts: overrides for the internal-coordinate finite-difference defaults
        :type internal_fd_opts: dict | None
        :param cartesian_fd_opts: overrides for the Cartesian finite-difference defaults
        :type cartesian_fd_opts: dict | None
        :return: None
        :rtype: None
        """

        self._coords = CoordinateSet(coords, MolecularCartesianCoordinateSystem(masses, coords))
        self.registered_converters = None
        if isinstance(internals, CoordinateSet):
            self._int_spec = None
            self._ints = internals
        else:
            self._int_spec = self.canonicalize_internal_coordinate_spec(internals)
            self._ints = None
        self._jacobians = self._get_jacobian_storage()
        # self._int_spec = internals
        self._frame = None
        self._tr_modes = None
        self._internal_fd_opts={} if internal_fd_opts is None else internal_fd_opts
        self._cart_fd_opts={} if cartesian_fd_opts is None else cartesian_fd_opts

    def get_direct_converter(self, target):
        """
        **LLM Docstring**

        Provide a converter from this embedding's coordinate system directly to plain (non-molecular) 3D Cartesian coordinates, if `target` is Cartesian-compatible.

        :param target: the coordinate system being converted to
        :type target: object
        :return: a `MolecularCartesianToRegularCartesianConverter`, or `None` if `target` isn't Cartesian-compatible
        :rtype: MolecularCartesianToRegularCartesianConverter | None
        """
        if target.is_compatible(CartesianCoordinates3D):
            carts = MolecularCartesianToRegularCartesianConverter(self.coords.system)
            return carts
    def get_inverse_converter(self, target):
        """
        **LLM Docstring**

        Provide a converter from plain (non-molecular) 3D Cartesian coordinates into this embedding's coordinate system, if `target` is Cartesian-compatible.

        :param target: the coordinate system being converted from
        :type target: object
        :return: a `RegularCartesianToMolecularCartesianConverter`, or `None` if `target` isn't Cartesian-compatible
        :rtype: RegularCartesianToMolecularCartesianConverter | None
        """
        if target.is_compatible(CartesianCoordinates3D):
            carts = RegularCartesianToMolecularCartesianConverter(self.coords.system)
            return carts

    def __del__(self):
        """
        **LLM Docstring**

        Deregister any coordinate converters this embedding registered, via `cleanup`, when the object is garbage collected.

        :return: None
        :rtype: None
        """
        self.cleanup()

    def cleanup(self):
        """
        **LLM Docstring**

        Deregister every converter previously registered via `register`, if any.

        :return: None
        :rtype: None
        """
        if hasattr(self, 'registered_converters') and self.registered_converters is not None:
            for d in self.registered_converters:
                d.deregister()

    def register(self):
        """
        **LLM Docstring**

        Register the Cartesian coordinate converters for this embedding's coordinate system with the global converter registry, if not already registered.

        :return: None
        :rtype: None
        """
        if not self.registered_converters:
            self.registered_converters = (
                MolecularCartesianToRegularCartesianConverter(self._coords.system),
                RegularCartesianToMolecularCartesianConverter(self._coords.system)
            )
            for r in self.registered_converters:
                r.register()
    @property
    def coords(self):
        """
        **LLM Docstring**

        Property getter/setter for the Cartesian coordinates. The getter registers the coordinate converters (via `register`) before returning them. The setter accepts a raw array or an already-systemed `CoordinateSet`, invalidates the Jacobian cache and inertial-frame cache, and marks the converters as needing re-registration if the coordinate system changed.

        :param coords: (setter only) the new Cartesian coordinates
        :type coords: np.ndarray | CoordinateSet
        :return: (getter) the Cartesian coordinates
        :rtype: CoordinateSet
        """
        self.register()
        return self._coords
    @coords.setter
    def coords(self, coords):
        """
        **LLM Docstring**

        Property getter/setter for the Cartesian coordinates. The getter registers the coordinate converters (via `register`) before returning them. The setter accepts a raw array or an already-systemed `CoordinateSet`, invalidates the Jacobian cache and inertial-frame cache, and marks the converters as needing re-registration if the coordinate system changed.

        :param coords: (setter only) the new Cartesian coordinates
        :type coords: np.ndarray | CoordinateSet
        :return: (getter) the Cartesian coordinates
        :rtype: CoordinateSet
        """
        if hasattr(coords, "system"):
            sys = coords.system
        else:
            sys = self._coords.system
            coords = CoordinateSet(coords, sys)
        if sys is not self._coords.system:
            self._regged = False
        self._jacobians = self._get_jacobian_storage()
        self._frame = None
        self._coords = coords
    @property
    def masses(self):
        """
        **LLM Docstring**

        The atomic masses, taken from the Cartesian coordinate system.

        :return: the atomic masses
        :rtype: np.ndarray
        """
        return self._coords.system.masses

    @staticmethod
    def _wrap_conv(f):
        """
        **LLM Docstring**

        Wrap a user-supplied coordinate-conversion function so it always returns a uniform `(values, opts)` pair: if the function returns a bare array, an empty options dict is paired with it; if it returns `(values, opts)`, the caller's original keyword arguments are merged underneath the returned `opts`.

        :param f: the conversion function to wrap, or `None`
        :type f: callable | None
        :return: the wrapped conversion function, or `None` if `f` is `None`
        :rtype: callable | None
        """
        if f is None:
            return f

        def wrapped(*args, **kwargs):
            """
            **LLM Docstring**

            Call the wrapped conversion function and normalize its result to a `(values, opts)` pair, merging the caller's keyword arguments underneath any options the function itself returned.

            :param args: positional arguments forwarded to the wrapped function
            :type args: tuple
            :param kwargs: keyword arguments forwarded to the wrapped function, and merged into the returned options
            :type kwargs: dict
            :return: `(values, opts)`
            :rtype: tuple
            """
            vals = f(*args, **kwargs)
            if not isinstance(vals, np.ndarray):
                vals, opts = vals
            else:
                opts = {}
            opts = dict(kwargs, **opts)
            return vals, opts

        return wrapped

    @classmethod
    def canonicalize_internal_coordinate_spec(cls, spec):
        """
        **LLM Docstring**

        Normalize the many accepted forms of an internal-coordinate specification (an options dict with `'zmatrix'`/`'specs'`/`'conversion'` keys, a bare callable conversion function, a raw Z-matrix-like array, or a list of generic-internal-coordinate specs) into the single canonical dict form (`'specs'`, `'zmatrix'`, `'conversion'`, `'inverse'`, `'converter_options'`) used internally, wrapping any conversion callables via `_wrap_conv` and filling in default embedding/jacobian-prep converter options for Z-matrices.

        :param spec: the internal-coordinate specification to canonicalize
        :type spec: dict | callable | Iterable | None
        :return: the canonicalized specification dict, or `None`/the original value if `spec` is `None`
        :rtype: dict | None
        :raises ValueError: if a Z-matrix-like spec doesn't have exactly 4 columns
        """
        if spec is not None:
            if hasattr(spec, 'items'):
                specs = spec.get('specs')
                if specs is None:
                    try:
                        zmatrix = spec['zmatrix']
                    except KeyError:
                        zmatrix = None
                    else:
                        zmatrix = np.asanyarray(zmatrix).astype(int)
                        if zmatrix.shape[1] != 4:
                            raise ValueError("can't understand Z-matrix {}".format(zmatrix))
                    spec['zmatrix'] = zmatrix
                    spec['conversion'] = cls._wrap_conv(spec.get('conversion'))
                    spec['inverse'] = cls._wrap_conv(spec.get('inverse'))
                    converter_options = spec.get('converter_options')
                    if converter_options is None:
                        converter_options = {}
                    if 'embedding_coords' not in converter_options:
                        if spec['zmatrix'] is not None:
                            converter_options['embedding_coords'] = MolecularZMatrixCoordinateSystem.embedding_coords
                    if 'jacobian_prep' not in converter_options:
                        if spec['zmatrix'] is not None:
                            converter_options['jacobian_prep'] = ZMatrixCoordinates.jacobian_prep_coordinates
                    spec['converter_options'] = converter_options
                else:
                    spec['zmatrix'] = None
                    spec['conversion'] = cls._wrap_conv(spec.get('conversion'))
                    spec['inverse'] = cls._wrap_conv(spec.get('inverse'))
            elif callable(spec):
                zmatrix = None
                specs = None
                conversion = cls._wrap_conv(spec)
                spec = {'specs': specs, 'zmatrix': zmatrix, 'conversion': conversion, 'inverse': None, 'converter_options': {}}
            else:
                if any(len(s) < 4 or hasattr(s, 'items') for s in spec):
                    conversion = None
                    specs = spec
                    zmatrix = None
                else:
                    conversion = None
                    specs = None
                    zmatrix = np.asanyarray(spec).astype(int)
                    if zmatrix.shape[1] != 4:
                        raise ValueError("can't understand Z-matrix {}".format(zmatrix))
                spec = {'specs': specs, 'zmatrix': zmatrix, 'conversion': conversion, 'inverse': None, 'converter_options': {}}
        return spec

    @property
    def internals(self):
        """
        **LLM Docstring**

        Property getter/setter for the raw (canonicalized) internal-coordinate specification. The setter re-canonicalizes the given specification and invalidates any already-computed internal coordinates.

        :param internals: (setter only) the new internal-coordinate specification, in any form accepted by `canonicalize_internal_coordinate_spec`
        :type internals: object
        :return: (getter) the canonicalized specification dict, or `None` if none is set
        :rtype: dict | None
        """
        if self._int_spec is not None:
            return self._int_spec

    @internals.setter
    def internals(self, internals):
        """
        **LLM Docstring**

        Property getter/setter for the raw (canonicalized) internal-coordinate specification. The setter re-canonicalizes the given specification and invalidates any already-computed internal coordinates.

        :param internals: (setter only) the new internal-coordinate specification, in any form accepted by `canonicalize_internal_coordinate_spec`
        :type internals: object
        :return: (getter) the canonicalized specification dict, or `None` if none is set
        :rtype: dict | None
        """
        self._int_spec = self.canonicalize_internal_coordinate_spec(internals)
        self._ints = None

    @property
    def zmatrix(self):
        """
        **LLM Docstring**

        Property getter/setter for just the Z-matrix ordering array out of the internal-coordinate specification. The setter validates the Z-matrix shape, builds a fresh specification if none exists yet, and invalidates any already-computed internal coordinates.

        :param zmat: (setter only) the new Z-matrix ordering array
        :type zmat: np.ndarray | None
        :return: (getter) the stored Z-matrix array, or `None`
        :rtype: np.ndarray | None
        :raises ValueError: if `zmat` doesn't have exactly 4 columns
        """
        if self._int_spec is not None:
            return self._int_spec['zmatrix']

    @zmatrix.setter
    def zmatrix(self, zmat):
        """
        **LLM Docstring**

        Property getter/setter for just the Z-matrix ordering array out of the internal-coordinate specification. The setter validates the Z-matrix shape, builds a fresh specification if none exists yet, and invalidates any already-computed internal coordinates.

        :param zmat: (setter only) the new Z-matrix ordering array
        :type zmat: np.ndarray | None
        :return: (getter) the stored Z-matrix array, or `None`
        :rtype: np.ndarray | None
        :raises ValueError: if `zmat` doesn't have exactly 4 columns
        """
        if zmat is not None:
            zmat = np.asanyarray(zmat).astype(int)
            if zmat.shape[1] != 4:
                raise ValueError("can't understand Z-matrix {}".format(zmat))
        if self._int_spec is None:
            self._int_spec = self.canonicalize_internal_coordinate_spec(zmat)
        else:
            self._int_spec['zmatrix'] = zmat
        self._ints = None

    @classmethod
    def convert_to_internals(cls, coords, masses, spec):
        """
        **LLM Docstring**

        Build the internal-coordinate `CoordinateSet` described by `spec`: constructs (and registers converters for) a generic-internal, Z-matrix, or iterative-Z-matrix coordinate system as appropriate, converts `coords` into it, layers on any extra custom `conversion`/`inverse` via a `CompositeCoordinateSystem`, and returns the resulting coordinates together with the (possibly updated, e.g. with redundant-transformation info) spec.

        :param coords: the Cartesian coordinates to convert
        :type coords: CoordinateSet
        :param masses: the atomic masses
        :type masses: np.ndarray
        :param spec: the canonicalized internal-coordinate specification (as produced by `canonicalize_internal_coordinate_spec`)
        :type spec: dict
        :return: `(coords, spec)` -- the internal coordinates and the (possibly updated) specification
        :rtype: tuple[CoordinateSet, dict]
        """
        spec = spec.copy()  # don't mutate user data
        opts = spec.copy()
        specs = opts.pop('specs', None)
        zmatrix = opts.pop('zmatrix', None)
        conversion = opts.pop('conversion', None)
        inverse = opts.pop('inverse', None)
        jacobian = opts.pop('jacobian', None)
        inverse_jacobian = opts.pop('inverse_jacobian', None)
        if specs is not None:
            ints = MolecularGenericInternalCoordinateSystem(masses, coords, specs=specs, **opts)
            ints.registered_converters = (
                MolecularCartesianToGICConverter(coords.system, ints),
                MolecularGICToCartesianConverter(ints, coords.system)
            )
            # for conv in ints.registered_converters:
            #     conv.register()
            coords = coords.convert(ints, reference_internals=spec.get('reference_internals'))
            cops = coords.converter_options
            for k in [
                'redundant_transformation',
                'redundant_inverse',
                'reference_internals',
                'reference_coordinates'
            ]:
                if k in cops:
                    spec[k] = cops[k]
        elif zmatrix is not None:
            iterative = opts.pop('iterative', False)
            if iterative:
                zms = MolecularIZCoordinateSystem(masses, coords, ordering=zmatrix, **opts)
                zms.registered_converters = (
                    MolecularCartesianToIZConverter(coords.system, zms),
                    MolecularIZToCartesianConverter(zms, coords.system),
                    MolecularIZToRegularIZConverter(zms),
                    RegularIZToMolecularIZConverter(zms)
                )
                # for conv in zms.registered_converters:
                #     conv.register()
            else:
                zms = MolecularZMatrixCoordinateSystem(masses, coords, ordering=zmatrix, **opts)
                zms.registered_converters = (
                    MolecularCartesianToZMatrixConverter(coords.system, zms),
                    MolecularZMatrixToCartesianConverter(zms, coords.system),
                    MolecularZMatrixToRegularZMatrixConverter(zms),
                    RegularZMatrixToMolecularZMatrixConverter(zms)
                )
                # for conv in zms.registered_converters:
                #     conv.register()
            coords = coords.convert(zms)
        if conversion is not None:
            conv = CompositeCoordinateSystem.register(
                coords.system,
                conversion,
                inverse_conversion=inverse,
                jacobian=jacobian,
                inverse_jacobian=inverse_jacobian
            )
            coords = coords.convert(conv)

        return coords, spec

    def internal_coordinates_from_spec(self, spec:dict):
        """
        **LLM Docstring**

        Build the internal coordinates for this embedding's current Cartesian coordinates and masses from a given specification, via `convert_to_internals`.

        :param spec: the canonicalized internal-coordinate specification
        :type spec: dict
        :return: `(coords, spec)`, as returned by `convert_to_internals`
        :rtype: tuple[CoordinateSet, dict]
        """
        return self.convert_to_internals(
            self._coords,
            self.masses,
            spec
        )

    @property
    def internal_coordinates(self):
        """
        **LLM Docstring**

        Property getter/setter for the internal coordinates. The getter lazily computes them from the stored specification (via `internal_coordinates_from_spec`) the first time they're needed. The setter requires an already-built `CoordinateSet`.

        :param ics: (setter only) the new internal coordinates
        :type ics: CoordinateSet
        :return: (getter) the internal coordinates, or `None` if no internal-coordinate specification is set
        :rtype: CoordinateSet | None
        :raises ValueError: if the setter is given something that isn't a `CoordinateSet`
        """
        if self._ints is None and (
                self._int_spec is not None
                and any(self._int_spec.get(k) is not None for k in {'zmatrix', 'conversion', 'specs'})
        ):
            self._ints, self._int_spec = self.internal_coordinates_from_spec(self._int_spec)
        return self._ints

    @internal_coordinates.setter
    def internal_coordinates(self, ics):
        """
        **LLM Docstring**

        Property getter/setter for the internal coordinates. The getter lazily computes them from the stored specification (via `internal_coordinates_from_spec`) the first time they're needed. The setter requires an already-built `CoordinateSet`.

        :param ics: (setter only) the new internal coordinates
        :type ics: CoordinateSet
        :return: (getter) the internal coordinates, or `None` if no internal-coordinate specification is set
        :rtype: CoordinateSet | None
        :raises ValueError: if the setter is given something that isn't a `CoordinateSet`
        """
        if not isinstance(ics, CoordinateSet):
            raise ValueError("{} must be a {} to be valid internal coordinates".format(
                ics, CoordinateSet.__name__
            ))
        self._ints = ics

    def strip_embedding_coordinates(self, coords):
        """
        **LLM Docstring**

        Drop the fixed embedding coordinates (e.g. the 6 translation/rotation degrees of freedom implied by the Z-matrix embedding) from a coordinate array or list of derivative tensors, if the underlying internal-coordinate system defines any.

        :param coords: the coordinates (or list of derivative tensors) to strip
        :type coords: np.ndarray | list[np.ndarray]
        :return: the coordinates with embedding coordinates removed, or unchanged if there are none to strip or they're already stripped
        :rtype: np.ndarray | list[np.ndarray]
        """
        embedding_coords = self._get_embedding_coords()
        if embedding_coords is not None:
            nc = 3 * len(self.masses)
            if isinstance(coords, list):
                if coords[0].shape[-1] == nc - len(embedding_coords): return coords
                good_coords = np.setdiff1d(np.arange(nc), embedding_coords)
                coords = [
                    c[..., good_coords]
                    for c in coords
                ]
            else:
                if coords.shape[-1] == nc - len(embedding_coords): return coords
                good_coords = np.setdiff1d(np.arange(nc), embedding_coords)
                coords = np.reshape(coords, coords.shape[:-2] + (-1,))
                coords = coords[..., good_coords]
        return coords

    def strip_derivative_embedding(self, derivs):
        """
        **LLM Docstring**

        Drop the fixed embedding coordinates from every axis of each tensor in a list of Cartesian-derivative tensors, if the underlying internal-coordinate system defines any.

        :param derivs: the list of Cartesian-derivative tensors (order-`n` tensor at index `n-1`) to strip
        :type derivs: list[np.ndarray]
        :return: the derivative tensors with embedding coordinates removed from every relevant axis, or unchanged if there are none to strip or they're already stripped
        :rtype: list[np.ndarray]
        """
        embedding_coords = self._get_embedding_coords()
        if embedding_coords is not None:
            nc = 3 * len(self.masses)
            if derivs[0].shape[-2] == nc - len(embedding_coords): return derivs
            good_coords = np.setdiff1d(np.arange(nc), embedding_coords)
            new_derivs = []
            for n,t in enumerate(derivs):
                idx = (...,) + np.ix_(*((good_coords,) * (n + 1))) + (slice(None),)
                new_derivs.append(t[idx])
            derivs = new_derivs
        return derivs

    def restore_embedding_coordinates(self, coords):
        """
        **LLM Docstring**

        Reinsert the fixed embedding coordinates (filled in from the reference internal coordinates) back into a stripped coordinate array or list, undoing `strip_embedding_coordinates`.

        :param coords: the stripped coordinates (or list) to restore
        :type coords: np.ndarray | list[np.ndarray]
        :return: the coordinates with embedding coordinates reinserted, or unchanged if there are none to restore or they're already present
        :rtype: np.ndarray | list[np.ndarray]
        """
        embedding_coords = self._get_embedding_coords()
        if embedding_coords is not None:
            nc = 3 * len(self.masses)
            if isinstance(coords, list):
                #TODO: add more robust check for derivatives
                #TODO: this is probably invalid and will need to be overwritten by downstream users
                if coords[0].shape[-1] == nc: return coords
                new_coords = []
                good_coords = np.setdiff1d(np.arange(nc), embedding_coords)
                for c in coords:
                    new_c = np.zeros(c.shape[:-1] + (nc,), dtype=c.dtype)
                    new_c[..., good_coords] = c
                    new_coords.append(c)
                coords = new_coords
            else:
                if coords.shape[-1] == nc: return coords
                coord_shape = coords.shape[:-1]
                fc_shape = self.internal_coordinates.shape
                new_coords = type(self.internal_coordinates)(
                    np.zeros(coord_shape + fc_shape, dtype=self.internal_coordinates.dtype),
                    self.internal_coordinates.system,
                    converter_options=self.internal_coordinates.converter_options
                )
                flat_coords = new_coords.reshape(coord_shape + (-1,))
                good_coords = np.setdiff1d(np.arange(np.prod(fc_shape, dtype=int)), embedding_coords)
                flat_coords[..., good_coords] = coords
                flat_coords[..., embedding_coords] = self.internal_coordinates.flatten()[embedding_coords,][np.newaxis]
                coords = new_coords
        return coords

    def restore_derivative_embedding(self, derivs):
        """
        **LLM Docstring**

        Reinsert zeroed-out placeholder entries for the fixed embedding coordinates back into every axis of each tensor in a list of stripped Cartesian-derivative tensors, undoing `strip_derivative_embedding`.

        :param derivs: the stripped list of Cartesian-derivative tensors to restore
        :type derivs: list[np.ndarray]
        :return: the derivative tensors with embedding-coordinate axes reinserted (as zeros), or unchanged if there are none to restore or they're already present
        :rtype: list[np.ndarray]
        """
        embedding_coords = self._get_embedding_coords()
        if embedding_coords is not None:
            nc = 3 * len(self.masses)
            if derivs[0].shape[-2] == nc: return derivs
            good_coords = np.setdiff1d(np.arange(nc), embedding_coords)
            coord_shape = derivs[0].shape[:-2]
            new_derivs = []
            for n,t in enumerate(derivs):
                new_d = np.zeros(
                    coord_shape + (nc,) * (n+2),
                    dtype=self.internal_coordinates.dtype
                )
                idx = (...,) + np.ix_(*((good_coords,) * (n + 1))) + (slice(None),)
                new_d[idx] = t
                new_derivs.append(new_d)
            derivs = new_derivs
        return derivs

    def get_internals(self, *, coords=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the internal coordinates, either the cached ones for this embedding's own geometry or freshly computed ones for an alternate set of Cartesian `coords`, optionally stripping the fixed embedding coordinates.

        :param coords: alternate Cartesian coordinates to convert instead of using the cached internal coordinates
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the internal coordinates, or `None` if none are defined
        :rtype: CoordinateSet | None
        """
        if coords is None:
            ics = self.internal_coordinates
        else:
            ccoords = type(self.coords)(coords, self.coords.system)
            ics, _ = self.convert_to_internals(ccoords, self.masses, self.internals)
        if ics is None:
            return None
        if strip_embedding:
            ics = self.strip_embedding_coordinates(ics)
        return ics
    def get_cartesians(self, *, coords=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch Cartesian coordinates, either this embedding's own cached ones or the Cartesian coordinates corresponding to a given set of internal `coords`, optionally restoring any stripped embedding coordinates first.

        :param coords: internal-coordinate values to convert to Cartesians instead of returning the cached Cartesian coordinates
        :type coords: np.ndarray | None
        :param strip_embedding: whether `coords` has had its embedding coordinates stripped and needs them restored before conversion
        :type strip_embedding: bool
        :return: the Cartesian coordinates
        :rtype: CoordinateSet
        """
        if coords is None:
            return self.coords
        else:
            if strip_embedding:
                coords = self.restore_embedding_coordinates(coords)
            else:
                coords = type(self.internal_coordinates)(coords,
                                                         self.internal_coordinates.system,
                                                         converter_options=self.internal_coordinates.converter_options
                                                         )
            return coords.convert(self.coords.system)


    @property
    def redundant_internal_transformation(self):
        """
        **LLM Docstring**

        The redundant-to-non-redundant transformation matrix associated with the current internal coordinates, if the internal-coordinate system used a redundant coordinate generator.

        :return: the redundant transformation, or `None` if not applicable
        :rtype: np.ndarray | None
        """
        return self.internal_coordinates.converter_options.get('redundant_transformation')

    @classmethod
    def _get_jacobian_storage(cls):
        """
        **LLM Docstring**

        Build a fresh, empty nested dict used to cache computed Jacobian tensors by coordinate-system pair and (for internals) by re-embedding mode.

        :return: the empty Jacobian cache structure
        :rtype: dict
        """
        return {
            'internals': {"default":[], "reembed":[]},
            'fast-internals': {"default":[], "reembed":[], "reembed-strip":[]},
            'cartesian': []
        }
    internal_fd_defaults=dict(
        strip_dummies=False,
        stencil=None,
        mesh_spacing=1.0e-3,
        all_numerical=None,
        reembed=True,
        planar_ref_tolerance=None,
        parallelizer=None
    )
    def _get_internal_fd_opts(self, **opts):
        """
        **LLM Docstring**

        Merge the class-level `internal_fd_defaults`, this instance's `_internal_fd_opts` overrides, and any explicitly passed `opts`, with later sources taking precedence, to resolve the finite-difference options used for internal-coordinate Jacobians.

        :param opts: explicit per-call overrides, taking highest precedence
        :type opts: dict
        :return: the merged finite-difference options
        :rtype: dict
        """
        return dict(
            dict(
                self.internal_fd_defaults,
                **self._internal_fd_opts
            ),
            **opts
        )
    def _get_int_jacobs(self,
                        jacs,
                        coords=None,
                        strip_embedding=False,
                        **fd_opts
                        ):
        """
        Gets the specified dX/dRs

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """

        if coords is None:
            intcds = self.internal_coordinates
            ccoords = self.coords
        else:
            # intcds = self.internal_coordinates
            ccoords = type(self.coords)(coords, self.coords.system)
            intcds, _ = self.convert_to_internals(ccoords, self.masses, self.internals)
        # intcds = self.internal_coordinates
        # ccoords = self.coords
        carts = ccoords.system
        internals = intcds.system

        generics = 'GenericInternals' in internals.name
        zmatrix = 'ZMatrix' in internals.name

        if "analytic_derivative_order" in fd_opts:
            fd_opts['analytic_deriv_order'] = fd_opts.pop('analytic_derivative_order')
        fd_opts['analytic_deriv_order'] = fd_opts.get('analytic_deriv_order', -1)
        fd_opts = self._get_internal_fd_opts(**fd_opts)
        (
            all_numerical, reembed, planar_ref_tolerance,
            strip_dummies, mesh_spacing, stencil, parallelizer,
            analytic_deriv_order
        ) = [
            fd_opts.pop(k) for k in
            (
                "all_numerical", "reembed", "planar_ref_tolerance",
                "strip_dummies", "mesh_spacing", "stencil", "parallelizer",
                "analytic_deriv_order"
            )
        ]

        if all_numerical is None:
            all_numerical = zmatrix and (analytic_deriv_order is not None and analytic_deriv_order == 0)

        converter_options = (
            dict(
                reembed=reembed,
                planar_ref_tolerance=planar_ref_tolerance,
                strip_dummies=strip_dummies,
                strip_embedding=strip_embedding
            )
                if zmatrix else
            {}
        )

        if isinstance(jacs, int):
            jacs = list(range(1, jacs + 1))
        max_jac = max(jacs)

        if coords is None:

            exist_jacs = self._jacobians['internals']
            if reembed:
                exist_jacs = exist_jacs["reembed"]
            else:
                exist_jacs = exist_jacs["default"]
            # {"default": [], "reembed": [], "fast": [], "generic": []}
            need_jacs = [x + 1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
            if len(need_jacs) > 0:
                if generics:
                    exist_jacs[:] = [
                        x #.squeeze() for x in
                        for x in
                        intcds.jacobian(carts, jacs,
                                        # odd behaves better
                                        mesh_spacing=mesh_spacing,
                                        stencil=stencil,
                                        all_numerical=all_numerical,
                                        converter_options=converter_options,
                                        analytic_deriv_order=analytic_deriv_order,
                                        **fd_opts
                                        )
                    ]
                    # self._jacobians['internals'] = exist_jacs
                else:
                    stencil = (max(need_jacs) + 2 + (1 + max(need_jacs)) % 2) if stencil is None else stencil
                    # odd behaves better
                    with Parallelizer.lookup(parallelizer) as par:
                        new_jacs = [
                            x#.squeeze() if isinstance(x, np.ndarray) else x
                            for x in intcds.jacobian(carts, need_jacs,
                                                     # odd behaves better
                                                     mesh_spacing=mesh_spacing,
                                                     stencil=stencil,
                                                     all_numerical=all_numerical,
                                                     converter_options=converter_options,
                                                     parallelizer=par,
                                                     analytic_deriv_order=analytic_deriv_order,
                                                     **fd_opts
                                                     )
                        ]
                    for j, v in zip(need_jacs, new_jacs):
                        for d in range(j - len(exist_jacs)):
                            exist_jacs.append(None)
                        exist_jacs[j - 1] = v

            return [exist_jacs[j - 1] for j in jacs]
        else:
            if generics:
                new_jacs = [
                    x for x in #.squeeze() for x in
                    intcds.jacobian(carts, jacs,
                                    # odd behaves better
                                    mesh_spacing=mesh_spacing,
                                    stencil=stencil,
                                    all_numerical=all_numerical,
                                    converter_options=converter_options,
                                    analytic_deriv_order=analytic_deriv_order,
                                    **fd_opts
                                    )
                ]
                # self._jacobians['internals'] = exist_jacs
            else:
                stencil = (max_jac + 2 + (1 + max_jac) % 2) if stencil is None else stencil
                # odd behaves better
                with Parallelizer.lookup(parallelizer) as par:
                    new_jacs = [
                        # x.squeeze() if isinstance(x, np.ndarray) else x
                        x
                        for x in intcds.jacobian(carts, max_jac,
                                                 # odd behaves better
                                                 mesh_spacing=mesh_spacing,
                                                 stencil=stencil,
                                                 all_numerical=all_numerical,
                                                 converter_options=converter_options,
                                                 parallelizer=par,
                                                 analytic_deriv_order=analytic_deriv_order,
                                                 **fd_opts
                                                 )
                    ]
            return new_jacs

    cart_fd_defaults=dict(
        strip_dummies=False,
        stencil=None,
        mesh_spacing=1.0e-3,
        all_numerical=None,
        parallelizer=None,
        analytic_deriv_order=-1
    )
    def _get_cart_fd_opts(self, **opts):
        """
        **LLM Docstring**

        Merge the class-level `cart_fd_defaults`, this instance's `_cart_fd_opts` overrides, and any explicitly passed `opts`, with later sources taking precedence, to resolve the finite-difference options used for Cartesian-coordinate Jacobians.

        :param opts: explicit per-call overrides, taking highest precedence
        :type opts: dict
        :return: the merged finite-difference options
        :rtype: dict
        """
        return dict(
            dict(
                self.cart_fd_defaults,
                **self._cart_fd_opts
            ),
            **opts
        )
    def _get_cart_jacobs(self, jacs,
                         coords=None,
                         strip_embedding=False,
                         **fd_opts
                         ):
        """
        Gets the specified dR/dXs

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        if coords is None:
            # intcds = self.internal_coordinates
            ccoords = self.coords
        else:
            # intcds = self.internal_coordinates
            ccoords = type(self.coords)(coords, self.coords.system)
        carts = ccoords.system
        internals = self.internal_coordinates.system
        generics = 'GenericInternals' in internals.name
        zmatrix = 'ZMatrix' in internals.name

        if "analytic_derivative_order" in fd_opts:
            fd_opts['analytic_deriv_order'] = fd_opts.pop('analytic_derivative_order')
        fd_opts = self._get_cart_fd_opts(**fd_opts)
        (
            all_numerical, strip_dummies, mesh_spacing, stencil, parallelizer,
            analytic_deriv_order
        ) = [
            fd_opts.pop(k) for k in
            (
                "all_numerical", "strip_dummies", "mesh_spacing", "stencil", "parallelizer",
                "analytic_deriv_order"
            )
        ]
        if zmatrix and not strip_embedding:
            analytic_deriv_order = 0

        if all_numerical is None:
            all_numerical = zmatrix and (analytic_deriv_order is not None and analytic_deriv_order == 0)

        converter_options = (
            dict(strip_dummies=strip_dummies)
                if zmatrix else
            {}
        )

        if isinstance(jacs, int):
            jacs = list(range(1, jacs + 1))

        if coords is None:
            exist_jacs = self._jacobians['cartesian']
            max_jac = max(jacs)
            need_jacs = [x + 1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
            if len(need_jacs) > 0:
                if generics:
                    exist_jacs[:] = [
                        x.squeeze() for x in
                        ccoords.jacobian(internals,
                                         list(range(1, max_jac+1)),
                                         # mesh_spacing=mesh_spacing,
                                         # stencil=stencil,
                                         # all_numerical=all_numerical,
                                         analytic_deriv_order=analytic_deriv_order,
                                         converter_options=converter_options,
                                         **fd_opts
                                         )
                    ]
                else:
                    stencil = (max(need_jacs) + 2 + (1 + max(need_jacs)) % 2) if stencil is None else stencil
                    # odd behaves better
                    with Parallelizer.lookup(parallelizer) as par:
                        exist_jacs[:] = [
                            x.squeeze() if isinstance(x, np.ndarray) else x
                            for x in ccoords.jacobian(internals,
                                                      order=list(range(1, max_jac + 1)),
                                                      mesh_spacing=mesh_spacing,
                                                      stencil=stencil,
                                                      all_numerical=all_numerical,
                                                      converter_options=converter_options,
                                                      analytic_deriv_order=analytic_deriv_order,
                                                      parallelizer=par,
                                                      **fd_opts
                                                      )
                        ]

            return [exist_jacs[j - 1] for j in jacs]
        else:
            max_jac = max(jacs)
            if generics:
                new_jacs = [
                    x.squeeze() for x in
                    ccoords.jacobian(internals,
                                     list(range(1, max_jac + 1)),
                                     # mesh_spacing=mesh_spacing,
                                     # stencil=stencil,
                                     # all_numerical=all_numerical,
                                     analytic_deriv_order=analytic_deriv_order,
                                     converter_options=converter_options
                                     )
                ]
            else:
                stencil = (max_jac + 2 + (1 + max_jac) % 2) if stencil is None else stencil
                # odd behaves better
                with Parallelizer.lookup(parallelizer) as par:
                    new_jacs = [
                        x.squeeze() if isinstance(x, np.ndarray) else x
                        for x in ccoords.jacobian(internals,
                                                  order=list(range(1, max_jac + 1)),
                                                  mesh_spacing=mesh_spacing,
                                                  stencil=stencil,
                                                  all_numerical=all_numerical,
                                                  analytic_deriv_order=analytic_deriv_order,
                                                  converter_options=converter_options,
                                                  parallelizer=par
                                                  )
                        ]
            return new_jacs

    @property
    def embedding_coords(self):
        """
        **LLM Docstring**

        The indices of the internal-coordinate system's fixed embedding coordinates (translation/rotation degrees of freedom), if any are defined.

        :return: the embedding-coordinate indices, or `None`
        :rtype: np.ndarray | None
        """
        return self._get_embedding_coords()
    def _get_embedding_coords(self):
        """
        **LLM Docstring**

        Look up the embedding-coordinate indices from the internal-coordinate system, trying its `embedding_coords` attribute first and falling back to its `converter_options['embedding_coords']`.

        :return: the embedding-coordinate indices, or `None` if neither source defines them
        :rtype: np.ndarray | None
        """
        try:
            embedding = self.internal_coordinates.system.embedding_coords
        except AttributeError:
            try:
                embedding = self.internal_coordinates.system.converter_options['embedding_coords']
            except KeyError:
                embedding = None
        return embedding

    cartesian_by_internals_method = 'fast'
    def get_cartesians_by_internals(self, order=None, strip_embedding=False,
                                    reembed=True, method=None, coords=None,
                                    **fd_opts
                                    ):
        """
        **LLM Docstring**

        Compute the Cartesians-by-internals Jacobian expansion up to the requested `order`, either via the fast route (inverting the internals-by-Cartesians Jacobian through a translation/rotation-invariant reduction, with caching) or the classic finite-difference/analytic route, depending on `method` (auto-selected based on the internal-coordinate system type) and whether `reembed`/`strip_embedding`/explicit `coords` are requested.

        :param order: the highest derivative order to compute; if `None`, returns whatever is already cached
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param reembed: whether to use the Eckart-reembedded (translation/rotation-invariant) formulation, for the `'fast'` method
        :type reembed: bool
        :param method: which computation strategy to use (`'fast'` or `'classic'`); auto-selected if `None`
        :type method: str | None
        :param coords: alternate Cartesian coordinates to compute the Jacobian at, instead of this embedding's own geometry
        :type coords: np.ndarray | None
        :param fd_opts: extra finite-difference options forwarded to the underlying Jacobian computation
        :type fd_opts: dict
        :return: the Cartesians-by-internals Jacobian tensors, one per order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached/computed orders are available than requested (classic route)
        """
        if method is None:
            int_sys = self.internal_coordinates.system
            if "GenericInternals" in int_sys.name:
                method = 'classic'
            else:
                method = self.cartesian_by_internals_method

        if method == 'fast':
            if reembed:
                if coords is None:
                    if strip_embedding:
                        key = 'reembed-strip'
                    else:
                        key = 'reembed'
                    fast_ints = self._jacobians['fast-internals'][key]
                else:
                    fast_ints = []
                if len(fast_ints) < order:
                    L_base, L_inv = self.get_translation_rotation_invariant_transformation(
                        strip_embedding=strip_embedding,
                        mass_weighted=False,
                        coords=coords
                    )
                    jacs_1 = self.get_internals_by_cartesians(order, strip_embedding=strip_embedding, coords=coords)
                    new_tf = nput.tensor_reexpand([L_inv], jacs_1, axes=[-1, -2], order=order)
                    inverse_tf = nput.inverse_transformation(new_tf, order, allow_pseudoinverse=False)
                    fast_ints[:] = [
                        nput.vec_tensordot(j, L_inv, axes=[-1, -2], shared=L_base.ndim-2)
                        for j in inverse_tf
                    ]
                return fast_ints
            else:
                wtf = self.get_internals_by_cartesians(order, strip_embedding=False, coords=coords) # faster to just do these derivs.
                base = nput.inverse_transformation(wtf, order, allow_pseudoinverse=True)
        else:
            if coords is None:
                base = (
                    self._get_int_jacobs(order, reembed=reembed, strip_embedding=strip_embedding, **fd_opts)
                        if order is not None else
                    self._jacobians['internals']["default" if not reembed else "reembed"]
                )
                if order is not None:
                    if len(base) < order:
                        raise ValueError("insufficient {} (have {} but expected {})".format(
                            'CartesiansByInternals',
                            len(base),
                            order
                        ))
                    base = base[:order]
            else:
                coords = np.asanyarray(coords)
                if order is None: order = 1
                base = self._get_int_jacobs(order, reembed=reembed, strip_embedding=strip_embedding, coords=coords, **fd_opts)

            # raise Exception(fd_opts.get('use_direct_expansions', True), base[0].shape)
            if not dev.is_list_like(fd_opts.get('use_direct_expansions', True)):
                _ = []
                if coords is not None:
                    sh = coords.shape[:-2]
                else:
                    sh = self.coords.shape[:-2]
                nc = 3 * len(self.masses)
                nr = len(sh)
                n = -2 if base[0].shape[-2:] == (len(self.masses), 3) else -1
                for i, b in enumerate(base):
                    rem = round(
                        np.power(np.prod(b.shape[nr:n], dtype=int), 1 / (i+1))
                    )
                    b = b.reshape(sh + (rem,) * (i+1) + (nc,))
                    _.append(b)
                base = _

                embedding_coords = self._get_embedding_coords() if strip_embedding else None
                if embedding_coords is not None and strip_embedding and base[0].shape[-2] == (3 * len(self.masses)):
                    good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
                    base = [t[np.ix_(*((good_coords,) * (t.ndim - 1)))] for t in base]
        return base

    def get_internals_by_cartesians(self, order=None, strip_embedding=False, coords=None, **opts):
        """
        **LLM Docstring**

        Compute the internals-by-Cartesians Jacobian expansion up to the requested `order`, via finite difference/analytic derivatives (through `_get_cart_jacobs`), reshaping the results to `(..., ncart, ncart, ..., nint)`-style tensors and optionally stripping the fixed embedding coordinates.

        :param order: the highest derivative order to compute; if `None`, returns whatever is already cached
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param coords: alternate Cartesian coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param opts: extra finite-difference options forwarded to `_get_cart_jacobs`
        :type opts: dict
        :return: the internals-by-Cartesians Jacobian tensors, one per order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached/computed orders are available than requested
        """
        if coords is not None:
            coords = np.asanyarray(coords)
            if order is None: order = 1
        base = (
                   self._get_cart_jacobs(order, coords=coords, strip_embedding=strip_embedding, **opts)
                        if order is not None else
                   self._jacobians['cartesian']
        )
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'InternalsByCartesians',
                    len(base),
                    order
                ))
            base = base[:order]

        _ = []
        if coords is not None:
            sh = coords.shape[:-2]
        else:
            sh = self.coords.shape[:-2]
        nc = 3 * len(self.masses)
        for i, b in enumerate(base):
            b = b.reshape(sh + (nc,) * (i + 1) + (-1,))
            _.append(b)
        base = _

        if strip_embedding:
            base = self.strip_embedding_coordinates(base)
        return base

    def embed_coords(self, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=True):
        """
        **LLM Docstring**

        Eckart-embed a set of Cartesian coordinates onto this embedding's reference geometry.

        :param coords: the coordinates to embed
        :type coords: np.ndarray
        :param sel: subset of atoms to use for the embedding fit
        :type sel: Iterable[int] | None
        :param in_paf: whether to embed into the principal-axis frame
        :type in_paf: bool
        :param planar_ref_tolerance: tolerance for detecting a (near-)planar reference structure
        :type planar_ref_tolerance: float | None
        :param proper_rotation: whether to restrict the embedding to proper (determinant +1) rotations
        :type proper_rotation: bool
        :return: the Eckart-embedded coordinates
        :rtype: np.ndarray
        """
        return StructuralProperties.get_eckart_embedded_coords(
            self.masses,
            self.coords,
            coords,
            sel=sel,
            in_paf=in_paf,
            planar_ref_tolerance=planar_ref_tolerance,
            proper_rotation=proper_rotation
        )

    def unembed_derivs(self, coords, derivs, sel=None, in_paf=False, planar_ref_tolerance=None):
        """
        **LLM Docstring**

        Undo an Eckart embedding's rotation on a set of Cartesian derivative tensors, transforming them back by the combination of the embedding's axis frame and rotation.

        :param coords: the (embedded) coordinates the derivatives were computed at
        :type coords: np.ndarray
        :param derivs: the Cartesian derivative tensors to un-rotate
        :type derivs: list[np.ndarray]
        :param sel: subset of atoms used for the embedding fit
        :type sel: Iterable[int] | None
        :param in_paf: whether the embedding used the principal-axis frame
        :type in_paf: bool
        :param planar_ref_tolerance: tolerance for detecting a (near-)planar reference structure
        :type planar_ref_tolerance: float | None
        :return: the un-rotated derivative tensors
        :rtype: list[np.ndarray]
        """
        emb = nput.eckart_embedding(
            self.coords,
            coords,
            self.masses,
            sel=sel,
            in_paf=in_paf,
            planar_ref_tolerance=planar_ref_tolerance,
            transform_coordinates=False
        )

        tfs = emb.coord_data.axes @ emb.rotations
        return nput.transform_cartesian_derivatives(
            derivs,
            tfs
        )


    @property
    def inertia_tensor(self):
        """
        **LLM Docstring**

        The molecule's inertia tensor at its current Cartesian coordinates.

        :return: the inertia tensor
        :rtype: np.ndarray
        """
        return StructuralProperties.get_prop_inertia_tensors(
            self.coords,
            self.masses
        )

    @property
    def inertial_frame(self):
        """
        Provides the inertial axis frame

        :return:
        :rtype:
        """

        if self._frame is None:
            # Need to put B in Hartree?
            #  I've got moments of inertia in amu * bohr^2 at the moment
            #  So we convert (amu * bohr^2) to (m_e * bohr^2) since hb^2/(m_e bohr^2) == E_h
            mom_i, eigs = StructuralProperties.get_prop_moments_of_inertia(
                self.coords,
                self.masses
            )
            B_e = 1 / (2 * mom_i)  # * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass"))
            self._inert_frame = B_e, eigs

        return self._inert_frame

    def inertial_frame_derivatives(self):
        """
        **LLM Docstring**

        The first and second derivatives of the inertia tensor with respect to mass-weighted Cartesian displacements.

        :return: `[I0Y, I0YY]`, as returned by `StructuralProperties.get_prop_inertial_frame_derivatives`
        :rtype: list[np.ndarray]
        """

        return StructuralProperties.get_prop_inertial_frame_derivatives(
                self.coords,
                self.masses
            )

    @property
    def translation_rotation_modes(self):
        """
        **LLM Docstring**

        The (cached) translation and rotation eigenvectors of the molecule at its current geometry.

        :return: the translation/rotation eigenvectors
        :rtype: tuple
        """
        if self._tr_modes is None:
            self._tr_modes = StructuralProperties.get_prop_translation_rotation_eigenvectors(
                self.coords,
                self.masses
            )

        return self._tr_modes


    def get_translation_rotation_invariant_transformation(self,
                                                          order=0,
                                                          mass_weighted=True,
                                                          strip_embedding=True,
                                                          coords=None):
        """
        **LLM Docstring**

        Build the transformation (and its inverse) that projects out the translational and rotational degrees of freedom from a set of Cartesian coordinates.

        :param order: the derivative order of the transformation to build
        :type order: int
        :param mass_weighted: whether the transformation should act on mass-weighted coordinates
        :type mass_weighted: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param coords: alternate Cartesian coordinates to build the transformation at, instead of this embedding's own geometry
        :type coords: np.ndarray | None
        :return: the translation/rotation-invariant transformation and its inverse
        :rtype: tuple
        """
        return nput.translation_rotation_invariant_transformation(
                self.coords if coords is None else coords,
                masses=self.masses,
                mass_weighted=mass_weighted,
                strip_embedding=strip_embedding
            )

class ModeEmbedding:
    """
    Provides a specialization on a `MoleculaEmbedding` to express all properties
    in terms of the attendant normal modes
    """
    modes: NormalModes
    def __init__(self,
                 embedding: MolecularEmbedding,
                 modes,  #:#NormalModesManager,
                 mass_weight=None,
                 dimensionless=None,
                 masses=None
                 ):
        """
        **LLM Docstring**

        Set up a mode-aware specialization of a `MolecularEmbedding`: resolves whatever form `modes` was given in (a manager, a `MolecularVibrations`, or a raw normal-modes object) down to a plain normal-modes object, optionally converting it to a dimensionless or mass-weighted basis, and records whether the resulting modes are mass-weighted.

        :param embedding: the underlying molecular embedding to specialize
        :type embedding: MolecularEmbedding
        :param modes: the normal modes (or a manager/vibrations object wrapping them) to express properties in terms of
        :type modes: object
        :param mass_weight: whether to convert the modes to a mass-weighted basis
        :type mass_weight: bool | None
        :param dimensionless: whether to convert the modes to a dimensionless basis
        :type dimensionless: bool | None
        :param masses: masses to use for the mass-weighting/dimensionless conversions; defaults to the embedding's own masses
        :type masses: np.ndarray | None
        :return: None
        :rtype: None
        """

        self.embedding = embedding
        self.masses = masses
        if hasattr(modes, 'get_modes'):
            modes = modes.get_modes(quiet=True, allow_compute=False)
        if hasattr(modes, 'modes'):
            modes = modes.modes
        if hasattr(modes, 'basis') and hasattr(modes.basis, 'to_new_modes'):
            modes = modes.basis.to_new_modes()
        if modes is not None:
            if dimensionless:
                modes = modes.make_dimensionless(masses=masses)
            elif mass_weight:
                modes = modes.make_mass_weighted(masses=masses)

        if modes is not None:
            self.mass_weighted = modes.mass_weighted
        else:
            self.mass_weighted = True
        self.modes = modes

    def mw_conversion(self, strip_dummies=None):
        """
        **LLM Docstring**

        Build the diagonal mass-weighting matrix (`sqrt(mass)` per Cartesian coordinate) used to convert plain Cartesian derivatives into mass-weighted ones.

        :param strip_dummies: whether to exclude dummy (non-positive-mass) atoms from the mass vector
        :type strip_dummies: bool | None
        :return: the diagonal mass-weighting matrix
        :rtype: np.ndarray
        """
        masses = self.masses
        if masses is None:
            masses = self.embedding.masses
        if strip_dummies:
            masses = masses[masses > 0]
        mvec = np.broadcast_to(
                np.asanyarray(masses)[:, np.newaxis],
                (len(masses), 3)
            ).flatten()
        return np.diag(np.sign(mvec) * np.sqrt(np.abs(mvec)))
    def mw_inverse(self, strip_dummies=None):
        """
        **LLM Docstring**

        Build the diagonal inverse-mass-weighting matrix (`1/sqrt(mass)` per Cartesian coordinate) used to convert mass-weighted Cartesian derivatives back into plain ones.

        :param strip_dummies: whether to exclude dummy (non-positive-mass) atoms from the mass vector
        :type strip_dummies: bool | None
        :return: the diagonal inverse-mass-weighting matrix
        :rtype: np.ndarray
        """
        masses = self.masses
        if masses is None:
            masses = self.embedding.masses
        if strip_dummies:
            masses = masses[masses > 0]
        mvec = np.broadcast_to(
                np.asanyarray(masses)[:, np.newaxis],
                (len(masses), 3)
            ).flatten()
        return np.diag(np.sign(mvec) / np.sqrt(np.abs(mvec)))

    def get_mw_cartesians_by_internals(self, order=None, mass_weighted=None, coords=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-internals Jacobian expansion, converted to mass-weighted Cartesians if `mass_weighted` (or `self.mass_weighted` by default) is set.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param mass_weighted: whether to mass-weight the result; defaults to `self.mass_weighted`
        :type mass_weighted: bool | None
        :param coords: alternate coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates
        :type strip_embedding: bool
        :return: the (optionally mass-weighted) Cartesians-by-internals Jacobian tensors
        :rtype: list[np.ndarray]
        """
        RX = self.embedding.get_cartesians_by_internals(
            order=order,
            coords=coords,
            strip_embedding=strip_embedding
        )
        if mass_weighted is None:
            mass_weighted = self.mass_weighted
        if mass_weighted:
            XY = self.mw_conversion()
            RX = [
                np.tensordot(tf_X, XY, axes=[-1, 0])
                for tf_X in RX
            ]
        return RX
    def get_internals_by_mw_cartesians(self, order=None, mass_weighted=None, coords=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the internals-by-Cartesians Jacobian expansion, converted to be with respect to mass-weighted Cartesians if `mass_weighted` (or `self.mass_weighted` by default) is set.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param mass_weighted: whether to mass-weight the result; defaults to `self.mass_weighted`
        :type mass_weighted: bool | None
        :param coords: alternate coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates
        :type strip_embedding: bool
        :return: the (optionally mass-weighted) internals-by-Cartesians Jacobian tensors
        :rtype: list[np.ndarray]
        """
        XR = self.embedding.get_internals_by_cartesians(
                order=order,
                coords=coords,
                strip_embedding=strip_embedding
            )
        if mass_weighted is None:
            mass_weighted = self.mass_weighted
        if mass_weighted:
            YX = self.mw_inverse()
            YR = []
            if coords is not None:
                pad_dim = coords.ndim - 2
            else:
                pad_dim = 0
            for X_tf in XR:
                for d in range(X_tf.ndim-(1 + pad_dim)):
                    X_tf = np.tensordot(YX, X_tf, axes=[1, d+pad_dim])
                    for j in range(pad_dim):
                        X_tf = np.moveaxis(X_tf, 1+j, 0)
                YR.append(X_tf)
        else:
            YR = XR
        return YR

    def get_internals_by_cartesians(self, order=None, coords=None, strip_embedding=True):
        """
        expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians

        :param order:
        :param strip_embedding:
        :return:
        """
        if self.embedding.internals is None:
            if self.modes is None:
                raise NotImplementedError("not sure what's most consistent for just...plain Cartesians")
            return [self.modes.modes_by_coords]
        else:
            if self.modes is None or self.modes.is_cartesian:
                if self.mass_weighted:
                    YR = self.get_internals_by_mw_cartesians(
                        order=order,
                        coords=coords,
                        strip_embedding=strip_embedding
                    )
                else:
                    YR = self.embedding.get_internals_by_cartesians(
                        order=order,
                        coords=coords,
                        strip_embedding=strip_embedding
                    )
                if self.modes is not None:
                    YQ = self.modes.modes_by_coords
                    if self.mass_weighted:
                        RY = self.get_mw_cartesians_by_internals(
                            order=1,
                            strip_embedding=strip_embedding
                        )[0]
                    else:
                        RY = self.embedding.get_cartesians_by_internals(
                            order=1,
                            strip_embedding=strip_embedding
                        )[0]
                    RQ = np.tensordot(RY, YQ, axes=[-1, -2])
                    padding = RQ.ndim-2
                    for j in range(padding):
                        RQ = np.moveaxis(RQ, padding+j, 0)
                    YR = nput.tensor_reexpand(YR, [RQ], order=order)
                return YR
            else:
                tens = self.embedding.get_internals_by_cartesians(order, coords=coords)
                return nput.tensor_reexpand(
                    tens,
                    [self.modes.modes_by_coords],
                    order=order if order is not None else len(tens)
                )

    def get_cartesians_by_internals(self, order=None, coords=None, strip_embedding=True):
        """
        expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians

        :param order:
        :param strip_embedding:
        :return:
        """
        if self.embedding.internals is None:
            if self.modes is None:
                raise NotImplementedError("not sure what's most consistent for just...plain Cartesians")
            return [self.modes.coords_by_modes]
        else:
            if self.modes is None or self.modes.is_cartesian:
                if self.modes is not None: strip_embedding = True
                if self.mass_weighted:
                    RY = self.get_mw_cartesians_by_internals(
                        order=order,
                        strip_embedding=strip_embedding,
                        coords=coords
                    )
                else:
                    RY = self.embedding.get_cartesians_by_internals(
                        order=order,
                        strip_embedding=strip_embedding,
                        coords=coords
                    )
                if self.modes is not None:
                    QY = self.modes.coords_by_modes
                    if self.mass_weighted:
                        YR = self.get_internals_by_mw_cartesians(
                            order=1,
                            strip_embedding=strip_embedding,
                            coords=coords
                        )[0]
                    else:
                        YR = self.embedding.get_internals_by_cartesians(
                            order=1,
                            strip_embedding=strip_embedding,
                            coords=coords
                        )[0]
                    QR = np.tensordot(QY, YR, axes=[-1, -2])
                    padding = QR.ndim-2
                    for j in range(padding):
                        QR = np.moveaxis(QR, padding+j, 0)
                    RY = nput.tensor_reexpand([QR], RY, order=order)
                return RY
            else:
                tens = self.embedding.get_cartesians_by_internals(order, coords=coords),
                return nput.tensor_reexpand(
                    [self.modes.coords_by_modes],
                    tens,
                    order=order if order is not None else len(tens)
                )

    def get_inertia_tensor_expansion(self, order=None, strip_embedding=True):
        """
        **LLM Docstring**

        Compute the Taylor expansion of the inertia tensor in terms of this embedding's coordinates (internal coordinates or normal modes), by re-expanding the inertial-frame derivatives through the Cartesians-by-internals Jacobian.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the underlying Jacobian
        :type strip_embedding: bool
        :return: `[I0] + [derivative terms...]`, the inertia tensor and its derivatives
        :rtype: list[np.ndarray]
        """
        YI0 = self.embedding.inertial_frame_derivatives()
        QY = self.get_cartesians_by_internals(order=order, strip_embedding=strip_embedding)
        return [self.embedding.inertia_tensor] + nput.tensor_reexpand(QY, YI0, order=order)

    def get_inertial_frame(self):
        """
        **LLM Docstring**

        The molecule's inertial (principal-axis) frame, delegated to the underlying embedding.

        :return: the inertial frame, as returned by `MolecularEmbedding.inertial_frame`
        :rtype: tuple
        """
        return self.embedding.inertial_frame

    def get_modes_by_coords(self, mass_weighted=None, frequency_scaled=None):
        """
        **LLM Docstring**

        The modes-by-coordinates transformation matrix, optionally adjusting the mass-weighting/frequency-scaling convention of the modes first, and (if the modes are Cartesian but expressed relative to an internal-coordinate embedding) re-expressing them in terms of internal coordinates.

        :param mass_weighted: `True`/`False` to force mass-weighting on/off before extracting the matrix; `None` to leave the modes' current convention
        :type mass_weighted: bool | None
        :param frequency_scaled: `True`/`False` to adjust frequency scaling before extracting the matrix (both branches currently call `remove_frequency_scaling`); `None` to leave it unchanged
        :type frequency_scaled: bool | None
        :return: the modes-by-coordinates matrix, or `None` if no modes are set
        :rtype: np.ndarray | None
        """
        if self.modes is None:
            return None
        else:
            clean_modes = self.modes
            if mass_weighted is True:
                clean_modes = clean_modes.make_mass_weighted()
            elif mass_weighted is False:
                clean_modes = clean_modes.remove_mass_weighting()
            if frequency_scaled is True:
                clean_modes = clean_modes.remove_frequency_scaling()
            elif frequency_scaled is False:
                clean_modes = clean_modes.remove_frequency_scaling()
            if self.modes.is_cartesian and self.embedding.internals is not None:
                exp = self.get_mw_cartesians_by_internals(1, mass_weighted=mass_weighted)[0]
                return exp @ clean_modes.modes_by_coords
            else:
                return clean_modes.modes_by_coords
    def get_coords_by_modes(self, mass_weighted=None, frequency_scaled=None):
        """
        **LLM Docstring**

        The coordinates-by-modes transformation matrix, optionally adjusting the mass-weighting/frequency-scaling convention of the modes first, and (if the modes are Cartesian but expressed relative to an internal-coordinate embedding) re-expressing them in terms of internal coordinates.

        :param mass_weighted: `True`/`False` to force mass-weighting on/off before extracting the matrix; `None` to leave the modes' current convention
        :type mass_weighted: bool | None
        :param frequency_scaled: `True`/`False` to adjust frequency scaling before extracting the matrix (both branches currently call `remove_frequency_scaling`); `None` to leave it unchanged
        :type frequency_scaled: bool | None
        :return: the coordinates-by-modes matrix, or `None` if no modes are set
        :rtype: np.ndarray | None
        """
        if self.modes is None:
            return None
        else:
            clean_modes = self.modes
            if mass_weighted is True:
                clean_modes = clean_modes.make_mass_weighted()
            elif mass_weighted is False:
                clean_modes = clean_modes.remove_mass_weighting()
            if frequency_scaled is True:
                clean_modes = clean_modes.remove_frequency_scaling()
            elif frequency_scaled is False:
                clean_modes = clean_modes.remove_frequency_scaling()
            if self.modes.is_cartesian and self.embedding.internals is not None:
                exp = self.get_internals_by_mw_cartesians(1, mass_weighted=mass_weighted)[0]
                return clean_modes.coords_by_modes @ exp
            else:
                return clean_modes.coords_by_modes

def _get_best_axes(first_pos, axes):
    """
    Determine the best pair of inertial axes so that we don't get large-scale breakdowns from the choice of embedding

    :param first_pos:
    :type first_pos:
    :param axes:
    :type axes:
    :return:
    :rtype:
    """

    if axes.ndim > 2:
        axes = axes[..., (0, 1), :]
        ax_choice = (0, 1)
        ax_names = ["A", "B"]
    else:
        fp_norm = np.linalg.norm(first_pos, axis=-1)
        if fp_norm > 1.0e-10:  # not chilling at the origin...
            first_pos = first_pos / fp_norm
            # check if it lies along an axis or is perpendicular to an axis
            a_proj = np.dot(first_pos, axes[0])
            b_proj = np.dot(first_pos, axes[1])
            c_proj = np.dot(first_pos, axes[2])
            if np.abs(b_proj) < .05: # lies in the A/C plane
                if np.abs(a_proj) > .95:
                    ax_choice = (1, 2)
                    ax_names = ["B", "C"]
                else:
                    ax_choice = (0, 1)
                    ax_names = ["A", "B"]
            elif np.abs(c_proj) < .05: # lies in the A/B plane
                if np.abs(a_proj) > .95:
                    ax_choice = (1, 2)
                    ax_names = ["B", "C"]
                else:
                    ax_choice = (0, 2)
                    ax_names = ["A", "C"]
            elif np.abs(a_proj) < .05:  # lies in the B/C plane
                if np.abs(b_proj) > .95:
                    ax_choice = (0, 2)
                    ax_names = ["A", "C"]
                else:
                    ax_choice = (0, 1)
                    ax_names = ["A", "B"]
            else: # not in any of the planes so no issues
                ax_choice = (0, 1)
                ax_names = ["A", "B"]

        else:
            ax_choice = (0, 1)
            ax_names = ["A", "B"]
        axes = axes[ax_choice,]
    return axes, ax_names, ax_choice

class MolecularGenericInternalCoordinateSystem(GenericInternalCoordinateSystem):
    """
    Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding
    """
    name = "MolecularGenericInternals"
    class PassThroughRedundantGenerator:
        def __init__(self, redundant_transformation, redundant_inverse=None,
                     masses=None,
                     untransformed_coordinates=None,
                     relocalize=False,
                     **opts
                     ):
            """
            **LLM Docstring**

            Store an already-known redundant-coordinate transformation so it can be handed back unchanged (rather than computed) when `get_redundant_transformation` is called.

            :param redundant_transformation: the redundant-to-non-redundant transformation matrix to pass through
            :type redundant_transformation: np.ndarray
            :param redundant_inverse: the corresponding inverse transformation, if available
            :type redundant_inverse: np.ndarray | None
            :param masses: the atomic masses (stored as a 1-tuple due to a trailing comma in the assignment)
            :type masses: np.ndarray | None
            :param untransformed_coordinates: coordinates to leave untransformed
            :type untransformed_coordinates: object | None
            :param relocalize: whether the redundant coordinates should be relocalized
            :type relocalize: bool
            :param opts: extra options, stored but not otherwise used
            :type opts: dict
            :return: None
            :rtype: None
            """
            self.tf = redundant_transformation
            self.inv = redundant_inverse
            self.masses = masses,
            self.untransformed_coordinates = untransformed_coordinates
            self.relocalize = relocalize
            self.opts = opts
        def get_redundant_transformation(self, base_expansions, **ignored):
            """
            **LLM Docstring**

            Return the stored redundant transformation unchanged, along with the base derivative expansions re-expanded through it.

            :param base_expansions: the derivative expansions to re-express through the redundant transformation
            :type base_expansions: list[np.ndarray]
            :param ignored: any other arguments, accepted but not used
            :type ignored: dict
            :return: `(self.tf, reexpanded_expansions)`
            :rtype: tuple[np.ndarray, list[np.ndarray]]
            """
            return self.tf, nput.tensor_reexpand(base_expansions,
                                                 [self.tf],
                                                 order=len(base_expansions)
                                                 )
    def __init__(self, masses, coords, /, specs, converter_options=None,
                 redundant=False,
                 relocalize=False,
                 untransformed_coordinates=None,
                 redundant_transformation=None,
                 redundant_inverse=None,
                 angle_ordering='ijk',
                 **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        self.masses = masses
        self.coords = coords
        if converter_options is None:
            converter_options = dict(opts, reference_coordinates=coords, masses=masses, angle_ordering=angle_ordering, specs=specs)
            opts = {}
        if redundant and 'redundant_generator' not in converter_options:
            if redundant_transformation is not None:
                redundant_generator = self.PassThroughRedundantGenerator(
                    redundant_transformation,
                    redundant_inverse=redundant_inverse,
                    masses=self.masses,
                    untransformed_coordinates=untransformed_coordinates,
                    angle_ordering=angle_ordering,
                    relocalize=relocalize
                )
            else:
                redundant_generator = RedundantCoordinateGenerator(
                    specs,
                    masses=self.masses,
                    untransformed_coordinates=untransformed_coordinates,
                    angle_ordering=angle_ordering,
                    relocalize=relocalize
                )
            converter_options['redundant_generator'] = converter_options.get('redundant_generator', redundant_generator)
        nats = None#len(specs)
        super().__init__(converter_options=converter_options, dimension=(nats,), coordinate_shape=(nats,), opts=opts)

class MolecularZMatrixCoordinateSystem(ZMatrixCoordinateSystem):
    """
    Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding
    """
    name = "MolecularZMatrix"
    embedding_coords = [0, 1, 2, 4, 5, 8]
    def __init__(self, masses, coords, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        self.masses = masses
        self.coords = coords
        self.dummy_positions = [i for i, m in enumerate(masses) if m < 0]
        self.com = StructuralProperties.get_prop_center_of_mass(coords, masses)
        self.inertial_axes = StructuralProperties.get_prop_moments_of_inertia(coords, masses)[1]
        if converter_options is None:
            converter_options = {}
        converter_options = dict(converter_options, **opts)
        nats = len(masses)
        super().__init__(converter_options=converter_options, dimension=(nats, 3), coordinate_shape=(nats, 3), opts=opts)
        self.set_embedding()
    @property
    def origins(self):
        """
        **LLM Docstring**

        The Z-matrix embedding's origin points (typically the reference center of mass), from `converter_options['origins']`.

        :return: the origin points
        :rtype: np.ndarray
        """
        return self.converter_options['origins']
    @property
    def axes(self):
        """
        **LLM Docstring**

        The Z-matrix embedding's reference axes (typically two principal axes), from `converter_options['axes']`.

        :return: the reference axes
        :rtype: np.ndarray
        """
        return self.converter_options['axes']

    def get_direct_converter(self, target):
        """
        **LLM Docstring**

        Provide a converter from this molecular Z-matrix system directly to the plain (non-molecular) `ZMatrix` coordinate system, if `target` is one.

        :param target: the coordinate system being converted to
        :type target: object
        :return: a `MolecularZMatrixToRegularZMatrixConverter`, or `None` if `target` isn't a `ZMatrix` system
        :rtype: MolecularZMatrixToRegularZMatrixConverter | None
        """
        if target.name == 'ZMatrix':
            zmat = MolecularZMatrixToRegularZMatrixConverter(self.coords.system)
            return zmat
    def get_inverse_converter(self, target):
        """
        **LLM Docstring**

        Provide a converter from the plain `ZMatrix` coordinate system into this molecular Z-matrix system, if `target` is one.

        :param target: the coordinate system being converted from
        :type target: object
        :return: a `RegularZMatrixToMolecularZMatrixConverter`, or `None` if `target` isn't a `ZMatrix` system
        :rtype: RegularZMatrixToMolecularZMatrixConverter | None
        """
        if target.name == 'ZMatrix':
            zmat = RegularZMatrixToMolecularZMatrixConverter(self.coords.system)
            return zmat

    def pre_convert_to(self, system, opts=None):
        """
        **LLM Docstring**

        Re-establish the embedding options (via `set_embedding`) before delegating to the base class's `pre_convert_to`, preserving the existing atom `ordering` if the caller supplied its own options dict.

        :param system: the coordinate system being converted to
        :type system: object
        :param opts: explicit conversion options to use instead of `self.converter_options`
        :type opts: dict | None
        :return: the resolved conversion options
        :rtype: dict
        """
        self.set_embedding()
        if opts is None:
            opts = self.converter_options
        else:
            opts = opts | {'ordering':self.converter_options['ordering']}
        opts = super().pre_convert_to(system, opts)
        return opts

    def set_embedding(self):
        """
        **LLM Docstring**

        (Re)compute and store this Z-matrix system's embedding options -- the reference origin, reference axes (chosen via `_get_best_axes` to avoid ill-conditioned choices), axis labels, masses, dummy-atom positions, and reference coordinates -- based on the current center of mass and inertial axes, fixing up the Z-matrix `ordering`'s first three rows to reference the embedding's dummy origin/axis points if an ordering is present.

        :return: None
        :rtype: None
        """
        com = self.com
        axes = self.inertial_axes
        converter_options = self.converter_options
        if 'ordering' in converter_options:
            ordering = np.array(converter_options['ordering'], dtype=int).reshape(-1, 4)
            ordering[0, 1] = -3; ordering[0, 2] = -1; ordering[0, 3] = -2
            ordering[1, 2] = -1; ordering[1, 3] = -2
            if len(ordering) > 2:
                ordering[2, 3] = -2
            converter_options['ordering'] = ordering
            first = ordering[0, 0]
        else:
            first = 0

        first_pos = self.coords[..., first, :]
        axes, ax_names, ax_choice = _get_best_axes(first_pos, axes)

        if 'origins' not in converter_options or converter_options['origins'] is None:
            converter_options['origins'] = com
        if 'axes' not in converter_options or converter_options['axes'] is None:
            converter_options['axes'] = axes
        converter_options['axes_labels'] = ax_names
        converter_options['axes_choice'] = ax_choice
        converter_options['masses'] = self.masses
        converter_options['dummy_positions'] = self.dummy_positions
        converter_options['ref_coords'] = self.coords

    def jacobian(self,
                 coords,
                 *args,
                 reembed=None,
                 strip_dummies=None,
                 converter_options=None,
                 **kwargs
                 ):
        """
        **LLM Docstring**

        Compute the Jacobian of this Z-matrix system with respect to Cartesian coordinates, handling batched/multi-frame inputs, optional dummy-atom stripping, and temporarily overriding the `reembed` converter option for the duration of the call.

        :param coords: the Cartesian coordinates to compute the Jacobian at
        :type coords: np.ndarray
        :param args: extra positional arguments forwarded to the base class's `jacobian`
        :type args: tuple
        :param reembed: whether to re-embed (Eckart-align) during the Jacobian calculation; falls back to the converter options, then defaults to `True`
        :type reembed: bool | None
        :param strip_dummies: whether to exclude dummy-atom coordinates from the Jacobian; falls back to the converter options, then defaults to `False`
        :type strip_dummies: bool | None
        :param converter_options: extra converter options merged with `self.converter_options`
        :type converter_options: dict | None
        :param kwargs: extra keyword arguments forwarded to the base class's `jacobian`
        :type kwargs: dict
        :return: the computed Jacobian tensor(s)
        :rtype: list[np.ndarray]
        """

        if converter_options is None:
            converter_options = {}
        merged_convert_options = dict(self.converter_options, **converter_options)
        try:
            remb = merged_convert_options['reembed'] if reembed is None else reembed
        except KeyError:
            remb = None

        try:
            strip_dummies = merged_convert_options['strip_dummies'] if strip_dummies is None else strip_dummies
        except KeyError:
            strip_dummies = False

        if strip_dummies:
            dummies = self.dummy_positions
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(self.masses),
                dummies
            )
        else:
            main_excludes = None

        if coords.ndim > 3:
            og_shape = coords.shape
            coords = coords.reshape((-1,) + coords.shape[-2:])
        else:
            og_shape = (coords.shape[0],) if coords.ndim == 3 else ()


        try:
            self.converter_options['reembed'] = True if remb is None else remb
            jacs = super().jacobian(coords, *args, converter_options=converter_options, **kwargs)
            if isinstance(jacs, np.ndarray):
                jacs = [jacs]
            raw_jacs = []
            if main_excludes is not None:
                main_excludes = np.delete(
                    np.arange(3*len(self.masses)),
                    nput.block_broadcast_indices(main_excludes, 3)
                )
            for j in jacs:
                ## this was to group shapes up in pairs of two since the
                ## finite difference added too much structure
                # skip_dim = coords.ndim - 2
                # if skip_dim > 0:
                #     j = np.moveaxis(j, -3, 0) # skip_dim == 1 by construction so only need 1 move...
                # ext_dim = j.ndim - 2 - skip_dim
                # shp = og_shape + sum(
                #     ((j.shape[i] // 3, 3) for i in range(skip_dim, ext_dim+skip_dim)),
                #     ()
                # ) + j.shape[-2:]
                # j = j.reshape(shp)
                if main_excludes is not None:
                    j = np.take(j, main_excludes, axis=-1)
                raw_jacs.append(j)
            jacs = raw_jacs
            return jacs
        finally:
            if remb is not None:
                self.converter_options['reembed'] = remb

class MolecularCartesianCoordinateSystem(CartesianCoordinateSystem):
    """
    Mirrors the standard Cartesian coordinate system in _almost_ all regards, but forces an embedding
    """
    name= "MolecularCartesians"
    def __init__(self, masses, coords, dummy_positions=None, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        masses = np.asanyarray(masses)
        coords = np.asanyarray(coords)
        self.masses = masses
        self.coords = coords
        self.dummy_positions = [i for i,m in enumerate(masses) if m < 0] if dummy_positions is None else dummy_positions
        self.com = StructuralProperties.get_prop_center_of_mass(coords, masses)
        self.axes = StructuralProperties.get_prop_moments_of_inertia(coords, masses)[1]
        nats = len(masses)
        if converter_options is None or len(converter_options) == 0:
            converter_options = opts
            opts = {}
        super().__init__(converter_options=converter_options, dimension=(nats, 3), coordinate_shape=(nats, 3), **opts)
    def to_state(self, serializer=None):
        """
        **LLM Docstring**

        Serialize this coordinate system's state, adding the masses, coordinates, and dummy-atom positions on top of whatever the base class's `to_state` produces.

        :param serializer: the serializer to use, forwarded to the base class
        :type serializer: object | None
        :return: the serialized state dict
        :rtype: dict
        """
        base_data = super().to_state(serializer=serializer)
        base_data['masses'] = self.masses
        base_data['coords'] = self.coords
        base_data['dummy_positions'] = self.dummy_positions
        return base_data
    @classmethod
    def from_state(cls, data, serializer=None):
        """
        **LLM Docstring**

        Reconstruct a `MolecularCartesianCoordinateSystem` from a previously serialized state dict.

        :param data: the serialized state, as produced by `to_state`
        :type data: dict
        :param serializer: the serializer to use, accepted for interface consistency but not used in this method's body
        :type serializer: object | None
        :return: the reconstructed coordinate system
        :rtype: MolecularCartesianCoordinateSystem
        """
        # dim = data.pop('dimension', None)
        # coordinate_shape = data.pop('coordinate_shape', None)
        return cls(
            data['masses'],
            data['coords'],
            dummy_positions=data['dummy_positions'],
            converter_options=data['converter_options']
        )
    def pre_convert_to(self, system, opts=None):
        """
        **LLM Docstring**

        Ensure the masses are up to date in `converter_options`, and, if converting to a Z-matrix-family system, re-establish the embedding options (via `set_embedding`) before delegating to the base class's `pre_convert_to`.

        :param system: the coordinate system being converted to
        :type system: object
        :param opts: explicit conversion options to use instead of `self.converter_options`
        :type opts: dict | None
        :return: the resolved conversion options
        :rtype: dict
        """
        self.converter_options['masses'] = self.masses
        if 'ZMatrix' in system.name:
            self.set_embedding()
            if opts is None:
                opts = self.converter_options
            else:
                opts = self.converter_options | opts
            opts = super().pre_convert_to(system, opts)
        return opts

    def set_embedding(self):
        """
        Sets up the embedding options...
        :return:
        :rtype:
        """
        com = self.com
        axes = self.axes
        converter_options = self.converter_options
        if 'ordering' in converter_options:
            ordering = np.array(converter_options['ordering'], dtype=int).reshape(-1, 4)
            ordering[0, 1] = -3; ordering[0, 2] = -1; ordering[0, 3] = -2
            ordering[1, 2] = -1; ordering[1, 3] = -2
            if len(ordering) > 2:
                ordering[2, 3] = -2
            converter_options['ordering'] = ordering
            first = ordering[0, 0]
        else:
            first = 0

        first_pos = self.coords[first]
        axes, ax_names, ax_choice = _get_best_axes(first_pos, axes)

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['axes_labels'] = ax_names
        converter_options['axes_choice'] = ax_choice
        converter_options['masses'] = self.masses
        converter_options['ref_coords'] = self.coords
        converter_options['dummy_positions'] = self.dummy_positions

    def jacobian(self,
                 coords,
                 system,
                 order=None,
                 strip_dummies=None,
                 converter_options=None,
                 analytic_deriv_order=None,
                 **kwargs
                 ):
        """
        **LLM Docstring**

        Compute the Jacobian of these Cartesian coordinates with respect to `system`, resolving the analytic-derivative order (defaulting to purely numerical for Z-matrix targets) and optionally excluding dummy-atom coordinates.

        :param coords: the coordinates to compute the Jacobian at
        :type coords: np.ndarray
        :param system: the target coordinate system
        :type system: object
        :param order: the derivative order(s) to compute
        :type order: int | list[int] | None
        :param strip_dummies: whether to exclude dummy-atom coordinates; falls back to the converter options, then defaults to `False`
        :type strip_dummies: bool | None
        :param converter_options: extra converter options merged with `self.converter_options`
        :type converter_options: dict | None
        :param analytic_deriv_order: the order up to which to compute the Jacobian analytically rather than numerically; falls back to the converter options, then defaults based on whether `system` is a Z-matrix
        :type analytic_deriv_order: int | None
        :param kwargs: extra keyword arguments forwarded to the base class's `jacobian`
        :type kwargs: dict
        :return: the computed Jacobian tensor(s)
        :rtype: list[np.ndarray]
        """

        zmat_conv = 'ZMatrix' in system.name

        if converter_options is None:
            converter_options = {}
        merged_convert_options = dict(self.converter_options, **converter_options)
        try:
            strip_dummies = merged_convert_options['strip_dummies'] if strip_dummies is None else strip_dummies
        except KeyError:
            strip_dummies = False

        try:
            analytic_deriv_order = (
                merged_convert_options['analytic_deriv_order']
                    if analytic_deriv_order is None else
                analytic_deriv_order
            )
        except KeyError:
            if zmat_conv:
                analytic_deriv_order = 0
            else:
                # if return_derivs is None:
                #     return_derivs = order
                analytic_deriv_order = None #return_derivs if nput.is_numeric(return_derivs) else max(return_derivs)

        if strip_dummies:
            dummies = self.dummy_positions
            if len(dummies) == 0:
                dummies = None
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(len(self.masses)),
                dummies
            )
        else:
            main_excludes = None
        jacs = super().jacobian(coords, system,
                                order=order,
                                analytic_deriv_order=analytic_deriv_order,
                                converter_options=converter_options,
                                **kwargs)
        if isinstance(jacs, np.ndarray):
            jacs = [jacs]

        if analytic_deriv_order < 0:
            analytic_deriv_order = len(jacs)

        if zmat_conv:
            raw_jacs = jacs
            # for n,j in enumerate(jacs): # this expects a full filling of the jacobians which maybe I need to not expect...
            #     print("???", j.shape)
            #     baseline = 2*analytic_deriv_order + len(coords.shape)
            #     ext_dim = j.ndim - baseline
            #     shp = sum(
            #         ((j.shape[i] // 3, 3) for i in range(ext_dim)),
            #         ()
            #     ) + j.shape[-baseline:]
            #     j = j.reshape(shp)
            #     if dummies is not None:
            #         for i in range(ext_dim):
            #             j = np.take(j, main_excludes, axis=2*i)
            #         for i in range(analytic_deriv_order):
            #             j = np.take(j, main_excludes, axis=-2*(i+2))
            #
            #     base_shape = coords.shape[:-2]
            #     if base_shape != j.shape[:len(base_shape)]:
            #         j = np.moveaxis(j, -3, 0)
            #
            #     raw_jacs.append(j)
            #     print("???...", j.shape)
            jacs = raw_jacs
        return jacs

class MolecularCartesianToZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """
    def __init__(self, cart_system, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(cart_system, zmat_system)` type pair this converter handles.

        :param cart_system: the source Cartesian coordinate system
        :type cart_system: object
        :param zmat_system: the target Z-matrix coordinate system
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (cart_system, zmat_system)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(cart_system, zmat_system)`
        :rtype: tuple
        """
        return self._types
    def convert(self, coords, *,
                     masses,
                     dummy_positions, origins=None, axes=None, ordering=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        :param coords:
        :type coords: CoordinateSet
        :param molecule:
        :type molecule:
        :param origins:
        :type origins:
        :param axes:
        :type axes:
        :param ordering:
        :type ordering:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """

        if origins is not None:
            if np.linalg.norm(coords[0] - origins) < nput.Options.norm_zero_threshold:
                origins = (
                        origins
                        + nput.vec_normalize(np.random.uniform(size=origins.shape)) * 2 * nput.Options.norm_zero_threshold
                )

        zmcs, opts = self.convert_many(np.array([coords]),
                                       masses=masses,
                                       dummy_positions=dummy_positions,
                                       origins=origins, axes=axes, ordering=ordering, **kwargs)
        zmcs = zmcs[0]

        if 'derivs' in opts:
            derivs = opts['derivs']
            reshaped_derivs = [None] * len(derivs)
            for i, v in enumerate(derivs):
                reshaped_derivs[i] = v[0]
            opts['derivs'] = reshaped_derivs

        return zmcs, opts

    base_cartesian_type = CartesianCoordinates3D
    base_internal_type = ZMatrixCoordinates
    def convert_many(self,
                     coords,*,
                     masses,
                     dummy_positions,
                     origins=None, axes=None,
                     ordering=None,
                     strip_embedding=True,
                     strip_dummies=False,
                     return_derivs=None,
                     **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding

        :param coords: coordinates in Cartesians to convert
        :type coords: np.ndarray
        :param molecule:
        :type molecule: MolecularEmbedding
        :param origins: the origin for each individual structure
        :type origins: np.ndarray
        :param axes: the axes for each structure
        :type axes: np.ndarray
        :param ordering: the Z-matrix ordering spec
        :type ordering:
        :param strip_embedding: whether to strip the embedding coordinates
        :type strip_embedding:
        :param strip_dummies: whether to strip all dummy coordinates
        :type strip_dummies:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """
        # return_derivs=False

        n_sys = coords.shape[0]
        n_coords = coords.shape[1]
        n_atoms = len(masses)

        # we add three dummy atoms at the origins and along the axes before doing the conversion
        if origins.ndim == 1:
            origins = np.broadcast_to(origins[np.newaxis, np.newaxis], (n_sys, 1, 3))
        elif origins.ndim == 2:
            origins = origins[:, np.newaxis, :]
        elif origins.ndim > 2:
            origins = np.reshape(origins, (-1, 1, 3))
        if axes.ndim == 2:
            axes = np.broadcast_to(axes[np.newaxis], (n_sys, 2, 3))
        elif axes.ndim > 2:
            axes = np.reshape(axes, (-1, 2, 3))
        if origins.shape[0] != n_sys:
            if n_sys % origins.shape[0] != 0:
                raise ValueError("inconsistent shapes; origins shape {} but coords shape {}".format(
                    origins.shape,
                    coords.shape
                ))
            num_coords = n_sys // origins.shape[0]
            origins = np.broadcast_to(origins[:, np.newaxis, :, :], (origins.shape[0], num_coords) + origins.shape[1:])
            origins = origins.reshape((n_sys,) + origins.shape[2:])
        if axes.shape[0] != n_sys:
            if n_sys % axes.shape[0] != 0:
                raise ValueError("inconsistent shapes; axes shape {} but coords shape {}".format(
                    axes.shape,
                    coords.shape
                ))
            num_coords = n_sys // axes.shape[0]
            axes = np.broadcast_to(axes[:, np.newaxis, :, :], (axes.shape[0], num_coords) + axes.shape[1:])
            axes = axes.reshape((n_sys,) + axes.shape[2:])
        coords = np.concatenate([origins, origins+axes, coords], axis=1)
        if ordering is not None:
            ordering = np.array(ordering, dtype=int).reshape(-1, 4)
            ordering[0, 1] = -3; ordering[0, 2] = -1; ordering[0, 3] = -2
            ordering[1, 2] = -1; ordering[1, 3] = -2
            if len(ordering) > 2:
                ordering[2, 3] = -2
            ordering = ordering + 3
            ordering = np.concatenate([ [[0, -1, -2, -3], [1, 0, -1, -2], [2, 0, 1, -1]], ordering])
        res = CoordinateSet(coords, self.base_cartesian_type).convert(self.base_internal_type,
                                                                      ordering=ordering,
                                                                      origins=origins,
                                                                      axes=axes,
                                                                      return_derivs=return_derivs,
                                                                      masses=masses,
                                                                      **kwargs
                                                                      )

        if isinstance(res, tuple):
            zmcs, opts = res
        else:
            zmcs = res
            opts=res.converter_options
        opts['ordering'] = opts['ordering'][3:] - 3
        if masses is not None:
            opts['masses'] = masses
        # zmcs = zmcs[:, 2:]
        if strip_dummies:
            dummies = [0, 1, 2] + [x+3 for x in dummy_positions] # add on axes
        elif strip_embedding:
            dummies = [0, 1, 2]
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(n_atoms + 3),
                dummies
            )
            sub_excludes = main_excludes - 1 # drop one fewer terms to drop I think...
            if 'derivs' in opts:
                derivs = opts['derivs']
                reshaped_derivs = [None] * len(derivs)
                deriv_excludes = np.arange(3, n_atoms + 3)
                for i, v in enumerate(derivs):
                    # drop all terms relating to the embedding of the embedding
                    start_dim = v.ndim - 2*(i+2)
                    for j in range(start_dim, v.ndim-2, 2):
                        v = np.take(v, deriv_excludes, axis=j)
                    v = np.take(v, sub_excludes, axis=-2)
                    reshaped_derivs[i] = v

                opts['derivs'] = reshaped_derivs

            zmcs = zmcs[..., sub_excludes, :]
        return zmcs, opts
# MolecularCartesianToZMatrixConverter = MolecularCartesianToZMatrixConverter()
# MolecularCartesianToZMatrixConverter.register(CoordinateSystemConverters)
class MolecularCartesianToRegularCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, cart_system, **opts):
        """
        **LLM Docstring**

        Store the `(CartesianCoordinates3D, cart_system)` type pair this converter handles.

        :param cart_system: the molecular Cartesian coordinate system being converted from
        :type cart_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (CartesianCoordinates3D, cart_system)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(CartesianCoordinates3D, cart_system)`
        :rtype: tuple
        """
        return self._types
    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a molecular Cartesian system and a plain Cartesian system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        return coords, kwargs

class RegularCartesianToMolecularCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, cart_system, **opts):
        """
        **LLM Docstring**

        Store the `(cart_system, CartesianCoordinates3D)` type pair this converter handles.

        :param cart_system: the molecular Cartesian coordinate system being converted to
        :type cart_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (cart_system, CartesianCoordinates3D)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(cart_system, CartesianCoordinates3D)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a plain Cartesian system and a molecular Cartesian system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        return coords, kwargs

class MolecularZMatrixToCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """
    def __init__(self, zmat_system, cart_system, **opts):
        """
        **LLM Docstring**

        Store the `(zmat_system, cart_system)` type pair this converter handles.

        :param zmat_system: the source Z-matrix coordinate system
        :type zmat_system: object
        :param cart_system: the target Cartesian coordinate system
        :type cart_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (zmat_system, cart_system)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(zmat_system, cart_system)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Convert a single frame of Z-matrix coordinates to Cartesians by delegating to `convert_many` on a length-1 batch and unwrapping the result (including the per-frame `'derivs'` entries, if present).

        :param coords: the single-frame Z-matrix coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments forwarded to `convert_many`
        :type kw: dict
        :return: `(cartesian_coords, opts)` for the single frame
        :rtype: tuple
        """
        total_points, opts = self.convert_many(coords[np.newaxis], **kw)
        if 'derivs' in opts:
            opts['derivs'] = [o[0] for o in opts['derivs']]
        return total_points[0], opts

    base_cartesian_type = CartesianCoordinates3D
    base_internal_type = ZMatrixCoordinates
    def convert_many(self, coords, *,
                     masses, dummy_positions, ref_coords,
                     origins=None, axes=None, ordering=None,
                     reembed=False, axes_choice=None, return_derivs=None,
                     strip_dummies=False,
                     strip_embedding=True,
                     planar_ref_tolerance=None,
                     embedding_masses=None,
                     spec=None,
                     **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, attempting to preserve the embedding
        """
        # from .Molecule import Molecule
        # return_derivs = False

        n_sys = coords.shape[0]
        n_coords = coords.shape[1]
        n_atoms = len(masses)
        if embedding_masses is None:
            embedding_masses = masses
        if n_coords != n_atoms + 2:
            # means we already added the embedding
            if n_coords != n_atoms:
                raise ValueError('Embedding unclear when num_coords ({}) < num_atoms ({})'.format(
                    n_coords,
                    n_atoms
                ))

            # TODO:
            # these need broadcasting when
            # doing jacobian calcs but this is currently pretty
            # hacky
            # print(axes.shape, n_sys)
            if axes.ndim == 2:
                axes = axes[np.newaxis]
            if axes.shape[0] > 1 and axes.shape[0] < n_sys:
                # print((n_sys // axes.shape[0], axes.shape[0]) + axes.shape)
                axes = np.reshape(
                    np.broadcast_to(
                        axes[np.newaxis],
                        (n_sys // axes.shape[0],) + axes.shape
                    ),
                    (n_sys, ) + axes.shape[1:]
                )

            if origins.ndim == 1:
                origins = origins[np.newaxis]
            # print(axes.shape, n_sys)
            if origins.shape[0] > 1 and origins.shape[0] < n_sys:
                # print((n_sys // axes.shape[0], axes.shape[0]) + axes.shape)
                origins = np.reshape(
                    np.broadcast_to(
                        origins[np.newaxis],
                        (n_sys // origins.shape[0],) + origins.shape
                    ),
                    (n_sys, ) + origins.shape[1:]
                )

            x_ax = axes[..., 0, :]
            y_ax = axes[..., 1, :]
            extra_norms0 = nput.vec_norms(x_ax)
            extra_norms1 = nput.vec_norms(y_ax)
            extra_angles, _ = nput.vec_angles(x_ax, y_ax)

            extra_coords = np.zeros((n_sys, 2, 3))
            extra_coords[..., 0, 0] = extra_norms0
            extra_coords[..., 1, 0] = extra_norms1
            extra_coords[..., 1, 1] = extra_angles

            coords = np.concatenate([extra_coords, coords], axis=-2)
            if ordering is not None:
                ordering = np.array(ordering, dtype=int).reshape(-1, 4)
                ordering = ordering + 3
                ordering = np.concatenate([ [[0, -1, -2, -3], [1, 0, -1, -2], [2, 0, 1, -1]], ordering])

            # print(ordering[15:19])
            # raise Exception(
            #     [i for i,c in enumerate(McUtils.Coordinerds.extract_zmatrix_internals(ordering))
            #      if c is None]
            # )
            # if spec is None:
            if return_derivs:
                internals = McUtils.Coordinerds.extract_zmatrix_internals(ordering)
                if not strip_embedding:
                    spec = McUtils.Coordinerds.InternalSpec(
                        internals[3:],
                        ungraphed_internals=[
                            internals[3+i]
                            for i in MolecularZMatrixCoordinateSystem.embedding_coords
                        ],
                        masses=np.concatenate([[1e6]*3, embedding_masses])
                    )
                else:
                    spec = McUtils.Coordinerds.InternalSpec(
                        [
                            internals[3 + i]
                            for i in np.setdiff1d(
                                np.arange(len(internals) - 3),
                                MolecularZMatrixCoordinateSystem.embedding_coords
                            )
                        ],
                        masses=np.concatenate([[1e6]*3, embedding_masses])
                    )
            # else:
            #     internals = McUtils.Coordinerds.extract_zmatrix_internals(ordering)
            #     reindexing = np.arange(len(masses)) + 3
            #     spec = McUtils.Coordinerds.InternalSpec(
            #         internals[:6] + [
            #             c.reindex(reindexing)
            #             for c in spec.coords
            #         ],
            #         ungraphed_internals=internals[:6]
            #     )

        refuse_derivs = reembed and coords.squeeze().ndim != 2
        res = CoordinateSet(coords, self.base_internal_type).convert(self.base_cartesian_type,
                                                                     ordering=ordering,
                                                                     origins=origins,
                                                                     axes=axes,
                                                                     return_derivs=return_derivs and not refuse_derivs,
                                                                     masses=embedding_masses,
                                                                     spec=spec,
                                                                     **kwargs)

        if isinstance(res, tuple):
            carts, opts = res
        else:
            carts = res
            opts = res.converter_options

        if reembed:
            embed_carts = carts[..., 3:, :]
            reembed = not (
                    carts.squeeze().ndim == 2 and
                    np.allclose(ref_coords, embed_carts, atol=1.0e-5)
            ) # agree to like a ten thousandth of an angstrom
            if reembed:
                # if not return_derivs:
                embed_carts = StructuralProperties.get_eckart_embedded_coords(
                    masses,
                    ref_coords,
                    embed_carts,
                    planar_ref_tolerance=planar_ref_tolerance
                )
                carts = np.concatenate([
                    carts[..., :3, :],
                    embed_carts
                    ],
                axis=-2
                )
                # else:
                #     print(coords.shape, masses.shape)
                #     inert_coords, coord_coms, coord_axes = StructuralProperties.get_prop_principle_axis_data(
                #         coords,
                #         masses
                #     ).principle_axis_data
                #     if axes_choice is None:
                #         axes_choice = (0, 1)
                #     guh = self.convert_many(coords,
                #                             origins=coord_coms,
                #                             axes=coord_axes[:, axes_choice],
                #                             masses=masses,
                #                             ref_coords=ref_coords,
                #                             dummy_positions=dummy_positions,
                #                             reembed=False,
                #                             ordering=ordering,
                #                             return_derivs=return_derivs,
                #                             axes_choice=axes_choice,
                #                             **kwargs
                #                             )
                #     return guh

        opts['origins'] = origins
        opts['axes'] = axes
        if ordering is not None:
            opts['ordering'] = ordering[3:] - 3
        if strip_dummies:
            dummies = [0, 1, 2] + [x + 3 for x in dummy_positions]
        else:
            dummies = [0, 1, 2]
        # expansion_terms = kwargs.get('use_direct_expansions')
        main_excludes = np.setdiff1d(
            np.arange(len(masses) + 3),
            dummies
        )
        if 'derivs' in opts:
            sub_excludes = nput.block_broadcast_indices(main_excludes, 3)
            opts['derivs'] = [v[..., sub_excludes] for v in opts['derivs']]
        carts = carts[..., main_excludes, :]

        return carts, opts

class MolecularZMatrixToRegularZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """
    def __init__(self, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(zmat_system, ZMatrixCoordinateSystem)` type pair this converter handles.

        :param zmat_system: the molecular Z-matrix coordinate system being converted from
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (zmat_system, ZMatrixCoordinateSystem)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(zmat_system, ZMatrixCoordinateSystem)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a molecular Z-matrix system and a plain Z-matrix system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a molecular Z-matrix system and a plain Z-matrix system share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        return coords, kwargs

class RegularZMatrixToMolecularZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(ZMatrixCoordinateSystem, zmat_system)` type pair this converter handles.

        :param zmat_system: the molecular Z-matrix coordinate system being converted to
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (ZMatrixCoordinateSystem, zmat_system)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(ZMatrixCoordinateSystem, zmat_system)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a plain Z-matrix system and a molecular Z-matrix system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a plain Z-matrix system and a molecular Z-matrix system share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        return coords, kwargs

class MolecularCartesianToGICConverter(CartesianToGICSystemConverter):
    """
    ...
    """
    def __init__(self, cart_system, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(cart_system, zmat_system)` type pair this converter handles, and precompute which internal-coordinate specs are periodic (angle-like specs with more than 2 atoms) for later use in handling periodic wraparound.

        :param cart_system: the source Cartesian coordinate system
        :type cart_system: object
        :param zmat_system: the target generic-internal coordinate system
        :type zmat_system: object
        :param opts: extra options forwarded to the base converter
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (cart_system, zmat_system)
        self._periodics = np.where([len(x) > 2 for x in zmat_system.converter_options['specs']])[0]
        super().__init__(**opts)

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(cart_system, zmat_system)`
        :rtype: tuple
        """
        return self._types

    def convert_many(self, coords, *,
                     order=0,
                     return_derivs=None,
                     redundant_generator:RedundantCoordinateGenerator=None,
                     reference_internals=None,
                     redundant_transformation=None,
                     handle_periodicity=True,
                     **kw):
        """
        We'll implement this by having the ordering arg wrap around in coords?
        """
        internals, opts = super().convert_many(coords, order=order, return_derivs=return_derivs, **kw)
        if redundant_generator is not None or redundant_transformation is not None:
            #TODO: merge the redundant handling into the lower order stuff...
            if redundant_transformation is not None:
                red_tf = redundant_transformation
                red_exp = nput.tensor_reexpand(opts['derivs'], [red_tf], order=len(opts['derivs']))
            else:
                red_tf, red_exp = redundant_generator.get_redundant_transformation(
                    opts['derivs'],
                    untransformed_coordinates=redundant_generator.untransformed_coordinates,
                    masses=redundant_generator.masses,
                    relocalize=redundant_generator.relocalize
                )
            if reference_internals is None:
                opts['reference_internals'] = internals
                internals = np.zeros(internals.shape[:-1] + red_tf.shape[-1:], dtype=internals.dtype)
            else:
                disp = internals - reference_internals
                if handle_periodicity:
                    # if the disp changed by > np.pi, assume it wrapped
                    periodic_disps = disp[..., self._periodics]
                    big_changes = np.where(np.abs(periodic_disps) > np.pi)
                    if len(big_changes[0]) > 0:
                        mods = periodic_disps[big_changes] % np.pi
                        disp[..., self._periodics][big_changes] = mods

                internals = nput.vec_tensordot(
                    red_tf,
                    internals - reference_internals,
                    axes=[-2, 0]
                )
            # internals = internals @ red_tf
            opts['derivs'] = red_exp
            opts['redundant_transformation'] = red_tf
            opts['redundant_inverse'] = np.moveaxis(red_tf, -1, -2)

        return internals, opts
# MolecularCartesianToZMatrixConverter = MolecularCartesianToZMatrixConverter()
# MolecularCartesianToZMatrixConverter.register(CoordinateSystemConverters)

class MolecularGICToCartesianConverter(GICSystemToCartesianConverter):
    """
    ...
    """
    def __init__(self, cart_system, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(cart_system, zmat_system)` type pair this converter handles.

        :param cart_system: the target Cartesian coordinate system
        :type cart_system: object
        :param zmat_system: the source generic-internal coordinate system
        :type zmat_system: object
        :param opts: extra options forwarded to the base converter
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (cart_system, zmat_system)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(cart_system, zmat_system)`
        :rtype: tuple
        """
        return self._types
    def convert_many(self,
                     coords, *,
                     redundant_transformation=None,
                     redundant_inverse=None,
                     redundant_generator=None,
                     reference_internals=None,
                     **kw):
        """
        **LLM Docstring**

        Convert a batch of generic-internal coordinates to Cartesians, forwarding any redundant transformation (and its transpose as the corresponding inverse) as the `transformations` argument to the base converter.

        :param coords: the batch of internal coordinates to convert
        :type coords: np.ndarray
        :param redundant_transformation: the redundant-to-non-redundant transformation used to build the internal coordinates, if any
        :type redundant_transformation: np.ndarray | None
        :param redundant_inverse: accepted but not directly used (recomputed from `redundant_transformation` if needed)
        :type redundant_inverse: np.ndarray | None
        :param redundant_generator: accepted for interface consistency but not used in this method's body
        :type redundant_generator: object | None
        :param reference_internals: the reference internal coordinates the redundant transformation is defined relative to
        :type reference_internals: np.ndarray | None
        :param kw: extra keyword arguments forwarded to the base converter's `convert_many`
        :type kw: dict
        :return: `(carts, opts)`
        :rtype: tuple
        """
        # if redundant_inverse is None and redundant_transformation is not None:
        #     redundant_inverse = np.moveaxis(redundant_transformation, -1, -2)
        carts, opts = super().convert_many(coords,
                                           transformations=(
                                               [[redundant_transformation.T], [redundant_transformation]]
                                                 if redundant_transformation is not None else
                                               None
                                           ),
                                           reference_internals=reference_internals,
                                           **kw)
        # if redundant_transformation is not None and 'derivs' in opts:
        #     opts['derivs'] = nput.tensor_reexpand([redundant_inverse], opts['derivs'])

        return carts, opts

class RegularGICToMolecularGICConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, **opts):
        """
        **LLM Docstring**

        Store the `(GenericInternalCoordinates, MolecularGenericInternalCoordinateSystem)` type pair this converter handles.

        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (GenericInternalCoordinates, MolecularGenericInternalCoordinateSystem)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(GenericInternalCoordinates, MolecularGenericInternalCoordinateSystem)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a plain generic-internal system and a molecular generic-internal system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a plain generic-internal system and a molecular generic-internal system share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        return coords, kwargs

class MolecularGICConverterToRegularGIC(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, **opts):
        """
        **LLM Docstring**

        Store the `(MolecularGenericInternalCoordinateSystem, GenericInternalCoordinates)` type pair this converter handles.

        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (MolecularGenericInternalCoordinateSystem, GenericInternalCoordinates)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(MolecularGenericInternalCoordinateSystem, GenericInternalCoordinates)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a molecular generic-internal system and a plain generic-internal system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a molecular generic-internal system and a plain generic-internal system share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        return coords, kwargs

class MolecularIZCoordinateSystem(MolecularZMatrixCoordinateSystem):
    name = "MolecularIZMatrix"

class MolecularCartesianToIZConverter(MolecularCartesianToZMatrixConverter):
    """
    ...
    """
    base_internal_type = McUtils.Coordinerds.IterativeZMatrixCoordinates

class MolecularIZToCartesianConverter(MolecularZMatrixToCartesianConverter):
    """
    ...
    """
    base_internal_type = McUtils.Coordinerds.IterativeZMatrixCoordinates
    def convert_many(self, coords, *,
                     masses, dummy_positions, ref_coords,
                     origins=None, axes=None, ordering=None,
                     reembed=False, axes_choice=None, return_derivs=None,
                     strip_dummies=False,
                     strip_embedding=True,
                     planar_ref_tolerance=None,
                     embedding_masses=None,
                     fixed_atoms=None,
                     fixed_coords=None,
                     **kwargs):
        """
        **LLM Docstring**

        Convert a batch of iterative-Z-matrix coordinates to Cartesians, filling in sensible defaults for the extra dummy embedding masses/fixed atoms/fixed coordinates used by the iterative scheme before delegating to the base `MolecularZMatrixToCartesianConverter.convert_many`.

        :param coords: the batch of iterative-Z-matrix coordinates to convert
        :type coords: np.ndarray
        :param masses: the atomic masses
        :type masses: np.ndarray
        :param dummy_positions: indices of dummy atoms
        :type dummy_positions: list[int]
        :param ref_coords: the reference Cartesian coordinates for the embedding
        :type ref_coords: np.ndarray
        :param origins: the embedding origin points
        :type origins: np.ndarray | None
        :param axes: the embedding reference axes
        :type axes: np.ndarray | None
        :param ordering: the Z-matrix atom ordering
        :type ordering: np.ndarray | None
        :param reembed: whether to re-embed (Eckart-align) during conversion
        :type reembed: bool
        :param axes_choice: which axis pair was chosen for the embedding
        :type axes_choice: tuple | None
        :param return_derivs: which derivative orders to return
        :type return_derivs: object | None
        :param strip_dummies: whether to exclude dummy-atom coordinates
        :type strip_dummies: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates
        :type strip_embedding: bool
        :param planar_ref_tolerance: tolerance for detecting a (near-)planar reference structure
        :type planar_ref_tolerance: float | None
        :param embedding_masses: masses to use for the 3 extra dummy embedding atoms, plus the real atoms; defaults to very heavy dummy masses
        :type embedding_masses: np.ndarray | None
        :param fixed_atoms: indices treated as fixed for the iterative embedding; defaults to the first 3 (dummy) atoms
        :type fixed_atoms: list[int] | None
        :param fixed_coords: coordinate indices treated as fixed; defaults to the standard embedding coordinates offset by the 3 dummy atoms
        :type fixed_coords: list[int] | None
        :param kwargs: extra keyword arguments forwarded to the base converter
        :type kwargs: dict
        :return: `(carts, opts)`, with `opts['masses']` set to the original (non-dummy-augmented) masses
        :rtype: tuple
        """
        if embedding_masses is None:
            embedding_masses = np.concatenate([[1e8] * 3, masses])
        if fixed_atoms is None:
            fixed_atoms = [0, 1, 2]
        if fixed_coords is None:
            fixed_coords = [x+3 for x in MolecularZMatrixCoordinateSystem.embedding_coords]
        carts, opts = super().convert_many(
            coords,
            masses=masses, dummy_positions=dummy_positions, ref_coords=ref_coords,
            origins=origins, axes=axes, ordering=ordering,
            reembed=reembed, axes_choice=axes_choice, return_derivs=return_derivs,
            strip_dummies=strip_dummies,
            strip_embedding=strip_embedding,
            planar_ref_tolerance=planar_ref_tolerance,
            embedding_masses=embedding_masses,
            fixed_atoms=fixed_atoms,
            fixed_coords=fixed_coords,
            **kwargs
        )
        opts['masses'] = masses

        return carts, opts


class RegularIZToMolecularIZConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(IterativeZMatrixCoordinates, zmat_system)` type pair this converter handles.

        :param zmat_system: the molecular iterative-Z-matrix coordinate system being converted to
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (McUtils.Coordinerds.IterativeZMatrixCoordinates, zmat_system)
        super().__init__(**opts)
    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(IterativeZMatrixCoordinates, zmat_system)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a plain iterative-Z-matrix system and a molecular one share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a plain iterative-Z-matrix system and a molecular one share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        return coords, kwargs

class MolecularIZToRegularIZConverter(CoordinateSystemConverter):
    """
    ...
    """


    def __init__(self, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(zmat_system, IterativeZMatrixCoordinates)` type pair this converter handles.

        :param zmat_system: the molecular iterative-Z-matrix coordinate system being converted from
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self._types = (zmat_system, McUtils.Coordinerds.IterativeZMatrixCoordinates)
        super().__init__(**opts)

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(zmat_system, IterativeZMatrixCoordinates)`
        :rtype: tuple
        """
        return self._types

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a molecular iterative-Z-matrix system and a plain one share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a molecular iterative-Z-matrix system and a plain one share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        return coords, kwargs