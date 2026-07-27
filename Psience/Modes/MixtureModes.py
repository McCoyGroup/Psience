"""
Provides support for handling modes that arise from
"""

import enum
import numpy as np, scipy.linalg as slag, collections
from McUtils.Data import AtomData, UnitsData
import numpy as np, scipy
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
                 mass_weighted=False,
                 frequency_scaled=False,
                 g_matrix=None,
                 name=None,
                 ):
        """
        **LLM Docstring**

        Build a `MixtureModes` object -- a `CoordinateSystem` expressing coordinates as a (possibly non-orthogonal) linear rotation of some base coordinate system, with associated frequencies, masses, and mass-weighting/frequency-scaling/G-matrix bookkeeping.

        :param basis: the base coordinate system the modes are expressed in
        :type basis: CoordinateSystem
        :param coeffs: the mode coefficient (coordinates-by-modes) matrix, or a list of such matrices/tensors representing higher-order terms of a nonlinear coordinate expansion (the first entry is used as the linear matrix)
        :type coeffs: np.ndarray | list[np.ndarray]
        :param freqs: the vibrational frequencies associated with each mode
        :type freqs: np.ndarray | None
        :param origin: the reference geometry the modes are expanded about
        :type origin: np.ndarray | None
        :param masses: the atomic masses
        :type masses: np.ndarray | None
        :param inverse: the modes-by-coordinates (inverse) matrix
        :type inverse: np.ndarray | None
        :param mass_weighted: whether the modes are mass-weighted
        :type mass_weighted: bool
        :param frequency_scaled: whether the modes are frequency-scaled (dimensionless)
        :type frequency_scaled: bool
        :param g_matrix: the G-matrix associated with the modes
        :type g_matrix: np.ndarray | None
        :param name: a name for the mode set; defaults to the class's `name` attribute
        :type name: str | None
        :return: None
        :rtype: None
        """
        if (
                isinstance(coeffs, np.ndarray) or
                (not nput.is_numeric(coeffs[0]) and nput.is_numeric(coeffs[0][0]))
        ):
            coeffs = [coeffs]
        coeffs = [np.asanyarray(c) for c in coeffs]
        full_coeffs = coeffs
        coeffs = coeffs[0]
        if inverse is not None:
            inverse = np.asanyarray(inverse)
        if masses is not None:
            masses = np.asanyarray(masses)
        if freqs is not None:
            freqs = np.asanyarray(freqs)
        if g_matrix is not None:
            g_matrix = np.asanyarray(g_matrix)

        super().__init__(
            matrix=coeffs,
            inverse=inverse,
            name=self.name if name is None else name,
            basis=basis,
            dimension=(coeffs.shape[1],),
            origin=origin
        )
        self.freqs = freqs
        self.masses = masses
        self.mass_weighted = mass_weighted
        self.frequency_scaled = frequency_scaled
        self.g_matrix = g_matrix
        self._extended_coeffs = full_coeffs
        self._inverse_coeffs = None


    def to_state(self, serializer=None):
        """
        **LLM Docstring**

        Serialize this mode set's essential data (basis, matrix, inverse, frequencies, masses, mass-weighting/frequency-scaling flags, G-matrix) into a plain dict.

        :param serializer: the serializer used to recursively serialize the `basis` object
        :type serializer: object
        :return: the serialized state dict
        :rtype: dict
        """
        return {
            'basis':serializer.serialize(self.basis),
            'matrix':self.matrix,
            'inverse':self.inverse,
            'freqs':self.freqs,
            'masses':self.masses,
            'mass_weighted':self.mass_weighted,
            'frequency_scaled':self.frequency_scaled,
            'g_matrix':self.g_matrix
        }
    @classmethod
    def from_state(cls, data, serializer=None):
        """
        **LLM Docstring**

        Reconstruct a `MixtureModes` from a previously serialized state dict.

        :param data: the serialized state, as produced by `to_state`
        :type data: dict
        :param serializer: the serializer used to recursively deserialize the `basis` object
        :type serializer: object
        :return: the reconstructed mode set
        :rtype: MixtureModes
        """
        return cls(
            serializer.deserialize(data['basis']),
            data['matrix'],
            inverse=data['inverse'],
            freqs=data['freqs'],
            masses=data['masses'],
            mass_weighted=data['mass_weighted'],
            frequency_scaled=data['frequency_scaled'],
            g_matrix=data['g_matrix']
        )

    @classmethod
    def prep_modes(cls, modes):
        """
        **LLM Docstring**

        Coerce a variety of mode representations (an already-built `MixtureModes`-like object, a dict of matrix/inverse/basis/freqs, an object exposing `.matrix`, or a raw coefficient array) into a proper `MixtureModes` instance, inferring a sensible default basis (Cartesian or generic internal) if none is given.

        :param modes: the mode data to coerce
        :type modes: object | dict | np.ndarray
        :return: the constructed (or passed-through) mode set
        :rtype: MixtureModes
        :raises ValueError: if neither a matrix nor an inverse can be determined from `modes`
        """
        if isinstance(modes, cls) or all(hasattr(modes, b) for b in ['basis', 'matrix', 'inverse', 'freqs']):
            return modes

        matrix = None
        inverse = None
        basis = None
        freqs = None
        if isinstance(modes, dict):
            opts = modes.copy()
            for k in ["matrix", "modes"]:
                if k in opts:
                    matrix = opts[k]
                    del opts[k]
        elif hasattr(modes, 'matrix'):
            matrix = modes.matrix
            opts = {}

            for k in ['inverse', 'basis', 'freqs']:
                if hasattr(modes, k): opts[k] = getattr(modes, k)
        else:
            matrix = np.asanyarray(modes)
            opts = {}

        if 'inverse' in opts:
            inverse = opts['inverse']
            del opts['inverse']
        if 'basis' in opts:
            basis = opts['basis']
            del opts['basis']
        if 'freqs' in opts:
            freqs = opts['freqs']
            del opts['freqs']

        if matrix is None and inverse is None:
            raise ValueError(f"can't prep {cls.__name__} without matrix or inverse")

        if basis is None:
            if matrix.shape[0] % 3 == 0:
                basis = CartesianCoordinateSystem3D
            else:
                basis = InternalCoordinateSystem(dimension=(None, matrix.shape[0]))

        return cls(
            basis,
            matrix,
            inverse=inverse,
            freqs=freqs,
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
        elif isinstance(item, slice):
            ...
            # item = np.arange(self.matrix.shape[1])[item]
        elif not isinstance(item[0], (int, np.integer)):
            item = tuple(item[0])

        sub_modes = self.matrix[:, item]
        inv = self._inv
        if inv is not None:
            inv = inv[item, :]
        freq = self.freqs[item,]
        return self.modify(
            matrix=sub_modes,
            freqs=freq,
            inverse=inv
        )

    def modify(self,
               matrix=None,
               *,
               freqs=None,
               origin=None,
               masses=None,
               inverse=None,
               name=None,
               mass_weighted=None,
               frequency_scaled=None,
               g_matrix=None
               ):
        """
        **LLM Docstring**

        Build a new `MixtureModes` (of the same concrete subclass) with the given fields overridden, defaulting each unspecified field to this object's own current value.

        :param matrix: replacement modes-by-coordinates matrix; defaults to `self.matrix`
        :type matrix: np.ndarray | None
        :param freqs: replacement frequencies; defaults to `self.freqs`
        :type freqs: np.ndarray | None
        :param origin: replacement reference geometry; defaults to `self.origin`
        :type origin: np.ndarray | None
        :param masses: replacement masses; defaults to `self.masses`
        :type masses: np.ndarray | None
        :param inverse: replacement coords-by-modes matrix; defaults to `self.inverse`
        :type inverse: np.ndarray | None
        :param name: replacement name; defaults to `self.name`
        :type name: str | None
        :param mass_weighted: replacement mass-weighting flag; defaults to `self.mass_weighted`
        :type mass_weighted: bool | None
        :param frequency_scaled: replacement frequency-scaling flag; defaults to `self.frequency_scaled`
        :type frequency_scaled: bool | None
        :param g_matrix: replacement G-matrix; defaults to `self.g_matrix`
        :type g_matrix: np.ndarray | None
        :return: the new, modified mode set
        :rtype: MixtureModes
        """
        return type(self)(
            self.basis,
            self.matrix if matrix is None else matrix,
            freqs=self.freqs if freqs is None else freqs,
            origin=self.origin if origin is None else origin,
            masses=self.masses if masses is None else masses,
            inverse=self.inverse if inverse is None else inverse,
            name=self.name if name is None else name,
            mass_weighted=self.mass_weighted if mass_weighted is None else mass_weighted,
            frequency_scaled=self.frequency_scaled if frequency_scaled is None else frequency_scaled,
            g_matrix=self.g_matrix if g_matrix is None else g_matrix,
        )

    def rotate(self, rot, in_place=False):
        """
        **LLM Docstring**

        Intended to rotate the mode set by a given rotation. Not implemented -- the intended semantics were judged too ambiguous.

        :param rot: the rotation to apply
        :type rot: object
        :param in_place: whether to apply the rotation in place
        :type in_place: bool
        :return: never returns
        :rtype: None
        :raises NotImplementedError: always
        """
        raise NotImplementedError("too confusing...")

    def transform(self, tf, inv=None, origin=None):
        """
        **LLM Docstring**

        Intended to apply a linear transformation (and its inverse) to the mode set. Not implemented -- immediately raises before reaching the (dead) implementation code below it, which is retained but unreachable.

        :param tf: the transformation matrix
        :type tf: np.ndarray
        :param inv: the inverse transformation; computed from `tf` if `tf` is square and `inv` isn't given (in the unreachable code)
        :type inv: np.ndarray | None
        :param origin: the reference geometry to transform; defaults to `self.origin` (in the unreachable code)
        :type origin: np.ndarray | None
        :return: never returns
        :rtype: None
        :raises NotImplementedError: always, before the (dead) implementation runs
        :raises ValueError: (in the unreachable code) if `inv` isn't given and `tf` isn't square
        """
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
        """
        **LLM Docstring**

        Whether the modes are expressed relative to a Cartesian (as opposed to internal-coordinate) origin, inferred from the origin's dimensionality.

        :return: `True` if the origin is 2-dimensional (`(natoms, 3)`)
        :rtype: bool
        """
        return self.origin.ndim == 2

    def embed_coords(self, carts):
        """
        **LLM Docstring**

        Convert a batch of Cartesian displacement structures into mode coordinates, by subtracting the reference origin and projecting through the coords-by-modes (inverse) matrix.

        :param carts: the Cartesian structures to embed
        :type carts: np.ndarray
        :return: the mode coordinates
        :rtype: np.ndarray
        """
        flat_carts = (carts - self.origin[np.newaxis]).reshape((len(carts), -1))
        return (flat_carts[:, np.newaxis, :] @ self.inverse.T[np.newaxis]).reshape(
            flat_carts.shape[0],
            self.matrix.shape[1]
        )
    def unembed_coords(self, mode_coords):
        """
        **LLM Docstring**

        Convert a batch of mode coordinates back into Cartesian structures, by projecting through the modes-by-coordinates matrix and adding back the reference origin.

        :param mode_coords: the mode coordinates to unembed
        :type mode_coords: np.ndarray
        :return: the Cartesian structures
        :rtype: np.ndarray
        """
        origin = self.origin
        carts = (mode_coords[:, np.newaxis, :] @ self.matrix.T[np.newaxis]).reshape(
            mode_coords.shape[:1] + origin.shape
        )
        carts = carts + origin[np.newaxis]
        return carts

    @property
    def total_transformation(self):
        """
        **LLM Docstring**

        The full (possibly multi-term, nonlinear) forward coordinate-expansion tensors this mode set was constructed with.

        :return: the list of expansion-order tensors, starting with the linear modes-by-coordinates matrix
        :rtype: list[np.ndarray]
        """
        return self._extended_coeffs
    @property
    def inverse_transformation(self):
        """
        **LLM Docstring**

        The (lazily computed and cached) inverse of `total_transformation`, i.e. the Taylor-series expansion mapping mode coordinates back to base coordinates.

        :return: the inverse expansion-order tensors
        :rtype: list[np.ndarray]
        """
        if self._inverse_coeffs is None:
            self._inverse_coeffs = nput.inverse_transformation(self.total_transformation,
                                                               len(self.total_transformation),
                                                               reverse_expansion=[self.inverse])
        return self._inverse_coeffs
    def embed_derivs(self, derivs):
        """
        **LLM Docstring**

        Re-express a set of derivative tensors (with respect to the base coordinates) in terms of mode coordinates, by re-expanding through `total_transformation`.

        :param derivs: the derivative tensors to re-express
        :type derivs: list[np.ndarray]
        :return: the re-expressed derivative tensors
        :rtype: list[np.ndarray]
        """
        return nput.tensor_reexpand(self.total_transformation, derivs)
    def unembed_derivs(self, derivs):
        """
        **LLM Docstring**

        Re-express a set of derivative tensors (with respect to mode coordinates) back in terms of the base coordinates, by re-expanding through `inverse_transformation`.

        :param derivs: the derivative tensors to re-express
        :type derivs: list[np.ndarray]
        :return: the re-expressed derivative tensors
        :rtype: list[np.ndarray]
        """
        return nput.tensor_reexpand(self.inverse_transformation, derivs)

    @property
    def is_cartesian(self):
        """
        **LLM Docstring**

        Whether the modes are expressed over a Cartesian coordinate space, either inferred from the mode-matrix row count matching `3 * len(masses)` (if masses are known) or from the base coordinate system's name.

        :return: whether the modes are Cartesian-basis
        :rtype: bool
        """
        if self.masses is not None:
            return self.matrix.shape[0] // 3 == len(self.masses)
        else:
            return 'Cartesian' in self.basis.name
    @property
    def coords_by_modes(self):
        """
        **LLM Docstring**

        The coordinates-by-modes (inverse) transformation matrix.

        :return: `self.inverse`
        :rtype: np.ndarray
        """
        return self.inverse
    @property
    def modes_by_coords(self):
        """
        **LLM Docstring**

        The modes-by-coordinates transformation matrix.

        :return: `self.matrix`
        :rtype: np.ndarray
        """
        return self.matrix

    def _eval_G(self, masses):
        """
        **LLM Docstring**

        Compute the coordinate-space G-matrix (inverse effective-mass metric) for this mode set: for Cartesian modes, the diagonal inverse-mass matrix; for internal-coordinate modes, derived from the (already mass-weighted) inverse mode matrix.

        :param masses: the atomic masses to use (only used in the Cartesian branch)
        :type masses: np.ndarray
        :return: the G-matrix
        :rtype: np.ndarray
        :raises NotImplementedError: for non-mass-weighted internal-coordinate modes, which aren't supported
        """
        if not self.is_cartesian:  # a hack
            if self.mass_weighted:
                G = self.inverse.T @ self.inverse
            else:
                # if self.origin is None:
                #     raise ValueError("can't get mass weighting matrix (G^-1/2) without structure")
                # G = self.modes_by_coords @
                raise NotImplementedError("non-mass-weighted internal G-matrix not supported")
        else:
            G = np.diag(1 / np.repeat(masses, 3))
        return G

    def _get_gmatrix(self, masses=None):
        """
        **LLM Docstring**

        Resolve (and cache, if computed fresh for the default masses) the G-matrix and its fractional powers, computing it via `_eval_G` if not already cached and no alternate `masses` are given.

        :param masses: alternate masses to compute a fresh (uncached) G-matrix for; if `None`, uses (and caches into) `self.g_matrix`/`self.masses`
        :type masses: np.ndarray | None
        :return: `(masses, g12, gi12)` -- the masses used, and the G-matrix's square root and inverse square root
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
        """
        if masses is None:
            G = self.g_matrix
        else:
            G = None
        if G is None:
            if masses is None:
                masses = self.masses
                G = self.g_matrix = self._eval_G(self.masses)
            else:
                G = self._eval_G(self.masses)
        g12 = nput.fractional_power(G, 1/2)
        gi12 = nput.fractional_power(G, -1/2)
        return masses, g12, gi12

    @classmethod
    def _local_tf(cls, f, g):
        """
        **LLM Docstring**

        Compute the diagonal rescaling transformation that converts between a set of diagonal force constants `f` and G-matrix elements `g`, used to build "local mode" (Duschinsky-diagonal) representations.

        :param f: the diagonal force-constant values
        :type f: np.ndarray
        :param g: the diagonal G-matrix values
        :type g: np.ndarray
        :return: the diagonal rescaling matrix
        :rtype: np.ndarray
        """
        g1 = np.diag(g)
        f1 = np.diag(f)
        s = np.sign(g1) / np.sign(f1)
        return np.diag(s * np.power(np.abs(g1) / np.abs(f1), 1 / 4))
    @classmethod
    def compute_local_transformations(cls, f, g):
        """
        **LLM Docstring**

        Compute the pair of diagonal rescaling transformations (for the force-constant matrix and its inverse-G counterpart) that define the "local mode" representation from diagonal `f`/`g` values.

        :param f: the diagonal force-constant values
        :type f: np.ndarray
        :param g: the diagonal G-matrix values
        :type g: np.ndarray
        :return: `[local_tf(f, g), local_tf(g, f)]`
        :rtype: list[np.ndarray]
        """
        return [
            cls._local_tf(f, g),
            cls._local_tf(g, f)
        ]

    @classmethod
    def compute_local_hessian(cls, f, g):
        """
        **LLM Docstring**

        Compute the rescaled ("local mode") force-constant matrix from full `f`/`g` matrices, using the diagonal rescaling from `_local_tf`.

        :param f: the force-constant (Hessian) matrix
        :type f: np.ndarray
        :param g: the G-matrix
        :type g: np.ndarray
        :return: the rescaled local Hessian
        :rtype: np.ndarray
        """
        a = cls._local_tf(f, g)
        return a @ f @ a

    @classmethod
    def compute_local_gmatrix(cls, f, g):
        """
        **LLM Docstring**

        Compute the rescaled ("local mode") G-matrix from full `f`/`g` matrices, using the diagonal rescaling from `_local_tf`.

        :param f: the force-constant (Hessian) matrix
        :type f: np.ndarray
        :param g: the G-matrix
        :type g: np.ndarray
        :return: the rescaled local G-matrix
        :rtype: np.ndarray
        """
        ai = cls._local_tf(g, f)
        return ai @ g @ ai

    def compute_hessian(self, system='modes'):
        """
        **LLM Docstring**

        Compute the (signed-frequency-squared) Hessian matrix in either the mode basis or the underlying coordinate basis, by reconstructing it from the stored frequencies and the mode transformation.

        :param system: which basis to compute the Hessian in, `'modes'` or `'coords'`
        :type system: str
        :return: the Hessian matrix
        :rtype: np.ndarray
        :raises ValueError: if `system` is neither `'modes'` nor `'coords'`
        """
        if system == 'modes':
            pinv = self.coords_by_modes @ self.modes_by_coords
            return pinv @ np.diag(np.sign(self.freqs) * (self.freqs ** 2)) @ pinv.T
        elif system == 'coords':
            pinv = self.modes_by_coords
            return pinv @ np.diag(np.sign(self.freqs) * (self.freqs ** 2)) @ pinv.T
        else:
            raise ValueError(f'unknown system for normal modes "{system}", valid are "modes", "coords"')

    def compute_gmatrix(self, system='modes', return_fractional=False):
        """
        **LLM Docstring**

        Compute the G-matrix in either the mode basis or the underlying coordinate basis, using whichever information is available (mass-weighting, a stored G-matrix, or Cartesian masses).

        :param system: which basis to compute the G-matrix in, `'modes'` or `'coords'`
        :type system: str
        :param return_fractional: whether to also return the G-matrix's square root and inverse square root
        :type return_fractional: bool
        :return: the G-matrix, or `(g, g12, gi12)` if `return_fractional` is set; `None` if no G-matrix information is available in the `'modes'` case
        :rtype: np.ndarray | tuple | None
        :raises ValueError: if `system` is neither `'modes'` nor `'coords'`
        """
        if system == 'modes':
            if self.mass_weighted:
                g = self.modes_by_coords.T @ self.modes_by_coords
            elif self.g_matrix is not None:
                g = self.modes_by_coords.T @ self.g_matrix @ self.modes_by_coords
            elif self.is_cartesian:
                g = self.modes_by_coords.T @ np.diag(np.repeat(1/self.masses, 3))  @ self.modes_by_coords
            else:
                return None
            if return_fractional:
                g12 = nput.fractional_power(g, 1/2)
                gi12 = nput.fractional_power(g, -1/2)
                return g, g12, gi12
            else:
                return g
                # raise NotImplementedError("non-mass-weighted internal G-matrix not supported")
        elif system == 'coords':
            g, g12, gi12 = self._get_gmatrix()
            if return_fractional:
                return g, g12, gi12
            else:
                return g
        else:
            raise ValueError(f'unknown system for normal modes "{system}", valid are "modes", "coords"')


    zero_freq_cutoff = 1e-8
    @classmethod
    def _zero_projected_modes(cls, f, g):
        """
        **LLM Docstring**

        Solve the generalized eigenvalue problem `f v = freq^2 g v`, falling back to first projecting out the null space of `g` (via its own eigendecomposition) if the direct generalized eigensolve fails due to `g` being singular.

        :param f: the force-constant matrix
        :type f: np.ndarray
        :param g: the G-matrix
        :type g: np.ndarray
        :return: `(freqs, modes)` -- the signed square-root frequencies and the corresponding eigenvectors
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        try:
            freqs2, modes = scipy.linalg.eigh(f, g, type=3)
        except np.linalg.LinAlgError:
            g_vals, q = np.linalg.eigh(g)
            mask = np.where(np.abs(g_vals) > cls.zero_freq_cutoff)[0]
            g_tf = q[:, mask]
            f_tf = g_tf.T
            g = g_tf @ g @ g_tf.T
            f= f_tf.T @ f @ f_tf
            freqs2, modes = scipy.linalg.eigh(f, g, type=3)
            modes = modes @ f_tf

        return np.sign(freqs2) * np.sqrt(np.abs(freqs2)), modes

    def compute_freqs(self):
        """
        **LLM Docstring**

        Recompute the vibrational frequencies for this mode set from its coordinate-basis Hessian and G-matrix, via a generalized eigenvalue solve.

        :return: the signed square-root frequencies
        :rtype: np.ndarray
        """
        # self = self.remove_frequency_scaling().remove_mass_weighting()
        f = self.compute_hessian()
        g = self.compute_gmatrix()

        freqs2, modes = scipy.linalg.eigh(f, g, type=3)
        return np.sign(freqs2) * np.sqrt(np.abs(freqs2))

        return np.sign(freqs2) * np.sqrt(np.abs(freqs2))
    @property
    def local_hessian(self):
        """
        **LLM Docstring**

        The "local mode" (diagonally rescaled) force-constant matrix for this mode set.

        :return: the local Hessian
        :rtype: np.ndarray
        """
        # self = self.remove_frequency_scaling()
        f = self.compute_hessian()
        g = self.compute_gmatrix()
        return self.compute_local_hessian(f, g)

    @property
    def local_gmatrix(self):
        """
        **LLM Docstring**

        The "local mode" (diagonally rescaled) G-matrix for this mode set.

        :return: the local G-matrix
        :rtype: np.ndarray
        """
        # self = self.remove_frequency_scaling()
        f = self.compute_hessian()
        g = self.compute_gmatrix()
        return self.compute_local_gmatrix(f, g)

    @property
    def local_freqs(self):
        """
        **LLM Docstring**

        The diagonal entries of the local-mode Hessian, giving an approximate per-mode "local" force constant/frequency.

        :return: the local-mode diagonal values
        :rtype: np.ndarray
        """
        return np.diag(self.local_hessian)

    @property
    def local_mode_transformations(self):
        """
        **LLM Docstring**

        The pair of diagonal rescaling transformations mapping between this mode set's coordinate-basis Hessian/G-matrix and their local-mode counterparts.

        :return: `[hessian_scaling, gmatrix_scaling]`
        :rtype: list[np.ndarray]
        """
        f = self.compute_hessian()
        g = self.compute_gmatrix()
        return self.compute_local_transformations(f, g)


    def get_nearest_mode_transform(self,
                                   alternate_modes:np.ndarray,
                                   mass_weighted=None,
                                   atoms=None,
                                   maximum_similarity=True,
                                   unitarize=True,
                                   masses=None
                                   ):
        """
        **LLM Docstring**

        Find the transformation (and its inverse) that best aligns this mode set with an externally supplied set of `alternate_modes`, either via a maximum-similarity (orthogonal Procrustes-style) alignment or via direct projection (optionally unitarized), with optional mass-(un)weighting and atom-restriction beforehand.

        :param alternate_modes: the external mode matrix to align/project against
        :type alternate_modes: np.ndarray
        :param mass_weighted: if given, coerces this mode set to (or from) mass-weighted form before comparing
        :type mass_weighted: bool | None
        :param atoms: restrict the comparison to a subset of atoms via an atom projector (Cartesian modes only)
        :type atoms: Iterable[int] | None
        :param maximum_similarity: whether to use a maximum-similarity (orthogonal Procrustes) alignment rather than direct projection
        :type maximum_similarity: bool
        :param unitarize: whether to unitarize the resulting transformation (direct-projection path only)
        :type unitarize: bool
        :param masses: masses to use for the mass-(un)weighting conversion
        :type masses: np.ndarray | None
        :return: `(tf, inv)` -- the forward and inverse transformation matrices
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        if mass_weighted is not None:
            if not mass_weighted:
                if masses is None:
                    masses = self.masses
                    modes = self.remove_mass_weighting()
                else:
                    modes = self.remove_mass_weighting(masses)
            else:
                if masses is None:
                    masses = self.masses
                    modes = self.make_mass_weighted()
                else:
                    modes = self.make_mass_weighted(masses)
        else:
            modes = self

        inv = modes.coords_by_modes
        modes = modes.modes_by_coords
        if self.is_cartesian and atoms is not None:
            proj = self._atom_projector(modes.shape[0], atoms)
            # if mass_weighted:
            #     gi12 = np.diag(np.repeat(1/np.sqrt(masses), 3))
            #     proj = nput.find_basis(gi12 @ proj @ gi12)
            modes = proj @ modes
            alternate_modes = proj @ alternate_modes

        if maximum_similarity:
            tf = nput.maximum_similarity_transformation(modes, alternate_modes, apply_transformation=False)
            return tf, tf.T
        else:
            tf = inv @ alternate_modes
            if unitarize:
                tf = nput.unitarize_transformation(tf)
                inv = tf.T
            else:
                inv = np.linalg.pinv(tf)
            return tf, inv

    ModeData = collections.namedtuple("ModeData", ['freqs', 'modes', 'inverse'])
    default_zero_freq_cutoff = 1.0e-4  # 20 wavenumbers...

    @classmethod
    def get_normal_modes(cls,
                         f_matrix,
                         mass_spec,
                         # mass_units="AtomicMassUnits",
                         remove_transrot=True,
                         dimensionless=False,
                         mass_weighted=None,
                         zero_freq_cutoff=None,
                         return_gmatrix=False,
                         projector=None
                         ):
        """
        **LLM Docstring**

        Compute ordinary harmonic normal modes from a force-constant matrix and mass specification (masses, element symbols, or a full G-matrix), via a generalized eigenvalue solve on the mass-weighted Hessian, optionally removing near-zero-frequency (translation/rotation) modes, applying a custom constraint projector, mass-weighting and/or frequency-scaling (dimensionless-izing) the result, and supporting batched inputs (including "ragged" batches where a different number of modes survives per entry).

        :param f_matrix: the force-constant (Hessian) matrix/matrices
        :type f_matrix: np.ndarray
        :param mass_spec: the atomic masses (or element symbols), or an already-built G-matrix/mass matrix
        :type mass_spec: np.ndarray | list[str]
        :param remove_transrot: whether to discard near-zero-frequency modes
        :type remove_transrot: bool
        :param dimensionless: whether the resulting modes should be frequency-scaled (dimensionless)
        :type dimensionless: bool
        :param mass_weighted: whether the resulting modes should be mass-weighted; defaults to `dimensionless` when `mass_spec` is a bare mass vector
        :type mass_weighted: bool | None
        :param zero_freq_cutoff: the frequency cutoff below which modes are discarded; defaults to `cls.default_zero_freq_cutoff`
        :type zero_freq_cutoff: float | None
        :param return_gmatrix: whether to also return the (broadcast) mass/G-matrix used
        :type return_gmatrix: bool
        :param projector: an additional constraint projector to apply to the mass-weighted Hessian before diagonalizing
        :type projector: np.ndarray | None
        :return: the mode data (`ModeData(freqs, modes, inverse)`), or `(mode_data, mass_spec)` if `return_gmatrix` is set
        :rtype: object | tuple
        """

        f_matrix = np.asanyarray(f_matrix)
        base_shape = f_matrix.shape[:-2]
        f_matrix = f_matrix.reshape((-1,) + f_matrix.shape[-2:])
        if isinstance(mass_spec[0], str):  # atoms were supplied
            mass_spec = np.array([AtomData[a, "Mass"] for a in mass_spec]) * UnitsData.convert(
                "AtomicMassUnits",
                "AtomicUnitOfMass"
            )
        mass_spec = np.asanyarray(mass_spec)

        if mass_spec.ndim == 1:
            if mass_weighted is None: mass_weighted = dimensionless  # generally do the right thing for Cartesians
            mass_spec = np.broadcast_to(mass_spec[:, np.newaxis], (len(mass_spec), 3)).flatten()
            mass_spec = np.diag(1 / mass_spec)

        mass_spec = np.reshape(mass_spec, (-1,) + mass_spec.shape[-2:])
        if mass_spec.shape[0] < f_matrix.shape[0]:
            mass_spec = np.broadcast_to(mass_spec, f_matrix.shape)

        gi12 = nput.fractional_power(mass_spec, -1 / 2)
        g12 = nput.fractional_power(mass_spec, 1 / 2)
        if projector is None:
            if f_matrix.shape[0] == 1:
                freq2, modes = slag.eigh(f_matrix[0], mass_spec[0], type=3)
                V = gi12[0] @ modes
                V = V[np.newaxis]
                freq2 = freq2[np.newaxis]
                modes = modes[np.newaxis]
            else:
                mw_F = g12 @ f_matrix @ np.moveaxis(g12, -1, -2)
                freq2, V = np.linalg.eigh(mw_F)
                modes = g12 @ V
        else:
            pG12 = projector @ g12
            mw_F = pG12 @ f_matrix @ np.moveaxis(pG12, -1, -2)
            freq2, V = np.linalg.eigh(mw_F)
            modes = g12 @ V

        # modes = g12 @ V
        # therefore, inv = V.T @ gi12
        # and modes[:, nz] = g12 @ V[:, nz]; inv[nz, :] = (V.T)[nz, :] @ gi12 (this is only the _left_ inverse)

        blocks = False
        if remove_transrot:
            if zero_freq_cutoff is None:
                zero_freq_cutoff = cls.default_zero_freq_cutoff
            if zero_freq_cutoff > 0:
                nonzero = np.abs(freq2) >= zero_freq_cutoff ** 2
            else:
                nonzero = np.full(freq2.shape, True)

            if gi12.shape[0] != modes.shape[0]:
                gi12 = np.broadcast_to(gi12, f_matrix.shape)

            # need to check if shape is good
            nz_counts = np.sum(nonzero, axis=1)
            blocks = len(freq2) == 1 or np.sum(np.abs(np.diff(nz_counts))) != 0
            freq2_ = []
            modes_ = []
            inv_ = []
            for nz2, f2, m2, v2, gi in zip(nonzero, freq2, modes, V, gi12):
                freq2_.append(f2[nz2])
                modes_.append(m2[:, nz2])
                inv_.append(v2.T[nz2, :] @ gi)
            freq2 = freq2_
            modes = modes_
            inv = inv_
            del freq2_
            del modes_
            del inv_
            if not blocks:
                freq2 = np.array(freq2)
                modes = np.array(modes)
                inv = np.array(inv)

            # else:
            #     if nz_counts[0] == 0:
            #         inv = np.moveaxis(V, -1, -2) @ gi12
            #     else:
            #         nz_spec = [np.where(nz)[0] for nz in nonzero]
            #         freq2 = freq2[nonzero].reshape((modes.shape[0], nz_counts[0]))
            #         print(nz_spec, modes.shape)
            #         modes = np.moveaxis(
            #             nput.vector_take(
            #                 np.moveaxis(modes, -1, -2),
            #                 nz_spec,
            #                 shared=1
            #             ),
            #             -1, -2
            #         )
            #         inv = np.moveaxis(
            #                 nput.vector_take(
            #                     np.moveaxis(V, 1, 0),
            #                     nz_spec,
            #                     shared=1
            #                 )
            #         print(modes.shape, inv.shape)
        else:
            inv = np.moveaxis(V, -1, -2) @ gi12
        # else:
        #     if np.linalg.det(modes) == 1 and np.allclose(modes.T @ modes, np.eyelen(modes)):
        #         inv = modes.T
        #     else:
        #         inv = np.linalg.inv(modes)

        if blocks:
            freqs = [np.sign(f2) * np.sqrt(np.abs(f2)) for f2 in freq2]
        else:
            freqs = np.sign(freq2) * np.sqrt(np.abs(freq2))

        if mass_weighted:
            if blocks:
                modes_ = []
                inv_ = []
                for m, i, g, gi in zip(modes, inv, g12, gi12):
                    modes_.append(gi @ m)
                    inv_.append(i @ g)
                modes = modes_
                inv = inv_
            else:
                # gi12 = slag.fractional_matrix_power(mass_spec, -1 / 2)
                # g12 = slag.fractional_matrix_power(mass_spec, 1 / 2)
                inv = inv @ g12
                modes = gi12 @ modes
        if dimensionless:
            if blocks:
                modes_ = []
                inv_ = []
                for f, m, i in zip(freqs, modes, inv):
                    weighting = np.sqrt(np.abs(f))
                    modes_.append(m / weighting[np.newaxis, :])
                    inv_.append(i * weighting[:, np.newaxis])
                modes = modes_
                inv = inv_
            else:
                weighting = np.sqrt(np.abs(freqs))
                modes = modes / weighting[:, np.newaxis, :]
                inv = inv * weighting[:, :, np.newaxis]
        # if mode == 'reasonable':
        # modes, inv = inv.T, modes.T

        if blocks:
            modes_ = []
            inv_ = []
            for f, m, i in zip(freqs, modes, inv):
                sorting = np.argsort(f)
                modes_.append(m[:, sorting])
                inv_.append(i[sorting, :])
            modes = modes_
            inv = inv_
            modes, inv = [np.moveaxis(i, -1, -2) for i in inv], [np.moveaxis(m, -1, -2) for m in modes]
            # TODO: iterative reshaping
            if len(base_shape) == 0:
                freqs = freqs[0]
                modes = modes[0]
                inv = inv[0]
        else:
            modes, inv = np.moveaxis(inv, -1, -2), np.moveaxis(modes, -1, -2)
            freqs = np.reshape(freqs, base_shape + freqs.shape[1:])
            modes = np.reshape(modes, base_shape + modes.shape[1:])
            inv = np.reshape(inv, base_shape + inv.shape[1:])

        mode_data = cls.ModeData(freqs, modes, inv)
        if return_gmatrix:
            mass_spec = mass_spec.reshape(base_shape + mass_spec.shape[-2:])
            return mode_data, mass_spec
        else:
            return mode_data

    localization_type = 'ned'
    localization_zero_freq_cutoff = 99 / 219474.56
    def get_projected_localized_mode_transformation(self,
                                                    projectors,
                                                    masses=None, origin=None,
                                                    localization_type=None,
                                                    allow_mode_mixing=False,
                                                    maximum_similarity=False,
                                                    unitarize=True,
                                                    zero_freq_cutoff=None,
                                                    orthogonal_projection=False,
                                                    project_zero_gmatrix_modes=True,
                                                    project_zero_gmatrix_cutoff=1e-8,
                                                    atoms=None #ignored but tedious to have around
                                                    ):
        """
        **LLM Docstring**

        Shared driver behind the various coordinate/atom/fragment-localized-mode constructors: diagonalizes the (optionally translation/rotation-projected) mass-weighted Hessian restricted to (or blocked by) a set of projectors, either mixing all projected subspaces together (`allow_mode_mixing`) or diagonalizing each independently and concatenating the results, then finds the transformation aligning the current modes onto the resulting localized modes via `get_nearest_mode_transform`.

        :param projectors: the projection matrix (or matrices) defining the coordinate subspace(s) to localize within
        :type projectors: list[np.ndarray]
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry to use for the translation/rotation-invariant transformation (`'direct'` localization type only)
        :type origin: np.ndarray | None
        :param localization_type: `'ned'`/other for the standard projected-Hessian approach, or `'direct'` to first factor out translation/rotation via an explicit invariant transformation (Cartesian modes only)
        :type localization_type: str | None
        :param allow_mode_mixing: whether to mix all projected subspaces together into a single diagonalization rather than treating each independently
        :type allow_mode_mixing: bool
        :param maximum_similarity: forwarded to `get_nearest_mode_transform`
        :type maximum_similarity: bool
        :param unitarize: forwarded to `get_nearest_mode_transform`
        :type unitarize: bool
        :param zero_freq_cutoff: the frequency cutoff used when diagonalizing each subspace; defaults to `self.localization_zero_freq_cutoff`
        :type zero_freq_cutoff: float | None
        :param orthogonal_projection: whether non-square entries in `projectors` should be treated as bases for an orthogonal (rather than oblique) projection
        :type orthogonal_projection: bool
        :param project_zero_gmatrix_modes: accepted for interface consistency but not used in this method's body
        :type project_zero_gmatrix_modes: bool
        :param project_zero_gmatrix_cutoff: accepted for interface consistency but not used in this method's body
        :type project_zero_gmatrix_cutoff: float
        :param atoms: accepted for interface consistency with sibling localizers but not used in this method's body
        :type atoms: object | None
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        if zero_freq_cutoff is None:
            zero_freq_cutoff = self.localization_zero_freq_cutoff
        if localization_type is None:
            localization_type = self.localization_type

        if masses is None:
            masses = self.masses
            mw = self.make_mass_weighted()
        else:
            mw = self.make_mass_weighted(masses=masses)

        f = mw.compute_hessian('coords')

        if self.is_cartesian and localization_type == 'direct':
            if origin is None:
                origin = mw.remove_mass_weighting().origin
            tr_tf, tr_inv = nput.translation_rotation_invariant_transformation(
                origin.reshape(-1, 3),
                masses,
                mass_weighted=True
            )
        else:
            tr_tf = tr_inv = np.eye(mw.matrix.shape[0])

        projectors = [
            (
                nput.orthogonal_projection_matrix(p)
                    if orthogonal_projection else
                nput.projection_matrix(p)
            )
                if p.shape[0] != p.shape[1] else
            np.asanyarray(p)
            for p in projectors
        ]

        g_matrix = tr_tf.T @ tr_tf
        if allow_mode_mixing:
            proj = np.sum(projectors, axis=0)
            # proj = proj @ g12
            # proj = self._atom_projector(len(masses), atoms)
            freqs, modes, inv = self.get_normal_modes(
                tr_inv @ proj @ f @ proj.T @ tr_inv.T,
                g_matrix,
                remove_transrot=True,
                mass_weighted=True,
                zero_freq_cutoff=zero_freq_cutoff
            )
            # freqs1, _, _ = self.get_normal_modes(
            #      f,
            #     g_matrix,
            #     remove_transrot=True,
            #     mass_weighted=True,
            #     zero_freq_cutoff=zero_freq_cutoff
            # )
            # raise Exception(freqs * 219474.63, freqs1 * 219474.63)
        else:
            # freqs = []
            modes = []
            # inv = []
            for proj in projectors:
                f_matrix = tr_inv @ proj @ f @ proj.T @ tr_inv.T
                # import McUtils.Formatters as mfmt
                # print(
                #     mfmt.TableFormatter("{:.3f}").format(
                #         f_matrix * 219474.56
                #     )
                # )
                sub_freqs, sub_modes, sub_inv = self.get_normal_modes(
                    f_matrix,
                    g_matrix,
                    remove_transrot=True,
                    mass_weighted=True,
                    zero_freq_cutoff=zero_freq_cutoff
                )
                # print(sub_freqs * 219474.56)

                # freqs.append(sub_freqs)
                modes.append(sub_modes)
                # inv.append(sub_inv)
            # freqs = np.concatenate(freqs)
            modes = np.concatenate(modes, axis=1)
            # inv = np.concatenate(inv, axis=0)

        if localization_type == 'direct':
            modes = tr_tf @ modes

        return mw.get_nearest_mode_transform(
            modes,
            mass_weighted=True,
            maximum_similarity=maximum_similarity,
            unitarize=unitarize
        )

    def get_atom_localized_mode_transformation(self,
                                               atoms,
                                               masses=None, origin=None,
                                               localization_type='ned',
                                               allow_mode_mixing=False,
                                               maximum_similarity=False,
                                               orthogonal_projection=False,
                                               unitarize=True,
                                               zero_freq_cutoff=None
                                               ):
        """
        **LLM Docstring**

        Build a localized-mode transformation that concentrates vibrational character onto (or, if `orthogonal_projection`, away from) a specified set of atoms, via `get_projected_localized_mode_transformation` with per-atom (or combined) atom-selection projectors.

        :param atoms: the atom index (or indices) to localize modes onto
        :type atoms: int | Iterable[int]
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry, forwarded to `get_projected_localized_mode_transformation`
        :type origin: np.ndarray | None
        :param localization_type: the localization scheme to use
        :type localization_type: str
        :param allow_mode_mixing: whether to build one combined projector spanning all `atoms` rather than one per atom
        :type allow_mode_mixing: bool
        :param maximum_similarity: forwarded to `get_projected_localized_mode_transformation`
        :type maximum_similarity: bool
        :param orthogonal_projection: whether to localize onto the complement of `atoms` (all other atoms) instead of `atoms` itself
        :type orthogonal_projection: bool
        :param unitarize: forwarded to `get_projected_localized_mode_transformation`
        :type unitarize: bool
        :param zero_freq_cutoff: forwarded to `get_projected_localized_mode_transformation`
        :type zero_freq_cutoff: float | None
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        if nput.is_numeric(atoms):
            atoms = [atoms]
        if masses is None:
            m = self.masses
        else:
            m = masses
        nats = len(m)
        if orthogonal_projection:
            atoms = np.setdiff1d(np.arange(nats), atoms)

        proj = (
            [self._atom_projector(nats, atoms)]
                if allow_mode_mixing else
            [self._atom_projector(nats, a) for a in atoms]
        )
        # gi12 = np.diag(np.repeat(1 / np.sqrt(m), 3))
        # proj = nput.find_basis(gi12 @ proj @ gi12)

        return self.get_projected_localized_mode_transformation(
            proj,
            origin=origin,
            masses=masses,
            localization_type=localization_type,
            maximum_similarity=maximum_similarity,
            unitarize=unitarize,
            allow_mode_mixing=allow_mode_mixing,
            zero_freq_cutoff=zero_freq_cutoff
        )

    def get_fragment_localized_mode_transformation(self,
                                                   fragment,
                                                   masses=None, origin=None,
                                                   localization_type='ned',
                                                   allow_mode_mixing=True,
                                                   maximum_similarity=False,
                                                   orthogonal_projection=False,
                                                   unitarize=True,
                                                   zero_freq_cutoff=None,
                                                   **etc
                                                   ):
        """
        **LLM Docstring**

        Build a localized-mode transformation for a molecular fragment (a set of atoms), via `get_atom_localized_mode_transformation`.

        :param fragment: the atom indices making up the fragment
        :type fragment: Iterable[int]
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry
        :type origin: np.ndarray | None
        :param localization_type: the localization scheme to use
        :type localization_type: str
        :param allow_mode_mixing: whether to build one combined projector for the whole fragment
        :type allow_mode_mixing: bool
        :param maximum_similarity: forwarded through to the underlying localizer
        :type maximum_similarity: bool
        :param orthogonal_projection: whether to localize onto the complement of the fragment instead
        :type orthogonal_projection: bool
        :param unitarize: forwarded through to the underlying localizer
        :type unitarize: bool
        :param zero_freq_cutoff: forwarded through to the underlying localizer
        :type zero_freq_cutoff: float | None
        :param etc: extra options forwarded to `get_atom_localized_mode_transformation`
        :type etc: dict
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        return self.get_atom_localized_mode_transformation(
            fragment,
            masses=masses, origin=origin,
            localization_type=localization_type,
            allow_mode_mixing=allow_mode_mixing,
            maximum_similarity=maximum_similarity,
            orthogonal_projection=orthogonal_projection,
            unitarize=unitarize,
            zero_freq_cutoff=zero_freq_cutoff,
            **etc
        )

    def get_coordinate_projected_localized_mode_transformation(self,
                                                               coordinate_constraints,
                                                               atoms=None,
                                                               masses=None, origin=None,
                                                               localization_type='ned',
                                                               allow_mode_mixing=False,
                                                               maximum_similarity=False,
                                                               orthogonal_projection=True,
                                                               unitarize=True
                                                               ):
        """
        **LLM Docstring**

        Build a localized-mode transformation that concentrates vibrational character along a set of internal-coordinate directions (bond/angle/dihedral bases for Cartesian modes, or by coordinate index for internal-coordinate modes), optionally restricted to a subset of atoms, via `get_projected_localized_mode_transformation`.

        :param coordinate_constraints: the internal-coordinate specification(s) (for Cartesian modes) or coordinate index/indices (for internal-coordinate modes) to localize along
        :type coordinate_constraints: object | Iterable[int]
        :param atoms: restrict the resulting projector(s) to this subset of atoms (Cartesian modes only)
        :type atoms: Iterable[int] | None
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry, forwarded to `get_projected_localized_mode_transformation`
        :type origin: np.ndarray | None
        :param localization_type: the localization scheme to use
        :type localization_type: str
        :param allow_mode_mixing: whether to combine all coordinate projectors into one before localizing
        :type allow_mode_mixing: bool
        :param maximum_similarity: forwarded to `get_projected_localized_mode_transformation`
        :type maximum_similarity: bool
        :param orthogonal_projection: whether the coordinate bases define orthogonal (rather than oblique) projections
        :type orthogonal_projection: bool
        :param unitarize: forwarded to `get_projected_localized_mode_transformation`
        :type unitarize: bool
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        :raises ValueError: if `coordinate_constraints` isn't recognized as coordinate indices for internal-coordinate modes
        """


        if self.is_cartesian:
            if nput.is_numeric(coordinate_constraints[0]):
                coordinate_constraints = [coordinate_constraints]
            nmw = self.remove_mass_weighting()
            if masses is None:
                m = self.masses
            else:
                m = masses
            basis, _, _ = nput.internal_basis(nmw.origin.reshape((-1, 3)), coordinate_constraints, masses=m)

            g12 = np.diag(np.repeat(1/np.sqrt(m), 3))
            projections = [
                nput.projection_matrix(g12 @ b)
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(g12 @ b)
                for b in basis
            ]

            if allow_mode_mixing and orthogonal_projection:
                proj = projections[0]
                for p in projections[1:]:
                    proj = p @ proj @ p
                projections = [proj]

            if atoms is not None:
                if nput.is_numeric(atoms):
                    atoms = [atoms]
                nats = len(m)
                a_proj = self._atom_projector(nats, atoms)
                projections = [
                    a_proj @ proj @ a_proj
                    for proj in projections
                ]
        else:
            if nput.is_numeric(coordinate_constraints):
                coordinate_constraints = [coordinate_constraints]
            if not (
                # Fast check
                isinstance(coordinate_constraints, np.ndarray)
                and np.issubdtype(coordinate_constraints.dtype, np.integer)
            ) and not all(nput.is_int(c) for c in coordinate_constraints):
                raise ValueError(f"internal modes can only be constrained by coordinate index (got {coordinate_constraints})")
            basis = np.zeros((len(coordinate_constraints),self.coords_by_modes.shape[1] ), dtype=float)
            coordinate_constraints = np.asanyarray(coordinate_constraints)
            for n,i in enumerate(coordinate_constraints):
                basis[n, i] = 1


            # g12 = nput.fractional_power(self.g_matrix, 1/2)
            # projections = [
            #     nput.projection_matrix(g12 @ b[:, np.newaxis])
            #         if not orthogonal_projection else
            #     nput.orthogonal_projection_matrix(g12 @ b[:, np.newaxis])
            #     for b in basis
            # ]

            projections = [
                nput.projection_matrix(b[:, np.newaxis])
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(b[:, np.newaxis])
                for b in basis
            ]

        # if orthogonal_projection:
        #     projections = [proj @ gi12 for proj in projections]
        # else:
        #     projections = [proj @ gi12 for proj in projections]


        return self.get_projected_localized_mode_transformation(
            projections,
            origin=origin,
            masses=masses,
            localization_type=localization_type,
            maximum_similarity=maximum_similarity,
            unitarize=unitarize,
            allow_mode_mixing=allow_mode_mixing
        )

    def get_internal_localized_mode_transformation(
            self,
            expansion_coordinates: "Iterable[Iterable[int]|dict]",
            fixed_atoms=None,
            mass_weighted=False,
            project_transrot=True,
            atoms=None,
            maximum_similarity=False,
            orthogonal_projection=False,
            projection=False,
            allow_mode_mixing=False,
            unitarize=True,
            origin=None,
            masses=None,
            localization_type='ned'
    ):
        """
        **LLM Docstring**

        Build a localized-mode transformation aligned with a set of internal-coordinate displacement directions (computed at the reference geometry), either by projecting the Hessian onto (or away from) those directions and re-diagonalizing (`projection`/`orthogonal_projection`), or by directly finding the nearest-mode alignment to the raw coordinate-derivative directions (optionally mass-weighted and/or translation/rotation-projected).

        :param expansion_coordinates: the internal coordinate(s) whose Cartesian derivative directions define the localization target
        :type expansion_coordinates: Iterable[Iterable[int]] | dict
        :param fixed_atoms: atoms to hold fixed when computing the internal-coordinate derivative tensors
        :type fixed_atoms: Iterable[int] | None
        :param mass_weighted: whether to mass-weight the coordinate-derivative directions
        :type mass_weighted: bool
        :param project_transrot: whether to project out translational/rotational components from the derivative directions (direct-alignment path only)
        :type project_transrot: bool
        :param atoms: restrict the projector(s)/alignment to a subset of atoms
        :type atoms: Iterable[int] | None
        :param maximum_similarity: forwarded to the underlying localizer
        :type maximum_similarity: bool
        :param orthogonal_projection: whether to use orthogonal (rather than oblique) projections in the projection-based path
        :type orthogonal_projection: bool
        :param projection: whether to use the projection-and-rediagonalize approach rather than direct nearest-mode alignment
        :type projection: bool
        :param allow_mode_mixing: whether to combine all coordinate directions into a single projected subspace (projection-based path)
        :type allow_mode_mixing: bool
        :param unitarize: forwarded to the underlying localizer
        :type unitarize: bool
        :param origin: reference geometry to use
        :type origin: np.ndarray | None
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param localization_type: the localization scheme to use for the projection-based path
        :type localization_type: str
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        # from .ObliqueModes import ObliqueModeGenerator

        nmw = self.remove_mass_weighting()
        base_derivs = nput.internal_coordinate_tensors(
            nmw.origin,
            expansion_coordinates,
            fixed_atoms=fixed_atoms
        )

        if projection or orthogonal_projection:
            if masses is None:
                m = self.masses
            else:
                m = masses
            base_derivs = np.diag(np.repeat(1 / np.sqrt(m), 3)) @ base_derivs

            if origin is None:
                origin = self.origin
            if masses is None:
                m = self.masses
            else:
                m = masses
            projector = nput.translation_rotation_projector(
                origin.reshape(-1, 3),
                m,
                mass_weighted=mass_weighted
            )
            base_derivs = np.moveaxis(projector, -1, -2) @ base_derivs

            # g12 = np.diag(np.repeat(np.sqrt(m), 3))
            gi12 = np.diag(np.repeat(1 / np.sqrt(m), 3))
            if allow_mode_mixing:
                projections = [
                    nput.projection_matrix(gi12 @ base_derivs)
                        if not orthogonal_projection else
                    nput.orthogonal_projection_matrix(gi12 @ base_derivs)
                ]
            else:
                projections = [
                    nput.projection_matrix(gi12 @ b[:, np.newaxis])
                        if not orthogonal_projection else
                    nput.orthogonal_projection_matrix(gi12 @ b[:, np.newaxis])
                    for b in base_derivs.T
                ]

            # if orthogonal_projection:
            #     projections = [proj @ gi12 for proj in projections]
            # else:
            #     projections = [proj @ gi12 for proj in projections]

            if atoms is not None:
                if nput.is_numeric(atoms):
                    atoms = [atoms]
                nats = len(m)
                a_proj = self._atom_projector(nats, atoms)
                projections = [
                    a_proj @ proj @ a_proj
                    for proj in projections
                ]

            return self.get_projected_localized_mode_transformation(
                projections,
                origin=origin,
                masses=masses,
                localization_type=localization_type,
                maximum_similarity=maximum_similarity,
                unitarize=unitarize,
                allow_mode_mixing=allow_mode_mixing
            )

        else:
            # if project_transrot: mass_weighted = True

            if mass_weighted:
                if masses is None:
                    m = self.masses
                else:
                    m = masses
                base_derivs = np.diag(np.repeat(1 / np.sqrt(m), 3)) @ base_derivs

            if project_transrot:
                if origin is None:
                    origin = self.origin
                if masses is None:
                    m = self.masses
                else:
                    m = masses
                projector = nput.translation_rotation_projector(
                    origin.reshape(-1, 3),
                    m,
                    mass_weighted=mass_weighted,
                    direction='reverse'
                )
                base_derivs = projector @ base_derivs

            return self.get_nearest_mode_transform(
                base_derivs,
                atoms=atoms,
                mass_weighted=mass_weighted,
                maximum_similarity=maximum_similarity,
                unitarize=unitarize
            )

    def get_displacement_localized_mode_transformation(self,
                                                       mode_blocks=None,
                                                       atoms=None,
                                                       mass_weighted=True,
                                                       unitarize=True,
                                                       **maximizer_opts
                                                       ):
        """
        **LLM Docstring**

        Build a unitary localized-mode transformation by iteratively (Jacobi-style) rotating blocks of modes to maximize a displacement-localization criterion (via `nput.jacobi_maximize`/`nput.displacement_localizing_rotation_generator`), optionally restricted to a subset of atoms and/or applied independently within specified mode blocks.

        :param mode_blocks: groups of mode indices to localize independently; a single flat list of indices localizes all modes together, and `None` localizes all modes as one block
        :type mode_blocks: Iterable[int] | list[Iterable[int]] | None
        :param atoms: restrict the localization criterion to a subset of atoms
        :type atoms: Iterable[int] | None
        :param mass_weighted: whether to localize the mass-weighted (rather than un-mass-weighted) mode matrix
        :type mass_weighted: bool
        :param unitarize: must be `True`; the Jacobi-rotation scheme only produces unitary transformations
        :type unitarize: bool
        :param maximizer_opts: extra options forwarded to `nput.jacobi_maximize`
        :type maximizer_opts: dict
        :return: `(tf_inv, tf_inv.T)`, the (unitary) localizing transformation's inverse and its transpose (used as the forward transformation)
        :rtype: tuple[np.ndarray, np.ndarray]
        :raises ValueError: if `unitarize` is `False`
        """
        if not unitarize:
            raise ValueError("Jacobi iterations only give a unitary transformation")

        if mass_weighted:
            modes = self.make_mass_weighted().matrix
        else:
            modes = self.remove_mass_weighting().matrix

        if atoms is not None:
            modes = self._atom_projector(modes.shape[0], atoms) @ modes
        if mode_blocks is None:
            mode_blocks = [np.arange(modes.shape[1])]
        elif nput.is_numeric(mode_blocks[0]):
            mode_blocks = [mode_blocks]

        tf_inv = np.zeros((modes.shape[1], modes.shape[1]), dtype=float)
        for block in mode_blocks:
            _, sub_inv, _ = nput.jacobi_maximize(
                modes[:, block],
                nput.displacement_localizing_rotation_generator,
                **maximizer_opts
            )
            tf_inv[np.ix_(block, block)] = sub_inv

        return tf_inv, tf_inv.T

    def get_mass_scaled_mode_transformation(self,
                                            mass_scaling,
                                            *,
                                            atoms,
                                            localization_cutoff=.8,
                                            num_modes=None,
                                            project_transrot=False,
                                            unitarize=True,
                                            **diag_opts
                                            ):
        """
        **LLM Docstring**

        Build a localized-mode transformation by artificially scaling the masses of a set of atoms (making them effectively heavier or lighter) and re-diagonalizing, then selecting the resulting modes that are most concentrated on those atoms (by projected-displacement norm), optionally applying a localization cutoff or a fixed mode count.

        :param mass_scaling: the scale factor (or per-atom scale factors) to apply to the selected atoms' masses
        :type mass_scaling: float | np.ndarray
        :param atoms: the atom(s) whose masses should be scaled
        :type atoms: int | Iterable[int]
        :param localization_cutoff: minimum normalized displacement-on-atoms required for a mode to be selected; if `None`, the top `num_modes` (or `3*len(atoms)`) modes by that metric are taken instead
        :type localization_cutoff: float | None
        :param num_modes: the number of modes to select; defaults to `3 * len(atoms)`
        :type num_modes: int | None
        :param project_transrot: whether to project out translation/rotation using the mass-scaled masses before diagonalizing
        :type project_transrot: bool
        :param unitarize: accepted for interface consistency with sibling localizers but not used in this method's body
        :type unitarize: bool
        :param diag_opts: extra options forwarded to `NormalModes.get_normal_modes`
        :type diag_opts: dict
        :return: `(tf, tf.T)`, the selected mode transformation and its transpose, or `(None, None)` if no modes clear `localization_cutoff`
        :rtype: tuple[np.ndarray, np.ndarray] | tuple[None, None]
        """
        from .NormalModes import NormalModes
        self = self.make_mass_weighted()

        f = self.compute_hessian('coords')

        if not nput.is_numeric(mass_scaling):
            mass_scaling = np.asanyarray(mass_scaling)

        if nput.is_int(atoms):
            atoms = [atoms]
        atoms = tuple(atoms)

        scaled_masses = np.array(self.masses)
        scaled_masses[atoms,] *= mass_scaling

        if project_transrot:
            proj = nput.translation_rotation_projector(
                self.remove_mass_weighting().origin,
                masses=scaled_masses,
                mass_weighted=True
            )
        else:
            proj = None
        freqs0, q1, _ = NormalModes.get_normal_modes(f, scaled_masses,
                                                     projector=proj,
                                                     mass_weighted=True,
                                                     **diag_opts
                                                     )

        atom_inds = np.concatenate([
            n * 3 + np.arange(3)
            for n in atoms
        ])
        atom_proj = q1[atom_inds, :]
        max_local_modes = np.linalg.norm(atom_proj, axis=0)

        max_num_modes = len(atoms) * 3
        if localization_cutoff is None:
            num_modes = max_num_modes if num_modes is None else num_modes
            mode_pos = np.argsort(-max_local_modes)[:num_modes]
        else:
            mode_pos = np.where(max_local_modes > localization_cutoff)
            if len(mode_pos) == 0 or len(mode_pos[0]) == 0:
                return None, None

            mode_vals = max_local_modes[mode_pos]
            num_modes = max_num_modes if num_modes is None else num_modes
            ord = np.argsort(-mode_vals)[:num_modes]
            mode_pos = mode_pos[0][ord,]

        tf = self.coords_by_modes@q1[:, mode_pos]
        ord2 = np.argsort(np.argmax(np.abs(tf), axis=0))
        tf = tf[:, ord2]

        return tf, tf.T


    class LocalizationMethods(enum.Enum):
        MaximumSimilarity = 'target_modes'
        AtomLocalized = 'atoms'
        FragmentLocalized = 'fragment'
        DisplacmentMinimized = 'displacements'
        Internals = 'coordinates'
        CoordinateConstraints = 'constraints'
        Projected = 'projections'
        MassScaled = 'mass_scaling'

    @property
    def localizer_dispatch(self):
        """
        **LLM Docstring**

        The mapping from `LocalizationMethods` value to the `(constructor_method, primary_argument_name)` pair used by `localize` to dispatch a localization request to the right underlying method.

        :return: the method-name-to-`(constructor, arg_name)` mapping
        :rtype: dict
        """
        return {
            self.LocalizationMethods.MaximumSimilarity.value:(self.get_nearest_mode_transform, 'target_modes'),
            self.LocalizationMethods.AtomLocalized.value:(self.get_atom_localized_mode_transformation, 'atoms'),
            self.LocalizationMethods.FragmentLocalized.value:(self.get_fragment_localized_mode_transformation, 'fragment'),
            self.LocalizationMethods.DisplacmentMinimized.value:(self.get_displacement_localized_mode_transformation, 'mode_blocks'),
            self.LocalizationMethods.Internals.value:(self.get_internal_localized_mode_transformation, 'internals'),
            self.LocalizationMethods.CoordinateConstraints.value:(self.get_coordinate_projected_localized_mode_transformation, 'coordinate_constraints'),
            self.LocalizationMethods.Projected.value:(self.get_projected_localized_mode_transformation, 'projections'),
            self.LocalizationMethods.MassScaled.value:(self.get_mass_scaled_mode_transformation, 'mass_scaling'),
        }

    localization_options = (
        "method",
        "atoms",
        "fragment",
        "masses",
        "target_modes",
        "internals",
        "mode_blocks",
        "coordinate_constraints",
        "projections",
        "reorthogonalize",
        "mass_scaling",
        "unitarize",
        "localization_cutoff",
        "num_modes",
        "project_transrot"
        "mass_weighted",
        "maximum_similarity"
    )
    def localize(self,
                 method=None,
                 *,
                 atoms=None,
                 fragment=None,
                 target_modes=None,
                 internals=None,
                 mode_blocks=None,
                 coordinate_constraints=None,
                 projections=None,
                 reorthogonalize=None,
                 mass_scaling=None,
                 unitarize=True,
                 allow_mode_mixing=False,
                 project_zero_gmatrix_modes=None,
                 project_zero_gmatrix_cutoff=1e-8,
                 **opts
                 ):
        """
        **LLM Docstring**

        Top-level entry point for building a `LocalizedModes` object: infers which localization method to use from whichever keyword argument was supplied (if `method` isn't given explicitly), dispatches to the corresponding `get_*_localized_mode_transformation` method, optionally projects out (and, if mixing is allowed, re-diagonalizes within) any modes with a near-zero effective G-matrix eigenvalue, optionally re-orthogonalizes the result, and wraps the resulting transformation in a `LocalizedModes` object.

        :param method: which localization method to use; inferred from the other arguments if not given
        :type method: str | LocalizationMethods | None
        :param atoms: atom(s) to localize onto, for the `'atoms'` method (or as a restriction for others)
        :type atoms: int | Iterable[int] | None
        :param fragment: fragment atoms to localize onto, for the `'fragment'` method
        :type fragment: Iterable[int] | None
        :param target_modes: external modes to align to, for the `'target_modes'` method
        :type target_modes: np.ndarray | None
        :param internals: internal coordinates to localize along, for the `'coordinates'` method
        :type internals: object | None
        :param mode_blocks: mode index groupings, for the `'displacements'` method
        :type mode_blocks: object | None
        :param coordinate_constraints: coordinate constraints, for the `'constraints'` method
        :type coordinate_constraints: object | None
        :param projections: explicit projectors, for the `'projections'` method
        :type projections: list[np.ndarray] | None
        :param reorthogonalize: whether to re-orthogonalize the localized modes' mass-weighted representation after localization; defaults to `not unitarize`
        :type reorthogonalize: bool | None
        :param mass_scaling: mass-scaling factor(s), for the `'mass_scaling'` method
        :type mass_scaling: float | np.ndarray | None
        :param unitarize: whether the underlying localizer should produce a unitary transformation
        :type unitarize: bool
        :param allow_mode_mixing: whether the underlying localizer is allowed to mix modes across different projected subspaces
        :type allow_mode_mixing: bool
        :param project_zero_gmatrix_modes: whether to project out (and handle) modes with a near-zero effective G-matrix eigenvalue after localizing; defaults to `allow_mode_mixing`
        :type project_zero_gmatrix_modes: bool | None
        :param project_zero_gmatrix_cutoff: the G-matrix eigenvalue magnitude below which a mode is treated as zero
        :type project_zero_gmatrix_cutoff: float
        :param opts: extra options forwarded to the dispatched localizer method
        :type opts: dict
        :return: the constructed `LocalizedModes` object, or `None` if the underlying localizer returned no transformation
        :rtype: LocalizedModes | None
        """

        from .LocalizedModes import LocalizedModes

        if method is None:
            if target_modes is not None:
                method = 'target_modes'
            elif internals is not None:
                method = 'coordinates'
            elif mode_blocks is not None:
                method = 'displacements'
            elif mass_scaling is not None:
                method = 'mass_scaling'
            elif atoms is not None:
                method = 'atoms'
            elif fragment is not None:
                method = 'fragment'
            elif coordinate_constraints is not None:
                method = 'constraints'
            elif projections is not None:
                method = 'projections'
            else:
                method = 'displacements'

        if isinstance(method, str):
            method = self.LocalizationMethods(method)

        if hasattr(method, 'value'):
            method = self.localizer_dispatch.get(method.value, method)

        args = ()
        try:
            method, arg_names = method
        except TypeError:
            ...
        else:
            if isinstance(arg_names, str):
                arg_names = [arg_names]
            all_kw = dict(
                atoms=atoms,
                target_modes=target_modes,
                internals=internals,
                mode_blocks=mode_blocks,
                reorthogonalize=reorthogonalize,
                coordinate_constraints=coordinate_constraints,
                unitarize=unitarize,
                projections=projections,
                mass_scaling=mass_scaling
            )
            args = tuple(all_kw.get(k) for k in arg_names)
            if 'atoms' not in arg_names:
                opts['atoms'] = atoms # taken by all of the localizers
        tf, inv = method(*args,
                         unitarize=unitarize,
                         allow_mode_mixing=allow_mode_mixing,
                         **opts)
        if tf is None: return None
        # inverse = tf.T @ self.inverse # assumes a unitary localization

        if project_zero_gmatrix_modes is None:
            project_zero_gmatrix_modes = allow_mode_mixing
        if project_zero_gmatrix_modes:
            g0 = self.compute_gmatrix()
            g = tf.T @ g0 @ tf
            g_vals, q = np.linalg.eigh(g)
            mask = np.where(np.abs(g_vals) > project_zero_gmatrix_cutoff)[0]
            g_tf = q[:, mask]
            # g_tf = nput.maximum_similarity_transformation(g_tf, np.eye(q.shape[0]), apply_transformation=True)
            # print(mat.shape, g_tf.shape, tf.shape, q_g.shape)
            tf = tf @ g_tf
            inv = g_tf.T @ inv

            if allow_mode_mixing:
                f = self.compute_hessian()
                g = g0

                g = tf.T @ g @ tf
                f = inv @ f @ inv.T

                _, modes = scipy.linalg.eigh(f, g, type=3)
                tf = tf @ modes
                inv = modes.T @ inv
                # return np.sign(freqs2) * np.sqrt(np.abs(freqs2))

        if reorthogonalize is None:
            reorthogonalize = not unitarize
        if reorthogonalize:
            modes = self.make_mass_weighted().matrix @ tf
            g = modes.T @ modes
            tf = tf @ nput.fractional_power(g, -1 / 2)
            inv = nput.fractional_power(g, 1 / 2) @ inv

        modes = LocalizedModes(self, tf, inverse=inv)
        return modes

    def make_mass_weighted(self, masses=None):
        """
        **LLM Docstring**

        Build a mass-weighted version of this mode set by rescaling the mode/inverse matrices and origin through the G-matrix's square root/inverse-square-root; returns `self` unchanged if already mass-weighted.

        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :return: the mass-weighted mode set
        :rtype: MixtureModes
        """
        if self.mass_weighted: return self
        masses, g12, gi12 = self._get_gmatrix(masses=masses)
        L = g12 @ self.modes_by_coords
        Linv = self.coords_by_modes @ gi12
        origin = (self.origin.flatten()[np.newaxis, :] @ gi12).reshape(self.origin.shape)

        return self.modify(L,
                           inverse=Linv,
                           masses=masses,
                           origin=origin,
                           mass_weighted=True
                           )
    def remove_mass_weighting(self, masses=None):
        """
        **LLM Docstring**

        Build a non-mass-weighted version of this mode set by undoing the mass-weighting rescaling; returns `self` unchanged if not already mass-weighted.

        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :return: the non-mass-weighted mode set
        :rtype: MixtureModes
        """
        if not self.mass_weighted: return self
        masses, g12, gi12 = self._get_gmatrix(masses=masses)
        L = gi12 @ self.modes_by_coords
        Linv = self.coords_by_modes @ g12
        origin = (self.origin.flatten()[np.newaxis, :] @ g12).reshape(self.origin.shape)

        return self.modify(L,
                           inverse=Linv,
                           masses=masses,
                           origin=origin,
                           mass_weighted=False
                           )

    def _frequency_scaling(self, freqs=None):
        """
        **LLM Docstring**

        Compute the per-mode scaling factor (`sqrt(|freq|)`) used to convert between frequency-scaled and non-frequency-scaled mode bases.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.freqs`
        :type freqs: np.ndarray | None
        :return: `(freqs, conv)` -- the frequencies used and their corresponding scaling factors
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        # L = self.matrix.shape.T
        if freqs is None:
            freqs = self.freqs
        conv = np.sqrt(np.abs(freqs))
        return freqs, conv

    def make_frequency_scaled(self, freqs=None):
        """
        **LLM Docstring**

        Build a frequency-scaled (dimensionless) version of this mode set by rescaling the mode/inverse matrices by the per-mode frequency-scaling factor; returns `self` unchanged if already frequency-scaled.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.freqs`
        :type freqs: np.ndarray | None
        :return: the frequency-scaled mode set
        :rtype: MixtureModes
        """
        if self.frequency_scaled: return self
        freqs, conv = self._frequency_scaling(freqs=freqs)
        L = self.matrix * conv[np.newaxis, :]
        Linv = self.inverse / conv[:, np.newaxis] # del_Q X

        return self.modify(L,
                           inverse=Linv,
                           freqs=freqs,
                           frequency_scaled=True
                           )

    def remove_frequency_scaling(self, freqs=None):
        """
        **LLM Docstring**

        Build a non-frequency-scaled version of this mode set by undoing the frequency-scaling rescaling; returns `self` unchanged if not already frequency-scaled.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.freqs`
        :type freqs: np.ndarray | None
        :return: the non-frequency-scaled mode set
        :rtype: MixtureModes
        """
        if not self.frequency_scaled: return self
        freqs, conv = self._frequency_scaling(freqs=freqs)
        L = self.matrix / conv[np.newaxis, :]
        Linv = self.inverse * conv[:, np.newaxis] # del_Q X

        return self.modify(L,
                           inverse=Linv,
                           freqs=freqs,
                           frequency_scaled=False
                           )

    def make_dimensionless(self, freqs=None, masses=None):
        """
        **LLM Docstring**

        Build a fully dimensionless version of this mode set: mass-weighted and frequency-scaled.

        :param freqs: the frequencies to compute the frequency scaling from
        :type freqs: np.ndarray | None
        :param masses: the masses to compute the mass-weighting from
        :type masses: np.ndarray | None
        :return: the dimensionless mode set
        :rtype: MixtureModes
        """
        new = self
        if not new.mass_weighted: new = new.make_mass_weighted(masses=masses)
        if not new.frequency_scaled: new = new.make_frequency_scaled(freqs=freqs)

        return new

    def make_dimensioned(self, freqs=None, masses=None):
        """
        **LLM Docstring**

        Build a fully dimensioned version of this mode set: not mass-weighted and not frequency-scaled.

        :param freqs: the frequencies to compute the frequency-descaling from
        :type freqs: np.ndarray | None
        :param masses: the masses to compute the mass-unweighting from
        :type masses: np.ndarray | None
        :return: the dimensioned mode set
        :rtype: MixtureModes
        """
        new = self
        if not new.mass_weighted: new = new.remove_mass_weighting(masses=masses)
        if not new.frequency_scaled: new = new.remove_frequency_scaling(freqs=freqs)
        return new

    # def take_submodes(self, pos, axes=(0, 1, 2)):
    #     if self.is_cartesian:
    #         pos = np.asanyarray(pos)
    #         axes = np.asanyarray(axes)
    #         pos = (axes[np.newaxis, :] + pos[:, np.newaxis]*len(axes)).flatten()
    #     atom_disps = self.matrix[pos]
    #     f_base_sub = atom_disps.T @ np.diag(self.freqs**2) @ atom_disps
    #     return type(self).from_fg(
    #         self.basis.take_subbasis(pos),
    #         f_base_sub,
    #         self.g_matrix[np.ix_(pos, pos)],
    #         ...
    #     )

    def apply_projection(self, proj, project_transrot=True, masses=None, origin=None):
        """
        **LLM Docstring**

        Apply an arbitrary projection matrix to the mode/inverse matrices (optionally first combining it with a translation/rotation projector), returning a new mode set restricted to (or excluding) the projected subspace.

        :param proj: the projection matrix to apply
        :type proj: np.ndarray
        :param project_transrot: whether to additionally combine `proj` with a translation/rotation projector before applying
        :type project_transrot: bool
        :param masses: masses to use for the translation/rotation projector and to store on the result
        :type masses: np.ndarray | None
        :param origin: reference geometry to use for the translation/rotation projector and to store on the result
        :type origin: np.ndarray | None
        :return: the projected mode set
        :rtype: MixtureModes
        """
        if project_transrot:
            if masses is None:
                m = self.masses
            else:
                m = masses
            if origin is None:
                o = self.origin
            else:
                o = origin

            tr_proj = nput.translation_rotation_projector(
                np.asanyarray(o).reshape((-1, 3)),
                masses=m,
                mass_weighted=self.mass_weighted
            )
            proj = tr_proj @ proj @ tr_proj



        mbc = proj @ self.modes_by_coords
        cbm = self.coords_by_modes @ proj

        return self.modify(
            matrix=mbc,
            inverse=cbm,
            masses=masses,
            origin=origin
        )

    @classmethod
    def _atom_projector(cls, n, i, orthogonal_projection=False):
        """
        **LLM Docstring**

        Build a diagonal `3n x 3n` projection matrix selecting (or, if `orthogonal_projection`, excluding) the Cartesian coordinates belonging to a given set of atom indices.

        :param n: the total number of atoms
        :type n: int
        :param i: the atom index (or indices) to select
        :type i: int | Iterable[int]
        :param orthogonal_projection: whether to build the complementary projector (excluding `i` instead of selecting it)
        :type orthogonal_projection: bool
        :return: the diagonal atom-selection projection matrix
        :rtype: np.ndarray
        """
        if nput.is_numeric(i):
            i = [i]
        z = np.zeros((3 * n, 3 * n))
        for i in i:
            x = np.arange(3 * i, 3 * (i + 1))
            z[x, x] = 1
        if orthogonal_projection:
            z = np.eye(3 * n) - z
        return z
    def apply_constraints(self,
                          coordinate_constraints,
                          atoms=None,
                          masses=None,
                          origin=None,
                          orthogonal_projection=True,
                          ):
        """
        **LLM Docstring**

        Build a new mode set with the given internal coordinates constrained (projected out or onto, depending on `orthogonal_projection`), optionally restricted to a subset of atoms, via `apply_projection`.

        :param coordinate_constraints: the internal coordinate specification(s) to constrain
        :type coordinate_constraints: object
        :param atoms: restrict the constraint projector(s) to a subset of atoms
        :type atoms: Iterable[int] | None
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry to use instead of `self.origin` (un-mass-weighted)
        :type origin: np.ndarray | None
        :param orthogonal_projection: whether the constraint bases define orthogonal (rather than oblique) projections, and whether multiple constraints are combined by sequential orthogonal projection (`True`) or simple summation (`False`)
        :type orthogonal_projection: bool
        :return: the constrained mode set
        :rtype: MixtureModes
        """

        if nput.is_numeric(coordinate_constraints[0]):
            coordinate_constraints = [coordinate_constraints]

        if masses is not None:
            m = masses
        else:
            m = self.masses
        if origin is not None:
            o = origin
        else:
            o = self.remove_mass_weighting().origin
        basis, _, _ = nput.internal_basis(
            np.asanyarray(o).reshape((-1, 3)),
            coordinate_constraints,
            masses=m,
            project_transrot=False
            )

        gi12 = np.diag(np.repeat(1 / np.sqrt(m), 3))
        if self.mass_weighted:
            projections = [
                nput.projection_matrix(gi12 @ b)
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(gi12 @ b)
                for b in basis
            ]
        else:
            projections = [
                nput.projection_matrix(b)
                    if not orthogonal_projection else
                nput.orthogonal_projection_matrix(b)
                for b in basis
            ]

        if atoms is not None:
            if nput.is_numeric(atoms):
                atoms = [atoms]
            nats = len(self.masses)
            a_proj = self._atom_projector(nats, atoms)
            projections = [
                a_proj @ proj @ a_proj
                for proj in projections
            ]

        if orthogonal_projection:
            proj = projections[0]
            for p in projections[1:]:
                proj = p @ proj @ p
        else:
            proj = np.sum(projections, axis=0)

        new = self.apply_projection(proj, masses=m, origin=o, project_transrot=False)
        # if self.mass_weighted:
        #     new = new.make_mass_weighted(masses=m)

        return new

    def apply_transformation(self, tf, **opts):
        """
        **LLM Docstring**

        Apply an arbitrary linear transformation to this mode set, returning the result as a `LocalizedModes` object.

        :param tf: the transformation to apply, in any form accepted by `LocalizedModes`
        :type tf: np.ndarray | tuple
        :param opts: extra options forwarded to the `LocalizedModes` constructor
        :type opts: dict
        :return: the transformed modes
        :rtype: LocalizedModes
        """
        from .LocalizedModes import LocalizedModes

        return LocalizedModes(self, tf, **opts)

    def make_oblique(self):
        """
        **LLM Docstring**

        Build an "oblique" mode representation (see `ObliqueModeGenerator`) that makes the mode transformation as close to orthogonal as possible while still diagonalizing the mode-basis Hessian/G-matrix.

        :return: the oblique-transformed modes
        :rtype: LocalizedModes
        """
        from .ObliqueModes import ObliqueModeGenerator

        f = self.compute_hessian()
        g = self.compute_gmatrix()

        _, _, u, ui = ObliqueModeGenerator(f, g).run()
        return self.apply_transformation((u, ui))

