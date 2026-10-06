"""Choose one Z-matrix for interpolating atom-mapped reaction endpoints.

The discrete search delegates graph internals and constrained Z-matrix assembly
to Molecule. Formed/broken bonds are mandatory. Other coordinates are ranked by
normalized endpoint changes and tried in small subsets. Every candidate is
scored in one fixed Cartesian frame, without per-image alignment.
"""

from dataclasses import dataclass
import itertools

import numpy as np
import McUtils.Numputils as nput
import McUtils.Coordinerds as coordops
from McUtils.Data import AtomData, UnitsData

__all__ = [
    'ReactionCoordinateSearch', 'ReactionCoordinateSearchResult',
    'ReactionCoordinateCandidate', 'ReactionCoordinateChange',
    'ZMatrixReactionInterpolator'
]


@dataclass(frozen=True)
class ReactionCoordinateChange:
    coordinate: tuple
    reactant_value: float
    product_value: float
    change: float
    normalized_change: float


@dataclass
class ReactionCoordinateCandidate:
    zmatrix: np.ndarray
    required_coordinates: tuple
    extra_coordinates: tuple
    images: np.ndarray
    displacement: float
    steric_penalty: float
    additional_penalty: float
    score: float
    source: str
    root: object
    interpolator: object


@dataclass
class ReactionCoordinateSearchResult:
    best: ReactionCoordinateCandidate
    candidates: tuple
    changes: tuple
    mandatory_coordinates: tuple
    rejected: tuple
    parameters: np.ndarray
    budget_exhausted: bool
    reactant: object = None
    product: object = None

    @property
    def zmatrix(self):
        return self.best.zmatrix

    @property
    def images(self):
        return self.best.images

    def interpolate(self, parameters):
        return self.best.interpolator(parameters)

    def get_optimizer(self, energy_evaluator, **options):
        """Connect the selected interpolation to McUtils' native optimizers."""
        from .ReactionPathOptimizer import ReactionPathOptimizer
        return ReactionPathOptimizer.from_coordinate_search(self, energy_evaluator, **options)

    def get_profile_generator(self, energy_evaluator, *, method='neb', **options):
        """Return a ProfileGenerator-compatible optimizer for this initial path."""
        return self.get_optimizer(energy_evaluator, method=method, **options)


class ZMatrixReactionInterpolator:
    """A single ordering and external frame for the whole path, in Bohr/radians.

    Three fixed reference points encode the six external pose coordinates of
    the first three atoms; they are removed from the returned Cartesians. They
    allow both endpoints to retain their supplied embedding. No image is fitted
    or rotated after conversion. Molecular bond/angle/dihedral references all
    come from the supplied Z-matrix.

    Torsions follow their shortest periodic arc. A custom ``interpolation`` is
    called as ``f(parameters, start, end)`` with the end torsions unwrapped;
    it must return shape ``(len(parameters), natoms + 2, 3)`` and preserve the
    endpoints and reference points. Linear interpolation is the default.

    Distances, bend angles, reference planes, and Cartesian round trips are
    checked on evaluation. Continuous nonsingularity is certified for the
    default linear path by conservative interval bounds, or a candidate is
    rejected. Custom paths are checked only at their evaluated images.
    """

    def __init__(self, reactant, product, zmatrix, *, origin=None, axes=None,
                 interpolation=None, singularity_tolerance=1.e-7,
                 roundtrip_tolerance=1.e-7, certificate_depth=12):
        self.endpoints = np.array([reactant, product], dtype=float)
        if (self.endpoints.ndim != 3 or self.endpoints.shape[-1] != 3
                or self.endpoints.shape[1] == 0 or not np.all(np.isfinite(self.endpoints))):
            raise ValueError('Endpoints must be finite and have matching (natoms, 3) shapes')
        self.natoms = self.endpoints.shape[1]
        ordering_array = np.asarray(zmatrix)
        if not np.issubdtype(ordering_array.dtype, np.integer):
            raise ValueError('Z-matrix indices must be integers')
        self.zmatrix = ordering_array.astype(int, copy=True)
        if self.zmatrix.shape != (self.natoms, 4):
            raise ValueError('A full (natoms, 4) Z-matrix is required')
        if sorted(self.zmatrix[:, 0].tolist()) != list(range(self.natoms)):
            raise ValueError('Z-matrix must contain every atom exactly once')
        self.singularity_tolerance = float(singularity_tolerance)
        self.roundtrip_tolerance = float(roundtrip_tolerance)
        if (not np.isfinite(self.singularity_tolerance) or self.singularity_tolerance <= 0
                or not np.isfinite(self.roundtrip_tolerance) or self.roundtrip_tolerance <= 0):
            raise ValueError('Conversion tolerances must be finite and positive')
        self.certificate_depth = ReactionCoordinateSearch._integer(certificate_depth, 'certificate_depth')
        self.interpolation = interpolation
        self.certified = False
        self.certificate_intervals = 0
        default_origin, default_axes = self._frame(self.endpoints)
        self.origin = np.array(default_origin if origin is None else origin, dtype=float)
        self.axes = np.array(default_axes if axes is None else axes, dtype=float)
        if (self.origin.shape != (3,) or self.axes.shape != (2, 3)
                or not np.all(np.isfinite(self.origin)) or not np.all(np.isfinite(self.axes))):
            raise ValueError('Embedding requires a finite origin (3,) and axes (2, 3)')
        if not np.allclose(self.axes @ self.axes.T, np.eye(2), atol=1.e-10):
            raise ValueError('Embedding axes must be orthonormal')
        # The first physical atoms use external references for their pose only.
        ordering = self.zmatrix.copy()
        ordering[0, 1:] = [-3, -1, -2]
        if self.natoms > 1:
            ordering[1, 2:] = [-1, -2]
        if self.natoms > 2:
            ordering[2, 3] = -2
        self.ordering = np.concatenate([
            [[0, -1, -2, -3], [1, 0, -1, -2], [2, 0, 1, -1]], ordering + 3
        ])
        coordops.validate_zmatrix(self.ordering, raise_exception=True)
        self.references = np.array([self.origin, self.origin + self.axes[0],
                                    self.origin + self.axes[1]])
        carts = np.concatenate([np.broadcast_to(self.references, (2, 3, 3)),
                                self.endpoints], axis=1)
        values = coordops.CoordinateSet(carts, coordops.CartesianCoordinates3D).convert(
            coordops.ZMatrixCoordinates, ordering=self.ordering
        )
        self.start, self.end = np.asarray(values).copy()
        self.end[:, 2] = self.start[:, 2] + self._periodic(self.end[:, 2] - self.start[:, 2])
        self.options = dict(ordering=self.ordering, origins=self.origin, axes=self.axes)
        check = self(np.array([0., 1.]))
        if not np.allclose(check, self.endpoints, rtol=0, atol=self.roundtrip_tolerance):
            raise ValueError('Z-matrix endpoint conversion does not preserve the embedding')
        # General reference triples can become collinear between sampled images.
        # Certification is deliberately conservative rather than silently accepting them.
        if interpolation is None:
            self._certify_linear_path()

    @staticmethod
    def _periodic(delta):
        return (delta + np.pi) % (2 * np.pi) - np.pi

    @staticmethod
    def _frame(endpoints):
        coords = endpoints[0]
        center = coords.mean(axis=0)
        disp = coords - center
        a = disp[np.argmax(np.linalg.norm(disp, axis=1))]
        norm = np.linalg.norm(a)
        a = a / norm if norm > 1.e-10 else np.array([1., 0., 0.])
        perpendicular = disp - np.outer(disp @ a, a)
        b = perpendicular[np.argmax(np.linalg.norm(perpendicular, axis=1))]
        norm = np.linalg.norm(b)
        if norm <= 1.e-10:
            b = np.eye(3)[np.argmin(np.abs(a))]
            b -= np.dot(b, a) * a
            norm = np.linalg.norm(b)
        b /= norm
        c = np.cross(a, b)
        span = max(float(np.max(np.linalg.norm(endpoints - center, axis=-1))), 1.)
        return center - span * (2 * a + 3 * b + 5 * c), np.array([a, b])

    def _values(self, parameters):
        if self.interpolation is None:
            return self.start + parameters[:, None, None] * (self.end - self.start)
        values = np.asarray(self.interpolation(parameters.copy(), self.start.copy(), self.end.copy()), dtype=float)
        if values.shape != (len(parameters), self.natoms + 2, 3):
            raise ValueError('Interpolation returned an invalid internal-coordinate shape')
        return values

    def _convert(self, values):
        return np.asarray(coordops.CoordinateSet(values, coordops.ZMatrixCoordinates,
                                                 self.options).convert(coordops.CartesianCoordinates3D))

    def _validate(self, values, carts):
        tol = self.singularity_tolerance
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(carts)):
            raise ValueError('Nonfinite interpolated coordinates')
        if np.any(values[..., 0] <= tol):
            raise ValueError('Zero Z-matrix distance')
        bends = values[:, 1:, 1]
        if np.any(bends <= tol) or np.any(bends >= np.pi - tol):
            raise ValueError('Singular Z-matrix bend')
        if not np.allclose(carts[:, :3], self.references, rtol=0, atol=self.roundtrip_tolerance):
            raise ValueError('Interpolation moved the external reference frame')
        for row in self.ordering[3:]:
            j, k, l = row[1:]
            a, b = carts[:, j] - carts[:, k], carts[:, l] - carts[:, k]
            sine = np.linalg.norm(np.cross(a, b), axis=-1) / np.maximum(
                np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1), 1.e-30)
            if np.any(sine <= tol):
                raise ValueError('Collinear Z-matrix reference plane')
        back = np.asarray(coordops.CoordinateSet(carts, coordops.CartesianCoordinates3D).convert(
            coordops.ZMatrixCoordinates, ordering=self.ordering))
        error = back - values
        error[..., 2] = self._periodic(error[..., 2])
        if np.max(np.abs(error)) > self.roundtrip_tolerance:
            raise ValueError('Interpolated Z-matrix does not round-trip')

    def __call__(self, parameters):
        parameters = np.asarray(parameters, dtype=float)
        scalar = parameters.ndim == 0
        parameters = np.atleast_1d(parameters)
        if parameters.ndim != 1 or len(parameters) == 0 or not np.all(np.isfinite(parameters)) or np.any((parameters < 0) | (parameters > 1)):
            raise ValueError('Interpolation parameters must lie in [0, 1]')
        values = self._values(parameters)
        carts = self._convert(values)
        self._validate(values, carts)
        for index, target in ((0., self.endpoints[0]), (1., self.endpoints[1])):
            if not np.allclose(carts[parameters == index, 3:], target, rtol=0,
                               atol=self.roundtrip_tolerance):
                raise ValueError('Interpolation does not preserve the endpoints')
        return carts[0, 3:] if scalar else carts[:, 3:]

    def _certify_linear_path(self):
        """Enclose the forward construction on adaptively subdivided t intervals.

        Outward-rounded boxes bound every reference vector and its cross
        product. Strictly positive lower bounds prove that no normalization
        or reference-plane construction becomes singular anywhere in [0, 1].
        An unresolved interval is a rejection, even if its sampled images work.
        """
        def box(lo, hi):
            return np.array([np.nextafter(lo, -np.inf), np.nextafter(hi, np.inf)])

        def add(a, b):
            return box(a[0] + b[0], a[1] + b[1])

        def negate(a):
            return np.array([-a[1], -a[0]])

        def subtract(a, b):
            return add(a, negate(b))

        def multiply(a, b):
            products = np.array([a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]])
            return box(np.min(products, axis=0), np.max(products, axis=0))

        def norm(a):
            near = np.where((a[0] <= 0) & (a[1] >= 0), 0., np.minimum(np.abs(a[0]), np.abs(a[1])))
            far = np.maximum(np.abs(a[0]), np.abs(a[1]))
            squares = multiply(box(near, far), box(near, far))
            squared = np.array([0., 0.])
            for component in squares.T:
                squared = add(squared, component)
            return box(np.sqrt(max(0., squared[0])), np.sqrt(max(0., squared[1])))

        def normalize(a, length):
            if length[0] <= self.singularity_tolerance:
                raise ValueError('Unresolved interval vector normalization')
            reciprocal = box(1 / length[1], 1 / length[0])
            return multiply(a, reciprocal[:, None])

        def cross(a, b):
            indices, other = [1, 2, 0], [2, 0, 1]
            return subtract(multiply(a[:, indices], b[:, other]), multiply(a[:, other], b[:, indices]))

        def trig(interval, cosine=False):
            # Include every extremum, not just the interval's endpoints.
            lo, hi = interval
            if cosine:
                lo, hi = lo + np.pi / 2, hi + np.pi / 2
                lo, hi = np.nextafter(lo, -np.inf), np.nextafter(hi, np.inf)
            values = [np.sin(lo), np.sin(hi)]
            first = int(np.ceil((lo - np.pi / 2) / np.pi))
            last = int(np.floor((hi - np.pi / 2) / np.pi))
            values.extend(1. if k % 2 == 0 else -1. for k in range(first, last + 1))
            return box(min(values), max(values))

        def enclose(left, right):
            delta = subtract(box(self.end, self.end), box(self.start, self.start))
            q = add(box(self.start, self.start), multiply(box(left, right)[:, None, None], delta))
            points = {i: box(p, p) for i, p in enumerate(self.references)}
            for index, row in enumerate(self.ordering[3:]):
                atom, j, k, l = row
                v, u = subtract(points[k], points[j]), subtract(points[l], points[k])
                nv, nu = norm(v), norm(u)
                perpendicular = cross(v, u)
                nc = norm(perpendicular)
                if (nv[0] <= self.singularity_tolerance or nu[0] <= self.singularity_tolerance
                        or nc[0] <= self.singularity_tolerance * nv[1] * nu[1]):
                    raise ValueError('Unresolved interval reference plane')
                e1, e3 = normalize(v, nv), normalize(perpendicular, nc)
                e2 = cross(e3, e1)
                radius, bend, torsion = (q[:, index + 2, component] for component in range(3))
                if (radius[0] <= self.singularity_tolerance or bend[0] <= self.singularity_tolerance
                        or bend[1] >= np.pi - self.singularity_tolerance):
                    raise ValueError('Singular interval internal coordinate')
                direction = add(multiply(trig(bend, cosine=True)[:, None], e1),
                                multiply(trig(bend)[:, None],
                                         subtract(multiply(trig(torsion, cosine=True)[:, None], e2),
                                                  multiply(trig(torsion)[:, None], e3))))
                points[atom] = add(points[j], multiply(radius[:, None], direction))
                if not np.all(np.isfinite(points[atom])):
                    raise ValueError('Nonfinite interval enclosure')

        pending = [(0., 1., 0)]
        intervals = 0
        while pending:
            left, right, depth = pending.pop()
            try:
                enclose(left, right)
                intervals += 1
            except ValueError as error:
                if depth >= self.certificate_depth:
                    raise ValueError('Could not certify a nonsingular Z-matrix on the continuous path') from error
                middle = (left + right) / 2
                pending.extend([(left, middle, depth + 1), (middle, right, depth + 1)])
        self.certified = True
        self.certificate_intervals = intervals


class ReactionCoordinateSearch:
    """Bounded sparse search over constrained bond-graph Z-matrices.

    Inputs are aligned Molecules with identical atom ordering, in BohrRadius.
    ``from_alignment`` accepts the Molecule-valued ConformerAlignmentResult.
    Arrays are not accepted here because Molecule supplies the topology toolkit.

    Formed/broken bonds and user ``required_coordinates`` are mandatory. The
    candidate pool combines both endpoints' unpruned bond-graph internals and
    internals of their union augmented with all reactive-atom pairs. Virtual
    nearest-fragment links and endpoint Z-matrix coordinates supply scalar
    fragment coordinates instead of orientation dictionaries. Stretches are ranked by abs(delta)/length_scale;
    bends/torsions by abs(delta)/angle_scale, with periodic torsion differences.

    Beam search tries at most ``max_extra_coordinates`` from the highest
    ``max_changing_coordinates`` eligible changes. Both endpoints and their
    union supply assembly seeds, with a connected union seed for spectators.
    ``roots`` defaults to None and up to three
    reactive atom roots. Required coordinates are checked after assembly.
    The metric is displacement + steric_weight*steric_penalty + image_penalty.
    Near ties within ``score_tolerance`` favor fewer supplemental constraints.
    This is a bounded heuristic, not a claim of global optimality.

    ``pair_metric(left, right)`` can replace Cartesian RMSD, e.g. with a geodesic
    distance. All coordinate and metric lengths use Bohr. ``interpolation`` is
    forwarded to ZMatrixReactionInterpolator for a parametric internal path.
    The steric term integrates sum(exp(decay*(1-distance/contact))) over the
    path parameter. Contacts are scaled sums of covalent radii. Pairs bonded
    in either endpoint are excluded; ``steric_exclusion_depth`` also permits
    excluding their graph neighbors. ``image_penalty(images)`` is an optional
    per-image contribution in metric units, available for a future attraction.
    No Lennard-Jones attraction is enabled by default.
    """

    def __init__(self, react, prod, *, required_coordinates=(), reactive_atoms=None,
                 candidate_coordinates=(), parameters=None, num_images=21,
                 length_scale=None, angle_scale=1., change_threshold=1.e-3,
                 max_changing_coordinates=12, max_extra_coordinates=3,
                 beam_width=4, max_candidates=128, max_assemblies=2048,
                 roots=None, score_tolerance=None, pair_metric=None,
                 atom_weights=None, interpolation=None, image_penalty=None,
                 steric_weight=None, steric_decay=8., steric_contact_scaling=1.,
                 steric_exclusion_depth=1, covalent_radii=None,
                 embedding_origin=None, embedding_axes=None,
                 singularity_tolerance=1.e-7, roundtrip_tolerance=1.e-7,
                 certificate_depth=12):
        if any(not all(hasattr(m, name) for name in ('coords', 'atoms', 'get_bond_zmatrix', 'get_bond_graph_internals'))
               for m in (react, prod)):
            raise TypeError('Coordinate search requires two Molecule conformers')
        self.react, self.prod = react, prod
        self.endpoints = np.array([react.coords, prod.coords], dtype=float)
        if (self.endpoints.ndim != 3 or self.endpoints.shape[2] != 3
                or self.endpoints.shape[1] == 0 or not np.all(np.isfinite(self.endpoints))):
            raise ValueError('Endpoint coordinates must be finite and have the same (natoms, 3) shape')
        self.natoms = self.endpoints.shape[1]
        identities = [tuple((AtomData[a, 'Number'], AtomData[a, 'MassNumber']) for a in m.atoms)
                      for m in (react, prod)]
        if identities[0] != identities[1]:
            raise ValueError('Endpoint atoms must agree index by index')
        self.bonds = tuple({coordops.canonicalize_internal(b[:2]): b[2] if len(b) > 2 else 1
                            for b in m.bonds} for m in (react, prod))
        added, removed = react.edge_graph.graph_difference(prod.edge_graph)
        self.formed_bonds = tuple(sorted((int(i), int(j)) for i, j in added if i < j))
        self.broken_bonds = tuple(sorted((int(i), int(j)) for i, j in removed if i < j))
        self.mandatory = tuple(sorted(set(self.formed_bonds + self.broken_bonds)
                                      | set(self._coordinates(required_coordinates)), key=lambda c: (len(c), c)))
        changed = set(self.mandatory) | {b for b in self.bonds[0].keys() & self.bonds[1].keys()
                                        if self.bonds[0][b] != self.bonds[1][b]}
        self.reactive_atoms = sorted(set(a for c in changed for a in c)
                                     | set(() if reactive_atoms is None else reactive_atoms))
        if any(isinstance(a, (bool, np.bool_)) or not isinstance(a, (int, np.integer))
               or not 0 <= a < self.natoms for a in self.reactive_atoms):
            raise ValueError('Reactive atoms must be valid atom indices')
        union_bonds = set(self.bonds[0]) | set(self.bonds[1])
        self.union = react.modify(bonds=sorted(union_bonds))
        # Virtual links make scalar inter-fragment coordinates available even
        # for spectator complexes and small fragment-assembly toolkit gaps.
        connected_bonds = set(union_bonds)
        fragments = self.union.edge_graph.get_fragments()
        for fragment in fragments[1:]:
            left, right = np.meshgrid(fragments[0], fragment, indexing='ij')
            lengths = np.mean(np.linalg.norm(self.endpoints[:, left] - self.endpoints[:, right], axis=-1), axis=0)
            index = np.unravel_index(np.argmin(lengths), lengths.shape)
            connected_bonds.add(coordops.canonicalize_internal((int(left[index]), int(right[index]))))
        connected = react.modify(bonds=sorted(connected_bonds))
        self.sources = (('reactant', react), ('product', prod), ('union', self.union))
        if connected_bonds != union_bonds:
            self.sources += (('connected_union', connected),)
        self.augmented = react.modify(bonds=sorted(connected_bonds | set(itertools.combinations(self.reactive_atoms, 2))))
        self.supplied_candidates = self._coordinates(candidate_coordinates)
        self.parameters = np.linspace(0., 1., self._integer(num_images, 'num_images', minimum=2)) if parameters is None else np.asarray(parameters, dtype=float).copy()
        if (self.parameters.ndim != 1 or len(self.parameters) < 2
                or not np.all(np.isfinite(self.parameters)) or self.parameters[0] != 0 or self.parameters[-1] != 1
                or np.any(np.diff(self.parameters) <= 0)):
            raise ValueError('Parameters must increase strictly from 0 to 1')
        angstrom = UnitsData.convert('Angstroms', 'BohrRadius')
        self.length_scale = self._positive(angstrom if length_scale is None else length_scale, 'length_scale')
        self.angle_scale = self._positive(angle_scale, 'angle_scale')
        self.change_threshold = self._positive(change_threshold, 'change_threshold', zero=True)
        self.max_changing = self._integer(max_changing_coordinates, 'max_changing_coordinates', minimum=0)
        self.max_extra = self._integer(max_extra_coordinates, 'max_extra_coordinates', minimum=0)
        self.beam_width = self._integer(beam_width, 'beam_width')
        self.max_candidates = self._integer(max_candidates, 'max_candidates')
        self.max_assemblies = self._integer(max_assemblies, 'max_assemblies')
        self.roots = tuple([None] + self.reactive_atoms[:3] if roots is None else roots)
        if not self.roots or any(r is not None and (isinstance(r, bool) or not isinstance(r, (int, np.integer)) or not 0 <= r < self.natoms) for r in self.roots):
            raise ValueError('Roots must be atom indices or None')
        self.score_tolerance = self._positive(1.e-4 * angstrom if score_tolerance is None else score_tolerance, 'score_tolerance', zero=True)
        self.pair_metric = pair_metric
        self.weights = np.ones(self.natoms) if atom_weights is None else np.asarray(atom_weights, dtype=float)
        if self.weights.shape != (self.natoms,) or not np.all(np.isfinite(self.weights)) or np.any(self.weights <= 0):
            raise ValueError('Atom weights must be positive and finite')
        self.image_penalty = image_penalty
        self.steric_weight = self._positive(angstrom if steric_weight is None else steric_weight, 'steric_weight', zero=True)
        self.steric_decay = self._positive(steric_decay, 'steric_decay')
        self.steric_scaling = self._positive(steric_contact_scaling, 'steric_contact_scaling')
        if self.steric_decay > 100:
            raise ValueError('steric_decay must not exceed 100 to keep contact penalties finite')
        self.radii = (np.array([AtomData[a, 'CovalentRadius'] for a in react.atoms]) * angstrom
                      if covalent_radii is None else np.asarray(covalent_radii, dtype=float))
        if self.radii.shape != (self.natoms,) or not np.all(np.isfinite(self.radii)) or np.any(self.radii <= 0):
            raise ValueError('Covalent radii must be finite positive Bohr values per atom')
        depth = self._integer(steric_exclusion_depth, 'steric_exclusion_depth', minimum=0)
        distances = self.union.edge_graph.get_distances()
        self.steric_pairs = np.array([(i, j) for i, j in itertools.combinations(range(self.natoms), 2)
                                      if distances[i, j] > depth], dtype=int).reshape(-1, 2)
        self.interpolator_options = dict(origin=embedding_origin, axes=embedding_axes,
                                         interpolation=interpolation,
                                         singularity_tolerance=self._positive(singularity_tolerance, 'singularity_tolerance'),
                                         roundtrip_tolerance=self._positive(roundtrip_tolerance, 'roundtrip_tolerance'),
                                         certificate_depth=self._integer(certificate_depth, 'certificate_depth'))

    @classmethod
    def from_alignment(cls, result, **options):
        return cls(result.reactant, result.product, **options)

    @staticmethod
    def _positive(value, name, zero=False):
        value = float(value)
        if not np.isfinite(value) or (value < 0 if zero else value <= 0):
            raise ValueError(f'{name} must be finite and {"nonnegative" if zero else "positive"}')
        return value

    @staticmethod
    def _integer(value, name, minimum=1):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
            raise ValueError(f'{name} must be an integer >= {minimum}')
        return int(value)

    def _coordinates(self, specs):
        result = []
        for spec in specs:
            spec = tuple(spec)
            if len(spec) not in (2, 3, 4) or any(isinstance(a, (bool, np.bool_)) or not isinstance(a, (int, np.integer)) or not 0 <= a < self.natoms for a in spec):
                raise ValueError('Coordinates must contain two, three, or four valid atom indices')
            spec = coordops.canonicalize_internal(spec)
            if spec is None:
                raise ValueError('Coordinate indices must be distinct')
            if spec not in result:
                result.append(spec)
        return tuple(result)

    def coordinate_changes(self):
        pool = set(self.mandatory) | set(self.supplied_candidates)
        for molecule in (self.react, self.prod, self.augmented):
            pool.update(self._coordinates(molecule.get_bond_graph_internals(include_fragments=False, pruning=False)))
        for _, molecule in self.sources[:2]:
            try:
                pool.update(coordops.extract_zmatrix_internals(molecule.get_bond_zmatrix()))
            except (ValueError, IndexError, KeyError):
                # Small/light-atom graphs are not supported by every backbone builder.
                pass
        changes = []
        for spec in sorted(pool, key=lambda c: (len(c), c)):
            points = [self.endpoints[:, a] for a in spec]
            if len(spec) == 2:
                values = nput.pts_norms(*points)
            elif len(spec) == 3:
                values = nput.pts_angles(*points)[0]
            else:
                values = nput.pts_dihedrals(*points)
            delta = values[1] - values[0]
            if len(spec) == 4:
                delta = ZMatrixReactionInterpolator._periodic(delta)
            scale = self.length_scale if len(spec) == 2 else self.angle_scale
            if np.all(np.isfinite(values)):
                changes.append(ReactionCoordinateChange(spec, float(values[0]), float(values[1]),
                                                         float(abs(delta)), float(abs(delta) / scale)))
        return tuple(sorted(changes, key=lambda c: (-c.normalized_change, len(c.coordinate), c.coordinate)))

    def score_images(self, images):
        images = np.asarray(images, dtype=float)
        if images.shape != (len(self.parameters), self.natoms, 3) or not np.all(np.isfinite(images)):
            raise ValueError('Images must be finite and match the sampling grid')
        if self.pair_metric is None:
            steps = np.sqrt(np.sum(self.weights * np.sum(np.diff(images, axis=0)**2, axis=-1), axis=-1) / np.sum(self.weights))
        else:
            steps = np.array([self.pair_metric(left.copy(), right.copy()) for left, right in zip(images[:-1], images[1:])], dtype=float)
        if steps.shape != (len(images) - 1,) or not np.all(np.isfinite(steps)) or np.any(steps < 0):
            raise ValueError('Pair metric must return finite nonnegative scalar distances')
        steric = np.zeros(len(images))
        if self.steric_weight and len(self.steric_pairs):
            i, j = self.steric_pairs.T
            distances = np.linalg.norm(images[:, i] - images[:, j], axis=-1)
            contacts = self.steric_scaling * (self.radii[i] + self.radii[j])
            steric = np.sum(np.exp(self.steric_decay * (1 - distances / contacts)), axis=-1)
        extra = np.zeros(len(images)) if self.image_penalty is None else np.asarray(self.image_penalty(images.copy()), dtype=float)
        if extra.shape != (len(images),) or not np.all(np.isfinite(extra)):
            raise ValueError('Image penalty must return a finite scalar per image')
        # Integrals preserve the weight when the sampling grid is refined.
        dt = np.diff(self.parameters)
        steric_integral = float(np.sum(dt * (steric[:-1] + steric[1:]) / 2))
        additional = float(np.sum(dt * (extra[:-1] + extra[1:]) / 2))
        displacement = float(np.sum(steps))
        score = displacement + self.steric_weight * steric_integral + additional
        if not np.isfinite(score):
            raise ValueError('Nonfinite path score')
        return displacement, steric_integral, additional, score

    def _prefer(self, candidate, incumbent):
        if incumbent is None or candidate.score < incumbent.score - self.score_tolerance:
            return True
        return (abs(candidate.score - incumbent.score) <= self.score_tolerance
                and len(candidate.extra_coordinates) < len(incumbent.extra_coordinates))

    def search(self):
        changes = self.coordinate_changes()
        optional = tuple(c.coordinate for c in changes if c.coordinate not in self.mandatory
                         and c.normalized_change >= self.change_threshold)[:self.max_changing]
        ranks = {c: i for i, c in enumerate(optional)}
        def subset_rank(subset):
            return tuple(sorted(ranks[c] for c in subset))
        candidates, rejected, seen_matrices, seen_subsets = [], [], {}, set()
        assemblies, exhausted, best = 0, False, None

        def evaluate(extra):
            nonlocal assemblies, exhausted, best
            required = self.mandatory + tuple(extra)
            produced = []
            for name, molecule in self.sources:
                for root in self.roots:
                    if assemblies >= self.max_assemblies or len(candidates) >= self.max_candidates:
                        exhausted = True
                        return produced
                    assemblies += 1
                    try:
                        orderings = [molecule.get_bond_zmatrix(required_coordinates=required or None, root=root)]
                    except (ValueError, IndexError, KeyError) as error:
                        if self.natoms <= 3:
                            # A toolkit chain fallback covers small backbone-builder gaps.
                            orderings = [coordops.chain_zmatrix(p) for p in itertools.permutations(range(self.natoms))]
                        else:
                            rejected.append((required, name, root, str(error)))
                            continue
                    for ordering in orderings:
                        try:
                            ordering = np.asarray(ordering, dtype=int)
                            represented = set(coordops.extract_zmatrix_internals(ordering))
                            missing = set(required) - represented
                            if missing:
                                raise ValueError(f'Missing required coordinates: {sorted(missing)}')
                            key = tuple(ordering.ravel())
                            if key in seen_matrices:
                                cached = seen_matrices[key]
                                if cached is None:
                                    continue
                                candidate = ReactionCoordinateCandidate(ordering.copy(), required, tuple(extra), cached.images,
                                    cached.displacement, cached.steric_penalty, cached.additional_penalty, cached.score,
                                    name, root, cached.interpolator)
                            else:
                                if len(candidates) >= self.max_candidates:
                                    exhausted = True
                                    return produced
                                seen_matrices[key] = None
                                interpolator = ZMatrixReactionInterpolator(*self.endpoints, ordering, **self.interpolator_options)
                                images = interpolator(self.parameters)
                                displacement, steric, additional, score = self.score_images(images)
                                candidate = ReactionCoordinateCandidate(ordering.copy(), required, tuple(extra), images,
                                    displacement, steric, additional, score, name, root, interpolator)
                                seen_matrices[key] = candidate
                                candidates.append(candidate)
                            produced.append(candidate)
                            if self._prefer(candidate, best):
                                best = candidate
                        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
                            rejected.append((required, name, root, str(error)))
            return produced

        beam = [()]
        for depth in range(min(self.max_extra, len(optional)) + 1):
            states = []
            for subset in beam:
                if subset in seen_subsets:
                    continue
                seen_subsets.add(subset)
                produced = evaluate(subset)
                if produced:
                    states.append((min(c.score for c in produced), subset))
                if exhausted:
                    break
            if exhausted or depth == self.max_extra:
                break
            # Failed assemblies can become feasible after another required coordinate.
            parents = [s for _, s in sorted(states, key=lambda x: (x[0], subset_rank(x[1])))[:self.beam_width]] or beam[:self.beam_width]
            beam = sorted({tuple(sorted(set(s) | {c}, key=lambda c: (len(c), c)))
                           for s in parents for c in optional if c not in s}, key=subset_rank)
        # Try removing the least-changing supplements from the winner. The beam
        # may not have visited these subsets on its way to the winning ordering.
        reduced = True
        while best is not None and best.extra_coordinates and reduced and not exhausted:
            reduced = False
            current = best.extra_coordinates
            for remove in sorted(current, key=lambda c: ranks[c], reverse=True):
                subset = tuple(c for c in current if c != remove)
                if subset not in seen_subsets:
                    seen_subsets.add(subset)
                    evaluate(subset)
                if len(best.extra_coordinates) < len(current):
                    reduced = True
                    break
                if exhausted:
                    break
        if best is None:
            error = ValueError(f'No nonsingular Z-matrix satisfied all mandatory coordinates {self.mandatory}')
            error.rejected = tuple(rejected)
            raise error
        return ReactionCoordinateSearchResult(best, tuple(candidates), changes, self.mandatory,
                                               tuple(rejected), self.parameters.copy(), exhausted,
                                               self.react.modify(coords=self.endpoints[0].copy()),
                                               self.prod.modify(coords=self.endpoints[1].copy()))
