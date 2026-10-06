"""Rigid fragment preconditioning of atom-mapped reaction endpoints.

Implements stages S1--S5 of Robertson and Habershon, J. Comput. Chem.
42, 761--770 (2021), https://doi.org/10.1002/jcc.26495. This is an
independent implementation of the published geometric construction.

The appendix's Eq. (7) has an ambiguous repulsion sign and uses ``c_ab``
in the force but ``phi*c_ab`` in its switch. Here the repulsion points
away from the other atom and vanishes at ``phi*c_ab``. Empty sums are
zero; point fragments use direct translational forces when their zero
lever arms would otherwise suppress every force. Their molecular sphere
radius is their covalent radius. Exact exponential-map rotations replace
finite Cartesian steps after rigid projection, preserving all internal
distances. Optional seeded perturbations break symmetric stationary points.
If force stepping exhausts its budget, a local least-squares solve over rigid
fragment poses can polish the projected force residual. The paper prescribes
neither numerical step sizes nor stopping tolerances; both are configurable.

No interpolation or internal-coordinate selection is performed here.
"""

from dataclasses import dataclass
import warnings

import numpy as np
from scipy.optimize import least_squares
import McUtils.Numputils as nput
from McUtils.Data import AtomData, UnitsData
from McUtils.Graphs import EdgeGraph

__all__ = [
    "ConformerAlignment", "ConformerAlignmentResult", "ConformerAlignmentStage"
]


@dataclass
class ConformerAlignmentStage:
    """One stage's coordinates and residual (in internal Angstrom units)."""
    name: str
    coordinates: np.ndarray
    iterations: int = 0
    residual: float = 0.0
    converged: bool = True
    polish_evaluations: int = 0


@dataclass
class ConformerAlignmentResult:
    """Aligned copies, preserving atom order; arrays remain arrays.

    Coordinates, including stage snapshots, have the input coordinate units.
    ``converged`` refers to force relaxation, not a physical energy minimum.
    """
    reactant: object
    product: object
    reactant_fragments: tuple
    product_fragments: tuple
    formed_bonds: tuple
    broken_bonds: tuple
    stages: tuple
    coordinate_units: str

    @property
    def coordinates(self):
        return np.array([
            np.asarray(getattr(self.reactant, 'coords', self.reactant)),
            np.asarray(getattr(self.product, 'coords', self.product))
        ])

    @property
    def converged(self):
        return all(s.converged for s in self.stages)


class ConformerAlignment:
    """Prepare unimolecular or bimolecular reactant/product conformers.

    ``react`` and ``prod`` are Molecule objects, sequences of one or two
    Molecules, or arrays of shape ``(natoms, 3)``. Sequences are concatenated
    in the supplied order. Atom identities MUST match index by index between
    endpoints; no atom mapping, permutation, or reflection is attempted.
    Fragment membership may differ between endpoints (1->2, 2->1, 2->2).

    For arrays, supply ``atoms`` (or ``covalent_radii``) and both bond lists.
    Bond lists contain ``(i, j)`` or ``(i, j, order)`` records; empty lists
    describe unbonded atoms. Molecule bonds are used unless overridden.
    Connected components determine fragments unless explicit partitions are
    supplied. Each endpoint must have one or two fragments. Explicit
    partitions must cover every atom once and cannot cut a supplied bond.

    ``coordinate_units`` defaults to BohrRadius, as in Molecule. Molecule
    coordinates must be in BohrRadius. Supplied radii and all length-valued
    options use coordinate_units. The calculation is normalized internally
    to Angstroms so changing units does not change the result. Defaults for
    max_displacement and endpoint_separation are 0.05 and 2.0 Angstroms.
    endpoint_separation is an outward displacement PER FRAGMENT at the end;
    set it to zero to retain the S5 encounter geometries.

    The stage constants follow Table 1. ``covalent_radius_scaling`` defines
    c_ab = scaling*(radius_a + radius_b)/2 in Eq. (7); contact distances
    are 1.5*c_ab for reactive atom sets and 2*c_ab otherwise. The paper does
    not specify the proportionality constant, so it is exposed here.

    ``max_iterations`` is a positive integer or a dict with S2/S4/S5 keys.
    ``polish_max_evaluations`` separately limits a local rigid-pose
    least-squares polish if S4/S5 force stepping exhausts its budget; zero
    disables this extension. max_displacement/max_rotation bound force
    stepping, not the trial poses used during polishing.
    The convergence tolerance is the maximum projected atom-force norm in
    the internally normalized force field, not an energy gradient tolerance.
    ``alignment_weights`` optionally weights the final whole-image Kabsch
    fit; the default uses equal atom weights. The fragment algorithm uses
    geometric centres regardless of these weights, as in the paper.

    ``random_seed`` and ``perturbation`` control optional early rigid noise;
    the default is deterministic without noise. ``perturbation_steps`` caps
    its duration. ``align`` warns on incomplete force relaxation; strict=True
    raises instead. Stage snapshots and residuals remain available for review.

    Examples
    --------
    >>> result = ConformerAlignment(react, prod).align()
    >>> aligned_react, aligned_prod = result.reactant, result.product
    >>> result = ConformerAlignment([react_a, react_b], prod).align()
    """

    stage_constants = {
        'S2': (1.0, 0.0, 0.0, 0.0, 0.0),
        'S4': (50.0, 1.0, 1.0, 0.0, 0.0),
        'S5': (0.0, 1.0, 3.0, 50.0, 2.0)
    }

    def __init__(self, react, prod, *,
                 atoms=None, reactant_bonds=None, product_bonds=None,
                 reactant_fragments=None, product_fragments=None,
                 covalent_radii=None, coordinate_units='BohrRadius',
                 covalent_radius_scaling=1.0, alignment_weights=None,
                 endpoint_separation=None, max_displacement=None,
                 max_rotation=0.05, step_size=0.5,
                 max_iterations=3000, tolerance=1.e-5,
                 polish_max_evaluations=100,
                 random_seed=0, perturbation=0.0, perturbation_steps=100):
        self.react = self._prepare_conformer(react)
        self.prod = self._prepare_conformer(prod)
        self.to_angstrom = UnitsData.convert(coordinate_units, 'Angstroms')
        self.coordinate_units = coordinate_units
        for obj in (self.react, self.prod):
            if hasattr(obj, 'coords') and coordinate_units != 'BohrRadius':
                raise ValueError('Molecule coordinates must use BohrRadius')
        endpoints = [self._coordinates(self.react), self._coordinates(self.prod)]
        if endpoints[0].shape != endpoints[1].shape:
            raise ValueError('Endpoint coordinates must have the same shape')
        self.initial = np.array(endpoints) * self.to_angstrom
        self.natoms = self.initial.shape[1]

        atom_lists = [getattr(obj, 'atoms', None) for obj in (self.react, self.prod)]
        if atoms is not None:
            atom_lists = [atoms if a is None else a for a in atom_lists]
        identities = [None if a is None else self._atom_identities(a) for a in atom_lists]
        if all(a is not None for a in identities) and identities[0] != identities[1]:
            raise ValueError('Atom identities must agree index by index')
        if atoms is not None and any(
            a is not None and a != self._atom_identities(atoms) for a in identities
        ):
            raise ValueError('Explicit atoms disagree with conformer atoms')
        if any(a is not None and len(a) != self.natoms for a in identities):
            raise ValueError('Atom lists must have natoms entries')
        if covalent_radii is None:
            atom_list = next((a for a in atom_lists if a is not None), None)
            if atom_list is None:
                raise ValueError('Supply atoms or covalent_radii for array inputs')
            self.radii = np.array([AtomData[a, 'CovalentRadius'] for a in atom_list])
        else:
            self.radii = np.asarray(covalent_radii, dtype=float) * self.to_angstrom
        if (self.radii.shape != (self.natoms,)
                or not np.all(np.isfinite(self.radii)) or np.any(self.radii <= 0)):
            raise ValueError('Covalent radii must be finite positive values per atom')

        bond_lists = [reactant_bonds, product_bonds]
        for i, obj in enumerate((self.react, self.prod)):
            if bond_lists[i] is None:
                bond_lists[i] = getattr(obj, 'bonds', None)
            if bond_lists[i] is None:
                raise ValueError('Both endpoint bond lists are required')
        self.bonds = tuple(self._bond_map(b) for b in bond_lists)
        # Connectivity graphs have unit weights: changing a bond order does
        # not create or remove a bond. Keep orders separately for reactive atoms.
        self.graphs = tuple(
            EdgeGraph(np.arange(self.natoms), np.array(list(b), dtype=int).reshape(-1, 2))
            for b in self.bonds
        )
        self.fragments = tuple(
            self._partition(f, graph)
            for f, graph in zip((reactant_fragments, product_fragments), self.graphs)
        )
        self.owners = tuple(self._owners(f) for f in self.fragments)
        added, removed = self.graphs[0].graph_difference(self.graphs[1])
        # Undirected adjacency differences contain both directions of each bond.
        self.formed_bonds = tuple(sorted((int(i), int(j)) for i, j in added if i < j))
        self.broken_bonds = tuple(sorted((int(i), int(j)) for i, j in removed if i < j))
        changed = set(self.formed_bonds) | set(self.broken_bonds) | {
            b for b in self.bonds[0].keys() & self.bonds[1].keys()
            if self.bonds[0][b] != self.bonds[1][b]
        }
        self.reactive_atoms = set(a for b in changed for a in b)
        self.reactive_sets = tuple(self._reactive_sets(i) for i in range(2))
        self.shared_sets = tuple(tuple(
            np.intersect1d(r, p) for p in self.fragments[1]
        ) for r in self.fragments[0])
        self.molecular_radii = tuple(tuple(
            self._molecular_radius(self.initial[i], f) for f in fragments
        ) for i, fragments in enumerate(self.fragments))

        self.radius_scaling = self._positive(covalent_radius_scaling, 'covalent_radius_scaling')
        self.max_displacement = (0.05 if max_displacement is None else
                                 self._positive(max_displacement, 'max_displacement') * self.to_angstrom)
        self.endpoint_separation = (2.0 if endpoint_separation is None else
                                    self._positive(endpoint_separation, 'endpoint_separation', zero=True)
                                    * self.to_angstrom)
        self.max_rotation = self._positive(max_rotation, 'max_rotation')
        self.step_size = self._positive(step_size, 'step_size')
        self.tolerance = self._positive(tolerance, 'tolerance')
        self.polish_max_evaluations = self._integer(
            polish_max_evaluations, 'polish_max_evaluations', zero=True
        )
        self.perturbation = self._positive(perturbation, 'perturbation', zero=True)
        self.perturbation_steps = self._integer(perturbation_steps, 'perturbation_steps', zero=True)
        self.random_seed = random_seed
        if isinstance(max_iterations, dict):
            if set(max_iterations) - {'S2', 'S4', 'S5'}:
                raise ValueError('Iteration limits only accept S2, S4 and S5 keys')
            self.iterations = {s: self._integer(max_iterations.get(s, 3000), 'max_iterations')
                               for s in self.stage_constants}
        else:
            self.iterations = {s: self._integer(max_iterations, 'max_iterations')
                               for s in self.stage_constants}
        self.weights = (np.ones(self.natoms) if alignment_weights is None else
                        np.asarray(alignment_weights, dtype=float).copy())
        if (self.weights.shape != (self.natoms,) or not np.all(np.isfinite(self.weights))
                or np.any(self.weights <= 0)):
            raise ValueError('alignment_weights must be finite positive values per atom')

    @staticmethod
    def _positive(value, name, zero=False):
        value = float(value)
        if not np.isfinite(value) or (value < 0 if zero else value <= 0):
            raise ValueError(f'{name} must be finite and {"nonnegative" if zero else "positive"}')
        return value

    @staticmethod
    def _integer(value, name, zero=False):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise ValueError(f'{name} must be an integer')
        if value < (0 if zero else 1):
            raise ValueError(f'{name} is out of range')
        return int(value)

    @staticmethod
    def _prepare_conformer(obj):
        if isinstance(obj, (tuple, list)) and obj and all(hasattr(m, 'coords') for m in obj):
            if len(obj) > 2:
                raise ValueError('Only one or two molecules per endpoint are supported')
            if len(obj) == 1:
                return obj[0]
            # A lazy import avoids a Reactions initialization dependency cycle.
            from .Reaction import Reaction
            combined = Reaction.concenate_molecules(obj)
            charges = [m.charge for m in obj]
            if all(c is not None for c in charges):
                combined = combined.modify(charge=sum(charges))
            return combined
        return obj

    @staticmethod
    def _coordinates(obj):
        coords = np.asarray(getattr(obj, 'coords', obj), dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 3 or len(coords) == 0:
            raise ValueError('Each conformer must have shape (natoms, 3) with natoms > 0')
        if not np.all(np.isfinite(coords)):
            raise ValueError('Conformer coordinates must be finite')
        return coords.copy()

    @staticmethod
    def _atom_identities(atoms):
        return tuple((AtomData[a, 'Number'], AtomData[a, 'MassNumber']) for a in atoms)

    def _bond_map(self, bonds):
        result = {}
        for bond in bonds:
            if len(bond) < 2:
                raise ValueError('Bonds must contain two atom indices')
            i, j = bond[:2]
            if any(isinstance(a, (bool, np.bool_)) or not isinstance(a, (int, np.integer))
                   or a < 0 or a >= self.natoms for a in (i, j)) or i == j:
                raise ValueError('Invalid bond atom indices')
            pair = tuple(sorted((int(i), int(j))))
            order = bond[2] if len(bond) > 2 else 1
            if pair in result and result[pair] != order:
                raise ValueError('Conflicting duplicate bonds')
            result[pair] = order
        return result

    def _partition(self, fragments, graph):
        if fragments is None:
            fragments = graph.get_fragments()
        raw = [np.asarray(f) for f in fragments]
        if len(raw) not in (1, 2) or any(
            f.ndim != 1 or f.size == 0 or not np.issubdtype(f.dtype, np.integer)
            for f in raw
        ):
            raise ValueError('Provide one or two nonempty integer fragment index lists')
        if sorted(np.concatenate(raw).tolist()) != list(range(self.natoms)):
            raise ValueError('Fragments must partition all atom indices exactly once')
        fragments = tuple(np.sort(f).astype(int) for f in raw)
        owners = self._owners(fragments)
        if any(owners[i] != owners[j] for i, j in graph.edges):
            raise ValueError('A fragment partition cannot cut an endpoint bond')
        return fragments

    def _owners(self, fragments):
        owners = np.empty(self.natoms, dtype=int)
        for i, fragment in enumerate(fragments):
            owners[fragment] = i
        return owners

    def _reactive_sets(self, image):
        sets = [[set() for _ in self.fragments[image]] for _ in self.fragments[image]]
        changes = self.formed_bonds if image == 0 else self.broken_bonds
        for a, b in changes:
            n, m = self.owners[image][[a, b]]
            if n != m:
                sets[n][m].add(a)
                sets[m][n].add(b)
        return tuple(tuple(np.array(sorted(s), dtype=int) for s in row) for row in sets)

    def _molecular_radius(self, coords, fragment):
        distances = np.linalg.norm(coords[fragment] - coords[fragment].mean(axis=0), axis=1)
        h = distances.mean() + 2 * distances.std()
        return self.radii[fragment].max() if h < 1.e-12 else h

    @staticmethod
    def _centers(coords, fragments):
        return np.array([coords[f].mean(axis=0) for f in fragments])

    @staticmethod
    def _direction(v, fallback=None):
        norm = np.linalg.norm(v)
        if norm > 1.e-12:
            return v / norm
        return np.array([1., 0., 0.]) if fallback is None else fallback

    @staticmethod
    def _rotation_from_vector(vector):
        """Exponential-map rotation for row vectors, including zero motion."""
        angle = np.linalg.norm(vector)
        if angle == 0:
            return np.eye(3)
        # Normalize explicitly so small optimizer steps retain their axes.
        return nput.rotation_matrix(vector / angle, angle).T

    @classmethod
    def _rotation(cls, source, target, weights=None):
        """Proper least-squares rotation of row vectors, with no recentering."""
        source, target = np.atleast_2d(source), np.atleast_2d(target)
        covariance = source.T @ (target if weights is None else weights[:, None] * target)
        if np.linalg.norm(covariance) < 1.e-14:
            return np.eye(3)
        u, singular, vt = np.linalg.svd(covariance)
        if singular[1] <= 1.e-12 * singular[0]:
            # A single shared atom or collinear subset leaves a free twist.
            # Choose the smallest rotation instead of arbitrary SVD null axes.
            a, b = u[:, 0], vt[0]
            axis = np.cross(a, b)
            norm = np.linalg.norm(axis)
            cosine = np.clip(np.dot(a, b), -1., 1.)
            if norm < 1.e-12:
                if cosine >= 0:
                    return np.eye(3)
                basis = np.eye(3)[np.argmin(np.abs(a))]
                axis = cls._direction(np.cross(a, basis))
                angle = np.pi
            else:
                axis /= norm
                angle = np.arctan2(norm, cosine)
            return nput.rotation_matrix(axis, angle).T
        sign = 1. if np.linalg.det(u @ vt) >= 0 else -1.
        return u @ np.diag([1., 1., sign]) @ vt

    def _partner_center(self, coords, image, fragment):
        total = np.zeros(3)
        for n, row in enumerate(self.reactive_sets[image]):
            indices = row[fragment]
            if len(indices):
                total += coords[indices].mean(axis=0)
        return total / len(self.fragments[image])

    def _position(self, coords):
        # S1: sequential reactant placement, then product placement from B/C.
        for f in self.fragments[0]:
            coords[0, f] -= coords[0, f].mean(axis=0)
        for m, f in enumerate(self.fragments[0]):
            coords[0, f] += self._partner_center(coords[0], 0, m)
        for m, f in enumerate(self.fragments[1]):
            beta, sigma, nb = np.zeros(3), np.zeros(3), 0
            for row in self.shared_sets:
                shared = row[m]
                if len(shared):
                    sigma += coords[0, shared].mean(axis=0)
                    reactive = np.array([a for a in shared if a in self.reactive_atoms], dtype=int)
                    if len(reactive):
                        beta += coords[0, reactive].mean(axis=0)
                        nb += 1
            sigma /= len(self.fragments[0])
            beta = beta / len(self.fragments[0]) if nb else sigma
            coords[1, f] += (beta + sigma) / 2 - coords[1, f].mean(axis=0)

    def _orient(self, coords):
        # S3: reactive sites face their partners; shared subsets correlate P/R.
        centers = self._centers(coords[0], self.fragments[0])
        targets = [self._partner_center(coords[0], 0, n) - centers[n]
                   for n in range(len(self.fragments[0]))]
        for n, f in enumerate(self.fragments[0]):
            reactive = np.unique(np.concatenate(self.reactive_sets[0][n]))
            if len(reactive):
                vector = coords[0, reactive].mean(axis=0) - centers[n]
                if np.linalg.norm(vector) > 1.e-12 and np.linalg.norm(targets[n]) > 1.e-12:
                    rot = self._rotation(vector, targets[n])
                    coords[0, f] = (coords[0, f] - centers[n]) @ rot + centers[n]
        centers = self._centers(coords[0], self.fragments[0])
        for m, f in enumerate(self.fragments[1]):
            center = coords[1, f].mean(axis=0)
            source = coords[1, f] - center
            target, weight = np.zeros_like(source), 0
            for n, row in enumerate(self.shared_sets):
                shared = row[m]
                if len(shared):
                    rot = self._rotation(coords[1, shared] - center,
                                         coords[0, shared] - centers[n])
                    target += len(shared)**2 * (source @ rot)
                    weight += len(shared)**2
            coords[1, f] = source @ self._rotation(source, target / weight) + center

    def _molecular_force(self, coords, image, n, interimage=False):
        f = self.fragments[image][n]
        center = coords[image, f].mean(axis=0)
        other_image = 1 - image if interimage else image
        force, count = np.zeros(3), 0
        for m, other in enumerate(self.fragments[other_image]):
            if (not interimage and n == m) or (interimage and np.intersect1d(f, other).size):
                continue
            d = coords[other_image, other].mean(axis=0) - center
            radius = self.molecular_radii[image][n] + self.molecular_radii[other_image][m]
            distance = np.linalg.norm(d)
            if distance < radius:
                # An antisymmetric direction resolves exactly coincident spheres.
                sign = 1 if (image, n) < (other_image, m) else -1
                force += (distance - radius) * self._direction(d, np.array([sign, 0., 0.]))
                count += 1
        return force / (3 * len(f) * count) if count else force

    @staticmethod
    def _pair_vectors(source, target, center, weights=None, denominator=1):
        """Eqs. (4)/(6): ALL pairs between the specified atom subsets."""
        displacement = target[None, :, :] - source[:, None, :]
        distances = np.linalg.norm(displacement, axis=-1)
        direction = np.divide(displacement, distances[..., None],
                              out=np.zeros_like(displacement), where=distances[..., None] > 1.e-12)
        lever = source[:, None, :] - center
        scale = np.abs(np.sum(displacement * lever, axis=-1))
        # Sites at the fragment centre also need a translational force.
        scale = np.where(np.linalg.norm(lever, axis=-1) < 1.e-12, distances, scale)
        wt = np.ones((len(source), 1, 1)) if weights is None else weights[:, None, None]
        translation = np.sum(wt * scale[..., None] * direction, axis=(0, 1))
        rotation = np.sum(wt * np.cross(displacement, lever), axis=(0, 1))
        return translation / denominator, rotation / denominator

    def _attraction(self, coords, image, n, interimage=False):
        f = self.fragments[image][n]
        center = coords[image, f].mean(axis=0)
        translation, rotation, pairs = np.zeros(3), np.zeros(3), 0
        if interimage:
            for m in range(len(self.fragments[1 - image])):
                shared = self.shared_sets[n][m] if image == 0 else self.shared_sets[m][n]
                if len(shared):
                    weights = np.array([1. if a in self.reactive_atoms else .5 for a in shared])
                    t, r = self._pair_vectors(coords[image, shared], coords[1 - image, shared], center, weights)
                    translation += t
                    rotation += r
                    pairs += len(shared)**2
        else:
            for m in range(len(self.fragments[image])):
                a, b = self.reactive_sets[image][n][m], self.reactive_sets[image][m][n]
                if len(a) and len(b):
                    t, r = self._pair_vectors(coords[image, a], coords[image, b], center)
                    translation += t
                    rotation += r
                    pairs += len(a) * len(b)
        if pairs:
            translation /= 3 * len(f) * pairs
            rotation /= 3 * len(f) * pairs
        lever = coords[image, f] - center
        return translation - np.cross(rotation, lever)

    def _atomic_force(self, coords, image, n):
        # Eq. (7), with an explicitly outward force and phi*c_ab target.
        f = self.fragments[image][n]
        center = coords[image, f].mean(axis=0)
        translation, rotation, count = np.zeros(3), np.zeros(3), 0
        for m, other in enumerate(self.fragments[image]):
            if m == n:
                continue
            d = coords[image, other][None, :, :] - coords[image, f][:, None, :]
            distances = np.linalg.norm(d, axis=-1)
            direction = np.divide(d, distances[..., None], out=np.zeros_like(d),
                                  where=distances[..., None] > 1.e-12)
            direction[distances <= 1.e-12] = [1 if n < m else -1, 0., 0.]
            reactive = np.isin(f, self.reactive_sets[image][n][m])
            phi = np.where(reactive, 1.5, 2.)[:, None]
            contacts = phi * self.radius_scaling * (self.radii[f, None] + self.radii[other][None, :]) / 2
            active = distances < contacts
            y = np.minimum(distances - contacts, 0)[..., None] * direction
            ynorm = np.linalg.norm(y, axis=-1)
            unit_y = np.divide(y, ynorm[..., None], out=np.zeros_like(y), where=ynorm[..., None] > 1.e-12)
            lever = coords[image, f][:, None, :] - center
            scale = np.abs(np.sum(y * lever, axis=-1))
            scale = np.where(np.linalg.norm(lever, axis=-1) < 1.e-12, ynorm, scale)
            translation += np.sum(scale[..., None] * unit_y, axis=(0, 1))
            rotation += np.sum(np.cross(y, lever), axis=(0, 1))
            count += np.count_nonzero(active)
        if count:
            translation /= len(f)**2 * count
            rotation /= len(f)**2 * count
        return translation - np.cross(rotation, coords[image, f] - center)

    def _forces(self, coords, stage):
        molecular, reactive, shared, spectator, atomic = self.stage_constants[stage]
        forces = np.zeros_like(coords)
        for image, fragments in enumerate(self.fragments):
            for n, f in enumerate(fragments):
                if molecular:
                    forces[image, f] += molecular * self._molecular_force(coords, image, n)
                if reactive:
                    forces[image, f] += reactive * self._attraction(coords, image, n)
                if shared:
                    forces[image, f] += shared * self._attraction(coords, image, n, interimage=True)
                if spectator:
                    forces[image, f] += spectator * self._molecular_force(coords, image, n, interimage=True)
                if atomic:
                    forces[image, f] += atomic * self._atomic_force(coords, image, n)
        return forces

    def _rigid_velocities(self, coords, forces):
        """Block form of the appendix's rigid orthogonal projection.

        The pseudoinverse handles the missing rotation axes of linear/point
        fragments, without constructing a dense 6*nfragments projector.
        """
        velocities, projected = [], np.zeros_like(coords)
        for image, fragments in enumerate(self.fragments):
            row = []
            for f in fragments:
                lever = coords[image, f] - coords[image, f].mean(axis=0)
                translation = forces[image, f].mean(axis=0)
                inertia = np.eye(3) * np.sum(lever**2) - lever.T @ lever
                torque = np.sum(np.cross(lever, forces[image, f]), axis=0)
                angular = np.linalg.pinv(inertia, rcond=1.e-12) @ torque
                projected[image, f] = translation + np.cross(angular, lever)
                row.append((translation, angular))
            velocities.append(row)
        return velocities, projected

    def _move(self, coords, velocities, step):
        for image, fragments in enumerate(self.fragments):
            for f, (translation, angular) in zip(fragments, velocities[image]):
                center = coords[image, f].mean(axis=0)
                coords[image, f] = ((coords[image, f] - center)
                                    @ self._rotation_from_vector(step * angular)
                                    + center + step * translation)

    def _remove_common_motion(self, coords, velocities, projected):
        """Remove rigid motion of the combined R/P system, a free gauge."""
        lever = coords - coords.mean(axis=(0, 1))
        translation = projected.mean(axis=(0, 1))
        flat = lever.reshape(-1, 3)
        inertia = np.eye(3) * np.sum(flat**2) - flat.T @ flat
        torque = np.sum(np.cross(lever, projected), axis=(0, 1))
        angular = np.linalg.pinv(inertia, rcond=1.e-12) @ torque
        projected -= translation + np.cross(angular, lever)
        for image, fragments in enumerate(self.fragments):
            for f, (t, r) in zip(fragments, velocities[image]):
                t -= translation + np.cross(angular, lever[image, f].mean(axis=0))
                r -= angular
        return projected

    def _relax(self, coords, stage, rng):
        previous, step = None, self.step_size
        for iteration in range(self.iterations[stage] + 1):
            forces = self._forces(coords, stage)
            velocities, projected = self._rigid_velocities(coords, forces)
            projected = self._remove_common_motion(coords, velocities, projected)
            residual = float(np.max(np.linalg.norm(projected, axis=-1)))
            noisy = (stage != 'S2' and self.perturbation > 0
                     and iteration < min(self.perturbation_steps, self.iterations[stage]))
            if not np.isfinite(residual):
                raise FloatingPointError(f'Nonfinite alignment force in {stage}')
            if (residual <= self.tolerance and not noisy) or iteration == self.iterations[stage]:
                evaluations = 0
                if residual > self.tolerance and stage != 'S2' and self.polish_max_evaluations:
                    polished, new_residual, evaluations = self._polish(coords, stage)
                    if new_residual < residual:
                        coords[:] = polished
                        residual = new_residual
                return ConformerAlignmentStage(stage, coords.copy() / self.to_angstrom,
                                                iteration, residual, residual <= self.tolerance,
                                                evaluations)
            if previous is not None:
                if np.sum(previous * projected) < 0:
                    step = max(step * .5, self.step_size * 1.e-4)
                else:
                    step = min(step * 1.05, self.step_size)
            previous = projected
            if noisy:
                for row in velocities:
                    for translation, angular in row:
                        translation += self.perturbation * rng.normal(size=3)
                        angular += self.perturbation * rng.normal(size=3)
            max_motion = max(
                np.linalg.norm(t) + np.linalg.norm(r) * np.max(np.linalg.norm(
                    coords[image, f] - coords[image, f].mean(axis=0), axis=-1
                ))
                for image, fragments in enumerate(self.fragments)
                for f, (t, r) in zip(fragments, velocities[image])
            )
            max_angular = max(np.linalg.norm(r) for row in velocities for _, r in row)
            move_step = min(step, self.max_displacement / max(max_motion, 1.e-14),
                            self.max_rotation / max(max_angular, 1.e-14))
            self._move(coords, velocities, move_step)
            # Remove a common translation only; relative endpoint positions remain.
            coords -= coords.mean(axis=(0, 1))

    def _polish(self, coords, stage):
        """Resolve slow/stalled force relaxation within exact rigid poses."""
        base = coords.copy()
        parts = [(image, f) for image, fragments in enumerate(self.fragments) for f in fragments]

        def decode(pose):
            trial = base.copy()
            for vector, (image, f) in zip(pose.reshape(-1, 6), parts):
                center = base[image, f].mean(axis=0)
                trial[image, f] = ((base[image, f] - center)
                                    @ self._rotation_from_vector(vector[3:])
                                    + center + vector[:3])
            return trial

        def objective(pose):
            trial = decode(pose)
            velocities, projected = self._rigid_velocities(trial, self._forces(trial, stage))
            return self._remove_common_motion(trial, velocities, projected).ravel()

        solution = least_squares(objective, np.zeros(6 * len(parts)),
                                max_nfev=self.polish_max_evaluations,
                                ftol=1.e-10, xtol=1.e-10, gtol=1.e-10)
        residual = float(np.max(np.linalg.norm(solution.fun.reshape(2, self.natoms, 3), axis=-1)))
        # Solver success alone does not establish convergence of the force field.
        return decode(solution.x), residual, solution.nfev

    def _finalize(self, coords):
        if self.endpoint_separation:
            for image, fragments in enumerate(self.fragments):
                if len(fragments) == 2:
                    centers = self._centers(coords[image], fragments)
                    outward = self._direction(centers[0] - centers[1])
                    coords[image, fragments[0]] += self.endpoint_separation * outward
                    coords[image, fragments[1]] -= self.endpoint_separation * outward
        # One final proper whole-image alignment; no interpolation is needed.
        coords -= np.average(coords, axis=1, weights=self.weights)[:, None, :]
        coords[1] = coords[1] @ self._rotation(coords[1], coords[0], self.weights)

    def align(self, *, strict=False):
        """Run the five stages on copies and return diagnostics with endpoints."""
        coords = self.initial.copy()
        rng = np.random.default_rng(self.random_seed)
        self._position(coords)
        stages = [ConformerAlignmentStage('S1', coords.copy() / self.to_angstrom)]
        stages.append(self._relax(coords, 'S2', rng))
        self._orient(coords)
        stages.append(ConformerAlignmentStage('S3', coords.copy() / self.to_angstrom))
        # Unimolecular endpoints have no fragment-placement freedom. The proper
        # Kabsch fit is exact; applying Eq. (6) cannot improve this rigid fit.
        if all(len(f) == 1 for f in self.fragments):
            for stage in ('S4', 'S5'):
                stages.append(ConformerAlignmentStage(stage, coords.copy() / self.to_angstrom))
            self._finalize(coords)
        else:
            for stage in ('S4', 'S5'):
                stages.append(self._relax(coords, stage, rng))
            self._finalize(coords)
        outputs = []
        for image, original in enumerate((self.react, self.prod)):
            new_coords = coords[image] / self.to_angstrom
            outputs.append(original.modify(coords=new_coords) if hasattr(original, 'modify') else new_coords.copy())
        result = ConformerAlignmentResult(*outputs,
                    tuple(tuple(f.tolist()) for f in self.fragments[0]),
                    tuple(tuple(f.tolist()) for f in self.fragments[1]),
                    self.formed_bonds, self.broken_bonds, tuple(stages), self.coordinate_units)
        if not result.converged:
            failed = ', '.join(s.name for s in stages if not s.converged)
            message = f'Conformer alignment did not converge in {failed}; inspect result.stages'
            if strict:
                error = RuntimeError(message)
                error.result = result
                raise error
            warnings.warn(message, RuntimeWarning, stacklevel=2)
        return result
