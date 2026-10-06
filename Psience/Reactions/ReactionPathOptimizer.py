"""Adapters from a selected reaction interpolation to native McUtils optimizers.

McUtils exposes its optimizer framework through Numputils.Optimization. These
hooks use its chain runner and NEB/string step finders, and its single-point
runner and eigenvalue-following step finder, without implementing another
optimizer. Cartesian coordinates are in Bohr, energies in Hartree, and forces
and Hessians in Hartree/Bohr and Hartree/Bohr**2 respectively.
"""

from dataclasses import dataclass
import warnings

import numpy as np
import McUtils.Numputils as nput
from McUtils.Data import UnitsData

from .ProfileGenerator import ProfileGenerator, InterpolatingProfileGenerator

__all__ = [
    'ReactionPathOptimizer', 'ReactionOptimizationResult',
    'EigenvectorFollowingProfileGenerator'
]


@dataclass
class ReactionOptimizationResult:
    method: str
    coordinates: np.ndarray
    images: tuple
    energies: np.ndarray
    converged: bool
    optimizer_converged: bool
    force_error: float
    errors: np.ndarray
    iterations: np.ndarray
    image_numbers: object = None
    eigenvalues: object = None
    saddle_order: object = None
    climbing_image: object = None
    trajectory: tuple = ()
    coordinate_search: object = None

    @property
    def highest_energy_image(self):
        return self.images[int(np.argmax(self.energies))]

    @property
    def transition_state(self):
        """The refined saddle candidate; inspect converged before using it."""
        if self.method != 'eigenvector-following':
            raise ValueError('Refine a path image with eigenvector following first')
        return self.images[0]


class _MolecularChainProjection:
    """Apply the same rigid-motion convention to native tangents and steps."""

    def __init__(self, *args, cartesian_projector, **options):
        self.cartesian_projector = cartesian_projector
        super().__init__(*args, **options)

    def get_tangent(self, guess, mask, cur, prev, next):
        tangent = super().get_tangent(guess, mask, cur, prev, next)
        projector = self.cartesian_projector(guess[:, cur])
        tangent = (projector @ tangent[..., None])[..., 0]
        norm = np.linalg.norm(tangent, axis=-1, keepdims=True)
        return tangent / np.maximum(norm, 1.e-12)

    def __call__(self, guess, mask, projector=None, **options):
        cur = mask[1][0]
        molecular_projector = self.cartesian_projector(guess[:, cur])
        # Native chain finders split mixed climbing/nonclimbing batches.
        # Project the reassembled outputs so each subset retains its own shape.
        step, gradient = super().__call__(guess, mask, projector=projector, **options)
        return ((molecular_projector @ step[..., None])[..., 0],
                (molecular_projector @ gradient[..., None])[..., 0])


class _MolecularNEBStepFinder(_MolecularChainProjection, nput.NudgedElasticBandStepFinder):
    pass


class _MolecularStringStepFinder(_MolecularChainProjection, nput.StringMethodStepFinder):
    pass


@ProfileGenerator.register('reaction-path')
class ReactionPathOptimizer(InterpolatingProfileGenerator):
    """Connect an initial fixed-frame interpolation to a physical potential.

``from_coordinate_search`` retains the selected globally valid Z-matrix and
its sampled parameter grid. Optimizers subsequently relax Cartesian images;
the initial interpolation certificate does not certify the relaxed path.
No per-image re-embedding is performed, and chain endpoints remain fixed.

``energy_evaluator`` follows Molecule's EnergyEvaluator protocol: evaluate
accepts coordinates in evaluator.distance_units and returns expansion terms
already normalized to Hartree/Bohr**order. Only inputs are converted here.
The coordinate-selection metric and its steric penalty are not a PES.

``get_callbacks`` and ``get_step_finder`` are hooks for direct use with the
McUtils runners. Callbacks accept flattened batches and an optional mask.
``optimize_path`` and ``find_transition_state`` add endpoint handling and
convergence diagnostics; ``generate`` implements the ProfileGenerator API.

For an isolated molecule, remove_rigid_motion removes overall translation
and rotation in the Euclidean Cartesian metric, using Numputils projectors
with unit masses. Turn it off for an external field or a model with physical
Cartesian degrees of freedom. Eigenvector following stabilizes excluded
rigid modes above the physical spectrum by gauge_curvature and checks saddle
order in the remaining physical subspace. target_mode indices refer to that spectrum.
For molecular chains, tangents and steps use that same projector. Spring
distances and string redistribution still use the fixed-frame Cartesian arc.
"""

    def __init__(self, reactant_complex, product_complex, *, energy_evaluator,
                 method='neb', coordinate_interpolator=None, num_images=21,
                 initial_image_positions=None, remove_rigid_motion=True,
                 gauge_curvature=1., saddle_eigenvalue_tolerance=1.e-7,
                 internals=None):
        endpoints = np.asarray([reactant_complex.coords, product_complex.coords], dtype=float)
        if (endpoints.ndim != 3 or endpoints.shape[-1] != 3 or endpoints.shape[1] == 0
                or not np.all(np.isfinite(endpoints))):
            raise ValueError('Endpoints must be finite conformers with equal (natoms, 3) shapes')
        if tuple(reactant_complex.atoms) != tuple(product_complex.atoms):
            raise ValueError('Endpoint atoms must agree index by index')
        self.natoms = endpoints.shape[1]
        self.ndim = 3 * self.natoms
        self.endpoints = endpoints.copy()
        self.method = self._method(method)
        self.remove_rigid_motion = bool(remove_rigid_motion)
        if not np.isfinite(gauge_curvature) or gauge_curvature <= 0:
            raise ValueError('gauge_curvature must be positive and finite')
        if not np.isfinite(saddle_eigenvalue_tolerance) or saddle_eigenvalue_tolerance < 0:
            raise ValueError('saddle_eigenvalue_tolerance must be nonnegative and finite')
        self.gauge_curvature = float(gauge_curvature)
        self.saddle_eigenvalue_tolerance = float(saddle_eigenvalue_tolerance)
        self._energy_evaluator = reactant_complex.get_energy_evaluator(energy_evaluator)
        self.coordinate_search = None
        self.last_result = None
        self.internals = internals
        if internals is not None:
            if coordinate_interpolator is not None:
                raise ValueError('Supply a Z-matrix or coordinate_interpolator, not both')
            from .ReactionCoordinateSearch import ZMatrixReactionInterpolator
            coordinate_interpolator = ZMatrixReactionInterpolator(endpoints[0], endpoints[1], internals)
        if coordinate_interpolator is None:
            # A Cartesian fallback also preserves the supplied embedding.
            def coordinate_interpolator(t):
                t = np.asarray(t, dtype=float)
                return endpoints[0] + t[..., None, None] * (endpoints[1] - endpoints[0])
        super().__init__(reactant_complex, product_complex,
                         coordinate_interpolator=coordinate_interpolator,
                         num_images=num_images, initial_image_positions=initial_image_positions)

    @classmethod
    def from_coordinate_search(cls, result, energy_evaluator, **options):
        if result.reactant is None or result.product is None:
            raise ValueError('The coordinate search result must retain its endpoint Molecules')
        options.setdefault('num_images', len(result.parameters))
        options.setdefault('initial_image_positions', result.parameters.copy())
        optimizer = cls(result.reactant, result.product, energy_evaluator=energy_evaluator,
                        coordinate_interpolator=result.best.interpolator, **options)
        optimizer.coordinate_search = result
        return optimizer

    @staticmethod
    def _method(method):
        aliases = {'neb': 'neb', 'string': 'string', 'growing-string': 'string',
                   'eigenvector-following': 'eigenvector-following',
                   'eigenvalue-following': 'eigenvector-following'}
        try:
            return aliases[method]
        except (KeyError, TypeError):
            raise ValueError(f'Unknown reaction optimizer {method!r}') from None

    @property
    def energy_evaluator(self):
        return self._energy_evaluator

    def _points(self, guess):
        guess = np.asarray(guess, dtype=float)
        if guess.ndim < 1 or guess.shape[-1] != self.ndim or not np.all(np.isfinite(guess)):
            raise ValueError(f'Optimizer coordinates must be finite with final dimension {self.ndim}')
        return guess

    def _evaluate(self, guess, order):
        guess = self._points(guess)
        coords = guess.reshape(-1, self.natoms, 3)
        coords = coords * UnitsData.convert('BohrRadius', self.energy_evaluator.distance_units)
        term = np.asarray(self.energy_evaluator.evaluate(coords, order=[order])[0], dtype=float)
        shape = guess.shape[:-1] + (self.ndim,) * order
        if term.size != int(np.prod(shape, dtype=int)) or not np.all(np.isfinite(term)):
            raise ValueError(f'Energy evaluator returned an invalid order-{order} term')
        return term.reshape(shape)

    def _projector(self, guess):
        guess = self._points(guess)
        if not self.remove_rigid_motion:
            return np.broadcast_to(np.eye(self.ndim), guess.shape[:-1] + (self.ndim, self.ndim))
        # Separate evaluations allow batches containing both linear and
        # nonlinear conformers, whose rigid-mode counts differ.
        projectors = [nput.translation_rotation_projector(c, orthonormal=True)
                      for c in guess.reshape(-1, self.natoms, 3)]
        return np.asarray(projectors).reshape(guess.shape[:-1] + (self.ndim, self.ndim))

    def _potential(self, guess, mask=None):
        return self._evaluate(guess, 0)

    def _gradient(self, guess, mask=None):
        gradient = self._evaluate(guess, 1)
        if self.remove_rigid_motion:
            gradient = (self._projector(guess) @ gradient[..., None])[..., 0]
        return gradient

    def _hessian(self, guess, mask=None):
        hessian = self._evaluate(guess, 2)
        hessian = (hessian + hessian.swapaxes(-1, -2)) / 2
        if self.remove_rigid_motion:
            projector = self._projector(guess)
            hessian = projector @ hessian @ projector
            gauge = self.gauge_curvature + np.max(np.abs(np.linalg.eigvalsh(hessian)), axis=-1)
            hessian = hessian + gauge[..., None, None] * (np.eye(self.ndim) - projector)
        return hessian

    def get_callbacks(self):
        """Return value, gradient, Hessian callbacks for flattened Bohr guesses."""
        return self._potential, self._gradient, self._hessian

    def get_step_finder(self, method=None, **options):
        method = self._method(self.method if method is None else method)
        value, gradient, hessian = self.get_callbacks()
        if method == 'neb':
            if self.remove_rigid_motion:
                return _MolecularNEBStepFinder(value, gradient, cartesian_projector=self._projector, **options)
            return nput.NudgedElasticBandStepFinder(value, gradient, **options)
        if method == 'string':
            if self.remove_rigid_motion:
                return _MolecularStringStepFinder(value, gradient, cartesian_projector=self._projector, **options)
            return nput.StringMethodStepFinder(value, gradient, **options)
        return nput.EigenvalueFollowingStepFinder(value, gradient, hessian, **options)

    def _initial_images(self, initial_images=None, num_images=None, initial_image_positions=None):
        if initial_images is None:
            initial_images = super().generate(num_images=num_images,
                                              initial_image_positions=initial_image_positions)
        coords = np.asarray([i.coords if hasattr(i, 'coords') else i for i in initial_images], dtype=float)
        if (coords.ndim != 3 or coords.shape[1:] != (self.natoms, 3) or len(coords) < 3
                or not np.all(np.isfinite(coords))):
            raise ValueError('A path needs at least three finite images with shape (natoms, 3)')
        if not np.allclose(coords[[0, -1]], self.endpoints, rtol=0, atol=1.e-7):
            raise ValueError('Initial path endpoints must match the supplied conformers')
        coords = coords.copy()
        coords[[0, -1]] = self.endpoints
        return coords

    @staticmethod
    def _error(gradient, use_max=True):
        if np.size(gradient) == 0:
            return 0.
        error = np.max(np.abs(gradient), axis=-1) if use_max else np.linalg.norm(gradient, axis=-1)
        return float(np.max(error))

    def _path_residual(self, finder, coords, fixed_images, climbing_image, use_max):
        chain = coords.reshape(1, len(coords), self.ndim)
        gradients = self._gradient(chain)[0]
        residuals = []
        for j in range(1, len(coords) - 1):
            if j in fixed_images:
                continue
            tangent = finder.get_tangent(chain, np.array([0]), j, j - 1, j + 1)[0]
            grad = gradients[j]
            if j == climbing_image:
                grad = grad - 2 * np.dot(grad, tangent) * tangent
            else:
                grad = grad - np.dot(grad, tangent) * tangent
                if isinstance(finder, nput.StringMethodStepFinder):
                    spring = 0
                else:
                    k = finder.spring_constants
                    k = k if np.ndim(k) == 0 else k[j]
                    spring = k * (np.linalg.norm(chain[0, j] - chain[0, j - 1])
                                  - np.linalg.norm(chain[0, j + 1] - chain[0, j])) * tangent
                grad = grad + spring
            residuals.append(grad)
        return self._error(np.asarray(residuals), use_max)

    def _finish(self, result, strict):
        self.last_result = result
        if not result.converged:
            message = (f'{result.method} did not converge: final force error {result.force_error:.3g}'
                       + (f', saddle order {result.saddle_order}' if result.saddle_order is not None else ''))
            if strict:
                error = RuntimeError(message)
                error.result = result
                raise error
            warnings.warn(message, RuntimeWarning, stacklevel=3)
        return result

    def optimize_path(self, method=None, *, initial_images=None, num_images=None,
                      initial_image_positions=None, step_finder_options=None,
                      optimizer_settings=None, climb=False, climbing_image=None,
                      reparametrizer=None, strict=False, **options):
        """Relax a fixed-endpoint Cartesian path using NEB or the string method.

Step-finder options (step_size, spring_constants, ...) and runner options
(tol, max_iterations, max_displacement_norm, ...) are separate dictionaries.
NEB climbing selects the highest initial interior image unless an index is
given. The string method uses native arc-length redistribution. For a saddle,
refine its highest-energy image with find_transition_state.
"""
        method = self._method(self.method if method is None else method)
        if method not in ('neb', 'string'):
            raise ValueError('optimize_path requires neb or string')
        coords = self._initial_images(initial_images, num_images, initial_image_positions)
        settings = dict(optimizer_settings or {}, **options)
        if (settings.get('periodic', False) or settings.get('reembed', False)
                or settings.get('unitary', False) or settings.get('generate_rotation', False)):
            raise ValueError('Reaction paths require fixed endpoints in their supplied Cartesian frame')
        if settings.get('orthogonal_projection_generator') is not None:
            raise ValueError('Rigid-motion projection is supplied by the reaction callbacks')
        settings['reembed'] = False
        fixed = settings.pop('fixed_images', None)
        fixed = sorted(set(() if fixed is None else fixed) | {0, len(coords) - 1})
        if any(isinstance(j, (bool, np.bool_)) or not isinstance(j, (int, np.integer))
               or not 0 <= j < len(coords) for j in fixed):
            raise ValueError('fixed_images must contain valid image indices')
        if climbing_image is not None:
            climb = True
        if climb:
            if method != 'neb':
                raise ValueError('Use NEB for climbing, or refine a string image with eigenvector following')
            if climbing_image is None:
                free = [j for j in range(1, len(coords) - 1) if j not in fixed]
                if not free:
                    raise ValueError('Climbing requires a free interior image')
                energies = self._potential(coords.reshape(len(coords), -1))
                climbing_image = free[int(np.argmax(energies[free]))]
            if (isinstance(climbing_image, (bool, np.bool_))
                    or not isinstance(climbing_image, (int, np.integer))
                    or not 0 < climbing_image < len(coords) - 1 or climbing_image in fixed):
                raise ValueError('climbing_image must be a free interior image index')
        # Each image owns any quasi-Newton/conjugate-gradient history used by
        # a caller-supplied inner step finder. Native runners also accept a
        # singleton, which would share that state between different images.
        finders = [self.get_step_finder(method, **(step_finder_options or {})) for _ in coords]
        if method == 'string' and reparametrizer is None:
            reparametrizer = nput.InterpolatingReparametrizer()
        settings.setdefault('logger', False)
        raw = nput.iterative_chain_minimize(
            coords.reshape(len(coords), -1), finders, fixed_images=fixed,
            climb=climb, climbing_nodes=[climbing_image] if climb else None,
            reparametrizer=reparametrizer, **settings)
        (payload, image_numbers), native_converged, (errors, iterations) = raw
        if settings.get('return_trajectory', False):
            optimized, trajectory = payload
            trajectory = tuple(t.reshape(-1, self.natoms, 3) for t in trajectory)
        else:
            optimized, trajectory = payload, ()
        optimized = np.asarray(optimized).reshape(-1, self.natoms, 3)
        if not np.array_equal(optimized[fixed], coords[fixed]):
            raise ValueError('The reparametrizer moved a fixed image')
        force_error = self._path_residual(finders[0], optimized, fixed, climbing_image,
                                          settings.get('use_max_for_error', True))
        tol = settings.get('tol', 1.e-8)
        result = ReactionOptimizationResult(
            method, optimized, tuple((self.products if j == len(optimized) - 1 else self.reactants)
                                     .modify(coords=c.copy()) for j, c in enumerate(optimized)),
            self._potential(optimized.reshape(len(optimized), -1)),
            bool(native_converged and force_error < tol), bool(native_converged), force_error,
            np.asarray(errors), np.asarray(iterations), image_numbers=image_numbers,
            climbing_image=climbing_image, trajectory=trajectory, coordinate_search=self.coordinate_search)
        return self._finish(result, strict)

    def _physical_modes(self, guess):
        projector = self._projector(guess)
        values, vectors = np.linalg.eigh(projector)
        basis = vectors[:, values > .5]
        hessian = self._evaluate(guess, 2)
        hessian = (hessian + hessian.T) / 2
        values, vectors = np.linalg.eigh(basis.T @ hessian @ basis)
        return values, basis @ vectors

    def find_transition_state(self, *, initial_guess=None, initial_images=None,
                              target_mode=None, step_finder_options=None,
                              optimizer_settings=None, strict=False, **options):
        """Refine a candidate with native eigenvector following and check its index.

Without a guess, start at the highest interior image of the last optimized
path, or of the selected interpolation. Its neighbor tangent seeds the mode.
With an explicit guess and no target, follow the lowest physical Hessian mode.
A zero-gradient minimum is reported as unconverged even if the runner stops.
"""
        tangent = None
        if initial_guess is None:
            if initial_images is None and self.last_result is not None and self.last_result.method in ('neb', 'string'):
                initial_images = self.last_result.coordinates
            coords = self._initial_images(initial_images)
            energies = self._potential(coords.reshape(len(coords), -1))
            j = int(np.argmax(energies[1:-1])) + 1
            initial_guess = coords[j]
            tangent = (coords[j + 1] - coords[j - 1]).reshape(-1)
        guess = np.asarray(initial_guess.coords if hasattr(initial_guess, 'coords') else initial_guess, dtype=float)
        if guess.shape not in ((self.natoms, 3), (self.ndim,)):
            raise ValueError('A TS guess must be a single conformer with the endpoint atom ordering')
        guess = self._points(guess.reshape(-1)).copy()
        values, modes = self._physical_modes(guess)
        if len(values) == 0:
            raise ValueError('There are no physical degrees of freedom after removing rigid motion')
        finder_options = dict(step_finder_options or {})
        if target_mode is None:
            target_mode = finder_options.pop('target_mode', None)
        elif 'target_mode' in finder_options:
            raise ValueError('Supply target_mode only once')
        use_path_tangent = target_mode is None and tangent is not None
        if target_mode is None:
            target_mode = tangent if tangent is not None else modes[:, 0]
        if isinstance(target_mode, (int, np.integer)) and not isinstance(target_mode, (bool, np.bool_)):
            if not 0 <= target_mode < len(values):
                raise ValueError('target_mode index must refer to the physical Hessian spectrum')
            target_mode = modes[:, target_mode]
        else:
            target_mode = np.asarray(target_mode, dtype=float)
            if target_mode.shape not in ((self.natoms, 3), (self.ndim,)):
                raise ValueError('target_mode must be a physical mode index or Cartesian vector')
            target_mode = self._projector(guess) @ target_mode.reshape(-1)
        norm = np.linalg.norm(target_mode)
        if not np.isfinite(norm) or norm < 1.e-12:
            if use_path_tangent and np.isfinite(norm):
                target_mode, norm = modes[:, 0], 1.
            else:
                raise ValueError('target_mode has no finite physical displacement')
        finder = self.get_step_finder('eigenvector-following', target_mode=target_mode / norm, **finder_options)
        settings = dict(optimizer_settings or {}, **options)
        if settings.get('track_best', False):
            raise ValueError('Energy-minimum tracking is unsuitable for a saddle search')
        if settings.get('unitary', False) or settings.get('generate_rotation', False):
            raise ValueError('TS refinement uses Cartesian coordinates')
        if self.remove_rigid_motion:
            if settings.get('orthogonal_projection_generator') is not None:
                raise ValueError('Rigid-motion projection is supplied by the reaction callbacks')
            settings['orthogonal_projection_generator'] = self._projector
            settings.setdefault('prevent_oscillations', False)
        settings.setdefault('logger', False)
        settings.setdefault('max_displacement_norm', .2)
        raw = nput.iterative_step_minimize(guess, finder, **settings)
        if settings.get('return_trajectory', False):
            raw, trajectory = raw
            trajectory = tuple(t.reshape(-1, self.natoms, 3) for _, t in trajectory)
        else:
            trajectory = ()
        optimized, native_converged, (errors, iterations) = raw
        optimized = np.asarray(optimized).reshape(-1)
        eigenvalues, _ = self._physical_modes(optimized)
        saddle_order = int(np.count_nonzero(eigenvalues < -self.saddle_eigenvalue_tolerance))
        force_error = self._error(self._gradient(optimized), settings.get('use_max_for_error', True))
        tol = settings.get('tol', 1.e-8)
        coords = optimized.reshape(1, self.natoms, 3)
        result = ReactionOptimizationResult(
            'eigenvector-following', coords,
            (self.reactants.modify(coords=coords[0].copy()),),
            self._potential(optimized).reshape(1),
            bool(native_converged and force_error < tol and saddle_order == 1),
            bool(native_converged), force_error, np.asarray(errors), np.asarray(iterations),
            eigenvalues=eigenvalues, saddle_order=saddle_order,
            trajectory=trajectory, coordinate_search=self.coordinate_search)
        return self._finish(result, strict)

    def optimize(self, method=None, **options):
        method = self._method(self.method if method is None else method)
        if method == 'eigenvector-following':
            return self.find_transition_state(**options)
        return self.optimize_path(method, **options)

    def generate(self, *, method=None, return_result=False, **options):
        """Return optimized Molecules, or the complete diagnostics on request."""
        result = self.optimize(method, **options)
        return result if return_result else list(result.images)

    def evaluate_profile_energies(self, profile, **options):
        coords = np.asarray([p.coords for p in profile])
        return self._potential(coords.reshape(len(coords), -1))

    def evaluate_profile_distances(self, profile, normalize=True):
        coords = np.asarray([p.coords for p in profile]).reshape(len(profile), -1)
        distances = np.concatenate(([0.], np.cumsum(np.linalg.norm(np.diff(coords, axis=0), axis=-1))))
        if normalize and distances[-1] > 0:
            distances = distances / distances[-1]
        return distances


@ProfileGenerator.register('eigenvalue-following')
@ProfileGenerator.register('eigenvector-following')
class EigenvectorFollowingProfileGenerator(ReactionPathOptimizer):
    """ProfileGenerator.resolve hook for a single native TS refinement."""

    def __init__(self, *args, method='eigenvector-following', **options):
        super().__init__(*args, method=method, **options)
