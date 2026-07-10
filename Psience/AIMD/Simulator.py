import collections, functools, numpy as np
import abc
from McUtils.Data import AtomData, UnitsData
from McUtils.Zachary import FiniteDifferenceDerivative, RBFDInterpolator
import McUtils.Numputils as nput

from ..Molecools import Molecule, MolecularZMatrixCoordinateSystem
from ..Molecools.Properties import PropertyManager

__all__ = [
    "AIMDSimulator",
    "PairwisePotential",
    "MDThermostat",
    "MDStepPredictor"
]


# ============================================================
# Thermostats
# ============================================================

class MDThermostat:
    """
    Base class for MD thermostats. Subclasses register themselves under
    `name` so they can be constructed from a plain string or a dict spec.
    """
    name = None
    _registry = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name is not None:
            MDThermostat._registry[cls.name] = cls

    def __init__(self, target_temperature=None, **opts):
        if target_temperature is None:
            raise ValueError(f"{type(self).__name__} requires `target_temperature`")
        self.target_temperature = target_temperature
        self.opts = opts  # unconsumed kwargs kept around for introspection/debugging

    def apply_thermostat(self, positions, velocities, sim=None):
        """
        Adjust (positions, velocities) for one integration step.
        `sim` is the parent AIMDSimulator, passed for access to dt, masses,
        instantaneous temperature, etc. Must return (positions, velocities).
        """
        raise NotImplementedError("subclasses must implement `apply_thermostat`")

    @classmethod
    def from_spec(cls, spec, **default_opts):
        """
        Build a thermostat from:
          - None                 -> None
          - MDThermostat instance -> returned as-is
          - str                  -> registry lookup, using `default_opts`
          - dict with 'type' key -> registry lookup, dict overrides `default_opts`
        """
        if spec is None:
            return None
        if isinstance(spec, MDThermostat):
            return spec

        if isinstance(spec, str):
            kind = spec
            opts = dict(default_opts)
        elif isinstance(spec, dict):
            spec = dict(spec)  # don't mutate caller's dict
            if 'type' not in spec:
                raise ValueError("thermostat dict spec must include a 'type' key")
            kind = spec.pop('type')
            opts = dict(default_opts)
            opts.update(spec)  # explicit dict entries win over defaults
        else:
            raise TypeError(f"can't build a thermostat from {spec!r}")

        if kind not in cls._registry:
            raise ValueError(
                f"unknown thermostat type '{kind}'; registered types: {sorted(cls._registry)}"
            )
        return cls._registry[kind](**opts)


class RescaleThermostat(MDThermostat):
    """Crude direct velocity rescaling to the target temperature every step."""
    name = 'rescale'

    def apply_thermostat(self, positions, velocities, sim=None):
        T = sim._temperature(velocities)
        scale = np.sqrt(self.target_temperature / T)
        scale = np.reshape(scale, scale.shape + (1,) * (velocities.ndim - np.ndim(scale)))
        return positions, velocities * scale


class BerendsenThermostat(MDThermostat):
    """Weak-coupling thermostat. Good for equilibration, not rigorous NVT sampling."""
    name = 'berendsen'

    def __init__(self, target_temperature=None, tau=None, **opts):
        super().__init__(target_temperature=target_temperature, **opts)
        self.tau = tau  # falls back to 100*dt if unset, resolved at call time

    def apply_thermostat(self, positions, velocities, sim=None):
        tau = self.tau if self.tau is not None else 100 * sim.dt
        T = sim._temperature(velocities)
        lam = np.sqrt(1 + (sim.dt / tau) * (self.target_temperature / T - 1))
        lam = np.reshape(lam, lam.shape + (1,) * (velocities.ndim - np.ndim(lam)))
        return positions, velocities * lam


class LangevinThermostat(MDThermostat):
    """BBK-style Langevin thermostat -- proper canonical sampling."""
    name = 'langevin'

    def __init__(self, target_temperature=None, friction=0.01, seed=None, **opts):
        super().__init__(target_temperature=target_temperature, **opts)
        self.friction = friction
        self._rng = np.random.default_rng(seed)

    def apply_thermostat(self, positions, velocities, sim=None):
        gamma = self.friction
        mass = sim._mass
        sigma = np.sqrt(2 * gamma * self.target_temperature * UnitsData.convert("Kelvins", "Hartrees") * mass / sim.dt)
        rand = self._rng.normal(size=velocities.shape)
        velocities = velocities + sim.dt * (-gamma * velocities) + np.sqrt(sim.dt) * sigma * rand / mass
        return positions, velocities

class MDStepPredictor(metaclass=abc.ABCMeta):
    """
    Base class for MD step integrators. Subclasses register themselves under
    `name` so they can be constructed from a plain string or a dict spec,
    exactly like MDThermostat.
    """
    name = None
    _registry = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name is not None:
            MDStepPredictor._registry[cls.name] = cls

    def __init__(self, **opts):
        self.opts = opts

    @abc.abstractmethod
    def predict_step(self, positions, velocities, forces, sim=None):
        """
        Advance one step. `forces` is the force at the *current* positions
        (possibly None, in which case the predictor should evaluate it via
        `sim.get_forces(positions)`). `sim` is the parent AIMDSimulator, for
        access to `dt`, `_mass`, `get_forces`, `_recenter`, etc.
        Must return (new_positions, new_velocities, new_forces).
        """
        raise NotImplementedError("subclasses must implement `predict_step`")

    @classmethod
    def from_spec(cls, spec, **default_opts):
        """
        Build a step predictor from:
          - None                  -> default ('velocity-verlet')
          - MDStepPredictor instance -> returned as-is
          - str                   -> registry lookup, using `default_opts`
          - dict with 'type' key  -> registry lookup, dict overrides `default_opts`
        """
        if spec is None:
            spec = 'velocity-verlet'
        if isinstance(spec, MDStepPredictor):
            return spec

        if isinstance(spec, str):
            kind = spec
            opts = dict(default_opts)
        elif isinstance(spec, dict):
            spec = dict(spec)
            if 'type' not in spec:
                raise ValueError("step predictor dict spec must include a 'type' key")
            kind = spec.pop('type')
            opts = dict(default_opts)
            opts.update(spec)
        else:
            raise TypeError(f"can't build a step predictor from {spec!r}")

        if kind not in cls._registry:
            raise ValueError(
                f"unknown step predictor type '{kind}'; registered types: {sorted(cls._registry)}"
            )
        return cls._registry[kind](**opts)


class VelocityVerletStepPredictor(MDStepPredictor):
    """Standard velocity-Verlet -- symplectic, time-reversible, 1 force eval/step."""
    name = 'velocity-verlet'

    def predict_step(self, positions, velocities, forces, sim=None):
        if forces is None:
            forces = sim.get_forces(positions)

        dt, mass = sim.dt, sim._mass
        new_positions = positions + velocities * dt + forces / (2 * mass) * dt ** 2
        new_positions = sim._recenter(new_positions)

        new_forces = sim.get_forces(new_positions)
        new_velocities = velocities + dt * (forces + new_forces) / (2 * mass)

        return new_positions, new_velocities, new_forces


class PositionVerletStepPredictor(MDStepPredictor):
    """
    Position-Verlet: velocity-centered dual of velocity-Verlet. Same order of
    accuracy and symplectic properties, different bookkeeping (updates
    velocity at the half-step, position twice per step). Costs the same one
    force eval/step as velocity-Verlet; included mostly for completeness /
    cases where you specifically want position updates centered on the force
    evaluation point.
    """
    name = 'position-verlet'

    def predict_step(self, positions, velocities, forces, sim=None):
        if forces is None:
            forces = sim.get_forces(positions)

        dt, mass = sim.dt, sim._mass
        mid_positions = sim._recenter(positions + velocities * dt / 2)

        mid_forces = sim.get_forces(mid_positions)
        new_velocities = velocities + mid_forces / mass * dt

        new_positions = sim._recenter(mid_positions + new_velocities * dt / 2)
        new_forces = mid_forces  # force at the midpoint stands in for this step

        return new_positions, new_velocities, new_forces


class BeemanStepPredictor(MDStepPredictor):
    """
    Beeman's algorithm: uses the previous step's forces in addition to the
    current ones for a more accurate velocity estimate. NOT symplectic --
    energy drifts more than velocity-Verlet over long trajectories -- but
    gives better instantaneous velocities/temperature, which can matter if
    you're feeding those into a sensitive thermostat. Needs to persist the
    prior force between calls, so state lives on the instance.
    """
    name = 'beeman'

    def __init__(self, **opts):
        super().__init__(**opts)
        self._f_prev = None  # f(t - dt); bootstrapped to f(t) on the first call

    def predict_step(self, positions, velocities, forces, sim=None):
        if forces is None:
            forces = sim.get_forces(positions)
        f_prev = self._f_prev if self._f_prev is not None else forces

        dt, mass = sim.dt, sim._mass
        new_positions = positions + velocities * dt + (4 * forces - f_prev) * dt ** 2 / (6 * mass)
        new_positions = sim._recenter(new_positions)

        new_forces = sim.get_forces(new_positions)
        new_velocities = velocities + (2 * new_forces + 5 * forces - f_prev) * dt / (6 * mass)

        self._f_prev = forces
        return new_positions, new_velocities, new_forces

class InternalCoordinateStepPredictor(MDStepPredictor):
    """
    As before, but with `include_curvature=True` this adds the
    -1/2 qdot^T (dG/dq) qdot term arising from the curvilinear metric,
    using order-2 Jacobians (d^2q/dx^2) already available via the same
    `.jacobian(...)` call used for the B-matrix.

    NOTE ON SHAPES: this assumes the same square-ish (n_cartesian ==
    n_internal) convention your existing get_forces() reshape uses. If your
    `self._internals` setup uses a genuinely redundant or reduced internal
    set, the reshapes below will need adjusting to match however
    `ccs(coords).jacobian(ics, order=2, ...)` actually shapes its output --
    I haven't been able to verify that against a live run, so treat the
    tensor-shape bookkeeping here as a first draft to validate, not gospel.
    A good validation step: finite-difference G(q) numerically and compare
    to `_dG_dq` before trusting this in production.
    """
    name = 'internal-verlet'

    def __init__(self,
                 include_curvature=False, **opts):
        super().__init__(**opts)
        self.include_curvature = include_curvature

    def _get_internals(self, sim, coords):
        ccs, ics = sim._internals
        icoords = ccs.convert_coords(coords, ics)
        if not isinstance(icoords, np.ndarray):
            icoords = icoords[0]
        return icoords

    def _get_B(self, sim, coords):
        ccs, ics = sim._internals
        B = ccs(coords).jacobian(ics, order=1, all_numerical=True)[0]
        ncs = np.prod(coords.shape[-2:], dtype=int)
        return np.reshape(np.moveaxis(B, 1, 0), (len(coords), -1, ncs))

    def _get_dBdx(self, sim, coords):
        # second Cartesian derivative of the internals: d^2 q_m / dx_j dx_l
        # shape -> (n_traj, n_internal, n_cartesian, n_cartesian)
        ccs, ics = sim._internals
        d2q = ccs(coords).jacobian(ics, order=2, all_numerical=True)[0]
        ncs = np.prod(coords.shape[-2:], dtype=int)
        n_traj = len(coords)
        # mirror the order=1 moveaxis/reshape convention, extended one axis --
        # VERIFY this against your actual jacobian() output shape before trusting it
        d2q = np.moveaxis(d2q, 2, 0)  # bring n_traj to front, if it sits at axis 2 like order=1's axis 1
        return np.reshape(d2q, (n_traj, -1, ncs, ncs))

    def _mass_weighted_pinv(self, sim, B):
        m_inv = 1 / sim._mass.reshape(len(B), -1)
        Minv_Bt = B.transpose(0, 2, 1) * m_inv[:, :, np.newaxis]
        G = B @ Minv_Bt
        G_inv = np.linalg.pinv(G)
        return Minv_Bt @ G_inv, G, G_inv

    def _dG_dq(self, sim, coords, B, A):
        # dB_mj/dq_k = sum_l  d2q[m,j,l] * A[l,k]
        C = self._get_dBdx(sim, coords)             # (n_traj, n_int_m, n_cart_j, n_cart_l)
        dB_dq = np.einsum('tmjl,tlk->tmjk', C, A)    # (n_traj, n_int_m, n_cart_j, n_int_k)

        m_inv = 1 / sim._mass.reshape(len(B), -1)    # (n_traj, n_cart_j)
        # dG_mn/dq_k = sum_j [ dB_mj/dq_k * m_inv_j * B_nj  +  B_mj * m_inv_j * dB_nj/dq_k ]
        term1 = np.einsum('tmjk,tj,tnj->tmnk', dB_dq, m_inv, B)
        term2 = np.einsum('tmj,tj,tnjk->tmnk', B, m_inv, dB_dq)
        return term1 + term2                          # (n_traj, n_int_m, n_int_n, n_int_k)

    def _curvature(self, dG_dq, qdot):
        # curvature_m = -1/2 * sum_{n,k} qdot_n * dG_mn/dq_k * qdot_k
        return -0.5 * np.einsum('tmnk,tn,tk->tm', dG_dq, qdot, qdot)

    def _back_transform(self, sim, x0, q_target):
        x = x0
        for _ in range(self.max_iters):
            residual = q_target - self._get_internals(sim, x)
            if np.max(np.abs(residual)) < self.back_transform_tol:
                break
            B = self._get_B(sim, x)
            A, _, _ = self._mass_weighted_pinv(sim, B)
            x = sim._recenter(x + (A @ residual.reshape(len(x), -1, 1)).reshape(x.shape))
        return x

    def _internal_accel(self, sim, coords, velocities, forces):
        B = self._get_B(sim, coords)
        A, G, G_inv = self._mass_weighted_pinv(sim, B)
        q = self._get_internals(sim, coords)
        qdot = (B @ velocities.reshape(len(coords), -1, 1)).reshape(q.shape)
        f_int = (B @ (forces / sim._mass).reshape(len(coords), -1, 1)).reshape(q.shape)

        accel = f_int
        if self.include_curvature:
            dG_dq = self._dG_dq(sim, coords, B, A)
            accel = accel + self._curvature(dG_dq, qdot)

        return q, qdot, f_int, accel, B, A

    def predict_step(self, positions, velocities, forces, sim=None):
        if forces is None:
            forces = sim.get_forces(positions)

        q, qdot, f_int, accel, B, A = self._internal_accel(sim, positions, velocities, forces)

        dt = sim.dt
        q_new = q + qdot * dt + accel / 2 * dt ** 2
        new_positions = self._back_transform(sim, positions, q_new)

        new_forces = sim.get_forces(new_positions)
        _, qdot_new_placeholder, f_int_new, accel_new, B_new, A_new = self._internal_accel(
            sim, new_positions, velocities, new_forces  # velocities here is a placeholder for qdot estimate below
        )
        qdot_new = qdot + (accel + accel_new) / 2 * dt

        new_velocities = (A_new @ qdot_new.reshape(len(positions), -1, 1)).reshape(positions.shape)

        return new_positions, new_velocities, new_forces

class AIMDSimulator:
    __props__ = (
        "atomic_structures",
        "internals",
        "velocities",
        "track_kinetic_energy",
        "track_velocities",
        "timestep",
        "sampling_rate",
        "thermostat",
        "target_temperature"
    )
    def __init__(self,
                 atoms,
                 coordinates,
                 force_function,
                 atomic_structures=None,
                 internals=None,
                 velocities=0,
                 track_kinetic_energy=False,
                 track_velocities=False,
                 timestep=.1,
                 sampling_rate=1,
                 step_predictor=None,  # None | str | dict | MDStepPredictor instance
                 thermostat=None,           # None | str | dict | MDThermostat instance
                 target_temperature=None,   # Kelvin; used as default if `thermostat` doesn't set its own
                 remove_rotation=False
                 ):
        self.masses = np.array([
                AtomData[a, "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass")
            if isinstance(a, str) else a for a in atoms
        ])
        coords = np.asanyarray(coordinates)
        if coords.ndim == 1:
            coords = coords[np.newaxis]

        if atomic_structures is None:
            atomic_structures = (
                coords.shape[-1] == 3
                    and
                coords.shape[-2] == len(self.masses)
            )
        self._atomic_structs = atomic_structures
        if self._atomic_structs:
            if coords.ndim == 2: coords = coords[np.newaxis]
            self._mass = np.expand_dims(self.masses, -1)
            for _ in range(coords.ndim - 2):
                self._mass = np.expand_dims(self._mass, 0)
            # self._mass = self.masses[np.newaxis, :,  np.newaxis]
        else: # regular coordinates
            self._mass = self.masses[np.newaxis]
            for _ in range(coords.ndim - 2):
                self._mass = np.expand_dims(self._mass, 0)
        self.coords = coords
        # for _ in range(coordinates.ndim - 2):
        #     self._mass = np.expand_dims(self._mass, -1)
        if isinstance(velocities, (int, float, np.integer, np.floating)):
            velocities = np.full_like(self.coords, velocities)
        velocities = np.asanyarray(velocities)
        self.velocities = velocities
        if internals is not None:
            if not isinstance(internals, MolecularZMatrixCoordinateSystem):
                base_mol = Molecule(
                    ['H']*len(self.masses),
                    coords=self.coords,
                    masses=self.masses,
                    internals=internals
                )
                internals = [base_mol.coords.system, base_mol.internal_coordinates.system]
            else:
                internals = [internals.molecule.coords.system, internals]
        self._internals = internals
        self.force_function = force_function
        self._prev_forces = None

        self.trajectory = collections.deque()
        self.trajectory.append(self.coords)
        if track_kinetic_energy:
            self.kinetic_energies = collections.deque()
            self.kinetic_energies.append(
                self._ke(self._mass, self.velocities)
            )
        else:
            self.kinetic_energies = None
        if track_velocities:
            self.velocity_deque = collections.deque()
            self.velocity_deque.append(self.velocities)
        else:
            self.velocity_deque = None
        self.steps = 0
        self.dt = timestep
        self.sampling_rate = sampling_rate

        self.thermostat = MDThermostat.from_spec(thermostat, target_temperature=target_temperature)
        self.target_temperature = (
            self.thermostat.target_temperature if self.thermostat is not None else target_temperature
        )
        self.step_predictor = MDStepPredictor.from_spec(step_predictor)

        if self._atomic_structs:
            n_dof = 3 * len(self.masses) - 3  # COM removed
            if remove_rotation:
                n_dof -= 3
        else:
            n_dof = self.coords.shape[-1] * self.coords.shape[-2] if self.coords.ndim > 2 else self.coords.shape[-1]
        self._n_dof = max(n_dof, 1)

    def _temperature(self, vels):
        ke = self._ke(self._mass, vels)
        return 2 * ke / (self._n_dof * UnitsData.convert("Kelvins", "Hartrees"))

    @classmethod
    def mode_energies_to_velocities(cls, modes, masses, energy_splits, inverse=None):

        masses = np.asanyarray(masses)
        e_part = np.asanyarray(energy_splits)
        smol = e_part.ndim == 1
        if smol: e_part = e_part[np.newaxis]

        # if not np.allclose(np.linalg.norm(modes, axis=0), np.ones(modes.shape[1])):  # i.e. came out of a generalized eigenvalue run
        #     modes = modes.reshape(-1, 3, modes.shape[1]) * np.sqrt(masses[:, np.newaxis, np.newaxis])
        #     modes = modes.reshape(-1, modes.shape[-1])

        if inverse is None:
            if not np.allclose(np.linalg.norm(modes, axis=0), np.ones(modes.shape[1])):
                raise ValueError("modes aren't normalized -> non-unitary transformation")
            inverse = modes.T

        tripmass = np.repeat(masses, 3).flatten()
        g = inverse.T @ np.diag(tripmass) @ inverse

        gv, Q = np.linalg.eigh(g)
        sorting = np.argsort(np.argmax(Q ** 2, axis=0))
        gv = gv[sorting]  # maximum similarity to OG vectors
        Q = Q[:, sorting]

        axes = (Q.T @ modes).reshape(modes.shape[0], -1, 3)
        v_part = np.sign(e_part) * np.sqrt(2 * np.abs(e_part)) # / gv[np.newaxis])
        vels = np.sum(
            axes[np.newaxis] *
                v_part[:, :, np.newaxis, np.newaxis],
            axis=1
        )

        return vels

    #TODO: add internal coordinate propagation following the Krimm paper using the B-matrix
    #      or not...maybe I just need a faster Jacobian?
    def get_forces(self, coords):
        if self._internals is not None:
            ccs, ics = self._internals
            icoords = ccs.convert_coords(coords, ics)
            if not isinstance(icoords, np.ndarray):
                icoords = icoords[0] # converter info
            forces = self.force_function(icoords)
            jacs = ccs(coords).jacobian(ics, order=1, all_numerical=True)[0] #dRdX?
            ncs = np.prod(coords.shape[-2:], dtype=int)
            jacs = np.reshape(np.moveaxis(jacs, 1, 0), (len(coords), ncs, ncs))
            forces = jacs@forces.reshape((len(coords), ncs, 1))
            forces = forces.reshape(coords.shape)
            # raise Exception(forces.shape)
        else:
            forces = self.force_function(coords)
        return forces

    def _recenter(self, coords):
        if self._atomic_structs:
            com = np.tensordot(self.masses, coords, axes=[0, -2]) / np.sum(self._mass)
            coords = coords - com[:, np.newaxis, :]
        return coords

    def step(self):
        forces = self._prev_forces
        coords, vels, forces_new = self.step_predictor.predict_step(
            self.coords, self.velocities, forces, sim=self
        )

        if self.thermostat is not None:
            coords, vels = self.thermostat.apply_thermostat(coords, vels, sim=self)

        self._prev_forces = forces_new
        self.velocities = vels
        self.coords = coords
        self.steps += 1

        return coords, vels, forces_new

    @staticmethod
    def _ke(_mass, v):
        # _mass and v are 2D b.c. atomic structs
        return 1 / 2 * np.sum(np.sum(_mass * v ** 2, axis=-1), axis=-1)
    def propagate(self, num_steps=1):

        for _ in range(num_steps):
            c, v, f = self.step()
            if self.steps % self.sampling_rate == 0:
                if self.kinetic_energies is not None:
                    self.kinetic_energies.append(self._ke(self._mass, v))
                self.trajectory.append(c)
                if self.velocity_deque is not None:
                    self.velocity_deque.append(v)

        return self.trajectory

    def build_interpolation(self, energy_function, interpolation_order=2,
                            equilibration_steps=None,
                            interpolator_class=None, eckart_embed=True, reference=None, **interpolator_options):

        traj = np.array(self.trajectory)
        if equilibration_steps is not None:
            traj = traj[equilibration_steps:]
        traj = traj.reshape((-1,) + traj.shape[-2:])

        vals = [energy_function(traj)]

        # if clustering_criterion > ...:
        #     ...

        # TODO: add clustering radius + energy cutoff within disks

        if interpolation_order > 0:
            vals.append(-self.force_function(traj))
        if interpolation_order > 1:
            cshape = (-1,) + traj.shape[-2:]
            npts = np.prod(traj.shape[-2:], dtype=int)
            grad = lambda x:-self.force_function(x.reshape(cshape)).reshape(x.shape)
            hess_fun = FiniteDifferenceDerivative(grad, function_shape=(npts, npts))
            hess = np.moveaxis(
                hess_fun.derivatives(traj.reshape(-1, npts)).derivative_tensor([1])[0],
                1, 0
            )
            vals.append(hess)
        if interpolation_order > 2:
            raise ValueError("can only do order 2 interps for now because I am laaazy")

        if eckart_embed:
            if isinstance(reference, Molecule):
                ref = reference
            elif reference is None:
                ref_pos = np.argmin(vals[0])
                ref = Molecule(
                    ["H"] * len(self.masses),
                    traj[ref_pos],
                    masses=self.masses,
                )
            else:
                ref = Molecule(
                    ["H"] * len(self.masses),
                    reference,
                    masses=self.masses,
                )
            ref = ref.get_embedded_molecule(load_properties=False)

            # all_crds = []
            # all_ders = [[] for _ in range(interpolation_order)]

            rots, _, (pax_traj, _, pax_rots) = ref.get_embedding_data(traj)
            traj = pax_traj @ np.swapaxes(rots, -2, -1)

            if interpolation_order > 0:
                npts = np.prod(traj.shape[-2:], dtype=int)
                cshape = traj.shape[-2:]

                # new_derivs = [vals[1] @ pax_rots @ np.swapaxes(rots, -2, -1)]
                # raise Exception(new_derivs[0])
                new_derivs = PropertyManager._transform_derivatives(
                    [v.reshape((-1,) + (npts,)*(i+1)) for i,v in enumerate(vals[1:])],
                    rots @ pax_rots
                )
                vals = [vals[0]] + [v.reshape((-1,) + cshape*(i+1)) for i,v in enumerate(new_derivs)]


            # for crd_ders in zip(traj, *vals[1:]): #TODO: speed this up...
            #     c = crd_ders[0]
            #     mol = Molecule(
            #         ["H"] * len(self.masses),
            #         c,
            #         masses=self.masses,
            #         potential_derivatives=crd_ders[1:]
            #     ).get_embedded_molecule(ref, load_properties=False)
            #     all_crds.append(mol.coords)
            #     pes = mol.potential_derivatives
            #     for i, p in enumerate(pes):
            #         all_ders[i].append(p)

            # traj = np.array(all_crds)

        if interpolator_class is None:
            interpolator_class = RBFDInterpolator

        return interpolator_class(traj, *vals, **interpolator_options)

    def extract_trajectory(self, flatten=True, embed=True, extract_velocities=None):
        ref_struct = None
        if isinstance(embed, (np.ndarray, list, tuple)):
            ref_struct = np.asanyarray(embed)
            embed = True
        base_coords = np.array(self.trajectory)
        if flatten:
            base_coords = base_coords.reshape((-1,)+base_coords.shape[-2:])
        if extract_velocities is None:
            extract_velocities = self.velocity_deque is not None
        if extract_velocities:
            velocities = np.array(self.velocity_deque)
        else:
            velocities = None
        if embed:
            if ref_struct is None:
                if not flatten:
                    ref_crds = base_coords.reshape((-1,)+base_coords.shape[-2:])
                else:
                    ref_crds = base_coords
                if self.kinetic_energies is not None:
                    # assuming symplectic integration...
                    kes = np.array(self.kinetic_energies).flatten()
                    max_pos = np.argmax(kes)
                    ref_struct = base_coords[max_pos]
                else:
                    ref_struct = ref_crds[0]
            mol = Molecule(["H"]*len(self.masses), ref_struct, masses=self.masses)
            if extract_velocities:
                embedding = mol.get_embedding_data(base_coords)
                coords = np.tensordot(
                    nput.vec_tensordot(
                        embedding.coord_data.coords,
                        embedding.rotations,
                        shared=1,
                        axes=[-1, 2]
                    ),
                    embedding.reference_data.axes,
                    axes=[-1, 1]
                )
                velocities = np.tensordot(
                    nput.vec_tensordot(
                        np.reshape(velocities, coords.shape),  # dx/dt has same embedding
                        embedding.rotations,
                        shared=1,
                        axes=[-1, 2]
                    ),
                    embedding.reference_data.axes,
                    axes=[-1, 1]
                )
                coords = np.reshape(coords, base_coords.shape)
                velocities = np.reshape(velocities, base_coords.shape)
            else:
                coords = mol.embed_coords(base_coords)
        else:
            coords = base_coords
        if extract_velocities:
            return coords, velocities
        else:
            return coords



class PairwisePotential:
    def __init__(self, fun, deriv=None):
        self.fun = fun
        if deriv is None:
            if hasattr(deriv, 'deriv'):
                fun.deriv()
            else:
                deriv = self.numerical_deriv(fun)
        self.deriv = deriv
        self._inds = None # to cache indices when applied to same-shape data

    @classmethod
    def numerical_deriv(cls, fun, step_size=.01):
        @functools.wraps(fun)
        def deriv(dists):
            return (fun(dists + step_size) - fun(dists - step_size)) / (2 * step_size)
        return deriv

    def _get_inds_cached(self, coords):
        if self._inds is not None and len(coords) != self._inds[0]:
            self._inds = None
        if self._inds is None:
            row_inds, col_inds = np.triu_indices(len(coords), k=1)
            _, row_map = np.unique(row_inds, return_index=True) # One element less than the coords
            sorting = np.argsort(col_inds)
            _, col_map = np.unique(col_inds, return_index=True) # One element less than the coords
            self._inds = [len(coords), (row_inds, col_inds, row_map, sorting, col_map)]
        return self._inds[1]

    def eval(self, coords):
        row_inds, col_inds, row_map, col_sort, col_map = self._get_inds_cached(coords)
        pot = self.fun

        # compute unsigned force for each pair
        diffs = coords[row_inds] - coords[col_inds]
        dists = np.linalg.norm(diffs, axis=1)
        pot_vals = pot(dists)

        # print(row_map)

        # now reduce over different chunks for the positive contribs...
        pos_chunks = np.add.reduceat(
            pot_vals,
            row_map
        )
        neg_chunks = np.add.reduceat(
            pot_vals[col_sort],
            col_map
        )

        return np.concatenate(
            [
                pos_chunks[:1],
                pos_chunks[1:] + neg_chunks[:-1],
                pos_chunks[-1:]
            ],
            axis=0
        )

    def forces(self, coords):
        row_inds, col_inds, row_map, col_sort, col_map = self._get_inds_cached(coords)
        pot_deriv = self.deriv

        # compute unsigned force for each pair
        diffs = coords[row_inds] - coords[col_inds]
        dists = np.linalg.norm(diffs, axis=1)
        base_force = pot_deriv(dists)
        normals = diffs / dists[:, np.newaxis]
        force_list = -normals * base_force[:, np.newaxis]

        # now reduce over different chunks for the positive contribs...
        pos_chunks = np.add.reduceat(
            force_list,
            row_map
        )
        neg_chunks = np.add.reduceat(
            force_list[col_sort],
            col_map
        )

        return np.concatenate(
            [
                pos_chunks[:1],
                pos_chunks[1:] - neg_chunks[:-1],
                pos_chunks[-1:]
            ],
            axis=0
        )

    def __call__(self, coords):
        return self.eval(coords)
