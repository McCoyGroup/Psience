"""Cartesian RMSD metadynamics regressions; no electronic-structure backend needed.

The reference energy uses Horn's quaternion eigenvalue method, independent of
the production alignment. Forces are checked by central energy differences.
Run normally with unittest after merging, or use the adjacent draft runner.
"""
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose

import Psience.AIMD as aimd
from Psience.AIMD.Simulator import AIMDSimulator, RMSDBiasPotential


REF = np.array([[0., 0., 0.], [1.4, .2, -.1], [.1, 1.1, .3], [.2, -.1, 1.5]])
QUERY = REF + np.array([[.03, -.05, .02], [.2, .04, .06], [-.1, .15, -.03], [.04, -.07, .12]])
MASSES = np.array([12., 1., 16., 2.])


def reference_rmsd_squared(x, h, weights=None):
    """Proper-rotation least squares via the largest quaternion eigenvalue."""
    w = np.ones(len(x)) if weights is None else np.asarray(weights, dtype=float)
    w = w / w.sum()
    x = x - np.sum(w[:, None] * x, axis=0)
    h = h - np.sum(w[:, None] * h, axis=0)
    s = (w[:, None] * x).T @ h
    xx, xy, xz, yx, yy, yz, zx, zy, zz = s.ravel()
    horn = np.array([
        [xx+yy+zz, yz-zy, zx-xz, xy-yx],
        [yz-zy, xx-yy-zz, xy+yx, zx+xz],
        [zx-xz, xy+yx, -xx+yy-zz, yz+zy],
        [xy-yx, zx+xz, yz+zy, -xx-yy+zz],
    ])
    d2 = np.sum(w[:, None] * (x*x + h*h)) - 2 * np.linalg.eigvalsh(horn)[-1]
    return max(float(d2), 0.)


def reference_energy(x, history, height=.03, width=.4, weights=None):
    return sum(height * np.exp(-reference_rmsd_squared(x, h, weights) / (2*width**2))
               for h in history)


def numerical_force(x, history, height=.03, width=.4, weights=None):
    force = np.zeros_like(x, dtype=float)
    eps = 2.e-5
    for index in np.ndindex(x.shape):
        plus, minus = x.copy(), x.copy()
        plus[index] += eps
        minus[index] -= eps
        force[index] = -(reference_energy(plus, history, height, width, weights)
                         - reference_energy(minus, history, height, width, weights)) / (2*eps)
    return force


def rotation():
    axis = np.array([.3, -.7, .5])
    axis /= np.linalg.norm(axis)
    skew = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
    theta = 1.17
    return np.eye(3) + np.sin(theta)*skew + (1-np.cos(theta))*(skew @ skew)


class RMSDBiasRegressionTests(unittest.TestCase):
    def make_bias(self, **kwargs):
        opts = dict(k=5, height=.03, width=.4, mass_weighted=False)
        opts.update(kwargs)
        return RMSDBiasPotential(np.zeros_like, **opts)

    def test_public_export(self):
        self.assertIs(getattr(aimd, 'RMSDBiasPotential', None), RMSDBiasPotential)

    def test_empty_history_preserves_base_force_and_input(self):
        coords = REF[None].copy()
        saved = coords.copy()
        base = np.full(coords.shape, .7)
        bias = RMSDBiasPotential(lambda c: base, k=2)
        result = bias(coords)
        assert_allclose(result, base)
        assert_allclose(coords, saved)
        self.assertFalse(np.shares_memory(result, base))

    def test_second_call_and_repeated_identical_structure(self):
        bias = self.make_bias()
        for _ in range(8):
            assert_allclose(bias(REF[None]), 0., atol=2.e-15)

    def test_force_is_negative_energy_gradient(self):
        bias = self.make_bias()
        expected = numerical_force(QUERY, [REF])
        actual = bias._bias_force(QUERY, [REF], bias.get_weights(QUERY[None]))
        assert_allclose(actual, expected, atol=2.e-10, rtol=2.e-7)

    def test_mass_weighted_force_is_negative_energy_gradient(self):
        bias = self.make_bias(masses=MASSES, mass_weighted=True)
        expected = numerical_force(QUERY, [REF], weights=MASSES)
        actual = bias._bias_force(QUERY, [REF], bias.get_weights(QUERY[None]))
        assert_allclose(actual, expected, atol=2.e-10, rtol=2.e-7)

    def test_mass_units_do_not_change_bias(self):
        b1 = self.make_bias(masses=MASSES, mass_weighted=True)
        b2 = self.make_bias(masses=MASSES*1822.888, mass_weighted=True)
        f1 = b1._bias_force(QUERY, [REF], b1.get_weights(QUERY[None]))
        f2 = b2._bias_force(QUERY, [REF], b2.get_weights(QUERY[None]))
        assert_allclose(f1, f2, atol=1.e-14)

    def test_supplied_masses_do_not_change_unweighted_bias(self):
        bias = self.make_bias(masses=MASSES)
        actual = bias._bias_force(QUERY, [REF], bias.get_weights(QUERY[None]))
        assert_allclose(actual, numerical_force(QUERY, [REF]), atol=2.e-10, rtol=2.e-7)

    def test_default_matches_unweighted_cartesian_rmsd(self):
        bias = RMSDBiasPotential(np.zeros_like, masses=MASSES)
        self.assertFalse(bias.mass_weighted)

    def test_rigid_translation_and_rotation_have_zero_force(self):
        bias = self.make_bias()
        x = REF @ rotation() + [3., -2., .8]
        actual = bias._bias_force(x, [REF], bias.get_weights(x[None]))
        assert_allclose(actual, 0., atol=2.e-14)

    def test_force_rotates_with_query_frame(self):
        bias = self.make_bias()
        rot = rotation()
        w = bias.get_weights(QUERY[None])
        f = bias._bias_force(QUERY, [REF], w)
        transformed = bias._bias_force(QUERY @ rot + [2., .4, -.7], [REF], w)
        assert_allclose(transformed, f @ rot, atol=2.e-14)

    def test_reference_frame_does_not_change_force(self):
        bias = self.make_bias()
        w = bias.get_weights(QUERY[None])
        actual = bias._bias_force(QUERY, [REF @ rotation() + [1., 2., 3.]], w)
        assert_allclose(actual, numerical_force(QUERY, [REF]), atol=2.e-10, rtol=2.e-7)

    def test_bias_adds_no_net_force_or_torque(self):
        for weighted in (False, True):
            with self.subTest(weighted=weighted):
                bias = self.make_bias(masses=MASSES, mass_weighted=weighted)
                f = bias._bias_force(QUERY, [REF], bias.get_weights(QUERY[None]))
                assert_allclose(f.sum(axis=0), 0., atol=2.e-14)
                assert_allclose(np.cross(QUERY, f).sum(axis=0), 0., atol=2.e-14)

    def test_planar_geometry_force(self):
        ref = np.array([[0., 0., 0.], [1.7, 0., 0.], [-.5, 1.6, 0.]])
        query = ref.copy()
        query[1, 0] += .12
        bias = self.make_bias()
        actual = bias._bias_force(query, [ref], bias.get_weights(query[None]))
        assert_allclose(actual, numerical_force(query, [ref]), atol=2.e-10, rtol=2.e-7)

    def test_diatomic_bias_repels_and_uses_per_atom_rmsd(self):
        ref = np.array([[-.5, 0., 0.], [.5, 0., 0.]])
        query = ref * 1.2
        bias = self.make_bias()
        f = bias._bias_force(query, [ref], bias.get_weights(query[None]))
        # RMSD = .1, and each atom carries weight 1/2.
        expected = 2*bias.alpha * .03*np.exp(-bias.alpha*.1**2) * .1/2
        assert_allclose(f[:, 0], [-expected, expected], atol=1.e-14)
        assert_allclose(f[:, 1:], 0., atol=1.e-14)

    def test_single_atom_has_no_conformational_bias_force(self):
        bias = self.make_bias()
        w = bias.get_weights(np.zeros((1, 1, 3)))
        assert_allclose(bias._bias_force(np.array([[4., 2., 1.]]), [np.zeros((1, 3))], w), 0.)

    def test_chiral_reflection_is_not_an_allowed_alignment(self):
        query = REF * [-1., 1., 1.]
        bias = self.make_bias()
        expected = numerical_force(query, [REF])
        self.assertGreater(np.linalg.norm(expected), 1.e-4)
        assert_allclose(bias._bias_force(query, [REF], bias.get_weights(query[None])),
                        expected, atol=2.e-10, rtol=2.e-7)

    def test_multi_hill_force_is_additive(self):
        bias = self.make_bias()
        refs = [REF, REF*1.05, QUERY*.94]
        actual = bias._bias_force(QUERY, refs, bias.get_weights(QUERY[None]))
        assert_allclose(actual, numerical_force(QUERY, refs), atol=5.e-10, rtol=3.e-7)

    def test_independent_history_stores_every_walker(self):
        coords = np.array([REF, QUERY, REF*1.1])
        bias = self.make_bias(k=2)
        # Isolate deposition from alignment failures on the unpatched code.
        with patch.object(bias, '_bias_force', side_effect=lambda x, h, w: np.zeros_like(x)):
            bias(coords)
            assert_allclose(bias._histories[0], coords)

    def test_independent_ring_buffer_wraps_when_walkers_differ_from_k(self):
        bias = self.make_bias(k=2)
        snapshots = [np.array([REF*(1+j*.02), QUERY*(1+j*.03), REF*(1+j*.04)]) for j in range(5)]
        with patch.object(bias, '_bias_force', side_effect=lambda x, h, w: np.zeros_like(x)):
            for coords in snapshots:
                bias(coords)
        assert_allclose(bias._histories[0], snapshots[-1])
        assert_allclose(bias._histories[1], snapshots[-2])

    def test_shared_history_pools_all_walkers(self):
        coords = np.array([REF, QUERY, REF*1.1])
        bias = self.make_bias(k=5, shared_history=True)
        bias(coords)
        assert_allclose(bias._histories[:3], coords)

    def test_shared_ring_capacity_is_number_of_structures(self):
        bias = self.make_bias(k=2, shared_history=True)
        coords = np.array([REF, QUERY, REF*1.1])
        with patch.object(bias, '_bias_force', side_effect=lambda x, h, w: np.zeros_like(x)):
            bias(coords)
        assert_allclose(bias._histories[0], coords[2])
        assert_allclose(bias._histories[1], coords[1])

    def test_stride_counts_calls_in_standalone_wrapper(self):
        bias = self.make_bias(stride=3)
        for _ in range(2):
            assert_allclose(bias(REF[None]), 0.)
        bias(REF[None])
        actual = bias(QUERY[None])[0]
        assert_allclose(actual, numerical_force(QUERY, [REF]), atol=2.e-10, rtol=2.e-7)

    def test_reset_can_reuse_instance_with_different_atom_count(self):
        bias = self.make_bias()
        bias(REF[None])
        bias.reset()
        self.assertEqual(len(bias.get_weights(REF[None, :3])), 3)
        assert_allclose(bias(REF[None, :3]), 0.)

    def test_integer_coordinates_get_floating_point_forces(self):
        bias = self.make_bias()
        ref = np.array([[-1, 0, 0], [1, 0, 0]])
        query = ref * 2
        bias(ref[None])
        f = bias(query[None])
        self.assertTrue(np.issubdtype(f.dtype, np.floating))
        self.assertGreater(f[0, 1, 0], 0.)

    def test_invalid_parameters_are_rejected(self):
        invalid = [dict(k=0), dict(k=-1), dict(k=1.5), dict(k=True),
                   dict(stride=0), dict(stride=1.2), dict(stride=True),
                   dict(width=0), dict(width=-.1), dict(width=np.nan),
                   dict(height=-.1), dict(height=np.inf),
                   dict(masses=[1., -1.]), dict(masses=[1., np.nan])]
        for opts in invalid:
            with self.subTest(opts=opts), self.assertRaises(ValueError):
                self.make_bias(**opts)


class RMSDBiasEvaluationTests(unittest.TestCase):
    def make_bias(self, **kwargs):
        opts = dict(k=5, height=.03, width=.4)
        opts.update(kwargs)
        return RMSDBiasPotential(np.zeros_like, **opts)

    def test_bias_energy_matches_independent_quaternion_reference(self):
        for weighted in (False, True):
            bias = self.make_bias(masses=MASSES, mass_weighted=weighted)
            bias.deposit(REF[None])
            expected = reference_energy(QUERY, [REF], weights=MASSES if weighted else None)
            assert_allclose(bias.bias_energy(QUERY[None]), [expected], atol=2.e-14)

    def test_energy_and_force_queries_do_not_deposit(self):
        bias = self.make_bias()
        bias.deposit(REF[None])
        history = bias._histories.copy()
        count, pointer = bias._call_count, bias._history_pointer
        for _ in range(4):
            bias.forces(QUERY[None])
            bias.bias_energy(QUERY[None])
        assert_allclose(bias._histories, history)
        self.assertEqual((bias._call_count, bias._history_pointer), (count, pointer))

    def test_deposited_geometry_is_copied(self):
        bias = self.make_bias()
        coords = REF[None].copy()
        bias.deposit(coords)
        coords[:] = 20.
        assert_allclose(bias.bias_energy(REF[None]), [.03], atol=1.e-14)

    def test_independent_walkers_match_separate_instances(self):
        bias = self.make_bias(k=2)
        singles = [self.make_bias(k=2), self.make_bias(k=2)]
        for step in range(6):
            coords = np.array([REF*(1+step*.03), QUERY*(1-step*.015)])
            expected = np.concatenate([b(c[None]) for b, c in zip(singles, coords)])
            assert_allclose(bias(coords), expected, atol=1.e-14)

    def test_shared_bias_includes_other_walkers(self):
        bias = self.make_bias(shared_history=True)
        bias.deposit(np.array([REF, QUERY]))
        query = REF*1.07
        expected = numerical_force(query, [REF, QUERY])
        assert_allclose(bias.forces(query[None])[0], expected, atol=5.e-10, rtol=3.e-7)

    def test_shape_errors_do_not_advance_or_write_history(self):
        bias = self.make_bias(masses=MASSES)
        bias.deposit(REF[None])
        for coords in (REF, np.zeros((0, 4, 3)), np.zeros((1, 3, 3)),
                       np.full((1, 4, 3), np.nan), np.zeros((2, 4, 3))):
            with self.subTest(shape=coords.shape), self.assertRaises(ValueError):
                bias(coords)
        self.assertEqual(bias._call_count, 0)
        assert_allclose(bias._histories[0], REF[None])

    def test_bad_base_force_shape_does_not_deposit(self):
        bias = RMSDBiasPotential(lambda c: np.zeros(3))
        with self.assertRaises(ValueError):
            bias(REF[None])
        self.assertEqual(bias._call_count, 0)
        self.assertIsNone(bias._histories)


class AIMDMetadynamicsIntegrationTests(unittest.TestCase):
    def simulation(self, stride=2, k=8, predictor='velocity-verlet', shared=False):
        bias = RMSDBiasPotential(lambda c: -.02*c, k=k, height=.01, width=.3,
                                 stride=stride, shared_history=shared)
        coords = np.array([REF, QUERY]) if shared else REF[None].copy()
        velocities = np.broadcast_to((QUERY-REF)*.03, coords.shape).copy()
        sim = AIMDSimulator(np.ones(4), coords, bias, velocities=velocities,
                            timestep=.03, step_predictor=predictor,
                            track_velocities=True, track_kinetic_energy=True)
        return sim, bias

    def test_force_probe_does_not_change_deposition_schedule(self):
        sim, bias = self.simulation()
        for _ in range(4):
            sim.get_forces(sim.coords)
        self.assertEqual(bias._call_count, 0)
        sim.propagate(3)
        self.assertEqual(bias._call_count, 3)
        self.assertEqual(bias._history_pointer, 1)
        assert_allclose(bias._histories[0], sim.trajectory[2])

    def test_schedule_counts_completed_steps_for_each_cartesian_predictor(self):
        for predictor in ('velocity-verlet', 'position-verlet', 'beeman'):
            with self.subTest(predictor=predictor):
                sim, bias = self.simulation(predictor=predictor)
                sim.propagate(5)
                self.assertEqual(bias._call_count, 5)
                self.assertEqual(bias._history_pointer, 2)
                self.assertEqual(len(sim.trajectory), 6)
                self.assertEqual(len(sim.velocity_deque), 6)
                self.assertEqual(len(sim.kinetic_energies), 6)

    def test_split_propagation_matches_single_run(self):
        whole, bwhole = self.simulation(k=2)
        split, bsplit = self.simulation(k=2)
        whole.propagate(12)
        split.propagate(5)
        split.propagate(7)
        assert_allclose(split.coords, whole.coords, atol=1.e-14)
        assert_allclose(split.velocities, whole.velocities, atol=1.e-14)
        assert_allclose(bsplit._histories, bwhole._histories, atol=1.e-14)

    def test_cached_force_is_invalidated_after_deposition_and_eviction(self):
        for shared in (False, True):
            sim, bias = self.simulation(stride=1, k=2, shared=shared)
            sim.propagate(4)
            self.assertIsNone(sim._prev_forces)
            fresh = sim.get_forces(sim.coords)
            expected = sim._recenter(sim.coords + sim.dt*sim.velocities
                                     + sim.dt**2*fresh/(2*sim._mass))
            sim.step()
            assert_allclose(sim.coords, expected, atol=1.e-14)

    def test_external_deposit_invalidates_cached_force(self):
        sim, bias = self.simulation(stride=10)
        sim.step()
        bias.deposit((QUERY*1.1)[None])
        force = sim.get_forces(sim.coords)
        expected = sim._recenter(sim.coords + sim.dt*sim.velocities + sim.dt**2*force/(2*sim._mass))
        sim.step()
        assert_allclose(sim.coords, expected, atol=1.e-14)

    def test_reset_invalidates_cached_force(self):
        sim, bias = self.simulation(stride=10)
        bias.deposit((QUERY*1.1)[None])
        sim.step()
        bias.reset()
        force = sim.get_forces(sim.coords)
        expected = sim._recenter(sim.coords + sim.dt*sim.velocities + sim.dt**2*force/(2*sim._mass))
        sim.step()
        assert_allclose(sim.coords, expected, atol=1.e-14)

    def test_deposition_is_independent_of_trajectory_sampling_rate(self):
        full, bfull = self.simulation()
        sparse, bsparse = self.simulation()
        sparse.sampling_rate = 3
        full.propagate(7)
        sparse.propagate(7)
        assert_allclose(full.coords, sparse.coords, atol=1.e-14)
        assert_allclose(bfull._histories, bsparse._histories, atol=1.e-14)
        self.assertEqual(bsparse._call_count, 7)
        self.assertEqual(len(sparse.trajectory), 3)

    def test_internal_coordinate_callback_is_rejected_for_cartesian_bias(self):
        bias = RMSDBiasPotential(np.zeros_like)
        with self.assertRaisesRegex(ValueError, 'Cartesian'):
            AIMDSimulator(np.ones(4), REF, bias, internals={'zmatrix': [[0, -1, -1, -1]]*4})

    def test_zero_height_matches_unbiased_trajectory(self):
        base_force = lambda c: -.02*c
        velocities = (QUERY-REF)[None]*.03
        plain = AIMDSimulator(np.ones(4), REF, base_force, velocities=velocities, timestep=.03)
        biased = AIMDSimulator(np.ones(4), REF,
                               RMSDBiasPotential(base_force, k=2, height=0., stride=2),
                               velocities=velocities, timestep=.03)
        plain.propagate(12)
        biased.propagate(12)
        assert_allclose(biased.coords, plain.coords, atol=1.e-14)
        assert_allclose(biased.velocities, plain.velocities, atol=1.e-14)

    def test_zero_height_with_seeded_langevin_matches_unbiased(self):
        base_force = lambda c: -.02*c
        opts = dict(velocities=(QUERY-REF)[None]*.03, timestep=.03,
                    thermostat={'type': 'langevin', 'friction': .1, 'seed': 43}, target_temperature=300.)
        plain = AIMDSimulator(np.ones(4), REF, base_force, **opts)
        biased = AIMDSimulator(np.ones(4), REF,
                               RMSDBiasPotential(base_force, k=2, height=0., stride=2), **opts)
        plain.propagate(12)
        biased.propagate(12)
        assert_allclose(biased.coords, plain.coords, atol=1.e-14)
        assert_allclose(biased.velocities, plain.velocities, atol=1.e-14)

    def test_metadynamics_escapes_a_double_well_with_sub_barrier_initial_energy(self):
        # Toy diatomic PES V(r) = .64 (r-1)^2 (r-2)^2 has minima at 1 and
        # 2, with a .04 barrier at 1.5. Initial kinetic energy is only .0016.
        # This checks actual propagation and exploration, beyond force shape.
        def forces(c):
            diff = c[:, 1]-c[:, 0]
            r = np.linalg.norm(diff, axis=-1)
            grad = .64*(r-1)*(r-2)*(2*r-3)
            f = grad[:, None]*diff/r[:, None]
            return np.stack([f, -f], axis=1)

        maxima = []
        for height in (0., .01):
            bias = RMSDBiasPotential(forces, k=30, height=height, width=.12, stride=10)
            sim = AIMDSimulator([1., 1.], [[-.5, 0., 0.], [.5, 0., 0.]], bias,
                                velocities=np.array([[[-.04, 0., 0.], [.04, 0., 0.]]]), timestep=.02)
            trajectory = np.array(sim.propagate(600))
            self.assertTrue(np.all(np.isfinite(trajectory)))
            maxima.append(np.linalg.norm(trajectory[:, :, 1]-trajectory[:, :, 0], axis=-1).max())
        self.assertLess(maxima[0], 1.2)
        self.assertGreater(maxima[1], 1.8)

    def test_unbiased_interpolation_does_not_deposit_or_mix_bias_derivatives(self):
        sim, bias = self.simulation(k=2)
        sim.propagate(4)
        hist = bias._histories.copy()
        count = bias._call_count

        class CaptureInterpolator:
            def __init__(self, coords, *values, **opts):
                self.coords, self.values = coords, values

        interp = sim.build_interpolation(lambda c: .01*np.sum(c*c, axis=(-2, -1)),
                                         interpolation_order=2, eckart_embed=False,
                                         interpolator_class=CaptureInterpolator)
        assert_allclose(interp.values[1], .02*interp.coords, atol=1.e-12)
        assert_allclose(np.reshape(interp.values[2], (-1, 12, 12)),
                        np.broadcast_to(.02*np.eye(12), (len(interp.coords), 12, 12)), atol=2.e-9)
        assert_allclose(bias._histories, hist)
        self.assertEqual(bias._call_count, count)

    def test_frozen_bias_velocity_verlet_conserves_total_energy(self):
        ref = np.array([[-.5, 0., 0.], [.5, 0., 0.]])
        bias = RMSDBiasPotential(lambda c: -.1*c, height=.02, width=.3)
        bias.deposit(ref[None])
        # Passing the pure force method freezes the deposited potential.
        sim = AIMDSimulator([1., 1.], ref*1.15, bias.forces,
                            velocities=np.array([[[-.01, 0., 0.], [.01, 0., 0.]]]), timestep=.02)
        energy = []
        for _ in range(150):
            energy.append(float(.05*np.sum(sim.coords**2) + bias.bias_energy(sim.coords)[0]
                                + sim._ke(sim._mass, sim.velocities)[0]))
            sim.step()
        self.assertLess(np.ptp(energy), 2.e-7)
        self.assertEqual(bias._call_count, 0)


class MetadynamicsThermostatTests(unittest.TestCase):
    def test_langevin_noise_variance_follows_fluctuation_dissipation(self):
        from McUtils.Data import UnitsData
        # Check the distribution of one thermal kick in actual AIMD state.
        # It must satisfy Var(delta v_i) = 2 gamma k_B T dt / m_i.
        masses = [1., 2., 5., 10.]
        coords = np.broadcast_to(REF, (6000, 4, 3)).copy()
        sim = AIMDSimulator(masses, coords, np.zeros_like, timestep=.025,
                            thermostat={'type': 'langevin', 'friction': .2, 'seed': 867},
                            target_temperature=300.)
        _, velocities = sim.thermostat.apply_thermostat(sim.coords, sim.velocities, sim)
        expected = 2*.2*300.*UnitsData.convert('Kelvins', 'Hartrees')*.025 / np.array(masses)
        actual = velocities.var(axis=(0, 2))
        assert_allclose(actual/expected, 1., atol=.035)

    def test_langevin_kick_scales_with_square_root_of_timestep(self):
        kicks = []
        for timestep in (.025, .1):
            sim = AIMDSimulator(np.ones(4), REF, np.zeros_like, timestep=timestep,
                                thermostat={'type': 'langevin', 'friction': .2, 'seed': 43},
                                target_temperature=300.)
            _, v = sim.thermostat.apply_thermostat(sim.coords, sim.velocities, sim)
            kicks.append(v)
        assert_allclose(kicks[1], 2*kicks[0], atol=1.e-14)


if __name__ == '__main__':
    unittest.main()
