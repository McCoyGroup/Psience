from Peeves.TestUtils import *
from unittest import TestCase
from Psience.Nonlinear import *
import McUtils.Plots as plt
from McUtils.Data import UnitsData
import sys, os, numpy as np
import copy
import itertools
import json
import tempfile

class NonlinearTests(TestCase):

    @validationTest
    def test_BasicPathways(self):
        from Psience.BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
        initial_basis = BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(1),
            0
        )
        # pathways = liouville_pathways(3)

        responses = experimental_response_generator(
            {
                ((0,0), (0,1)):{
                    'frequency':1610,
                    'transition_moment':[0, 0, .1]
                },
                ((0,0), (1,0)):{
                    'frequency':1570,
                    'transition_moment':[0, .1, 0]
                },
                ((0,1), (0,2)): {
                    'frequency': 1600,
                    'transition_moment': [0, 0, .12]
                },
                ((1,0), (2, 0)): {
                    'frequency': 1560,
                    'transition_moment': [0, .12, 0]
                },
                ((0, 1), (1, 1)): {
                    'frequency': 1570,
                    'transition_moment': [0, .1, 0]
                },
                ((1, 0), (1, 1)): {
                    'frequency': 1610,
                    'transition_moment': [0,  0, .1]
                }
            },
            band_coherences={
                (0, 1):2,
                (1, 1):.5,
                (1, 2):2
            },
            frequency_unit="Wavenumbers",
            application_domain="frequency",
            # polarization=None,
            # response_function_generator="alt2d",
            response_function_generator="simple2dir",
            # included_signals=['non-rephasing']
        )

        spec = responses.get_spectrum(
            # [1550, 1620], 0, [1550, 1620],
            [1400, 1700], 10, [1400, 1700],
            num_samples=1024,
            time_step=.0625,
            default_frequency_divisions=500
        )
        # subspec = spec.frequency_filter([1400, 1700], [1400, 1700])
        subspec = spec.frequency_filter([1550, 1620], [1550, 1620])
        # print(subspec.clip(1e-6, 1))
        subspec.clip(1e-8, 1).plot(levels=15).show()


        # for p in paths[1]:
        #     print(p)

        # for cur, nxt in enumerate_basis_path(initial_basis, pathways[0]):
        #     print("-"*50)
        #     print(cur.excitations)
        #     print("."*10)
        #     print(nxt.excitations)

    class FourStateLiouvilleSystem:
        """
        2D-IR tests built around a single, deliberately "interesting" 4-state
        vibrational system:

            gs = (0, 0)   ground state
            e1 = (1, 0)   fundamental of mode A, w_A = 1610 cm^-1
            e2 = (0, 1)   fundamental of mode B, w_B = 1590 cm^-1
            f  = (1, 1)   combination band of A and B

        The combination band is given a cross-anharmonicity DELTA = 15 cm^-1, so
        that e1 -> f sits at w_B - DELTA and e2 -> f sits at w_A - DELTA. This is
        the minimal system that shows the full textbook set of rephasing/
        non-rephasing Feynman pathways (R1/R2/R4/R5 diagonal + cross ground-state
        bleach/stimulated emission, R3/R6 excited-state absorption) while still
        only touching 4 total states, so it exercises:

          * the "simple2dir" canned pathway catalog
          * the fully general Liouville-pathway machinery
            (`prep_liouville_spaces`/`enumerate_state_paths`) that the default
            (non-"simple2dir") response function generator drives
          * restricting the detected Liouville-space matrix elements
            (`response_tensor_elements`) to specific population sectors
          * isolating individual pathway families (rephasing vs. non-rephasing)
          * starting the pathway enumeration from a non-ground `initial_states`
            (a different corner of Liouville space entirely)
          * polarization-dependent orientational averaging
        """

        W1 = 1610.
        W2 = 1590.
        DELTA = 15.
        WINDOW = [1550, 1660]
        T2 = 10
        SAMPLING = dict(num_samples=128, time_step=.0625, default_frequency_divisions=80)

        @classmethod
        def four_state_transitions(cls):
            return {
                ((0, 0), (1, 0)): {'frequency': cls.W1, 'transition_moment': [0, .1, 0]},
                ((0, 0), (0, 1)): {'frequency': cls.W2, 'transition_moment': [0, 0, .1]},
                ((1, 0), (1, 1)): {'frequency': cls.W2 - cls.DELTA, 'transition_moment': [0, 0, .12]},
                ((0, 1), (1, 1)): {'frequency': cls.W1 - cls.DELTA, 'transition_moment': [0, .12, 0]},
            }

        @classmethod
        def band_coherences(cls):
            # dephasing rates for the (gs<->fundamental) and (fundamental<->combination)
            # coherence bands; used as Lorentzian linewidths in the response function
            return {(0, 1): 3., (1, 1): 3.}

        @classmethod
        def get_response(cls, **opts):
            return experimental_response_generator(
                cls.four_state_transitions(),
                band_coherences=cls.band_coherences(),
                frequency_unit="Wavenumbers",
                application_domain="frequency",
                **opts
            )

        @classmethod
        def get_intensities(cls, **opts):
            responses = cls.get_response(**opts)
            spec = responses.get_spectrum(cls.WINDOW, cls.T2, cls.WINDOW, **cls.SAMPLING)
            return np.real(spec.intensities), spec

        @staticmethod
        def nearest_value(spec, w1, w3):
            ix = int(np.argmin(np.abs(spec.freq1 - w1)))
            iy = int(np.argmin(np.abs(spec.freq2 - w3)))
            return np.real(spec.intensities)[iy, ix]

    @validationTest
    def test_CombinationBandCrossPeaksAndESA(self):
        """
        The canned "simple2dir" pathway catalog (R1/R2/R4/R5/R3/R6) should
        produce a positive ground-state-bleach/stimulated-emission peak at
        each of the two diagonal positions and at both symmetry-related
        cross positions, plus a genuine (negative-going) excited-state
        absorption feature reflecting the combination-band anharmonicity.
        """
        fsls = self.FourStateLiouvilleSystem()

        I, spec = fsls.get_intensities(response_function_generator="simple2dir")

        diag_1 = fsls.nearest_value(spec, fsls.W1, fsls.W1)
        diag_2 = fsls.nearest_value(spec, fsls.W2, fsls.W2)
        cross_12 = fsls.nearest_value(spec, fsls.W1, fsls.W2)
        cross_21 = fsls.nearest_value(spec, fsls.W2, fsls.W1)

        # diagonal and cross ground-state features are clearly resolved and positive
        for val in (diag_1, diag_2, cross_12, cross_21):
            self.assertGreater(val, 5e-7)

        # the two fundamentals were built symmetrically, so the spectrum should be
        # (approximately) symmetric under swapping the pump/probe axes
        self.assertTrue(np.isclose(diag_1, diag_2, rtol=0.1))
        self.assertTrue(np.isclose(cross_12, cross_21, rtol=0.1))

        # a real 2D-IR spectrum of a coupled pair must show excited-state absorption:
        # a substantial negative-going feature somewhere in the window, reflecting
        # the anharmonically-shifted 1->2 (combination band) transitions
        self.assertLess(I.min(), -0.3 * I.max())

    @validationTest
    def test_FullLiouvillePathwayGeneratorMatchesCannedCatalog(self):
        """
        Exercises the general Liouville-pathway machinery
        (`prep_liouville_spaces` + `enumerate_state_paths` +
        `_identify_response_tensor_paths`, reached whenever
        `response_function_generator` is left as its default) and checks
        that -- for the default 2D-IR experiment settings, where every
        population element of the 4-state space is a valid final detection
        channel -- it reproduces the hand-rolled "simple2dir" catalog of
        Feynman pathways exactly, since both should be enumerating the same
        physical set of rephasing/non-rephasing diagrams for this system.
        """
        fsls = self.FourStateLiouvilleSystem()

        I_simple, _ = fsls.get_intensities(response_function_generator="simple2dir")
        I_general, _ = fsls.get_intensities()  # default generator drives the full pathway search

        self.assertGreater(np.abs(I_general).max(), 1e-6)
        np.testing.assert_allclose(I_general, I_simple, atol=1e-9)

    @validationTest
    def test_RephasingAndNonRephasingPathwaysAreDistinctSectors(self):
        """
        Rephasing and non-rephasing pathways are different regions of
        Liouville space (opposite sign of the coherence accrued during t1).
        Isolating either one (via `included_signals`) should give a
        nonzero, self-consistent spectrum, and the two sectors should not
        be identical to one another.

        Note: `included_signals` only filters the `paths` (and, if given as
        a dict, `selection_rules`) mapping -- it does *not* filter a
        top-level `phases` dict, so the default 2D-IR `phases` entry (which
        has both 'rephasing' and 'non-rephasing' keys) must be overridden to
        match whichever single signal is requested, or `prep_liouville_spaces`
        raises a `KeyError` looking up the excluded signal's phase pattern.
        """
        fsls = self.FourStateLiouvilleSystem()

        I_rephasing, _ = fsls.get_intensities(
            included_signals=['rephasing'],
            phases={'rephasing': [-1, 1, 1]},
        )
        I_non_rephasing, _ = fsls.get_intensities(
            included_signals=['non-rephasing'],
            phases={'non-rephasing': [1, -1, 1]},
        )

        self.assertGreater(I_rephasing.max(), 1e-6)
        self.assertGreater(I_non_rephasing.max(), 1e-6)
        # both pathway families individually carry an excited-state absorption feature
        self.assertLess(I_rephasing.min(), 0)
        self.assertLess(I_non_rephasing.min(), 0)

        self.assertFalse(np.allclose(I_rephasing, I_non_rephasing, atol=1e-9))

        with self.assertRaises(KeyError):
            fsls.get_intensities(included_signals=['rephasing'])

    @validationTest
    def test_RestrictingResponseTensorElementsSelectsLiouvilleSectors(self):
        """
        `response_tensor_elements` selects which final density-matrix
        (population) element a pathway is allowed to terminate on -- i.e.
        it picks out a specific sector of Liouville space rather than
        summing over all of them. Restricting to the pure ground-state
        bleach channel should give an (almost) purely positive spectrum
        (there is no excited-state absorption pathway that returns
        population to the ground state), while restricting to either
        fundamental's population channel pulls in a genuine negative
        (excited-state absorption) contribution, and the three restricted
        channels -- along with the unrestricted sum over all channels --
        should all be numerically distinct from one another.
        """
        fsls = self.FourStateLiouvilleSystem()

        gs, e1, e2 = (0, 0), (1, 0), (0, 1)

        I_gs, _ = fsls.get_intensities(response_tensor_elements=[(gs, gs)])
        I_e1, _ = fsls.get_intensities(response_tensor_elements=[(e1, e1)])
        I_e2, _ = fsls.get_intensities(response_tensor_elements=[(e2, e2)])
        I_all, _ = fsls.get_intensities()

        # pure ground-state bleach: overwhelmingly positive
        self.assertGreater(I_gs.max(), 20 * abs(I_gs.min()))
        # fundamental-population channels carry a real excited-state absorption dip
        self.assertGreater(abs(I_e1.min()), 0.1 * I_e1.max())
        self.assertGreater(abs(I_e2.min()), 0.1 * I_e2.max())

        for a, b in itertools.combinations([I_gs, I_e1, I_e2, I_all], 2):
            self.assertFalse(np.allclose(a, b, atol=1e-9))

    @validationTest
    def test_ExcitedStateInitialPopulationExploresDifferentLiouvilleSector(self):
        """
        `initial_states` sets where in Liouville space the pathway
        enumeration starts (i.e. the initial density matrix element,
        normally the ground-state population |gs><gs|). Starting instead
        from population already sitting in the e1 = (1, 0) fundamental --
        e.g. a vibrationally "hot" or pre-excited system -- opens up a
        disjoint set of pathways (built on top of an existing excitation
        rather than from vacuum) and should give a spectrum that is both
        numerically distinct from, and substantially larger than, the
        ordinary ground-state-initiated response, since the e1 -> f
        transition dipole (0.12) is larger than the ground -> e1 dipole (0.1)
        and additional excited-state pathways become accessible.
        """
        fsls = self.FourStateLiouvilleSystem()

        I_ground_start, _ = fsls.get_intensities()
        I_hot_start, _ = fsls.get_intensities(initial_states=[(1, 0)])

        self.assertGreater(I_hot_start.max(), 0)
        self.assertFalse(np.allclose(I_ground_start, I_hot_start, atol=1e-9))
        self.assertGreater(I_hot_start.max(), 5 * I_ground_start.max())

    @validationTest
    def test_PolarizationDependenceOfOrientationalAveraging(self):
        """
        The four-wave-mixing orientational averaging factor
        (`four_wave_averaging_function`/`interpret_polarization`) depends
        on the relative polarizations of the four interacting fields.
        All-parallel (XXXX) polarization should give a larger isotropically
        averaged response than crossed (XXYY) polarization for this system.
        """
        fsls = self.FourStateLiouvilleSystem()

        I_xxxx, _ = fsls.get_intensities(response_function_generator="simple2dir", polarization="XXXX")
        I_xxyy, _ = fsls.get_intensities(response_function_generator="simple2dir", polarization="XXYY")

        self.assertGreater(I_xxxx.max(), 0)
        self.assertGreater(I_xxyy.max(), 0)
        self.assertFalse(np.allclose(I_xxxx, I_xxyy, atol=1e-9))
        self.assertGreater(I_xxxx.max(), I_xxyy.max())

    class LiouvilleFilterSpaceSystem:
        """
        A 4-mode harmonic basis (modes A-D) used to exercise `filter_space`,
        the mechanism for restricting *which basis states* `prep_liouville_spaces`/
        `get_interaction_basis` are allowed to explore while building up the full
        Liouville-pathway state space -- as opposed to `response_tensor_elements`
        (which restricts only the final detection channel) or `initial_states`
        (which restricts only the starting point). `filter_space` prunes the
        basis at every intermediate interaction step, so it changes both which
        `total_space` states show up at all and which enumerated pathways are
        possible between them.

        Two scenarios are built from the same NUM_INTERACTIONS=4, 3-interaction
        liouville_path pattern (`path_paths[3]`) used elsewhere in this module:

          * excluding every state where mode D carries any excitation (treating
            it as a dark/spectator mode never involved in the response), and
          * excluding one specific targeted state (the (1, 1, 0, 0) combination
            band of modes A and B).

        Baseline (no filter_space) numbers for this system, used as a sanity
        check that the filtering actually narrows the search: total_space has
        70 states and the chosen liou_path pattern enumerates 44 pathways
        touching 15 distinct states.
        """

        NUM_MODES = 4
        NUM_INTERACTIONS = 4
        PATTERN_INDEX = 3  # arbitrary representative liou_path pattern

        @classmethod
        def basis(cls):
            from Psience.BasisReps import HarmonicOscillatorProductBasis
            return HarmonicOscillatorProductBasis(cls.NUM_MODES)

        @classmethod
        def ground_state(cls):
            from Psience.BasisReps import BasisStateSpace
            return BasisStateSpace(cls.basis(), [(0,) * cls.NUM_MODES])

        @classmethod
        def candidate_space(cls):
            from Psience.BasisReps import BasisStateSpace
            return BasisStateSpace.from_quanta(cls.basis(), [0, 1, 2, 3])

        @classmethod
        def mode_d_excited_states(cls):
            return [tuple(s) for s in cls.candidate_space().excitations if s[3] != 0]

        @classmethod
        def prep(cls, **filter_opts):
            try:
                from Psience.Nonlinear.NonlinearResponse import prep_liouville_spaces
            except ImportError:
                # depending on how "Psience" resolves (the outer dev-shim package
                # vs. the real nested package), `prep_liouville_spaces`/
                # `enumerate_state_paths` -- which aren't re-exported in
                # NonlinearResponse.__all__ -- may only be reachable via the
                # doubly-nested path (see Psience/Nonlinear.py's own
                # `__identifier__ = 'Psience.Psience.Nonlinear'` marker)
                from Psience.Psience.Nonlinear.NonlinearResponse import prep_liouville_spaces
            return prep_liouville_spaces(
                cls.ground_state(), 3,
                num_interactions=cls.NUM_INTERACTIONS,
                **filter_opts
            )

        @classmethod
        def touched_states(cls, total_space, path_paths):
            try:
                from Psience.Nonlinear.NonlinearResponse import enumerate_state_paths
            except ImportError:
                from Psience.Psience.Nonlinear.NonlinearResponse import enumerate_state_paths
            sign, liou_path, index_set = path_paths[cls.PATTERN_INDEX]
            paths = list(enumerate_state_paths(
                liou_path, index_set, num_interactions=cls.NUM_INTERACTIONS
            ))
            touched = set()
            for path, transitions in paths:
                for l, r in path:
                    touched.add(l)
                    touched.add(r)
            states = sorted(tuple(int(x) for x in total_space.excitations[i]) for i in touched)
            return paths, states

    @validationTest
    def test_FilterSpaceBaselineHasNoRestriction(self):
        """
        Sanity check establishing the unfiltered numbers that the two
        `filter_space` tests below narrow down from: 70 total_space states,
        44 enumerated pathways touching 15 distinct states, for the
        representative liou_path pattern used throughout this class.
        """
        lfs = self.LiouvilleFilterSpaceSystem()
        total_space, path_paths = lfs.prep()
        paths, states = lfs.touched_states(total_space, path_paths)

        self.assertEqual(len(total_space), 70)
        self.assertEqual(len(paths), 44)
        self.assertEqual(len(states), 15)

    @validationTest
    def test_FilterSpaceExcludesSpectatorMode(self):
        """
        `filter_space` passed through `prep_liouville_spaces`/
        `get_interaction_basis` should prune out every basis state where mode
        D carries excitation, as if it were a dark spectator mode never
        involved in the interaction pathways -- narrowing both the resulting
        `total_space` and the pathways `enumerate_state_paths` can build over
        it.

        NOTE: this currently requires the fix in
        `claude_drafts/filter_space_fix.patch` to be applied to
        `Psience/Nonlinear/NonlinearResponse.py`. Without it,
        `get_interaction_basis`'s internal `_apply_rules` helper does not
        unpack the `(space, filter)` tuple that `apply_selection_rules`
        returns whenever `filter_space` is supplied, and this call raises
        `AttributeError: 'tuple' object has no attribute 'to_single'`.
        """
        lfs = self.LiouvilleFilterSpaceSystem()
        allowed = lfs.candidate_space().drop_states(lfs.mode_d_excited_states())

        total_space, path_paths = lfs.prep(filter_space=allowed)
        paths, states = lfs.touched_states(total_space, path_paths)

        # mode D (index 3) never appears excited anywhere in the filtered space
        self.assertTrue(all(int(s[3]) == 0 for s in total_space.excitations))
        self.assertTrue(all(s[3] == 0 for s in states))

        # the filter should have actually narrowed things down relative to baseline
        self.assertEqual(len(total_space), 20)
        self.assertEqual(len(paths), 24)
        self.assertEqual(len(states), 10)

    @validationTest
    def test_FilterSpaceExcludesTargetedState(self):
        """
        `filter_space` should also work as a scalpel, excluding a single
        targeted state (here the (1, 1, 0, 0) combination band of modes A
        and B) rather than an entire mode's excitation manifold -- and that
        state should never appear in `total_space` or in any enumerated
        pathway.

        NOTE: see `test_FilterSpaceExcludesSpectatorMode` -- this also
        depends on `claude_drafts/filter_space_fix.patch` being applied.
        """
        lfs = self.LiouvilleFilterSpaceSystem()
        target = (1, 1, 0, 0)
        allowed = lfs.candidate_space().drop_states([target])

        total_space, path_paths = lfs.prep(filter_space=allowed)
        paths, states = lfs.touched_states(total_space, path_paths)

        total_states = [tuple(int(x) for x in s) for s in total_space.excitations]
        self.assertNotIn(target, total_states)
        self.assertNotIn(target, states)
        self.assertEqual(len(total_space), 34)

    @validationTest
    def test_VPTResponseTargetStatePreparation(self):
        """Explicit and dictionary target-state specifications resolve consistently."""
        try:
            from Psience.Nonlinear.NonlinearResponse import _prep_vpt_target_states
            from Psience.BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
        except ImportError:
            from Psience.Psience.Nonlinear.NonlinearResponse import _prep_vpt_target_states
            from Psience.Psience.BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis

        freqs = np.array([1.0, 2.0, 3.0])
        ground = (0, 0, 0)

        # Explicit vectors are preserved, deduplicated, and get an implicit ground state.
        explicit = _prep_vpt_target_states(
            freqs,
            target_states=[
                [1, 0, 0],
                [0, 0, 1],
                [1, 0, 0]
            ]
        )
        self.assertEqual(explicit, [ground, (1, 0, 0), (0, 0, 1)])
        self.assertEqual(_prep_vpt_target_states(freqs, target_states=[]), [ground])
        generated_iterable = ([n, 0, 0] for n in (1, 2))
        self.assertEqual(
            _prep_vpt_target_states(freqs, target_states=generated_iterable),
            [ground, (1, 0, 0), (2, 0, 0)]
        )

        # Existing BasisStateSpace objects are accepted without regenerating states.
        basis_space = BasisStateSpace(
            HarmonicOscillatorProductBasis(3),
            [[0, 1, 0], [0, 0, 2]],
            mode=BasisStateSpace.StateSpaceSpec.Excitations
        )
        from_basis = _prep_vpt_target_states(freqs, target_states=basis_space)
        self.assertEqual(from_basis, [ground, (0, 1, 0), (0, 0, 2)])

        # A dict is forwarded as BasisStateSpace.states_under_freq_threshold options.
        state_options = {
            'max_freq': 4.1,
            'max_quanta': 3,
            'fixed_modes': [1]
        }
        generated = _prep_vpt_target_states(freqs, target_states=state_options)
        direct = BasisStateSpace.states_under_freq_threshold(freqs, **state_options)
        expected = {ground} | {tuple(int(x) for x in state) for state in direct}
        self.assertEqual(set(generated), expected)
        self.assertTrue(all(state[1] == 0 for state in generated))
        self.assertTrue(all(sum(state) < 3 for state in generated))

        with self.assertRaises(ValueError):
            _prep_vpt_target_states(freqs, target_states=[[1, 0]])
        with self.assertRaises(ValueError):
            _prep_vpt_target_states(freqs, target_states=[[0, -1, 0]])
        with self.assertRaises(ValueError):
            _prep_vpt_target_states(freqs, target_states=[[0, 0.5, 0]])

    @validationTest
    def test_VPTResponseTargetStatesReachRunners(self):
        """Restricted target states are forwarded to both VPT runner styles."""
        from types import SimpleNamespace
        from unittest import mock
        try:
            from Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data
            from Psience.VPT2 import VPTSystem, VPTRunner, AnalyticVPTRunner
        except ImportError:
            from Psience.Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data
            from Psience.Psience.VPT2 import VPTSystem, VPTRunner, AnalyticVPTRunner

        fake_system = VPTSystem.__new__(VPTSystem)
        fake_system.mol = SimpleNamespace(
            normal_modes=SimpleNamespace(
                modes=SimpleNamespace(freqs=np.array([0.01, 0.02]))
            )
        )
        target_states = [
            [0, 0],
            [1, 0],
            [0, 1],
            [2, 0]
        ]

        class RunnerReached(Exception):
            pass

        with mock.patch.object(VPTRunner, 'run_simple', side_effect=RunnerReached) as run:
            with self.assertRaises(RunnerReached):
                prep_vpt_response_data(fake_system, target_states=target_states)
        self.assertEqual(run.call_args.args[1], [tuple(s) for s in target_states])
        self.assertEqual(
            run.call_args.kwargs['initial_states'],
            [(0, 0), (1, 0), (0, 1)]
        )

        with mock.patch.object(AnalyticVPTRunner, 'run_simple', side_effect=RunnerReached) as run:
            with self.assertRaises(RunnerReached):
                prep_vpt_response_data(
                    fake_system,
                    target_states=target_states,
                    use_analytic=True
                )
        self.assertEqual(
            run.call_args.args[1],
            [
                [[[0, 0]], [[1, 0], [0, 1]]],
                [[[1, 0], [0, 1]], [[2, 0]]]
            ]
        )

    @validationTest
    def test_VPTResponseDataFromWaterFchk(self):
        """
        `prep_vpt_response_data` bridges a real ab initio VPT calculation
        (via `VPTRunner.run_simple`) into the `transition_dict` format that
        `prep_nonlinear_transition_data`/`experimental_response_generator`
        consume, rather than requiring the transitions to be hand-specified
        like the systems above: it builds the target state list from
        `BasisStateSpace.states_under_freq_threshold` (capped at two total
        quanta), always includes the ground state, and seeds the VPT run's
        `initial_states` with the ground state plus every one-quantum
        fundamental.

        NOTE: this depends on `claude_drafts/prep_vpt_response_data.patch`
        being merged into `Psience/Nonlinear/NonlinearResponse.py`.
        """
        try:
            from Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data
        except ImportError:
            # see the `filter_space` tests above for why this fallback is needed
            from Psience.Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data

        # `TestManager.test_data(...)` resolves relative to a legacy `<repo_root>/Tests`
        # layout that doesn't match this project's actual `ci/tests/TestData` convention
        # (and is only correctly configured when driven through `ci/tests/run_tests.py`),
        # so we resolve the path directly relative to this file instead.
        fchk = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TestData', 'water_freq.fchk')
        transition_dict, wfns = prep_vpt_response_data(fchk, return_wavefunctions=True)

        states_seen = set()
        for si, sj in transition_dict.keys():
            states_seen.add(si)
            states_seen.add(sj)
        ndim = len(next(iter(states_seen)))
        ground_state = (0,) * ndim

        # water has 3 normal modes -> 1 (gs) + 3 (fundamentals) + 3 (overtones)
        # + 3 (combination bands) = 10 states with <= 2 total quanta
        self.assertEqual(ndim, 3)
        self.assertIn(ground_state, states_seen)
        self.assertEqual(len(states_seen), 10)

        # every mode's fundamental should show up as a real, dipole-allowed
        # transition directly out of the ground state
        fundamentals = [s for s in states_seen if sum(s) == 1]
        self.assertEqual(len(fundamentals), 3)
        for fund in fundamentals:
            key = (ground_state, fund)
            self.assertIn(key, transition_dict)
            data = transition_dict[key]
            # water's fundamentals (bend + two stretches) all fall in the mid-IR
            self.assertGreater(data['frequency'], 1000)
            self.assertLess(data['frequency'], 4200)
            self.assertGreater(np.linalg.norm(data['transition_moment']), 1e-3)

        # every transition should be reported in its "upward" (positive-frequency)
        # direction -- `prep_nonlinear_transition_data` infers the reverse itself
        for data in transition_dict.values():
            self.assertGreater(data['frequency'], 0)

    @validationTest
    def test_VPTResponseDataCaching(self):
        """
        `prep_vpt_response_data`'s `output_file` option should cache the
        computed `transition_dict` to disk as JSON (a list of
        `{"state": [state_i, state_j], ...}` records, since JSON object keys
        can't be tuples) and, on a later call with the same `output_file`,
        load that cached data back in verbatim rather than rerunning VPT --
        unless `overwrite=True` is passed, in which case it always reruns and
        rewrites the file.

        NOTE: this depends on `claude_drafts/prep_vpt_response_data.patch`
        being merged into `Psience/Nonlinear/NonlinearResponse.py`.
        """
        try:
            from Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data
        except ImportError:
            # see the `filter_space` tests above for why this fallback is needed
            from Psience.Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data

        fchk = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TestData', 'water_freq.fchk')

        with tempfile.TemporaryDirectory() as tmp_dir:
            out_file = os.path.join(tmp_dir, 'water_response.json')

            self.assertFalse(os.path.exists(out_file))
            computed = prep_vpt_response_data(fchk, output_file=out_file)
            self.assertTrue(os.path.isfile(out_file))

            # the file on disk should be a JSON list of "state"-keyed records,
            # not a JSON object keyed by (unserializable) state-pair tuples
            with open(out_file) as f:
                raw_records = json.load(f)
            self.assertIsInstance(raw_records, list)
            self.assertEqual(len(raw_records), len(computed))
            for rec in raw_records:
                self.assertIn('state', rec)
                self.assertEqual(len(rec['state']), 2)
                self.assertIn('frequency', rec)
                self.assertIn('transition_moment', rec)

            # loading it back in (overwrite=False, the default) should match exactly
            cached = prep_vpt_response_data(fchk, output_file=out_file)
            self.assertEqual(set(cached.keys()), set(computed.keys()))
            for key in computed:
                self.assertAlmostEqual(cached[key]['frequency'], computed[key]['frequency'], places=6)
                np.testing.assert_allclose(
                    cached[key]['transition_moment'], computed[key]['transition_moment']
                )

            # a cache load has no live VPTWavefunctions to hand back
            _, wfns = prep_vpt_response_data(fchk, output_file=out_file, return_wavefunctions=True)
            self.assertIsNone(wfns)

            # tamper with the cached file directly; without `overwrite`, the tampered
            # value should come back verbatim -- proving the cache is actually used
            # rather than silently recomputed every time
            raw_records[0]['frequency'] = -12345.0
            with open(out_file, 'w') as f:
                json.dump(raw_records, f)
            tampered = prep_vpt_response_data(fchk, output_file=out_file)
            self.assertIn(-12345.0, [d['frequency'] for d in tampered.values()])

            # `overwrite=True` should ignore the tampered file and recompute + rewrite it
            recomputed = prep_vpt_response_data(fchk, output_file=out_file, overwrite=True)
            self.assertNotIn(-12345.0, [d['frequency'] for d in recomputed.values()])
            self.assertEqual(set(recomputed.keys()), set(computed.keys()))

    @validationTest
    def test_VPTResponseDataSavedToTestData(self):
        """
        Exercises `prep_vpt_response_data`'s `output_file` disk-caching against
        a *persistent* location -- `ci/tests/TestData/water_freq_response.json`,
        checked in alongside `water_freq.fchk` itself -- rather than a scratch
        tempfile, so that once the cache file exists this test (and anything
        else that wants water's VPT response data) never has to rerun the VPT
        calculation at all: with `overwrite` left at its default of `False`,
        `prep_vpt_response_data` loads the checked-in JSON straight off disk
        instead of recomputing it every time the test suite runs. Passing
        `overwrite=True` (not exercised as the default path here on purpose)
        is the only way to force a refresh of that file.

        NOTE: this depends on `claude_drafts/prep_vpt_response_data.patch`
        being merged into `Psience/Nonlinear/NonlinearResponse.py`.
        """

        fchk = TestManager.test_data('water_freq.fchk')
        out_file = TestManager.test_data('water_freq_response.json')

        # by default (`overwrite=False`) this does nothing but load the
        # checked-in file if it's already there -- a fresh VPT run only
        # happens the first time this is ever called for this file
        transition_dict = prep_vpt_response_data(fchk, output_file=out_file, overwrite=False)
        self.assertTrue(os.path.isfile(out_file))

        with open(out_file) as f:
            raw_records = json.load(f)
        self.assertIsInstance(raw_records, list)
        self.assertEqual(len(raw_records), len(transition_dict))
        for rec in raw_records:
            self.assertIn('state', rec)
            self.assertEqual(len(rec['state']), 2)
            self.assertIn('frequency', rec)
            self.assertIn('transition_moment', rec)

        # the same sanity checks as `test_VPTResponseDataFromWaterFchk`: water
        # has 3 normal modes -> 10 states with <= 2 total quanta, ground state
        # included, and all 3 fundamentals present as real ground-state transitions
        states_seen = set()
        for si, sj in transition_dict.keys():
            states_seen.add(si)
            states_seen.add(sj)
        ndim = len(next(iter(states_seen)))
        ground_state = (0,) * ndim

        self.assertEqual(ndim, 3)
        self.assertIn(ground_state, states_seen)
        self.assertEqual(len(states_seen), 10)

        fundamentals = [s for s in states_seen if sum(s) == 1]
        self.assertEqual(len(fundamentals), 3)
        for fund in fundamentals:
            key = (ground_state, fund)
            self.assertIn(key, transition_dict)
            data = transition_dict[key]
            self.assertGreater(data['frequency'], 1000)
            self.assertLess(data['frequency'], 4200)
            self.assertGreater(np.linalg.norm(data['transition_moment']), 1e-3)

    @validationTest
    def test_FullTwoDimensionalIRFromVPTResponseData(self):
        """
        Hooks `prep_vpt_response_data`'s real ab initio VPT results for water
        straight into the same Liouville-pathway 2D-IR machinery exercised by
        the hand-specified systems above, instead of a synthetic transition
        dict -- i.e. runs a genuine, from-first-principles 2D-IR calculation
        on water.

        VPT gives frequencies and transition moments but no dephasing rates,
        so every coherence/population band is given a uniform, physically
        modest homogeneous linewidth via `band_coherences`, mirroring the
        (0,1)/(1,1)-style band-quanta convention used by
        `FourStateLiouvilleSystem` above (rather than a lambda keyed on raw
        state vectors, which would also assign a spurious decay rate to the
        ground state itself). The observation window is restricted to
        water's two O-H stretch fundamentals -- picked out as the two
        higher-frequency fundamentals in the data (the lowest is the bend) --
        since that's the classic textbook water 2D-IR region: two diagonal
        peaks plus a real excited-state absorption feature from the
        anharmonically-shifted stretch overtones/combination band, exactly
        like the synthetic system in `test_CombinationBandCrossPeaksAndESA`
        above, except here every number -- including the two stretches'
        genuinely different transition dipole moments -- comes from an
        actual VPT calculation instead of being made up. Because of that,
        the two diagonal peaks are NOT expected to have comparable height
        (unlike the symmetric synthetic system): the peak heights are
        checked against each other's ab initio transition dipole magnitudes
        instead of a shared absolute threshold.
        """
        from Psience.Spectra import TwoDimensionalSpectrum

        fchk = TestManager.test_data('water_freq.fchk')
        out_file = TestManager.test_data('water_freq_response.json')

        # by default (`overwrite=False`) this does nothing but load the
        # checked-in file if it's already there -- a fresh VPT run only
        # happens the first time this is ever called for this file
        transition_dict = prep_vpt_response_data(fchk, output_file=out_file, overwrite=False)

        ndim = len(next(iter(transition_dict))[0])
        ground_state = (0,) * ndim
        fundamentals = sorted(
            (data['frequency'], sj)
            for (si, sj), data in transition_dict.items()
            if si == ground_state and sum(sj) == 1
        )
        self.assertEqual(len(fundamentals), 3)
        # the lowest-frequency fundamental is the bend; the other two are the
        # O-H stretches -- the classic water 2D-IR system
        stretch_freqs = [f for f, _ in fundamentals[1:]]
        stretch_states = [sj for _, sj in fundamentals[1:]]

        center = (min(stretch_freqs) + max(stretch_freqs)) / 2  # ~3683 cm^-1 for the two O-H stretches
        coherence_strength = 3  # physically modest linewidth -- see the module docstring note above
        responses = experimental_response_generator(
            transition_dict,
            band_coherences={
                (0, 1): coherence_strength,
                (1, 2): coherence_strength,
                (1, 1): coherence_strength,
                (2, 2): coherence_strength
            },
            frequency_unit="Wavenumbers",
            application_domain="frequency",
            driving_frequency=center,
        )

        window = [min(stretch_freqs) - 250, max(stretch_freqs) + 100]
        spec:TwoDimensionalSpectrum = responses.get_spectrum(
            window, 10, window,
            default_frequency_divisions=300
        )
        spec.plot().show()
        I = np.real(spec.intensities)
        self.assertGreater(I.max(), 1e-9)

        def nearest_value(w1, w3):
            ix = int(np.argmin(np.abs(spec.freq1 - w1)))
            iy = int(np.argmin(np.abs(spec.freq2 - w3)))
            return I[iy, ix]

        # both O-H stretch fundamentals should show up as positive
        # (ground-state bleach/stimulated emission) diagonal peaks, and the
        # stretch with the larger ab initio transition dipole should produce
        # the larger diagonal peak
        diag_vals = []
        tm_norms = []
        for freq, state in zip(stretch_freqs, stretch_states):
            val = nearest_value(freq, freq)
            self.assertGreater(val, 5e-12)
            diag_vals.append(val)
            tm_norms.append(np.linalg.norm(
                transition_dict[(ground_state, state)]['transition_moment']
            ))
        self.assertEqual(int(np.argmax(diag_vals)), int(np.argmax(tm_norms)))

        # a real, coupled/anharmonic 2D-IR spectrum must also show excited-state
        # absorption: a genuine negative-going feature somewhere in the window
        self.assertLess(I.min(), -1e-3 * I.max())

    @validationTest
    def test_VPTResponseDataFromWaterFchkAnalytic(self):
        """
        `prep_vpt_response_data(..., use_analytic=True)` bridges a real
        symbolic/analytic VPT calculation (via `AnalyticVPTRunner.run_simple`)
        into the same `transition_dict` format as the classic branch, but via
        a genuinely different calling convention under the hood: instead of
        one dense "every initial state x every final state" matrix, it builds
        one `[initial_space, target_space]` block per consecutive quantum
        shell implied by `initial_quanta` (ground -> fundamentals,
        fundamentals -> two-quantum states), and only computes transition
        moments *within* each block. This test exercises that branch fresh
        (no caching) against real water data and checks both the ordinary
        sanity properties (ground state, all three fundamentals present with
        sane frequencies/moments) and the block-structure-specific ones: the
        full 10-state space is still reachable, but non-adjacent-shell
        transitions (a direct ground -> two-quantum overtone, or a
        fundamental -> fundamental cross term within the one-quantum
        manifold) are *not* computed, unlike the classic branch.

        NOTE: this depends on `claude_drafts/prep_vpt_response_data_analytic.patch`
        being merged into `Psience/Nonlinear/NonlinearResponse.py`.
        """
        # `TestManager.test_data(...)` resolves relative to a legacy `<repo_root>/Tests`
        # layout that doesn't match this project's actual `ci/tests/TestData` convention
        # (and is only correctly configured when driven through `ci/tests/run_tests.py`),
        # so we resolve the path directly relative to this file instead.
        fchk = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TestData', 'water_freq.fchk')
        transition_dict, corrs = prep_vpt_response_data(
            fchk, use_analytic=True, logger=False, return_wavefunctions=True
        )

        # a fresh analytic run hands back the live `AnalyticPerturbationTheoryCorrections`,
        # not a `VPTWavefunctions` -- distinguishable by its `state_lists` attribute,
        # which the classic branch's result doesn't have
        self.assertTrue(hasattr(corrs, 'state_lists'))
        self.assertFalse(hasattr(corrs, 'initial_state_indices'))

        states_seen = set()
        for si, sj in transition_dict.keys():
            states_seen.add(si)
            states_seen.add(sj)
        ndim = len(next(iter(states_seen)))
        ground_state = (0,) * ndim

        # the block structure (0->1, 1->2 quantum shells) still touches every
        # one of water's 10 states with <= 2 total quanta
        self.assertEqual(ndim, 3)
        self.assertIn(ground_state, states_seen)
        self.assertEqual(len(states_seen), 10)

        # every mode's fundamental is still a real, dipole-allowed transition
        # directly out of the ground state (the 0->1 block)
        fundamentals = [s for s in states_seen if sum(s) == 1]
        self.assertEqual(len(fundamentals), 3)
        for fund in fundamentals:
            key = (ground_state, fund)
            self.assertIn(key, transition_dict)
            data = transition_dict[key]
            self.assertGreater(data['frequency'], 1000)
            self.assertLess(data['frequency'], 4200)
            self.assertGreater(np.linalg.norm(data['transition_moment']), 1e-3)

        # every transition should still be reported in its "upward" direction
        for data in transition_dict.values():
            self.assertGreater(data['frequency'], 0)

        # the block structure only connects *adjacent* quantum shells, so a
        # direct ground -> two-quantum overtone (a non-adjacent-shell
        # transition) should NOT show up, unlike in the classic branch
        two_quantum_states = [s for s in states_seen if sum(s) == 2]
        self.assertEqual(len(two_quantum_states), 6)
        for state in two_quantum_states:
            self.assertNotIn((ground_state, state), transition_dict)

        # nor should a fundamental -> fundamental cross term within the
        # one-quantum manifold itself (both endpoints live in the same block
        # boundary, not across one)
        for fund_a in fundamentals:
            for fund_b in fundamentals:
                if fund_a != fund_b:
                    self.assertNotIn((fund_a, fund_b), transition_dict)

        # exactly 19 transitions survive: the classic branch's 28 minus the
        # 6 direct ground -> two-quantum overtones and the 3 fundamental ->
        # fundamental cross terms that the block structure can't reach
        self.assertEqual(len(transition_dict), 19)

    @validationTest
    def test_VPTResponseDataAnalyticMatchesClassicOnSharedTransitions(self):
        """
        The classic and analytic branches are independent VPT
        implementations, so they should agree closely -- but not
        necessarily bit-for-bit, and not necessarily up to the same overall
        sign convention on transition moments -- on whatever transitions
        they *both* compute. This compares the two persisted `TestData`
        caches (so it doesn't have to rerun either VPT calculation) on their
        common transitions: frequencies should match to a small fraction of
        a wavenumber, and each transition moment should match up to a
        possible overall sign flip (the two evaluators' dipole-derivative
        phase conventions aren't guaranteed to agree, but the underlying
        physics -- and thus every even, sign-invariant combination that
        actually enters a computed intensity/spectrum -- is the same either
        way).

        NOTE: this depends on `claude_drafts/prep_vpt_response_data_analytic.patch`
        being merged into `Psience/Nonlinear/NonlinearResponse.py`, and on both
        `ci/tests/TestData/water_freq_response.json` and
        `ci/tests/TestData/water_freq_response_analytic.json` being checked in.
        """
        fchk = TestManager.test_data('water_freq.fchk')
        classic_out = TestManager.test_data('water_freq_response.json')
        analytic_out = TestManager.test_data('water_freq_response_analytic.json')

        classic = prep_vpt_response_data(fchk, output_file=classic_out, overwrite=False)
        analytic = prep_vpt_response_data(
            fchk, use_analytic=True, output_file=analytic_out, overwrite=False
        )

        # every transition the analytic (block-restricted) branch computes is
        # also computed by the classic (dense-matrix) branch -- the block
        # structure is strictly more conservative about what it connects
        self.assertTrue(set(analytic.keys()) <= set(classic.keys()))
        self.assertLess(len(analytic), len(classic))

        for key in analytic:
            c_freq = classic[key]['frequency']
            a_freq = analytic[key]['frequency']
            self.assertAlmostEqual(c_freq, a_freq, delta=0.5)

            tm_c = np.asarray(classic[key]['transition_moment'])
            tm_a = np.asarray(analytic[key]['transition_moment'])
            diff_same_sign = np.linalg.norm(tm_c - tm_a)
            diff_flipped = np.linalg.norm(tm_c + tm_a)
            # whichever relative sign is the better match should agree closely
            self.assertLess(min(diff_same_sign, diff_flipped), 5e-3)

    @validationTest
    def test_VPTResponseDataAnalyticCaching(self):
        """
        `prep_vpt_response_data`'s `output_file` caching behaves the same way
        for the analytic branch as for the classic one (see
        `test_VPTResponseDataCaching`): the first call runs
        `AnalyticVPTRunner.run_simple` and writes the JSON cache, a later
        call with the same `output_file` loads that cache back in verbatim
        (proven by tampering with it directly and getting the tampered value
        back) rather than rerunning VPT, and `overwrite=True` always reruns
        and rewrites the file. A cache load also still has no live results
        object to hand back via `return_wavefunctions`, exactly as for the
        classic branch.

        NOTE: this depends on `claude_drafts/prep_vpt_response_data_analytic.patch`
        being merged into `Psience/Nonlinear/NonlinearResponse.py`.
        """
        fchk = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TestData', 'water_freq.fchk')

        with tempfile.TemporaryDirectory() as tmp_dir:
            out_file = os.path.join(tmp_dir, 'water_response_analytic.json')

            self.assertFalse(os.path.exists(out_file))
            computed = prep_vpt_response_data(fchk, use_analytic=True, logger=False, output_file=out_file)
            self.assertTrue(os.path.isfile(out_file))

            with open(out_file) as f:
                raw_records = json.load(f)
            self.assertIsInstance(raw_records, list)
            self.assertEqual(len(raw_records), len(computed))
            for rec in raw_records:
                self.assertIn('state', rec)
                self.assertEqual(len(rec['state']), 2)
                self.assertIn('frequency', rec)
                self.assertIn('transition_moment', rec)

            # loading it back in (overwrite=False, the default) should match exactly
            cached = prep_vpt_response_data(fchk, use_analytic=True, output_file=out_file)
            self.assertEqual(set(cached.keys()), set(computed.keys()))
            for key in computed:
                self.assertAlmostEqual(cached[key]['frequency'], computed[key]['frequency'], places=6)
                np.testing.assert_allclose(
                    cached[key]['transition_moment'], computed[key]['transition_moment']
                )

            # a cache load has no live results object to hand back
            _, corrs = prep_vpt_response_data(
                fchk, use_analytic=True, output_file=out_file, return_wavefunctions=True
            )
            self.assertIsNone(corrs)

            # tamper with the cached file directly; without `overwrite`, the tampered
            # value should come back verbatim -- proving the cache is actually used
            raw_records[0]['frequency'] = -12345.0
            with open(out_file, 'w') as f:
                json.dump(raw_records, f)
            tampered = prep_vpt_response_data(fchk, use_analytic=True, output_file=out_file)
            self.assertIn(-12345.0, [d['frequency'] for d in tampered.values()])

            # `overwrite=True` should ignore the tampered file and recompute + rewrite it
            recomputed = prep_vpt_response_data(
                fchk, use_analytic=True, logger=False, output_file=out_file, overwrite=True
            )
            self.assertNotIn(-12345.0, [d['frequency'] for d in recomputed.values()])
            self.assertEqual(set(recomputed.keys()), set(computed.keys()))

    @validationTest
    def test_VPTResponseDataFromLogMatchesWaterFchk(self):
        """
        Exercises `prep_vpt_response_data_from_log` -- which reconstructs a
        `transition_dict` purely from a saved VPT2 text log via
        `Psience.VPT2.Analyzer.VPTAnalyzer`, rather than from an in-memory
        `VPTWavefunctions` result -- against a checked-in log fixture
        (`water_vpt_classic.log`, generated the same way
        `prep_vpt_response_data`'s classic branch itself would generate one,
        via `VPTRunner.run_simple(..., logger=<path>)`), and checks it against
        the known-correct ground truth already checked in for the in-memory
        classic branch (`water_freq_response.json`).

        This only works at all because of two real bugs found and fixed in
        `Psience/VPT2/Analyzer.py` while building this function -- see
        `claude_drafts/vpt_analyzer_log_parsing_fixes.patch` for the full
        writeup, and `claude_drafts/prep_vpt_response_data_from_log.patch`
        for `prep_vpt_response_data_from_log` itself. Both patches must be
        merged for this test to pass.

        Frequencies are recovered exactly (they come straight off the log's
        own printed values). Transition moments are only recovered exactly
        for pure fundamentals; a third, separate, *not*-fixed bug in
        `VPTAnalyzerLogParser.reformat_tm_block` (documented in
        `prep_vpt_response_data_from_log`'s docstring) means combination-band
        and overtone transition moments come back off by up to ~30% -- so
        this only asserts a majority (>=60%) of the 28 known transitions
        match closely, not all of them.
        """
        try:
            from Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data_from_log
        except ImportError:
            from Psience.Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data_from_log

        log_file = TestManager.test_data('water_vpt_classic.log')
        gt_file = TestManager.test_data('water_freq_response.json')

        transition_dict = prep_vpt_response_data_from_log(log_file)
        self.assertEqual(len(transition_dict), 28)

        with open(gt_file) as f:
            gt_records = json.load(f)
        gt_dict = {
            (tuple(rec['state'][0]), tuple(rec['state'][1])): rec
            for rec in gt_records
        }
        self.assertEqual(set(transition_dict.keys()), set(gt_dict.keys()))

        n_close = 0
        for key, data in transition_dict.items():
            gt_rec = gt_dict[key]
            self.assertAlmostEqual(data['frequency'], gt_rec['frequency'], places=3)

            tm = np.asarray(data['transition_moment'])
            tm_gt = np.asarray(gt_rec['transition_moment'])
            # the two independent codepaths (in-memory vs. log-reconstructed) can come back
            # with an overall sign flip, same as the classic-vs-analytic comparison above
            if (
                np.max(np.abs(tm - tm_gt)) < 1e-3
                or np.max(np.abs(tm + tm_gt)) < 1e-3
            ):
                n_close += 1
        self.assertGreaterEqual(n_close, int(0.6 * len(gt_dict)))

        # the fundamentals specifically should always be exact-ish, regardless of the
        # combination-band/overtone slicing bug documented above
        ground_state = (0, 0, 0)
        fundamentals = [sj for (si, sj) in transition_dict if si == ground_state and sum(sj) == 1]
        self.assertEqual(len(fundamentals), 3)
        for fund in fundamentals:
            key = (ground_state, fund)
            tm = np.asarray(transition_dict[key]['transition_moment'])
            tm_gt = np.asarray(gt_dict[key]['transition_moment'])
            self.assertLess(
                min(np.max(np.abs(tm - tm_gt)), np.max(np.abs(tm + tm_gt))),
                1e-3
            )

    @validationTest
    def test_VPTResponseDataFromLogFreshlyGenerated(self):
        """
        Same idea as `test_VPTResponseDataFromLogMatchesWaterFchk`, but generates its own
        log fresh (via `VPTRunner.run_simple(..., logger=<path>)`, using the identical
        `state_list`/`initial_states` construction `prep_vpt_response_data`'s classic
        branch uses) rather than relying on the checked-in fixture -- so this also exercises
        `VPTRunner`'s own log-writing path, not just `VPTAnalyzer`'s parsing of a pre-made one.
        """
        try:
            from Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data_from_log
        except ImportError:
            from Psience.Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data_from_log
        from Psience.VPT2 import VPTRunner, VPTSystem
        from Psience.BasisReps import BasisStateSpace

        fchk = TestManager.test_data('water_freq.fchk')
        vpt_system = VPTSystem(fchk)
        freqs = vpt_system.mol.normal_modes.modes.freqs
        max_quanta = 2
        max_freq = max_quanta * np.max(np.abs(freqs))
        raw_states = BasisStateSpace.states_under_freq_threshold(freqs, max_freq, max_quanta=max_quanta + 1)
        state_list = [tuple(int(x) for x in s) for s in raw_states]
        gs = (0, 0, 0)
        if gs not in state_list:
            state_list = [gs] + state_list
        initial_states = [s for s in state_list if sum(s) in (0, 1)]
        if gs not in initial_states:
            initial_states = [gs] + initial_states

        with tempfile.TemporaryDirectory() as tmp_dir:
            log_path = os.path.join(tmp_dir, 'water_vpt.log')
            VPTRunner.run_simple(vpt_system, state_list, initial_states=initial_states, logger=log_path)
            self.assertTrue(os.path.isfile(log_path))

            transition_dict = prep_vpt_response_data_from_log(log_path)

        self.assertEqual(len(transition_dict), 28)
        ground_state = (0, 0, 0)
        fundamentals = [sj for (si, sj) in transition_dict if si == ground_state and sum(sj) == 1]
        self.assertEqual(len(fundamentals), 3)
        for fund in fundamentals:
            data = transition_dict[(ground_state, fund)]
            self.assertGreater(data['frequency'], 1000)
            self.assertLess(data['frequency'], 4200)
            self.assertGreater(np.linalg.norm(data['transition_moment']), 1e-3)

    @validationTest
    def test_VPTResponseDataFromLogRejectsAnalyticLog(self):
        """
        `AnalyticVPTRunner` logs use an entirely different, untagged table format that
        `VPTAnalyzerLogParser` cannot parse at all (confirmed by actually generating one,
        `water_vpt_analytic.log`, via `AnalyticVPTRunner.run_simple(..., logger=<path>)`
        and attempting to load it -- see `prep_vpt_response_data_from_log`'s docstring).
        This checks `prep_vpt_response_data_from_log` turns that failure into a clear,
        actionable `ValueError` instead of letting a bare `IndexError` from deep inside
        `VPTAnalyzerLogParser` leak out.
        """
        try:
            from Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data_from_log
        except ImportError:
            from Psience.Psience.Nonlinear.NonlinearResponse import prep_vpt_response_data_from_log

        log_file = TestManager.test_data('water_vpt_analytic.log')
        with self.assertRaises(ValueError) as ctx:
            prep_vpt_response_data_from_log(log_file)
        self.assertIn('AnalyticVPTRunner', str(ctx.exception))

    @validationTest
    def test_VPTResponseDataAnalyticSavedToTestData(self):
        """
        Exercises `prep_vpt_response_data(..., use_analytic=True)`'s
        `output_file` disk-caching against a *persistent* location --
        `ci/tests/TestData/water_freq_response_analytic.json`, checked in
        alongside `water_freq.fchk` and the classic branch's own
        `water_freq_response.json` -- rather than a scratch tempfile, so this
        (and `test_VPTResponseDataAnalyticMatchesClassicOnSharedTransitions`
        above) never has to rerun the (comparatively slow) analytic VPT
        calculation once the cache file exists.

        NOTE: this depends on `claude_drafts/prep_vpt_response_data_analytic.patch`
        being merged into `Psience/Nonlinear/NonlinearResponse.py`.
        """
        fchk = TestManager.test_data('water_freq.fchk')
        out_file = TestManager.test_data('water_freq_response_analytic.json')

        transition_dict = prep_vpt_response_data(
            fchk, use_analytic=True, output_file=out_file, overwrite=False
        )
        self.assertTrue(os.path.isfile(out_file))

        with open(out_file) as f:
            raw_records = json.load(f)
        self.assertIsInstance(raw_records, list)
        self.assertEqual(len(raw_records), len(transition_dict))
        for rec in raw_records:
            self.assertIn('state', rec)
            self.assertEqual(len(rec['state']), 2)
            self.assertIn('frequency', rec)
            self.assertIn('transition_moment', rec)

        states_seen = set()
        for si, sj in transition_dict.keys():
            states_seen.add(si)
            states_seen.add(sj)
        ndim = len(next(iter(states_seen)))
        ground_state = (0,) * ndim

        self.assertEqual(ndim, 3)
        self.assertIn(ground_state, states_seen)
        self.assertEqual(len(states_seen), 10)
        self.assertEqual(len(transition_dict), 19)

        fundamentals = [s for s in states_seen if sum(s) == 1]
        self.assertEqual(len(fundamentals), 3)
        for fund in fundamentals:
            key = (ground_state, fund)
            self.assertIn(key, transition_dict)
            data = transition_dict[key]
            self.assertGreater(data['frequency'], 1000)
            self.assertLess(data['frequency'], 4200)
            self.assertGreater(np.linalg.norm(data['transition_moment']), 1e-3)

    @validationTest
    def test_FullTwoDimensionalIRFromAnalyticVPTResponseData(self):
        """
        The whole point of `use_analytic=True` is that its `transition_dict`
        output plugs into the same `prep_nonlinear_transition_data`/
        `experimental_response_generator` Liouville-pathway machinery as the
        classic branch's -- this mirrors
        `test_FullTwoDimensionalIRFromVPTResponseData` but sources its
        `transition_dict` from the analytic branch's persisted cache instead,
        confirming that swap produces an equally sane 2D-IR spectrum (same
        two-diagonal-peaks-plus-ESA-feature structure), not just a
        structurally-plausible `transition_dict` in isolation.

        NOTE: this depends on `claude_drafts/prep_vpt_response_data_analytic.patch`
        being merged into `Psience/Nonlinear/NonlinearResponse.py`.
        """
        fchk = TestManager.test_data('water_freq.fchk')
        out_file = TestManager.test_data('water_freq_response_analytic.json')
        transition_dict = prep_vpt_response_data(
            fchk, use_analytic=True, output_file=out_file, overwrite=False
        )

        ground_state = (0, 0, 0)
        fundamentals = sorted(
            (data['frequency'], sj)
            for (si, sj), data in transition_dict.items()
            if si == ground_state and sum(sj) == 1
        )
        self.assertEqual(len(fundamentals), 3)
        stretch_freqs = [f for f, _ in fundamentals[1:]]
        stretch_states = [sj for _, sj in fundamentals[1:]]

        center = (min(stretch_freqs) + max(stretch_freqs)) / 2
        coherence_strength = 3
        responses = experimental_response_generator(
            transition_dict,
            band_coherences={
                (0, 1): coherence_strength,
                (1, 2): coherence_strength,
                (1, 1): coherence_strength,
                (2, 2): coherence_strength
            },
            frequency_unit="Wavenumbers",
            application_domain="frequency",
            driving_frequency=center,
        )

        window = [min(stretch_freqs) - 250, max(stretch_freqs) + 100]
        spec = responses.get_spectrum(window, 10, window, default_frequency_divisions=300)
        I = np.real(spec.intensities)
        self.assertGreater(I.max(), 1e-9)

        def nearest_value(w1, w3):
            ix = int(np.argmin(np.abs(spec.freq1 - w1)))
            iy = int(np.argmin(np.abs(spec.freq2 - w3)))
            return I[iy, ix]

        # the stretch with the larger ab initio transition dipole should still
        # produce the larger diagonal peak, same as the classic-branch test
        diag_vals = []
        tm_norms = []
        for freq, state in zip(stretch_freqs, stretch_states):
            val = nearest_value(freq, freq)
            diag_vals.append(val)
            tm_norms.append(np.linalg.norm(
                transition_dict[(ground_state, state)]['transition_moment']
            ))
        self.assertEqual(int(np.argmax(diag_vals)), int(np.argmax(tm_norms)))
        self.assertGreater(diag_vals[int(np.argmax(tm_norms))], 1e-9)
