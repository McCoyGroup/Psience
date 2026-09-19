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

    @debugTest
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
