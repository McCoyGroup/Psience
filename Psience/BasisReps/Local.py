import numpy as np
import scipy.linalg
from McUtils.Data import UnitsData
import McUtils.Devutils as dev
import McUtils.Numputils as nput
import McUtils.Coordinerds as coordops
from .StateSpaces import AbstractStateSpace, BasisStateSpace, BraKetSpace
from .Representations import ExpansionRepresentation
from .HarmonicOscillator import HarmonicOscillatorProductBasis
from ..Wavefun import MatrixWavefunctions
from . import Util as util

__all__ = [
    "LocalHarmonicModel",
    "CustomLocalHarmonicModel",
    "TaborCHModel"
]

class LocalHarmonicModel:
    def __init__(self, f, g,
                 internals=None,
                 anharmonic_scalings=None,
                 anharmonic_couplings=None,
                 anharmonic_shifts=None,
                 freqs=None,
                 modes=None,
                 dipole_derivatives=None,
                 **operator_settings):
        self.f = np.asanyarray(f)
        self.g = np.asanyarray(g)
        if freqs is None:
            freqs = scipy.linalg.eigvalsh(self.f, self.g, type=3)
            freqs = np.sign(freqs) * np.sqrt(np.abs(freqs))
        self.freqs = np.asanyarray(freqs)
        self.local_freqs = np.diag(f)
        self.basis = HarmonicOscillatorProductBasis(len(f))
        self.modes = modes
        self.rep = 1 / 2 * (
                self.basis.representation(
                    "p", "p",
                    coeffs=self.g,
                    name="T",
                    **operator_settings
                ) +
                self.basis.representation(
                    "x", "x",
                    coeffs=self.f,
                    axes=[[0, 1], [0, 1]],
                    name="V",
                    **operator_settings
                )
        )
        # self.dip_reps: list[ExpansionRepresentation] | None
        if dipole_derivatives is not None:
            dipole_derivatives = [
                self.basis.representation(
                    "x",
                    coeffs=d,
                    axes=[[0], [0]],
                    name="M",
                    **operator_settings
                )
                for d in dipole_derivatives
            ]
        self.dip_reps = dipole_derivatives
        self.internals = internals
        self.scalings = anharmonic_scalings
        self.couplings = anharmonic_couplings
        self.shifts = anharmonic_shifts

    def get_hamiltonian(self, states, coupled_space=None, remove_zpe=True, **rep_opts):
        states = self.prep_states(states)
        if coupled_space is None:
            if isinstance(states, BraKetSpace):
                coupled_space = states
                states = states.bras.take_unique()
            else:
                coupled_space = self.rep.get_transformed_space(states, filter_space=states)[0]
        ham = self.rep.get_representation_matrix(
            coupled_space,
            states,
            zero_element_warning=False,
            **rep_opts
        ).asarray()

        if remove_zpe:
            zpe = np.sum(np.abs(self.freqs)) * 1/2
            ham = ham - np.eye(len(ham)) * zpe

        if (
                self.scalings is not None
                or self.couplings is not None
                or self.shifts is not None
        ):
            if self.internals is not None:
                ham = util.modify_internal_hamiltonian(
                    ham,
                    self.internals,
                    states=states.excitations,
                    scaling_types=self.scalings,
                    coupling_types=self.couplings,
                    shift_types=self.shifts,
                )
            else:
                ham = util.modify_hamiltonian(
                    ham,
                    [tuple(s) for s in states.excitations],
                    scaling_types=self.scalings,
                    coupling_types=self.couplings,
                    shift_types=self.shifts
                )

        return ham

    def get_dipoles(self, states, coupled_space=None, **rep_opts):
        if self.dip_reps is None:
            raise ValueError("`dipole_derivatives` must be supplied to compute dipole matrices")
        states = self.prep_states(states)
        if coupled_space is None:
            if isinstance(states, BraKetSpace):
                coupled_space = states
                states = states.bras.take_unique()
            else:
                coupled_space = self.dip_reps[0].get_transformed_space(states, filter_space=states)[0]
        return np.moveaxis(
            np.array([
                r.get_representation_matrix(
                    coupled_space,
                    states,
                    zero_element_warning=False,
                    **rep_opts
                ).asarray()
                for r in self.dip_reps
            ]),
            0, -1
        )

    def prep_states(self, states):
        if not isinstance(states, AbstractStateSpace):
            if isinstance(states, dict):
                states = states.copy()
                thresh = states.pop("max_freq")
                states = BasisStateSpace.states_under_freq_threshold(
                    self.local_freqs,
                    thresh,
                    **states
                )
                if not isinstance(states, AbstractStateSpace):
                    states = BasisStateSpace(self.basis, states)
                    states = states.as_sorted()
            elif nput.is_int(states):
                states = BasisStateSpace.from_quanta(self.basis, states)
            else:
                if (
                        nput.is_numeric_array_like(states)
                    and np.asanyarray(states).ndim == 2
                ):
                    states = BasisStateSpace(self.basis, states)
                else:
                    s0 = self.prep_states(states[0])
                    for s1 in states[1:]:
                        s0 = s0.union(self.prep_states(s1))
                    states = s0

        return states
    def get_wavefunctions(self, states, coupled_space=None, remove_zpe=True, include_dipoles=True,
                          expansion_opts=None,
                          **wfn_opts
                          ):
        if expansion_opts is None:
            expansion_opts = {}

        states = self.prep_states(states)

        ham = self.get_hamiltonian(states, coupled_space=coupled_space, remove_zpe=remove_zpe, **expansion_opts)
        engs, coeffs = np.linalg.eigh(ham)
        ground_state = BasisStateSpace.from_quanta(self.basis, 0)
        gs_loc = states.find(ground_state, missing_val=-1)[0]
        if gs_loc == -1:
            if remove_zpe:
                zpe = 0
            else:
                zpe = np.sum(np.abs(self.freqs)) * 1/2
            engs = np.concatenate([[zpe], engs])
            coeffs = np.pad(coeffs, [[1, 0], [1, 0]])
            coeffs[0, 0] = 1
            states = ground_state.union(states)
        elif gs_loc > 0:
            engs = np.concatenate([engs[(gs_loc,)], engs[:gs_loc], engs[gs_loc+1:]])
            coeffs = np.concatenate([coeffs[:, (gs_loc,)], engs[:, :gs_loc], engs[:, gs_loc+1:]], axis=1)


        if include_dipoles:
            dips = self.get_dipoles(states, coupled_space=coupled_space, **expansion_opts)
        else:
            dips = None

        return LocalHarmonicWavefunctions(self,
                                          engs, coeffs,
                                          basis=states,
                                          dipole_matrix=dips,
                                          hamiltonian=ham,
                                          **wfn_opts
                                          )

    @classmethod
    def _match_mode_label(cls, test, freq, label):
        if label is None: return test is None
        if len(test) == 2 and (
                not (test[0] is None or isinstance(test[0], str))
        ):
            if nput.is_numeric(test[0]):
                test_freqs = test
                test = None
            else:
                test_freqs, test = test
        else:
            test_freqs = None

        if test_freqs is not None:
            lf, mF = test_freqs
            if freq < lf or freq > mF: return False

        if test is None: return True

        if isinstance(test, str):
            test = [test]
        n_test = len(test)
        if len(label) < n_test: return False
        return all(
            t is None or (t == l if isinstance(t, str) else l in t)
            for t,l in zip(test, label[-n_test:])
        )

    @classmethod
    def localize_internal_modes(cls,
                                nms, internals,
                                localization_mode_spaces=None,
                                mode_labels=None,
                                use_nonnegative_modes=True
                                ):

        if use_nonnegative_modes:
            neg_modes = np.where(nms.freqs < 1e-4 / UnitsData.hartrees_to_wavenumbers)
            if len(neg_modes) > 0 and len(neg_modes[0]) > 0:
                neg_modes = neg_modes[0]
                rem = np.setdiff1d(np.arange(len(nms.freqs)), neg_modes)
                nms = nms[rem]
                if mode_labels is not None:
                    if hasattr(mode_labels, 'items'):
                        mode_labels = dev.dict_take(mode_labels, rem)
                    else:
                        mode_labels = [mode_labels[r] for r in rem]

        if isinstance(internals, dict):

            if mode_labels is not None:
                # for i, (freq, lab) in enumerate(zip(reversed(nms.freqs), reversed(mode_labels))):
                #     print(
                #         "Mode {} ({}): {:.0f} {}".format(i + 1,
                #                                          len(nms.freqs) - (i+1),
                #                                          freq * UnitsData.hartrees_to_wavenumbers,
                #                                          "mixed"
                #                                          if lab.type is None else
                #                                          lab.type
                #                                          )
                #     )

                internals = internals.copy()
                mod_internals = {}
                for coord, c in internals.items():
                    if (
                            len(c) == 2
                            and not (c[0] is None or isinstance(c[0], str))
                    ):
                        mod_internals[coord] = (
                            c[0],
                            [
                                i for i,(f,l) in enumerate(zip(nms.freqs, mode_labels))
                                if cls._match_mode_label(c[1], f, l)
                            ]
                        )
                internals.update(mod_internals)

                if localization_mode_spaces is not None:
                    localization_mode_spaces = {
                        k: (
                            lm
                                if all(nput.is_int(i) for i in lm) else
                            [
                                i
                                    for i, (f,l) in enumerate(zip(nms.freqs, mode_labels))
                                if cls._match_mode_label(lm, f, l.type if l is not None else l)
                            ]
                        )
                        for k,lm in localization_mode_spaces.items()
                    }

            dim = len(nms.freqs)
            spaces = {}
            subints = {}
            subspaces_used = False
            for n, (coord, c) in enumerate(internals.items()):
                if (
                        len(c) == 2
                        and all(nput.is_int(i) for i in c[1])
                ):
                    subspaces_used = True
                    key = tuple(sorted(c[1]))
                    if key not in spaces:
                        spaces[key] = {}
                    spaces[key][coord] = n
                    subints[coord] = c[0]
                elif localization_mode_spaces is not None:
                    space = localization_mode_spaces.get(coord)
                    if space is not None:
                        subspaces_used = True
                        key = tuple(space)
                        if key not in spaces:
                            spaces[key] = {}
                        spaces[key][coord] = n
                        subints[coord] = c
                    elif isinstance(c, str):
                        space = localization_mode_spaces.get(c)
                        if space is not None:
                            subspaces_used = True
                            key = tuple(space)
                            if key not in spaces:
                                spaces[key] = {}
                            spaces[key][coord] = n
                            subints[coord] = c
                        else:
                            key = tuple(range(dim))
                            if key not in spaces:
                                spaces[key] = {}
                            spaces[key][coord] = n
                            subints[coord] = c
                    else:
                        for i in range(len(c)):
                            sub = tuple(c[i:])
                            space = localization_mode_spaces.get(sub)
                            if space is not None:
                                subspaces_used = True
                                key = tuple(space)
                                if key not in spaces:
                                    spaces[key] = {}
                                spaces[key][coord] = n
                                subints[coord] = c
                                break
                        else:
                            key = tuple(range(dim))
                            if key not in spaces:
                                spaces[key] = {}
                            spaces[key][coord] = n
                            subints[coord] = c

                else:
                    key = tuple(range(dim))
                    if key not in spaces:
                        spaces[key] = {}
                    spaces[key][coord] = n
                    subints[coord] = c

            # noneg_selection = None
            # if use_nonnegative_modes:
            #     noneg_selection = np.where(nms.freqs < 1e-4 / UnitsData.hartrees_to_wavenumbers)
            #     if len(noneg_selection) > 0 and len(noneg_selection[0]) > 0:
            #         subspaces_used = True
            #         noneg_selection= noneg_selection[0]


            internals = subints
            if subspaces_used:
                ord = []
                tranfs = []
                for space, ints in spaces.items():
                    space = np.array([
                        dim + s if s < 0 else s
                        for s in space
                    ])
                    # if noneg_selection is not None:
                    #     space = np.setdiff1d(space, noneg_selection)

                    subnms = nms[space]
                    sublocs = subnms.localize(internals=ints)
                    tf = np.zeros((dim, len(ints)))
                    inv = np.zeros((len(ints), dim))
                    t_sub, inv_sub = sublocs.localizing_transformation
                    tf[space, :] = t_sub
                    inv[:, space] = inv_sub
                    tranfs.append((tf, inv))
                    ord.extend(ints.values())

                tf = np.concatenate([t for t, i in tranfs], axis=1)
                inv = np.concatenate([i for t, i in tranfs], axis=0)
                ord = np.argsort(ord)
                tf = tf[:, ord]
                inv = inv[ord, :]
                loc_modes = nms.apply_transformation((tf, inv))
            else:
                loc_modes = nms.localize(internals=internals)
        else:
            # noneg_selection = None
            # if use_nonnegative_modes:
            #     noneg_selection = np.where(nms.freqs < 1e-4 / UnitsData.hartrees_to_wavenumbers)
            #     if len(noneg_selection) > 0 and len(noneg_selection[0]) > 0:
            #         noneg_selection = noneg_selection[0]
            #     else:
            #         noneg_selection = None
            #
            # if noneg_selection is not None:
            #     dim = len(nms.freqs)
            #     subnms = nms[noneg_selection]
            #     sublocs = subnms.localize(internals=internals)
            #     tf = np.zeros((dim, len(internals)))
            #     inv = np.zeros((len(internals), dim))
            #     t_sub, inv_sub = sublocs.localizing_transformation
            #     tf[noneg_selection, :] = t_sub
            #     inv[:, noneg_selection] = inv_sub
            #     loc_modes = nms.apply_transformation((tf, inv))
            # else:
            loc_modes = nms.localize(internals=internals)

        return loc_modes

    @classmethod
    def from_modes(cls,
                   nms,
                   internals=None,
                   mode_labels=None,
                   localization_mode_spaces=None,
                   oblique=True,
                   include_complement=False,
                   dipole_derivatives=None,
                   localize=True,
                   use_nonnegative_modes=True,
                   **opts):
        from ..Modes import NormalModes
        nms:NormalModes

        nms = nms.remove_frequency_scaling().remove_mass_weighting()
        if localize:
            loc_modes = cls.localize_internal_modes(
                nms,
                internals,
                localization_mode_spaces=localization_mode_spaces,
                mode_labels=mode_labels,
                use_nonnegative_modes=use_nonnegative_modes
            )

            if oblique:
                loc_modes = loc_modes.make_oblique()
            if include_complement:
                loc_modes = loc_modes.get_complement(concatenate=True)
                if internals is not None:
                    k = len(internals)
                    M = len(loc_modes.freqs)
                    if k < M:
                        if isinstance(internals, dict):
                            internals = internals.copy()
                            for i in range(k, M):
                                internals[i] = coordops.mode_label(
                                    None,
                                    None,
                                    None,
                                    "complement"
                                )
                        else:
                            internals = list(internals) + list(range(k, M))
        else:
            loc_modes = nms

        if dipole_derivatives is not None:
            # freq_conv = np.sqrt(loc_modes.compute_freqs())
            L = loc_modes.coords_by_modes
            if not oblique:
                a = np.diag(1 / np.sqrt(loc_modes.local_freqs))
                L = a @ L
            dipole_derivatives = [ L @ d for d in dipole_derivatives]


        # loc_modes = loc_modes.make_mass_weighted()
        f = loc_modes.local_hessian
        g = loc_modes.local_gmatrix


        return cls(
            f, g,
            internals=internals,
            dipole_derivatives=dipole_derivatives,
            modes=loc_modes,
            **opts
        )

    @classmethod
    def state(cls, nquanta_or_state, *tags, indices=None):
        if isinstance(nquanta_or_state, util.product_state):
            return nquanta_or_state
        nquanta = nquanta_or_state
        if nput.is_int(nquanta):
            if all(nput.is_int(q) for q in tags):
                return (nquanta,) + tags
        elif isinstance(nquanta, str):
            tags = (nquanta,) + tags
            nquanta = None
        elif nquanta is not None:
            if (
                    all(q is None or isinstance(q, str) for q in nquanta)
                    and len(tags) == 0
            ):
                return cls.state(None, *nquanta, indices=indices)
            elif nput.is_int(nquanta[0]):
                return cls.state(*nquanta, *tags, indices=indices)
            elif len(tags) == 0:
                base_data = [cls.state(n) for n in nquanta]
                indices = tuple(b.indices[0] if b.indices is not None else None for b in base_data)
                nquanta = tuple(b.quanta[0] if b.quanta is not None else None for b in base_data)
                tags = tuple(b.tags[0] if b.tags is not None else None for b in base_data)
            else:
                base_data = [cls.state(nquanta)] + [cls.state(t) for t in tags]
                indices = tuple(b.indices[0] if b.indices is not None else None for b in base_data)
                nquanta = tuple(b.quanta[0] if b.quanta is not None else None for b in base_data)
                tags = tuple(b.tags[0] if b.tags is not None else None for b in base_data)

        if isinstance(tags, str):
            tags = ((tags,),)
        elif all(q is None or isinstance(q, str) for q in tags):
            tags = (tags,)
        if indices is not None:
            if nput.is_int(indices):
                indices = (indices,)
            elif all(i is None for i in indices):
                indices = None
        if nquanta is not None:
            if nput.is_int(nquanta):
                nquanta = (nquanta,)
            elif all(i is None for i in nquanta):
                nquanta = None
        return util.product_state(indices, nquanta, tags)

    @classmethod
    def state_pair(cls, shared_atoms, state_1, state_2=None):
        if not (
                shared_atoms is None
                or nput.is_int(shared_atoms)
                or all(nput.is_int(a) for a in shared_atoms)
        ) and state_2 is None:
            state_1, state_2 = shared_atoms, state_1
            shared_atoms = None
        if shared_atoms is not None:
            return (
                shared_atoms,
                cls.state(state_1),
                cls.state(state_2),
            )
        else:
            return (
                cls.state(state_1),
                cls.state(state_2),
            )

    @classmethod
    def from_molecule(cls,
                      mol,
                      modes=None,
                      mode_labels=None,
                      internals=None,
                      coordinate_filter=None,
                      allowed_coordinate_types=None,
                      excluded_coordinate_types=None,
                      allowed_ring_types=None,
                      excluded_ring_types=None,
                      allowed_group_types=None,
                      excluded_group_types=None,
                      include_stretches=True,
                      include_bends=True,
                      include_dihedrals=False,
                      dipole_derivatives=None,
                      include_dipole=True,
                      internal_coordinate_sorting=None,
                      prune_excess_internals=False,
                      **opts):
        from ..Molecools import Molecule
        mol: Molecule
        if modes is None:
            modes = mol.get_normal_modes(use_internals=False)


        if internals is None:
            internals = mol.get_labeled_internals(
                coordinate_filter=coordinate_filter,
                allowed_coordinate_types=allowed_coordinate_types,
                excluded_coordinate_types=excluded_coordinate_types,
                allowed_ring_types=allowed_ring_types,
                excluded_ring_types=excluded_ring_types,
                allowed_group_types=allowed_group_types,
                excluded_group_types=excluded_group_types,
                include_stretches=include_stretches,
                include_bends=include_bends,
                include_dihedrals=include_dihedrals,
                coordinate_sorting=internal_coordinate_sorting,
                pruning=prune_excess_internals
            )

        if mode_labels is True:
            mode_labels = mol.get_mode_labels(
                # internals,
                modes=modes,
                pruning=prune_excess_internals,
                use_redundants=not prune_excess_internals
            )

        if dipole_derivatives is None and include_dipole:
            if mol.dipole_derivatives is not None:
                dipole_derivatives = mol.dipole_derivatives[1].T

        return cls.from_modes(
            modes,
            internals=internals,
            mode_labels=mode_labels,
            dipole_derivatives=dipole_derivatives,
            **opts
        )

class LocalHarmonicWavefunctions(MatrixWavefunctions):
    def __init__(self, model, energies, coefficients, **etc):
        super().__init__(energies, coefficients, **etc)
        self.model = model

    def get_modification_dict(self,
                              *,
                              model=None,
                              energies=None,
                              wavefunctions=None,
                              **etc
                              ):
        return dict(
            dict(model=model if model is not None else self.model),
            **super().get_modification_dict(
                energies=energies,
                wavefunctions=wavefunctions,
                **etc
            )
        )

    def modify(self,
               *,
               model=None,
               energies=None,
               wavefunctions=None,
               **etc
               ):
        mod_dict = self.get_modification_dict(
            model=model,
            energies=energies,
            wavefunctions=wavefunctions,
            **etc
        )
        return type(self)(
            mod_dict.pop('model'),
            mod_dict.pop('energies'),
            mod_dict.pop('wavefunctions'),
            **mod_dict
        )

class CustomLocalHarmonicModel(LocalHarmonicModel):
    default_localization_settings = {}
    default_molecule_settings = {}
    default_scalings = {}
    default_couplings = {}
    default_shifts = {}

    def __init__(self,
                 f, g,
                 internals=None,
                 anharmonic_scalings=None,
                 anharmonic_couplings=None,
                 anharmonic_shifts=None,
                 dipole_derivatives=None,
                 **etc
                 ):
        if anharmonic_scalings is None:
            anharmonic_scalings = {}
        anharmonic_scalings = dev.merge_dicts(self.default_scalings, anharmonic_scalings)
        if anharmonic_couplings is None:
            anharmonic_couplings = {}
        anharmonic_couplings = dev.merge_dicts(self.default_couplings, anharmonic_couplings)
        if anharmonic_shifts is None:
            anharmonic_shifts = {}
        anharmonic_shifts = dev.merge_dicts(self.default_shifts, anharmonic_shifts)

        super().__init__(
            f, g,
            internals=internals,
            anharmonic_scalings=anharmonic_scalings,
            anharmonic_couplings=anharmonic_couplings,
            anharmonic_shifts=anharmonic_shifts,
            dipole_derivatives=dipole_derivatives,
            **etc
        )

    @classmethod
    def from_modes(cls,
                   nms,
                   internals=None,
                   dipole_derivatives=None,
                   **opts):
        opts = dev.merge_dicts(cls.default_localization_settings, opts)
        return super().from_modes(
            nms,
            internals=internals,
            dipole_derivatives=dipole_derivatives,
            **opts
        )

    @classmethod
    def from_molecule(cls,
                      mol,
                      modes=None,
                      internals=None,
                      dipole_derivatives=None,
                      **opts):
        opts = dev.merge_dicts(cls.default_molecule_settings, opts)
        return super().from_molecule(
            mol,
            modes=modes,
            internals=internals,
            dipole_derivatives=dipole_derivatives,
            **opts
        )


class TaborCHModel(CustomLocalHarmonicModel):
    default_localization_settings = {"oblique":False}
    default_molecule_settings = dict(
        allowed_coordinate_types={'CH', 'HCH', 'NH', 'OH'},
        excluded_ring_types={'benzene'}
    )
    default_scalings = {
        LocalHarmonicModel.state("NH", "stretch"): 0.961,
        LocalHarmonicModel.state("methyl", "CH", "stretch"): 0.961,
        LocalHarmonicModel.state("ethyl", "CH", "stretch"): 0.961,
        LocalHarmonicModel.state("HCH", "bend"): 0.975,  # fundamentals and overtones
        LocalHarmonicModel.state(("HCH", "bend"), ("HCH", "bend")): 0.975,  # combination band
    }
    default_couplings = {
        # CH to adjacent bend overtone
        LocalHarmonicModel.state_pair(
            2,  # shared atoms
            ("CH", "stretch"),  # stretch fundamental
            (2, ("HCH", "bend"))  # bend overtone
        ): 22 / UnitsData.hartrees_to_wavenumbers,

        # CH3 to opposite bend combination
        LocalHarmonicModel.state_pair(
            ((2, 1),),  # shared atoms
            ("CH", "stretch"),  # stretch fundamental
            (("HCH", "bend"), ("HCH", "bend"))  # bend overtone
        ): 5.6 / UnitsData.hartrees_to_wavenumbers,
        LocalHarmonicModel.state_pair(
            ((1, 2),),  # shared atoms
            ("CH", "stretch"),  # stretch fundamental
            (("HCH", "bend"), ("HCH", "bend"))  # bend overtone
        ): 5.6 / UnitsData.hartrees_to_wavenumbers,

        # CH3 to adjacent bend combination
        LocalHarmonicModel.state_pair(
            ((2, 2),),  # shared atoms
            ("CH", "stretch"),  # stretch fundamental
            (("HCH", "bend"), ("HCH", "bend"))  # bend overtone
        ): 1.5 / UnitsData.hartrees_to_wavenumbers,
        LocalHarmonicModel.state_pair(
            ((2, 2),),  # shared atoms
            ("CH", "stretch"),  # stretch fundamental
            (("HCH", "bend"), ("HCH", "bend"))  # bend overtone
        ): 1.5 / UnitsData.hartrees_to_wavenumbers
    }
    default_shifts = {
        # internal-matching localization overshoots mildly
        ("ethyl", "CH", "stretch"): -4 / UnitsData.hartrees_to_wavenumbers,
        # internal-matching localization requires less bend overtone shifting
        ("HCH", "bend"): -2 * 9.6 / UnitsData.hartrees_to_wavenumbers
    }