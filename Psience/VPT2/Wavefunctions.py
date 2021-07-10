"""
Provides classes to support wave functions coming from VPT calculations
"""

import numpy as np, itertools as ip, time, enum

from McUtils.Numputils import SparseArray
from McUtils.Data import UnitsData

from ..BasisReps import BasisStateSpace, SimpleProductBasis, ExpansionWavefunctions, ExpansionWavefunction

from .Terms import DipoleTerms
from .Hamiltonian import PerturbationTheoryCorrections

__all__ = [
    'PerturbationTheoryWavefunction',
    'PerturbationTheoryWavefunctions',
]


class PerturbationTheoryWavefunction(ExpansionWavefunction):
    """
    These things are fed the first and second order corrections
    """

    def __init__(self, mol, basis, corrections):
        """
        :param mol: the molecule the wavefunction is for
        :type mol: Molecule
        :param basis: the basis the expansion is being done in
        :type basis: RepresentationBasis
        :param corrections: the corrections to the terms
        :type corrections: PerturbationTheoryCorrections
        """
        self.mol = mol
        self.corrs = corrections
        self.rep_basis = basis
        super().__init__(self.corrs.energies, self.corrs.wavefunctions, None)

    @property
    def order(self):
        return self.corrs.order

    def expectation(self, operator, other):
        return NotImplemented

    @property
    def zero_order_energy(self):
        return self.corrs.energies[0]


class PerturbationTheoryWavefunctions(ExpansionWavefunctions):
    """
    These things are fed the first and second order corrections
    """

    def __init__(self, mol, basis, corrections,
                 modes=None,
                 mode_selection=None,
                 logger=None,
                 operator_settings=None
                 ):
        """
        :param mol: the molecule the wavefunction is for
        :type mol: Molecule
        :param basis: the basis the expansion is being done in
        :type basis: SimpleProductBasis
        :param corrections: the corrections to the terms
        :type corrections: PerturbationTheoryCorrections
        """
        self.mol = mol
        self.corrs = corrections
        self.modes = modes
        self.mode_selection = mode_selection
        self.rep_basis = basis  # temporary hack until I decided how to merge the idea of
        # AnalyticWavefunctions with a RepresentationBasis
        self._tm_dat = None
        self._dipole_terms = None
        self._dipole_partitioning = self.DipolePartitioningMethod.Standard
        self.logger = logger
        self.operator_settings = operator_settings if operator_settings is not None else {}
        self._order = None
        super().__init__(
            self.corrs.energies,
            self.corrs.wfn_corrections,
            None
        )

    def energies_to_order(self, order):
        return np.sum(self.corrs.energy_corrs[:, :order+1], axis=1)

    @property
    def order(self):
        if self._order is None:
            return self.corrs.order
        else:
            return self._order

    def expectation(self, operator, other):
        return NotImplemented

    @property
    def zero_order_energies(self):
        return self.corrs.energy_corrs[:, 0]

    def get_M0(self, mu_0):
        return self.rep_basis.representation(name='M(0)', coeffs=mu_0,
                                             **self.operator_settings
                                             )
    def get_M1(self, mu_1):
        return self.rep_basis.representation("x", name='M(1)', coeffs=mu_1,
                                             **self.operator_settings
                                             )
    def get_M2(self, mu_2):
        return 1 / 2 * self.rep_basis.representation("x", "x", name="M(2)", coeffs=mu_2,
                                                     **self.operator_settings
                                                     )
    def get_M3(self, mu_3):
        return 1 / 6 * self.rep_basis.representation("x", "x", "x", name="M(3)", coeffs=mu_3,
                                                     **self.operator_settings
                                                     )

    def _mu_representations(self,
                            a,
                            mu,
                            M,
                            space,
                            bra_spaces,
                            ket_spaces,
                            order,
                            k,
                            partitioning,
                            rep_inds
                            ):

        # define out reps based on partitioning style
        if partitioning is None:
            partitioning = self.dipole_partitioning
        elif not isinstance(partitioning, self.DipolePartitioningMethod):
            partitioning = self.DipolePartitioningMethod(partitioning)

        def get_states(m, filters):
            with self.logger.block(tag="getting coupled states for {}".format(m)):
                start = time.time()

                for i in range(order+1):
                    for j in range(order+1):
                        if i + j + k:
                            ...
                        filter = filters[j]


                inds, filter = bra_space.get_representation_brakets(
                    other=ket_space, filter=filter, return_filter=True,
                    selection_rules=m.selection_rules
                )  # no selection rules here
                end = time.time()
                self.logger.log_print('took {t:.3f}s...', t=end - start)
            return inds, filter

        filter = ket_space.to_single()
        if all(
                    isinstance(m, (np.ndarray, SparseArray)) and m.shape == (M, M)
                    or isinstance(m, (int, float, np.integer, np.floating)) and m == 0
                    for m in mu
        ): # we were passed representations to reuse
            mu_terms = mu
            if partitioning == self.DipolePartitioningMethod.Standard and len(mu_terms) < 4:
                dts = self.dipole_terms[a]
                mu_3 = dts[3]
                m3 = self.get_M3(mu_3)
                if rep_inds[3] is None:
                    m3_inds, filter = get_states(m3, filter)
                    with self.logger.block(tag="building {h}".format(h=m3)):
                        rep_inds[3] = m3_inds
                        rep_3 = m3.get_representation_matrix(rep_inds[3], self.corrs.total_basis)
                mu_terms.append(rep_3)
                # mu_terms = mu_terms + [rep_3]
        else:
            mu_0, mu_1, mu_2, mu_3 = mu

            m0 = self.get_M0(mu_0)
            if rep_inds[0] is None:
                m0_inds, filter = get_states(m0, filter)
                rep_inds[0] = m0_inds
            m1 = self.get_M1(mu_1)
            if rep_inds[1] is None:
                m1_inds, filter = get_states(m1, filter)
                rep_inds[1] = m1_inds
            m2 = self.get_M2(mu_2)
            if rep_inds[2] is None:
                m2_inds, filter = get_states(m2, filter)
                rep_inds[2] = m2_inds

            if partitioning == self.DipolePartitioningMethod.Intuitive:
                reps = [m0, m1, m2]
            elif partitioning == self.DipolePartitioningMethod.Standard:
                m3 = self.get_M3(mu_3)
                if rep_inds[3] is None:
                    m3_inds, filter = get_states(m3, filter)
                    rep_inds[3] = m3_inds
                reps = [m0, m1, m2, m3]
            else:
                raise ValueError("don't know how to interpret dipole partitioning {}".format(partitioning))

            # reps = [m1, m2, m3]
            mu_terms = [None] * len(reps)
            for i, h in enumerate(reps):
                m_pairs = rep_inds[i]
                if m_pairs is None:
                    raise ValueError("representation indices haven't filled enough to calculate {}".format(h))

                with self.logger.block(tag="building {}".format(h)):
                    start = time.time()
                    sub = h.get_representation_matrix(m_pairs, self.corrs.total_basis, zero_element_warning=False) # expect zeros in M(3)?
                    end = time.time()
                    self.logger.log_print('took {t:.3f}s...', t=end - start)
                # sub = generate_rep(h, m_pairs)
                mu_terms[i] = sub

        return mu_terms

    def _transition_moments(self,
                            mu_x, mu_y, mu_z,
                            correction_terms=None,
                            lower_states=None,
                            excited_states=None,
                            degenerate_transformation=None,
                            partitioning=None,
                            order=None
                            ):
        """
        Calculates the x, y, and z components of the transition moment between
        the wavefunctions stored.
        By default, only calculates moments involving the ground states.

        :param mu_x: dipole x components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_x: Iterable[np.ndarray]
        :param mu_y: dipole y components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_y: Iterable[np.ndarray]
        :param mu_z: dipole z components (1st, 2nd, 3rd derivatives in normal modes)
        :type mu_z: Iterable[np.ndarray]
        :param partitioning: whether to put linear and constant in mu_0 (standard) or linear in mu_1 (intuitive)
        :type partitioning: None | str
        :return:
        :rtype:
        """

        if order is None:
            order = self.order

        logger = self.logger
        with logger.block(tag="Calculating intensities:"):

            space = self.corrs.coupled_states
            # store this for handling the degenerate_transformation
            lower_states_input = (0,) if lower_states is None else lower_states
            upper_states_input = tuple(range(space.nstates)) if excited_states is None else excited_states

            corr_terms = self.corrs.wfn_corrections
            if degenerate_transformation is None:
                degenerate_transformation = self.corrs.degenerate_transf
            if degenerate_transformation is not None:
                deg_transf_lower = degenerate_transformation[lower_states_input, :]
                deg_transf_upper = degenerate_transformation[upper_states_input, :]

                corr_terms_lower = [deg_transf_lower.dot(corr_terms[i]) for i in range(order)]
                corr_terms_upper = [deg_transf_upper.dot(corr_terms[i]).T for i in range(order)]

                correction_terms = corr_terms_lower, corr_terms_upper

                bra_space = space.to_single().take_unique()
                ket_space = space.to_single().take_unique()

            if correction_terms is None:
                if lower_states is None:
                    low_spec = (0,)
                else:
                    if isinstance(lower_states, (int, np.integer)):
                        lower_states = (lower_states,)
                    low_spec = lower_states
                lower_states = space[(low_spec,)]

                if excited_states is None:
                    up_spec = tuple(range(space.nstates))
                    excited_states = space
                else:
                    if isinstance(excited_states, (int, np.integer)):
                        excited_states = (excited_states,)
                    up_spec = excited_states
                    excited_states = space[(up_spec,)]

                bra_space = lower_states if isinstance(lower_states,
                                                       BasisStateSpace) else lower_states.to_single().take_unique()
                ket_space = excited_states if isinstance(excited_states,
                                                         BasisStateSpace) else excited_states.to_single().take_unique()

                # M = len(space.indices)
                logger.log_print(
                    [
                        "lower/upper states: {l}/{u}"
                    ],
                    l=len(low_spec),
                    u=len(up_spec)
                )

                with logger.block(tag='taking correction subspaces:'):
                    start = time.time()
                    logger.log_print('getting ground space...')
                    corr_terms_lower = [corr_terms[i][low_spec, :] for i in range(order)]
                    logger.log_print(
                        [
                            "highest-order space: {s}",
                        ],
                        s=corr_terms_lower[-1]
                    )

                    logger.log_print('getting excited space...')
                    corr_terms_upper = [corr_terms[i][up_spec, :].T for i in range(order)]
                    logger.log_print(
                        [
                            "highest-order space: {s}",
                        ],
                        s=corr_terms_upper[-1]
                    )

                    end = time.time()
                    logger.log_print(
                        [
                            "took {t:.3f}s..."
                        ],

                        t=end - start
                    )
            else:
                corr_terms_lower, corr_terms_upper = correction_terms

            rep_inds = [None] * 4
            M = corr_terms[0].shape[1]
            mu_reps = []
            total_space = self.corrs.total_basis
            # corr_vecs = self.corrs.wfn_corrections#[..., M]
            transition_moment_components = np.zeros(
                (order, 3)).tolist()  # x, y, and z components of the 0th, 1st, and 2nd order stuff

            mu = [mu_x, mu_y, mu_z]

            # define out reps based on partitioning style
            if partitioning is None:
                partitioning = self.dipole_partitioning
            elif not isinstance(partitioning, self.DipolePartitioningMethod):
                partitioning = self.DipolePartitioningMethod(partitioning)

            logger.log_print(
                [
                    "dipole partitioning: {m}"
                ],
                m=partitioning.value
            )

            # raise Exception(corr_terms[1], corr_terms_lower[1], corr_terms_upper[1], low_spec, len(low_spec), len(up_spec))

            for a in range(3):  # x, y, and z
                with logger.block(tag="Getting M representations for axis {}:".format(a)):
                    start = time.time()
                    mu_terms = self._mu_representations(
                        a,
                        mu[a],
                        M,
                        total_space,
                        bra_space,
                        ket_space,
                        partitioning,
                        rep_inds
                    )
                    end = time.time()
                    logger.log_print(
                        [
                            "took {t:.3f}s..."
                        ],
                        t=end-start
                    )

                # print(">>>>", mu_terms)
                mu_reps.append(mu_terms)

                if partitioning == self.DipolePartitioningMethod.Intuitive:
                    mu_terms = [mu_terms[0], mu_terms[1], mu_terms[2]]
                    # print(mu_terms)
                elif partitioning == self.DipolePartitioningMethod.Standard:
                    mu_terms = [mu_terms[0] + mu_terms[1], mu_terms[2], mu_terms[3]]
                else:
                    raise ValueError("don't know how to interpret dipole partitioning {}".format(partitioning))

                # define out reps based on partitioning style
                with logger.block(tag="axis: {}".format(a)):
                    logger.log_print(
                        [
                            "non-zero dipole terms: {nt}"
                        ],
                        nt=len(mu_terms) - mu_terms.count(0)
                    )

                    cur_lw = np.get_printoptions()['linewidth']
                    try:
                        np.set_printoptions(linewidth=100000) #infinite line width basically...

                        with logger.block(tag="calculating corrections..."):
                            start = time.time()
                            for q in range(order):  # total quanta
                                # logger.log_print("calculating corrections at order {q}...", q=q)
                                terms = []
                                # should do this smarter
                                for i, j, k in ip.product(range(q + 1), range(q + 1), range(q + 1)):
                                    if i + j + k == q:
                                        if len(mu_terms) <= k or len(corr_terms) <= i or len(corr_terms) <= j:
                                            new = np.zeros((len(low_spec), len(up_spec)))
                                        else:
                                            m = mu_terms[k]
                                            if isinstance(m, (int, float, np.integer, np.floating)) and m == 0:
                                                # to make it easy to zero stuff out
                                                new = np.zeros((len(lower_states_input), len(upper_states_input)))
                                            else:
                                                c_lower = corr_terms_lower[i]
                                                c_upper = corr_terms_upper[j]
                                                num = c_lower.dot(m)
                                                new = num.dot(c_upper)
                                            if isinstance(new, SparseArray):
                                                new = new.asarray()
                                            new = new.reshape((len(lower_states_input), len(upper_states_input)))
                                        with logger.block(tag="<{i}|M({k})|{j}>".format(i=i, j=j, k=k)):
                                            logger.log_print(
                                                str(new).splitlines()
                                            )
                                        terms.append(new)
                                        # raise Exception(new.toarray())
                                # print(q, a)
                                transition_moment_components[q][a] = terms
                    finally:
                        np.set_printoptions(linewidth=cur_lw)

                        end = time.time()
                        logger.log_print(
                            "took {t:.3f}s",
                            t=end - start
                        )

            # we calculate it explicitly like this up front in case we want to use it later since the shape
            # can be a bit confusing ([0-order, 1-ord, 2-ord], [x, y, z])
            tmom = [
                sum(
                    sum(ip.chain(transition_moment_components[i][j]))
                    if not isinstance(transition_moment_components[i][j], (float, int, np.floating, np.integer))
                    else transition_moment_components[i][j]
                    for i in range(order)  # correction order
                ) for j in range(3)  # xyz
            ]

            # mu_reps already organized like x/y/z
            # mu_reps = [
            #     [m[a] for m in mu_reps]
            #     for a in range(3)
            # ]

        return [tmom, transition_moment_components, mu_reps, (corr_terms_lower, corr_terms_upper)]

    class TermHolder(tuple):
        """symbolic wrapper on Tuple so we can know that we've canonicalized some term"""
    @property
    def dipole_terms(self):
        if self._dipole_terms is None:
            self._dipole_terms = self.TermHolder(DipoleTerms(self.mol, modes=self.modes, mode_selection=self.mode_selection).get_terms())
        elif not isinstance(self._dipole_terms, self.TermHolder):
            self._dipole_terms = self.TermHolder(DipoleTerms(self.mol, modes=self.modes, mode_selection=self.mode_selection, derivatives=self._dipole_terms).get_terms())
        return self._dipole_terms

    @dipole_terms.setter
    def dipole_terms(self, v):
        self._dipole_terms = v
        self._tm_dat = None

    class DipolePartitioningMethod(enum.Enum):
        Standard="standard"
        Intuitive="intuitive"

    @property
    def dipole_partitioning(self):
        return self._dipole_partitioning

    @dipole_partitioning.setter
    def dipole_partitioning(self, p):
        part = self.DipolePartitioningMethod(p)
        if part is not self._dipole_partitioning:
            self._dipole_partitioning = part
            if self._tm_dat is not None:
                self._tm_dat[0] = None

    def _load_tm_dat(self):
        if self._tm_dat is None:
            self._tm_dat = self._transition_moments(*self.dipole_terms)
        elif self._tm_dat[0] is None:
            self._tm_dat = self._transition_moments(*self._tm_dat[2], correction_terms=self._tm_dat[3])

    @property
    def transition_moments(self):
        """
        Computes the transition moments between wavefunctions stored in the object

        :return:
        :rtype:
        """

        self._load_tm_dat()
        return self._tm_dat[0]

    @property
    def transition_moment_corrections(self):
        """
        Computes the transition moment corrections between wavefunctions stored in the object

        :return:
        :rtype:
        """

        self._load_tm_dat()
        return self._tm_dat[1]

    @property
    def oscillator_strengths(self):
        """
        Computes the oscillator strengths for transitions from the ground state to the other states

        :return:
        :rtype:
        """

        tms = self.transition_moments
        return self._oscillator_strengths(tms)

    def oscillator_strengths_to_order(self, order):
        """

        :param tms:
        :type tms:
        :return:
        :rtype:
        """

        return self._oscillator_strengths(self.transition_moments_to_order(order))

    def _oscillator_strengths(self, tms):

        gs_tms = np.array([tms[i][0] for i in range(3)]).T
        osc = np.linalg.norm(gs_tms, axis=1) ** 2
        return osc

    @property
    def intensities(self):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        return self._intensities(self.oscillator_strengths)

    def intensities_to_order(self, order):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        return self._intensities(self.oscillator_strengths_to_order(order), energy_order=order)

    def _intensities(self, oscs, energy_order=None):
        if energy_order is None:
            eng = self.energies
        else:
            eng = self.energies_to_order(energy_order)
        units = UnitsData.convert("OscillatorStrength", "KilometersPerMole")
        return units * (eng - eng[0]) * oscs

    @property
    def zero_order_intensities(self):
        """
        Computes the harmonic intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        eng = self.zero_order_energies
        tm = np.array([
            self.transition_moment_corrections[0][x][0][0] for x in range(3)
        ]).T
        osc = np.linalg.norm(tm, axis=1) ** 2
        units = 3.554206329390961e6
        return units * (eng - eng[0]) * osc

    def generate_intensity_breakdown(self, include_wavefunctions=True):
        """
        Generates a breakdown of the terms that contribute to the intensity
        Returns in a format that can be directly exported to JSON if desired.

        :return:
        :rtype:
        """
        from collections import OrderedDict

        engs = self.energies
        freqs = engs - engs[0]
        # raise Exception(self.corrs.coupled_states.to_single().excitations)
        terms = OrderedDict((
            ("frequencies", freqs.tolist()),
            ('axes', self.mol.inertial_axes.tolist()),
            ("states", self.corrs.states.excitations),
            ('breakdowns', OrderedDict())
        ))

        if include_wavefunctions:
            terms["wavefunctions"] = {
                "corrections": [x.toarray().tolist() for x in self.corrs.wfn_corrections],
                "coupled_states": self.corrs.total_basis.excitations
            }
        bds = terms['breakdowns']

        # dts = self.dipole_terms

        dipole_breakdowns = OrderedDict((
            ("Full",      (True, True, True, True)),
            ("Constant",  (True, False, False, False)),
            ("Linear",    (False, True, False, False)),
            ("Quadratic", (False, False, True, False)),
            ("Cubic",     (False, False, False, True)),
            ("Order0",    (True, True, False, False)),
            ("Order1",    (True, True, True, False))
            # ("Order2",  (True, True, True, True)),
        ))

        # wfn_terms.append(freqs.tolist())
        for key in dipole_breakdowns:
            if self.logger is not None:
                self.logger.log_print(
                    "-" * 50 + "{k} {m}" + "-" * 50,
                    k=key,
                    m=self.dipole_partitioning.value,
                    padding="",
                    newline=" "
                )
            if key == "Full":
                # gotta come first
                ints = self.intensities
            else:
                mu_reps = self._tm_dat[2]
                dts = dipole_breakdowns[key]
                new_mu = [
                    [ m[i] if d and len(m) > i else 0 for i, d in enumerate(dts)]
                    for m in mu_reps # x/y/z
                ]
                tms = self._transition_moments(*new_mu)
                oscs = self._oscillator_strengths(tms[0])
                ints = self._intensities(oscs)

            full_corrs = self.transition_moment_corrections
            # raise Exception([np.array(x).shape for x in full_corrs])
            corrs_keys = [y for x in ip.chain(
                [(i, j, k) for i, j, k in ip.product(range(q + 1), range(q + 1), range(q + 1)) if i + j + k == q]
                for q in range(len(full_corrs))
            ) for y in x]

            full_corrs = np.concatenate([np.array(x) for x in full_corrs], axis=1)
            # transpose so we have x-y-z as outermost dimension and correction as innermost
            full_corrs = full_corrs.transpose((0, 2, 3, 1)).tolist()

            bds[key] = OrderedDict((
                ("intensities", ints),
                ("corrections", {
                    "keys": corrs_keys,
                    "values": OrderedDict(zip(['X', 'Y', 'Z'], full_corrs))
                })
            ))

        return terms

    @classmethod
    def write_CSV_breakdown(cls, file, intensity_breakdown, padding=None):
        """
        Writes an intensity breakdown to a CSV by annoyingly flattening all the arrays

        :param file:
        :type file:
        :param intensity_breakdown:
        :type intensity_breakdown:
        :return:
        :rtype:
        """
        import csv

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        if padding is None:
            padding = []
        else:
            padding = list(padding)

        if isinstance(file, str):
            with open(file, "w+") as fp:
                cls.write_CSV_breakdown(fp, intensity_breakdown)
        else:
            writer = csv.writer(file)
            # write header
            ib = intensity_breakdown
            bds = ib['breakdowns']
            keys = ['State', 'Freq.'] + list(bds.keys())
            writer.writerows([padding + keys])

            # write main intensity table
            states = ib['states']
            freqs = h2w * np.array(ib['frequencies'])
            intensity_table = []
            for i, s in enumerate(states):
                f = freqs[i]

                ints = [b['intensities'][i] for b in bds.values()]
                intensity_table.append([s, f] + ints)
            writer.writerows([padding + x for x in intensity_table])

            for k, b in bds.items():

                corr_block = [[k]]
                corrs = b['corrections']
                corr_keys = ["<gs^({0})|mu^({2})|es^({1})>".format(*lab) for lab in corrs['keys']]
                for ax, vals in corrs['values'].items():
                    corr_block.append([ax])
                    corr_block.append([" "] + corr_keys)
                    for s, v in zip(states, vals[0]):
                        corr_block.append([s] + list(v))
                writer.writerows([padding + x for x in corr_block])

            if 'wavefunctions' in ib:

                coupled_state_blocks = [["Wavefunction Corrections"]]
                wfn_corrs = np.array(ib['wavefunctions']["corrections"])
                coupled_states = ib['wavefunctions']["coupled_states"]

                # raise Exception(coupled_states)
                cs_keys = list(coupled_states)
                for corr_block, state in zip(wfn_corrs.transpose((1, 0, 2)), states):
                    coupled_state_blocks.append([state] + cs_keys)
                    for i, b in enumerate(corr_block):
                        coupled_state_blocks.append(["|n^({})>".format(i)] + list(b))

                writer.writerows([padding + x for x in coupled_state_blocks])

    @classmethod
    def _format_energies_table(cls, states, zpe, freqs, real_fmt='{:>10.3f}', dash_fmt='{:>10}'):

        n_modes = len(states[0])

        zpe = np.asanyarray(zpe)
        freqs = np.asanyarray(freqs)

        ncols = freqs.shape[1] if len(freqs.shape) > 0 else 1
        states_fmt = '{:<1.0f} ' * n_modes
        padding = len(states_fmt.format(*[0 for _ in range(n_modes)]))
        bar = "  Harmonic Anharmonic"
        if zpe is not None:
            header = ["State" + " " * (padding - len("State")) + bar + " " + bar]
            header += [
                " " * padding
                + " " * (len("  Harmonic") - 3) + "ZPE"
                + " " + " " * (len("Anharmonic") - 3) + "ZPE"
                + " "
                + " " * (len("  Harmonic") - 9) + "Frequency"
                + " " + " " * (len("Anharmonic") - 9) + "Frequency"
            ]
            zpe_fmt = (real_fmt + " ") * ncols + (dash_fmt + " ") * ncols
            freq_fmt = (dash_fmt + " ") * ncols + (real_fmt + " ") * ncols
            lines = [
                        states_fmt.format(*[0] * n_modes) + zpe_fmt.format(*zpe, *["-"] * ncols)] + [
                        states_fmt.format(*s) + freq_fmt.format(*["-"] * ncols, *e)
                        for s, e in
                        zip(states[1:], freqs)
                    ]
        else:
            header = ["State" + " " * (padding - len("State")) + bar]
            header += [
                " " * padding
                + " " * (len("  Harmonic") - 9) + "Frequency"
                + " " + " " * (len("Anharmonic") - 9) + "Frequency"
            ]
            freq_fmt = ncols + (real_fmt + " ") * ncols
            lines = [
                states_fmt.format(*s) + freq_fmt.format(*e)
                for s, e in
                zip(states, freqs)
            ]

        return "\n".join([
            *header,
            *lines,
        ])


    def format_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>10.3f}', dash_fmt='{:>10}'):
        # simple utility function pulled from the unit tests

        if states is None:
            states = self.corrs.states.excitations

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        harm_engs = self.zero_order_energies * h2w
        engs = self.energies * h2w
        if zpe is None and np.all(states[0] == 0):
            zpe = [harm_engs[0], engs[0]]

        if freqs is None:
            harm_freq = harm_engs[1:] - harm_engs[0]
            anh_freqs = engs[1:] - engs[0]
            freqs = np.column_stack([harm_freq, anh_freqs])

        return self._format_energies_table(
            states,
            zpe,
            freqs,
            real_fmt=real_fmt,
            dash_fmt=dash_fmt
        )



    def format_intensities_table(self):
        # simple utility function pulled from the unit tests

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        states = self.corrs.states.excitations

        engs = h2w * self.energies
        freqs = engs - engs[0]
        ints = self.intensities

        harm_engs = h2w * self.zero_order_energies
        harm_freqs = harm_engs - harm_engs[0]
        harm_ints = self.zero_order_intensities

        n_modes = self.corrs.total_basis.ndim
        padding = np.max([len(str("0 " * n_modes)), 1]) + 1
        bar = "Frequency    Intensity"
        spacer = "    "
        report = (
                         " " * (padding + 2) + " " * (len("Frequency") - 2) + "Harmonic" + " " * (
                             len("    Intensity") + 2 - len("Harmonic")) + spacer + " " * (
                                     len("Frequency") - 2) + "Anharmonic\n"
                                                             "State" + " " * (padding - 3) + bar + spacer + bar + "\n"
                 ) + "\n".join(
            (
                (
                        "  " + ("{:<1.0f} " * n_modes) + " "
                        + "{:>9.2f} {:>12.4f}"
                        + spacer
                        + "{:>9.2f} {:>12.4f}"
                ).format(*s, hf, hi, f, i)
                for s, hf, hi, f, i in zip(states[1:], harm_freqs[1:], harm_ints[1:], freqs[1:], ints[1:])
            ))

        return report