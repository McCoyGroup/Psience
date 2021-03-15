"""
Provides classes to support wave functions coming from VPT calculations
"""

import numpy as np, itertools as ip, time, enum

from McUtils.Numputils import SparseArray

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

    def __init__(self, mol, basis, corrections, logger=None):
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
        self.rep_basis = basis  # temporary hack until I decided how to merge the idea of
        # AnalyticWavefunctions with a RepresentationBasis
        self._tm_dat = None
        self._dipole_terms = None
        self._dipole_partitioning = self.DipolePartitioningMethod.Standard
        self.logger = logger
        self._order = None
        super().__init__(
            self.corrs.energies,
            self.corrs.wfn_corrections,
            None
        )

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

    def _build_representation_matrix(self, h, m_pairs):
        """
        Constructs a representation of an operator in a coupled
        space...

        :param h:
        :type h:
        :param cs:
        :type cs:
        :return:
        :rtype:
        """
        logger = self.logger

        if len(m_pairs) > 0:
            logger.log_print(["coupled space dimension {d}"], d=len(m_pairs))
            sub = h[m_pairs]
            if isinstance(sub, SparseArray):
                sub = sub.asarray()
            SparseArray.clear_cache()
        else:
            logger.log_print('no states to couple!')
            sub = 0

        logger.log_print("constructing sparse representation...")

        N = len(self.corrs.total_basis)
        if isinstance(sub, (int, np.integer, np.floating, float)):
            if sub == 0:
                sub = SparseArray.empty((N, N), dtype=float)
            else:
                raise ValueError("Using a constant shift of {} will force Hamiltonians to be dense...".format(sub))
                sub = np.full((N, N), sub)
        else:
            # figure out the appropriate inds for this data in the sparse representation
            row_inds = self.corrs.total_basis.find(m_pairs.bras)
            col_inds = self.corrs.total_basis.find(m_pairs.kets)

            # upper triangle of indices
            up_tri = np.array([row_inds, col_inds]).T
            # lower triangle is made by transposition
            low_tri = np.array([col_inds, row_inds]).T
            # but now we need to remove the duplicates, because many sparse matrix implementations
            # will sum up any repeated elements
            full_inds = np.concatenate([up_tri, low_tri])
            full_dat = np.concatenate([sub, sub])

            _, idx = np.unique(full_inds, axis=0, return_index=True)
            sidx = np.sort(idx)
            full_inds = full_inds[sidx]
            full_dat = full_dat[sidx]
            sub = SparseArray.from_data((full_dat, full_inds.T), shape=(N, N))

        return sub

    @staticmethod
    def _generate_rep(h, m_pairs, logger, i, M, space):
            # m_pairs = rep_inds[i]

            # print(">", len(m_pairs[0]))
            if logger is not None:
                start = time.time()
                logger.log_print(
                    [
                        "calculating M^({})...",
                        "(coupled space size {})"
                    ],
                    i,
                    len(m_pairs[0])
                )
            # print(m_pairs)
            sub = h[m_pairs[0], m_pairs[1]]
            if isinstance(sub, SparseArray):
                sub = sub.asarray()  #
            SparseArray.clear_ravel_caches()
            if logger is not None:
                end = time.time()
                logger.log_print(
                    [
                        "took {}s"
                    ],
                    round(end - start, 3)
                )
            if logger is not None:
                logger.log_print(
                    "constructing sparse representation..."
                )

            if isinstance(sub, (int, np.integer, np.floating, float)):
                if sub == 0:
                    sub = SparseArray.empty((M, M), dtype=float)
                else:
                    raise ValueError(
                        "Using a constant shift of {} will force representations to be dense...".format(sub))
                    sub = np.full((N, N), sub)
            elif len(sub) == ():
                sub = SparseArray.empty((M, M), dtype=float)
            else:
                # figure out the appropriate inds for this data in the sparse representation
                row_inds = space.find(m_pairs[0])
                col_inds = space.find(m_pairs[1])

                # upper triangle of indices
                up_tri = np.array([row_inds, col_inds]).T
                # lower triangle is made by transposition
                low_tri = np.array([col_inds, row_inds]).T
                # but now we need to remove the duplicates, because many sparse matrix implementations
                # will sum up any repeated elements
                full_inds = np.concatenate([up_tri, low_tri])
                full_dat = np.concatenate([sub, sub])

                _, idx = np.unique(full_inds, axis=0, return_index=True)
                sidx = np.sort(idx)
                full_inds = full_dat[sidx]
                full_dat = full_dat[sidx]
                sub = SparseArray((full_dat, full_inds.T), shape=(M, M))

            return sub

    def _mu_representations(self,
                            a,
                            mu,
                            M,
                            space,
                            bra_space,
                            ket_space,
                            partitioning,
                            rep_inds
                            ):

        # define out reps based on partitioning style
        if partitioning is None:
            partitioning = self.dipole_partitioning
        elif not isinstance(partitioning, self.DipolePartitioningMethod):
            partitioning = self.DipolePartitioningMethod(partitioning)

        logger = self.logger
        if all(
                    isinstance(m, (np.ndarray, SparseArray)) and m.shape == (M, M)
                    or isinstance(m, (int, float, np.integer, np.floating)) and m == 0
                    for m in mu
        ): # we were passed representations to reuse
            mu_terms = mu
            if partitioning == self.DipolePartitioningMethod.Standard and len(mu_terms) < 4:
                dts = self.dipole_terms[a]
                mu_3 = dts[3]
                m3 = 1 / 6 * self.rep_basis.representation("x", "x", "x", coeffs=mu_3)
                if rep_inds[3] is None:
                    m3_inds = bra_space.get_representation_brakets(
                        other=ket_space,
                        selection_rules=bra_space.basis.selection_rules("x", "x", "x")
                    )
                    rep_inds[3] = m3_inds
                with self.logger.block(tag="getting M^(3)"):
                    rep_3 = self._build_representation_matrix(m3, rep_inds[3])
                mu_terms.append(rep_3)
                # mu_terms = mu_terms + [rep_3]
        else:
            mu_0, mu_1, mu_2, mu_3 = mu

            m0 = self.rep_basis.representation(coeffs=mu_0)
            if rep_inds[0] is None:
                m0_inds = bra_space.get_representation_brakets(other=ket_space,
                                                               selection_rules=[[]])  # no selection rules here
                rep_inds[0] = m0_inds
            m1 = (
                    self.rep_basis.representation("x", coeffs=mu_1)
                    + self.rep_basis.representation(coeffs=mu_0)
            )
            if rep_inds[1] is None:
                m1_inds = bra_space.get_representation_brakets(
                    other=ket_space,
                    selection_rules=bra_space.basis.selection_rules("x")
                )
                rep_inds[1] = m1_inds
            m2 = 1 / 2 * self.rep_basis.representation("x", "x", coeffs=mu_2)
            if rep_inds[2] is None:
                m2_inds = bra_space.get_representation_brakets(
                    other=ket_space,
                    selection_rules=bra_space.basis.selection_rules("x", "x"))
                rep_inds[2] = m2_inds

            if partitioning == self.DipolePartitioningMethod.Intuitive:
                reps = [m0, m1, m2]
            elif partitioning == self.DipolePartitioningMethod.Standard:
                m3 = 1 / 6 * self.rep_basis.representation("x", "x", "x", coeffs=mu_3)
                if rep_inds[3] is None:
                    m3_inds = bra_space.get_representation_brakets(
                        other=ket_space,
                        selection_rules=bra_space.basis.selection_rules("x", "x", "x")
                    )
                    rep_inds[3] = m3_inds
                reps = [m0, m1, m2, m3]
            else:
                raise ValueError("don't know how to interpret dipole partitioning {}".format(partitioning))

            # reps = [m1, m2, m3]
            mu_terms = [None] * len(reps)
            for i, h in enumerate(reps):
                m_pairs = rep_inds[i]
                if m_pairs is None:
                    raise ValueError("representation indices haven't filled enough to calculate {}".format(i))

                with self.logger.block(tag="getting M^({})".format(i)):
                    sub = self._build_representation_matrix(h, m_pairs)
                # sub = generate_rep(h, m_pairs)
                mu_terms[i] = sub

        return mu_terms

    def _transition_moments(self,
                            mu_x, mu_y, mu_z,
                            lower_states=None,
                            excited_states=None,
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

            if lower_states is None:
                low_spec = (0,)
            else:
                if isinstance(lower_states, (int, np.integer)):
                    lower_states = (lower_states,)
                low_spec = lower_states
            lower_states = space[low_spec]

            if excited_states is None:
                up_spec = tuple(range(space.nstates))
                excited_states = space
            else:
                if isinstance(excited_states, (int, np.integer)):
                    excited_states = (excited_states,)
                up_spec = excited_states
                excited_states = space[up_spec]

            bra_space = lower_states if isinstance(lower_states, BasisStateSpace) else lower_states.to_single()
            ket_space = excited_states if isinstance(excited_states, BasisStateSpace) else excited_states.to_single()

            # M = len(space.indices)
            logger.log_print(
                [
                    "lower/upper states: {l}/{u}"
                ],
                l=len(low_spec),
                u=len(up_spec)
            )

            # corr_vecs = self.corrs.wfn_corrections#[..., M]
            transition_moment_components = np.zeros((3, 3)).tolist()  # x, y, and z components of the 0th, 1st, and 2nd order stuff

            mu = [mu_x, mu_y, mu_z]
            rep_inds = [None] * 4
            corr_terms = self.corrs.wfn_corrections
            M = corr_terms[0].shape[1]
            mu_reps = []
            total_space = self.corrs.total_basis

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

            for a in range(3):  # x, y, and z

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

                    with logger.block(tag="calculating corrections..."):
                        start = time.time()
                        for q in range(order):  # total quanta
                            terms = []
                            for i, j, k in ip.product(range(q + 1), range(q + 1), range(q + 1)):
                                if i + j + k == q:
                                    if len(mu_terms) <= k:
                                        new = np.zeros((len(low_spec), len(up_spec)))
                                    elif len(corr_terms) <= i or len(corr_terms) <= j:
                                        new = np.zeros((len(low_spec), len(up_spec)))
                                    else:
                                        m = mu_terms[k]
                                        if isinstance(m, (int, float, np.integer, np.floating)) and m == 0:
                                            # to make it easy to zero stuff out
                                            new = np.zeros((len(low_spec), len(up_spec)))
                                        else:
                                            c_lower = corr_terms[i][low_spec, :]
                                            c_upper = corr_terms[j][up_spec, :]
                                            num = c_lower.dot(m)
                                            new = num.dot(c_upper.T)
                                        if isinstance(new, SparseArray):
                                            new = new.asarray()
                                    terms.append(
                                        new.reshape((len(low_spec), len(up_spec))
                                                    ))
                                    # raise Exception(new.toarray())
                            # print(q, a)
                            transition_moment_components[q][a] = terms

                        end = time.time()
                        logger.log_print(
                            "took {t}s",
                            t=round(end - start, 3)
                        )

            # we calculate it explicitly like this up front in case we want to use it later since the shape
            # can be a bit confusing ([0-order, 1-ord, 2-ord], [x, y, z])
            tmom = [
                sum(
                    sum(ip.chain(transition_moment_components[i][j]))
                    if not isinstance(transition_moment_components[i][j], (float, int, np.floating, np.integer))
                    else transition_moment_components[i][j]
                    for i in range(3)  # correction order
                ) for j in range(3)  # xyz
            ]

            # mu_reps already organized like x/y/z
            # mu_reps = [
            #     [m[a] for m in mu_reps]
            #     for a in range(3)
            # ]

        return [tmom, transition_moment_components, mu_reps]

    @property
    def dipole_terms(self):
        if self._dipole_terms is None:
            self._dipole_terms = DipoleTerms(self.mol).get_terms()
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
            self._tm_dat = self._transition_moments(*self._tm_dat[2])

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

    def _intensities(self, oscs):
        eng = self.energies
        units = 3.554206329390961e6
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
            ("Cubic",     (False, False, False, True))
        ))

        # wfn_terms.append(freqs.tolist())
        for key in dipole_breakdowns:
            if self.logger is not None:
                self.logger.log_print(
                    "-" * 50 + "{} {}" + "-" * 50,
                    key,
                    self.dipole_partitioning,
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
        from McUtils.Data import UnitsData

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

    def format_energies_table(self):

        from McUtils.Data import UnitsData

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        states = self.corrs.states.excitations
        harm_engs = h2w * self.zero_order_energies
        engs = h2w * self.energies

        n_modes = self.corrs.states.ndim
        harm_freq = harm_engs[1:] - harm_engs[0]
        freqs = engs[1:] - engs[0]

        return "\n".join([
                             "State Energies:",
                             ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}").format(harm_engs[0], engs[0], "-",
                                                                                         "-")
                         ] +
                         [
                             ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}").format(*s, "-", "-", e1, e2)
                             for s, e1, e2 in
                             zip(states[1:], harm_freq, freqs)
                         ]
                         )
