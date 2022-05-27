"""
Provides classes to support wave functions coming from VPT calculations
"""

import numpy as np, itertools as ip, time, enum

from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Scaffolding import ParameterManager, Checkpointer, NullCheckpointer
from McUtils.Data import UnitsData

from ..BasisReps import BasisStateSpace, SelectionRuleStateSpace, BasisMultiStateSpace, SimpleProductBasis, ExpansionWavefunctions, ExpansionWavefunction, BraKetSpace

from .Terms import DipoleTerms
from .Solver import PerturbationTheoryCorrections

__all__ = [
    # 'PerturbationTheoryWavefunction',
    'PerturbationTheoryWavefunctions'
]

__reload_hook__ = ["..Wavefun"]


# class PerturbationTheoryWavefunction(ExpansionWavefunction):
#     """
#     These things are fed the first and second order corrections
#     """
#
#     def __init__(self, mol, basis, corrections):
#         """
#         :param mol: the molecule the wavefunction is for
#         :type mol: Molecule
#         :param basis: the basis the expansion is being done in
#         :type basis: RepresentationBasis
#         :param corrections: the corrections to the terms
#         :type corrections: PerturbationTheoryCorrections
#         """
#         self.mol = mol
#         self.corrs = corrections
#         self.rep_basis = basis
#         super().__init__(self.corrs.energies, self.corrs.wavefunctions, None)
#
#     @property
#     def order(self):
#         return self.corrs.order
#
#     def expectation(self, operator, other):
#         return NotImplemented
#
#     @property
#     def zero_order_energy(self):
#         return self.corrs.energies[0]

class PerturbationTheoryWavefunctions(ExpansionWavefunctions):
    """
    These things are fed the first and second order corrections
    """

    def __init__(self, mol, basis, corrections,
                 modes=None,
                 mode_selection=None,
                 logger=None,
                 checkpoint=None,
                 results=None,
                 operator_settings=None,
                 expansion_options=None,
                 degenerate_transformation_layout=None # need to see if this is really necessary...
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

        self.checkpointer = Checkpointer.build_canonical(checkpoint)
        if results is None:
            self.results = self.checkpointer
        else:
            self.results = Checkpointer.build_canonical(results)
            
        if expansion_options is None:
            expansion_options = {}
        self.expansion_options = expansion_options
        self.operator_settings = operator_settings if operator_settings is not None else {}
        self._order = None
        self.degenerate_transformation_layout="column" if degenerate_transformation_layout is None else degenerate_transformation_layout
        # super().__init__(
        #     self.corrs.energies,
        #     self.corrs.wfn_corrections,
        #     None
        # )
    @property
    def energies(self):
        return self.corrs.energies

    def to_state(self, serializer=None):
        keys = dict(
            corrections=self.corrs,
            molecule=self.mol,
            modes=self.modes,
            mode_selection=self.mode_selection,
            basis=self.rep_basis,
            expansion_options=self.expansion_options,
            operator_settings=self.operator_settings,
            logger=self.logger
        )
        return keys

    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(
            serializer.deserialize(data['molecule']),
            serializer.deserialize(data['basis']), serializer.deserialize(data['corrections']),
            modes=serializer.deserialize(data['modes']),
            mode_selection=serializer.deserialize(data['mode_selection']),
            logger=serializer.deserialize(data['logger']),
            operator_settings=serializer.deserialize(data['operator_settings']),
            expansion_options=serializer.deserialize(data['expansion_options'])
        )

    @property
    def degenerate_transformation(self):
        return self.corrs.degenerate_transf


    def energies_to_order(self, order=None):
        if order is None:
            order = self.corrs.order
        return np.sum(self.corrs.energy_corrs[:, :order+1], axis=1)
    @property
    def deperturbed_energies(self):
        return self.energies_to_order()
    def deperturbed_frequencies(self, order=None):
        eng = self.energies_to_order(order=order)
        return eng[1:] - eng[0]

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
                                             axes=[[], []],
                                             **self.operator_settings
                                             )
    def get_M1(self, mu_1):
        return self.rep_basis.representation("x", name='M(1)', coeffs=mu_1,
                                             axes=[[0], [1]],
                                             **self.operator_settings
                                             )
    def get_M2(self, mu_2):
        return 1 / 2 * self.rep_basis.representation("x", "x", name="M(2)",
                                                     coeffs=mu_2,
                                                     axes=[[0, 1], [1, 2]],
                                                     **self.operator_settings
                                                     )
    def get_M3(self, mu_3):
        return 1 / 6 * self.rep_basis.representation("x", "x", "x", name="M(3)", coeffs=mu_3,
                                                     axes=[[0, 1, 2], [1, 2, 3]],
                                                     **self.operator_settings
                                                     )

    def get_Mi(self, i, mu):
        return 1 / np.math.factorial(i) * self.rep_basis.representation(*(("x",)*i), name="M({})".format(i), coeffs=mu,
                                                     axes=[list(range(0, i)), list(range(1, i+1))],
                                                     **self.operator_settings
                                                     )

    # @profile
    def _get_rep_states(self, m, k, order, ket_spaces, bra_spaces):
            with self.logger.block(tag="getting coupled states for {}".format(m)):
                start = time.time()
                inds = None
                rules = m.selection_rules
                # bra_space = bra_spaces.to_single().take_unique()
                # ket_space = ket_spaces.to_single().take_unique()
                # inds, filter = bra_space.get_representation_brakets(
                #     other=ket_space, return_filter=True,
                #     selection_rules=m.selection_rules
                # )

                for i in range(order):
                    ket_space = None
                    for j in range(order):
                        if i + j + k <= (order-1):
                            # self.logger.log_print("<{i}|{k}|{j}>", i=i, j=j, k=k)
                            for nk, kets in enumerate(ket_spaces.spaces):
                                if ket_space is None:
                                    ket_space = kets[j]
                                else:
                                    ket_space = ket_space.union(kets[j])

                    bra_space = None
                    for nb, bras in enumerate(bra_spaces.spaces):
                        if ket_space is None:
                            bra_space = bras[i]
                        else:
                            bra_space = ket_space.union(bras[i])

                    # self.logger.log_print("<{bras}|{k}|{kets}>", bras=bra_space, kets=ket_space, k=k)

                    if ket_space is not None:
                        filt=None
                        new_inds, filt = bra_space.get_representation_brakets(
                            other=ket_space, filter=filt, return_filter=True,
                            selection_rules=rules
                        )  # no selection rules here
                        # filters[nk][j] = filter
                        if inds is None:
                            inds = new_inds #type: BraKetSpace
                        else:
                            inds = inds.concatenate(new_inds)
                end = time.time()

                self.logger.log_print([
                    'got {sp}',
                    'took {t:.3f}s...'
                    ],
                    sp=inds,
                    t=end - start
                )

            if inds is None:
                raise ValueError("Dipole operator {} (order {}) has empty representation at order {}".format(
                    m,
                    k,
                    order-1
                ))

            inds = inds.remove_duplicates(assume_symmetric=True)

            return inds

    # from memory_profiler import profile
    # @profile
    def _mu_representations(self,
                            mu,
                            M,
                            space,
                            bra_spaces,
                            ket_spaces,
                            order,
                            partitioning,
                            rep_inds,
                            allow_higher_dipole_terms=False
                            ):

        # define out reps based on partitioning style
        if partitioning is None:
            partitioning = self.dipole_partitioning
        elif not isinstance(partitioning, self.DipolePartitioningMethod):
            partitioning = self.DipolePartitioningMethod(partitioning)

        # filters = [[None] * (order+1) for _ in range(ket_spaces.nstates)]
        def get_states(m, k):
            return self._get_rep_states(m, k, order, ket_spaces, bra_spaces)

        if all(
                    isinstance(m, (np.ndarray, SparseArray)) and m.shape == (M, M)
                    or isinstance(m, (int, float, np.integer, np.floating)) and m == 0
                    for m in mu
        ): # we were passed representations to reuse
            mu_terms = mu
            if (
                    partitioning == self.DipolePartitioningMethod.Standard and
                    order > 1
                    and len(mu_terms) < 4
            ):
                mu_3 = [d[3] for d in self.dipole_terms]
                m3 = self.get_M3(mu_3)
                if rep_inds[3] is None:
                    m3_inds = get_states(m3, 2)
                    with self.logger.block(tag="building {h}".format(h=m3)):
                        rep_inds[3] = m3_inds
                        rep_3 = m3.get_representation_matrix(rep_inds[3], self.corrs.total_basis, remove_duplicates=False)
                mu_terms.append(rep_3)
                # mu_terms = mu_terms + [rep_3]
        else:
            int_part = partitioning == self.DipolePartitioningMethod.Intuitive


            mu_coeffs = []
            # mu_0 = np.array([m[0] for m in mu])
            # mu_coeffs.append(mu_0)
            # if (not int_part and order > 0) or (int_part and order > 1):
            #     mu_1 = np.array([m[1] for m in mu])
            # else:
            #     mu_1 = None
            # mu_coeffs.append(mu_1)
            # if (not int_part and order > 1) or (int_part and order > 2):
            #     mu_2 = np.array([m[2] for m in mu])
            # else:
            #     mu_2 = None
            # mu_coeffs.append(mu_2)
            # if (not int_part and order > 2) or (int_part and order > 3):
            #     mu_3 = np.array([m[3] for m in mu])
            # else:
            #     mu_3 = None
            # mu_coeffs.append(mu_3)

            n = self.corrs.states.ndim
            for i in range(order + (1 if not int_part else 0)):
                mu_sel = []
                if isinstance(mu, (int, np.integer)) and mu == 0:
                    mu_sel = [np.zeros((n,)*i), np.zeros((n,)*i), np.zeros((n,)*i)]
                else:
                    for m in mu:
                        if isinstance(m, (int, np.integer)) and m == 0:
                            mu_sel.append(np.zeros((n,)*i))
                        else:
                            try:
                                mu_sel.append(m[i])
                            except IndexError:
                                if allow_higher_dipole_terms:
                                    mu_sel.append(np.zeros((n,)*i))
                                else:
                                    raise ValueError("dipole expansion not provided to order {}".format(i))

                mu_coeffs.append(np.array(mu_sel))

            # if (not int_part and order > 3) or (int_part and order > 4):
            #     raise NotImplementedError("dipole representations up to order {} with {} partitioning not supported yet".format(
            #         order,
            #         self.DipolePartitioningMethod
            #     ))

            reps = []
            for i, mu in enumerate(mu_coeffs):
                m = self.get_Mi(i, mu)
                if partitioning == self.DipolePartitioningMethod.Intuitive:
                    inds = get_states(m, i)
                elif partitioning == self.DipolePartitioningMethod.Standard:
                    inds = get_states(m, max(0, i-1))
                else:
                    raise ValueError("don't know how to interpret dipole partitioning {}".format(partitioning))

                if i >= len(rep_inds):
                    for d in range(i-len(rep_inds) + 1):
                        rep_inds.append(None)
                rep_inds[i] = inds
                reps.append(m)

            # reps = [m1, m2, m3]
            mu_terms = [None] * len(reps)
            tb = self.corrs.total_basis
            # with self.corrs.disk_backed():

            for i, h in enumerate(reps):
                m_pairs = rep_inds[i]
                if m_pairs is None:
                    raise ValueError("representation indices haven't filled enough to calculate {}".format(h))

                with self.logger.block(tag="building {}".format(h)):
                    start = time.time()
                    sub = h.get_representation_matrix(m_pairs, tb
                                                      , zero_element_warning=False # expect zeros in M(3)?
                                                      # , remove_duplicates=False
                                                      )
                    end = time.time()
                    self.logger.log_print('took {t:.3f}s...', t=end - start)
                # sub = generate_rep(h, m_pairs)
                mu_terms[i] = sub

        return mu_terms

    # @profile
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
            exp_order = None if 'expansion_order' not in self.expansion_options else self.expansion_options['expansion_order']
            if exp_order is None:
                exp_order = self.order
            elif isinstance(exp_order, int):
                exp_order = exp_order + 1
            else:
                if 'dipole' in exp_order:
                    exp_order = exp_order['dipole'] + 1
                elif 'default' in exp_order:
                    exp_order = exp_order['default'] + 1
                else:
                    exp_order = self.order
            order = exp_order

        logger = self.logger
        with logger.block(tag="Calculating intensities:", printoptions={'linewidth':int(1e8)}):

            space = self.corrs.coupled_states
            # store this for handling the degenerate_transformation
            lower_states_input = (0,) if lower_states is None else lower_states
            upper_states_input = tuple(range(space.nstates)) if excited_states is None else excited_states

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

                bra_spaces = lower_states  # if isinstance(lower_states, BasisStateSpace) else lower_states.to_single().take_unique()
                ket_spaces = excited_states  # if isinstance(excited_states, BasisStateSpace) else excited_states.to_single().take_unique()

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
                    corr_terms_lower = [self.corrs.wfn_corrections[i][low_spec, :] for i in range(order)]
                    logger.log_print(
                        [
                            "highest-order space: {s}",
                        ],
                        s=corr_terms_lower[-1]
                    )

                    logger.log_print('getting excited space...')
                    corr_terms_upper = [self.corrs.wfn_corrections[i][up_spec, :].T for i in range(order)]
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
            M = self.corrs.wfn_corrections[0].shape[1]
            total_space = self.corrs.total_basis
            transition_moment_components = np.zeros((order, 3)).tolist()  # x, y, and z components of the 0th, 1st, and 2nd order stuff

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

            mu_reps = self._mu_representations(
                mu,
                M,
                total_space,
                bra_spaces,
                ket_spaces,
                order,
                partitioning,
                rep_inds,
                allow_higher_dipole_terms=False if 'allow_higher_dipole_terms' not in self.expansion_options else self.expansion_options['allow_higher_dipole_terms']
            )

            # raise Exception("...")
                # logger.log_print(
                #     [
                #         "took {t:.3f}s..."
                #     ],
                #     t=end-start
                # )

            # raise Exception(mu_terms)
            if degenerate_transformation is None:
                degenerate_transformation = self.corrs.degenerate_transf
            if degenerate_transformation is not None:
                    # raise  Exception("ugh")
                    # degenerate_transformation = degenerate_transformation.T

                if len(lower_states_input) > 1:
                    deg_transf_lower = degenerate_transformation[np.ix_(lower_states_input, lower_states_input)]
                    if self.degenerate_transformation_layout == "column":
                        deg_transf_lower = deg_transf_lower.T
                else:
                    deg_transf_lower = None
                if len(upper_states_input) > 1:
                    if self.degenerate_transformation_layout == "column":
                        deg_transf_upper = degenerate_transformation[:, upper_states_input]
                    else:
                        deg_transf_upper = degenerate_transformation[upper_states_input, :].T
                else:
                    deg_transf_upper = None
            else:
                deg_transf_lower = deg_transf_upper = None

            if degenerate_transformation is not None:
                transition_moment_components_deg = np.zeros((order, 3)).tolist()  # x, y, and z components of the 0th, 1st, and 2nd order stuff
            else:
                transition_moment_components_deg = None

            corr_terms = self.corrs.wfn_corrections
            for a in range(3):
                mu_basic_terms = [r[:, :, a] for r in mu_reps]
                if partitioning == self.DipolePartitioningMethod.Intuitive:
                    mu_terms = [mu_basic_terms[0]]
                    if order > 1:
                        mu_terms.append(mu_basic_terms[1])
                    if order > 2:
                        mu_terms.append(mu_basic_terms[2])
                elif partitioning == self.DipolePartitioningMethod.Standard:
                    mu_terms = [mu_basic_terms[0] + mu_basic_terms[1]]
                    if order > 1:
                        mu_terms.append(mu_basic_terms[2])
                    if order > 2:
                        mu_terms.append(mu_basic_terms[3])
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
                            # logger.log_print("calculating corrections at order {q}...", q=q)
                            terms_depert = []
                            terms_deg = []
                            # should do this smarter
                            for i, j, k in ip.product(range(q + 1), range(q + 1), range(q + 1)):
                                if i + j + k == q:
                                    if len(mu_terms) <= k or len(corr_terms) <= i or len(corr_terms) <= j:
                                        new = np.zeros((len(low_spec), len(up_spec)))
                                        new_deg = new
                                    else:
                                        m = mu_terms[k]
                                        if isinstance(m, (int, float, np.integer, np.floating)) and m == 0:
                                            # to make it easy to zero stuff out
                                            new = np.zeros((len(lower_states_input), len(upper_states_input)))
                                            new_deg = new
                                        else:
                                            c_lower = corr_terms_lower[i]
                                            c_upper = corr_terms_upper[j]
                                            num = c_lower.dot(m)
                                            new = num.dot(c_upper)
                                            if degenerate_transformation is not None:
                                                new_deg = new
                                                if deg_transf_lower is not None:
                                                    new_deg = deg_transf_lower.dot(new_deg)
                                                if deg_transf_upper is not None:
                                                    new_deg = new_deg.dot(deg_transf_upper)
                                            else:
                                                new_deg = None
                                        if isinstance(new, SparseArray):
                                            new = new.asarray()
                                        new = new.reshape((len(lower_states_input), len(upper_states_input)))
                                    terms_depert.append(new)

                                    if degenerate_transformation is not None:
                                        if isinstance(new_deg, SparseArray):
                                            new_deg = new_deg.asarray()
                                        new_deg = new_deg.reshape((len(lower_states_input), len(upper_states_input)))
                                        terms_deg.append(new_deg)
                                    # raise Exception(new.toarray())
                            # print(q, a)
                            transition_moment_components[q][a] = terms_depert
                            if degenerate_transformation is not None:
                                transition_moment_components_deg[q][a] = terms_deg
                            end = time.time()
                            logger.log_print(
                                "took {t:.3f}s",
                                t=end - start
                            )

            # we calculate it explicitly like this up front in case we want to use it later since the shape
            # can be a bit confusing ([0-order, 1-ord, 2-ord], [x, y, z])
            tmom = self._compute_tmom_to_order(transition_moment_components, order-1)

            if degenerate_transformation is not None:
                tmom_deg = self._compute_tmom_to_order(transition_moment_components_deg, order-1)
            else:
                tmom_deg = None

            # mu_reps already organized like x/y/z
            # mu_reps = [
            #     [m[a] for m in mu_reps]
            #     for a in range(3)
            # ]

        # print(
        #     (tmom, transition_moment_components),
        #     (tmom_deg, transition_moment_components_deg)
        #
        # )
        with self.checkpointer:
            if degenerate_transformation is not None:
                self.checkpointer["transition_moments"] = transition_moment_components_deg
                self.checkpointer["nondegenerate_transition_moments"] = transition_moment_components
            else:
                self.checkpointer["transition_moments"] = transition_moment_components

        return [(tmom, transition_moment_components), (tmom_deg, transition_moment_components_deg), mu_reps, (corr_terms_lower, corr_terms_upper)]

    @classmethod
    def _compute_tmom_to_order(cls, transition_moment_components, order):
        # we calculate it explicitly like this up front in case we want to use it later since the shape
        # can be a bit confusing ([0-order, 1-ord, 2-ord], [x, y, z])
        return [
            sum(
                sum(ip.chain(transition_moment_components[i][j]))
                if not isinstance(transition_moment_components[i][j], (float, int, np.floating, np.integer))
                else transition_moment_components[i][j]
                for i in range(order+1)  # correction order
            ) for j in range(3)  # xyz
        ]

    class TermHolder(tuple):
        """symbolic wrapper on Tuple so we can know that we've canonicalized some term"""
    @property
    def dipole_terms(self):
        if self._dipole_terms is None:
            if 'dipole_terms' in self.expansion_options:
                self._dipole_terms = self.TermHolder(self.expansion_options['dipole_terms'])
            else:
                self._dipole_terms = self.TermHolder(
                    DipoleTerms(self.mol, modes=self.modes, mode_selection=self.mode_selection,
                                **ParameterManager(self.expansion_options).filter(DipoleTerms)
                                ).get_terms()
                )
        elif not isinstance(self._dipole_terms, self.TermHolder):
            self._dipole_terms = self.TermHolder(
                DipoleTerms(self.mol, modes=self.modes, mode_selection=self.mode_selection,
                            dipole_derivatives=self._dipole_terms,
                            **ParameterManager(self.expansion_options).filter(DipoleTerms)
                            ).get_terms()
            )
        # print([x for x in self._dipole_terms])
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
        if self._tm_dat[1][0] is not None:
            return self._tm_dat[1][0]
        else:
            return self._tm_dat[0][0]

    @property
    def deperturbed_transition_moments(self):
        """
        Computes the transition moments between wavefunctions stored in the object

        :return:
        :rtype:
        """

        self._load_tm_dat()
        return self._tm_dat[0][0]

    @property
    def transition_moment_corrections(self):
        """
        Computes the transition moment corrections between wavefunctions stored in the object

        :return:
        :rtype:
        """

        self._load_tm_dat()
        if self._tm_dat[1][0] is not None:
            return self._tm_dat[1][1]
        else:
            return self._tm_dat[0][1]

    @property
    def deperturbed_transition_moment_corrections(self):
        """
        Computes the transition moment corrections between wavefunctions stored in the object

        :return:
        :rtype:
        """

        self._load_tm_dat()
        return self._tm_dat[0][1]

    @property
    def oscillator_strengths(self):
        """
        Computes the oscillator strengths for transitions from the ground state to the other states

        :return:
        :rtype:
        """

        tms = self.transition_moments
        return self._oscillator_strengths(tms)

    @property
    def deperturbed_oscillator_strengths(self):
        """
        Computes the oscillator strengths for transitions from the ground state to the other states

        :return:
        :rtype:
        """

        tms = self.deperturbed_transition_moments
        return self._oscillator_strengths(tms)

    def oscillator_strengths_to_order(self, order):
        """

        :param tms:
        :type tms:
        :return:
        :rtype:
        """

        return self._oscillator_strengths(self._compute_tmom_to_order(
            self.transition_moment_corrections,
            order
        ))

    def deperturbed_oscillator_strengths_to_order(self, order):
        """

        :param tms:
        :type tms:
        :return:
        :rtype:
        """

        return self._oscillator_strengths(self._compute_tmom_to_order(
            self.deperturbed_transition_moment_corrections,
            order
        ))

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

        freqs, ints = self._intensities(self.oscillator_strengths)
        with self.results: # TODO: this really isn't the right place for this...
            self.results['spectrum'] = (
                freqs,
                ints
            )

        return ints

    @property
    def deperturbed_intensities(self):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        return self._dep_intensities(self.deperturbed_oscillator_strengths)[1]

    def intensities_to_order(self, order):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        return self._intensities(self.oscillator_strengths_to_order(order), energy_order=order)[1]

    def deperturbed_intensities_to_order(self, order):
        """
        Computes the intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        return self._dep_intensities(self.deperturbed_oscillator_strengths_to_order(order), energy_order=order)[1]

    def _dep_intensities(self, oscs, energy_order=None):
        if energy_order is None:
            eng = self.deperturbed_energies
        else:
            eng = self.energies_to_order(energy_order)
        freqs = eng - eng[0]
        units = UnitsData.convert("OscillatorStrength", "KilometersPerMole")
        return freqs, units * freqs * oscs

    def _intensities(self, oscs, energy_order=None):
        if energy_order is None:
            eng = self.energies
        elif self.degenerate_transformation is not None:
            raise NotImplementedError("ugh")
        else:
            eng = self.energies_to_order(energy_order)
        freqs = eng - eng[0]
        units = UnitsData.convert("OscillatorStrength", "KilometersPerMole")
        return freqs, units * freqs * oscs

    @property
    def zero_order_intensities(self):
        """
        Computes the harmonic intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        eng = self.zero_order_energies
        ints = self.deperturbed_intensities_to_order(0)

        with self.results: #TODO: this is the wrong place for this
            self.results['zero_order_spectrum'] = (
                eng - eng[0],
                ints
            )

        return ints

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
                freqs, ints = self._intensities(oscs)

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
    def _format_energies_table(cls, states, zpe, freqs, real_fmt='{:>12.5f}'):

        n_modes = len(states[0])

        zpe = np.asanyarray(zpe)
        freqs = np.asanyarray(freqs)

        ncols = freqs.shape[1] if len(freqs.shape) > 0 else 1
        states_fmt = '{:<1.0f} ' * n_modes
        padding = len(states_fmt.format(*[0 for _ in range(n_modes)]))
        num_digits = len(real_fmt.format(0))
        dash_fmt = '{:>' + str(num_digits) + "}"
        htag = "Harmonic"; hpad = (num_digits - len(htag))
        atag = "Anharmonic"; apad = (num_digits - len(atag))
        bar = (
                (" " * hpad if hpad > 0 else "") + htag
                + " "
                + (" " * apad if apad > 0 else "") + atag
                )
        if zpe is not None:
            header = ["State" + " " * (padding - len("State")) + bar + " " + bar]
            header += [
                " " * padding
                + " " * (num_digits - 3) + "ZPE"
                + " " + " " * (num_digits - 3) + "ZPE"
                + " "
                + " " * (num_digits - 9) + "Frequency"
                + " " + " " * (num_digits - 9) + "Frequency"
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
                + " " * (num_digits - 9) + "Frequency"
                + " " + " " * (num_digits - 9) + "Frequency"
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


    def format_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'):
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
            real_fmt=real_fmt
        )

    def format_deperturbed_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'):
        # simple utility function pulled from the unit tests

        if states is None:
            states = self.corrs.states.excitations

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        harm_engs = self.zero_order_energies * h2w
        engs = self.deperturbed_energies * h2w
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
            real_fmt=real_fmt
        )

    def format_property_matrices(self, states, prop_corrs, real_fmt="{:>.8e}", padding_fmt='{:>16}'):

        # reshape data so it's one big matrix...
        order = len(prop_corrs)
        term_keys = []
        terms = []
        for q in range(order):  # total quanta
            for pq in prop_corrs[q]:
                if len(pq) > 1:
                    raise NotImplementedError(
                        "don't know what to do when more than one starting state (got {})".format(
                            len(pq)
                        )
                    )
                terms.append(pq[0])
            for i, j, k in ip.product(range(q + 1), range(q + 1), range(q + 1)):
                if i + j + k == q:
                    term_keys.append("<{i}|M({k})|{j}>".format(i=i, j=j, k=k))

        # print(len(terms), [len(x) for x in terms])
        terms = np.array(terms).T

        nels = len(padding_fmt.format(real_fmt.format(0)))
        header_fmt = "{:<" + str(nels) + "}"
        header = " ".join(header_fmt.format(k) for k in term_keys)

        n_modes = self.corrs.total_basis.ndim
        padding = np.max([len(str("0 " * n_modes)), 1]) + 1
        header = "State" + " " * (padding - 3) + header

        n_terms = len(term_keys)
        prop_mat = header + "\n" + "\n".join(
            "  " + ("{:<1.0f} " * n_modes).format(*s) +
            " ".join(padding_fmt.format(real_fmt.format(ll)) for ll in l)
            for s, l in zip(states, terms)
        )

        return prop_mat

    def format_dipole_contribs_tables(self):

        states = self.corrs.states.excitations
        corrs = self.transition_moment_corrections
        mats = []
        for a in range(3):
            prop = [x[a] for x in corrs]
            mats.append(self.format_property_matrices(states, prop))

        return mats

    def format_deperturbed_dipole_contribs_tables(self):

        states = self.corrs.states.excitations
        corrs = self.deperturbed_transition_moment_corrections
        mats = []
        for a in range(3):
            prop = [x[a] for x in corrs]
            mats.append(self.format_property_matrices(states, prop))

        return mats

    def format_energy_corrections_table(self, real_fmt="{:>12.5f}"):

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        states = self.corrs.states.excitations
        e_corrs = self.corrs.all_energy_corrs #ordered like (state, order, corr_list)
        order = len(e_corrs[0])
        term_keys = []
        for k in range(1, order):
            if e_corrs[0][k] is not None: # from odd order terms...
                for i in range(0, k):
                    term_keys.append("<0|dH({j})|{i}>".format(j=k-i, i=i))
        terms = h2w * np.array([np.concatenate([y for y in x if y is not None]) for x in e_corrs])

        nels = len(real_fmt.format(0))
        header_fmt = "{:<" + str(nels) + "}"
        header = " ".join(header_fmt.format(k) for k in term_keys)

        n_modes = self.corrs.total_basis.ndim
        padding = np.max([len(str("0 " * n_modes)), 1]) + 1
        header = "State" + " " * (padding - 3) + header

        prop_mat = header + "\n" + "\n".join(
            "  " + ("{:<1.0f} " * n_modes).format(*s) +
            " ".join(real_fmt.format(ll) for ll in l)
            for s, l in zip(states, terms)
        )

        return prop_mat

    def _format_intensities_table(self, states, freqs, ints, harm_freqs, harm_ints, real_fmt="{:>12.5f}"):
        # simple utility function pulled from the unit tests

        n_modes = self.corrs.total_basis.ndim
        padding = np.max([len(str("0 " * n_modes)), 1]) + 1

        num_digits = len(real_fmt.format(0))
        ftag = "Frequency"; fpad = (num_digits - len(ftag))
        itag = "Intensity"; ipad = (num_digits - len(itag))
        bar = (
            (" " * fpad if fpad > 0 else "") + ftag
            + " "
            + (" " * ipad if ipad > 0 else "") + itag
        )

        spacer = "    "
        report = (
                    " " * (padding + 2) + " " * (num_digits - 2) + "Harmonic" +
                    " " * (num_digits - len("Harmonic")) + spacer
                    + " " * (num_digits - 2) + "Anharmonic\n"
                    + "State" + " " * (padding - 3) + bar + spacer + bar + "\n"
                 ) + "\n".join(
            (
                (
                        "  " + ("{:<1.0f} " * n_modes) + " "
                        + real_fmt + " " + real_fmt
                        + spacer
                        + real_fmt + " " + real_fmt
                ).format(*s, hf, hi, f, i)
                for s, hf, hi, f, i in zip(states[1:], harm_freqs[1:], harm_ints[1:], freqs[1:], ints[1:])
            ))

        return report

    def format_intensities_table(self, real_fmt="{:>12.5f}"):
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        states = self.corrs.states.excitations

        engs = h2w * self.energies
        freqs = engs - engs[0]
        ints = self.intensities

        harm_engs = h2w * self.zero_order_energies
        harm_freqs = harm_engs - harm_engs[0]
        harm_ints = self.zero_order_intensities

        return self._format_intensities_table(states, freqs, ints, harm_freqs, harm_ints)

    def format_deperturbed_intensities_table(self, real_fmt="{:>12.5f}"):
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        states = self.corrs.states.excitations

        engs = h2w * self.deperturbed_energies
        freqs = engs - engs[0]
        ints = self.deperturbed_intensities

        harm_engs = h2w * self.zero_order_energies
        harm_freqs = harm_engs - harm_engs[0]
        harm_ints = self.zero_order_intensities

        return self._format_intensities_table(states, freqs, ints, harm_freqs, harm_ints, real_fmt=real_fmt)