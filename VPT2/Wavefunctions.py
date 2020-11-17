"""
Provides classes to support wave functions coming from VPT calculations
"""

import numpy as np, itertools as ip

from ..BasisReps import SimpleProductBasis, ExpansionWavefunctions, ExpansionWavefunction

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

    def __init__(self, mol, basis, corrections):
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
        self.dipole_partitioning = 'standard'
        super().__init__(
            self.corrs.energies,
            self.corrs.wfn_corrections,
            None
        )

    @property
    def order(self):
        return self.corrs.order

    def expectation(self, operator, other):
        return NotImplemented

    @property
    def zero_order_energies(self):
        return self.corrs.energy_corrs[:, 0]

    def _transition_moments(self, mu_x, mu_y, mu_z, partitioning=None):
        """
        Calculates the x, y, and z components of the
        transition moment between the wavefunctions stored

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

        M = self.corrs.coupled_states  # the coupled subspace space we're working in
        corr_vecs = self.corrs.wfn_corrections#[..., M]
        transition_moment_components = np.zeros((3, 3)).tolist()  # x, y, and z components of the 0th, 1st, and 2nd order stuff
        # raise Exception([mu_1[0].shape)
        mu = [mu_x, mu_y, mu_z]
        for a in range(3):  # x, y, and z
            mu_0, mu_1, mu_2, mu_3 = mu[a]
            if isinstance(mu[a][0], np.ndarray) and mu[a][0].shape == (M, M):
                mu_terms = mu
            else:
                if partitioning is None:
                    partitioning = self.dipole_partitioning
                if partitioning == 'intuitive':
                    m1 = self.rep_basis.representation(coeffs=mu_0)
                    m1 = m1[np.ix_(M, M)]
                    m2 = self.rep_basis.representation("x", coeffs=mu_1)
                    m2 = m2[np.ix_(M, M)]
                    m3 = 1 / 2 * self.rep_basis.representation("x", "x", coeffs=mu_2)
                    m3 = m3[np.ix_(M, M)]
                else:
                    m1 = self.rep_basis.representation("x", coeffs=mu_1) + self.rep_basis.representation(coeffs=mu_0)
                    m1 = m1[np.ix_(M, M)]
                    m2 = 1 / 2 * self.rep_basis.representation("x", "x", coeffs=mu_2)
                    m2 = m2[np.ix_(M, M)]
                    m3 = 1 / 6 * self.rep_basis.representation("x", "x", "x", coeffs=mu_3)
                    m3 = m3[np.ix_(M, M)]

                mu_terms = [m1, m2, m3]

            mu_terms = [m.todense() if not isinstance(m, (int, float, np.integer, np.floating, np.ndarray)) else m for m
                        in mu_terms]

            corr_terms = [corr_vecs[:, 0], corr_vecs[:, 1], corr_vecs[:, 2]]

            for q in range(len(mu_terms)):  # total quanta
                # print([(i, j, k) for i,j,k in ip.product(range(q+1), range(q+1), range(q+1)) if i+j+k==q])
                transition_moment_components[q][a] = [
                    np.dot(np.dot(corr_terms[i], mu_terms[k]), corr_terms[j].T)
                    for i, j, k in ip.product(range(q + 1), range(q + 1), range(q + 1)) if i + j + k == q
                ]

        # we calculate it explicitly like this up front in case we want to use it later since the shape
        # can be a bit confusing ([0-order, 1-ord, 2-ord], [x, y, z])
        tmom = [
            sum(
                sum(ip.chain(transition_moment_components[i][j]))
                for i in range(3)  # correction order
            ) for j in range(3)  # xyz
        ]
        # raise Exception(tmom)
        return tmom, transition_moment_components

    @property
    def dipole_terms(self):
        if self._dipole_terms is None:
            self._dipole_terms = DipoleTerms(self.mol).get_terms()
        return self._dipole_terms

    @dipole_terms.setter
    def dipole_terms(self, v):
        self._dipole_terms = v
        self._tm_dat = None

    @property
    def transition_moments(self):
        """
        Computes the transition moments between wavefunctions stored in the object

        :return:
        :rtype:
        """

        if self._tm_dat is None:
            self._tm_dat = self._transition_moments(*self.dipole_terms)
        return self._tm_dat[0]

    @property
    def transition_moment_corrections(self):
        """
        Computes the transition moment corrections between wavefunctions stored in the object

        :return:
        :rtype:
        """

        if self._tm_dat is None:
            self._tm_dat = self._transition_moments(*self.dipole_terms)
        return self._tm_dat[1]

    @property
    def oscillator_strengths(self):
        """
        Computes the oscillator strengths for transitions from the ground state to the other states

        :return:
        :rtype:
        """

        tms = self.transition_moments
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
        eng = self.energies
        units = 3.554206329390961e6
        return units * (eng - eng[0]) * self.oscillator_strengths

    @property
    def zero_order_intensities(self):
        """
        Computes the harmonic intensities for transitions from the ground state to the other states

        :return:
        :rtype:
        """
        eng = self.zero_order_energies
        tm = np.array([self.transition_moment_corrections[0][x][0][0] for x in range(3)]).T
        osc = np.linalg.norm(tm, axis=1) ** 2
        units = 3.554206329390961e6
        return units * (eng - eng[0]) * osc

    def generate_intensity_breakdown(self):
        """
        Generates a breakdown of the terms that contribute to the intensity
        Returns in a format that can be directly exported to JSON if desired.

        :return:
        :rtype:
        """
        from collections import OrderedDict

        engs = self.energies
        freqs = engs - engs[0]

        # harm_engs = self.zero_order_energies
        # harm_freqs = harm_engs - harm_engs[0]

        terms = OrderedDict((
            ("frequencies", freqs.tolist()),
            ('axes', self.mol.inertial_axes.tolist()),
            ("states", self.rep_basis.unravel_state_inds(self.corrs.states)),
            ('breakdowns', OrderedDict()),
            ("wavefunctions", {
                "corrections": self.corrs.wfn_corrections[:, :, self.corrs.coupled_states].tolist(),
                "coupled_states": self.rep_basis.unravel_state_inds(self.corrs.coupled_states)
            })
        ))
        bds = terms['breakdowns']

        dts = self.dipole_terms

        dipole_breakdowns = OrderedDict((
            ("Linear", [[d[0], d[1], np.zeros(d[2].shape), np.zeros(d[3].shape)] for d in dts]),
            ("Quadratic", [[d[0], d[1], d[2], np.zeros(d[3].shape)] for d in dts]),
            ("Cubic", [[d[0], d[1], np.zeros(d[2].shape), d[3]] for d in dts]),
            ("Full", dts)
        ))

        # wfn_terms.append(freqs.tolist())
        for key in dipole_breakdowns:
            dt = dipole_breakdowns[key]
            self.dipole_terms = dt
            ints = self.intensities

            # for q in range(len(mu_terms)): # total quanta
            #     # print([(i, j, k) for i,j,k in ip.product(range(q+1), range(q+1), range(q+1)) if i+j+k==q])
            #     transition_moment_components[q][a] = [
            #         np.dot(np.dot(corr_terms[i], mu_terms[k]), corr_terms[j].T)
            #         for i, j, k in ip.product(range(q+1), range(q+1), range(q+1)) if i+j+k==q
            #     ]

            full_corrs = self.transition_moment_corrections
            corrs_keys = [y for x in ip.chain(
                [(i, j, k) for i, j, k in ip.product(range(q + 1), range(q + 1), range(q + 1)) if i + j + k == q]
                for q in range(len(full_corrs))
            ) for y in x]

            full_corrs = np.concatenate([np.array(x) for x in full_corrs], axis=1)
            # transpose so we have x-y-z as outermost dimension and correction as innermost
            full_corrs = full_corrs.transpose((0, 2, 3, 1)).tolist()

            # harm_ints = self.zero_order_intensities
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
                corr_keys = ["<gs^({})|mu^({})|es^({})>".format(*lab) for lab in corrs['keys']]
                for ax, vals in corrs['values'].items():
                    corr_block.append([ax])
                    corr_block.append([" "] + corr_keys)
                    for s, v in zip(states, vals[0]):
                        corr_block.append([s] + list(v))
                writer.writerows([padding + x for x in corr_block])

            coupled_state_blocks = [["Wavefunction Corrections"]]
            wfn_corrs = np.array(ib['wavefunctions']["corrections"])
            coupled_states = ib['wavefunctions']["coupled_states"]
            num_corrs = wfn_corrs.shape[1]
            coupled_state_blocks.append(sum([ [s] + ["|n^({})>".format(i) for i in range(num_corrs)] for s in states ], []))
            for corr_block, state in zip(wfn_corrs.transpose((2, 0, 1)), coupled_states):
                row = []
                for corr in corr_block:
                    row.extend([state] + list(corr))
                coupled_state_blocks.append(row)

            writer.writerows([padding + x for x in coupled_state_blocks])

            # ints = ib[ey]
            # for

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
