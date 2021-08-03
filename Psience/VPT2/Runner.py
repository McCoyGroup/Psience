"""
A little package of utilities for setting up/running VPT jobs
"""

import numpy as np

from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis


class VPTSystem:
    def __init__(self, src):
        ...



class VPTRunner:
    """
    A helper class to make it easier to run jobs
    """

    def __init__(self,
                 system,
                 modes=None,
                 mode_selection=None,

                 ):

        self.system = system
        ...

    @staticmethod
    def get_states(n_quanta, n_modes, target_modes=None):
        whee = [np.flip(x) for x in BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(n_modes),
            range(n_quanta)
        ).excitations]
        if target_modes is not None:
            whee = [
                p for p in whee if sum(p) == 0 or any(p[i] > 0 for i in target_modes)
            ]
        return whee

    @staticmethod
    def get_degenerate_polyad_space(states, polyadic_pairs):
        # we build a graph of connected states by the polyadic rules
        polyadic_pairs = np.array(polyadic_pairs)
        states = np.array(states)
        poss_degs = [[] for _ in states]
        check_list = states.tolist()
        for n, s in enumerate(states):  # build list-of-lists structure
            for i, nt_spec in enumerate(polyadic_pairs):
                if np.all(s - nt_spec[0] >= 0):
                    new = (s - nt_spec[0] + nt_spec[1]).tolist()
                    if new not in check_list:
                        check_list.append(new)
                        poss_degs.append([])
                    # else:
                    #     poss_degs[idx].append(slist)
                    poss_degs[n].append(new)

        # from the populated lists build the real connection graph
        groups = [[] for _ in check_list]
        new_checks = []
        for i, s1 in enumerate(check_list):
            if s1 not in new_checks:
                new_checks.append(s1)
                groups[i].append(s1)
                groups[i].extend(poss_degs[i])
                for s2, p in zip(check_list[i + 1:], poss_degs[i + 1:]):
                    if s2 not in new_checks:
                        if s2 in poss_degs[i]:
                            new_checks.append(s2)
                            groups[i].extend(p)

        return [g for g in groups if len(g) > 1]

    def get_Hamiltonian(self):
        ...

    def get_wavefunctions(self):
        ...

    def check_potential_derivatives(self):
        ...

