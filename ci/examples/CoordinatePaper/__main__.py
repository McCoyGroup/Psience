import os.path

import numpy as np
from Psience.Molecools import Molecule
from Psience.Modes import NormalModes, LocalizedModes
from Psience.VPT2 import VPTAnalyzer
from McUtils.Formatters import TableFormatter
from McUtils.Profilers import Timer
from . import expansions
from . import modes
from . import vpt

tests = []

def inactive(test):
    ...
def runnable(test):
    tests.append(test)
    return test

@inactive
def Duschinksy():

    mol1, (locs1, nms1, rpnms1) = modes.get_gaussian_modes('b3lyp', 'methanol_vpt_3',
                                                           use_internals=True,
                                                           return_mol=True, return_status=False)
    mol2, (locs2, nms2, rpnms2) = modes.get_gaussian_modes('b3lyp', 'methanol_vpt_5',
                                                           use_internals=True,
                                                           return_mol=True, return_status=False)

    locs1: LocalizedModes
    locs2: LocalizedModes

    print("="*50)
    d_cart = (locs1.coords_by_modes @ locs2.modes_by_coords)
    print(np.round(np.linalg.norm(d_cart, axis=0), 3))
    print(
        TableFormatter('{:.3f}').format(
            d_cart**2
        )
    )

    tf1 = mol1.get_internals_by_cartesians(1)[0]
    inv1 = mol1.get_cartesians_by_internals(1)[0]
    tf2 = mol2.get_internals_by_cartesians(1)[0]
    inv2 = mol2.get_cartesians_by_internals(1)[0]

    print("="*50)
    d_int = (locs1.coords_by_modes @ tf1 @ inv2 @locs2.modes_by_coords)
    print(np.round(np.linalg.norm(d_int, axis=0), 3))
    print(
        TableFormatter('{:.3f}').format(
            d_int**2
        )
    )

@inactive
def AIMNetVPT():
    for k in [0, -30]:
        print(k, "-", "Full")
        vpt.run_aimnet_vpt(k)
        print(k, "-", "Full Internals")
        vpt.run_aimnet_vpt(k, use_internals=True)
        print(k, "-", "Partial")
        vpt.run_aimnet_vpt(k, mode_selection=[-4, -3, -2, -1])
        print(k, "-", "Partial Internals")
        vpt.run_aimnet_vpt(k, mode_selection=[-4, -3, -2, -1], use_internals=True)
        print(k, "-", "Partial")
        vpt.run_aimnet_vpt(k, mode_selection=[-7, -6, -5, -4, -3, -2, -1])
        print(k, "-", "Partial Internals")
        vpt.run_aimnet_vpt(k, mode_selection=[-7, -6, -5, -4, -3, -2, -1], use_internals=True)


@inactive
def B3LYPVPT():
    for k in [1, 3]:
        print(k, "-", "Full")
        vpt.run_gaussian_vpt('b3lyp', k)
        print(k, "-", "Full Internals")
        vpt.run_gaussian_vpt('b3lyp', k, use_internals=True)
        print(k, "-", "Partial")
        vpt.run_gaussian_vpt('b3lyp', k, mode_selection=[-4, -3, -2, -1])
        print(k, "-", "Partial Internals")
        vpt.run_gaussian_vpt('b3lyp', k, mode_selection=[-4, -3, -2, -1], use_internals=True)
        print(k, "-", "Partial")
        vpt.run_gaussian_vpt('b3lyp', k, mode_selection=[-7, -6, -5, -4, -3, -2, -1])
        print(k, "-", "Partial Internals")
        vpt.run_gaussian_vpt('b3lyp', k, mode_selection=[-7, -6, -5, -4, -3, -2, -1], use_internals=True)

@inactive
def wb97VPT():
    for k in [1, 3]:
        print(k, "-", "Full")
        vpt.run_wb97_vpt(k)
        print(k, "-", "Full Internals")
        vpt.run_wb97_vpt(k, use_internals=True)
        print(k, "-", "Partial")
        vpt.run_wb97_vpt(k, mode_selection=[-4, -3, -2, -1])
        print(k, "-", "Partial Internals")
        vpt.run_wb97_vpt(k, mode_selection=[-4, -3, -2, -1], use_internals=True)
        print(k, "-", "Partial")
        vpt.run_wb97_vpt(k, mode_selection=[-7, -6, -5, -4, -3, -2, -1])
        print(k, "-", "Partial Internals")
        vpt.run_wb97_vpt(k, mode_selection=[-7, -6, -5, -4, -3, -2, -1], use_internals=True)

@runnable
def NORPVPT():
        # vpt.run_aimnet_vpt(0, use_internals=True, use_reaction_path=False)
        # vpt.run_wb97_vpt(1, use_internals=True, use_reaction_path=False)
        vpt.run_b3lyp_vpt(1, use_internals=False, use_reaction_path=False)
        vpt.run_b3lyp_vpt(1, use_internals=True, use_reaction_path=False)

@runnable
def analyzeStretchResults():
    def fprint(a):
        a = np.asanyarray(a)
        if a.ndim == 1: a = a[np.newaxis]
        print(TableFormatter("{:>4.0f}").format(a))

    def print_freq_comps_info(tag, log_getter, key, **opts):
        print("="*20, tag, "="*20)
        base = log_getter(key, use_reaction_path=False, **opts)
        if os.path.exists(base):
            a0 = VPTAnalyzer(base)
        else:
            a0 = None
        a1 = VPTAnalyzer(log_getter(key, **opts))
        a2 = VPTAnalyzer(log_getter(key, mode_selection=[-7, -6, -5, -4, -3, -2, -1], **opts))
        if a0 is not None:
            e0 = a0.deperturbed_energies[1][1:8]
            fprint(e0)
        else:
            e0 = None
        e1 = a1.deperturbed_energies[1][1:8]
        e2 = a2.deperturbed_energies[1][1:8]
        fprint(e1)
        fprint(e2)
        if e0 is not None:
            fprint(e0 - e2)
        fprint(e0 - e1)
        fprint(e1 - e2)

    # print_freq_comps_info(
    #     "AIMNet Cartesian",
    #     vpt.get_aimnet_logfile,
    #     0
    # )
    #
    # print_freq_comps_info(
    #     "AIMNet Internals",
    #     vpt.get_aimnet_logfile,
    #     0,
    #     use_internals=True
    # )
    #
    # print_freq_comps_info(
    #     "WB97 Cartesians",
    #     vpt.get_wb97_logfile,
    #     1
    # )
    #
    # print_freq_comps_info(
    #     "WB97 Internals",
    #     vpt.get_wb97_logfile,
    #     1,
    #     use_internals=True
    # )

    print_freq_comps_info(
        "B3LYP Cartesians",
        vpt.get_b3lyp_logfile,
        1
    )

    print_freq_comps_info(
        "B3LYP Internals",
        vpt.get_b3lyp_logfile,
        1,
        use_internals=True
    )


if __name__ == '__main__':
    print(expansions.methanol_gaussian_zmatrix)

    # for t in tests:
    #     with Timer(t.__name__, file=None):
    #         print("="*50)
    #         print(t.__name__, ":")
    #         t()
