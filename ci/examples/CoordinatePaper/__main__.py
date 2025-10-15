import os.path

import numpy as np
from Psience.Molecools import Molecule
from Psience.Modes import NormalModes, LocalizedModes
from Psience.VPT2 import VPTAnalyzer
from McUtils.Formatters import TableFormatter
from McUtils.Data import UnitsData
import McUtils.Formatters as mfmt
import McUtils.Numputils as nput
import McUtils.Coordinerds as coordops
import McUtils.Plots as plt
from McUtils.Profilers import Timer
from . import expansions
from . import modes
from . import vpt
from . import util
from . import analysis

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
def setup_mace_expansions():
    # from ase.optimize import BFGS
    # import scipy.optimize._optimize
    # new = expansions.methanol_mace
    # print(new.calculate_energy(order=1)[1]
    #       * UnitsData.convert("Hartrees", "ElectronVolts")
    #       / UnitsData.convert("BohrRadius", "Angstroms")
    #       )
    # new = expansions.methanol_mace.optimize(
    #     max_iterations=15,
    #     max_displacement=1e-12,
    #     tol=1e-12,
    #     # logger=True,
    #     track_best=True,
    #     # method='gradient-descent'
    #     # method='conjugate-gradient',
    #     # method='quasi-newton',
    #     optimizer_settings={
    #         # 'restart_parameter': None
    #     },
    #     # line_search=True
    #     # restart_interval=100
    #     # mode='scipy',
    #     # method='cg'
    #     # method='ase'
    # )
    # print(np.linalg.norm(new.calculate_energy(order=1)[1]) * UnitsData.hartrees_to_wavenumbers)
    # print(new.coords.tolist())
    # # new.plot().show()
    # return
    for i in [-50]:#range(0, -70, -10):
        with Timer(tag=f"{i}"):
            base_struct = expansions.get_mace_structure(i, 'optimized', overwrite=True)
            harmonic_expansion = expansions.get_mace_expansion(i, 'optimized', return_harmonic=True, overwrite=True)
            mol = expansions.methanol_mace.modify(coords=base_struct)
            rpnms, stat = mol.get_reaction_path_modes(
                potential_derivatives=harmonic_expansion[1:],
                return_status=True,
                zero_gradient_cutoff=25 / UnitsData.hartrees_to_wavenumbers
            )
            print(harmonic_expansion[1])
            print(f"MACE-OFF ({i}):", stat, rpnms.freqs * UnitsData.hartrees_to_wavenumbers)
            # return
            # structs = expansions.get_aimnet_expansion(i, 'optimized',
            #                                           chk_pattern='checkpoints/test_aimnet_{key}.json'
            # )
            # print(f"AIMNet2 ({i}):", structs[0][0].freqs * UnitsData.hartrees_to_wavenumbers)
            # structs = expansions.get_mace_expansion(i, 'optimized', overwrite=True)
            # print(f"MACE-OFF ({i}):", structs[0][0].freqs * UnitsData.hartrees_to_wavenumbers)

@inactive
def analytic_expansions():
    woof = expansions.get_aimnet_expansion(0, 'optimized', analytic_derivative_order=3)
    print(np.array(woof[1][4]).shape)

@runnable
def AIMNetVPT():
    degs = True
    # vpt.run_aimnet_vpt(0, use_degeneracies=False, use_internals=False, use_reaction_path=False)
    # vpt.run_aimnet_vpt(0, use_degeneracies=False, use_internals=True, use_reaction_path=False)
    step_size = None#.25
    import warnings
    warnings.filterwarnings('error')
    # overwrite = True
    # for k in range(0, -70, -10):
    #     print(f"OPTIMIZING: {k}")
    #     # meth, meth_int = expansions.get_aimnet_methanol()
    #     base_struct = expansions.get_aimnet_structure(k, 'optimized', overwrite=overwrite)
    #     harmonic_expansion = expansions.get_aimnet_expansion(k, 'optimized', overwrite=overwrite)
    overwrite = False
    for k in range(0, -70, -10):
        # overwrite = k < -15
        # meth, meth_int = expansions.get_aimnet_methanol()
        # base_struct = expansions.get_aimnet_structure(k, 'optimized', overwrite=overwrite)
        # harmonic_expansion = expansions.get_aimnet_expansion(k, 'optimized', overwrite=overwrite)
        # wtf1a = nput.tensor_reexpand(
        #     meth_int.modify(coords=base_struct).get_cartesians_by_internals(2),
        #     harmonic_expansion[1:]
        # )
        # print(TableFormatter("{:.1f}").format([wtf1a[0] * UnitsData.hartrees_to_wavenumbers]))
        # print(TableFormatter("{:.1f}").format([np.asanyarray(harmonic_expansion[1]) * UnitsData.hartrees_to_wavenumbers]))
        # # print(harmonic_expansion[1])
        # continue
        # print(harmonic_expansion[0])
        # print(harmonic_expansion[1] * 219474.63)
        # mol = expansions.methanol.modify(energy_evaluator='aimnet2', coords=base_struct)
        # print(coordops.zmatrix_unit_convert(
        #     mol.modify(internals=expansions.methanol_zmatrix).internal_coordinates,
        #     UnitsData.bohr_to_angstroms,
        #     rad2deg=True
        # ))
        # rpnms, stat = mol.get_reaction_path_modes(
        #     potential_derivatives=harmonic_expansion[1:],
        #     return_status=True,
        #     zero_gradient_cutoff=25 / UnitsData.hartrees_to_wavenumbers
        # )
        # print(f"AIMNet2 ({k}):", stat, rpnms.freqs * UnitsData.hartrees_to_wavenumbers)
        # # print(f"MACE-OFF ({i}):", stat, rpnms.freqs * UnitsData.hartrees_to_wavenumbers)
        # return
        print(k, "-", "Full")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs, step_size=step_size, overwrite=overwrite)
        print(k, "-", "Full Internals")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs, use_internals=True, step_size=step_size, overwrite=overwrite)
        for spec in [
            [-4, -3, -2, -1],
            [-7, -6, -5, -4, -3, -2, -1],
            [-11, -10, -8, -1],
            [-8, -1],
            [-11, -1],
            [-11, -4, -3, -2, -1],
            [-11, -7, -6, -5, -4, -3, -2, -1],
            [-1],
        ]:
            print(k, "-", "Partial", spec)
            vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=spec, overwrite=overwrite)
            print(k, "-", "Partial Internals", spec)
            vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=spec, use_internals=True, overwrite=overwrite)

@runnable
def MACEVPT():
    degs = True
    # vpt.run_aimnet_vpt(0, use_degeneracies=False, use_internals=False, use_reaction_path=False)
    # vpt.run_aimnet_vpt(0, use_degeneracies=False, use_internals=True, use_reaction_path=False)
    step_size = None#.25
    for k in range(0, -70, -10):
        print(k, "-", "Full")
        vpt.run_mace_vpt(k, use_degeneracies=degs, step_size=step_size)
        print(k, "-", "Full Internals")
        vpt.run_mace_vpt(k, use_degeneracies=degs, use_internals=True, step_size=step_size)
        for spec in [
            [-4, -3, -2, -1],
            [-7, -6, -5, -4, -3, -2, -1],
            [-11, -10, -8, -1],
            [-8, -1],
            [-11, -1],
            [-11, -4, -3, -2, -1],
            [-11, -7, -6, -5, -4, -3, -2, -1],
            [-1],
        ]:
            print(k, "-", "Partial", spec)
            vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=spec)
            print(k, "-", "Partial Internals", spec)
            vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=spec, use_internals=True)

@inactive
def AIMNetVPTMore():
    degs = True
    # vpt.run_aimnet_vpt(0, use_degeneracies=False, use_internals=False, use_reaction_path=False)
    # vpt.run_aimnet_vpt(0, use_degeneracies=False, use_internals=True, use_reaction_path=False)
    for k in range(-53, -70, -10):
        print(k, "-", "Full")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs)
        print(k, "-", "Full Internals")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs, use_internals=True)
        for spec in [
            [-4, -3, -2, -1],
            [-7, -6, -5, -4, -3, -2, -1],
            [-11, -10, -8, -1],
            [-8, -1],
            [-11, -1],
            [-11, -4, -3, -2, -1],
            [-11, -7, -6, -5, -4, -3, -2, -1],
            [-1],
        ]:
            print(k, "-", "Partial")
            vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=spec)
            print(k, "-", "Partial Internals")
            vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=spec, use_internals=True)


@inactive
def B3LYPVPT():
    degs = True
    vpt.run_b3lyp_vpt(1, use_degeneracies=degs, use_internals=False, use_reaction_path=False)
    vpt.run_b3lyp_vpt(1, use_degeneracies=degs, use_internals=True, use_reaction_path=False)
    for k in [1, 2, 3, 4, 5, 6, 7]:
        print(k, "-", "Full")
        vpt.run_gaussian_vpt('b3lyp', k, use_degeneracies=degs)
        print(k, "-", "Full Internals")
        vpt.run_gaussian_vpt('b3lyp', k, use_degeneracies=degs, use_internals=True)
        for spec in [
            [-4, -3, -2, -1],
            [-7, -6, -5, -4, -3, -2, -1],
            [-11, -10, -8, -1],
            [-8, -1],
            [-11, -1],
            [-11, -4, -3, -2, -1],
            [-11, -7, -6, -5, -4, -3, -2, -1],
            [-1],
        ]:
            print(k, "-", "Partial", spec)
            vpt.run_gaussian_vpt('b3lyp', k, use_degeneracies=degs, mode_selection=spec)
            print(k, "-", "Partial Internals", spec)
            vpt.run_gaussian_vpt('b3lyp', k, use_degeneracies=degs, mode_selection=spec, use_internals=True)
        # print(k, "-", "OH")
        # vpt.run_gaussian_vpt('b3lyp', k, use_degeneracies=degs, mode_selection=[-1])
        # print(k, "-", "OH Internals")
        # vpt.run_gaussian_vpt('b3lyp', k, use_degeneracies=degs, mode_selection=[-1], use_internals=True)

@inactive
def wb97VPT():
    degs = True
    vpt.run_wb97_vpt(1, use_internals=False, use_reaction_path=False)
    vpt.run_wb97_vpt(1, use_internals=True, use_reaction_path=False)
    for k in [
        1,
        2, 3,
        4,
        5, 6,
        7
    ]:
        print(k, "-", "Full")
        vpt.run_wb97_vpt(k, use_degeneracies=degs)
        print(k, "-", "Full Internals")
        vpt.run_wb97_vpt(k, use_degeneracies=degs, use_internals=True)
        for spec in [
            [-4, -3, -2, -1],
            [-7, -6, -5, -4, -3, -2, -1],
            [-11, -10, -8, -1],
            [-8, -1],
            [-11, -1],
            [-11, -4, -3, -2, -1],
            [-11, -7, -6, -5, -4, -3, -2, -1],
            [-1],
        ]:
            print(k, "-", "Partial", spec)
            vpt.run_wb97_vpt(k, use_degeneracies=degs, mode_selection=spec)
            print(k, "-", "Partial Internals", spec)
            vpt.run_wb97_vpt(k, use_degeneracies=degs, mode_selection=spec, use_internals=True)

@inactive
def NORPVPT():
        vpt.run_aimnet_vpt(0, use_internals=True, use_degeneracies=False, use_reaction_path=True)
        vpt.run_aimnet_vpt(0, use_internals=False, use_degeneracies=False, use_reaction_path=True)
        vpt.run_wb97_vpt(1, use_internals=False, use_degeneracies=False, use_reaction_path=True)
        vpt.run_wb97_vpt(1, use_internals=True, use_degeneracies=False, use_reaction_path=True)
        vpt.run_b3lyp_vpt(1, use_internals=False, use_degeneracies=False, use_reaction_path=True)
        vpt.run_b3lyp_vpt(1, use_internals=True, use_degeneracies=False, use_reaction_path=True)
        ...



@inactive
def analyzeStretchResults():
    b3lyp_key = 1
    wb97_key = 1
    aimnet_key = 0
    distortion_key = ""

    print()
    base_opts = dict(use_degeneracies=False, use_reaction_path=False)
    print(
        analysis.make_freq_comp_tables(
            [
                ["OH", "CH", "CH", "CH", "HCH", "HCH", "HCH", "HOCH"],
                [1, 2, 3, 4, 5, 6, 7, 12]
            ],
            [
                ["B3LYP", "wB97X-D3", "AIMNet"],
                ["Cart.", "Int.", "Cart.", "Int.", "Cart.", "Int."],
            ],
            [
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=False, **base_opts),
                vpt.get_b3lyp_logfile(b3lyp_key,  use_internals=True, **base_opts),
                vpt.get_wb97_logfile(wb97_key,  use_internals=False, **base_opts),
                vpt.get_wb97_logfile(wb97_key,  use_internals=True, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key,  use_internals=False, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key,  use_internals=True, **base_opts)
            ],
            use_tex=True,
            caption="internal vs cartesian comparison without reaction path projection, no degeneracies, obviously the same",
            label='lab:int_cart_norp_nodeg'
            # use_zero_order=True
        )
    )

    print()
    base_opts = dict(use_degeneracies=True, use_reaction_path=False)
    print(
        analysis.make_freq_comp_tables(
            [
                ["OH", "CH", "CH", "CH", "HCH", "HCH", "HCH", "HOCH"],
                [1, 2, 3, 4, 5, 6, 7, 12]
            ],
            [
                ["B3LYP", "wB97X-D3", "AIMNet"],
                ["Cart.", "Int.", "Cart.", "Int.", "Cart.", "Int."],
            ],
            [
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=False, **base_opts),
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=True, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=False, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=True, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=False, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=True, **base_opts)
            ],
            use_tex=True,
            caption="internal vs cartesian comparison without reaction path projection, with degeneracies, not the same",
            label='lab:int_cart_norp_deg',
            use_deperturbed=False
            # use_zero_order=True
        )
    )

    print()
    base_opts = dict(use_degeneracies=True, use_reaction_path=True)
    print(
        analysis.make_freq_comp_tables(
            [
                ["OH", "CH", "CH", "CH", "HCH", "HCH", "HCH"],
                [1, 2, 3, 4, 5, 6, 7]
            ],
            [
                ["B3LYP", "wB97X-D3", "AIMNet"],
                ["Cart.", "Int.", "Cart.", "Int.", "Cart.", "Int."],
            ],
            [
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=False, **base_opts),
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=True, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=False, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=True, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=False, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=True, **base_opts)
            ],
            use_tex=True,
            caption="internal vs cartesian comparison with reaction path projection and degeneracies" + distortion_key,
            label='lab:int_cart_rp_deg'
            # use_zero_order=True
        )
    )

    print()
    base_opts = dict(use_degeneracies=True, use_reaction_path=True, mode_selection=[-7, -6, -5, -4, -3, -2, -1])
    print(
        analysis.make_freq_comp_tables(
            [
                ["OH", "CH", "CH", "CH", "HCH", "HCH", "HCH"],
                [1, 2, 3, 4, 5, 6, 7]
            ],
            [
                ["B3LYP", "wB97X-D3", "AIMNet"],
                ["Cart.", "Int.", "Cart.", "Int.", "Cart.", "Int."],
            ],
            [
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=False, **base_opts),
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=True, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=False, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=True, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=False, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=True, **base_opts)
            ],
            use_tex=True,
            caption="internal vs cartesian comparison with reaction path projection, stretch-bend subspace" + distortion_key,
            label='lab:int_cart_stretch_bend'
            # use_zero_order=True
        )
    )

    print()
    base_opts = dict(use_degeneracies=True, use_reaction_path=True, mode_selection=[-4, -3, -2, -1])
    print(
        analysis.make_freq_comp_tables(
            [
                ["OH", "CH", "CH", "CH"],
                [1, 2, 3, 4]
            ],
            [
                ["B3LYP", "wB97X-D3", "AIMNet"],
                ["Cart.", "Int.", "Cart.", "Int.", "Cart.", "Int."],
            ],
            [
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=False, **base_opts),
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=True, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=False, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=True, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=False, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=True, **base_opts)
            ],
            use_tex=True,
            caption="internal vs cartesian comparison with reaction path projection, stretch subspace" + distortion_key,
            label='lab:int_cart_stretches'
            # use_zero_order=True
        )
    )


    b3lyp_key = 3
    wb97_key = 3
    aimnet_key = -30
    distortion_key = " at $\\tau=30\\circ$"

    print()
    base_opts = dict(use_degeneracies=True, use_reaction_path=True)
    print(
        analysis.make_freq_comp_tables(
            [
                ["OH", "CH", "CH", "CH", "HCH", "HCH", "HCH"],
                [1, 2, 3, 4, 5, 6, 7]
            ],
            [
                ["B3LYP", "wB97X-D3", "AIMNet"],
                ["Cart.", "Int.", "Cart.", "Int.", "Cart.", "Int."],
            ],
            [
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=False, **base_opts),
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=True, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=False, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=True, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=False, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=True, **base_opts)
            ],
            use_tex=True,
            caption="internal vs cartesian comparison with reaction path projection and degeneracies" + distortion_key,
            label='lab:int_cart_norp_deg'
            # use_zero_order=True
        )
    )

    print()
    base_opts = dict(use_degeneracies=True, use_reaction_path=True, mode_selection=[-7, -6, -5, -4, -3, -2, -1])
    print(
        analysis.make_freq_comp_tables(
            [
                ["OH", "CH", "CH", "CH", "HCH", "HCH", "HCH"],
                [1, 2, 3, 4, 5, 6, 7]
            ],
            [
                ["B3LYP", "wB97X-D3", "AIMNet"],
                ["Cart.", "Int.", "Cart.", "Int.", "Cart.", "Int."],
            ],
            [
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=False, **base_opts),
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=True, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=False, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=True, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=False, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=True, **base_opts)
            ],
            use_tex=True,
            caption="internal vs cartesian comparison with reaction path projection, stretch-bend subspace" + distortion_key,
            label='lab:int_cart_stretch_bend'
            # use_zero_order=True
        )
    )

    print()
    base_opts = dict(use_degeneracies=True, use_reaction_path=True, mode_selection=[-4, -3, -2, -1])
    print(
        analysis.make_freq_comp_tables(
            [
                ["OH", "CH", "CH", "CH"],
                [1, 2, 3, 4]
            ],
            [
                ["B3LYP", "wB97X-D3", "AIMNet"],
                ["Cart.", "Int.", "Cart.", "Int.", "Cart.", "Int."],
            ],
            [
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=False, **base_opts),
                vpt.get_b3lyp_logfile(b3lyp_key, use_internals=True, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=False, **base_opts),
                vpt.get_wb97_logfile(wb97_key, use_internals=True, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=False, **base_opts),
                vpt.get_aimnet_logfile(aimnet_key, use_internals=True, **base_opts)
            ],
            use_tex=True,
            caption="internal vs cartesian comparison with reaction path projection, stretch subspace" + distortion_key,
            label='lab:int_cart_stretches'
            # use_zero_order=True
        )
    )

    return

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
        "AIMNet Cartesians",
        vpt.get_aimnet_logfile,
        0
    )

    print_freq_comps_info(
        "AIMNet Internals",
        vpt.get_aimnet_logfile,
        0,
        use_internals=True
    )

    print_freq_comps_info(
        "WB97 Cartesians",
        vpt.get_wb97_logfile,
        1
    )

    print_freq_comps_info(
        "WB97 Internals",
        vpt.get_wb97_logfile,
        1,
        use_internals=True
    )

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

@inactive
def plot_OH_stretch_surface():
    eng_array = []
    spec = list(range(-60, 10, 10))
    for k in spec:
        stretch_log = vpt.get_aimnet_logfile(k,
                                             use_degeneracies=False,
                                             # mode_selection=[-7, -6, -5, -4, -3, -2, -1],
                                             use_internals=False)
        energies = VPTAnalyzer(stretch_log).zero_order_energies[1]
        eng_array.append(energies[1:4])

    g_matrix = ...

    fig = None
    for e_list in np.array(eng_array).T[:1]:
        fig = plt.Plot(*util.symmetrize_potential(spec, e_list), figure=fig)
    fig.show()

@inactive
def get_timings():
    import time
    # anth = Molecule.from_string('anthracene')
    # print(anth.atoms)
    # print(anth.coords.tolist())

    anth = Molecule(
        ('C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
        [[0.43459809047069276, 2.6947451134350238, 0.5528980961258969],
         [2.520675463818493, 1.1675342601855154, 0.24241312524169933],
         [4.923480398080228, 2.1500885928527023, 0.44682842213341023],
         [6.977204717144367, 0.5904344629023636, 0.12967926483081624],
         [6.594754633724675, -1.9348663395912933, -0.3884626606942818],
         [4.174351071959761, -2.8826343951329685, -0.5857689271872423],
         [2.1087305270311694, -1.3528545366522018, -0.2747515334415861],
         [-0.3148694632808378, -2.2896810925617905, -0.4697789471444097],
         [-2.4430404427727113, -0.7812545199788724, -0.16318704502444908],
         [-4.852547106310639, -1.7330836248552046, -0.3612610497362914],
         [-6.9411803237457566, -0.20442243344539196, -0.05046827680719325],
         [-6.501740332497076, 2.312928728793216, 0.4661066639580801],
         [-4.101836698319057, 3.2855946728605163, 0.6684381902950978],
         [-1.9887771604883142, 1.7651290860727675, 0.3593489528453505],
         [5.120868589905046, 4.151100146685542, 0.8572719910865615],
         [8.888588217196185, 1.3318669533675938, 0.2840940110725707],
         [8.211394448672712, -3.1677259052687283, -0.6391565214029209],
         [3.888684109420547, -4.871249406030155, -0.9938461283850493],
         [-0.6861635380524018, -4.288174754199079, -0.8798628702851997],
         [-5.166146074027914, -3.6836911256623033, -0.7614626029433617],
         [-8.809015740505894, -1.016743221315567, -0.2192997841238593],
         [-8.209323389290347, 3.473826817783778, 0.7020195168277881],
         [-3.8286899981332465, 5.283132519754418, 1.078208112758665]]
    )
    internals = anth.get_bond_zmatrix()
    iterations = 1
    order = 3
    with Timer("numerical", file=None, number=iterations):
        for _ in range(iterations):
            ints_by_carts_numerical = anth.modify(internals=internals).get_cartesians_by_internals(order,
                                                                                                   strip_embedding=True,
                                                                                                   reembed=True,
                                                                                                   method='classic',
                                                                                                   analytic_deriv_order=0,
                                                                                                   stencil=5
                                                                                                   )

    with Timer("analytic", file=None, number=iterations):
        for _ in range(iterations):
            analytic = anth.modify(internals=internals).get_cartesians_by_internals(order,
                                                                                    strip_embedding=True,
                                                                                    reembed=True,
                                                                                    method='fast'
                                                                                    )

    # print(
    #     np.round(analytic[1] - ints_by_carts_numerical[1], 8)
    # )

    # with Timer("direct", file=None, number=iterations):
    #     for _ in range(iterations):
    #         direct = nput.internal_coordinate_tensors(
    #             anth.coords,
    #             coordops.extract_zmatrix_internals(internals),
    #             order=order,
    #             return_inverse=True
    #         )


if __name__ == '__main__':
    # print(expansions.methanol_gaussian_zmatrix)

    for t in tests:
        with Timer(t.__name__, file=None):
            print("="*50)
            print(t.__name__, ":")
            t()
