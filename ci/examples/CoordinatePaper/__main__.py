import os.path

import numpy as np
from Psience.Molecools import Molecule
from Psience.Modes import NormalModes, LocalizedModes
from Psience.VPT2 import VPTAnalyzer
from McUtils.Formatters import TableFormatter
import McUtils.Formatters as mfmt
import McUtils.Plots as plt
from McUtils.Profilers import Timer
from . import expansions
from . import modes
from . import vpt
from . import util

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
    degs = False
    vpt.run_aimnet_vpt(0, use_degeneracies=False, use_internals=False, use_reaction_path=False)
    vpt.run_aimnet_vpt(0, use_degeneracies=False, use_internals=True, use_reaction_path=False)
    for k in range(0, -70, -10):
        print(k, "-", "Full")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs)
        print(k, "-", "Full Internals")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs, use_internals=True)
        print(k, "-", "Partial")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=[-4, -3, -2, -1])
        print(k, "-", "Partial Internals")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=[-4, -3, -2, -1], use_internals=True)
        print(k, "-", "Partial")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=[-7, -6, -5, -4, -3, -2, -1])
        print(k, "-", "Partial Internals")
        vpt.run_aimnet_vpt(k, use_degeneracies=degs, mode_selection=[-7, -6, -5, -4, -3, -2, -1], use_internals=True)


@inactive
def B3LYPVPT():
    vpt.run_b3lyp_vpt(1, use_internals=False, use_reaction_path=False)
    vpt.run_b3lyp_vpt(1, use_internals=True, use_reaction_path=False)
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
    vpt.run_wb97_vpt(1, use_internals=False, use_reaction_path=False)
    vpt.run_wb97_vpt(1, use_internals=True, use_reaction_path=False)
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

@inactive
def NORPVPT():
        vpt.run_aimnet_vpt(0, use_internals=True, use_degeneracies=False, use_reaction_path=True)
        vpt.run_aimnet_vpt(0, use_internals=False, use_degeneracies=False, use_reaction_path=True)
        vpt.run_wb97_vpt(1, use_internals=False, use_degeneracies=False, use_reaction_path=True)
        vpt.run_wb97_vpt(1, use_internals=True, use_degeneracies=False, use_reaction_path=True)
        vpt.run_b3lyp_vpt(1, use_internals=False, use_degeneracies=False, use_reaction_path=True)
        vpt.run_b3lyp_vpt(1, use_internals=True, use_degeneracies=False, use_reaction_path=True)
        ...


def fprint(tag, a=None):
    if a is None: a, tag = tag, None
    a = np.asanyarray(a)
    if a.ndim == 1: a = a[np.newaxis]
    if tag is not None:
        print(tag, TableFormatter("{:>4.0f}").format(a))
    else:
        print(TableFormatter("{:>4.0f}").format(a))

def make_freq_comp_tables(
        states,
        tags,
        log_files,
        header_spans=None,
        use_zero_order=False,
        use_deperturbed=True,
        use_tex=False,
        **etc
):
    if len(states) == 2 and isinstance(tags[1], (tuple, list, np.ndarray)):
        states, inds = states
    else:
        inds = list(range(1, len(states)+1))

    data_sets = []
    for lf in log_files:
        a1 = VPTAnalyzer(lf)
        if use_zero_order:
            e1 = a1.zero_order_energies[1][inds,]
        elif use_deperturbed:
            e1 = a1.deperturbed_energies[1][inds,]
        else:
            e1 = a1.energies[1][inds,]
        data_sets.append(e1)

    if isinstance(tags[0], str):
        tags = ["States"] + list(tags)
        header_spans = None
    else:
        if header_spans is None:
            header_spans = [
                [1] + [len(tags[1]) // len(tags[0])] * len(tags[0]),
                [1] + [1] * len(tags[1])
            ]
        tags = [
            [""] + list(tags[0]),
            ["States"] + tags[1]
        ]

    if use_tex:
        # def format_mat(lhs, m, label=None, digits=2):
        #     f_og = m.astype(object)
        #     f_og[np.tril_indices_from(f_og, -1)] = ""
        #     fsym = TeX.bold(lhs).as_expr()
        #     TeX.Writer.real_digits = digits
        #     fexpr = fsym.Eq(TeX.Matrix(f_og))
        #     return TeX.Equation(fexpr, label=label).format_tex()
        return mfmt.TeX.Table(
            tags,
            [
                [s] + a.tolist()
                for s, a in zip(
                    states,
                    np.array(data_sets).T
                )
            ],
            number_format="{:4>.0f}",
            resizeable=True,
            header_spans=header_spans,
            **etc
        ).format_tex()
    else:
        return TableFormatter(
            column_formats=[""] + ["{:>4.0f}"] * len(tags),
            headers=tags,
            header_spans=header_spans,
            column_join=" | ",
            column_alignments=["^", ">", ">", ">"]
        ).format([
            [s] + a.tolist()
            for s,a in zip(
                states,
                np.array(data_sets).T
            )
        ])


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
        e0 = a0.deperturbed_energies[1][1:13]
        fprint("Full:", e0)
    else:
        e0 = None
    e1 = a1.deperturbed_energies[1][1:13]
    e2 = a2.deperturbed_energies[1][1:13]
    fprint("Proj:", e1)
    fprint("Subs:", e2)
    # if e0 is not None:
    #     fprint(e0 - e1)
    #     fprint(e0 - e2)
    # fprint(e1 - e2)

@inactive
def analyzeStretchResults():
    b3lyp_key = 1
    wb97_key = 1
    aimnet_key = 0
    distortion_key = ""

    print()
    base_opts = dict(use_degeneracies=False, use_reaction_path=False)
    print(
        make_freq_comp_tables(
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
        make_freq_comp_tables(
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
        make_freq_comp_tables(
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
        make_freq_comp_tables(
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
        make_freq_comp_tables(
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
        make_freq_comp_tables(
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
        make_freq_comp_tables(
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
        make_freq_comp_tables(
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

@runnable
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


if __name__ == '__main__':
    # print(expansions.methanol_gaussian_zmatrix)

    for t in tests:
        with Timer(t.__name__, file=None):
            print("="*50)
            print(t.__name__, ":")
            t()
