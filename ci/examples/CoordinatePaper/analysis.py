import os.path
import numpy as np
import typing
import McUtils.Plots as plt
from McUtils.Formatters import TableFormatter
import McUtils.Devutils as dev
import McUtils.Numputils as nput
import McUtils.Formatters as mfmt
from Psience.VPT2 import VPTAnalyzer

from .vpt import get_log_generator, LevelsOfTheory
from .modes import get_cpmo3_modes

_analyzer_cache = {}
def get_energies(log_file:'str|typing.Callable', *logfile_args,
                 energy_type='auto',
                 cache=None,
                 step_size=None,
                 **logfile_opts):
    if isinstance(log_file, LevelsOfTheory) or (isinstance(log_file, str) and not os.path.isfile(log_file)):
        log_gen = LevelsOfTheory(log_file)
        if any(
                log_gen == x
                for x in [LevelsOfTheory.AIMNet2, LevelsOfTheory.AIMNet2Old, LevelsOfTheory.MACE]
        ):
            logfile_opts['step_size'] = step_size
        log_file = get_log_generator(log_file)
    if not isinstance(log_file, str):
        log_file = log_file(
            *logfile_args,
            **logfile_opts
        )
    if cache is None:
        cache = _analyzer_cache
    if cache is False:
        analyzer = VPTAnalyzer(log_file)
    else:
        if log_file in cache:
            analyzer = cache[log_file]
        else:
            analyzer = VPTAnalyzer(log_file)
            cache[log_file] = analyzer
    return (
        analyzer.deperturbed_energies
            if energy_type == 'deperturbed' else
        analyzer.zero_order_energies
            if energy_type == 'harmonic' else
        analyzer.energies
    )

_job_props_ = ("mode_selection", "step_size")
def get_plot_data(k, log_gen,
                  key_function=None,
                  plot_vals=None,
                  zpe_correct=False,
                  anharmonic_only=False,
                  step_size=None,
                  **job_opts
                  ):
    if isinstance(log_gen, (LevelsOfTheory, str)):
        log_gen = LevelsOfTheory(log_gen)
        if key_function is None:
            if any(
                    log_gen == x
                    for x in [LevelsOfTheory.AIMNet2, LevelsOfTheory.AIMNet2Old, LevelsOfTheory.MACE]
            ):
                key_function = lambda i : -i
                job_opts['step_size'] = step_size
            else:
                key_function = lambda i :( i//10) + 1
        if plot_vals is None:
            # if log_gen == LevelsOfTheory.wB97:
            #     plot_vals = (0, 20, 40, 50)
            # else:
            plot_vals = (0, 10, 20, 30, 40, 50, 60)

    if zpe_correct:
        base_vals = [
            get_energies(
                log_gen,
                key_function(i),
                **job_opts
            )[1][0]
            for i in plot_vals
        ]
    elif anharmonic_only:
        joops = job_opts.copy()
        joops['energy_type'] = 'harmonic'
        base_vals = [
            -get_energies(
                log_gen,
                key_function(i),
                **joops
            )[1][k]
            for i in plot_vals
        ]
    else:
        base_vals = [0] * len(plot_vals)
    return (
        plot_vals,
        [
            get_energies(
                log_gen,
                key_function(i),
                **job_opts
            )[1][k] + b
            for i,b in zip(plot_vals, base_vals)
        ]
    )

def energy_plot(k, log_gen, key_function,
                plot_vals=(0, 10, 20, 30, 40, 50, 60),
                figure=None, plot_styles=None,
                energy_type='deperturbed',
                mode_selection=None,
                step_size=None,
                use_degeneracies=True,
                use_reaction_path=True,
                use_internals=False,
                zpe_correct=False,
                anharmonic_only=False,
                plot_type=None,
                **opts):
    if plot_styles is None:
        plot_styles = {}
    if plot_type is None:
        plot_type = plt.Plot
    job_opts = dict(
        mode_selection=mode_selection,
        use_internals=use_internals,
        use_degeneracies=use_degeneracies,
        energy_type=energy_type,
        zpe_correct=zpe_correct,
        anharmonic_only=anharmonic_only,
        step_size=step_size,
        use_reaction_path=use_reaction_path
    )
    return plot_type(
        *get_plot_data(
            k, log_gen, key_function,
            plot_vals=plot_vals,
            **job_opts
        ),
        figure=figure,
        **opts,
        **plot_styles
    )

def energy_comp_plot(
        k, log_gen, key_function,
        figure=None,
        **opts
):
    figure = energy_plot(k, log_gen, key_function, figure=figure,
                         label='All', **opts)
    energy_plot(k, log_gen, key_function, mode_selection=[-7, -6, -5, -4, -3, -2, -1],
                label='Stretch-Bend', figure=figure, **opts)
    energy_plot(k, log_gen, key_function, mode_selection=[-4, -3, -2, -1], figure=figure,
                label='Stretch',
                **opts)
    energy_plot(k, log_gen, key_function, mode_selection=[-1], figure=figure,
                label='Internals',
                **opts)
    return figure

def plot_aimnet_energies(k, **opts):
    return energy_plot(
        k,
        'aimnet',
        lambda i: -i,
        **opts
    )
def aimnet_comp_plot(k, **opts):
    return energy_comp_plot(
        k,
        'aimnet',
        lambda i :-i,
        **opts
    )

def plot_aimnet_old_energies(k, **opts):
    return energy_plot(
        k,
        'aimnet-old',
        lambda i: -i,
        **opts
    )
def aimnet_old_comp_plot(k, **opts):
    return energy_comp_plot(
        k,
        'aimnet-old',
        lambda i :-i,
        **opts
    )

def plot_mace_energies(k, **opts):
    return energy_plot(
        k,
        'mace',
        lambda i: -i,
        **opts
    )
def mace_comp_plot(k, **opts):
    return energy_comp_plot(
        k,
        'mace',
        lambda i :-i,
        **opts
    )

def plot_b3lyp_energies(k, **opts):
    return energy_plot(
        k,
        'b3lyp',
        lambda i :( i//10) + 1,
        **opts
    )

def b3lyp_comp_plot(k, **opts):
    return energy_comp_plot(
        k,
        'b3lyp',
        lambda i :( i//10) + 1,
        **opts
    )

def plot_wb97_energies(k, **opts):
    return energy_plot(
        k,
        'wb97',
        lambda i :( i//10) + 1,
        # plot_vals=plot_vals,
        **opts
    )

def wb97_comp_plot(k, **opts):
    return energy_comp_plot(
        k,
        'wb97',
        lambda i :( i//10) + 1,
        **opts
    )

default_lot_styles = dict(
    plot_legend=True,
    legend_style=dict(borderaxespad=0, frameon=False, ncol=3),
    # plot_range=[[0, 60], [2800, 3000]],
    image_size=800,
    aspect_ratio=1/1.6,
    use_internals=True,
    energy_type='deperturbed'
)
default_b3lyp_styles = dict(
    linestyle='dashed'
)
default_wb97_styles = dict(
    linestyle='dotted'
)
default_aimnet_styles = dict(
    # linestyle='dotted'
)


def lot_subspace_comp_plot(
        key,
        *,
        b3lyp_styles=None,
        wb97_styles=None,
        aimnet_styles=None,
        prefix="",
        **base_opts
):

    if b3lyp_styles is None:
        b3lyp_styles = {}
    b3lyp_styles = dict(default_b3lyp_styles, **b3lyp_styles)
    if wb97_styles is None:
        wb97_styles = {}
    wb97_styles = dict(default_wb97_styles, **wb97_styles)
    if aimnet_styles is None:
        aimnet_styles = {}
    aimnet_styles = dict(default_aimnet_styles, **aimnet_styles)

    base_opts = dict(default_lot_styles, **base_opts)
    b3lyp_styles = {k: b3lyp_styles[k]
                    for k in b3lyp_styles.keys() - base_opts.keys()}
    wb97_styles = {k: wb97_styles[k]
                   for k in wb97_styles.keys() - base_opts.keys()}
    aimnet_styles = {k: aimnet_styles[k]
                     for k in aimnet_styles.keys() - base_opts.keys()}

    f1 = plot_b3lyp_energies(key,
                             label=prefix + 'B3LYP - All',
                             **b3lyp_styles,
                             **base_opts)
    plot_b3lyp_energies(key,
                        label=prefix + 'B3LYP - Stretch Bend',
                        mode_selection=[-7, -6, -5, -4, -3, -2, -1],
                        figure=f1,
                        **b3lyp_styles,
                        **base_opts)
    plot_b3lyp_energies(key,
                        label=prefix + 'B3LYP - Stretch',
                        mode_selection=[-4, -3, -2, -1],
                        figure=f1,
                        **b3lyp_styles,
                        **base_opts)

    plot_wb97_energies(key,
                       label=prefix + '$\\omega$B97X-D3 All',
                       figure=f1,
                       **wb97_styles,
                       **base_opts)
    plot_wb97_energies(key,
                       label=prefix + '$\\omega$B97X-D3 Stretch-Bend',
                       mode_selection=[-7, -6, -5, -4, -3, -2, -1],
                       figure=f1,
                       **wb97_styles,
                       **base_opts)
    plot_wb97_energies(key,
                       label=prefix + '$\\omega$B97X-D3 Stretch',
                       mode_selection=[-4, -3, -2, -1],
                       figure=f1,
                       **wb97_styles,
                       **base_opts)

    plot_aimnet_energies(key,
                         label=prefix + 'AIMNet2 All',
                         figure=f1,
                         **aimnet_styles,
                         **base_opts)
    plot_aimnet_energies(key,
                         label=prefix + 'AIMNet2 Stretch-Bend',
                         mode_selection=[-7, -6, -5, -4, -3, -2, -1],
                         figure=f1,
                         **aimnet_styles,
                         **base_opts)
    plot_aimnet_energies(key,
                         label=prefix + 'AIMNet2 Stretch',
                         mode_selection=[-4, -3, -2, -1],
                         figure=f1,
                         **aimnet_styles,
                         **base_opts)

    # f1
    # f1.savefig('comp_plot_stretches_stretch_only.svg')
    return f1


def lot_comp_plot(
        key,
        *,
        b3lyp_styles=None,
        wb97_styles=None,
        aimnet_styles=None,
        prefix="",
        figure=None,
        **base_opts
):

    if b3lyp_styles is None:
        b3lyp_styles = {}
    b3lyp_styles = dict(default_b3lyp_styles, **b3lyp_styles)
    if wb97_styles is None:
        wb97_styles = {}
    wb97_styles = dict(default_wb97_styles, **wb97_styles)
    if aimnet_styles is None:
        aimnet_styles = {}
    aimnet_styles = dict(default_aimnet_styles, **aimnet_styles)

    base_opts = dict(default_lot_styles, **base_opts)
    b3lyp_styles = {k: b3lyp_styles[k]
                    for k in b3lyp_styles.keys() - base_opts.keys()}
    wb97_styles = {k: wb97_styles[k]
                   for k in wb97_styles.keys() - base_opts.keys()}
    aimnet_styles = {k: aimnet_styles[k]
                     for k in aimnet_styles.keys() - base_opts.keys()}


    figure = plot_aimnet_energies(key,
                         label=prefix + 'AIMNet2',
                         figure=figure,
                         **aimnet_styles,
                         **base_opts)
    plot_b3lyp_energies(key,
                                 label=prefix + 'B3LYP',
                                 figure=figure,
                                 **b3lyp_styles,
                                 **base_opts)

    plot_wb97_energies(key,
                       label=prefix + '$\\omega$B97X-D3',
                       figure=figure,
                       **wb97_styles,
                       **base_opts)

    return figure


def subspace_comp_plot(
        plot_function,
        key,
        *,
        include_single_space=False,
        prefix='',
        figure=None,
        **base_opts
):

    base_opts = dict(default_lot_styles, **base_opts)

    figure = plot_function(key,
                           label=prefix + 'All',
                           figure=figure,
                           **base_opts)

    plot_function(key,
                  label=prefix + 'Stretch-Bend',
                  figure=figure,
                  mode_selection=[-7, -6, -5, -4, -3, -2, -1],
                  **base_opts)

    plot_function(key,
                  label=prefix + 'Stretch',
                  figure=figure,
                  mode_selection=[-4, -3, -2, -1],
                  **base_opts)
    if include_single_space:
        plot_function(key,
                      label=prefix + 'OH',
                      figure=figure,
                      mode_selection=[-1],
                      **base_opts)

    return figure


def energy_diffence_plot(k, log_gen1, log_gen2,
                         key_function1=None,
                         key_function2=None,
                         plot_vals=None,
                         figure=None, plot_styles=None,
                         energy_type='deperturbed',
                         mode_selection=None,
                         step_size=None,
                         use_reaction_path=True,
                         use_degeneracies=True,
                         use_internals=False,
                         zpe_correct=False,
                         anharmonic_only=False,
                         plot_type=None,
                         **opts):
    if plot_styles is None:
        plot_styles = {}
    if plot_type is None:
        plot_type = plt.Plot
    job_opts = dict(
        mode_selection=mode_selection,
        use_internals=use_internals,
        use_degeneracies=use_degeneracies,
        energy_type=energy_type,
        zpe_correct=zpe_correct,
        anharmonic_only=anharmonic_only,
        step_size=step_size,
        use_reaction_path=use_reaction_path
    )
    plot_vals, d1 = get_plot_data(
        k, log_gen1,
        key_function=key_function1,
        plot_vals=plot_vals,
        **job_opts
    )
    _, d2 = get_plot_data(
        k, log_gen2,
        key_function=key_function2,
        plot_vals=plot_vals,
        **job_opts
    )
    return plot_type(
        plot_vals,
        [e1 - e2 for e1, e2 in zip(d1, d2)],
        figure=figure,
        **opts,
        **plot_styles
    )

plotter_map = {
    LevelsOfTheory.AIMNet2:plot_aimnet_energies,
    LevelsOfTheory.AIMNet2Old:plot_aimnet_old_energies,
    LevelsOfTheory.MACE:plot_mace_energies,
    LevelsOfTheory.B3LYP:plot_b3lyp_energies,
    LevelsOfTheory.wB97:plot_wb97_energies,
}
def resolve_plotter(lot):
    lot = LevelsOfTheory(lot)
    return plotter_map[lot]

label_map = {
    LevelsOfTheory.AIMNet2: 'AIMNet2',
    LevelsOfTheory.AIMNet2Old: 'AIMNet2',
    LevelsOfTheory.MACE: 'MACE-OFF',
    LevelsOfTheory.B3LYP: "B3LYP",
    LevelsOfTheory.wB97: "$\\omega$B97X-D3",
}
label_formats = {
    'lot':'{lot_label}',
    'freq':'{freq_type}{state}',
    'subspace':'{subspace}'
}
freq_types = {
    'deperturbed':'$\\nu$',
    'degenerate':'$\\nu$',
    'harmonic':'$\\omega$'
}
def resolve_freq_type_label(energy_type, anharmonic_only):
    freq_label = freq_types.get(energy_type, energy_type)
    if anharmonic_only:
        freq_label = (
            ("$\\Delta" + freq_label[1:])
                    if freq_label[0] == "$" else
            ("$\\Delta$" + freq_label)
        )
    return freq_label
oh_stretch_pos = [10]
ch_stretch_pos = [9, 8, 7]
ch_bend_pos = [6, 5, 4]
misc_pos = [3, 2, 1, 0]
def resolve_subspace_label(mode_selection):
    if mode_selection is None:
        return 'All'
    else:
        sel = [(12 + s) if s < 0 else s for s in mode_selection]
        oh = oh_stretch_pos in sel
        all_ch = all(c in sel for c in ch_stretch_pos)
        all_bend = all(c in sel for c in ch_bend_pos)
        all_misc = all(c in sel for c in misc_pos)
        if (oh and all_ch and all_bend and all_misc): return 'All'

        any_misc = any(c in sel for c in misc_pos)
        label = []
        if (oh and all_ch):
            label.append("Stretch")
        elif oh:
            label.append("OH Stretch")
        elif all_ch:
            label.append("CH Stretch")
        if all_bend:
            label.append("Bend")
        if any_misc:
            label.append("Extra")
        return "+".join(label)
def resolve_label(*, lot, state, label=None,
                  label_format=None,
                  label_format_join=' - ',
                  energy_type='deperturbed',
                  anharmonic_only=False,
                  mode_selection=None,
                  **opts):
    if label is True:
        if label_format is None:
            label_format = ['lot', 'freq']
            if mode_selection is not None:
                label_format.append('subspace')

    if label is not True: return label

    if not isinstance(label_format, str):
        label_format = label_format_join.join(
            label_formats.get(l, l)
            for l in label_format
        )

    return label_format.format(
        lot_label=label_map.get(LevelsOfTheory(lot)),
        state=state,
        freq_type=resolve_freq_type_label(energy_type, anharmonic_only),
        subspace=resolve_subspace_label(mode_selection)
    )


def plot_generic(*, lot, state, label=None, **opts):
    plotter = resolve_plotter(lot)
    return plotter(state,
                   label=resolve_label(lot=lot, state=state, label=label, **opts),
                   **opts)

default_multi_plot_styles = dict(
    style_list={
        'color': (
            ['001C7F', '017517', '8C0900']
        )
    },
    plot_range=[[0, 60],  None],
    aspect_ratio=1 / 1.6,
    axes_labels=['$\\Delta\\tau$ ($^\\circ$)', 'Frequency (cm$^{-1}$)'],
    plot_legend=True,
    image_size=800,
    legend_style={'ncol': 3, 'frameon': False}
)
def plot_multi(
        *plot_specs: dict,
        energy_type='deperturbed',
        mode_selection=None,
        step_size=None,
        use_degeneracies=None,
        use_reaction_path=False,
        use_internals=None,
        zpe_correct=None,
        anharmonic_only=None,
        label=True,
        lot_styles=None,
        figure=None,
        **global_settings
):
    global_settings = dict(
        default_multi_plot_styles,
        **global_settings
    )
    if lot_styles is None:
        lot_styles = {}
    for i,f in enumerate(plot_specs):
        base_opts = dict(
            energy_type=energy_type,
            mode_selection=mode_selection,
            step_size=step_size,
            use_degeneracies=use_degeneracies,
            use_reaction_path=use_reaction_path,
            use_internals=use_internals,
            zpe_correct=zpe_correct,
            anharmonic_only=anharmonic_only,
            label=label
        )
        f = dict(
            dict(
                {
                    k: v
                    for k, v in base_opts.items()
                    if v is not None
                },
                **lot_styles.get(f['lot'], {})
            ),
            **f
        )

        f['figure'] = f.get('figure', figure)
        s = f['state']
        if dev.is_number(s):
            if i == 0:
                figure = plot_generic(**f, **global_settings)
            else:
                _ = plot_generic(**f) #TDB if I prefer to update the object each iteration or not
        else:
            for j,r in enumerate(s):
                if i == 0 and j == 0:
                    figure = plot_generic(**dict(f, state=r), **global_settings)
                    if f['figure'] is None:
                        f['figure'] = figure
                else:
                    _ = plot_generic(**dict(f, state=r))  # TDB if I prefer to update the object each iteration or not
    return figure

def plot_comb(**styles):
    f1 = plots.plot_mace_energies(2,
                                  )
    plots.plot_mace_energies(3, figure=f1, energy_type='deperturbed',
                             label='MACE - $\\omega_2$',
                             use_reaction_path=False)  # , plot_range=[-190, -180])
    plots.plot_mace_energies(4, figure=f1, energy_type='deperturbed',
                             label='MACE - $\\omega_4$',
                             use_reaction_path=False)  # , plot_range=[-190, -180])

    plots.plot_b3lyp_energies(2, figure=f1, linestyle='dashed', energy_type='deperturbed',
                              label='B3LYP - $\\omega_2$')
    plots.plot_b3lyp_energies(3, figure=f1, linestyle='dashed', energy_type='deperturbed',
                              label='B3LYP - $\\omega_3$')  # , plot_range=[-190, -180])
    plots.plot_b3lyp_energies(4, figure=f1, linestyle='dashed', energy_type='deperturbed',
                              label='B3LYP - $\\omega_4$')  # , plot_range=[-190, -180])

    plots.plot_wb97_energies(2, figure=f1, linestyle='dotted', energy_type='deperturbed',
                             label='$\\omega$B97X-D3 - $\\omega_2$')
    plots.plot_wb97_energies(3, figure=f1, linestyle='dotted', energy_type='deperturbed',
                             label='$\\omega$B97X-D3 - $\\omega_3$')
    plots.plot_wb97_energies(4, figure=f1, linestyle='dotted', energy_type='deperturbed',
                             label='$\\omega$B97X-D3 - $\\omega_4$',
                             image_size=800)
    # f1.savefig('CH_harmonic_freq.svg')
    f1



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

def setup_cpmo3_hessian(struct_id, return_structure=False, **opts):
    mol, (modes, _) = get_cpmo3_modes(struct_id, return_structure=True, **opts)
    if hasattr(modes, 'coords_by_modes'):
        exp = [modes.coords_by_modes]
    else:
        exp = modes[1]

    f = nput.tensor_reexpand(
        exp,
        mol.potential_derivatives
    )[1]
    g = nput.metric_tensor(modes[0], masses=mol.atomic_masses)
    if return_structure:
        return mol, (f, g)
    else:
        return f,g
