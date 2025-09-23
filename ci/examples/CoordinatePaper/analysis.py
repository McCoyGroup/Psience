import os.path
import numpy as np
import typing
import McUtils.Plots as plt
from McUtils.Formatters import TableFormatter
import McUtils.Numputils as nput
import McUtils.Formatters as mfmt
from Psience.VPT2 import VPTAnalyzer

from .vpt import get_log_generator, LevelsOfTheory
from .modes import get_cpmo3_modes

_analyzer_cache = {}
def get_energies(log_file:'str|typing.Callable', *logfile_args,
                 energy_type='auto',
                 cache=None,
                 **logfile_opts):
    if isinstance(log_file, LevelsOfTheory) or (isinstance(log_file, str) and not os.path.isfile(log_file)):
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

def get_plot_data(k, log_gen,
                  key_function=None,
                  plot_vals=None,
                  zpe_correct=False,
                  anharmonic_only=False,
                  **job_opts
                  ):
    if isinstance(log_gen, (LevelsOfTheory, str)):
        log_gen = LevelsOfTheory(log_gen)
        if key_function is None:
            if log_gen == LevelsOfTheory.AIMNet2:
                key_function = lambda i : -i
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
        anharmonic_only=anharmonic_only
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
        anharmonic_only=anharmonic_only
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
