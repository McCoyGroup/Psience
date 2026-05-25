import numpy as np

import McUtils.Plots as plt
import McUtils.Numputils as nput

__all__ = [
    "plot_energy_levels"
]


def plot_energy_levels(energy_list, figure: plt.Graphics = None,
                       x_list=None,
                       bar_spacing=.2,
                       bar_width=None,
                       plot_range=None,
                       color='black',
                       ticks=None,
                       labels=None,
                       bar_styles=None,
                       connect=False,
                       connection_style=None,
                       scaled_x=False,
                       **styles):

    if isinstance(energy_list, dict): # replace with collections.Mapping...
        labels = list(energy_list.keys())
        energy_list = list(energy_list.values())

    if isinstance(energy_list[0], dict):
        if 'x' in energy_list[0]:
            x_list = [e['x'] for e in energy_list]
        if bar_styles is None:
            bar_styles = [None] * len(energy_list)
        _ = []
        for el, b in zip(energy_list, bar_styles):
            if b is None: b = {}
            b = dict(el, **b)
            b.pop('x', None)
            b.pop('y')
            _.append(b)
        bar_styles = _
        energy_list = [e['y'] for e in energy_list]

    if plot_range is None and figure is not None:
        plot_range = figure.plot_range

    if x_list is None:
        neng = len(energy_list)
        x_list = np.linspace(0, 1, neng)
    if x_list is not None and not scaled_x:
        x_list = np.asanyarray(x_list)
        x_min = np.min(x_list)
        x_max = np.max(x_list)
        if bar_width is None:
            d = np.abs(np.diff(np.sort(np.array(x_list).flatten())))
            d[d < 1e-8] = np.max(2) * 2
            bar_width = np.min(d)
        x_range = [x_min - bar_width / 2, x_max + bar_width / 2]
    elif plot_range is not None and not nput.is_numeric(plot_range[0]):
        scaled_x = True
        x_range = plot_range[0]
    else:
        x_range = [0, 1]
    if nput.is_numeric(x_list[0]):
        x_list = np.asanyarray(x_list)
        # x_min = np.min(x_list)
        # x_max = np.max(x_list)
        if bar_width is None:
            d = np.abs(np.diff(np.sort(np.array(x_list).flatten())))
            d[d < 1e-8] = np.max(2) * 2
            bar_width = np.min(d)
        pad = bar_width / 2
        s = bar_spacing * bar_width
        x_list = [[x - pad + s, x + pad - s] for x in x_list]

    x_list = np.asanyarray(x_list)
    x_min = np.min(x_list)
    x_max = np.max(x_list)
    if bar_width is None:
        d = np.abs(np.diff(np.sort(np.array(x_list).flatten())))
        d[d < 1e-8] = np.max(2) * 2
        bar_width = np.min(d)
    if scaled_x:
        cur_range = [x_min - bar_width/2, x_max + bar_spacing/2]
        x_list = np.reshape(nput.vec_rescale(x_list.flatten(), x_range, cur_range=cur_range), x_list.shape)


    if plot_range is None or nput.is_numeric(plot_range[0]):
        plot_range = [x_range, plot_range]
    elif plot_range is False:
        plot_range = None

    if ticks is False:
        ticks = None
    elif ticks is None or ticks is True or (len(ticks) > 0 and nput.is_numeric(ticks[0])):
        if labels is None:
            x_ticks = []
        else:
            tick_pos = [(a+b)/2 for a,b in x_list]
            x_ticks = (tick_pos, dict(labels=labels))
        ticks = [x_ticks, ticks]
    elif labels is not None:
        x_ticks, y_ticks = ticks
        if len(x_ticks) > 0 and nput.is_numeric(x_ticks[0]):
            x_ticks = (x_ticks, dict(labels=labels))
        ticks = [x_ticks, y_ticks]
    elif ticks is False:
        ticks = None

    if bar_styles is None:
        bar_styles = {}
    if isinstance(bar_styles, dict):
        bar_styles = [bar_styles] * len(energy_list)
    if color is None or isinstance(color, str):
        color = [color] * len(energy_list)

    if connect:
        _ = []
        for (a,b), (c,d) in zip(x_list[:-1], x_list[1:]):
            if b > c and a < b:
                b, a = a, b
            _.append([a, b])
        if len(x_list) > 1:
            (a, b), (c, d) = x_list[-1], x_list[-2]
            if b < c and a < b:
                b, a = a, b
            _.append([a, b])
        else:
            _.append(x_list[-1])
        x_list = _


    for i, (f, x, c, s) in enumerate(zip(energy_list, x_list, color, bar_styles)):
        style_dict = dict(styles, **s)
        c = style_dict.pop('color', c)
        figure = plt.HorizontalLinePlot(
            x,
            f,
            figure=figure,
            color=c,
            plot_range=plot_range,
            ticks=ticks,
            **style_dict
        )
    if connect:
        if connection_style is None:
            connection_style = {'linestyle':'dashed'}
        dat = list(zip(energy_list, x_list, color, bar_styles))
        for i, ((f, x, c, s), (f2, x2, c2, s2)) in enumerate(zip(dat[:-1], dat[1:])):
            style_dict = styles | connection_style | s | s2
            c = style_dict.pop('color', c)
            for fp in zip(f, f2):
                plt.Plot(
                    [x[-1], x2[0] if abs(x2[0] - x[-1]) < abs(x2[1] - x[-1]) else x2[1]],
                    fp,
                    figure=figure,
                    color=c,
                    **style_dict
                )
    return figure