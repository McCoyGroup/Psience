import numpy as np

import McUtils.Plots as plt
import McUtils.Numputils as nput

__all__ = [
    "plot_energy_levels"
]

def plot_energy_levels(energy_list, figure: plt.Graphics = None,
                       x_list=None,
                       bar_spacing=.025,
                       bar_width=None,
                       plot_range=None,
                       color='black',
                       ticks=None,
                       labels=None,
                       bar_styles=None,
                       **styles):

    if isinstance(energy_list, dict): # replace with collections.Mapping...
        labels = energy_list.keys()
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
    if plot_range is not None and not nput.is_numeric(plot_range[0]):
        x_range = plot_range[0]
    else:
        if x_list is not None:
            x_list = np.asanyarray(x_list)
            x_min = np.min(x_list)
            x_max = np.max(x_list)
            diff = (x_max - x_min) / (1 - 2*bar_spacing)
            x_range = [x_min - bar_spacing * diff, x_max + bar_spacing * diff]
        else:
            x_range = [0, 1]
    if x_list is None:
        neng = len(energy_list)
        x_list = [[i / neng + bar_spacing, (i + 1) / neng - bar_spacing] for i in range(neng)]
    elif nput.is_numeric(x_list[0]):
        x_list = np.asanyarray(x_list)
        x_min = np.min(x_list)
        x_max = np.max(x_list)
        if bar_width is None:
            bar_width = np.min(np.diff(np.sort(x_list)))
        pad = bar_width / 2
        x_range = [x_min - pad, x_max + pad]
        s = bar_spacing * (x_range[-1] - x_range[0])
        x_list = [[x - pad + s, x + pad - s] for x in x_list]

    x_list = np.asanyarray(x_list)
    x_min = np.min(x_list)
    x_max = np.max(x_list)
    diff = (x_max - x_min) / (1 - 2 * bar_spacing)
    cur_range = [x_min - bar_spacing * diff, x_max + bar_spacing * diff]
    x_list = np.reshape(nput.vec_rescale(x_list.flatten(), x_range, cur_range=cur_range), x_list.shape)

    if plot_range is None or nput.is_numeric(plot_range[0]):
        plot_range = [x_range, plot_range]

    if ticks is None or nput.is_numeric(ticks[0]):
        if labels is None:
            x_ticks = []
        else:
            tick_pos = [(a+b)/2 for a,b in x_list]
            x_ticks = (tick_pos, dict(labels=labels))
        ticks = [x_ticks, ticks]
    elif labels is not None:
        x_ticks, y_ticks = ticks
        if nput.is_numeric(x_ticks[0]):
            x_ticks = (x_ticks, dict(labels=labels))
        ticks = [x_ticks, y_ticks]

    if bar_styles is None:
        bar_styles = {}
    if isinstance(bar_styles, dict):
        bar_styles = [bar_styles] * len(energy_list)
    if color is None or isinstance(color, str):
        color = [color] * len(energy_list)

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
    return figure