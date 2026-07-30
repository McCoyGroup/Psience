import numpy as np
import McUtils.Plots as plt
import McUtils.Numputils as nput
__all__ = ['plot_energy_levels']

def plot_energy_levels(energy_list, figure: plt.Graphics=None, x_list=None, bar_spacing=0.2, bar_width=None, plot_range=None, color='black', ticks=None, labels=None, primary_label=None, bar_styles=None, connect=False, connection_style=None, scaled_x=False, **styles):
    ...