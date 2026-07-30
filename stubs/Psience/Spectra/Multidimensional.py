import numpy as np, json, abc, collections
import McUtils.Plots as plt
from McUtils.Data import UnitsData
import McUtils.Numputils as nput
__all__ = ['TwoDimensionalSpectrum']

class TwoDimensionalSpectrum:
    """
    Base class to support spectral operation
    """

    def __init__(self, freq1, freq2, intensities, **meta):
        """
        :param frequencies: frequency list
        :type frequencies: np.ndarray
        :param intensities: intensity list
        :type intensities: np.ndarray
        :param meta: metadata
        :type meta:
        """
        ...

    def take_subspectrum(self, sample_x, sample_y):
        """
        Takes a subset of frequencies/intensities specified by `pos`

        :param pos:
        :type pos:
        :return:
        :rtype:
        """
        ...

    def frequency_filter(self, freq_span_x, freq_span_y):
        """
        **LLM Docstring**

        Restrict the spectrum to a rectangular window of frequency values along both axes, via `take_subspectrum`.

        :param freq_span_x: the `(min, max)` frequency bounds to keep along the `freq1` axis
        :type freq_span_x: tuple[float, float]
        :param freq_span_y: the `(min, max)` frequency bounds to keep along the `freq2` axis
        :type freq_span_y: tuple[float, float]
        :return: the frequency-restricted subspectrum
        :rtype: TwoDimensionalSpectrum
        """
        ...

    def intensity_filter(self, int_min, int_max):
        """
        **LLM Docstring**

        Restrict the spectrum to the rectangular index range spanning the points whose intensity falls within `[int_min, int_max]`. Note the intensity test uses the bitwise `&` operator rather than `and`, which for numeric `int_min <= self.intensities` (a scalar comparison producing a bare bool) does not behave as intended for filtering an array by both bounds.

        :param int_min: the minimum intensity to include
        :type int_min: float
        :param int_max: the maximum intensity to include
        :type int_max: float
        :return: the intensity-restricted subspectrum
        :rtype: TwoDimensionalSpectrum
        """
        ...

    def clip(self, int_min, int_max, clip_abs=True):
        """
        **LLM Docstring**

        Clip the spectrum's intensities to a `[int_min, int_max]` range, either clipping the absolute magnitude while preserving sign (and zeroing out points below `int_min` entirely) or clipping the raw signed values directly.

        :param int_min: the minimum intensity magnitude (or raw value, if `clip_abs` is `False`) to keep
        :type int_min: float
        :param int_max: the maximum intensity magnitude (or raw value) to clip to
        :type int_max: float
        :param clip_abs: whether to clip based on absolute magnitude (preserving sign) rather than the raw signed values
        :type clip_abs: bool
        :return: the clipped spectrum
        :rtype: TwoDimensionalSpectrum
        """
        ...
    default_styles = {'cmap': 'RdBu', 'levels': 10}
    default_line_style = {'colors': 'black'}

    def plot(self, plot_filled=True, contour_line_style=None, figure=None, symmetric_range=True, remove_baseline=True, vmin=None, vmax=None, levels=None, **opts):
        """
        **LLM Docstring**

        Render the 2D spectrum as a filled and/or line contour plot, optionally removing the median baseline from the intensities first and using a symmetric (about zero) color/level range by default.

        :param plot_filled: whether to draw a filled contour plot (via `ContourPlot`)
        :type plot_filled: bool
        :param contour_line_style: style options for an overlaid contour-line plot; `True`/`None` for default line styling, `False` to omit the line plot entirely (only meaningful together with `plot_filled`)
        :type contour_line_style: dict | bool | None
        :param figure: an existing figure to draw into
        :type figure: object | None
        :param symmetric_range: whether to force the color/level range to be symmetric about zero
        :type symmetric_range: bool
        :param remove_baseline: whether to subtract the median intensity before plotting
        :type remove_baseline: bool
        :param vmin: the minimum value for the color scale; derived from the data if not given
        :type vmin: float | None
        :param vmax: the maximum value for the color scale; derived from the data if not given
        :type vmax: float | None
        :param levels: the number of contour levels, or explicit level values
        :type levels: int | Iterable[float] | None
        :param opts: extra plotting options forwarded to `ContourPlot`/`ContourLinePlot`
        :type opts: dict
        :return: the resulting figure
        :rtype: object
        """
        ...