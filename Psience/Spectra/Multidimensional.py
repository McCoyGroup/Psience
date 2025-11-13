import numpy as np, json, abc, collections
import McUtils.Plots as plt
from McUtils.Data import UnitsData
import McUtils.Numputils as nput

__all__ = [
    "TwoDimensionalSpectrum"
]

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
        self.freq1 = freq1
        self.freq2 = freq2
        self.intensities = intensities
        self.meta = meta

    def take_subspectrum(self, sample_x, sample_y):
        """
        Takes a subset of frequencies/intensities specified by `pos`

        :param pos:
        :type pos:
        :return:
        :rtype:
        """

        return type(self)(
            self.freq1[sample_x],
            self.freq2[sample_y],
            self.intensities[np.ix_(sample_y, sample_x)],
            **self.meta
        )

    def frequency_filter(self, freq_span_x, freq_span_y):
        fmin_x, fmax_x = freq_span_x
        fmin_y, fmax_y = freq_span_y
        return self.take_subspectrum(
            np.where(np.logical_and(fmin_x <= self.freq1, self.freq1 <= fmax_x))[0],
            np.where(np.logical_and(fmin_y <= self.freq2, self.freq2 <= fmax_y))[0],
        )

    def intensity_filter(self, int_min, int_max):
        pos = np.where(int_min <= self.intensities & self.intensities <= int_max)
        x_min, x_max = np.min(pos[0]), np.max(pos[0])
        y_min, y_max = np.min(pos[1]), np.max(pos[1])
        return self.take_subspectrum(
            np.arange(x_min, x_max+1),
            np.arange(y_min, y_max+1)
        )

    def clip(self, int_min, int_max, clip_abs=True):
        if clip_abs:
            abs_I = np.abs(self.intensities)
            mask = abs_I > int_min
            signs = np.sign(self.intensities) * mask
            return type(self)(
                self.freq1,
                self.freq2,
                signs * np.clip(abs_I, int_min, int_max),
                **self.meta
            )
        else:
            return type(self)(
                self.freq1,
                self.freq2,
                np.clip(self.intensities, int_min, int_max),
                **self.meta
            )

    default_styles = {
        "cmap":'RdBu',
        "levels":10,
    }
    default_line_style = {'colors':'black'}
    def plot(self,
             plot_filled=True,
             contour_line_style=None,
             figure=None,
             symmetric_range=True,
             remove_baseline=True,
             vmin=None,
             vmax=None,
             levels=None,
             **opts
             ):
        # print(self.freq1.shape, self.freq2.shape, self.intensities.shape)
        ints = self.intensities
        if remove_baseline:
            ints = ints - np.median(ints)
        if symmetric_range:
            if vmin is None:
                if vmax is None:
                    vmax = np.max(np.abs(ints))
                vmin = -vmax
            elif vmax is None:
                vmax = -vmin
            if levels is not None and nput.is_int(levels):
                levels = np.linspace(vmin, vmax, levels)
                levels = levels[np.abs(levels) > 1e-12]
        opts.update(
            levels=levels,
            vmin=vmin,
            vmax=vmax
        )
        # print(vmax, vmin)
        if plot_filled:
            base_opts = dict(self.default_styles, **opts)
            figure = plt.ContourPlot(
                self.freq1,
                self.freq2,
                ints,
                figure=figure,
                **base_opts
            )
            if contour_line_style is None or contour_line_style is True:
                contour_line_style = {}
            if contour_line_style is not False:
                base_opts.pop('colors', None)
                base_opts.pop('cmap', None)
                plt.ContourLinePlot(
                    self.freq1,
                    self.freq2,
                    self.intensities,
                    figure=figure,
                    **collections.ChainMap(contour_line_style, base_opts, self.default_line_style)
                )
        else:
            if contour_line_style is None or contour_line_style is True:
                contour_line_style = {}
            figure = plt.ContourLinePlot(
                self.freq1,
                self.freq2,
                self.intensities,
                figure=figure,
                **collections.ChainMap(contour_line_style, opts, self.default_styles, self.default_line_style)
            )

        return figure