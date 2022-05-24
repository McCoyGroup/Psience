"""
Provides a general base spectrum class that can be extended to new fancy spectral forms
"""

import numpy as np, json, abc
import McUtils.Plots as plt

__all__ = [
    "BaseSpectrum",
    "DiscreteSpectrum",
    "ContinuousSpectrum",
    "BroadenedSpectrum"
]

# TODO: should add in support for various flavors of unit conversion
#       like oscillator strengths to km mol^-1 units and whatnot
class BaseSpectrum:
    """
    Base class to support spectral operation
    """
    def __init__(self, frequencies, intensities, **meta):
        """
        :param frequencies: frequency list
        :type frequencies: np.ndarray
        :param intensities: intensity list
        :type intensities: np.ndarray
        :param meta: metadata
        :type meta:
        """
        self.frequencies = frequencies
        self.intensities = intensities
        self.meta = meta

    def take_subspectrum(self, pos):
        """
        Takes a subset of frequencies/intensities specified by `pos`

        :param pos:
        :type pos:
        :return:
        :rtype:
        """

        return type(self)(self.frequencies[pos], self.intensities[pos], **self.meta)

    def __getitem__(self, item):

        f = self.frequencies[item]
        if isinstance(f, (int, float, np.integer, np.floating)):
            return (f, self.intensities[item])
        else:
            return self.take_subspectrum(item)

    def frequency_filter(self, freq_min, freq_max):
        """
        Filters by frequencies >= `freq_min` and <= `freq_max`

        :param freq_min: min frequency
        :type freq_min: float
        :param freq_max: max frequency
        :type freq_max: float
        :return: subspectrum
        :rtype: BaseSpectrum
        """

        freqs = self.frequencies
        return self.take_subspectrum(np.where(np.logical_and(freq_min <= freqs, freqs <= freq_max))[0])

    def intensity_filter(self, int_min, int_max):
        """
        Filters by intensities >= `int_min` and <= `int_max`

        :param int_min: min intensity
        :type int_min: float
        :param int_max: max intensity
        :type int_max: float
        :return: subspectrum
        :rtype: BaseSpectrum
        """

        ints = self.intensities
        return self.take_subspectrum(np.where(np.logical_and(int_min <= ints, ints <= int_max))[0])


    def save(self, file):
        """
        Saves the spectrum in JSON format

        :param file: str | file-like object
        :type file:
        :return:
        :rtype:
        """

        dump = json.dumps({
            "frequencies": self.frequencies,
            "intensities": self.intensities,
            "metadata": self.meta
        })
        if isinstance(file, str):
            with open(file, "w") as f:
                f.write(dump)
        else:
            file.write(dump)

    @classmethod
    def load(cls, file):
        """
        Saves a spectrum from a JSON file

        :param file: str | file-like object
        :type file:
        :return:
        :rtype:
        """

        if isinstance(file, str):
            with open(file, 'r') as f:
                dump = json.load(f)
        else:
            dump = json.load(file)

        return cls(dump['frequencies'], dump['intensities'], **dump['metadata'])

    @abc.abstractmethod
    def plot(self, **opts):
        """
        A stub so that subclasses can implement their own `plot` methods

        :param opts: plotting options to be fed through to whatever the plotting function uses
        :type opts:
        :return:
        :rtype:
        """
        raise NotImplementedError("{}: unclear how to plot spectrum".format(
            type(self).__name__
        ))


class DiscreteSpectrum(BaseSpectrum):
    """
    Concrete implementation of `BaseSpectrum` that exists
    solely to allow for plotting and broadening.
    """

    def plot(self,
             figure=None,
             plot_style=None,
             **opts
             ):
        """
        Plots a spectrum using `McUtils.Plots.StickSpectrum`

        :param figure: figure to plot the spectrum on
        :type figure: None | McUtils.Plots.Graphics
        :param opts: any of the many, many options supported by `McUtils.Plots.Graphics`
        :type opts:
        :return:
        :rtype:
        """

        # suuuuper straightforward
        return plt.StickPlot(self.frequencies, self.intensities,
                             figure=figure,
                             plot_style=plot_style,
                             **opts
                             )

    def broaden(self, broadening_type="gaussian", breadth=10):
        """
        Applies a broadening to the spectrum

        :param broadening_type:
        :type broadening_type:
        :param breadth:
        :type breadth:
        :return:
        :rtype:
        """
        return BroadenedSpectrum(self.frequencies, self.intensities,
                                 broadening_type=broadening_type, breadth=breadth,
                                 **self.meta
                                 )

class ContinuousSpectrum(BaseSpectrum):
    """
    Concrete implementation of `BaseSpectrum` that exists
    solely to allow for plotting & maybe some day for interchange with like experimental formats
    """

    def plot(self,
             figure=None,
             filled=False,
             plot_style=None,
             **opts
             ):
        """
        Plots a spectrum using `McUtils.Plots.Plot`

        :param figure: figure to plot the spectrum on
        :type figure: None | McUtils.Plots.Graphics
        :param opts: any of the many, many options supported by `McUtils.Plots.Graphics`
        :type opts:
        :return:
        :rtype:
        """

        # suuuuper straightforward
        main = plt.Plot(self.frequencies, self.intensities,
                         figure=figure,
                         plot_style=plot_style,
                         **opts
                         )
        if filled:
            plt.FilledPlot(self.frequencies, self.intensities,
                     figure=main,
                     plot_style=plot_style,
                     **opts
                     )
        return main

class BroadenedSpectrum(BaseSpectrum):
    """
    A stick spectrum with associated broadening function
    """

    def __init__(self, frequencies, intensities, broadening_type="gaussian", breadth=10, **meta):
        """
        :param frequencies:
        :type frequencies:
        :param intensities:
        :type intensities:
        :param broadening_type: the type of broadening to apply (can be any function)
        :type broadening_type: "gaussian" | "lorentzian" | function
        :param breadth: the breadth or list of breads for the peaks in the spectrum
        :type breadth:
        :param meta:
        :type meta:
        """

        super().__init__(frequencies, intensities, **meta)

        self.broadening_type = broadening_type
        self.breadth = breadth

    def _eval_gauss_broadening(self, pts, height, center, breadth, target_zero=.1, adjust_width=True):
        """
        Evaluates a Gaussian centered around `center` with breadth `breadth` at `pts`

        :param pts:
        :type pts: np.ndarray
        :param center:
        :type center: float
        :param breadth:
        :type breadth: float
        :return:
        :rtype:
        """
        from scipy.stats import norm

        if adjust_width:
            h = height
            if h < 1:
                h = 1
            z = target_zero / h
            breadth = np.sqrt(breadth**2/(2*np.log(1/z))) # chosen to make the pdf(center-breadth) == target
        bd = norm(loc=center, scale=breadth)
        return height * bd.pdf(pts) * (np.sqrt(2 * np.pi) * breadth) # remove the normalization

    def _eval_lorentz_broadening(self, pts, height, center, breadth):
        """
        Evaluates a Lorentzian centered around `center` with breadth `breadth` at `pts`

        :param pts:
        :type pts: np.ndarray
        :param center:
        :type center: float
        :param breadth:
        :type breadth: float
        :return:
        :rtype:
        """
        from scipy.stats import cauchy

        bd = cauchy(loc=center, scale=breadth).pdf
        # if height < 1:
        #     height = 1
        return height * bd(pts) * (np.pi * breadth )  # remove the normalization

    def _get_pts(self, step_size=.5, freq_min=None, freq_max=None, stddevs=5, adjust_width=True):
        """
        Evaluates the points needed to plot the broadened spectrum

        :param step_size:
        :type step_size:
        :param freq_min:
        :type freq_min:
        :param freq_max:
        :type freq_max:
        :return:
        :rtype:
        """

        freqs = self.frequencies
        ints = self.intensities

        if isinstance(self.breadth, (int, float, np.integer, np.floating)):
            breadths = [self.breadth] * len(freqs)
        else:
            breadths = self.breadth

        if freq_min is None:
            min_pos = np.argmin(freqs)
            freq_min = freqs[min_pos] - stddevs*breadths[min_pos]
        if freq_max is None:
            max_pos = np.argmax(freqs)
            freq_max = freqs[max_pos] + stddevs*breadths[max_pos]

        bt = self.broadening_type

        if isinstance(bt, str):
            if bt.lower() == 'gaussian':
                bt = self._eval_gauss_broadening
            elif bt.lower() == 'lorentzian':
                bt = self._eval_lorentz_broadening
            else:
                raise ValueError("{}.{}: don't know how to handle broadening type '{}'".format(
                    type(self).__name__,
                    'plot',
                    bt.lower()
                ))

        freq_vals = np.arange(freq_min, freq_max, step_size)
        vals = np.sum([bt(freq_vals, i, c, b, adjust_width=adjust_width) for i, c, b in zip(ints, freqs, breadths)], axis=0)

        return freq_vals, vals

    def plot(self, step_size=.5, freq_min=None, freq_max=None, figure=None, plot_style=None, filled=False, adjust_width=True, **opts):
        """
        Applies the broadening then plots it using `McUtils.Plots.Plot`

        :param step_size: step size to use when getting evaluation points for evaluating the broadening
        :type step_size:
        :param freq_min: min freq for evaluation
        :type freq_min:
        :param freq_max: max freq for evaluation
        :type freq_max:
        :param figure:
        :type figure:
        :param plot_style:
        :type plot_style:
        :param opts:
        :type opts:
        :return:
        :rtype:
        """

        freqs, ints = self._get_pts(step_size=step_size, freq_min=freq_min, freq_max=freq_max, adjust_width=adjust_width)
        if plot_style is None:
            plot_style = {}
        if filled:
            main = plt.FilledPlot(freqs, ints,
                     figure=figure,
                     plot_style=plot_style,
                     **opts
                     )
            plt.Plot(freqs, ints,
                        figure=main,
                        plot_style=plt.Plot.filter_options(plot_style),
                        **plt.Plot.filter_options(opts)
                        )
        else:
            main = plt.Plot(freqs, ints,
                     figure=figure,
                     plot_style=plot_style,
                     **opts
                     )
        return main