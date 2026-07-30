"""
Provides a general base spectrum class that can be extended to new fancy spectral forms
"""
import numpy as np, json, abc
import McUtils.Plots as plt
from McUtils.Data import UnitsData
__all__ = ['BaseSpectrum', 'DiscreteSpectrum', 'ContinuousSpectrum', 'BroadenedSpectrum', 'MultiSpectrum']

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
        ...

    def take_subspectrum(self, pos):
        """
        Takes a subset of frequencies/intensities specified by `pos`

        :param pos:
        :type pos:
        :return:
        :rtype:
        """
        ...

    def __getitem__(self, item):
        """
        **LLM Docstring**

        Index into the spectrum: a scalar index returns the corresponding `(frequency, intensity)` pair, while any other index (slice, array, etc.) returns a subspectrum via `take_subspectrum`.

        :param item: the index or slice to apply
        :type item: int | slice | np.ndarray
        :return: a `(frequency, intensity)` pair for a scalar index, or a subspectrum otherwise
        :rtype: tuple[float, float] | BaseSpectrum
        """
        ...

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
        ...

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
        ...

    def save(self, file):
        """
        Saves the spectrum in JSON format

        :param file: str | file-like object
        :type file:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def load(cls, file):
        """
        Saves a spectrum from a JSON file

        :param file: str | file-like object
        :type file:
        :return:
        :rtype:
        """
        ...

    @abc.abstractmethod
    def plot(self, **opts):
        """
        A stub so that subclasses can implement their own `plot` methods

        :param opts: plotting options to be fed through to whatever the plotting function uses
        :type opts:
        :return:
        :rtype:
        """
        ...

class DiscreteSpectrum(BaseSpectrum):
    """
    Concrete implementation of `BaseSpectrum` that exists
    solely to allow for plotting and broadening.
    """

    @classmethod
    def from_transition_moments(cls, frequencies, transition_moments, **meta):
        """
        Assumes frequencies and transition moments in a.u.

        :param frequencies:
        :type frequencies:
        :param transition_moments:
        :type transition_moments:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def from_raman_moments(cls, frequencies, transition_polarizabilities, pump_frequency=0, **meta):
        """
        **LLM Docstring**

        Build a discrete Raman spectrum from transition frequencies and polarizability transition moments (in atomic units), converting frequencies to wavenumbers and computing intensities from the transition-polarizability magnitudes scaled by the (pump-shifted) frequency to the fourth power.

        :param frequencies: the transition frequencies, in Hartrees
        :type frequencies: np.ndarray
        :param transition_polarizabilities: the polarizability transition-moment tensors
        :type transition_polarizabilities: np.ndarray
        :param pump_frequency: the pump laser frequency to add before the frequency^4 scaling
        :type pump_frequency: float
        :param meta: extra metadata stored on the spectrum
        :type meta: dict
        :return: the constructed Raman spectrum
        :rtype: DiscreteSpectrum
        """
        ...

    def normalize(self, which=None):
        """
        **LLM Docstring**

        Build a copy of the spectrum with intensities normalized by either the overall maximum intensity or a specific reference intensity.

        :param which: the index of a specific transition to normalize against; if `None`, normalizes by the maximum intensity
        :type which: int | None
        :return: the normalized spectrum
        :rtype: DiscreteSpectrum
        """
        ...

    def plot(self, figure=None, plot_style=None, **opts):
        """
        Plots a spectrum using `McUtils.Plots.StickSpectrum`

        :param figure: figure to plot the spectrum on
        :type figure: None | McUtils.Plots.Graphics
        :param opts: any of the many, many options supported by `McUtils.Plots.Graphics`
        :type opts:
        :return:
        :rtype:
        """
        ...

    def broaden(self, breadth=10, *, broadening_type='gaussian', noising_function=None):
        """
        Applies a broadening to the spectrum

        :param broadening_type:
        :type broadening_type:
        :param breadth:
        :type breadth:
        :return:
        :rtype:
        """
        ...

class ContinuousSpectrum(BaseSpectrum):
    """
    Concrete implementation of `BaseSpectrum` that exists
    solely to allow for plotting & maybe some day for interchange with like experimental formats
    """

    def plot(self, figure=None, filled=False, plot_style=None, **opts):
        """
        Plots a spectrum using `McUtils.Plots.Plot`

        :param figure: figure to plot the spectrum on
        :type figure: None | McUtils.Plots.Graphics
        :param opts: any of the many, many options supported by `McUtils.Plots.Graphics`
        :type opts:
        :return:
        :rtype:
        """
        ...

class BroadenedSpectrum(BaseSpectrum):
    """
    A stick spectrum with associated broadening function
    """

    def __init__(self, frequencies, intensities, broadening_type='gaussian', breadth=10, noising_function=None, **meta):
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
        ...

    def _eval_gauss_broadening(self, pts, height, center, breadth, target_zero=0.1, renormalize=True, adjust_width=True):
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
        ...

    def _eval_lorentz_broadening(self, pts, height, center, breadth, renormalize=False):
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
        ...

    def _get_pts(self, step_size=0.5, freq_min=None, freq_max=None, stddevs=5, adjust_width=True, renormalize=True):
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
        ...

    def plot(self, step_size=0.5, freq_min=None, freq_max=None, figure=None, plot_style=None, filled=False, adjust_width=True, renormalize=True, **opts):
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
        ...

class MultiSpectrum:
    """
    A wrapper for multiple spectra, really just for the plotting convenience
    """
    '\n        Base class to support spectral operation\n        '

    def __init__(self, spectra: 'Iterable[BaseSpectrum]', **meta):
        """
        :param frequencies: frequency list
        :type frequencies: np.ndarray
        :param intensities: intensity list
        :type intensities: np.ndarray
        :param meta: metadata
        :type meta:
        """
        ...

    def __getitem__(self, item):
        """
        **LLM Docstring**

        Index into the collection of spectra: a scalar index returns the corresponding `BaseSpectrum` directly, while any other index (slice, array, etc.) returns a new `MultiSpectrum` wrapping the selected subset.

        :param item: the index or slice to apply
        :type item: int | slice | np.ndarray
        :return: the selected spectrum, or a `MultiSpectrum` wrapping the selected subset
        :rtype: BaseSpectrum | MultiSpectrum
        """
        ...

    def frequency_filter(self, freq_min, freq_max):
        """
        Filters by frequencies >= `freq_min` and <= `freq_max`

        :param freq_min: min frequency
        :type freq_min: float
        :param freq_max: max frequency
        :type freq_max: float
        :return: subspectrum
        :rtype: MultiSpectrum
        """
        ...

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
        ...

    def plot(self, figure=None, **opts):
        """
        A just plots all the spectra on the same figure

        :param opts: plotting options to be fed through to whatever the plotting function uses
        :type opts:
        :return:
        :rtype:
        """
        ...