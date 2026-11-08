"""
Provides a general base spectrum class that can be extended to new fancy spectral forms
"""

import numpy as np, json, abc
import McUtils.Plots as plt
from McUtils.Data import UnitsData

__all__ = [
    "BaseSpectrum",
    "DiscreteSpectrum",
    "ContinuousSpectrum",
    "BroadenedSpectrum",
    "MultiSpectrum"
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
        """
        **LLM Docstring**

        Index into the spectrum: a scalar index returns the corresponding `(frequency, intensity)` pair, while any other index (slice, array, etc.) returns a subspectrum via `take_subspectrum`.

        :param item: the index or slice to apply
        :type item: int | slice | np.ndarray
        :return: a `(frequency, intensity)` pair for a scalar index, or a subspectrum otherwise
        :rtype: tuple[float, float] | BaseSpectrum
        """

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

        oscillator_strengths = np.linalg.norm(transition_moments, axis=-1) ** 2
        units = UnitsData.convert("OscillatorStrength", "KilometersPerMole")
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        meta.update(frequency_unit="Wavenumbers", intensity_unit="KilometersPerMole")
        return cls(
            frequencies * h2w,
            units * frequencies * oscillator_strengths,
            **meta
        )

    @staticmethod
    def _raman_invariants(polarizabilities):
        """Return squared isotropic and anisotropic Placzek invariants.

        The nonresonant electric-dipole approximation uses symmetric 3 x 3
        tensors. Complex symmetric tensors are supported using absolute squares.
        """
        pol = np.asarray(polarizabilities)
        if pol.ndim < 2 or pol.shape[-2:] != (3, 3):
            raise ValueError("Raman polarizabilities must have trailing shape (3, 3)")
        if not np.all(np.isfinite(pol)):
            raise ValueError("Raman polarizabilities must be finite")
        if not np.allclose(pol, np.swapaxes(pol, -2, -1)):
            raise ValueError("the Placzek approximation requires symmetric polarizabilities")
        diagonal = np.diagonal(pol, axis1=-2, axis2=-1)
        isotropic = np.abs(np.sum(diagonal, axis=-1) / 3) ** 2
        anisotropic = (
            np.abs(diagonal[..., 0] - diagonal[..., 1]) ** 2
            + np.abs(diagonal[..., 1] - diagonal[..., 2]) ** 2
            + np.abs(diagonal[..., 2] - diagonal[..., 0]) ** 2
            + 6 * (np.abs(pol[..., 0, 1]) ** 2
                   + np.abs(pol[..., 0, 2]) ** 2
                   + np.abs(pol[..., 1, 2]) ** 2)
        ) / 2
        return isotropic, anisotropic

    @classmethod
    def from_raman_moments(cls, frequencies, transition_polarizabilities, pump_frequency=0, **meta):
        """Build a nonresonant Raman spectrum from transition tensors in a.u.

        Frequencies and the optional pump photon energy are in Hartrees. Without
        a pump (None or zero), return transition strengths 45|a|^2 + 7 gamma^2
        in Angstroms**6. These include vibrational matrix elements; they are not
        the derivative-based harmonic activities in Angstroms**4/amu.

        A positive pump returns the differential *energy* scattering cross
        section in cm^2/sr: (2 pi nu_s)^4 (45|a|^2 + 7 gamma^2) / 45, with the
        polarizability volume and scattered wavenumber expressed in cm. The
        geometry is 90-degree detection, incident linear polarization normal to
        the scattering plane, and both scattered polarizations collected.
        Positive shifts are Stokes (pump minus shift); negative shifts describe
        anti-Stokes transitions. Closed channels have zero cross section. No
        thermal populations or frequency dependence of the input tensor is added.
        """
        frequencies = np.asarray(frequencies)
        if np.iscomplexobj(frequencies) or not np.all(np.isfinite(frequencies)):
            raise ValueError("Raman transition frequencies must be finite and real")
        isotropic, anisotropic = cls._raman_invariants(transition_polarizabilities)
        strengths = 45 * isotropic + 7 * anisotropic
        if strengths.shape != frequencies.shape:
            raise ValueError("each Raman transition frequency needs one polarizability tensor")
        if pump_frequency is not None and (
                not np.isscalar(pump_frequency) or np.iscomplexobj(pump_frequency)
                or not np.isfinite(pump_frequency) or pump_frequency < 0):
            raise ValueError("pump_frequency must be a finite nonnegative photon energy in Hartrees")
        if pump_frequency is None or pump_frequency == 0:
            intensities = strengths * UnitsData.convert("BohrRadius", "Angstroms") ** 6
            unit = "Angstroms**6"
            label = "Raman transition strength"
        else:
            scattered = np.maximum(pump_frequency - frequencies, 0)
            wavevectors = 2 * np.pi * scattered * UnitsData.convert("Hartrees", "Wavenumbers")
            intensities = wavevectors ** 4 * strengths * UnitsData.convert("BohrRadius", "Centimeters") ** 6 / 45
            unit = "Centimeters**2/Steradians"
            label = "Raman differential cross section"
            meta.update(pump_frequency=float(pump_frequency), pump_frequency_unit="Hartrees",
                        scattering_geometry="90 degrees; incident polarization perpendicular to scattering plane; unanalysed detection")
        meta.update(frequency_unit="Wavenumbers", intensity_unit=unit, intensity_label=label)
        return cls(frequencies * UnitsData.convert("Hartrees", "Wavenumbers"), intensities, **meta)

    def normalize(self, which=None):
        """
        **LLM Docstring**

        Build a copy of the spectrum with intensities normalized by either the overall maximum intensity or a specific reference intensity.

        :param which: the index of a specific transition to normalize against; if `None`, normalizes by the maximum intensity
        :type which: int | None
        :return: the normalized spectrum
        :rtype: DiscreteSpectrum
        """
        meta = dict(self.meta, intensity_unit=None, intensity_label="Relative intensity")
        return type(self)(
            self.frequencies,
            self.intensities / (np.max(self.intensities) if which is None else self.intensities[which]),
            **meta
        )

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

    def broaden(self, breadth=10, *, broadening_type="gaussian", noising_function=None,):
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
                                 noising_function=noising_function,
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

    def __init__(self, frequencies, intensities, broadening_type="gaussian", breadth=10,
                 noising_function=None,
                 **meta):
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
        self.noising_function = noising_function

    def _eval_gauss_broadening(self, pts, height, center, breadth, target_zero=.1,
                               renormalize=True, adjust_width=True):
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
        if renormalize:
            height = height * (np.sqrt(2 * np.pi) * breadth)  # remove the normalization
        return height * bd.pdf(pts)

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
        from scipy.stats import cauchy

        bd = cauchy(loc=center, scale=breadth).pdf
        # if height < 1:
        #     height = 1
        if renormalize:
            height = height * (np.pi * breadth)  # remove the normalization
        return height * bd(pts)

    def _get_pts(self, step_size=.5, freq_min=None, freq_max=None, stddevs=5,
                 adjust_width=True,
                 renormalize=True
                 ):
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
        vals = np.sum([
            bt(freq_vals, i, c, b, adjust_width=adjust_width, renormalize=renormalize)
            for i, c, b in zip(ints, freqs, breadths)
        ], axis=0)

        if self.noising_function is not None:
            noise = self.noising_function(freq_vals, vals)
            vals = noise + vals

        return freq_vals, vals

    def plot(self, step_size=.5, freq_min=None, freq_max=None, figure=None, plot_style=None, filled=False,
             adjust_width=True,
             renormalize=True,
             **opts
             ):
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

        freqs, ints = self._get_pts(step_size=step_size, freq_min=freq_min, freq_max=freq_max,
                                    adjust_width=adjust_width,
                                    renormalize=renormalize
                                    )
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

class MultiSpectrum:
    """
    A wrapper for multiple spectra, really just for the plotting convenience
    """
    """
        Base class to support spectral operation
        """

    def __init__(self, spectra:"Iterable[BaseSpectrum]", **meta):
        """
        :param frequencies: frequency list
        :type frequencies: np.ndarray
        :param intensities: intensity list
        :type intensities: np.ndarray
        :param meta: metadata
        :type meta:
        """
        self.specs = spectra
        self.meta = meta

    def __getitem__(self, item):
        """
        **LLM Docstring**

        Index into the collection of spectra: a scalar index returns the corresponding `BaseSpectrum` directly, while any other index (slice, array, etc.) returns a new `MultiSpectrum` wrapping the selected subset.

        :param item: the index or slice to apply
        :type item: int | slice | np.ndarray
        :return: the selected spectrum, or a `MultiSpectrum` wrapping the selected subset
        :rtype: BaseSpectrum | MultiSpectrum
        """

        specs = self.specs[item]
        if isinstance(specs, BaseSpectrum):
            return specs
        else:
            return type(self)(specs)

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

        return type(self)([
            s.frequency_filter(freq_min, freq_max)
            for s in self.specs
        ])

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

        return type(self)([
            s.intensity_filter(int_min, int_max)
            for s in self.specs
        ])

    def plot(self, figure=None, **opts):
        """
        A just plots all the spectra on the same figure

        :param opts: plotting options to be fed through to whatever the plotting function uses
        :type opts:
        :return:
        :rtype:
        """

        for p in self.specs:
            subfig = p.plot(figure=figure, **opts)
            if figure is None:
                figure = subfig
        return figure