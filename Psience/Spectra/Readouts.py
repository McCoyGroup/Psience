"""
Readout interface for `DiscreteSpectrum` (and so for `HarmonicSpectrum`), built on `McUtils.Readouts`.

``DiscreteSpectrum(...).to_readout()`` (or ``to_readout(spectrum)``) gives the sections::

    lines       a table of the lines: number, frequency, intensity
    stick       the stick spectrum as a figure (its numbers are the ``lines`` table)
    levels      the frequencies as a level ladder, one line per frequency (off by default)
    broadened   the spectrum broadened onto a grid, as a line plot (off by default; a long series)

A `Molecule` reaches this through its ``spectrum`` section (``mol.to_readout(include=["spectrum"])``), which is off by
default: it builds ``mol.get_harmonic_spectrum()`` and hands it to this interface, so the options below are given as
``spectrum={"stick": {...}}`` there.
"""

import numpy as np

from McUtils.Readouts import (ReadoutAdapter, readout_section, register_readout, ReadoutTable, TabularData, Column,
                              ReadoutSectionUnavailable)
from McUtils.Readouts.Builders import stick_chart, line_chart, level_chart
from .BaseSpectrum import DiscreteSpectrum

__all__ = ["DiscreteSpectrumReadoutInterface"]

_KINDS = {"ir": "IR", "raman": "Raman"}


@register_readout(DiscreteSpectrum)
class DiscreteSpectrumReadoutInterface(ReadoutAdapter):
    """
    **LLM Docstring**

    Readout sections for a discrete spectrum (frequencies in wavenumbers, intensities in km/mol unless the spectrum's
    ``meta`` gives ``intensity_unit``). ``meta["kind"]`` (``"ir"`` or ``"raman"``) only changes the captions. Values are
    converted to the readout's display units while the readout is built.
    """
    readout_id = "discrete_spectrum"
    readout_title = "Discrete spectrum"

    def __init__(self, spectrum, molecule=None):
        super().__init__(spectrum)
        self.molecule = molecule

    @property
    def spectrum(self):
        return self.obj

    @property
    def kind_label(self):
        kind = self.spectrum.meta.get("kind")
        return _KINDS.get(str(kind).lower(), str(kind)) if kind else None

    @property
    def intensity_unit(self):
        """
        The unit of the intensities: ``meta["intensity_unit"]`` when the spectrum states one (``None`` is a stated
        "no unit", as for a normalized spectrum), else km/mol (the IR convention).
        """
        return self.spectrum.meta.get("intensity_unit", "KilometersPerMole")

    @property
    def intensity_label(self):
        """``meta["intensity_label"]`` (``"Raman activity"``, ``"Relative intensity"``) or ``"Intensity"``."""
        return self.spectrum.meta.get("intensity_label") or "Intensity"

    @staticmethod
    def _fmt(values):
        """A format that keeps small or large intensities readable (a cross section is of order 1e-30)."""
        v = np.abs(np.asarray(values, dtype=float))
        v = v[np.isfinite(v) & (v > 0)]
        return "{:.4e}" if v.size and (v.max() < 1e-2 or v.max() > 1e6) else None

    def _name(self, what):
        k = self.kind_label
        return f"{k} {what}" if k else what.capitalize()

    def _lines(self, which=None, min_intensity=None):
        freqs = np.asarray(self.spectrum.frequencies, dtype=float)
        ints = np.asarray(self.spectrum.intensities, dtype=float)
        idx = np.arange(len(freqs))
        if which is not None:
            idx = np.array([int(w) % len(freqs) for w in ([which] if isinstance(which, (int, np.integer)) else which)])
        if min_intensity is not None and ints.size:
            idx = idx[ints[idx] >= float(min_intensity) * np.max(ints)]
        return idx, freqs[idx], ints[idx]

    @readout_section("lines", title="Lines", description="The spectrum's lines: frequency and intensity.")
    def readout_lines(self, ctx, which=None, min_intensity=None):
        """
        A row per line. ``which`` picks lines by (0-based) index; ``min_intensity`` keeps lines at least that fraction
        of the strongest line.
        """
        idx, f, i = self._lines(which, min_intensity)
        cols = [Column("line", idx + 1, label="Line", quantity="index"),
                ctx.column("frequency", f, quantity="frequency", unit="Wavenumbers", label="Frequency"),
                ctx.column("intensity", i, quantity="intensity", unit=self.intensity_unit, label=self.intensity_label,
                           fmt=self._fmt(i))]
        return ReadoutTable(TabularData(cols, name="table", title=self._name("lines")))

    @readout_section("stick", title="Stick spectrum", description="The stick spectrum: intensity against frequency.")
    def readout_stick(self, ctx, which=None, min_intensity=None, caption=None):
        """The lines as sticks. The plotted numbers are the ``frequency`` and ``intensity`` columns of its data."""
        _, f, i = self._lines(which, min_intensity)
        f, f_unit = ctx.convert(f, "frequency", "Wavenumbers")
        i, i_unit = ctx.convert(i, "intensity", self.intensity_unit) if self.intensity_unit else (i, None)
        return stick_chart(f, i, x_label="Frequency", y_label=self.intensity_label, x_unit=f_unit, y_unit=i_unit,
                           x_name="frequency", y_name="intensity", x_quantity="frequency", y_quantity="intensity",
                           caption=caption or f"{self._name('stick spectrum')}: intensity against frequency.",
                           id="stick", title=self._name("stick spectrum"))

    @readout_section("levels", title="Level ladder", default=False,
                     description="The frequencies drawn as a ladder of levels, one line per frequency.")
    def readout_levels(self, ctx, which=None, min_intensity=None, caption=None):
        """
        The frequencies as a ladder: a horizontal line at each frequency, numbered like the ``lines`` table. Intensities
        are not drawn; the data is the ``level`` number and the ``frequency``.
        """
        idx, f, _ = self._lines(which, min_intensity)
        f, f_unit = ctx.convert(f, "frequency", "Wavenumbers")
        return level_chart(f, labels=[str(i + 1) for i in idx], y_label="Frequency", y_unit=f_unit,
                           value_name="frequency",
                           caption=caption or f"{self._name('frequencies')} as a ladder of levels.",
                           id="levels", title="Level ladder")

    @readout_section("broadened", title="Broadened spectrum", default=False,
                     description="The spectrum broadened onto a frequency grid (a long series of points).")
    def readout_broadened(self, ctx, breadth=10, broadening_type="gaussian", step_size=1.0, freq_min=None,
                          freq_max=None, caption=None):
        """
        The lines broadened (``gaussian`` or ``lorentzian``, ``breadth`` in wavenumbers) and evaluated every
        ``step_size`` wavenumbers: one point per step, so hundreds to thousands of rows in its data.
        """
        from .BaseSpectrum import BroadenedSpectrum
        sp = self.spectrum
        if not len(sp.frequencies):
            raise ReadoutSectionUnavailable("the spectrum has no lines to broaden")
        b = BroadenedSpectrum(sp.frequencies, sp.intensities, broadening_type=broadening_type, breadth=breadth)
        x, y = b._get_pts(step_size=step_size, freq_min=freq_min, freq_max=freq_max)
        x, x_unit = ctx.convert(np.asarray(x, dtype=float), "frequency", "Wavenumbers")
        return line_chart(x, np.asarray(y, dtype=float), x_label="Frequency", y_label=self.intensity_label, x_unit=x_unit,
                          x_name="frequency", names=["intensity"], x_quantity="frequency",
                          caption=caption or f"{self._name('spectrum')} broadened by {broadening_type} peaks: "
                                             f"intensity against frequency.",
                          id="broadened", title=self._name("broadened spectrum"))
