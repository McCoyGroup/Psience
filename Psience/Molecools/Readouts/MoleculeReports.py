"""
The report template for `Molecule`: a fixed outline (structure, then vibrations), so a report on one molecule needs no
model to choose or order content. `MoleculeReadoutInterface.report_template` points here (imported lazily).
"""

from McUtils.Readouts.Reports import ReportTemplate

__all__ = ["MoleculeReportTemplate"]


class MoleculeReportTemplate(ReportTemplate):
    """
    :cvar sections: structure (``overview``, ``structure``) then vibrations (the frequency table only; animations and
        displacement vectors are separate sections a model or user can add)
    :cvar optional: sections that suit a longer report; not used by `ReportHarness` yet, listed for tools and callers
    """
    name = "molecule"
    sections = [
        ("Structure", ["overview", "structure"]),
        ("Vibrations", ["normal_modes.frequencies"]),
    ]
    optional = ["source", "rdkit", "cartesians", "internals", "bonds",
                "normal_modes.animations", "normal_modes.displacements",
                "spectrum.lines", "spectrum.stick", "spectrum.levels", "spectrum.broadened"]

    _captions = {
        "overview": "Formula, size, charge and energy of the molecule.",
        "structure": "The equilibrium geometry.",
        "cartesians": "Atomic masses and Cartesian coordinates.",
        "internals": "Bond lengths and angles of the molecule.",
        "bonds": "The bonds of the molecule with their orders and lengths.",
        "normal_modes.frequencies": "Harmonic vibrational frequencies of the normal modes.",
        "spectrum.lines": "The lines of the harmonic spectrum: frequency and intensity.",
        "spectrum.stick": "The harmonic IR stick spectrum: intensity against frequency.",
        "spectrum.levels": "The harmonic frequencies drawn as a ladder of levels.",
        "spectrum.broadened": "The harmonic IR spectrum broadened onto a frequency grid.",
    }

    def caption_for(self, handle, section, source):
        return self._captions.get(section)
