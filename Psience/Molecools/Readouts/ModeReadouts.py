"""
Readout interface for normal modes. Usually reached through `MoleculeReadoutInterface`'s
``normal_modes`` section, but it can be used on its own:
``NormalModesReadoutInterface(modes, molecule=mol).to_readout()``.
"""

import numpy as np

from McUtils.Jupyter import (
    ReadoutAdapter, readout_section, ReadoutSectionUnavailable,
    ReadoutTable, ReadoutScene, ReadoutGallery, ReadoutArray,
    TabularData, Column, ArrayData,
)

__all__ = [
    "NormalModesReadoutInterface",
]


class NormalModesReadoutInterface(ReadoutAdapter):
    """
    **LLM Docstring**

    Readout sections for a set of normal modes: a frequency (and, when dipole derivatives are
    available, intensity) table, a gallery of mode animations, and the Cartesian displacement
    vectors (exported, not displayed). Frequencies are converted from Hartrees to the readout's
    frequency unit while the readout is built.
    """
    readout_id = "normal_modes"
    readout_title = "Normal modes"

    def __init__(self, modes, molecule=None):
        super().__init__(modes)
        self.molecule = molecule

    @property
    def modes(self):
        return self.obj

    def get_molecule(self, ctx=None):
        mol = self.molecule
        if mol is None and ctx is not None:
            from ..Molecule import Molecule
            mol = ctx.find_parent(Molecule)
        return mol

    def _which(self, which):
        n = len(self.modes.freqs)
        if which is None:
            return list(range(n))
        if isinstance(which, (int, np.integer)):
            which = [which]
        return [int(w) % n for w in which]

    def _intensities(self):
        mol = self.molecule
        if mol is None or mol.dipole_derivatives is None:
            return None
        try:
            spec = mol.get_harmonic_spectrum()
        except Exception:
            return None
        ints = np.asarray(spec.intensities, dtype=float)
        return ints if len(ints) == len(self.modes.freqs) else None

    @readout_section("frequencies", title="Frequencies")
    def readout_frequencies(self, ctx, which=None, intensities=True):
        """Harmonic frequencies (and IR intensities, in km/mol, when available)."""
        idx = self._which(which)
        freqs = np.asarray(self.modes.freqs, dtype=float)[idx]
        cols = [
            Column("mode", np.array(idx) + 1, label="Mode", quantity="index"),
            ctx.column("frequency", freqs, quantity="frequency", unit="Hartrees", label="Frequency"),
        ]
        if intensities:
            ints = self._intensities()
            if ints is not None:
                cols.append(ctx.column("intensity", ints[idx], quantity="intensity", unit="KilometersPerMole",
                                       label="Intensity"))
        return ReadoutTable(TabularData(cols, name="table"))

    @readout_section("animations", title="Mode animations", available=lambda self: self.molecule is not None)
    def readout_animations(self, ctx, which=None, steps=12, extent=.5, duration=2.0, per_slide=None,
                           columns=None, shared_view=True, html_backend="x3d", animate_options=None):
        """
        One animated scene per mode, captioned with the mode's frequency. Scenes are built lazily
        with `Molecule.animate_mode`: ``backend='mesh3D'`` for PowerPoint (`.glb` + poster) and, with
        ``html_backend='x3d'``, the much lighter native X3D animation for HTML (``'mesh'`` reuses
        the mesh animation instead). ``shared_view`` draws every mode with one camera and scale.
        """
        mol = self.get_molecule(ctx)
        if mol is None:
            raise ReadoutSectionUnavailable("animations need the molecule the modes belong to")
        idx = self._which(which)
        freqs, freq_unit = ctx.convert(np.asarray(self.modes.freqs, dtype=float), "frequency", "Hartrees")
        animate_options = dict(animate_options or {})
        scenes = []
        for i in idx:
            def model(i=i):
                return mol.animate_mode(i, backend="mesh3D", steps=steps, extent=extent,
                                        animation_options={"animation_duration": duration}, **animate_options)
            html = None
            if html_backend == "x3d":
                def html(view=None, i=i):
                    opts = dict(animate_options)
                    if view is not None:
                        opts.setdefault("view_settings", view.x3d_viewpoint())
                    return mol.animate_mode(i, backend="x3d", steps=steps, extent=extent,
                                            animation_options={"animation_duration": duration}, **opts)
            caption = ctx.field("frequency", float(freqs[i]), quantity="frequency", label=f"Mode {i + 1}")
            caption.unit = freq_unit
            scenes.append(ReadoutScene(model=model, html=html, caption=caption, id=f"mode_{i + 1}", animated=True,
                                       meta={"mode": i, "steps": steps, "extent": extent, "duration": duration}))
        return ReadoutGallery(*scenes, id="gallery", shared_view=shared_view, per_slide=per_slide, columns=columns)

    @readout_section("displacements", title="Displacement vectors")
    def readout_displacements(self, ctx, which=None, display=False):
        """
        Normalized Cartesian displacement vectors, shape ``(mode, atom, xyz)``; exported with the
        data but not displayed unless ``display=True``.
        """
        idx = self._which(which)
        modes = self.modes.remove_mass_weighting()
        disp = np.asarray(modes.coords_by_modes, dtype=float)[idx]
        disp = disp / np.linalg.norm(disp, axis=1, keepdims=True)
        nat = disp.shape[1] // 3
        mol = self.molecule
        atoms = list(mol.atoms) if mol is not None else [str(a + 1) for a in range(nat)]
        arr = ArrayData("vectors", disp.reshape(len(idx), nat, 3), axes=("mode", "atom", "xyz"),
                        label="Normalized displacement vectors",
                        description="Cartesian displacements per mode (mass weighting removed), unit norm")
        return ReadoutArray(arr, display=display)
