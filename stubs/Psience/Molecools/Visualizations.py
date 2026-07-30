from __future__ import annotations
import itertools
import numpy as np
from McUtils.Coordinerds import CoordinateSet, CartesianCoordinates3D
from McUtils.Data import UnitsData, AtomData
import McUtils.Numputils as nput
import McUtils.Devutils as dev
import McUtils.Plots as plt
__all__ = ['MoleculePlotter', 'Graphics3DMoleculePlotter', 'JSMolMoleculePlotter', 'RDKitMoleculePlotter', 'JupyterMoleculePlotter']

class MoleculePlotter:
    """
    Base strategy: holds a molecule + resolved options, owns the shared plotting
    configuration and every data-derivation method, and dispatches to the right
    concrete plotter for the resolved mode.
    """
    highlight_styles = {'glow': 'green', 'color': 'white'}
    vector_style = {'color': 'black', 'radius': 0.1}
    principle_axes_style = [{'color': 'green'}, {'color': 'red'}, {'color': 'blue'}]
    draw_coords_style = {'line_color': 'black', 'line_thickness': 0.05}
    draw_coords_label_style = {'color': 'black', 'billboard': True}
    backend_options_resolution = {'mode_vectors': 'x3d', 'dipole': 'x3d', 'include_jsmol_script_interface': 'jsmol'}

    @staticmethod
    def _flat_color(i, a, styles):
        """
        **LLM Docstring**

        Style modifier used by the `'flat'` x3d subtheme: strips out shading by forcing the color to black and darkening the glow slightly, if a color was set.

        :param i: the atom/bond index (or index pair) the style applies to; unused in this modifier's body
        :type i: int | tuple
        :param a: the associated atom label(s); unused in this modifier's body
        :type a: str | tuple
        :param styles: the style dict to modify
        :type styles: dict
        :return: the (possibly updated) style dict
        :rtype: dict
        """
        ...
    subthemes = {}
    _mode_registry: dict = {}
    modes: tuple = ()

    def __init_subclass__(cls, **kwargs):
        """
        **LLM Docstring**

        Registers each concrete plotter subclass under every mode name it declares in its `modes` tuple, and merges the subclass's `subthemes` into the shared `plot_themes` dict for those modes.

        :param kwargs: forwarded to `super().__init_subclass__`
        :type kwargs: dict
        :return: None
        :rtype: None
        """
        ...

    def __init__(self, mol, geometries, *, mode, backend, atoms=None, masses=None, coords=None, bonds=None, **opts):
        """
        **LLM Docstring**

        Set up a plotter instance for a molecule (or molecule-like data) over one or more geometries: resolves the atom list, per-atom masses (converting to atomic units if given in amu), and stores the coordinates/bonds overrides and any extra plotting options.

        :param mol: the molecule being plotted
        :type mol: AbstractMolecule
        :param geometries: the geometry (or geometries) to render
        :type geometries: CoordinateSet | np.ndarray
        :param mode: the resolved display mode (e.g. `'x3d'`, `'jsmol'`, `'svg2d'`)
        :type mode: str
        :param backend: the resolved rendering backend
        :type backend: str
        :param atoms: atom labels/symbols to use instead of `mol.atoms`
        :type atoms: tuple[str] | None
        :param masses: atomic masses to use instead of those looked up from `atoms`
        :type masses: np.ndarray | None
        :param coords: coordinates to use instead of `mol.coords` (lazily resolved via the `coords` property)
        :type coords: CoordinateSet | None
        :param bonds: bonds to use instead of `mol.bonds` (lazily resolved via the `bonds` property)
        :type bonds: tuple[tuple] | None
        :param opts: additional plotting options stored for later use by `plot`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def coords(self):
        """
        **LLM Docstring**

        The coordinates to plot, lazily falling back to `self.mol.coords` if none were supplied explicitly.

        :return: the coordinates to render
        :rtype: CoordinateSet
        """
        ...

    @property
    def bonds(self):
        """
        **LLM Docstring**

        The bonds to plot, lazily falling back to `self.mol.bonds` if none were supplied explicitly.

        :return: the bonds to render
        :rtype: tuple[tuple]
        """
        ...

    @classmethod
    def _resolve_plot_mode(cls, mol, mode, backend, geometries, **ignored):
        """
        **LLM Docstring**

        Work out the display `mode` and rendering `backend` to use from whatever combination of `mode`/`backend`/other hints (an explicit backend alias, the presence of multiple geometries, or which backend-specific options were passed) was given, falling back to the molecule's own `display_mode` as a last resort.

        :param mol: the molecule being plotted, used for its default `display_mode`
        :type mol: AbstractMolecule
        :param mode: an explicitly requested display mode, if any
        :type mode: str | None
        :param backend: an explicitly requested backend, if any
        :type backend: str | None
        :param geometries: the geometries to be plotted, used to decide on a default backend if none is given
        :type geometries: tuple
        :param ignored: extra plotting options inspected only to detect backend-implying keys (via `backend_options_resolution`)
        :type ignored: dict
        :return: the resolved `(mode, backend)` pair
        :rtype: tuple[str, str]
        """
        ...

    @classmethod
    def _resolve_plot_theme(cls, mode, backend, theme, base_opts):
        """
        **LLM Docstring**

        Merge the appropriate theme dict (looked up by `(backend, mode)`, then `backend`, then `mode`) with the `'default'` theme, and layer it under `base_opts` so any option the caller didn't already specify gets the theme's value (dict-valued options are merged rather than replaced).

        :param mode: the display mode
        :type mode: str
        :param backend: the rendering backend
        :type backend: str
        :param theme: the theme name to select within the resolved theme set (e.g. `'default'`, `'simple'`, `'flat'`)
        :type theme: str
        :param base_opts: the options dict to fill in with theme defaults, modified in place
        :type base_opts: dict
        :return: `base_opts`, updated with theme defaults
        :rtype: dict
        """
        ...

    @classmethod
    def _default_plot_range(cls, geometries, pr, plot_range_padding, radii):
        """
        **LLM Docstring**

        Compute a default `[[xmin,xmax],[ymin,ymax],[zmin,zmax]]` plot range spanning all the geometries (if `pr` is not already given), then optionally pad it by `plot_range_padding` (a single number, `'auto'` meaning the largest atomic radius, or per-axis values/pairs).

        :param geometries: the geometry (or geometries) being plotted, used to compute the default bounding range
        :type geometries: np.ndarray
        :param pr: an explicit plot range to use instead of computing one from `geometries`
        :type pr: list | None
        :param plot_range_padding: extra padding to add around the range; `None` for no padding, `'auto'` to pad by the largest radius in `radii`, or a number/per-axis spec
        :type plot_range_padding: float | str | tuple | None
        :param radii: atomic radii, used to compute the padding amount when `plot_range_padding` is `'auto'`
        :type radii: np.ndarray
        :return: the (possibly padded) plot range
        :rtype: list
        """
        ...

    def _set_backend_figure_options(self, figure, mode, backend, include_save_buttons=None, dynamic_loading=None, recording_options=None, **ignored):
        """
        **LLM Docstring**

        Apply backend-specific figure-level options (save/record/view-settings buttons, dynamic loading, recording options for `x3d`; save buttons for `plotly3D`) to an already-built figure object, in place.

        :param figure: the figure object to configure
        :type figure: object
        :param mode: the display mode
        :type mode: str
        :param backend: the rendering backend
        :type backend: str
        :param include_save_buttons: whether to show export/record/view-settings buttons
        :type include_save_buttons: bool | None
        :param dynamic_loading: whether the figure should use dynamic (lazy) loading, `x3d` only
        :type dynamic_loading: bool | None
        :param recording_options: recording configuration to attach to the figure, `x3d` only
        :type recording_options: dict | None
        :param ignored: any other options, accepted but not used
        :type ignored: dict
        :return: None (the figure is modified in place)
        :rtype: None
        """
        ...

    @staticmethod
    def _flat_color(i, a, styles):
        """
        **LLM Docstring**

        Style modifier used by the `'flat'` x3d subtheme: strips out shading by forcing the color to black and darkening the glow slightly, if a color was set.

        :param i: the atom/bond index (or index pair) the style applies to; unused in this modifier's body
        :type i: int | tuple
        :param a: the associated atom label(s); unused in this modifier's body
        :type a: str | tuple
        :param styles: the style dict to modify
        :type styles: dict
        :return: the (possibly updated) style dict
        :rtype: dict
        """
        ...
    atom_color_updates = {}

    @classmethod
    def _resolve_atom_color(cls, atom_data):
        """
        **LLM Docstring**

        Look up the display color for an atom, preferring any override in `atom_color_updates` (keyed by element symbol) over the atom's default `IconColor`.

        :param atom_data: the atom-data record (as returned by `AtomData`) to color
        :type atom_data: dict
        :return: the resolved color
        :rtype: str
        """
        ...
    material_backends = {'mesh3D', 'x3d'}

    def _prep_display_atom_style(self, atom_style, highlight_atoms, *, backend, reflectiveness, highlight_styles):
        """
        **LLM Docstring**

        Normalize the `atom_style` specification (a bool, a dict keyed by index/element/string option, or a sequence) into a per-atom-index style dict, layering in reflectiveness options for the `x3d` backend and merging in `highlight_styles` for any atoms in `highlight_atoms`; also derives the resolved per-atom `colors`/`glows` lists (falling back to each atom's default color when not overridden).

        :param atom_style: the raw atom-style specification
        :type atom_style: dict | bool | Iterable | None
        :param highlight_atoms: atom indices to additionally style with `highlight_styles`
        :type highlight_atoms: Iterable[int] | None
        :param backend: the rendering backend, used to decide whether to add reflectiveness options
        :type backend: str
        :param reflectiveness: reflectiveness value (0-1-ish) to translate into `specularity`/`shininess` options on `x3d`
        :type reflectiveness: float | None
        :param highlight_styles: the style overrides to apply to highlighted atoms
        :type highlight_styles: dict
        :return: `(atom_style, highlight_atoms, colors, glows)` -- the normalized per-atom style dict, the (possibly unchanged) highlight list, and the resolved per-atom colors/glows
        :rtype: tuple[dict, Iterable | None, list, list]
        """
        ...

    def _prep_display_atom_text(self, atom_text, display_atom_numbers, label_style):
        """
        **LLM Docstring**

        Build the per-atom text-label specification: if `display_atom_numbers` is requested, builds a label dict per selected atom (merging in `label_style` and any per-atom overrides); otherwise passes `atom_text` through (defaulting to a list of `None`s).

        :param atom_text: an existing per-atom text specification, if provided directly
        :type atom_text: list | None
        :param display_atom_numbers: `True` to label every atom with its index, or an iterable/dict of atom indices (optionally mapped to per-atom style overrides) to label a subset
        :type display_atom_numbers: bool | Iterable[int] | dict
        :param label_style: base style applied to every generated atom-number label
        :type label_style: dict
        :return: a list, one entry per atom, of either `None` or a label dict (with at least a `'text'` key)
        :rtype: list
        """
        ...

    def _prep_display_bond_style(self, bond_style, highlight_bonds, *, backend, reflectiveness, highlight_atoms, highlight_styles, capped_bonds):
        """
        **LLM Docstring**

        Normalize the `bond_style` specification (a bool, a dict keyed by index-pair/element/string option, or a sequence) into a per-bond-pair style dict covering every atom pair, layering in `capped_bonds`/reflectiveness options and merging in `highlight_styles` for any bonds in `highlight_bonds` (defaulting to `highlight_atoms` if `highlight_bonds` isn't given).

        :param bond_style: the raw bond-style specification
        :type bond_style: dict | bool | Iterable | None
        :param highlight_bonds: bonds (as index pairs) to additionally style with `highlight_styles`
        :type highlight_bonds: Iterable | None
        :param backend: the rendering backend, used to decide whether to add reflectiveness options
        :type backend: str
        :param reflectiveness: reflectiveness value to translate into `specularity`/`shininess` options on `x3d`
        :type reflectiveness: float | None
        :param highlight_atoms: atom indices used as a fallback set of bonds to highlight if `highlight_bonds` is `None`
        :type highlight_atoms: Iterable[int] | None
        :param highlight_styles: the style overrides to apply to highlighted bonds
        :type highlight_styles: dict
        :param capped_bonds: whether bonds should be drawn with capped ends
        :type capped_bonds: bool
        :return: `(bond_style, highlight_bonds)` -- the normalized per-bond-pair style dict and the (possibly defaulted) highlight-bonds list
        :rtype: tuple[dict, Iterable | None]
        """
        ...

    def _get_atomic_radius(self, atom_data, radius_type=None):
        """
        **LLM Docstring**

        Resolve the display radius for a single atom, delegating to `self.mol._get_atomic_radius` if a molecule is attached, otherwise falling back to the atom's icon radius (or van der Waals radius if the icon radius is too small) or an explicitly named radius field.

        :param atom_data: the atom-data record to get a radius for
        :type atom_data: dict
        :param radius_type: the specific radius field to use (e.g. `"VanDerWaalsRadius"`); if `None`, uses the icon radius with a van der Waals fallback
        :type radius_type: str | None
        :return: the resolved atomic radius
        :rtype: float
        """
        ...

    def _get_atom_radii(self, atom_radii, atom_radius_scaling, radius_type):
        """
        **LLM Docstring**

        Resolve the final per-atom display radii by combining an `atom_radii` override (a single number, a dict keyed by index/element symbol, or `None`) with each atom's default radius (via `_get_atomic_radius`), then scaling by `atom_radius_scaling`.

        :param atom_radii: radius override(s); `None` to use each atom's default, a number for a uniform override, or a dict keyed by atom index or element symbol
        :type atom_radii: float | dict | None
        :param atom_radius_scaling: a uniform or per-atom scale factor applied to the resolved radii
        :type atom_radius_scaling: float | list
        :param radius_type: the radius field to use when falling back to `_get_atomic_radius`
        :type radius_type: str | None
        :return: the final per-atom radii
        :rtype: list[float]
        """
        ...

    def _prep_display_dipole(self, geometries, dipole, dipole_origin, units):
        """
        **LLM Docstring**

        Normalize the dipole vector and (optional) dipole origin for display: broadcasts a single vector/origin across all `geometries` if needed, and converts both from atomic units (Bohr) to the requested display `units`.

        :param geometries: the geometries being plotted, used to determine the broadcast length
        :type geometries: np.ndarray
        :param dipole: the dipole vector(s) to display, in atomic units
        :type dipole: np.ndarray | None
        :param dipole_origin: the origin point(s) for the dipole arrow, in atomic units
        :type dipole_origin: np.ndarray | None
        :param units: the display units to convert into (from `"BohrRadius"`); if `None`, no conversion is applied
        :type units: str | None
        :return: `(dipole, dipole_origin)`, normalized and unit-converted
        :rtype: tuple[np.ndarray | None, np.ndarray | None]
        """
        ...

    def _prep_principle_axes(self, geometries, units, principle_axes, principle_axes_origin, principle_axes_style):
        """
        **LLM Docstring**

        Resolve the principal-axes vectors to display: computes them from the moments of inertia if `principle_axes is True`, normalizes the per-axis style list, and broadcasts the axes/origin across all `geometries` as needed.

        :param geometries: the geometries being plotted, used for computing moments of inertia and for broadcasting
        :type geometries: np.ndarray
        :param units: the display units to convert the origin into (from `"BohrRadius"`)
        :type units: str | None
        :param principle_axes: `True` to compute the principal axes automatically, `False`/`None` to omit them, or explicit axis vectors
        :type principle_axes: bool | np.ndarray | None
        :param principle_axes_origin: the origin point(s) for the axis arrows
        :type principle_axes_origin: np.ndarray | None
        :param principle_axes_style: per-axis style override(s), merged onto the class-level `principle_axes_style` defaults
        :type principle_axes_style: dict | list | None
        :return: `(principle_axes, principle_axes_origin, principle_axes_style)`, all resolved/normalized
        :rtype: tuple[np.ndarray | None, np.ndarray | None, list | None]
        """
        ...

    def _prep_display_mode_vectors(self, geometries, units, mode_vectors, mode_vector_origins):
        """
        **LLM Docstring**

        Normalize the normal-mode displacement vectors and their origins for display: reshapes a flat vector into `(natoms, 3)`, broadcasts across all `geometries`, and converts to the requested display `units`.

        :param geometries: the geometries being plotted, used to determine the broadcast length
        :type geometries: np.ndarray
        :param units: the display units to convert into (from `"BohrRadius"`); if `None`, no conversion is applied
        :type units: str | None
        :param mode_vectors: the per-atom mode displacement vectors
        :type mode_vectors: np.ndarray | None
        :param mode_vector_origins: the per-atom origin points for the mode-vector arrows
        :type mode_vector_origins: np.ndarray | None
        :return: `(mode_vectors, mode_vector_origins)`, normalized and unit-converted
        :rtype: tuple[np.ndarray | None, np.ndarray | None]
        """
        ...

    def _prep_display_draw_coords(self, draw_coords, draw_coords_style):
        """
        **LLM Docstring**

        Normalize the `draw_coords` specification (an iterable of coordinate-index tuples, or an already-keyed style dict) into a dict keyed by coordinate tuple, and resolve the default `draw_coords_style` if none was given.

        :param draw_coords: the coordinates to draw as extra annotations (bonds/angles/dihedrals), given as an iterable of index tuples or a dict already mapping tuples to per-coordinate style overrides
        :type draw_coords: Iterable | dict | None
        :param draw_coords_style: the base style to apply to drawn coordinates; defaults to `self.draw_coords_style`
        :type draw_coords_style: dict | None
        :return: `(draw_coords, draw_coords_style)`, both normalized
        :rtype: tuple[dict | None, dict]
        """
        ...

    def _get_bond_primitives(self, geom, b, *, bond_list, bond_radius, radii, bond_center_radius_offset, multiple_bond_spacing, render_multiple_bonds, render_fractional_bonds, fractional_bond_offset, max_bond_orders, up_vector, cylinder_class, colors, glows, bond_style, theme_function, plotos, cylinder_options):
        """
        **LLM Docstring**

        Build the drawable primitives (one or more cylinder pairs, meeting at a midpoint) for a single bond, handling per-endpoint styling/coloring, radius-based endpoint offsetting, and -- for multiple/fractional bond orders -- the extra parallel cylinders offset perpendicular to the bond axis.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param b: the bond spec, an `(atom1, atom2[, order])` tuple
        :type b: tuple
        :param bond_list: the full list of bonds, used to find a reference atom for orienting multiple-bond offsets
        :type bond_list: list[tuple]
        :param bond_radius: the base cylinder radius for single bonds
        :type bond_radius: float
        :param radii: the per-atom display radii, used to offset bond endpoints away from atom centers
        :type radii: list[float]
        :param bond_center_radius_offset: how far to pull bond endpoints in from the atom surface; `'auto'`/a dict with `'padding'`/`'multi'` keys, a number, or `None`
        :type bond_center_radius_offset: str | dict | float | None
        :param multiple_bond_spacing: perpendicular spacing between parallel cylinders for double/triple bonds
        :type multiple_bond_spacing: float | None
        :param render_multiple_bonds: whether to draw double/triple bonds as multiple parallel cylinders rather than one
        :type render_multiple_bonds: bool
        :param render_fractional_bonds: whether to shorten cylinders proportionally for fractional bond orders
        :type render_fractional_bonds: bool
        :param fractional_bond_offset: scaling factor applied to the fractional-order shortening
        :type fractional_bond_offset: float
        :param max_bond_orders: optional mapping used to pad `bond_point_list` up to the number of cylinders implied by a maximum bond order for this pair
        :type max_bond_orders: dict | None
        :param up_vector: fallback reference vector for orienting multiple-bond offsets when no adjacent bond is found
        :type up_vector: Iterable[float] | None
        :param cylinder_class: the primitive class used to build each cylinder segment
        :type cylinder_class: type
        :param colors: per-atom colors
        :type colors: list
        :param glows: per-atom glow colors
        :type glows: list
        :param bond_style: per-bond/per-atom style overrides
        :type bond_style: dict
        :param theme_function: optional callable applied to each cylinder's resolved style before construction
        :type theme_function: callable | None
        :param plotos: base plot-wide options merged under everything else
        :type plotos: dict
        :param cylinder_options: theme-level cylinder options merged in
        :type cylinder_options: dict
        :return: the list of cylinder primitive objects making up this bond
        :rtype: list
        """
        ...

    def _get_bondlist_primitives(self, geom, bond_list, *, bond_radius, radii, bond_center_radius_offset, multiple_bond_spacing, render_multiple_bonds, render_fractional_bonds, fractional_bond_offset, max_bond_orders, up_vector, cylinder_class, colors, glows, bond_style, theme_function, plotos, cylinder_options):
        """
        **LLM Docstring**

        Build the drawable primitives for every bond in `bond_list` by calling `_get_bond_primitives` on each.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param bond_list: the bonds to draw
        :type bond_list: list[tuple]
        :param bond_radius: base bond cylinder radius
        :type bond_radius: float
        :param radii: per-atom display radii
        :type radii: list[float]
        :param bond_center_radius_offset: endpoint offset spec, forwarded to `_get_bond_primitives`
        :type bond_center_radius_offset: str | dict | float | None
        :param multiple_bond_spacing: spacing between parallel cylinders for multiple bonds
        :type multiple_bond_spacing: float | None
        :param render_multiple_bonds: whether to draw multiple parallel cylinders for double/triple bonds
        :type render_multiple_bonds: bool
        :param render_fractional_bonds: whether to shorten cylinders for fractional bond orders
        :type render_fractional_bonds: bool
        :param fractional_bond_offset: scaling factor for fractional-order shortening
        :type fractional_bond_offset: float
        :param max_bond_orders: optional mapping of maximum bond orders per pair
        :type max_bond_orders: dict | None
        :param up_vector: fallback orientation vector
        :type up_vector: Iterable[float] | None
        :param cylinder_class: the primitive class used for cylinder segments
        :type cylinder_class: type
        :param colors: per-atom colors
        :type colors: list
        :param glows: per-atom glow colors
        :type glows: list
        :param bond_style: per-bond/per-atom style overrides
        :type bond_style: dict
        :param theme_function: optional style-post-processing callable
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param cylinder_options: theme-level cylinder options
        :type cylinder_options: dict
        :return: a list, one entry per bond, each a list of cylinder primitives
        :rtype: list[list]
        """
        ...

    def _get_atom_primitives(self, geom, atoms, *, colors, radii, sphere_class, atom_style, theme_function, plotos, sphere_options):
        """
        **LLM Docstring**

        Build the drawable sphere primitives for a set of atoms, resolving each atom's final style (color, any per-atom `modifier` callable, theme function) before constructing it.

        :param geom: the atom coordinates to place spheres at
        :type geom: np.ndarray
        :param atoms: the atom indices being drawn (used only for iteration length via `zip`)
        :type atoms: Iterable
        :param colors: per-atom colors
        :type colors: list
        :param radii: per-atom display radii
        :type radii: list[float]
        :param sphere_class: the primitive class used to build each atom sphere
        :type sphere_class: type
        :param atom_style: per-atom style overrides, keyed by atom index
        :type atom_style: dict
        :param theme_function: optional callable applied to each sphere's resolved style before construction
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param sphere_options: theme-level sphere options
        :type sphere_options: dict
        :return: the list of sphere primitive objects, one per atom
        :rtype: list
        """
        ...

    def _get_dipole_primitives(self, geom, dip, *, dipole_origin, dipole_origin_mode, mode_vector_display_cutoff, arrow_class, theme_function, plotos, vector_style):
        """
        **LLM Docstring**

        Build the drawable arrow primitive for the dipole vector, anchored either at a fixed origin or at (an offset from) the molecule's center of mass, and skipped entirely if the dipole magnitude is below `mode_vector_display_cutoff`.

        :param geom: the current frame's Cartesian coordinates, used to compute the center of mass
        :type geom: np.ndarray
        :param dip: the dipole vector to draw
        :type dip: np.ndarray
        :param dipole_origin: an explicit origin (or offset, in `'shift'` mode) for the arrow
        :type dipole_origin: np.ndarray | None
        :param dipole_origin_mode: `'shift'` to add `dipole_origin` to the center of mass, otherwise use it as an absolute origin
        :type dipole_origin_mode: str
        :param mode_vector_display_cutoff: minimum vector norm required to draw the arrow at all
        :type mode_vector_display_cutoff: float
        :param arrow_class: the primitive class used to build the arrow
        :type arrow_class: type
        :param theme_function: optional callable applied to the arrow's resolved style before construction
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param vector_style: theme-level vector/arrow style options
        :type vector_style: dict
        :return: a list containing the dipole arrow primitive, or an empty list if below the display cutoff
        :rtype: list
        """
        ...

    def _get_pax_primitives(self, geom, pax: np.ndarray, *, principle_axes_origin, principle_axes_origin_mode, principle_axes_style, arrow_class, theme_function, plotos, vector_style):
        """
        **LLM Docstring**

        Build the drawable arrow primitives for the three principal axes, anchored either at a fixed origin or at (an offset from) the molecule's center of mass.

        :param geom: the current frame's Cartesian coordinates, used to compute the center of mass
        :type geom: np.ndarray
        :param pax: the `3x3` matrix of principal-axis column vectors
        :type pax: np.ndarray
        :param principle_axes_origin: an explicit origin (or offset, in `'shift'` mode) for the arrows
        :type principle_axes_origin: np.ndarray | None
        :param principle_axes_origin_mode: `'shift'` to add `principle_axes_origin` to the center of mass, otherwise use it as an absolute origin
        :type principle_axes_origin_mode: str
        :param principle_axes_style: per-axis style overrides
        :type principle_axes_style: list[dict]
        :param arrow_class: the primitive class used to build each arrow
        :type arrow_class: type
        :param theme_function: optional callable applied to each arrow's resolved style before construction
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param vector_style: theme-level vector/arrow style options
        :type vector_style: dict
        :return: the three principal-axis arrow primitives
        :rtype: list
        """
        ...

    def _get_draw_line_points(self, geom, k, v, radii):
        """
        **LLM Docstring**

        Resolve the two (and optional reference) points that define a drawn line/bond annotation: atom indices are looked up in `geom`, callables are evaluated on `geom`, and the endpoints are then pulled in by the corresponding atomic radii along the line direction.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param k: an `(xx, yy)` key, each either an atom index, a callable of `geom`, or an explicit point
        :type k: tuple
        :param v: the per-line style/options dict; consumes `'ref'` and `'offset'` from it
        :type v: dict
        :param radii: per-atom display radii, used to offset the line endpoints away from atom centers
        :type radii: list[float]
        :return: `(points, zz, normal)` -- the two (offset-adjusted) endpoint coordinates, the optional reference point (or `None`), and a normal vector for orienting a label
        :rtype: tuple[list[np.ndarray], np.ndarray | None, np.ndarray | list]
        """
        ...

    def _get_draw_line_label_props(self, xx, yy, zz, v, *, normal, default_label_style):
        """
        **LLM Docstring**

        Resolve the label text, position, and style for a drawn line annotation, if a `'label'` entry is present in `v`; computes a default offset perpendicular to the line unless one is given explicitly.

        :param xx: the first line endpoint
        :type xx: np.ndarray
        :param yy: the second line endpoint
        :type yy: np.ndarray
        :param zz: an optional reference point used to decide the label's normal/billboard behavior
        :type zz: np.ndarray | None
        :param v: the per-line style/options dict; consumes `'label'`, `'label_style'` from it
        :type v: dict
        :param normal: the fallback normal vector for the label if not billboarded
        :type normal: np.ndarray | list
        :param default_label_style: base style merged under any per-line `label_style` override
        :type default_label_style: dict
        :return: `(label_props, label_style)` where `label_props` is `(label, label_center)` or `None` if no label was requested
        :rtype: tuple[tuple | None, dict | None]
        """
        ...

    def _prep_draw_line(self, geom, key, props, *, default_label_style, line_class, theme_function, plotos, line_options, radii):
        """
        **LLM Docstring**

        Resolve everything needed to draw one line annotation: the (radius-adjusted) endpoints via `_get_draw_line_points`, the line color, and the label properties/style via `_get_draw_line_label_props`, plus the final merged drawing style.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param key: the `(xx, yy)` key identifying the line's endpoints
        :type key: tuple
        :param props: the per-line style/options dict (consumed destructively)
        :type props: dict
        :param default_label_style: base label style
        :type default_label_style: dict
        :param line_class: the primitive class that will draw the line (passed to `theme_function` only)
        :type line_class: type
        :param theme_function: optional callable applied to the resolved line style
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param line_options: theme-level line options
        :type line_options: dict
        :param radii: per-atom display radii
        :type radii: list[float]
        :return: `((xx, yy), label_props, sty, label_style)` -- the endpoints, label properties (or `None`), the resolved line style, and the label style
        :rtype: tuple
        """
        ...

    def _get_draw_coords_line(self, geom, key, props, *, radii, default_label_style, line_class, theme_function, plotos, line_options):
        """
        **LLM Docstring**

        Build the drawable primitives (a line, plus an optional label) for one `draw_coords` line entry, via `_prep_draw_line`.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param key: the `(xx, yy)` key identifying the line's endpoints
        :type key: tuple
        :param props: the per-line style/options dict
        :type props: dict
        :param radii: per-atom display radii
        :type radii: list[float]
        :param default_label_style: base label style
        :type default_label_style: dict
        :param line_class: the primitive class used to draw the line
        :type line_class: type
        :param theme_function: optional style-post-processing callable
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param line_options: theme-level line options
        :type line_options: dict
        :return: the line primitive, plus a text-label primitive if a label was requested
        :rtype: list
        """
        ...

    def _prep_draw_arc_points(self, geom, k, v, radii):
        """
        **LLM Docstring**

        Resolve the geometric parameters of a drawn arc (used for bond angles): looks up/evaluates the three key points, derives the two axis vectors (from explicit `axes` or from the point positions), pulls the endpoints in by any associated atomic radii, and determines the arc's radius and angle.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param k: the `(xx, yy, zz)` key, with `yy` as the angle vertex; each entry either an atom index, a callable, or an explicit point
        :type k: tuple
        :param v: the per-arc style/options dict; consumes `'axes'`, `'radius'`, `'angle'` from it
        :type v: dict
        :param radii: per-atom display radii, used to offset the arc endpoints away from atom centers
        :type radii: list[float]
        :return: `((xx, yy, zz), angle, axes, radius)` -- the (possibly offset-adjusted) points, the arc angle, the two unit axis vectors, and the arc radius
        :rtype: tuple
        """
        ...

    def _prep_draw_arc_label(self, yy, angle, axes, label, label_style, up_vector):
        """
        **LLM Docstring**

        Resolve the label text, position, normal, and billboard flag for a drawn arc label, placing it at the angular bisector of the arc by default.

        :param yy: the arc's vertex point
        :type yy: np.ndarray
        :param angle: the arc's angle
        :type angle: float
        :param axes: the two unit axis vectors defining the arc plane
        :type axes: list[np.ndarray]
        :param label: the label text
        :type label: str
        :param label_style: the per-arc label style dict; consumes `'offset'`, `'normal'`, `'billboard'`, `'offset_magnitude'` from it
        :type label_style: dict
        :param up_vector: fallback reference vector used if the axes are parallel (degenerate normal)
        :type up_vector: Iterable[float] | None
        :return: `(label, label_center, label_normal, label_billboard)`
        :rtype: tuple
        """
        ...

    def _prep_draw_arc(self, geom, key, props, *, up_vector, radii, default_label_style, disk_class, theme_function, plotos, disk_options):
        """
        **LLM Docstring**

        Resolve everything needed to draw one arc annotation: the arc's geometric parameters via `_prep_draw_arc_points`, the merged drawing style, and (if requested) the label properties via `_prep_draw_arc_label`.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param key: the `(xx, yy, zz)` key identifying the arc's points
        :type key: tuple
        :param props: the per-arc style/options dict (consumed destructively)
        :type props: dict
        :param up_vector: fallback reference vector for label placement
        :type up_vector: Iterable[float] | None
        :param radii: per-atom display radii
        :type radii: list[float]
        :param default_label_style: base label style
        :type default_label_style: dict
        :param disk_class: the primitive class that will draw the arc (passed to `theme_function` only)
        :type disk_class: type
        :param theme_function: optional callable applied to the resolved arc style
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param disk_options: theme-level disk/arc options
        :type disk_options: dict
        :return: `((yy, axes, angle, radius), label_props, sty, label_style)`
        :rtype: tuple
        """
        ...

    def _get_draw_coords_arc(self, geom, key, props, *, up_vector, radii, default_label_style, disk_class, theme_function, plotos, disk_options):
        """
        **LLM Docstring**

        Build the drawable primitives (an arc disk, plus an optional label) for one `draw_coords` angle entry, via `_prep_draw_arc`.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param key: the `(xx, yy, zz)` key identifying the arc's points
        :type key: tuple
        :param props: the per-arc style/options dict
        :type props: dict
        :param up_vector: fallback reference vector for label placement
        :type up_vector: Iterable[float] | None
        :param radii: per-atom display radii
        :type radii: list[float]
        :param default_label_style: base label style
        :type default_label_style: dict
        :param disk_class: the primitive class used to draw the arc
        :type disk_class: type
        :param theme_function: optional style-post-processing callable
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param disk_options: theme-level disk options
        :type disk_options: dict
        :return: the arc disk primitive, plus a text-label primitive if a label was requested
        :rtype: list
        """
        ...

    def _get_atom_text_primitives(self, geom, atom_text, *, radii, plotos):
        """
        **LLM Docstring**

        Build free-form text-label primitives (from the `atom_text` specification produced by `_prep_display_atom_text`), positioning each label near its atom (offset by the atom's radius unless an explicit position/offset is given).

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param atom_text: per-atom label specs (each `None`, a string, or a dict with at least a `'text'` key)
        :type atom_text: list
        :param radii: per-atom display radii, used for the default label offset
        :type radii: list[float]
        :param plotos: base plot-wide options merged into each label's style
        :type plotos: dict
        :return: the list of text-label primitives (skipping atoms with no text)
        :rtype: list
        """
        ...

    def _get_draw_coords_dihed(self, geom, key, props, *, default_label_style, disk_class, theme_function, plotos, disk_options):
        """
        **LLM Docstring**

        Build the drawable primitives (an arc disk plus an optional label) visualizing a dihedral angle defined by four points, computing the arc plane/axes/radius from the projected geometry rather than reusing the 3-point arc helpers.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param key: the `(xx, yy, zz, ll)` key identifying the four dihedral-defining points (atom indices, callables, or explicit points)
        :type key: tuple
        :param props: the per-arc style/options dict (consumed destructively); consumes `'angle'`, `'radius'`, `'label'`, `'label_style'`
        :type props: dict
        :param default_label_style: base label style
        :type default_label_style: dict
        :param disk_class: the primitive class used to draw the arc
        :type disk_class: type
        :param theme_function: optional callable applied to the resolved arc/label style
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param disk_options: theme-level disk options
        :type disk_options: dict
        :return: the arc disk primitive, plus a text-label primitive if a label was requested
        :rtype: list
        """
        ...

    def _get_mode_vector_primitives(self, geom, mode_vectors, *, mode_vector_display_cutoff, mode_vector_origins, mode_vector_origin_mode, arrow_class, theme_function, plotos, vector_style):
        """
        **LLM Docstring**

        Build the drawable arrow primitives for a set of per-atom normal-mode displacement vectors, anchored either at each atom's position or at (an offset from) it, and skipping any vector below `mode_vector_display_cutoff`.

        :param geom: the current frame's Cartesian coordinates
        :type geom: np.ndarray
        :param mode_vectors: the per-atom displacement vectors
        :type mode_vectors: np.ndarray
        :param mode_vector_display_cutoff: minimum vector norm required to draw an arrow for that atom
        :type mode_vector_display_cutoff: float
        :param mode_vector_origins: per-atom explicit origins (or offsets, in `'shift'` mode) for the arrows
        :type mode_vector_origins: np.ndarray | None
        :param mode_vector_origin_mode: `'shift'` to add `mode_vector_origins[j]` to the atom position, otherwise use it as an absolute origin
        :type mode_vector_origin_mode: str
        :param arrow_class: the primitive class used to build each arrow
        :type arrow_class: type
        :param theme_function: optional callable applied to each arrow's resolved style before construction
        :type theme_function: callable | None
        :param plotos: base plot-wide options
        :type plotos: dict
        :param vector_style: theme-level vector/arrow style options
        :type vector_style: dict
        :return: the list of mode-vector arrow primitives (one per atom whose displacement exceeds the cutoff)
        :rtype: list
        """
        ...

    @classmethod
    def plot_molecule(cls, mol, *geometries, figure=None, return_objects=False, bonds=None, bond_radius=None, atom_radius_scaling=None, atom_style=None, atom_radii=None, atom_text=None, display_atom_numbers=False, radius_type=None, bond_style=None, reconcile_bonds=True, capped_bonds=None, reflectiveness=None, vector_style=None, highlight_atoms=None, highlight_bonds=None, highlight_rings=None, highlight_styles=None, comparison_styles=None, animation_frame_styles=None, mode_vectors=None, mode_vector_origins=None, mode_vector_origin_mode='set', mode_vector_display_cutoff=0.01, principle_axes=None, principle_axes_origin=None, principle_axes_origin_mode='set', principle_axes_style=None, dipole=None, dipole_origin=None, dipole_origin_mode='set', render_multiple_bonds=None, render_fractional_bonds=None, fractional_bond_offset=None, bond_center_radius_offset=None, draw_coords=None, draw_coords_style=None, up_vector=None, multiple_bond_spacing=None, mode=None, backend=None, include_save_buttons=None, objects=False, graphics_class=None, cylinder_class=None, cylinder_options=None, sphere_class=None, sphere_options=None, arrow_class=None, arrow_options=None, line_class=None, line_options=None, disk_class=None, disk_options=None, animate=None, recording_options=None, animation_options=None, jsmol_load_script=None, include_jsmol_script_interface=None, dynamic_loading=None, units='Angstroms', label_style=None, theme='default', theme_function=None, plot_range_padding='auto', annotation_function=None, **plot_ops):
        """
        **LLM Docstring**

        Top-level entry point for plotting a molecule: gathers the full (very large) set of plotting options into one dict, resolves the display mode/backend and the concrete plotter class to dispatch to, normalizes the requested geometries into a `CoordinateSet`, resolves the theme, and finally delegates to the constructed plotter's `plot` method.

        :param mol: the molecule to plot
        :type mol: AbstractMolecule
        :param geometries: zero, one, or more geometries to render (zero uses `mol.coords`)
        :type geometries: np.ndarray
        :param figure: an existing figure to draw into, if continuing/overlaying a plot
        :type figure: object | None
        :param return_objects: whether to return the constructed graphics primitive objects alongside the figure
        :type return_objects: bool
        :param bonds: bonds to draw instead of `mol.bonds`; `False` to draw none
        :type bonds: tuple | bool | None
        :param mode: the requested display mode
        :type mode: str | None
        :param backend: the requested rendering backend
        :type backend: str | None
        :param theme: the theme name to select (e.g. `'default'`, `'simple'`, `'flat'`)
        :type theme: str
        :param plot_ops: every other plotting option (atom/bond styling, highlighting, mode vectors, dipole, principal axes, drawn coordinates, backend-specific figure/primitive classes and options, animation, and more) -- see the full parameter list in the source
        :type plot_ops: dict
        :return: the resulting figure (and, depending on `return_objects`/`objects`, the constructed graphics primitives)
        :rtype: object
        """
        ...

    def plot(self, **new_opts):
        """
        **LLM Docstring**

        Instance-level plotting entry point: merges any newly supplied options on top of the options stored at construction time, then dispatches to the concrete subclass's `plot_impl`.

        :param new_opts: options to override/add to `self.opts`
        :type new_opts: dict
        :return: whatever `plot_impl` returns for this plotter
        :rtype: object
        """
        ...

    def plot_impl(self, full_opts):
        """
        **LLM Docstring**

        Abstract rendering hook that concrete plotter subclasses must implement to actually produce a figure from the resolved options. Not implemented on the base class.

        :param full_opts: the fully merged/resolved plotting options
        :type full_opts: dict
        :return: never returns on the base class
        :rtype: object
        :raises NotImplementedError: always, on the base class
        """
        ...

class JupyterMoleculePlotter(MoleculePlotter):
    modes = ('jupyter',)

    def jupyter_viz(self):
        """
        **LLM Docstring**

        Build a `McUtils.Jupyter.MoleculeGraphics` widget for this molecule's atoms, Cartesian coordinates, and bonds.

        :return: the constructed Jupyter molecule-graphics widget
        :rtype: MoleculeGraphics
        """
        ...

    def plot_impl(self, cur_opts):
        """
        **LLM Docstring**

        Render this plotter's molecule as a Jupyter widget via `jupyter_viz`.

        :param cur_opts: the resolved plotting options, forwarded as keyword arguments to `jupyter_viz`
        :type cur_opts: dict
        :return: the constructed Jupyter molecule-graphics widget
        :rtype: MoleculeGraphics
        """
        ...

class JSMolMoleculePlotter(MoleculePlotter):
    modes = ('jsmol',)

    def jsmol_viz(self, xyz=None, animate=False, vibrate=False, script=None, include_script_interface=False, image_size=None, width=None, height=None, figure=None, **etc):
        """
        **LLM Docstring**

        Build a `JSMol.Applet` widget displaying this molecule, generating the XYZ block from the current coordinates if one isn't supplied, and resolving a combined width/height from `image_size` if given.

        :param xyz: an already-formatted XYZ block to display instead of generating one from `self.coords`
        :type xyz: str | None
        :param animate: whether the applet should animate between frames
        :type animate: bool
        :param vibrate: whether the applet should show a vibration animation
        :type vibrate: bool
        :param script: a Jmol script to load into the applet
        :type script: str | list[str] | None
        :param include_script_interface: whether to expose a script-input interface in the widget
        :type include_script_interface: bool
        :param image_size: a single size (applied to both width and height) or `(width, height)` pair
        :type image_size: float | tuple | None
        :param width: applet width, used if `image_size` doesn't supply one
        :type width: float | None
        :param height: applet height, used if `image_size` doesn't supply one
        :type height: float | None
        :param figure: accepted for interface consistency but not used in this method's body
        :type figure: object | None
        :param etc: additional keyword arguments forwarded to `JSMol.Applet`
        :type etc: dict
        :return: the constructed JSMol applet widget
        :rtype: object
        """
        ...

    def plot_impl(self, full_opts):
        """
        **LLM Docstring**

        Render this plotter's molecule as a JSMol applet: builds the JSMol-specific plotting options via `_prep_jsmol_plot_opts`, constructs the applet via `jsmol_viz`, and applies any shared backend figure options.

        :param full_opts: the fully merged/resolved plotting options
        :type full_opts: dict
        :return: the constructed JSMol applet figure
        :rtype: object
        """
        ...

    @classmethod
    def _jsmol_atom_sel_block(cls, atom_offset, atoms, *exprs):
        """
        **LLM Docstring**

        Build a Jmol script block that selects a set of atoms (by 0-based indices, offset by `atom_offset`, or an already-formatted selection string), applies one or more Jmol command expressions to the selection, and then deselects.

        :param atom_offset: offset added to integer atom indices to account for atoms already loaded into the figure
        :type atom_offset: int
        :param atoms: the atom indices to select, or an already-formatted Jmol selection expression string
        :type atoms: Iterable[int] | str
        :param exprs: one or more Jmol command strings to apply to the selection
        :type exprs: str
        :return: the list of Jmol script lines implementing the select/apply/deselect block
        :rtype: list[str]
        """
        ...

    def _jsmol_view_settings(self, up_vector=None, right_vector=None, view_vector=None, view_distance=None, view_matrix=None, view_center=None):
        """
        **LLM Docstring**

        Resolve a JSMol/X3D-style view specification (up/right/view vectors, distance, matrix, or center) into the settings dict used by the Jmol load script (rotation axis/angle in degrees, center, distance), via `McUtils.Plots.X3DScene.get_view_settings`.

        :param up_vector: the desired "up" direction for the view
        :type up_vector: Iterable[float] | None
        :param right_vector: the desired "right" direction for the view
        :type right_vector: Iterable[float] | None
        :param view_vector: the desired viewing direction
        :type view_vector: Iterable[float] | None
        :param view_distance: the desired viewing distance
        :type view_distance: float | None
        :param view_matrix: an explicit view/rotation matrix, if given directly
        :type view_matrix: np.ndarray | None
        :param view_center: the point the view should be centered on
        :type view_center: Iterable[float] | None
        :return: the resolved view-settings dict, with `view_angle` converted to degrees
        :rtype: dict
        """
        ...

    def _prep_jsmol_load_script(self, geom=None, background=None, atom_style=None, bond_style=None, use_default_radii=True, atom_radius_scaling=None, atom_radii=None, radius_type=None, bond_radius=None, highlight_atoms=None, highlight_bonds=None, highlight_rings=None, highlight_styles=None, display_atom_numbers=None, jsmol_load_script=None, atom_offset=0, use_default_bonds=True, draw_coords=None, draw_coords_style=None, reflectiveness=None, view_settings=None, **ignored):
        """
        **LLM Docstring**

        Build the full Jmol load-script line list for displaying this molecule with the requested styling: sets background, custom bonds (if `use_default_bonds` is `False`), atom/bond highlighting, per-atom/bond style overrides, custom atom radii/bond wireframe (if `use_default_radii` is `False`), atom-number labels, reflectiveness, extra drawn coordinates (bonds/angles as `draw` primitives, via `_prep_draw_line`/`_prep_draw_arc`), and camera/view settings.

        :param geom: the Cartesian coordinates (already unit-converted) to draw annotations relative to
        :type geom: np.ndarray | None
        :param background: the background color/spec to set
        :type background: str | None
        :param atom_style: per-atom style overrides, as produced by `_prep_display_atom_style`
        :type atom_style: dict | None
        :param bond_style: per-bond style overrides, as produced by `_prep_display_bond_style`
        :type bond_style: dict | None
        :param use_default_radii: whether to leave JSMol's default atom radii/bond wireframe alone rather than overriding them
        :type use_default_radii: bool
        :param atom_radius_scaling: uniform or per-atom radius scale factor
        :type atom_radius_scaling: float | list | None
        :param atom_radii: explicit atom radius override(s)
        :type atom_radii: float | dict | list | None
        :param radius_type: the radius field to fall back to when computing radii
        :type radius_type: str | None
        :param bond_radius: wireframe bond radius override
        :type bond_radius: float | None
        :param highlight_atoms: atom indices to highlight
        :type highlight_atoms: Iterable[int] | None
        :param highlight_bonds: bonds to highlight (their endpoint atoms are added to the highlight set)
        :type highlight_bonds: Iterable | None
        :param highlight_rings: rings (atom-index sequences) to highlight
        :type highlight_rings: Iterable | None
        :param highlight_styles: color/glow style to apply to highlighted atoms; defaults to `self.highlight_styles`
        :type highlight_styles: dict | None
        :param display_atom_numbers: `True` to label every atom with its index, or a subset of indices
        :type display_atom_numbers: bool | Iterable[int] | None
        :param jsmol_load_script: extra Jmol script text/lines to prepend
        :type jsmol_load_script: str | list[str] | None
        :param atom_offset: offset added to atom indices to account for atoms already loaded
        :type atom_offset: int
        :param use_default_bonds: whether to leave JSMol's automatic bonding alone rather than issuing explicit `connect` commands
        :type use_default_bonds: bool
        :param draw_coords: extra bond/angle annotations to draw, as produced by `_prep_display_draw_coords`
        :type draw_coords: dict | None
        :param draw_coords_style: base style for drawn-coordinate annotations
        :type draw_coords_style: dict | None
        :param reflectiveness: reflectiveness setting to apply (`True`/`False`/a 0-1 value)
        :type reflectiveness: bool | float | None
        :param view_settings: camera/view settings, as accepted by `_jsmol_view_settings`
        :type view_settings: dict | None
        :param ignored: any other options, accepted but not used
        :type ignored: dict
        :return: the list of Jmol script lines to load
        :rtype: list[str]
        """
        ...

    def _prep_jsmol_plot_opts(self, geometries=None, figure=None, atom_offset=0, extra_opts=None, script=None, jsmol_load_script=None, draw_coords=None, recording_options=None, dynamic_loading=None, include_jsmol_script_interface=None, include_save_buttons=None, **etc):
        """
        **LLM Docstring**

        Assemble the keyword arguments to pass to `jsmol_viz`: builds (or extends, if appending to an existing `figure`) the XYZ block, resolves whether to include the script interface, and (unless an explicit `script` was given) builds the full load script via `_prep_jsmol_load_script`.

        :param geometries: the geometry (or geometries) to render; only a single frame is supported
        :type geometries: np.ndarray | None
        :param figure: an existing JSMol figure to append this molecule's atoms to
        :type figure: object | None
        :param atom_offset: starting atom-index offset (incremented if appending to `figure`)
        :type atom_offset: int
        :param extra_opts: extra options bag; consumes `'xyz'`, `'use_default_radii'`, `'use_default_bonds'`, `'view_settings'`, `'background'`, `'include_script_interface'`
        :type extra_opts: dict | None
        :param script: an explicit Jmol script to use instead of building one
        :type script: str | list[str] | None
        :param jsmol_load_script: extra script lines to prepend when building the script
        :type jsmol_load_script: str | list[str] | None
        :param draw_coords: extra bond/angle annotations to draw
        :type draw_coords: dict | None
        :param recording_options: recording configuration forwarded to `jsmol_viz`
        :type recording_options: dict | None
        :param dynamic_loading: whether to use dynamic loading, forwarded to `jsmol_viz`
        :type dynamic_loading: bool | None
        :param include_jsmol_script_interface: fallback for `include_script_interface` if not given directly
        :type include_jsmol_script_interface: bool | None
        :param include_save_buttons: fallback source for `include_jsmol_script_interface` if that is also unset
        :type include_save_buttons: bool | None
        :param etc: remaining style/annotation options forwarded to `_prep_jsmol_load_script`
        :type etc: dict
        :return: the keyword-argument dict to pass to `jsmol_viz`
        :rtype: dict
        :raises ValueError: if `geometries` has more than 2 non-trivial dimensions (multiple frames aren't supported by this backend)
        """
        ...

class RDKitMoleculePlotter(MoleculePlotter):
    modes = ('rdkit', 'rdkit3d')

    def plot_impl(self, full_opts):
        """
        **LLM Docstring**

        Render this plotter's molecule via RDKit: dispatches to the molecule's `rdmol.plot` (for the `'rdkit3d'` mode) or `rdmol.draw` (for the flat `'rdkit'` mode) with the appropriately prepared options, then applies any shared backend figure options.

        :param full_opts: the fully merged/resolved plotting options
        :type full_opts: dict
        :return: the resulting RDKit figure
        :rtype: object
        """
        ...

    def _prep_rdkit_plot_opts(self, extra_opts=None, **etc):
        """
        **LLM Docstring**

        Prepare the options for the 3D RDKit plotting path (`rdmol.plot`); currently just passes through the extra plot options unchanged.

        :param extra_opts: free-form extra plotting options
        :type extra_opts: dict | None
        :param etc: all other resolved options, accepted but not used
        :type etc: dict
        :return: the (currently unmodified) extra options dict
        :rtype: dict
        """
        ...

    def _prep_rdkit_draw_opts(self, figure=None, atom_style=None, bond_style=None, atom_radii=None, atom_radius_scaling=None, radius_type=None, bond_radius=None, highlight_atoms=None, highlight_atom_radii=None, highlight_bonds=None, highlight_bond_radii=None, highlight_color=None, highlight_rings=None, highlight_styles=None, draw_coords=None, extra_opts=None, include_save_buttons=None, display_atom_numbers=None, label_style=None, **etc):
        """
        **LLM Docstring**

        Prepare the options for the flat 2D RDKit drawing path (`rdmol.draw`): resolves highlight colors from `highlight_styles` (blending glow/color if both given), computes explicit atom radii/bond radius unless RDKit's defaults are requested, and assembles the resulting keyword-argument dict.

        :param figure: an existing figure to draw into
        :type figure: object | None
        :param atom_style: per-atom style overrides (accepted but not directly used in the returned dict)
        :type atom_style: dict | None
        :param bond_style: per-bond style overrides (accepted but not directly used in the returned dict)
        :type bond_style: dict | None
        :param atom_radii: explicit atom radius override(s)
        :type atom_radii: float | dict | None
        :param atom_radius_scaling: uniform or per-atom radius scale factor
        :type atom_radius_scaling: float | list | None
        :param radius_type: the radius field to use when computing explicit radii
        :type radius_type: str | None
        :param bond_radius: explicit bond radius override
        :type bond_radius: float | None
        :param highlight_atoms: atom indices to highlight
        :type highlight_atoms: Iterable[int] | None
        :param highlight_atom_radii: per-highlighted-atom radius override
        :type highlight_atom_radii: dict | None
        :param highlight_bonds: bonds to highlight
        :type highlight_bonds: Iterable | None
        :param highlight_bond_radii: per-highlighted-bond radius override
        :type highlight_bond_radii: dict | None
        :param highlight_color: an explicit highlight color, overriding the derived one
        :type highlight_color: object | None
        :param highlight_rings: rings to highlight
        :type highlight_rings: Iterable | None
        :param highlight_styles: color/glow style used to derive `highlight_color` if not given directly; defaults to `self.highlight_styles`
        :type highlight_styles: dict | None
        :param draw_coords: extra bond/angle annotations to draw
        :type draw_coords: dict | None
        :param extra_opts: extra options bag; consumes `'use_default_radii'`
        :type extra_opts: dict | None
        :param include_save_buttons: whether to include save buttons in the figure
        :type include_save_buttons: bool | None
        :param display_atom_numbers: whether/which atoms to label with their index
        :type display_atom_numbers: bool | Iterable[int] | None
        :param label_style: base label style
        :type label_style: dict | None
        :param etc: any other options, accepted but not used
        :type etc: dict
        :return: the keyword-argument dict to pass to `rdmol.draw`
        :rtype: dict
        """
        ...

class Graphics3DMoleculePlotter(MoleculePlotter):
    """
    The main path (`x3d` / `matplotlib3D` / `plotly3D` / `svg3D`). This is the
    original `plot` body from the `extra_opts` split onward, carried over as-is;
    `self.<helper>` resolves on the plotter and `self.<mol attr>` is forwarded to
    the molecule, so nothing in the loop needed rewriting.
    """
    modes = ('x3d', 'matplotlib3D', 'plotly3D', 'svg3D', 'fast')

    def plot_impl(self, full_opts):
        """
        **LLM Docstring**

        The main 3D rendering path (used for `x3d`/`matplotlib3D`/`plotly3D`/`svg3D`): splits out `Graphics3D`-level figure options from the per-primitive plotting options, resolves every display option (atom/bond styling and radii, highlighting, dipole, principal axes, mode vectors, drawn coordinates, plot range), builds the corresponding `Sphere`/`Cylinder`/`Arrow`/`Line`/`Disk` primitives for each requested geometry frame via the various `_get_*_primitives`/`_prep_draw_*` helpers, renders them into a (possibly newly constructed) `Graphics3D` figure, and returns the figure (optionally alongside the constructed objects).

        :param full_opts: the fully merged/resolved plotting options (the complete keyword-argument set accepted by `plot_molecule`)
        :type full_opts: dict
        :return: the resulting figure, or (depending on `objects`/`return_objects`) the figure together with the constructed atom/bond/label objects
        :rtype: object
        """
        ...

class SVG2DMoleculePlotter(MoleculePlotter):
    """
    Flat 2D structure drawing on a `Graphics` (2D) SVG backend -- the plotter form
    of the old `Molecule.plot_2d`.

    Atoms are rendered as element-symbol **text** (skeletal-formula style), not as
    disks. By default carbons and hydrogens are left implicit -- pass
    ``draw_carbons=True`` / ``draw_hydrogens=True`` to label them. Where an atom is
    unlabeled its bonds run all the way to the vertex center; where it is labeled
    the bonds are trimmed back (by the atom radius) to leave room for the glyph. An
    explicit ``atom_text`` entry forces a label (and its text) even on an otherwise
    omitted carbon/hydrogen.

    Reuses the base display-prep methods (`_prep_display_atom_style`,
    `_prep_display_bond_style`, `_prep_display_atom_text`, `_get_atom_radii`), so
    coloring and highlighting behave like the 3D path. What's 2D-specific lives
    here: the embedding (`_embed_2d`), the bond-length layout (`_apply_bond_lengths`),
    and the flat `Line`/`Text` builders.

    Options beyond the shared `plot` surface arrive via `**plot_ops`:
        pose, principal_axis_order, bond_lengths, bond_layout_root,
        trim_bonds, half_colored_bonds, draw_carbons, draw_hydrogens,
        text_class, image_size.
    (`disk_class` / `disk_options` are accepted but unused now that atoms are text.)
    """
    modes = ('svg2d',)
    atom_color_updates = {'C': 'black'}
    highlight_styles = {'glow': '#ffc449'}

    def _embed_2d(self, coords, masses, pose=None, principal_axis_order=(0, 1)):
        """
        **LLM Docstring**

        Project a 3D geometry into 2D for the flat skeletal-style depiction: by default, projects onto two of the molecule's principal axes (about its center of mass); alternatively accepts a custom `pose` -- a callable, an already-2D coordinate array, a `(2,3)`/`(3,3)` projection/rotation matrix, or a full `(N,3)` alternate geometry to project instead.

        :param coords: the 3D Cartesian coordinates to embed
        :type coords: np.ndarray
        :param masses: the atomic masses, used to compute the center of mass and moments of inertia
        :type masses: np.ndarray
        :param pose: `None` for the default principal-axis projection, or a custom pose specification (callable, `(N,2)` array, `(2,3)`/`(3,3)` matrix, or `(N,3)` array)
        :type pose: callable | np.ndarray | None
        :param principal_axis_order: which two principal axes (by index) to use as the 2D `x`/`y` axes in the default projection
        :type principal_axis_order: tuple[int, int]
        :return: the 2D embedded coordinates, one row per atom
        :rtype: np.ndarray
        :raises ValueError: if `pose` is an array whose shape doesn't match any of the recognized conventions
        """
        ...

    def _apply_bond_lengths(self, xy, bond_list, bond_lengths, root=0):
        """
        **LLM Docstring**

        Adjust the 2D embedded coordinates so that bonded atom pairs sit at specified target bond lengths, propagating the adjustment outward from `root` via a breadth-first traversal of the bond graph (each newly visited atom is moved along its existing bond direction to hit the target length).

        :param xy: the initial 2D coordinates to adjust
        :type xy: np.ndarray
        :param bond_list: the bonds defining the molecular graph to traverse
        :type bond_list: list[tuple]
        :param bond_lengths: `None` to leave `xy` unchanged, a single number for a uniform target length, or a dict (keyed by index pair, element-symbol pair, or `'default'`) giving per-bond-type target lengths
        :type bond_lengths: float | dict | None
        :param root: the atom index to start the breadth-first propagation from
        :type root: int
        :return: the adjusted 2D coordinates
        :rtype: np.ndarray
        :raises ValueError: if `bond_lengths` is a value type that can't be interpreted
        """
        ...

    def _handle_overlaps_repulsions(self, xy, offset_radius=0.5, max_adjustment=1, max_iterations=10):
        """
        **LLM Docstring**

        Iteratively nudge overlapping 2D atom positions apart: repeatedly finds atom pairs closer than `offset_radius` (that haven't already been displaced beyond `max_adjustment`) and pushes each apart along their connecting direction, up to `max_iterations` passes.

        :param xy: the 2D coordinates to adjust
        :type xy: np.ndarray
        :param offset_radius: the minimum allowed distance between any two atoms
        :type offset_radius: float
        :param max_adjustment: the maximum total displacement allowed for any single atom
        :type max_adjustment: float
        :param max_iterations: the maximum number of repulsion passes to run
        :type max_iterations: int
        :return: the adjusted 2D coordinates
        :rtype: np.ndarray
        """
        ...

    def _plot_range_2d(self, xy, radii, plot_range_padding):
        """
        **LLM Docstring**

        Compute the 2D plot range spanning the given coordinates, padded either by a fixed amount, by `1.5x` the largest atomic radius (`'auto'`), or not at all.

        :param xy: the 2D coordinates to bound
        :type xy: np.ndarray
        :param radii: per-atom display radii, used when `plot_range_padding` is `'auto'`
        :type radii: np.ndarray
        :param plot_range_padding: `None` for no padding, `'auto'` for radius-based padding, or an explicit padding amount
        :type plot_range_padding: float | str | None
        :return: the `[[xmin,xmax],[ymin,ymax]]` plot range
        :rtype: list
        """
        ...

    def _atom_draw_flags(self, atom_filter, atom_style):
        """
        **LLM Docstring**

        Evaluate the atom-visibility filter for every atom, determining which atoms should be drawn with an explicit element-symbol label in the skeletal depiction.

        :param atom_filter: a callable `(element_symbol, atom_index, *, plotter, **atom_style_opts) -> bool | str` deciding whether/how to draw each atom
        :type atom_filter: callable
        :param atom_style: per-atom style dict, whose per-atom entries are passed as keyword arguments to `atom_filter`
        :type atom_style: dict
        :return: the per-atom filter results (booleans or filter-defined flag values, e.g. `'line'`)
        :rtype: list
        """
        ...

    def _get_atom_label_primitives_2d(self, xy, colors, atom_style, atom_text, labeled, *, text_class, label_style, plot_range, disk_class, radii, glows, theme_function, glow_radius, plotos, only_glows=False):
        """
        **LLM Docstring**

        Build the 2D text-label (and optional glow-disk) primitives for each atom: draws a glow disk behind labeled/glowing atoms if requested, then places the element-symbol (or overridden) text at each labeled atom's position with the resolved color/style.

        :param xy: the 2D atom coordinates
        :type xy: np.ndarray
        :param colors: per-atom colors
        :type colors: list
        :param atom_style: per-atom style overrides
        :type atom_style: dict
        :param atom_text: per-atom text overrides (from `_prep_display_atom_text`)
        :type atom_text: list
        :param labeled: per-atom flags indicating whether that atom should be labeled at all
        :type labeled: list[bool]
        :param text_class: the primitive class used to build each label
        :type text_class: type
        :param label_style: base label style
        :type label_style: dict
        :param plot_range: the current plot range, passed through to each label's style
        :type plot_range: list
        :param disk_class: the primitive class used to build glow disks
        :type disk_class: type
        :param radii: per-atom display radii
        :type radii: list[float]
        :param glows: per-atom glow colors
        :type glows: list
        :param theme_function: optional callable applied to each label's resolved style before construction
        :type theme_function: callable | None
        :param glow_radius: the radius to draw glow disks at (0 disables them)
        :type glow_radius: float
        :param plotos: base plot-wide options
        :type plotos: dict
        :param only_glows: if `True`, only build glow-disk primitives and skip the text labels
        :type only_glows: bool
        :return: the list of constructed label (and glow-disk) primitives
        :rtype: list
        """
        ...

    def _get_bond_primitives_2d(self, xy, bond_list, radii, colors, glows, labeled, *, line_class, bond_style, bond_radius, line_options, trim_bonds, render_multiple_bonds, multiple_bond_spacing, half_colored_bonds, theme_function, glow_radius, plotos, only_glows=False):
        """
        **LLM Docstring**

        Build the 2D line primitives for a set of bonds in the skeletal depiction: trims bond endpoints back from labeled atoms, draws one or more parallel offset lines depending on bond order (with optional fractional shortening), and optionally splits each line into two half-colored segments per endpoint atom (also handling a separate glow-only pass).

        :param xy: the 2D atom coordinates
        :type xy: np.ndarray
        :param bond_list: the bonds to draw, each an `(atom1, atom2, order)` tuple
        :type bond_list: list[tuple]
        :param radii: per-atom display radii, used for trimming bond endpoints
        :type radii: np.ndarray
        :param colors: per-atom colors
        :type colors: list
        :param glows: per-atom glow colors (or `None` entries to suppress glow rendering)
        :type glows: list
        :param labeled: per-atom flags indicating whether that atom is drawn with a label (used to decide whether to trim the bond endpoint)
        :type labeled: list[bool]
        :param line_class: the primitive class used to build each line segment
        :type line_class: type
        :param bond_style: per-bond/per-atom style overrides
        :type bond_style: dict
        :param bond_radius: default line width, if not already set via `line_options`
        :type bond_radius: float | None
        :param line_options: theme-level line options
        :type line_options: dict
        :param trim_bonds: whether to pull bond endpoints in from labeled atom positions
        :type trim_bonds: bool
        :param render_multiple_bonds: whether to draw multiple parallel lines for bond orders greater than 1
        :type render_multiple_bonds: bool
        :param multiple_bond_spacing: perpendicular spacing between parallel lines
        :type multiple_bond_spacing: float
        :param half_colored_bonds: whether to split each line into two segments colored by each endpoint atom when their colors/glows differ
        :type half_colored_bonds: bool
        :param theme_function: optional callable applied to each line's resolved style before construction
        :type theme_function: callable | None
        :param glow_radius: line width used for glow-highlighted segments
        :type glow_radius: float
        :param plotos: base plot-wide options
        :type plotos: dict
        :param only_glows: if `True`, only build glow-highlight line segments and skip the ordinary bond lines
        :type only_glows: bool
        :return: the list of constructed line primitives
        :rtype: list
        """
        ...

    @classmethod
    def _default_atom_filter(cls, draw_carbons, draw_hydrogens):
        """
        **LLM Docstring**

        Build the default atom-visibility filter used for the skeletal 2D depiction: carbons are hidden (drawn only as an implicit line vertex) unless `draw_carbons` is set, hydrogens are hidden entirely unless `draw_hydrogens` is set, and every other element is always shown.

        :param draw_carbons: whether carbon atoms should get an explicit element-symbol label
        :type draw_carbons: bool
        :param draw_hydrogens: whether hydrogen atoms should get an explicit element-symbol label
        :type draw_hydrogens: bool
        :return: the constructed atom filter callable
        :rtype: callable
        """
        ...

    def plot_impl(self, full_opts):
        """
        **LLM Docstring**

        The main 2D skeletal-depiction rendering path: splits out `Graphics`-level figure options, resolves all the 2D-specific options (pose, bond-length layout, overlap handling, atom filtering, etc.) alongside the shared display options, embeds each requested geometry frame into 2D (via `_embed_2d`/`_apply_bond_lengths`/`_handle_overlaps_repulsions` or a custom `layout_function`), builds the atom-label and bond-line primitives (plus optional glow passes and annotation-function extras) for each frame, and renders them into a (possibly newly constructed) 2D `Graphics` figure.

        :param full_opts: the fully merged/resolved plotting options
        :type full_opts: dict
        :return: the resulting figure, or (depending on `objects`/`return_objects`) the figure together with the constructed atom/bond/label objects, or a per-frame dict of constructed objects if `objects` is set
        :rtype: object
        """
        ...