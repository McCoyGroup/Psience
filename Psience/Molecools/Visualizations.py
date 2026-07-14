from __future__ import annotations

import itertools
import numpy as np

from McUtils.Coordinerds import CoordinateSet, CartesianCoordinates3D
from McUtils.Data import UnitsData, AtomData
import McUtils.Numputils as nput
import McUtils.Devutils as dev
import McUtils.Plots as plt

__all__ = [
    "MoleculePlotter",
    "Graphics3DMoleculePlotter",
    "JSMolMoleculePlotter",
    "RDKitMoleculePlotter",
    "JupyterMoleculePlotter",
]

class MoleculePlotter:
    """
    Base strategy: holds a molecule + resolved options, owns the shared plotting
    configuration and every data-derivation method, and dispatches to the right
    concrete plotter for the resolved mode.
    """

    # ==================================================================== #
    #  Plot configuration (moved off `Molecule`)
    # ==================================================================== #
    highlight_styles = {"glow": "green", "color": "white"}
    vector_style = {'color': 'black', 'radius': .1}
    principle_axes_style = [{'color': 'green'}, {'color': 'red'}, {'color': 'blue'}]
    draw_coords_style = {'line_color': 'black', 'line_thickness': .05}
    draw_coords_label_style = {'color': 'black', 'billboard': True}
    backend_options_resolution = {
        'mode_vectors': 'x3d',
        'dipole': 'x3d',
        'include_jsmol_script_interface': 'jsmol'
    }
    backend_aliases = {
        'x3d': ('x3d', 'x3d'),
        'jsmol': ('jsmol', 'jsmol'),
        'rdkit': ('rdkit', 'rdkit'),
        '2d': ('rdkit', 'rdkit'),
        'rdkit3d': ('rdkit3d', 'rdkit3d'),
        'matplotlib': ('fast', 'matplotlib3D'),
        'matplotlib3D': ('fast', 'matplotlib3D'),
        'matplotlib3d': ('fast', 'matplotlib3D'),
        'plotly': ('fast', 'plotly3D'),
        'plotly3D': ('fast', 'plotly3D'),
        'plotly3d': ('fast', 'plotly3D'),
        'svg3D': ('svg3D', 'svg3D'),
        'svg3d': ('svg3D', 'svg3D'),
        'svg': ('svg3D', 'svg3D'),
        'flat': ('svg3D', 'svg3D'),
        'svg2d': ('svg2d', 'svg2d'),
        'svg2D': ('svg2d', 'svg2d')
    }

    @staticmethod
    def _flat_color(i, a, styles):
        if 'color' in styles:
            styles = styles | {
                "color": "black",
                "glow": plt.ColorPalette.color_lighten(styles.get("glow", styles["color"]), -.05)
            }
        return styles

    plot_themes = {
        'default': {
            'bond_radius': .1,
            'bond_center_radius_offset': 'auto',
            'atom_radius_scaling': .25,
            'capped_bonds': False,
            'render_multiple_bonds': True,
            'render_fractional_bonds': True,
            'fractional_bond_offset': 1.1
        },
        "x3d": {
            "flat": {
                'atom_style': {"modifier": _flat_color},
                'bond_style': {"modifier": _flat_color},
                'theme_function': _flat_color
            }
        },
        'jsmol': {
            'default': {'extra_opts': {'use_default_bonds': False, 'use_default_radii': False}},
            'simple': {}
        },
        "matplotlib3D": {
            'default': {
                'multiple_bond_spacing': .15,
                'bond_center_radius_offset': {'padding': -.2, 'multi': -.5},
                'cylinder_options': {'edge_width': .05, 'edge_color': 'black', 'segments': 5},
                'sphere_options': {'edge_width': .025, 'edge_color': 'black'},
            },
            'simple': {
                'bond_radius': 0,
                'bond_style': {'color': 'none'},
                'multiple_bond_spacing': .1,
                'cylinder_options': {'edge_width': .05, 'edge_color': 'black', 'segments': 5},
                'sphere_options': {'edge_width': .025, 'edge_color': 'black'},
            }
        },
        "plotly3D": {
            'default': {
                'bond_radius': .15,
                'multiple_bond_spacing': .2,
                'bond_center_radius_offset': {'padding': .065, 'multi': -.05},
                'label_style': {'font_size': 25},
                'cylinder_options': {'edge_width': .1, 'edge_color': 'black', 'segments': 1,
                                     'hoverinfo': 'skip', 'hovertemplate': None},
                'sphere_options': {'edge_width': .05, 'edge_color': 'black',
                                   'hoverinfo': 'skip', 'hovertemplate': None},
            },
            'simple': {
                'bond_radius': 0,
                'bond_style': {'color': 'none'},
                'multiple_bond_spacing': .1,
                'label_style': {'font_size': 25},
                'cylinder_options': {'edge_width': .05, 'edge_color': 'black', 'segments': 1},
                'sphere_options': {'edge_width': .025, 'edge_color': 'black'}
            }
        },
        "svg3D": {
            'default': {
                'cylinder_options': {'stroke-width': ".01px", 'stroke': 'black'},
                'sphere_options': {'stroke-width': ".01px", 'stroke': 'black'},
                'line_options': {'stroke-width': ".05px", 'stroke': 'black'},
                'disk_options': {'stroke-width': ".05px", 'stroke': 'black'},
                'label_style': {'font_size': '.5px'}
            },
        },
        "rdkit": {
            "default": {
                "atom_radii": .15,
                "atom_radius_scaling": 1,
                "bond_radius": .2,
                'extra_opts': {'use_default_radii': False}
            }
        }
    }

    subthemes = {}

    # ==================================================================== #
    #  Registry + construction + molecule forwarding
    # ==================================================================== #
    _mode_registry: dict = {}
    modes: tuple = ()

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        for m in getattr(cls, 'modes', ()):
            MoleculePlotter._mode_registry[m] = cls
            MoleculePlotter.plot_themes[m] = cls.subthemes

    def __init__(self, mol, geometries, *, mode, backend,
                 atoms=None,
                 masses=None,
                 coords=None,
                 bonds=None,
                 **opts):
        self.mol = mol
        self.geometries = geometries
        self.mode = mode
        self.backend = backend
        self.opts = opts

        if atoms is None:
            atoms = self.mol.atoms
        self.atoms = atoms
        self._ats = [AtomData[a] for a in atoms]
        if masses is None:
            masses = [a["Mass"] for a in self._ats]
        self.masses = np.array(masses)
        if self.masses[0] < 500:
            self.atomic_masses = self.masses * UnitsData.convert("AtomicMassUnits", "ElectronMass")
        else:
            self.atomic_masses = self.masses
        self._coords = coords
        self._bonds = bonds

    @property
    def coords(self):
        if self._coords is None:
            self._coords = self.mol.coords
        return self._coords
    @property
    def bonds(self):
        if self._bonds is None:
            self._bonds = self.mol.bonds
        return self._bonds

    # def __getattr__(self, item):
    #     # Anything the plotter doesn't define (molecule data + not-yet-moved
    #     # helpers) is answered by the molecule. This is what lets the moved
    #     # method bodies -- and the `plot` body -- reference `self.<mol attr>`
    #     # unchanged.
    #     if item in ('mol', '__setstate__'):
    #         raise AttributeError(item)
    #     return getattr(self.mol, item)

    # ==================================================================== #
    #  Resolution + theming (moved off `Molecule`, keyed on the molecule)
    # ==================================================================== #
    @classmethod
    def _resolve_plot_mode(cls, mol, mode, backend, geometries, **ignored):
        if mode is None:
            if backend is not None:
                mode, backend = cls.backend_aliases.get(backend, (None, backend))
            else:
                if len(geometries) > 0:
                    backend = 'x3d'
                else:
                    for k, v in cls.backend_options_resolution.items():
                        if ignored.get(k) is not None:
                            backend = v
                            break
                mode = backend
        if mode is None:
            mode = mol.display_mode
        if backend is None:
            backend = mode
        return mode, backend

    @classmethod
    def _resolve_plot_theme(cls, mode, backend, theme, base_opts):
        if (backend, mode) in cls.plot_themes:
            theme_set = cls.plot_themes[(backend, mode)]
        elif backend in cls.plot_themes:
            theme_set = cls.plot_themes[backend]
        elif mode in cls.plot_themes:
            theme_set = cls.plot_themes[mode]
        else:
            theme_set = cls.plot_themes
        core_theme = cls.plot_themes.get('default', {}) | theme_set.get(theme, {})
        for c, v in core_theme.items():
            bop = base_opts.get(c)
            if bop is None:
                base_opts[c] = v
            elif isinstance(bop, dict):
                base_opts[c] = dev.merge_dicts(v, base_opts[c])
        return base_opts

    @classmethod
    def _default_plot_range(cls, geometries, pr, plot_range_padding, radii):
        if pr is None:
            min_any = np.min(geometries)
            max_any = np.max(geometries)
            pr = [[min_any, max_any], [min_any, max_any], [min_any, max_any]]
        if plot_range_padding is not None:
            if dev.str_is(plot_range_padding, 'auto'):
                plot_range_padding = np.max(radii)
            if nput.is_numeric(plot_range_padding):
                plot_range_padding = (plot_range_padding, plot_range_padding, plot_range_padding)
            x, y, z = plot_range_padding
            if nput.is_numeric(x):
                x = [x, x]
            if nput.is_numeric(y):
                y = [y, y]
            if nput.is_numeric(z):
                z = [z, z]
            px, py, pz = pr
            if nput.is_numeric(px):
                px = [px, px]
            px = [px[0] - x[0], px[1] + x[1]]
            if nput.is_numeric(py):
                py = [py, py]
            py = [py[0] - y[0], py[1] + y[1]]
            if nput.is_numeric(pz):
                pz = [pz, pz]
            pz = [pz[0] - z[0], pz[1] + z[1]]
            pr = (px, py, pz)
        return pr

    def _set_backend_figure_options(self,
                                    figure,
                                    mode,
                                    backend,
                                    include_save_buttons=None,
                                    dynamic_loading=None,
                                    recording_options=None,
                                    **ignored):
        if backend == 'x3d':
            if include_save_buttons is not None:
                figure.figure.include_export_button = include_save_buttons
                figure.figure.include_record_button = include_save_buttons
                figure.figure.include_view_settings_button = include_save_buttons

            if dynamic_loading is not None:
                figure.figure.dynamic_loading = dynamic_loading

            if recording_options is not None:
                figure.figure.recording_options = recording_options
        elif backend == 'plotly3D':
            if include_save_buttons is not None:
                figure.figure.include_save_buttons = True

    @staticmethod
    def _flat_color(i, a, styles):
        if 'color' in styles:
            styles = styles | {
                "color": "black",
                "glow": plt.ColorPalette.color_lighten(styles.get("glow", styles["color"]), -.05)
            }
        # if 'line_color' in styles:
        #     styles = styles | {
        #         "line_color": "line_color",
        #         "glow": styles.get("glow", styles["line_color"])
        #     }
        return styles

    def _prep_display_atom_style(self,
                                 atom_style,
                                 highlight_atoms,
                                 *,
                                 backend,
                                 reflectiveness,
                                 highlight_styles
                                 ):
        colors = [ at["IconColor"] for at in self._ats ]
        glows = [ None for at in self._ats ]
        if atom_style is None or atom_style is True:
            atom_style = {}
        elif atom_style is False:
            ...
        elif not isinstance(atom_style, dict):
            atom_style = {i: a for i, a in enumerate(atom_style)}

        if atom_style is not False:
            base_atom_style = {}
            if reflectiveness is not None and backend == 'x3d':
                base_atom_style.update({
                    'specularity': 'white',
                    'shininess': 100 * np.clip(1.1 - reflectiveness, 0, 1)
                })
            _atom_style = {i: {} for i in range(len(self._ats))}
            for k, v in atom_style.items():
                if isinstance(k, str):
                    base_atom_style[k] = v
                else:
                    _atom_style[k] = v
            for k, v in _atom_style.items():
                _atom_style[k] = base_atom_style | v
                colors[k] = v.get('color', colors[k])
                glows[k] = v.get('glow', glows[k])
            atom_style = _atom_style

        if highlight_atoms is not None:
            for k in highlight_atoms:
                if k not in atom_style: atom_style[k] = {}
                atom_style[k].update(highlight_styles)

        return atom_style, highlight_atoms, colors, glows

    def _prep_display_atom_text(self, atom_text, display_atom_numbers, label_style):
        if display_atom_numbers:
            if display_atom_numbers is True:
                display_atom_numbers = list(range(len(self._ats)))
            if not isinstance(display_atom_numbers, dict):
                display_atom_numbers = {i:{} for i in display_atom_numbers}
            atom_text = [
                {"text":i} | label_style | display_atom_numbers.get(i, {})
                    if i in display_atom_numbers else
                None
                for i in range(len(self._ats))
            ]
        if atom_text is None:
            atom_text = [None] * len(self._ats)
        return atom_text

    def _prep_display_bond_style(self,
                                 bond_style,
                                 highlight_bonds,
                                 *,
                                 backend,
                                 reflectiveness,
                                 highlight_atoms,
                                 highlight_styles,
                                 capped_bonds):
        if bond_style is None or bond_style is True:
            bond_style = {}
        elif bond_style is False:
            ...
        elif not isinstance(bond_style, dict):
            bond_style = {i:a for i,a in enumerate(bond_style)}

        if bond_style is not False:
            base_bond_style = {}
            if capped_bonds:
                base_bond_style['capped'] = True
            if reflectiveness is not None and backend == 'x3d':
                base_bond_style.update({
                    'specularity': 'white',
                    'shininess': 100 * np.clip(1.1 - reflectiveness, 0, 1)
                })
            _bond_style = {
                (i,j):{}
                for i,j in itertools.combinations(range(len(self._ats)), 2)
            }
            for k,v in bond_style.items():
                if isinstance(k, str):
                    base_bond_style[k] = v
                else:
                    _bond_style[k] = v
            for k,v in _bond_style.items():
                _bond_style[k] = dict(base_bond_style, **v)
            bond_style = _bond_style
            if highlight_bonds is None and highlight_atoms is not None:
                highlight_bonds = highlight_atoms
            if highlight_bonds is not None:
                for k in highlight_bonds:
                    if not nput.is_numeric(k): k = tuple(k)
                    if k not in bond_style: bond_style[k] = {}
                    bond_style[k].update(highlight_styles)
        return bond_style, highlight_bonds

    def _get_atomic_radius(self, atom_data, radius_type=None):
        if self.mol is not None:
            return self.mol._get_atomic_radius(atom_data, radius_type=radius_type)
        else:
            if radius_type is None:
                rad = atom_data["IconRadius"]
                if rad < .8:
                    rad = atom_data["VanDerWaalsRadius"]
            else:
                rad = atom_data[radius_type]
            return rad
    def _get_atom_radii(self, atom_radii, atom_radius_scaling, radius_type):
        if atom_radii is None:
            atom_radii = [None] * len(self._ats)
        elif nput.is_numeric(atom_radii):
            atom_radii = [atom_radii] * len(self._ats)
        elif dev.is_dict_like(atom_radii):
            atom_radii = [
                atom_radii.get(i,
                    atom_radii.get(a["ElementSymbol"])
                ) for i,a in enumerate(self._ats)
            ]
        atom_radii = [
            self._get_atomic_radius(at, radius_type)
                if c is None else
            c
            for c, at in
            zip(atom_radii, self._ats)
        ]
        if nput.is_numeric(atom_radius_scaling):
            atom_radius_scaling = [atom_radius_scaling] * len(atom_radii)
        radii = [ s * r for s, r in zip(atom_radius_scaling, atom_radii) ]
        return radii

    def _prep_display_dipole(self, geometries, dipole, dipole_origin, units):
        if dipole is not None:
            dipole = np.asanyarray(dipole)
            if dipole.ndim == 1:
                dipole = np.broadcast_to(dipole[np.newaxis], (len(geometries),) + dipole.shape)

            if units is not None:
                dipole = dipole * UnitsData.convert("BohrRadius", units)

            if dipole_origin is not None:
                dipole_origin = np.asanyarray(dipole_origin)
                if dipole_origin.ndim == 1:
                    dipole_origin = np.broadcast_to(
                        dipole_origin[np.newaxis],
                        (len(geometries),) + dipole_origin.shape
                    )

                if units is not None:
                    dipole_origin = dipole_origin * UnitsData.convert("BohrRadius", units)
        return dipole, dipole_origin

    def _prep_principle_axes(self, geometries, units, principle_axes, principle_axes_origin, principle_axes_style):
        if principle_axes is True:
            _, principle_axes = nput.moments_of_inertia(geometries, self.atomic_masses)
        elif principle_axes is False:
            principle_axes = None
        if principle_axes is not None:
            if principle_axes_style is None:
                principle_axes_style = {}
            if isinstance(principle_axes_style, dict):
                principle_axes_style = [principle_axes_style] * 3
            principle_axes_style = [
                dict(self.principle_axes_style[i], **principle_axes_style[i])
                for i in range(3)
            ]
            if principle_axes is not None:
                if principle_axes.ndim == 2:
                    principle_axes = np.broadcast_to(
                        principle_axes[np.newaxis],
                        (len(geometries),) + principle_axes.shape
                    )
                # if units is not None:
                #     principle_axes = principle_axes * UnitsData.convert("BohrRadius", units)

            if principle_axes_origin is not None:
                if principle_axes_origin.ndim == 1:
                    principle_axes_origin = np.broadcast_to(
                        principle_axes_origin[np.newaxis],
                        (len(geometries),) + principle_axes_origin.shape
                    )

                if units is not None:
                    principle_axes_origin = principle_axes_origin * UnitsData.convert("BohrRadius", units)
        return principle_axes, principle_axes_origin, principle_axes_style

    def _prep_display_mode_vectors(self, geometries, units, mode_vectors, mode_vector_origins):
        if mode_vectors is not None:
            mode_vectors = np.asanyarray(mode_vectors)
            if mode_vectors.ndim == 1:
                mode_vectors = np.reshape(mode_vectors, (-1, 3))
            if mode_vectors.ndim == 2:
                mode_vectors = np.broadcast_to(mode_vectors[np.newaxis], (len(geometries),) + mode_vectors.shape)

            if units is not None:
                mode_vectors = mode_vectors * UnitsData.convert("BohrRadius", units)

            if mode_vector_origins is not None:
                mode_vector_origins = np.asanyarray(mode_vector_origins)
                if mode_vector_origins.ndim == 1:
                    mode_vector_origins = np.reshape(mode_vector_origins, (-1, 3))
                if mode_vector_origins.ndim == 2:
                    mode_vector_origins = np.broadcast_to(mode_vector_origins[np.newaxis],
                                                          (len(geometries),) + mode_vector_origins.shape)

                if units is not None:
                    mode_vector_origins = mode_vector_origins * UnitsData.convert("BohrRadius", units)

        return mode_vectors, mode_vector_origins

    def _prep_display_draw_coords(self, draw_coords, draw_coords_style):
        if draw_coords is not None:
            if not isinstance(draw_coords, dict):
                draw_coords = {
                    tuple(a):{}
                    for a in draw_coords
                }
        if draw_coords_style is None:
            draw_coords_style = self.draw_coords_style
        return draw_coords, draw_coords_style

    def _get_bond_primitives(self,
                             geom,
                             b,
                             *,
                             bond_list,
                             bond_radius,
                             radii,
                             bond_center_radius_offset,
                             multiple_bond_spacing,
                             render_multiple_bonds,
                             render_fractional_bonds,
                             fractional_bond_offset,
                             max_bond_orders,
                             up_vector,
                             cylinder_class,
                             colors,
                             glows,
                             bond_style,
                             theme_function,
                             plotos,
                             cylinder_options
                             ):
        atom1 = b[0]
        atom2 = b[1]

        c1 = colors[atom1]
        c2 = colors[atom2]
        g1 = glows[atom1]
        g2 = glows[atom2]
        base_bstyle = dict(
            bond_style.get((atom2, atom1), {}),
            **bond_style.get((atom1, atom2), {})
        )
        b_sty_1 = dict(
            bond_style.get(atom1, {}),
            **base_bstyle
        )
        if b_sty_1.get('color') is None:
            b_sty_1['color'] = c1
        if b_sty_1.get('glow') is None and g1 is not None:
            b_sty_1['glow'] = g1
        b_sty_1 = b_sty_1.copy()
        modifier = b_sty_1.pop('modifier', None)
        if modifier is not None:
            b_sty_1 = modifier((atom1, atom2), (self.atoms[atom1], self.atoms[atom2]), b_sty_1)

        b_sty_2 = dict(
            bond_style.get(atom2, {}),
            **base_bstyle
        )
        if b_sty_2.get('color') is None:
            b_sty_2['color'] = c2
        if b_sty_2.get('glow') is None and g2 is not None:
            b_sty_2['glow'] = g2
        b_sty_2 = b_sty_2.copy()
        modifier = b_sty_2.pop('modifier', None)
        if modifier is not None:
            b_sty_2 = modifier((atom2, atom1), (self.atoms[atom2], self.atoms[atom1]), b_sty_2)

        p1 = geom[atom1]
        p2 = geom[atom2]
        disp_vector = p2 - p1
        nv, vn = nput.vec_normalize(disp_vector, return_norms=True)
        midpoint = ((p1 + nv * radii[atom1]) + (p2 - nv * radii[atom2])) / 2
        if bond_center_radius_offset is not None:
            if dev.str_is(bond_center_radius_offset, 'auto'):
                bond_center_radius_offset = {'padding': 0}
            if isinstance(bond_center_radius_offset, dict):
                pad = bond_center_radius_offset['padding']
                rz1 = (radii[atom1] ** 2 - (bond_radius + pad) ** 2)
                if rz1 > 0:
                    rz1 = np.sqrt(rz1)
                else:
                    rz1 = 0
                p1 = p1 + rz1 * nv
                rz2 = (radii[atom2] ** 2 - (bond_radius + pad) ** 2)
                if rz2 > 0:
                    rz2 = np.sqrt(rz2)
                else:
                    rz2 = 0
                p2 = p2 - rz2 * nv
            else:
                p1 = p1 + bond_center_radius_offset * radii[atom1] * nv
                p2 = p2 - bond_center_radius_offset * radii[atom2] * nv

        disp_vector = p2 - p1
        # if 'transparency' not in b_sty_1:
        #     b_sty_1['transparency'] = 0
        # if 'transparency' not in b_sty_2:
        #     b_sty_2['transparency'] = 0
        if b[2] == 0:
            bond_radius = 0
            bond_point_list = [
                [p1, p1, p1]
            ]
            # b_sty_2['color'] = 'white'
            # b_sty_1['color'] = 'white'
        elif not render_multiple_bonds or len(b) == 2 or b[2] <= 1 or b[2] > 3:
            bond_point_list = [
                [p1, p2, midpoint]
            ]
        else:
            if up_vector is None:
                up_vector = [0, 0, 1]

            for bb in bond_list:
                atom3 = bb[0]
                atom4 = bb[1]
                if atom3 in {atom1, atom2} and atom4 not in {atom1, atom2}:
                    u_vec = np.cross(disp_vector, geom[atom4] - geom[atom3])
                    break
                elif atom4 not in {atom1, atom2} and atom3 in {atom1, atom2}:
                    u_vec = np.cross(disp_vector, geom[atom3] - geom[atom4])
                    break
            else:
                u_vec = up_vector
            if multiple_bond_spacing is None:
                multiple_bond_spacing = bond_radius * 1.1

            axis = nput.vec_normalize(
                np.cross(disp_vector, u_vec)
            )

            if isinstance(bond_center_radius_offset, dict):
                pad = bond_center_radius_offset['padding']
                mpad = bond_center_radius_offset.get('multi', pad)
                rz1 = (radii[atom1] ** 2 - (bond_radius + pad) ** 2)
                if rz1 > 0:
                    rz1 = np.sqrt(rz1)
                    rz12 = (radii[atom1] ** 2 - (bond_radius + mpad + multiple_bond_spacing) ** 2)
                    if rz12 > 0:
                        rzd = rz1 - np.sqrt(rz12)
                    else:
                        rzd = rz1
                    p1 = p1 - rzd * nv

                rz2 = (radii[atom2] ** 2 - (bond_radius + pad) ** 2)
                if rz2 > 0:
                    rz2 = np.sqrt(rz2)
                    rz22 = (radii[atom2] ** 2 - (bond_radius + mpad + multiple_bond_spacing) ** 2)
                    if rz22 > 0:
                        rzd = rz2 - np.sqrt(rz22)
                    else:
                        rzd = rz2
                    p2 = p2 + rzd * nv

            if render_fractional_bonds:
                dr = (b[2] - int(b[2]))
                if dr > .05:
                    dr = dr * fractional_bond_offset
                    dr = (1 - dr) / 2
                    p10 = p1 + disp_vector * dr
                    p20 = p2 - disp_vector * dr
                else:
                    p10 = p1
                    p20 = p2
            else:
                p10 = p1
                p20 = p2

            if b[2] <= 2:

                p11 = p10 - axis * multiple_bond_spacing
                p21 = p20 - axis * multiple_bond_spacing
                mp1 = midpoint - axis * multiple_bond_spacing

                p12 = p1 + axis * multiple_bond_spacing
                p22 = p2 + axis * multiple_bond_spacing
                mp2 = midpoint + axis * multiple_bond_spacing
                bond_point_list = [
                    [p11, p21, mp1],
                    [p12, p22, mp2]
                ]
            else:
                dx = (multiple_bond_spacing)

                u_vec = nput.vec_normalize(u_vec)

                dv = axis * dx
                p11 = p10 - dv
                p21 = p20 - dv
                mp1 = midpoint - dv

                dv = (np.cos(2 * np.pi / 3) * axis + np.sin(2 * np.pi / 3) * u_vec) * dx
                p12 = p1 - dv
                p22 = p2 - dv
                mp2 = midpoint - dv

                dv = (np.cos(2 * np.pi / 3) * axis - np.sin(2 * np.pi / 3) * u_vec) * dx
                p13 = p1 - dv
                p23 = p2 - dv
                mp3 = midpoint - dv
                bond_point_list = [
                    [p11, p21, mp1],
                    [p12, p22, mp2],
                    [p13, p23, mp3],
                ]

        if max_bond_orders is not None:
            mbo = max_bond_orders.get((b[0], b[1]))
            target_bond_number = (
                1 if mbo > 3.01 else
                3 if mbo > 2.01 else
                2 if mbo > 1.01 else
                1
            )
            if len(bond_point_list) < target_bond_number:
                bond_point_list = bond_point_list + [bond_point_list[0]] * (target_bond_number - len(bond_point_list))

        bond_objs = []
        for pp1, pp2, mp in bond_point_list:
            sty1 = (plotos | cylinder_options | b_sty_1)
            if theme_function is not None:
                sty1 = theme_function((atom1, atom2), cylinder_class, sty1)
            cc1 = cylinder_class(
                pp1,
                mp,
                bond_radius,
                **sty1
            )
            sty2 = (plotos | cylinder_options | b_sty_2)
            if theme_function is not None:
                sty2 = theme_function((atom2, atom1), cylinder_class, sty2)
            cc2 = cylinder_class(
                mp,
                pp2,
                bond_radius,
                **sty2
            )
            bond_objs.extend([cc1, cc2])

        return bond_objs

    def _get_bondlist_primitives(self,
                                 geom,
                                 bond_list,
                                 *,
                                 bond_radius,
                                 radii,
                                 bond_center_radius_offset,
                                 multiple_bond_spacing,
                                 render_multiple_bonds,
                                 render_fractional_bonds,
                                 fractional_bond_offset,
                                 max_bond_orders,
                                 up_vector,
                                 cylinder_class,
                                 colors,
                                 glows,
                                 bond_style,
                                 theme_function,
                                 plotos,
                                 cylinder_options
                                 ):
        all_bonds = []
        for j, b in enumerate(bond_list):
            all_bonds.append(
                self._get_bond_primitives(
                    geom,
                    b,
                    bond_list=bond_list,
                    bond_radius=bond_radius,
                    radii=radii,
                    bond_center_radius_offset=bond_center_radius_offset,
                    multiple_bond_spacing=multiple_bond_spacing,
                    render_multiple_bonds=render_multiple_bonds,
                    render_fractional_bonds=render_fractional_bonds,
                    fractional_bond_offset=fractional_bond_offset,
                    max_bond_orders=max_bond_orders,
                    up_vector=up_vector,
                    cylinder_class=cylinder_class,
                    colors=colors,
                    glows=glows,
                    bond_style=bond_style,
                    theme_function=theme_function,
                    plotos=plotos,
                    cylinder_options=cylinder_options
                )
            )
        return all_bonds

    def _get_atom_primitives(self,
                             geom,
                             atoms,
                             *,
                             colors,
                             radii,
                             sphere_class,
                             atom_style,
                             theme_function,
                             plotos,
                             sphere_options
                             ):
        atom_objs = []
        for j, stuff in enumerate(zip(colors, radii, geom)):
            color, radius, coord = stuff
            a_sty = atom_style.get(j, {})
            if a_sty.get('color') is None:
                a_sty['color'] = color
            a_sty = a_sty.copy()
            modifier = a_sty.pop('modifier', None)
            if modifier is not None:
                a_sty = modifier(j, self.atoms[j], a_sty)

            asty = (plotos | sphere_options | a_sty)
            if theme_function is not None:
                asty = theme_function(j, sphere_class, asty)
            sphere = sphere_class(coord, radius, **asty)
            atom_objs.append(sphere)
        return atom_objs

    def _get_dipole_primitives(self,
                               geom,
                               dip,
                               *,
                               dipole_origin,
                               dipole_origin_mode,
                               mode_vector_display_cutoff,
                               arrow_class,
                               theme_function,
                               plotos,
                               vector_style
                               ):
        dipole_primitives = []
        if np.linalg.norm(dip) > mode_vector_display_cutoff:
            if dipole_origin is None or dipole_origin_mode == 'shift':
                com = np.tensordot(self.masses, geom, axes=[0, 0]) / np.sum(self.masses)
                if dipole_origin is not None:
                    com = com + dipole_origin
            else:
                com = dipole_origin
            sty = (plotos | vector_style)
            if theme_function is not None:
                sty = theme_function(None, arrow_class, sty)
            dipole_arrow = arrow_class(
                com,
                com + dip,
                **sty
            )
            dipole_primitives.append(dipole_arrow)
        return dipole_primitives

    def _get_pax_primitives(self,
                            geom,
                            pax:np.ndarray,
                            *,
                            principle_axes_origin,
                            principle_axes_origin_mode,
                            principle_axes_style,
                            arrow_class,
                            theme_function,
                            plotos,
                            vector_style
                            ):
        pax_objs = []
        if principle_axes_origin is None or principle_axes_origin_mode == 'shift':
            com = np.tensordot(self.masses, geom, axes=[0, 0]) / np.sum(self.masses)
            if principle_axes_origin is not None:
                com = com + principle_axes_origin
        else:
            com = principle_axes_origin
        for ax, sty in zip(pax.T, principle_axes_style):
            sty = (plotos | vector_style | sty)
            if theme_function is not None:
                sty = theme_function(None, arrow_class, sty)
            pax_arrow = arrow_class(
                com,
                com + ax,
                **sty
            )
            pax_objs.append(pax_arrow)
        return pax_objs

    def _get_draw_line_points(self, geom, k, v, radii):
        # draw bond
        # draw angle
        xx, yy = k
        rx = 0
        ry = 0
        if nput.is_int(xx):
            rx = radii[xx]
            xx = geom[xx]
        elif callable(xx):
            xx = xx(geom)
        if nput.is_int(yy):
            ry = radii[yy]
            yy = geom[yy]
        elif callable(yy):
            yy = yy(geom)

        zz = v.pop('ref', None)
        if zz is not None:
            if nput.is_int(zz):
                zz = geom[zz]
            elif callable(zz):
                zz = zz(geom)
            zz = np.asanyarray(zz)

        xx = np.asanyarray(xx)
        yy = np.asanyarray(yy)
        ax = nput.vec_normalize(xx - yy)
        xx = xx - rx * ax
        yy = yy + ry * ax
        if zz is not None:
            normal = nput.vec_crosses(xx - yy, zz - yy)
        else:
            normal = [0, 1, 0]
        offset = np.array(v.pop('offset', [0, 0, 0]))
        points = [xx + offset, yy + offset]
        return points, zz, normal

    def _get_draw_line_label_props(self, xx, yy, zz, v,
                                   *,
                                   normal,
                                   default_label_style
                                   ):
        label = v.pop('label', None)
        if label is not None:
            label_style = dict(default_label_style, **v.pop('label_style', {}))
            label_style['billboard'] = label_style.get('billboard',
                                                       zz is None and 'normal' not in label_style)
            if not label_style['billboard']:
                label_style['normal'] = label_style.get('normal', normal)
            offset_magnitude = label_style.pop('offset_magnitude', .2)
            offset_axis = nput.vec_crosses(xx - yy, label_style.get('normal', normal), normalize=True)
            label_offset = np.asanyarray(label_style.pop('offset', offset_magnitude * offset_axis))
            label_center = (xx + yy) / 2 + label_offset
            label_props = label, label_center
        else:
            label_props = None
            label_style = None
        return label_props, label_style

    def _prep_draw_line(self,
                        geom,
                        key,
                        props,
                        *,
                        default_label_style,
                        line_class,
                        theme_function,
                        plotos,
                        line_options,
                        radii
                        ):

        (xx, yy), zz, normal = self._get_draw_line_points(geom, key, props, radii)
        color = props.pop('line_color', props.pop('color', 'black'))

        label_props, label_style = self._get_draw_line_label_props(xx, yy, zz, props,
                                                                   normal=normal,
                                                                   default_label_style=default_label_style)

        sty = (plotos | line_options | dict(color=color) | props)
        if theme_function is not None:
            sty = theme_function(None, line_class, sty)

        return (xx, yy), label_props, sty, label_style
    def _get_draw_coords_line(self,
                              geom,
                              key,
                              props,
                              *,
                              radii,
                              default_label_style,
                              line_class,
                              theme_function,
                              plotos,
                              line_options
                              ):
        prims = []

        points, label_props, sty, label_style = self._prep_draw_line(
            geom,
            key,
            props,
            radii=radii,
            default_label_style=default_label_style,
            line_class=line_class,
            theme_function=theme_function,
            plotos=plotos,
            line_options=line_options
        )

        arc = line_class(points, **sty)
        prims.append(arc)

        if label_props is not None:
            label, label_center = label_props
            lab = plt.Text(
                label, label_center,
                **(plotos | label_style)
            )
            prims.append(lab)

        return prims

    def _prep_draw_arc_points(self,
                              geom,
                              k,
                              v,
                              radii):
        # draw angle
        xx, yy, zz = k
        rx, ry, rz = 0, 0, 0
        if nput.is_int(xx):
            rx = radii[xx]
            xx = geom[xx]
        elif callable(xx):
            xx = xx(geom)
        if nput.is_int(yy):
            ry = radii[yy]
            yy = geom[yy]
        elif callable(yy):
            yy = yy(geom)
        if nput.is_int(zz):
            rz = radii[zz]
            zz = geom[zz]
        elif callable(zz):
            zz = zz(geom)

        yy = np.asanyarray(yy)
        axes = v.pop('axes', None)
        if axes is None:
            xx = np.asanyarray(xx)
            zz = np.asanyarray(zz)
            axes = [xx - yy, zz - yy]
        else:
            axes = nput.vec_normalize(axes)
            xx = yy + axes[0] * np.linalg.norm(yy - xx)
            zz = yy + axes[1] * np.linalg.norm(yy - zz)

        radius = v.pop('radius', None)
        if radius is None:
            axes, norms = nput.vec_normalize(axes, return_norms=True)
            radius = min(norms)
        else:
            axes = nput.vec_normalize(axes)

        if rx > 0 or rz > 0:
            normal = np.cross(*axes)
            if rx > 0:
                perp_x = nput.vec_crosses(axes[0], normal, normalize=True)
                xx = xx - perp_x * rx
            if rz > 0:
                perp_z = nput.vec_crosses(normal, axes[1], normalize=True)
                zz = zz - perp_z * rz
            axes = nput.vec_normalize([xx - yy, zz - yy])

        angle = v.pop('angle', None)
        if angle is None:
            angle = nput.vec_angles(*axes, return_crosses=False)

        return (xx,yy,zz), angle, axes, radius

    def _prep_draw_arc_label(self,
                             yy, angle, axes,
                             label, label_style, up_vector):
        label_offset = np.asanyarray(label_style.pop('offset', [0, 0, 0]))
        label_normal = label_style.pop('normal', None)
        if label_normal is None:
            label_normal = nput.vec_crosses(*axes, normalize=True)
            if np.linalg.norm(label_normal) < 1e-8:
                if up_vector is None:
                    up_vector = [0, 0, 1]
                label_normal = nput.vec_crosses(axes[0], up_vector, normalize=True)
        label_billboard = np.asanyarray(label_style.pop('billboard', False))
        offset_magnitude = label_style.pop('offset_magnitude', .8)
        label_center = np.dot(
            axes[0] * offset_magnitude,
            nput.rotation_matrix(label_normal, -angle / 2)
        ) + yy
        return label, label_offset + label_center, label_normal, label_billboard

    def _prep_draw_arc(self,
                       geom, key, props,
                       *,
                       up_vector,
                       radii,
                       default_label_style,
                       disk_class,
                       theme_function,
                       plotos,
                       disk_options
                       ):

        v = props
        (xx, yy, zz), angle, axes, radius = self._prep_draw_arc_points(geom, key, props, radii)

        label = v.pop('label', None)
        label_style = dict(default_label_style, **v.pop('label_style', {}))
        sty = (plotos | disk_options | v)
        if theme_function is not None:
            sty = theme_function(None, disk_class, sty)

        if label is not None:
            label_props = self._prep_draw_arc_label(
                yy, angle, axes,
                label, label_style, up_vector
            )
        else:
            label_props = None

        return (yy, axes, angle, radius), label_props, sty, label_style

    def _get_draw_coords_arc(self,
                             geom,
                             key,
                             props,
                             *,
                             up_vector,
                             radii,
                             default_label_style,
                             disk_class,
                             theme_function,
                             plotos,
                             disk_options
                             ):
        prims = []

        arc_props, label_props, sty, label_style = self._prep_draw_arc(
            geom, key, props,
            up_vector=up_vector,
            radii=radii,
            default_label_style=default_label_style,
            disk_class=disk_class,
            theme_function=theme_function,
            plotos=plotos,
            disk_options=disk_options
        )

        yy, axes, angle, radius = arc_props
        arc = disk_class(
            yy,
            uv_axes=axes,
            angle=angle,
            radius=radius,
            **sty
        )
        prims.append(arc)

        if label_props is not None:
            label, center, normal, billboard = label_props
            lab = plt.Text(
                label,
                center,
                billboard=billboard,
                normal=normal,
                **(plotos | label_style)
            )
            prims.append(lab)


        return prims

    def _get_atom_text_primitives(self,
                                  geom,
                                  atom_text,
                                  *,
                                  radii,
                                  plotos):
        prims = []
        for a, r, t in zip(geom, radii, atom_text):
            if t is None: continue

            if not isinstance(t, dict):
                t = {"text": t}
            t = t.copy()
            text = t.pop('text')
            if nput.is_numeric(text):
                text = str(text)
            pos = t.pop('pos', a + t.pop('offset', np.array([r / 2, r / 2, r])))
            fs = t.pop('font_style', {})
            fs['size'] = fs.get('size', .5)
            t['font_style'] = fs
            t['billboard'] = t.get('billboard', True)
            t['color'] = t.get('color', 'black')
            lab = plt.Text(text, pos, **(plotos | t))
            prims.append(lab)
        return prims

    def _get_draw_coords_dihed(self,
                               geom,
                               key,
                               props,
                               *,
                               default_label_style,
                               disk_class,
                               theme_function,
                               plotos,
                               disk_options
                               ):
        prims = []
        k = key
        v = props

        # draw angle
        xx, yy, zz, ll = k
        if nput.is_int(xx):
            xx = geom[xx]
        elif callable(xx):
            xx = xx(geom)
        if nput.is_int(yy):
            yy = geom[yy]
        elif callable(yy):
            yy = yy(geom)
        if nput.is_int(zz):
            zz = geom[zz]
        elif callable(zz):
            zz = zz(geom)
        if nput.is_int(ll):
            ll = geom[ll]
        elif callable(ll):
            ll = ll(geom)

        xx = np.asanyarray(xx)
        yy = np.asanyarray(yy)
        zz = np.asanyarray(zz)
        ll = np.asanyarray(ll)

        normal = yy - zz
        ax1 = nput.project_out((xx - yy), normal[:, np.newaxis])
        ax2 = nput.project_out((ll - zz), normal[:, np.newaxis])
        angle = v.pop('angle', None)
        if angle is None:
            angle = nput.pts_dihedrals(xx, yy, zz, ll)

        c2 = (yy + zz) / 2
        if np.linalg.norm(ax1) < np.linalg.norm(ax2):
            c1 = c2 + ax1
        else:
            c1 = c2 + ax2
            normal = -normal
        axis = c1 - c2
        c3 = np.dot(
            axis,
            nput.rotation_matrix(normal, angle)
        ) + c2
        axes = [
            axis,
            c3 - c2,
        ]
        if np.dot(np.cross(*axes), normal) < 0:
            normal = -normal

        radius = v.pop('radius', None)
        if radius is None:
            radius = np.linalg.norm(axis)
            # angle=nput.vec_angles(*axes, return_crosses=False)

        label = v.pop('label', None)
        label_style = dict(default_label_style, **v.pop('label_style', {}))
        sty = (plotos | disk_options | v)
        if theme_function is not None:
            sty = theme_function(None, disk_class, sty)
        arc = disk_class(
            c2,
            uv_axes=axes,
            normal=normal,
            angle=angle,
            radius=radius,
            **sty
        )
        prims.append(arc)

        if label is not None:
            label_offset = np.asanyarray(label_style.pop('offset', [0, 0, 0]))
            label_normal = label_style.pop('normal', None)
            if label_normal is None:
                label_normal = normal
                # if np.linalg.norm(label_normal) < 1e-8:
                #     if up_vector is None:
                #         up_vector = [0, 0, 1]
                #     label_normal = nput.vec_crosses(axes[0], up_vector, normalize=True)
            label_billboard = np.asanyarray(label_style.pop('billboard', False))
            offset_magnitude = label_style.pop('offset_magnitude', .8)
            label_center = np.dot(
                axes[0] * offset_magnitude,
                nput.rotation_matrix(label_normal, -angle / 2)
            ) + c2
            lab = plt.Text(
                label,
                label_center + label_offset,
                billboard=label_billboard,
                normal=label_normal,
                **(plotos | label_style)
            )
            prims.append(lab)

        return prims

    def _get_mode_vector_primitives(self,
                                    geom,
                                    mode_vectors,
                                    *,
                                    mode_vector_display_cutoff,
                                    mode_vector_origins,
                                    mode_vector_origin_mode,
                                    arrow_class,
                                    theme_function,
                                    plotos,
                                    vector_style
                                    ):
        mode_arrows = []
        for j, v in enumerate(mode_vectors):
            if np.linalg.norm(v) > mode_vector_display_cutoff:
                if mode_vector_origins is None or mode_vector_origin_mode == 'shift':
                    com = geom[j]
                    if mode_vector_origins is not None:
                        com = com + mode_vector_origins[j]
                else:
                    com = mode_vector_origins[j]
                sty = (plotos | vector_style)
                if theme_function is not None:
                    sty = theme_function(None, arrow_class, sty)
                mode_arrow = arrow_class(
                    com,
                    com + v,
                    **sty
                )
                mode_arrows.append(mode_arrow)
        return mode_arrows

    # ==================================================================== #
    #  Dispatch -- the plot entry point `Molecule.plot` forwards to
    # ==================================================================== #
    @classmethod
    def plot_molecule(cls,
                      mol,
                      *geometries,
                      figure=None,
                      return_objects=False,
                      bonds=None,
                      bond_radius=None,
                      atom_radius_scaling=None,
                      atom_style=None,
                      atom_radii=None,
                      atom_text=None,
                      display_atom_numbers=False,
                      radius_type=None,
                      bond_style=None,
                      reconcile_bonds=True,
                      capped_bonds=None,
                      reflectiveness=None,
                      vector_style=None,
                      highlight_atoms=None,
                      highlight_bonds=None,
                      highlight_rings=None,
                      highlight_styles=None,
                      comparison_styles=None,
                      animation_frame_styles=None,
                      mode_vectors=None,
                      mode_vector_origins=None,
                      mode_vector_origin_mode='set',
                      mode_vector_display_cutoff=1e-2,
                      principle_axes=None,
                      principle_axes_origin=None,
                      principle_axes_origin_mode='set',
                      principle_axes_style=None,
                      dipole=None,
                      dipole_origin=None,
                      dipole_origin_mode='set',
                      render_multiple_bonds=None,
                      render_fractional_bonds=None,
                      fractional_bond_offset=None,
                      bond_center_radius_offset=None,
                      draw_coords=None,
                      draw_coords_style=None,
                      up_vector=None,
                      multiple_bond_spacing=None,
                      mode=None,
                      backend=None,
                      include_save_buttons=None,
                      objects=False,
                      graphics_class=None,
                      cylinder_class=None,
                      cylinder_options=None,
                      sphere_class=None,
                      sphere_options=None,
                      arrow_class=None,
                      arrow_options=None,
                      line_class=None,
                      line_options=None,
                      disk_class=None,
                      disk_options=None,
                      animate=None,
                      recording_options=None,
                      animation_options=None,
                      jsmol_load_script=None,
                      include_jsmol_script_interface=None,
                      dynamic_loading=None,
                      units="Angstroms",
                      label_style=None,
                      theme='default',
                      theme_function=None,
                      plot_range_padding='auto',
                      annotation_function=None,
                      **plot_ops
                      ):
        full_opts = dict(
            return_objects=return_objects,
            bonds=bonds,
            bond_radius=bond_radius,
            atom_radius_scaling=atom_radius_scaling,
            atom_style=atom_style,
            atom_radii=atom_radii,
            atom_text=atom_text,
            display_atom_numbers=display_atom_numbers,
            radius_type=radius_type,
            bond_style=bond_style,
            capped_bonds=capped_bonds,
            reflectiveness=reflectiveness,
            vector_style=vector_style,
            highlight_atoms=highlight_atoms,
            highlight_bonds=highlight_bonds,
            highlight_rings=highlight_rings,
            highlight_styles=highlight_styles,
            comparison_styles=comparison_styles,
            animation_frame_styles=animation_frame_styles,
            mode_vectors=mode_vectors,
            mode_vector_origins=mode_vector_origins,
            mode_vector_origin_mode=mode_vector_origin_mode,
            mode_vector_display_cutoff=mode_vector_display_cutoff,
            principle_axes=principle_axes,
            principle_axes_origin=principle_axes_origin,
            principle_axes_origin_mode=principle_axes_origin_mode,
            principle_axes_style=principle_axes_style,
            dipole=dipole,
            dipole_origin=dipole_origin,
            dipole_origin_mode=dipole_origin_mode,
            render_multiple_bonds=render_multiple_bonds,
            render_fractional_bonds=render_fractional_bonds,
            fractional_bond_offset=fractional_bond_offset,
            bond_center_radius_offset=bond_center_radius_offset,
            draw_coords=draw_coords,
            draw_coords_style=draw_coords_style,
            up_vector=up_vector,
            multiple_bond_spacing=multiple_bond_spacing,
            include_save_buttons=include_save_buttons,
            objects=objects,
            graphics_class=graphics_class,
            cylinder_class=cylinder_class,
            cylinder_options=cylinder_options,
            sphere_class=sphere_class,
            sphere_options=sphere_options,
            arrow_class=arrow_class,
            arrow_options=arrow_options,
            line_class=line_class,
            line_options=line_options,
            disk_class=disk_class,
            disk_options=disk_options,
            animate=animate,
            recording_options=recording_options,
            animation_options=animation_options,
            jsmol_load_script=jsmol_load_script,
            include_jsmol_script_interface=include_jsmol_script_interface,
            dynamic_loading=dynamic_loading,
            label_style=label_style,
            theme_function=theme_function,
            plot_range_padding=plot_range_padding,
            annotation_function=annotation_function,
            # carried through so the main-path body can pop them like the original
            # method's locals (they aren't in the original full_opts dict):
            reconcile_bonds=reconcile_bonds,
            units=units,
            extra_opts=plot_ops
        )

        mode, backend = cls._resolve_plot_mode(mol, mode, backend, geometries, figure=figure, **full_opts)
        plotter_cls = cls._mode_registry.get(mode, Graphics3DMoleculePlotter)
        
        if len(geometries) == 0:
            geometries = mol.coords
        elif len(geometries) == 1:
            geometries = CoordinateSet(np.asanyarray(geometries[0]), mol.coords.system)
        else:
            geometries = CoordinateSet([np.asanyarray(g) for g in geometries], mol.coords.system)
            
        plotter = plotter_cls(
            mol, 
            geometries,
            mode=mode,
            backend=backend
        )
        full_opts = plotter._resolve_plot_theme(mode, backend, theme, full_opts)
        
        return plotter.plot(
            figure=figure,
            **full_opts
        )

    def plot(self, **new_opts):
        cur_opts = self.opts.copy()
        cur_opts.update(new_opts)
        return self.plot_impl(cur_opts)

    def plot_impl(self, full_opts):
        raise NotImplementedError(f"{type(self).__name__} doesn't implement plot()")


class JupyterMoleculePlotter(MoleculePlotter):
    modes = ('jupyter',)

    def jupyter_viz(self):
        from McUtils.Jupyter import MoleculeGraphics

        return MoleculeGraphics(self.atoms,
                                np.ndarray.view(self.coords.convert(CartesianCoordinates3D)),
                                bonds=self.bonds
                                )

    def plot_impl(self, cur_opts):
        return self.jupyter_viz(**cur_opts)


class JSMolMoleculePlotter(MoleculePlotter):
    modes = ('jsmol',)

    def jsmol_viz(self,
                  xyz=None,
                  animate=False,
                  vibrate=False,
                  script=None,
                  include_script_interface=False,
                  image_size=None,
                  width=None,
                  height=None,
                  figure=None,
                  **etc
                  ):
        from McUtils.Jupyter import JSMol
        if xyz is None:
            xyz = self.mol._format_xyz(0,
                                   len(self.atoms),
                                   self.atoms,
                                   self.coords * UnitsData.convert("BohrRadius", "Angstroms")
                                   )
        if image_size is not None:
            if nput.is_numeric(image_size):
                w = h = image_size
            else:
                w,h = image_size
            if width is None:
                width = w
            if height is None:
                height = h
        opts = {
            k: v
            for k, v in dict(
                width=width,
                height=height
            ).items()
            if v is not None
        }
        return JSMol.Applet(xyz, animate=animate, vibrate=vibrate, load_script=script,
                            include_script_interface=include_script_interface,
                            **opts,
                            **etc
                            )

    def plot_impl(self, full_opts):
        figure = full_opts.pop('figure', None)
        plot_ops = self._prep_jsmol_plot_opts(geometries=self.geometries, figure=figure, **full_opts)
        figure = self.jsmol_viz(**plot_ops)  # forwarded to self.mol.jsmol_viz
        self._set_backend_figure_options(figure, self.mode, self.backend, **full_opts)
        return figure

    @classmethod
    def _jsmol_atom_sel_block(cls, atom_offset, atoms, *exprs):
        block = []
        if not isinstance(atoms, str):
            atoms = "atomno=[{}]".format(
                ",".join(
                    str(atom_offset + h + 1)
                    if nput.is_int(h) else
                    h
                    for h in atoms
                ))
        block.append(f'select {atoms}')
        block.extend(exprs)
        block.append('select none')
        return block

    def _jsmol_view_settings(
            self,
            up_vector=None, right_vector=None, view_vector=None, view_distance=None,
            view_matrix=None, view_center=None):
        if view_matrix is None and (
                view_vector is not None
                or right_vector is not None
                or up_vector is not None
        ):
            if view_vector is None:
                if (
                        up_vector is not None and right_vector is not None
                ):
                    view_vector = nput.vec_crosses(up_vector, right_vector, normalize=True)
                elif right_vector is not None:
                    view_vector = nput.vec_crosses([0, 0, 1], right_vector, normalize=True)
                elif up_vector is not None:
                    view_vector = nput.vec_crosses(up_vector, [0, 1, 0], normalize=True)

            if view_vector is not None:
                m = nput.rotation_matrix(
                    view_vector,
                    [1, 0, 0]
                )
            else:
                m = np.eye(3)

            if up_vector is None and right_vector is not None:
                if view_vector is None:
                    view_vector = [1, 0, 0]
                up_vector = nput.vec_normalize(
                    nput.vec_crosses(right_vector, view_vector)
                )
            elif up_vector is not None and view_vector is not None:
                up_vector = nput.vec_crosses(
                    view_vector,
                    nput.vec_crosses(view_vector, up_vector),
                    normalize=True
                )
            if up_vector is not None:
                m = m @ nput.rotation_matrix(
                    m.T @ up_vector,
                    [1, 0, 0]
                )
            view_matrix = m

        ang, cross = nput.extract_rotation_angle_axis(view_matrix)
        ang = np.rad2deg(ang)
        if view_vector is None:
            view_vector = view_matrix[:, -1]
        return {
            'view_angle': ang,
            'rotation_axis': cross,
            'view': view_vector,
            'center': view_center,
            'dist': view_distance
        }
    def _prep_jsmol_load_script(self,
                                geom=None,
                                background=None,
                                atom_style=None,
                                bond_style=None,
                                use_default_radii=True,
                                atom_radius_scaling=None,
                                atom_radii=None,
                                radius_type=None,
                                bond_radius=None,
                                highlight_atoms=None,
                                highlight_bonds=None,
                                highlight_rings=None,
                                highlight_styles=None,
                                display_atom_numbers=None,
                                jsmol_load_script=None,
                                atom_offset=0,
                                use_default_bonds=True,
                                draw_coords=None,
                                draw_coords_style=None,
                                reflectiveness=None,
                                view_settings=None,
                                **ignored
                                ):
        if jsmol_load_script is not None:
            bits = jsmol_load_script.split(";") if isinstance(jsmol_load_script, str) else list(jsmol_load_script)
        else:
            bits = []
        if background is not None:
            if isinstance(background, str) and background.startswith("#"):
                background = "[{}]".format(background.replace("#", "x"))
            bits.append(f'background {background}')
        if not use_default_bonds:
            for b in self.bonds:
                i,j = b[:2]
                if len(b) > 2:
                    t = b[2]
                else:
                    t = 1
                if not isinstance(t, str):
                    if t < .9:
                        t = "AROMATIC"
                    elif t < 1.1:
                        t = "SINGLE"
                    elif t < 1.9:
                        t = 'AROMATIC'#"AROMATICDOUBLE"
                    elif t < 2.1:
                        t = "DOUBLE"
                    elif t < 2.9:
                        t = "AROMATICDOUBLE"
                    else:
                        t = "TRIPLE"
                bits.append(
                    f"connect (atomno={i + atom_offset + 1}) (atomno={j + atom_offset + 1}) {t}"
                )
        if highlight_styles is None:
            highlight_styles = self.highlight_styles
        if highlight_bonds is not None:
            if highlight_atoms is None: highlight_atoms = []
            highlight_atoms = set(highlight_atoms)
            for b in highlight_bonds:
                highlight_atoms.update(b)
            highlight_atoms = list(highlight_atoms)
        if highlight_rings is not None:
            if highlight_atoms is None: highlight_atoms = []
            highlight_atoms = set(highlight_atoms)
            for b in highlight_rings:
                highlight_atoms.update(b)
            highlight_atoms = list(highlight_atoms)
        if highlight_atoms is not None:
            glow_color = highlight_styles.get('glow')
            if glow_color is not None:
                bits.extend(
                    self._jsmol_atom_sel_block(
                        atom_offset,
                        highlight_atoms,
                        f'color {glow_color}',
                        'halos on'
                    )
                )
            else:
                color = highlight_styles.get('color', 'green')
                bits.extend(
                    self._jsmol_atom_sel_block(
                        atom_offset,
                        highlight_atoms,
                        f'color {color}',
                    )
                )
        if atom_style is not None:
            default_styles = {
                k:s
                for k,s in atom_style.items()
                if isinstance(k, str)
            }
            all_styles = {
                'all':default_styles
            } | {
                k:s
                for k,s in atom_style.items()
                if not isinstance(k, str)
            }
            non_empty = False
            for k,s in all_styles.items():
                if len(s) > 0:
                    non_empty = True
                    if not isinstance(k, str):
                        if not nput.is_int(k):
                            k = 'atomno=[{}]'.format(",".join(str(i+atom_offset+1) for i in k))
                        else:
                            k = f'@{atom_offset + k + 1}'
                    bits.append(f'select {k}')
                    for kk, v in s.items():
                        bits.append(f'{kk} {v}')
            if non_empty:
                bits.append('select none')
        if bond_style is not None:
            default_styles = {
                k:s
                for k,s in bond_style.items()
                if isinstance(k, str)
            }
            all_styles = {
                'all':default_styles
            } | {
                k:s
                for k,s in bond_style.items()
                if not isinstance(k, str)
            }
            non_empty = False
            for k,s in all_styles.items():
                if len(s) > 0:
                    non_empty = True
                    if not isinstance(k, str):
                        k = 'atomno=[{}]'.format(",".join(str(i+1+atom_offset) for i in k))
                    bits.append(f'select {k}')
                    for kk, v in s.items():
                        bits.append(f'{kk} bonds {v}')
            if non_empty:
                bits.append('select none')
        if not use_default_radii:
            if atom_radii is not None:
                if nput.is_numeric(atom_radii):
                    atom_radii = [atom_radii] * len(self._ats)
                elif dev.is_dict_like(atom_radii):
                    atom_radii = [
                        atom_radii.get(i,
                            atom_radii.get(a["ElementSymbol"])
                        ) for i,a in enumerate(self._ats)
                    ]
                if atom_radius_scaling is None:
                    atom_radius_scaling = [1] * len(atom_radii)
                elif nput.is_numeric(atom_radius_scaling):
                    atom_radius_scaling = [atom_radius_scaling] * len(atom_radii)

                non_empty = False
                for i, (r,s) in enumerate(zip(atom_radii, atom_radius_scaling)):
                    if r is not None:
                        non_empty = True
                        if nput.is_numeric(r) and nput.is_numeric(s):
                            r = r * s
                        bits.extend(
                            self._jsmol_atom_sel_block(
                                atom_offset,
                                f'@{atom_offset + i + 1}',
                                f'spacefill {r}'
                            )
                        )
                if non_empty:
                    bits.append('select none')
            elif atom_radius_scaling is not None:
                if nput.is_numeric(atom_radius_scaling):
                    bits.extend(
                        self._jsmol_atom_sel_block(
                            atom_offset,
                            'all',
                            f'spacefill {100*atom_radius_scaling:.0f}%'
                        )
                    )
                else:
                    non_empty = False
                    for i,r in enumerate(atom_radius_scaling):
                        if r is not None:
                            non_empty=True
                            bits.extend(
                                self._jsmol_atom_sel_block(
                                    atom_offset,
                                    f'@{atom_offset + i+1}',
                                    f'spacefill {100 * r:.0f}%'
                                )
                            )
                    if non_empty:
                        bits.append('select none')
            if bond_radius is not None:
                bits.extend(
                    self._jsmol_atom_sel_block(
                        atom_offset,
                        'all',
                        f'wireframe {bond_radius}'
                    )
                )
        if display_atom_numbers:
            if display_atom_numbers is True:
                bits.extend(
                    self._jsmol_atom_sel_block(
                        atom_offset,
                        'all',
                        'label %i'
                    )
                )
            else:
                bits.extend(
                    self._jsmol_atom_sel_block(
                        atom_offset,
                        display_atom_numbers,
                        'label %i'
                    )
                )

        if reflectiveness is not None:
            if reflectiveness is True:
                bits.append('set specular on')
            elif reflectiveness is False:
                bits.append('set specular off')
            else:
                bits.append(f'set specular {100*reflectiveness:.0f}')

        if draw_coords is not None:
            radii = self._get_atom_radii(atom_radii, atom_radius_scaling, radius_type)
            draw_coords, draw_coords_style = self._prep_display_draw_coords(draw_coords, draw_coords_style)
            no_hover = False
            draw_font = None
            for n, (k, v) in enumerate(draw_coords.items()):
                if isinstance(v, str):
                    v = {'label': v}
                props = dict(draw_coords_style, **v)
                id = None
                style = None
                label_props = None
                if len(k) == 2:
                    id = f'bondline{n}'
                    points, label_props, style, label_style = self._prep_draw_line(
                        geom, k, props,
                        default_label_style=self.draw_coords_label_style,
                        line_class=None,
                        theme_function=None,
                        radii=radii,
                        plotos={},
                        line_options={}
                    )
                    points = " ".join(f"{{{x:.3f} {y:.3f} {z:.3f}}}" for x, y, z in points)
                    bits.append(f'draw ID {id} LINE {points}')
                elif len(k) == 3:
                    id = f'bondang{n}'
                    arc_props, label_props, style, label_style = self._prep_draw_arc(
                        geom, k, props,
                        up_vector=[0, 0, 1],
                        default_label_style=self.draw_coords_label_style,
                        disk_class=None,
                        theme_function=None,
                        radii=radii,
                        plotos={},
                        disk_options={}
                    )
                    if 'line_color' in style:
                        style['color'] = style.get('color', style['line_color'])
                    yy, axes, _, radius = arc_props
                    density = style.pop('angular_density', None)
                    normal, offset_angle, span_angle = nput.angle_arc_parameters(*axes)
                    points = nput.arc_points(yy, radius, offset_angle, span_angle, normal=normal, angular_density=density)
                    points = " ".join(f"{{{x:.3f} {y:.3f} {z:.3f}}}" for x,y,z in points)
                    bits.append(f'draw ID {id} CURVE {points}')
                else:
                    raise NotImplementedError(...)

                if id is not None and style is not None:
                    color = style.pop('color', None)
                    if color is not None:
                        bits.append(f'color ${id} {color}')

                if label_props is not None:
                    label, label_props = label_props[0], label_props[1:]
                    if label is not None:
                        # if not no_hover:
                        #     bits.append('set drawHover False')
                        #     no_hover = True
                        if nput.is_numeric(label_props[0]):
                            points = label_props
                        else:
                            points = label_props[0]
                        x,y,z = points
                        points = f"{{{x:.3f} {y:.3f} {z:.3f}}}"
                        title_command = f'TITLE "{label}"'
                        label_color = label_style.get('color', label_style.get('font_color'))
                        if label_color is not None:
                            title_command = f'{title_command} TITLE COLOR {label_color}'
                        font_command = ""
                        font_size = label_style.get('font_size')
                        if font_size is not None:
                            font_command = f'{font_command} {font_size}'.strip()
                            label_style['font_family'] = label_style.get('font_family', "sansserif")
                            label_style['font_type'] = label_style.get('font_type', "plain")
                        font_face = label_style.get('font_family')
                        if font_face is not None:
                            font_command = f'{font_command} {font_face}'.strip()
                        font_type = label_style.get('font_type')
                        if font_type is not None:
                            font_command = f'{font_command} {font_type}'.strip()
                        bits.append(f'draw ID {id}label LINE [{points} {points}] {title_command} DIAMETER .01')
                        if len(font_command) > 0:
                            draw_font = font_command
            if draw_font is not None:
                #JSMol is broken and its documentation is wrong, only one font style is supported
                bits.append(f'FONT draw {draw_font}')

        if view_settings is not None:
            vps = self._jsmol_view_settings(**view_settings)
            center = vps.pop('center')
            if center is not None:
                x, y, z = center
                bits.append(f"translate {{{x:.3f} {y:.3f} {z:.3f}}}")
            else:
                center = np.tensordot(self.masses, geom, axes=[0, 0]) / np.sum(self.masses)
                x, y, z = -center
                bits.append(f"translate {{{x:.3f} {y:.3f} {z:.3f}}}")
            center = np.asanyarray(center)
            ang = vps.pop('view_angle')
            ax = vps.pop('rotation_axis')
            x, y, z = center
            x2, y2, z2 = center + ax
            bits.append(f"rotate {{{x:.3f} {y:.3f} {z:.3f}}} {{{x2:.3f} {y2:.3f} {z2:.3f}}} {ang:.3f}")
            dist = vps.pop('dist')
            if dist is not None:
                view_geom = (geom - center[np.newaxis]) @ nput.rotation_matrix(ax, ang)
                x = view_geom[:, 0]
                x_range = np.max(x) - np.min(x)
                y = view_geom[:, 0]
                y_range = np.max(y) - np.min(y)
                default_dist = np.max([x_range, y_range])
                if use_default_radii:
                    if atom_radii is not None:
                        if nput.is_numeric(atom_radii):
                            atom_radii = [atom_radii] * len(self._ats)
                        elif dev.is_dict_like(atom_radii):
                            atom_radii = [
                                atom_radii.get(i,
                                               atom_radii.get(a["ElementSymbol"])
                                               ) for i, a in enumerate(self._ats)
                            ]
                    default_dist = default_dist + np.max(atom_radii)
                target_zoom = default_dist / np.sqrt(dist)
                bits.append(f"zoom {100*target_zoom:.0f}")

        return bits

    def _prep_jsmol_plot_opts(self, geometries=None, figure=None, atom_offset=0,
                              extra_opts=None, script=None, jsmol_load_script=None,
                              draw_coords=None, recording_options=None,
                              dynamic_loading=None, include_jsmol_script_interface=None,
                              include_save_buttons=None, **etc):
        if extra_opts is None:
            extra_opts = {}
        xyz = extra_opts.pop('xyz', None)
        prev_script = []
        if figure is not None:
            model_state = figure.model_file
            nats, base_xyz = model_state.split("\n", 1)
            nats = int(nats)
            prev_script = figure.load_script
            etc = dict(figure.attrs, **etc)
            atom_offset = atom_offset + nats
            xyz = f"{nats + len(self.atoms)}\n" + base_xyz + "\n" + self.mol._format_xyz(
                0, len(self.atoms), self.atoms,
                self.coords * UnitsData.convert("BohrRadius", "Angstroms")
            ).split("\n", 2)[-1]
        elif xyz is None and geometries is not None and len(geometries) > 0:
            geometries = np.asanyarray(geometries)
            if geometries.ndim == 3 and geometries.shape[0] == 1:
                geometries = geometries[0]
            if geometries.ndim > 2:
                raise ValueError(f"Jmol can't handle {geometries.shape}")
            xyz = self.mol._format_xyz(
                0, len(self.atoms), self.atoms,
                geometries * UnitsData.convert("BohrRadius", "Angstroms")
            )

        use_default_radii = extra_opts.pop('use_default_radii', True)
        use_default_bonds = extra_opts.pop('use_default_bonds', True)
        view_settings = extra_opts.pop('view_settings', None)
        background = extra_opts.get('background')

        include_script_interface = extra_opts.pop('include_script_interface', None)
        if include_script_interface is None:
            if include_jsmol_script_interface is None:
                include_jsmol_script_interface = include_save_buttons
            include_script_interface = include_jsmol_script_interface

        if script is None:
            if geometries is None:
                geometries = self.coords
            geometries = geometries * UnitsData.convert("BohrRadius", "Angstroms")
            script = prev_script + self._prep_jsmol_load_script(
                geom=geometries, atom_offset=atom_offset,
                jsmol_load_script=jsmol_load_script, background=background,
                use_default_radii=use_default_radii, use_default_bonds=use_default_bonds,
                draw_coords=draw_coords, view_settings=view_settings, **etc
            )

        return dict(
            xyz=xyz, script=script, recording_options=recording_options,
            dynamic_loading=dynamic_loading, autobond=use_default_bonds,
            include_script_interface=include_script_interface, **extra_opts
        )


class RDKitMoleculePlotter(MoleculePlotter):
    modes = ('rdkit', 'rdkit3d')

    def plot_impl(self, full_opts):
        figure = full_opts.pop('figure', None)
        if self.mode == 'rdkit3d':
            plot_opts = self._prep_rdkit_plot_opts(figure=figure, **full_opts)
            figure = self.mol.rdmol.plot(**plot_opts)   # forwarded to self.mol.rdmol
        else:
            draw_opts = self._prep_rdkit_draw_opts(figure=figure, **full_opts)
            figure = self.mol.rdmol.draw(**draw_opts)
        self._set_backend_figure_options(figure, self.mode, self.backend, **full_opts)
        return figure

    def _prep_rdkit_plot_opts(self, extra_opts=None, **etc):
        if extra_opts is None:
            extra_opts = {}
        return extra_opts

    def _prep_rdkit_draw_opts(self, figure=None, atom_style=None, bond_style=None,
                              atom_radii=None, atom_radius_scaling=None, radius_type=None,
                              bond_radius=None, highlight_atoms=None, highlight_atom_radii=None,
                              highlight_bonds=None, highlight_bond_radii=None, highlight_color=None,
                              highlight_rings=None, highlight_styles=None, draw_coords=None,
                              extra_opts=None, include_save_buttons=None,
                              display_atom_numbers=None, label_style=None, **etc):
        if extra_opts is None:
            extra_opts = {}
        use_default_radii = extra_opts.pop('use_default_radii', True)
        if highlight_styles is None:
            highlight_styles = self.highlight_styles

        highlight_atom_colors = None
        highlight_bond_colors = None
        if highlight_color is None and (highlight_atoms is not None or highlight_bonds is not None):
            glow = highlight_styles.get('glow')
            color = highlight_styles.get('color')
            if glow is not None:
                if color is None:
                    if isinstance(glow, str):
                        glow = np.array(plt.ColorPalette.parse_color_string(glow)) / 255
                    if len(glow) == 3:
                        glow = tuple(glow) + (.5,)
                    color = glow
                else:
                    if isinstance(glow, str):
                        glow = plt.ColorPalette.parse_color_string(glow)
                    if isinstance(color, str):
                        color = plt.ColorPalette.parse_color_string(color)
                    color = plt.prep_color(palette=[glow, color], blending=0.5, return_color_code=False)
                    color = np.array(color)
            highlight_color = color
        if highlight_atoms is not None:
            highlight_atom_colors = {}
            for a in highlight_atoms:
                highlight_atom_colors[a] = highlight_color
        if highlight_bonds is not None:
            highlight_bond_colors = {}
            for b in highlight_bonds:
                if not nput.is_int(b):
                    b = tuple(b)
                highlight_bond_colors[b] = highlight_color

        if not use_default_radii:
            atom_radii = self._get_atom_radii(atom_radii, atom_radius_scaling, radius_type)
        if use_default_radii:
            atom_radii = None
            bond_radius = None

        return dict(
            dict(
                figure=figure,
                highlight_atoms=highlight_atoms,
                highlight_bonds=highlight_bonds,
                highlight_rings=highlight_rings,
                highlight_atom_radii=highlight_atom_radii,
                highlight_bond_radii=highlight_bond_radii,
                atom_radii=atom_radii,
                bond_radius=bond_radius,
                highlight_atom_colors=highlight_atom_colors,
                highlight_bond_colors=highlight_bond_colors,
                include_save_buttons=include_save_buttons,
                display_atom_numbers=display_atom_numbers,
                draw_coords=draw_coords,
                label_style=label_style,
                highlight_bond_width_multiplier=None
            ),
            **extra_opts
        )


class Graphics3DMoleculePlotter(MoleculePlotter):
    """
    The main path (`x3d` / `matplotlib3D` / `plotly3D` / `svg3D`). This is the
    original `plot` body from the `extra_opts` split onward, carried over as-is;
    `self.<helper>` resolves on the plotter and `self.<mol attr>` is forwarded to
    the molecule, so nothing in the loop needed rewriting.
    """
    modes = ('x3d', 'matplotlib3D', 'plotly3D', 'svg3D', 'fast')

    def plot_impl(self, full_opts):
        from McUtils.Plots import Graphics3D, Sphere, Cylinder, Arrow, Line, Disk

        figure = full_opts.pop('figure', None)
        backend = self.backend
        mode = self.mode
        geometries = self.geometries

        # these two are ordinary locals in the original method; here they ride in
        # on full_opts, so pop them before the shared unpack below
        reconcile_bonds = full_opts.pop('reconcile_bonds', True)
        units = full_opts.pop('units', 'Angstroms')

        plot_ops = full_opts.pop('extra_opts')
        graphics_keys = Graphics3D.known_keys | Graphics3D.opt_keys | Graphics3D.figure_keys | Graphics3D.axes_keys
        graphics_opts = {k: plot_ops[k] for k in plot_ops.keys() & graphics_keys}
        plot_ops = {k: plot_ops[k] for k in plot_ops.keys() - graphics_keys}

        if 'background' not in graphics_opts:
            graphics_opts['background'] = 'transparent'

        full_opts_base = full_opts.copy()

        (
            return_objects,
            bonds,
            bond_radius,
            atom_radius_scaling,
            atom_style,
            atom_radii,
            atom_text,
            display_atom_numbers,
            radius_type,
            bond_style,
            capped_bonds,
            reflectiveness,
            vector_style,
            highlight_atoms,
            highlight_bonds,
            highlight_rings,
            highlight_styles,
            comparison_styles,
            animation_frame_styles,
            mode_vectors,
            mode_vector_origins,
            mode_vector_origin_mode,
            mode_vector_display_cutoff,
            principle_axes,
            principle_axes_origin,
            principle_axes_origin_mode,
            principle_axes_style,
            dipole,
            dipole_origin,
            dipole_origin_mode,
            render_multiple_bonds,
            render_fractional_bonds,
            fractional_bond_offset,
            bond_center_radius_offset,
            draw_coords,
            draw_coords_style,
            up_vector,
            multiple_bond_spacing,
            include_save_buttons,
            objects,
            graphics_class,
            cylinder_class,
            cylinder_options,
            sphere_class,
            sphere_options,
            arrow_class,
            arrow_options,
            line_class,
            line_options,
            disk_class,
            disk_options,
            animate,
            recording_options,
            animation_options,
            jsmol_load_script,
            include_jsmol_script_interface,
            dynamic_loading,
            label_style,
            plot_range_padding,
            theme_function,
            annotation_function
        ) = [
            full_opts.pop(f) for f in (  # allows for theme flexibility
                "return_objects",
                "bonds",
                "bond_radius",
                "atom_radius_scaling",
                "atom_style",
                "atom_radii",
                "atom_text",
                "display_atom_numbers",
                "radius_type",
                "bond_style",
                "capped_bonds",
                "reflectiveness",
                "vector_style",
                "highlight_atoms",
                "highlight_bonds",
                "highlight_rings",
                "highlight_styles",
                "comparison_styles",
                "animation_frame_styles",
                "mode_vectors",
                "mode_vector_origins",
                "mode_vector_origin_mode",
                "mode_vector_display_cutoff",
                "principle_axes",
                "principle_axes_origin",
                "principle_axes_origin_mode",
                "principle_axes_style",
                "dipole",
                "dipole_origin",
                "dipole_origin_mode",
                "render_multiple_bonds",
                "render_fractional_bonds",
                "fractional_bond_offset",
                "bond_center_radius_offset",
                "draw_coords",
                "draw_coords_style",
                "up_vector",
                "multiple_bond_spacing",
                "include_save_buttons",
                "objects",
                "graphics_class",
                "cylinder_class",
                "cylinder_options",
                "sphere_class",
                "sphere_options",
                "arrow_class",
                "arrow_options",
                "line_class",
                "line_options",
                "disk_class",
                "disk_options",
                "animate",
                "recording_options",
                "animation_options",
                "jsmol_load_script",
                "include_jsmol_script_interface",
                "dynamic_loading",
                "label_style",
                'plot_range_padding',
                "theme_function",
                "annotation_function"
            )
        ]
        if len(full_opts) > 0:
            raise ValueError(f"options unhandled: {full_opts}")

        if cylinder_class is None:
            cylinder_class = Cylinder
        if cylinder_options is None:
            cylinder_options = {}
        if sphere_class is None:
            sphere_class = Sphere
        if sphere_options is None:
            sphere_options = {}
        if arrow_class is None:
            arrow_class = Arrow
        if line_class is None:
            line_class = Line
        if line_options is None:
            line_options = {}
        if disk_class is None:
            disk_class = Disk
        if disk_options is None:
            disk_options = {}

        if backend == 'matplotlib3D':
            if mode != 'quality':
                plot_ops['rendering'] = plot_ops.get('rendering', 'flat')
                cylinder_options['edge_color'] = cylinder_options.get('edge_color', [.3] * 3)
                sphere_options['edge_color'] = sphere_options.get('edge_color', [.3] * 3)

        if units is not None:
            geometries = geometries * UnitsData.convert("BohrRadius", units)

        if geometries.ndim == 2:
            geometries = geometries[np.newaxis]
        if animate is None:
            animate = geometries.shape[0] > 1 and comparison_styles is None

        dipole, dipole_origin = self._prep_display_dipole(geometries, dipole, dipole_origin, units)
        principle_axes, principle_axes_origin, principle_axes_style = self._prep_principle_axes(
            geometries, units, principle_axes, principle_axes_origin, principle_axes_style
        )
        mode_vectors, mode_vector_origins = self._prep_display_mode_vectors(
            geometries, units,
            mode_vectors, mode_vector_origins
        )

        draw_coords, draw_coords_style = self._prep_display_draw_coords(
            draw_coords, draw_coords_style
        )

        geometries = geometries.convert(CartesianCoordinates3D)
        draw_bonds = bonds

        radii = self._get_atom_radii(atom_radii, atom_radius_scaling, radius_type)

        bonds = [None] * len(geometries)
        atoms = [None] * len(geometries)
        arrows = [None] * len(geometries)

        if vector_style is None:
            vector_style = {}
        vector_style = dict(self.vector_style, **vector_style)

        if label_style is None:
            label_style = {}

        atom_text = self._prep_display_atom_text(atom_text, display_atom_numbers, label_style)

        if highlight_styles is None:
            highlight_styles = self.highlight_styles

        if highlight_atoms is True:
            highlight_atoms = list(range(len(self.atoms)))
        if highlight_rings is not None:
            if highlight_atoms is None:
                highlight_atoms = []
            if highlight_bonds is None:
                highlight_bonds = []
            highlight_atoms = list(highlight_atoms)
            highlight_bonds = list(highlight_bonds)
            for r in highlight_rings:
                highlight_atoms.extend(r)
                highlight_bonds.extend(zip(
                    r, r[1:] + r[:1]
                ))

        atom_style, highlight_atoms, colors, glows = self._prep_display_atom_style(
            atom_style,
            highlight_atoms,
            backend=backend,
            reflectiveness=reflectiveness,
            highlight_styles=highlight_styles
        )

        og_highlight_bonds = list(highlight_bonds) if highlight_bonds is not None else None
        bond_style, highlight_bonds = self._prep_display_bond_style(
            bond_style,
            highlight_bonds,
            backend=backend,
            reflectiveness=reflectiveness,
            highlight_atoms=highlight_atoms,
            highlight_styles=highlight_styles,
            capped_bonds=capped_bonds
        )

        if draw_bonds is None:
            draw_bonds = self.bonds
        elif dev.str_is(draw_bonds, "recompute"):
            if units is not None:
                conv = UnitsData.convert(units, "BohrRadius")
            else:
                conv = 1
            draw_bonds = [
                self.mol.modify(coords=g, bonds=None).bonds
                for g in geometries * conv
            ]

        if draw_bonds is None or draw_bonds is False:
            draw_bonds = [None] * len(geometries)
        elif nput.is_int(draw_bonds[0][0]):
            draw_bonds = [draw_bonds] * len(geometries)

        max_bond_orders = {}
        if reconcile_bonds and any(b is not None for b in draw_bonds):
            for bond_set in draw_bonds:
                if bond_set is None: continue
                for b in bond_set:
                    if b is None: continue
                    if len(b) == 2:
                        b = (b[0], b[1], 1)
                    key = tuple(sorted(b[:2]))
                    max_bond_orders[key] = max(max_bond_orders.get(key, 0), b[2])
            new_bonds = []
            for bond_set in draw_bonds:
                subset = {}
                if bond_set is not None:
                    for b in bond_set:
                        if b is None: continue
                        if len(b) == 2:
                            b = (b[0], b[1], 1)
                        key = tuple(sorted(b[:2]))
                        subset[key] = b[2]
                new_bonds.append(subset)
            draw_bonds = new_bonds

        default_label_style = label_style | self.draw_coords_label_style
        radii = np.asanyarray(radii)

        # TODO: create overload-able function for this
        if backend == 'matplotlib3D':
            graphics_opts['box_ratios'] = graphics_opts.get('box_ratios', 'auto')
            graphics_opts['aspect_ratio'] = graphics_opts.get('aspect_ratio', 'equal')
            graphics_opts['frame'] = graphics_opts.get('frame', False)
            pr = graphics_opts.get('plot_range', None)
            graphics_opts['plot_range'] = self._default_plot_range(geometries, pr, plot_range_padding,
                                                                   atom_radius_scaling * radii)
            graphics_opts['autoscale'] = graphics_opts.get('autoscale', False)
        elif backend == 'plotly3D':
            graphics_opts['box_ratios'] = graphics_opts.get('box_ratios', 'cube')
            graphics_opts['frame'] = graphics_opts.get('frame', False)
            pr = graphics_opts.get('plot_range', None)
            graphics_opts['plot_range'] = self._default_plot_range(geometries, pr, plot_range_padding,
                                                                   atom_radius_scaling * radii)
            vs = graphics_opts.pop('view_settings', None)
            if vs is not None:
                view_dist = vs.pop('view_distance', None)
                if view_dist is not None:
                    pr = graphics_opts['plot_range']
                    # (xmin, xmax), (ymin, ymax), (zmin, zmax) = pr
                    # compute the portion of the field of view we actually want to sample
                    view_percents = view_dist / np.max([M - m for m, M in pr])
                    vs['view_distance'] = np.sqrt(view_percents)  # * np.sqrt(25 / 16 * 3)
                graphics_opts['view_settings'] = vs

        if figure is None:
            if graphics_class is None:
                graphics_class = Graphics3D
            figure = graphics_class(backend=backend, **graphics_opts)

        self._set_backend_figure_options(figure, mode, backend, **full_opts_base)

        if animation_frame_styles is not None:
            comparison_styles = animation_frame_styles
        if comparison_styles is None:
            comparison_styles = [None] * len(geometries)
        elif isinstance(comparison_styles, dict):
            comparison_styles = [None, comparison_styles]
        if len(comparison_styles) < len(geometries):
            comparison_styles = (list(comparison_styles) * len(geometries))[:len(geometries)]

        global_state = (
            atom_style,
            radii,
            atom_radius_scaling,
            atom_radii,
            radius_type,
            bond_style,
            bond_center_radius_offset,
            multiple_bond_spacing,
            render_multiple_bonds,
            render_fractional_bonds,
            fractional_bond_offset,
            up_vector,
            cylinder_class,
            colors,
            glows,
            display_atom_numbers,
            atom_text,
            label_style,
            theme_function,
            cylinder_options,
            sphere_class,
            sphere_options,
            draw_coords
        )

        for i, geom in enumerate(geometries):
            plotos = plot_ops.copy()

            (
                atom_style,
                radii,
                atom_radius_scaling,
                atom_radii,
                radius_type,
                bond_style,
                bond_center_radius_offset,
                multiple_bond_spacing,
                render_multiple_bonds,
                render_fractional_bonds,
                fractional_bond_offset,
                up_vector,
                cylinder_class,
                colors,
                glows,
                display_atom_numbers,
                atom_text,
                label_style,
                theme_function,
                cylinder_options,
                sphere_class,
                sphere_options,
                draw_coords
            ) = global_state

            substyle = comparison_styles[i]
            if substyle is not None:
                a_sty = substyle.pop('atom_style', atom_style)

                ha2 = substyle.pop('highlight_atoms', highlight_atoms)
                hs = substyle.pop('highlight_styles', highlight_styles)
                if ha2 is None: ha2 = []
                highlight_diffs = np.setdiff1d(highlight_atoms, ha2)
                a_sty = a_sty.copy()
                for d in highlight_diffs:
                    if d in a_sty:
                        subd = a_sty[d]
                        for k, v in highlight_styles.items():
                            if subd.get(k) == v:
                                subd.pop(k)

                atom_style, updata_ha, colors, glows = self._prep_display_atom_style(
                    a_sty,
                    ha2,
                    backend=backend,
                    reflectiveness=reflectiveness,
                    highlight_styles=hs,
                )

                atom_radius_scaling = substyle.pop('atom_radius_scaling', atom_radius_scaling)
                atom_radii = substyle.pop('atom_radii', atom_radii)
                radius_type = substyle.pop('radius_type', radius_type)
                radii = self._get_atom_radii(atom_radii, atom_radius_scaling, radius_type)

                atom_text = substyle.pop('atom_text', atom_text)
                display_atom_numbers = substyle.pop('display_atom_numbers', display_atom_numbers)
                label_style = substyle.pop('label_style', label_style)
                atom_text = self._prep_display_atom_text(atom_text, display_atom_numbers, label_style)

                bond_center_radius_offset = substyle.pop('bond_center_radius_offset', bond_center_radius_offset)
                multiple_bond_spacing = substyle.pop('multiple_bond_spacing', multiple_bond_spacing)
                render_multiple_bonds = substyle.pop('render_multiple_bonds', render_multiple_bonds)
                render_fractional_bonds = substyle.pop('render_fractional_bonds', render_fractional_bonds)
                fractional_bond_offset = substyle.pop('fractional_bond_offset', fractional_bond_offset)

                b_sty = substyle.pop('bond_style', bond_style)
                hb2 = substyle.pop('highlight_bonds', og_highlight_bonds)
                if b_sty is not None:
                    if hb2 is None: hb2 = []
                    if highlight_bonds is None: highlight_bonds = []
                    highlight_diffs = [h for h in highlight_bonds if h not in hb2]

                    b_sty = b_sty.copy()
                    for d in highlight_diffs:
                        if d in b_sty:
                            subd = b_sty[d]
                            for k, v in highlight_styles.items():
                                if subd.get(k) == v:
                                    subd.pop(k)
                    bond_style, highlight_bonds = self._prep_display_bond_style(
                        b_sty,
                        hb2,
                        backend=backend,
                        reflectiveness=reflectiveness,
                        highlight_atoms=highlight_atoms,
                        highlight_styles=highlight_styles,
                        capped_bonds=capped_bonds
                    )

                draw_coords = substyle.pop('draw_coords', draw_coords)
                theme_function = substyle.pop('theme_function', theme_function)
                sphere_class = substyle.pop('sphere_class', sphere_class)
                sphere_options = substyle.pop('sphere_options', sphere_options)
                cylinder_class = substyle.pop('cylinder_class', cylinder_class)
                cylinder_options = substyle.pop('cylinder_options', cylinder_options)

                plotos = plotos | substyle

            if backend in {'matplotlib3D', 'plotly3D'}:
                box_scalings = plotos.get('box_scalings')
                pr = graphics_opts.get('plot_range')
                if box_scalings is None:
                    if pr is None:
                        pr = [None, None, None]
                    rx, ry, rz = pr
                    if rx is None:
                        min_x = np.min(geom[:, 0] - radii)
                        max_x = np.max(geom[:, 0] + radii)
                    else:
                        min_x, max_x = rx
                    if ry is None:
                        min_y = np.min(geom[:, 1] - radii)
                        max_y = np.max(geom[:, 1] + radii)
                    else:
                        min_y, max_y = ry
                    if rz is None:
                        min_z = np.min(geom[:, 2] - radii)
                        max_z = np.max(geom[:, 2] + radii)
                    else:
                        min_z, max_z = rz
                    box_scalings = (
                            (np.max(figure.image_size) / 72) /
                            np.array([max_x - min_x, max_y - min_y, max_z - min_z])
                    )
                    plotos['box_scalings'] = box_scalings

            bond_list = draw_bonds[i]
            if bond_style is not False and bond_list is not None:
                bonds: list
                if reconcile_bonds:
                    # no reason to redo this so often but w/e
                    sublist = [
                        [i, j, bond_list.get((i, j), 0)]
                        for (i, j), t in max_bond_orders.items()
                    ]
                    bond_list = sorted(sublist, key=lambda b: (
                            self._ats[b[0]]["ElementSymbol"] in {"H", "X", "D"}
                            or
                            self._ats[b[1]]["ElementSymbol"] in {"H", "X", "D"}
                    ))
                else:
                    bond_list = sorted(bond_list, key=lambda b: (
                            self._ats[b[0]]["ElementSymbol"] in {"H", "X", "D"}
                            or
                            self._ats[b[1]]["ElementSymbol"] in {"H", "X", "D"}
                    ))
                bond_objs = self._get_bondlist_primitives(
                    geom,
                    bond_list,
                    bond_radius=bond_radius,
                    radii=radii,
                    bond_center_radius_offset=bond_center_radius_offset,
                    multiple_bond_spacing=multiple_bond_spacing,
                    render_multiple_bonds=render_multiple_bonds,
                    render_fractional_bonds=render_fractional_bonds,
                    fractional_bond_offset=fractional_bond_offset,
                    max_bond_orders=max_bond_orders,
                    up_vector=up_vector,
                    cylinder_class=cylinder_class,
                    colors=colors,
                    glows=glows,
                    bond_style=bond_style,
                    theme_function=theme_function,
                    plotos=plotos,
                    cylinder_options=cylinder_options
                )

                if not objects:
                    _ = []
                    for cset in bond_objs:
                        subbonds = []
                        for c in cset:
                            cyl = c.plot(figure)
                            if isinstance(cyl, (list, tuple)):
                                cyl = cyl[0]
                            subbonds.append(cyl)
                        _.append(subbonds)
                    bond_objs = _
                bonds[i] = bond_objs

            if atom_style is not False:
                atoms: list
                atoms[i] = [None] * len(geom)
                atom_objs = self._get_atom_primitives(
                    geom,
                    self._ats,
                    colors=colors,
                    radii=radii,
                    sphere_class=sphere_class,
                    atom_style=atom_style,
                    theme_function=theme_function,
                    plotos=plotos,
                    sphere_options=sphere_options
                )

                if not objects:
                    _ = []
                    for c in atom_objs:
                        cyl = c.plot(figure)
                        if isinstance(cyl, (list, tuple)):
                            cyl = cyl[0]
                        _.append(cyl)
                    atom_objs = _
                atoms[i] = atom_objs

            if dipole is not None:
                arrows: list
                if arrows[i] is None: arrows[i] = []
                dip = dipole[i]
                dip_obs = self._get_dipole_primitives(
                    geom,
                    dip,
                    dipole_origin=dipole_origin[i] if dipole_origin is not None else dipole_origin,
                    dipole_origin_mode=dipole_origin_mode,
                    mode_vector_display_cutoff=mode_vector_display_cutoff,
                    arrow_class=arrow_class,
                    theme_function=theme_function,
                    plotos=plotos,
                    vector_style=vector_style
                )
                if not objects:
                    _ = []
                    for c in dip_obs:
                        cyl = c.plot(figure)
                        if isinstance(cyl, (list, tuple)):
                            cyl = cyl[0]
                        _.append(cyl)
                    dip_obs = _
                arrows[i].extend(dip_obs)

            if principle_axes is not None:
                arrows: list
                if arrows[i] is None: arrows[i] = []
                pax_objs = self._get_pax_primitives(
                    geom,
                    principle_axes[i],
                    principle_axes_origin=principle_axes_origin,
                    principle_axes_origin_mode=principle_axes_origin_mode,
                    principle_axes_style=principle_axes_style,
                    arrow_class=arrow_class,
                    theme_function=theme_function,
                    plotos=plotos,
                    vector_style=vector_style
                )

                if not objects:
                    _ = []
                    for c in pax_objs:
                        cyl = c.plot(figure)
                        if isinstance(cyl, (list, tuple)):
                            cyl = cyl[0]
                        _.append(cyl)
                    pax_objs = _
                arrows[i].extend(pax_objs)

            if mode_vectors is not None:
                arrows: list
                if arrows[i] is None: arrows[i] = []
                mode_arrows = self._get_mode_vector_primitives(
                    geom,
                    mode_vectors[i],
                    mode_vector_display_cutoff=mode_vector_display_cutoff,
                    mode_vector_origins=mode_vector_origins[
                        i] if mode_vector_origins is not None else mode_vector_origins,
                    mode_vector_origin_mode=mode_vector_origin_mode,
                    arrow_class=arrow_class,
                    theme_function=theme_function,
                    plotos=plotos,
                    vector_style=vector_style
                )

                if not objects:
                    _ = []
                    for c in mode_arrows:
                        cyl = c.plot(figure)
                        if isinstance(cyl, (list, tuple)):
                            cyl = cyl[0]
                        _.append(cyl)
                    mode_arrows = _
                arrows[i].extend(mode_arrows)

            if draw_coords is not None:
                arrows: list
                if arrows[i] is None: arrows[i] = []
                for k, v in draw_coords.items():
                    if isinstance(v, str):
                        v = {'label': v}
                    v = dict(draw_coords_style, **v)
                    if len(k) == 2:
                        prims = self._get_draw_coords_line(
                            geom,
                            k,
                            v,
                            default_label_style=default_label_style,
                            radii=radii,
                            line_class=line_class,
                            theme_function=theme_function,
                            plotos=plotos,
                            line_options=line_options
                        )

                    elif len(k) == 3:
                        prims = self._get_draw_coords_arc(
                            geom,
                            k,
                            v,
                            up_vector=up_vector,
                            default_label_style=default_label_style,
                            radii=radii,
                            disk_class=disk_class,
                            theme_function=theme_function,
                            plotos=plotos,
                            disk_options=disk_options
                        )
                    else:
                        prims = self._get_draw_coords_dihed(
                            geom,
                            k,
                            v,
                            default_label_style=default_label_style,
                            disk_class=disk_class,
                            # radii=radii,
                            theme_function=theme_function,
                            plotos=plotos,
                            disk_options=disk_options
                        )

                    if not objects:
                        _ = []
                        for c in prims:
                            cyl = c.plot(figure)
                            if isinstance(cyl, (list, tuple)):
                                cyl = cyl[0]
                            _.append(cyl)
                        prims = _
                    arrows[i].extend(prims)

            if atom_text is not None:
                arrows: list
                if arrows[i] is None: arrows[i] = []
                prims = self._get_atom_text_primitives(
                    geom,
                    atom_text,
                    radii=radii,
                    plotos=plotos
                )
                if not objects:
                    _ = []
                    for c in prims:
                        cyl = c.plot(figure)
                        if isinstance(cyl, (list, tuple)):
                            cyl = cyl[0]
                        _.append(cyl)
                    prims = _
                arrows[i].extend(prims)

            if annotation_function is not None:
                arrows: list
                if arrows[i] is None: arrows[i] = []
                annotations = annotation_function(self, i, geom)
                if not objects:
                    _ = []
                    for c in annotations:
                        cyl = c.plot(figure)
                        if isinstance(cyl, (list, tuple)):
                            cyl = cyl[0]
                        _.append(cyl)
                    annotations = _
                arrows[i].extend(annotations)

        if animate:
            if animation_options is None: animation_options = {}
            figure = figure.animate_frames(
                [
                    (a if a is not None else [])
                    + (sum([list(b) for b in bl], []) if bl is not None else [])
                    + (ar if ar is not None else [])
                    for a, bl, ar in zip(atoms, bonds, arrows)
                ],
                **animation_options
            )

        if return_objects:
            return figure, atoms, bonds, arrows
        else:
            return figure


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
    subthemes = {
            'default': {
                'multiple_bond_spacing': .12,
                # 'bond_lengths': 1,
                # 'layout_function':'default',
                'disk_options': {'stroke-width': '.01px', 'stroke': 'black'},
                'line_options': {'stroke-width': '.06px', 'stroke': 'black'},
                'label_style': {'font_size': 6}
            }
        }

    # ------------------------------------------------------------------ #
    #  Embedding / poses
    # ------------------------------------------------------------------ #
    def _embed_2d(self, coords, masses, pose=None, principal_axis_order=(0, 1)):
        coords = np.asanyarray(coords, dtype=float)
        i, j = principal_axis_order

        def _pa_project(c):
            com = nput.center_of_mass(c, masses)
            _, axes = nput.moments_of_inertia(c, masses)  # ascending moment; cols=axes
            proj = (c - com[np.newaxis, :]) @ axes
            return np.stack([proj[:, i], proj[:, j]], axis=-1)

        if pose is None:
            return _pa_project(coords)
        if callable(pose):
            return np.asanyarray(pose(coords), dtype=float)

        arr = np.asanyarray(pose, dtype=float)
        n = len(coords)
        if arr.ndim == 2 and arr.shape == (n, 2):
            return arr
        if arr.ndim == 2 and arr.shape == (2, 3):
            com = nput.center_of_mass(coords, masses)
            return (coords - com[np.newaxis, :]) @ arr.T
        if arr.ndim == 2 and arr.shape == (3, 3):
            com = nput.center_of_mass(coords, masses)
            proj = (coords - com[np.newaxis, :]) @ arr
            return np.stack([proj[:, i], proj[:, j]], axis=-1)
        if arr.ndim == 2 and arr.shape == (n, 3):
            return _pa_project(arr)
        raise ValueError(
            f"can't interpret pose of shape {arr.shape} for {n} atoms; "
            "expected (N,2), (N,3), (2,3), (3,3), or a callable"
        )

    def _apply_bond_lengths(self, xy, bond_list, bond_lengths, root=0):
        if bond_lengths is None:
            return xy

        # TODO: write implementation for repulsions

        from collections import deque
        xy = np.array(xy, dtype=float)
        n = len(xy)

        def target_for(a, b):
            if nput.is_numeric(bond_lengths):
                return float(bond_lengths)
            if dev.is_dict_like(bond_lengths):
                key = tuple(sorted((a, b)))
                if key in bond_lengths:
                    return bond_lengths[key]
                ea = self._ats[a]["ElementSymbol"]
                eb = self._ats[b]["ElementSymbol"]
                for k in ((ea, eb), (eb, ea)):
                    if k in bond_lengths:
                        return bond_lengths[k]
                return bond_lengths.get('default', None)
            raise ValueError(f"can't interpret bond_lengths={bond_lengths!r}")

        adj = {k: [] for k in range(n)}
        for b in bond_list:
            a, c = int(b[0]), int(b[1])
            adj[a].append(c)
            adj[c].append(a)

        seen = set()
        for seed in [root] + list(range(n)):
            if seed in seen:
                continue
            seen.add(seed)
            q = deque([seed])
            while q:
                p = q.popleft()
                for c in adj[p]:
                    if c in seen:
                        continue
                    seen.add(c)
                    L = target_for(p, c)
                    if L is not None:
                        vv = xy[c] - xy[p]
                        nv = np.linalg.norm(vv)
                        if nv > 1e-8:
                            xy[c] = xy[p] + (vv / nv) * L
                    q.append(c)
        return xy

    def _handle_overlaps_repulsions(self, xy,
                                    offset_radius=.5,
                                    max_adjustment=1,
                                    max_iterations=10
                                    ):
        cur = xy
        xy = xy.copy()
        # r,c = np.triu_indices(len(xy), k=1)
        for _ in range(max_iterations):
            dm, diffs = nput.distance_matrix(xy, return_diffs=True)
            new_dms = np.linalg.norm(cur - xy, axis=-1)
            new_test = new_dms < max_adjustment
            filt_dm = new_test[:, np.newaxis] & new_test[np.newaxis, :]
            bad_pos = np.where((dm < offset_radius) & filt_dm)
            if len(bad_pos[0]) == 0:
                return xy
            m = bad_pos[1] > bad_pos[0]
            for i, j in zip(bad_pos[0][m], bad_pos[1][m]):
                v, d = nput.vec_normalize(xy[i] - xy[j], return_norms=True)
                di = min([offset_radius / 2, max_adjustment - new_dms[i]])
                dj = min([offset_radius / 2, max_adjustment - new_dms[j]])
                xy[i] += di * v
                xy[j] -= dj * v

        return xy

    # ------------------------------------------------------------------ #
    #  2D primitive builders (mirror the 3D style handling)
    # ------------------------------------------------------------------ #
    @staticmethod
    def _svg2d_clean_style(sty):
        # 2D Disk/Line don't take `glow`; fold it into a stroke outline instead
        sty = dict(sty)
        glow = sty.pop('glow', None)
        if glow is not None and 'stroke' not in sty:
            sty['stroke'] = glow
        return sty

    def _plot_range_2d(self, xy, radii, plot_range_padding):
        lo = np.min(xy, axis=0)
        hi = np.max(xy, axis=0)
        if plot_range_padding is None:
            pad = 0.0
        elif dev.str_is(plot_range_padding, 'auto'):
            pad = float(np.max(radii)) * 1.5
        else:
            pad = float(plot_range_padding)
        return [[lo[0] - pad, hi[0] + pad], [lo[1] - pad, hi[1] + pad]]

    def _atom_draw_flags(self, atom_filter, atom_style):
        # element-symbol visibility for a skeletal-style depiction: heteroatoms
        # always labeled, carbons/hydrogens only when explicitly requested
        # hydrogens = {"H", "D", "T"}
        flags = []
        for i,el in enumerate(self.atoms):
            flag = atom_filter(el, i, plotter=self, **atom_style[i])
            flags.append(flag)
        return flags

    def _get_atom_label_primitives_2d(self, xy, colors, atom_style, atom_text, labeled, *,
                                      text_class, label_style, plot_range, plotos):
        prims = []
        for j, (coord, color, lab) in enumerate(zip(xy, colors, labeled)):
            if not lab:
                continue
            # per-atom user text overrides the element symbol
            override = atom_text[j]
            extra = {}
            if override is not None:
                if isinstance(override, dict):
                    extra = dict(override)
                    text = extra.pop('text')
                else:
                    text = override
            else:
                text = self._ats[j]["ElementSymbol"]
            if nput.is_numeric(text):
                text = str(text)

            a_sty = dict(atom_style.get(j, {}))
            modifier = a_sty.pop('modifier', None)
            col = a_sty.pop('color', None) or color
            if modifier is not None:
                a_sty = modifier(j, self.atoms[j], dict(a_sty, color=col))
                col = a_sty.pop('color', col)
            a_sty = self._svg2d_clean_style(a_sty)

            # atom color wins over the base label color so labels read as CPK-ish;
            # a per-atom atom_style color still overrides via `col`
            sty = (label_style | {
                'color': col,
                'use_path':True,
                'invert':True,
                'anchor':(-.5, 1.25),
                'plot_range': plot_range
            } | a_sty | extra)
            pos = coord + np.asanyarray(sty.pop('offset', [0.0, 0.0]), dtype=float)
            prims.append(text_class(text, pos, **sty))
        return prims

    def _get_bond_primitives_2d(self, xy, bond_list, radii, colors, glows, labeled, *,
                                line_class, bond_style, bond_radius, line_options,
                                trim_bonds, render_multiple_bonds, multiple_bond_spacing,
                                half_colored_bonds, theme_function, plotos):
        prims = []
        base = dict(line_options)
        if bond_radius is not None and 'line_width' not in base and 'stroke-width' not in base:
            base['line_width'] = bond_radius
        for b in bond_list:
            i, j, order = b[0], b[1], b[2]
            pi, pj = xy[i], xy[j]
            d = pj - pi
            L = np.linalg.norm(d)
            if L < 1e-8:
                continue
            u = d / L
            perp = np.array([-u[1], u[0]])
            # trim back to leave room for a label; extend to the vertex center when
            # an endpoint atom is unlabeled (e.g. an omitted carbon)
            start = pi + u * radii[i] if (trim_bonds and labeled[i]) else pi
            end = pj - u * radii[j] if (trim_bonds and labeled[j]) else pj

            base_bstyle = dict(bond_style.get((j, i), {}), **bond_style.get((i, j), {}))
            b1 = dict(bond_style.get(i, {}), **base_bstyle)
            mod = b1.pop('modifier', None)
            if mod is not None:
                b1 = mod((i, j), (self.atoms[i], self.atoms[j]), b1)
            b1 = self._svg2d_clean_style(b1)
            b2 = dict(bond_style.get(j, {}), **base_bstyle)
            mod = b2.pop('modifier', None)
            if mod is not None:
                b2 = mod((j, i), (self.atoms[j], self.atoms[i]), b2)
            b2 = self._svg2d_clean_style(b2)

            n_lines = order if (render_multiple_bonds and order and order >= 1) else 1
            fraction = order - int(order)
            offsets = (np.arange(n_lines) - (n_lines - 1) / 2) * multiple_bond_spacing

            for fi,off in enumerate(offsets):
                s = start + perp * off
                e = end + perp * off
                if fi%2 == 1:
                    v = (e - s) * fraction/2
                    e = e - v
                    s = s + v
                if half_colored_bonds and b1.get('color') is None and b2.get('color') is None:
                    mid = (s + e) / 2
                    for pts, c, over in (([s, mid], colors[i], b1), ([mid, e], colors[j], b2)):
                        sty = (plotos | base | {'color': c} | over)
                        if theme_function is not None:
                            sty = theme_function((i, j), line_class, sty)
                        prims.append(line_class(np.array(pts), **sty))
                else:
                    sty = (plotos | base | b1)
                    if 'color' not in sty:
                        sty['color'] = colors[i]
                    if theme_function is not None:
                        sty = theme_function((i, j), line_class, sty)
                    prims.append(line_class(np.array([s, e]), **sty))
        return prims


    @classmethod
    def _default_atom_filter(cls, draw_carbons, draw_hydrogens):
        def filter(atom_type, i, *, plotter, **opts):
            if atom_type == 'C':
                if draw_carbons:
                    return True
                else:
                    return 'line'
            elif atom_type == 'H':
                if draw_hydrogens:
                    return True
                else:
                    return False
            else:
                return True
        return filter


    # ------------------------------------------------------------------ #
    #  Entry point
    # ------------------------------------------------------------------ #
    def plot_impl(self, full_opts):
        from McUtils.Plots import Graphics, Disk, Line, Text

        figure = full_opts.pop('figure', None)

        units = full_opts.pop('units', 'Angstroms')
        full_opts.pop('reconcile_bonds', None)
        plot_ops = full_opts.pop('extra_opts')

        # split 2D-Graphics figure options out of the free kwargs
        gk = (getattr(Graphics, 'known_keys', set()) | getattr(Graphics, 'opt_keys', set())
              | getattr(Graphics, 'figure_keys', set()) | getattr(Graphics, 'axes_keys', set()))
        graphics_opts = {k: plot_ops[k] for k in plot_ops.keys() & gk}
        plot_ops = {k: plot_ops[k] for k in plot_ops.keys() - gk}

        def _pop(k, default=None):
            return plot_ops.pop(k, full_opts.pop(k, default))

        # 2D-specific options (arrive via **plot_ops -> extra_opts)
        pose = _pop('pose', None)
        principal_axis_order = _pop('principal_axis_order', (0, 1))
        bond_lengths = _pop('bond_lengths', None)
        bond_layout_root = _pop('bond_layout_root', 0)
        trim_bonds = _pop('trim_bonds', True)
        half_colored_bonds = _pop('half_colored_bonds', True)
        draw_carbons = _pop('draw_carbons', False)
        draw_hydrogens = _pop('draw_hydrogens', False)
        atom_filter = _pop('atom_filter', None)
        text_class = _pop('text_class', None) or Text
        layout_function = _pop('layout_function')
        image_size = graphics_opts.pop('image_size', _pop('image_size', (400, 400)))
        offset_radius = _pop('offset_radius', .2)
        max_adjustment = _pop('max_adjustment', 1)
        max_adjustment_iterations = _pop('max_adjustment_iterations', 10)

        v = full_opts
        return_objects = v.pop('return_objects')
        objects = v.pop('objects')
        atom_style = v.pop('atom_style')
        atom_radii = v.pop('atom_radii')
        atom_radius_scaling = v.pop('atom_radius_scaling')
        radius_type = v.pop('radius_type')
        atom_text = v.pop('atom_text')
        display_atom_numbers = v.pop('display_atom_numbers')
        label_style = v.pop('label_style') or {}
        bond_style = v.pop('bond_style')
        bond_radius = v.pop('bond_radius')
        capped_bonds = v.pop('capped_bonds')
        reflectiveness = v.pop('reflectiveness')
        highlight_atoms = v.pop('highlight_atoms')
        highlight_bonds = v.pop('highlight_bonds')
        highlight_rings = v.pop('highlight_rings')
        highlight_styles = v.pop('highlight_styles')
        bonds = v.pop('bonds')
        render_multiple_bonds = v.pop('render_multiple_bonds')
        multiple_bond_spacing = v.pop('multiple_bond_spacing') or .12
        disk_class = v.pop('disk_class') or Disk
        disk_options = v.pop('disk_options') or {}
        line_class = v.pop('line_class') or Line
        line_options = v.pop('line_options') or {}
        theme_function = v.pop('theme_function')
        plot_range_padding = v.pop('plot_range_padding')
        annotation_function = v.pop('annotation_function')
        graphics_class = v.pop('graphics_class') or Graphics
        background = graphics_opts.pop('background', 'transparent')
        # (any remaining full_opts keys are 3D-only and ignored on this path)

        if atom_filter is None:
            atom_filter = self._default_atom_filter(draw_carbons=draw_carbons, draw_hydrogens=draw_hydrogens)

        masses = self.atomic_masses
        conv = UnitsData.convert("BohrRadius", units) if units is not None else 1.0

        if highlight_rings is not None:
            highlight_atoms = list(highlight_atoms or [])
            highlight_bonds = list(highlight_bonds or [])
            for r in highlight_rings:
                highlight_atoms.extend(r)
                highlight_bonds.extend(zip(r, r[1:] + r[:1]))

        atom_style, highlight_atoms, colors, glows = self._prep_display_atom_style(
            atom_style, highlight_atoms, backend='svg',
            reflectiveness=reflectiveness, highlight_styles=highlight_styles)
        bond_style, highlight_bonds = self._prep_display_bond_style(
            bond_style, highlight_bonds, backend='svg', reflectiveness=reflectiveness,
            highlight_atoms=highlight_atoms, highlight_styles=highlight_styles,
            capped_bonds=capped_bonds)
        atom_text = self._prep_display_atom_text(atom_text, display_atom_numbers, label_style)
        radii = np.asanyarray(self._get_atom_radii(atom_radii, atom_radius_scaling, radius_type))

        # which atoms get an element-symbol label; an explicit atom_text entry
        # forces a label even on an otherwise-omitted carbon/hydrogen
        draw_flags = self._atom_draw_flags(atom_filter, atom_style)
        labeled = [
            df or (atom_text[j] is not None)
            for j, df in enumerate(draw_flags)
        ]

        if bonds is None:
            draw_bonds = self.bonds
        elif bonds is False:
            draw_bonds = []
        else:
            draw_bonds = bonds
        draw_bonds = [
            b for b in draw_bonds
            if labeled[b[0]] and labeled[b[1]]
        ]

        # normalize geometries to a stack of frames, then embed each
        geom_arr = np.asanyarray(self.geometries)
        if geom_arr.ndim == 2:
            geom_arr = geom_arr[np.newaxis]
        geoms = [geom_arr[k] * conv for k in range(geom_arr.shape[0])]

        embeddings = []
        for g in geoms:
            if layout_function is not None:
                if dev.is_dict_like(layout_function):
                    layout_opts = layout_function.copy()
                    layout_function = layout_function.pop('layout')
                else:
                    layout_opts = {}
                if isinstance(layout_function, str):
                    coords = self.mol.edge_graph.layout(layout_function, **layout_opts)
                    c = np.array(list(coords.values()))
                    i = np.argsort([c[1] for c in coords])
                    xy = c[i] * 20 #TODO: cleanup scaling functions
                else:
                    xy = layout_function(self.mol, **layout_opts)
            else:
                xy = self._embed_2d(g, masses, pose=pose, principal_axis_order=principal_axis_order)
                xy = self._apply_bond_lengths(xy, draw_bonds, bond_lengths, root=bond_layout_root)
            xy = self._handle_overlaps_repulsions(xy,
                                                  offset_radius=offset_radius,
                                                  max_adjustment=max_adjustment,
                                                  max_iterations=max_adjustment_iterations)
            embeddings.append(xy)

        if figure is None:
            all_xy = np.concatenate(embeddings, axis=0)
            fig_opts = dict(
                backend='svg', image_size=image_size, aspect_ratio='equal',
                frame=False, background=background,
                plot_range=self._plot_range_2d(all_xy, radii, plot_range_padding),
            )
            fig_opts.update(graphics_opts)
            figure = graphics_class(**fig_opts)

        plot_range = figure.plot_range

        def _render(prims):
            if objects:
                return prims
            out = []
            for p in prims:
                art = p.plot(figure)
                if isinstance(art, (list, tuple)):
                    art = art[0]
                out.append(art)
            return out

        atoms_out = [None] * len(geoms)
        bonds_out = [None] * len(geoms)
        labels_out = [None] * len(geoms)

        labeled = [l is True for l in labeled]
        for i, xy in enumerate(embeddings):
            plotos = plot_ops.copy()
            bond_prims = self._get_bond_primitives_2d(
                xy, draw_bonds, radii, colors, glows, labeled,
                line_class=line_class, bond_style=bond_style, bond_radius=bond_radius,
                line_options=line_options, trim_bonds=trim_bonds,
                render_multiple_bonds=render_multiple_bonds,
                multiple_bond_spacing=multiple_bond_spacing,
                half_colored_bonds=half_colored_bonds,
                theme_function=theme_function, plotos=plotos)
            # atoms are drawn as element-symbol text (skeletal style); carbons and
            # hydrogens are only labeled when draw_carbons / draw_hydrogens is set
            atom_prims = self._get_atom_label_primitives_2d(
                xy, colors, atom_style, atom_text, labeled,
                text_class=text_class, label_style=label_style,
                plot_range=plot_range,
                plotos=plotos)
            extra_prims = []
            if annotation_function is not None:
                extra_prims = list(annotation_function(self.mol, i, xy))

            bonds_out[i] = _render(bond_prims)
            atoms_out[i] = _render(atom_prims)
            labels_out[i] = _render(extra_prims)

        if objects:
            packed = [
                {'atoms': a, 'bonds': b, 'labels': l}
                for a, b, l in zip(atoms_out, bonds_out, labels_out)
            ]
            return packed[0] if len(geoms) == 1 else packed
        if return_objects:
            return figure, atoms_out, bonds_out, labels_out
        return figure