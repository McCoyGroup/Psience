"""Periodic cells and asymmetric units backed by ordinary Molecule objects."""

import itertools
import copy
import os
import re

import numpy as np

from McUtils.Data import UnitsData
from McUtils.ExternalPrograms import CIFParser, CIFConverter
from McUtils.ExternalPrograms.Readouts.Crystals import expand_cif_sites
from McUtils.Jupyter.Readouts import cell_matrix, element_symbol
from McUtils.Plots import Line, Cylinder, Arrow, Text
from McUtils.Parsers import StringLineByLineReader

from .Molecule import Molecule

__all__ = ["PeriodicSystem"]


class _PeriodicCIFParser(CIFParser, StringLineByLineReader):
    """Extend the McUtils reader for uncertain numbers and modern symmetry tags."""

    def __init__(self, file, **options):
        # CIFParser's file stream uses a writable mmap; CIF import is read-only.
        if hasattr(file, "read"):
            text = file.read()
        else:
            with open(os.fspath(file), encoding="utf-8") as stream:
                text = stream.read()
        if isinstance(text, bytes):
            text = text.decode("utf-8")
        StringLineByLineReader.__init__(self, text, max_nesting_depth=1, **options)
        self.fields = None

    def check_tag(self, line, depth=0, active_tag=None, label=None, history=None):
        if line[:5].lower() == "data_":
            return self.LineReaderTags.BLOCK_START if active_tag is None else self.LineReaderTags.RESETTING_BLOCK_END
        if line[:5].lower() == "loop_":
            line = "loop_" + line[5:]
        return super().check_tag(line, depth, active_tag=active_tag, label=label, history=history)

    def _get_float(self, key, value):
        def number(v):
            v = re.sub(r"\(\d+\)", "", str(v).strip())
            return np.nan if v in (".", "?") else float(v)
        return number(value) if isinstance(value, str) else np.array([number(v) for v in value])

    def get_block_handlers(self):
        # Leave publication/other metadata as text; only structural fields need typing.
        handlers = {k: v for k, v in super().get_block_handlers().items()
                    if k.startswith(("cell_", "atom_site_fract_", "atom_site_occupancy"))}
        symmetry_handler = super().get_block_handlers()["symmetry_equiv_pos_as_xyz"]

        def symmetry(key, value):
            values = [value] if isinstance(value, str) else value
            return symmetry_handler(key, [re.sub(r"\s+", "", v.lower()) for v in values])
        handlers.update(symmetry_equiv_pos_as_xyz=symmetry, space_group_symop_operation_xyz=symmetry)
        return handlers

    def handle_block(self, label, block_data, join=True, depth=0):
        if label is not None:
            label = label.lower().replace(".", "_")
        elif block_data and all(isinstance(v, str) for v in block_data):
            # CIF loops are token streams: a logical row may span physical lines.
            keys, values = [], []
            for line in block_data:
                if not values and line.startswith("_"):
                    keys.append(line[1:].strip())
                else:
                    values.extend(self._tokenize(line))
            if keys:
                if len(values) % len(keys):
                    raise ValueError("CIF loop has an incomplete row")
                return {k: self.handle_block(k, values[i::len(keys)], join=False)
                        for i, k in enumerate(keys)}
        return super().handle_block(label, block_data, join=join, depth=depth)


class PeriodicSystem:
    """A molecule plus a periodic lattice, without inheriting from Molecule.

    ``lattice`` contains the three Cartesian lattice vectors as **rows**. Its
    lengths and ``origin`` use ``lattice_units`` (Bohr by default), while the
    stored molecule follows Molecule's usual Bohr convention. Internally all
    lengths are Bohr. Fractional coordinates obey ``xyz = origin + frac @ lattice``.

    Set ``symmetry_reduced=True`` when the molecule contains asymmetric-unit
    representatives instead of a realized cell. ``symmetry_operations`` is the
    complete list of fractional-coordinate 4x4 affine matrices, acting on column
    vectors. Identity is included automatically. Tolerances are fractional.
    ``periodicity`` enables translation/wrapping along a, b, and c independently;
    retain independent lattice rows for finite directions (e.g. slab vacuum).
    """

    def __init__(self, molecule, lattice, *, symmetry_reduced=False,
                 symmetry_operations=None, origin=None, lattice_units="BohrRadius",
                 space_group=None, site_labels=None, occupancies=None, metadata=None,
                 periodicity=(True, True, True)):
        if not isinstance(molecule, Molecule):
            raise TypeError("molecule must be a Molecule")
        if np.shape(molecule.coords) != (len(molecule.atoms), 3) or not np.all(np.isfinite(molecule.coords)):
            raise ValueError("a periodic system requires one finite (natoms, 3) geometry")
        scale = 1 if lattice_units is None else UnitsData.convert(lattice_units, "BohrRadius")
        lattice = np.array(lattice, dtype=float, copy=True) * scale
        if (lattice.shape != (3, 3) or not np.all(np.isfinite(lattice))
                or np.linalg.matrix_rank(lattice) != 3):
            raise ValueError("lattice must contain three finite, independent row vectors")
        origin = np.zeros(3) if origin is None else np.array(origin, dtype=float, copy=True) * scale
        if origin.shape != (3,) or not np.all(np.isfinite(origin)):
            raise ValueError("origin must be a finite three-vector")
        if isinstance(periodicity, (int, np.integer)) and not isinstance(periodicity, (bool, np.bool_)):
            if not 1 <= periodicity <= 3:
                raise ValueError("periodicity dimension must be 1, 2, or 3")
            periodicity = np.arange(3) < periodicity
        periodicity = np.asarray(periodicity)
        if periodicity.shape != (3,) or periodicity.dtype != np.dtype(bool):
            raise ValueError("periodicity must contain three booleans for a, b, and c")
        ops = np.eye(4)[None] if symmetry_operations is None else np.array(symmetry_operations, dtype=float)
        if ops.shape == (4, 4):
            ops = ops[None]
        if (ops.ndim != 3 or ops.shape[1:] != (4, 4) or len(ops) == 0
                or not np.all(np.isfinite(ops)) or not np.allclose(ops[:, 3], [0, 0, 0, 1])):
            raise ValueError("symmetry_operations must be finite affine matrices of shape (n, 4, 4)")
        rotations = ops[:, :3, :3]
        metric = lattice @ lattice.T
        if (not np.allclose(rotations, np.round(rotations), atol=1e-8, rtol=0)
                or not np.allclose(np.abs(np.linalg.det(rotations)), 1)
                or not np.allclose(rotations.transpose(0, 2, 1) @ metric @ rotations, metric)):
            raise ValueError("symmetry rotations must preserve the supplied lattice metric")
        if not np.allclose(rotations[:, ~periodicity][:, :, periodicity], 0, atol=1e-8, rtol=0):
            raise ValueError("symmetry operations must preserve the periodic directions")
        ops = np.concatenate([np.eye(4)[None], ops])
        ops[:, np.flatnonzero(periodicity), 3] %= 1
        _, unique = np.unique(np.round(ops.reshape(len(ops), -1), 8), axis=0, return_index=True)
        n = len(molecule.atoms)
        labels = list(site_labels) if site_labels is not None else [f"{a}{i + 1}" for i, a in enumerate(molecule.atoms)]
        occupancies = np.ones(n) if occupancies is None else np.array(occupancies, dtype=float, copy=True)
        if len(labels) != n or occupancies.shape != (n,) or not np.all(np.isfinite(occupancies)):
            raise ValueError("site labels and finite occupancies must have one entry per atom")
        if np.any((occupancies < 0) | (occupancies > 1)):
            raise ValueError("occupancies must lie between zero and one")
        self.molecule = molecule
        self.lattice = lattice
        self.origin = origin
        self.periodicity = periodicity.copy()
        self.symmetry_reduced = bool(symmetry_reduced)
        self.symmetry_operations = ops[np.sort(unique)].copy()
        self.space_group = dict(space_group or {})
        self.site_labels = labels
        self.occupancies = occupancies
        self.metadata = dict(metadata or {})

    @property
    def fractional_coordinates(self):
        return np.linalg.solve(self.lattice.T, (np.asarray(self.molecule.coords) - self.origin).T).T

    @property
    def cell_vertices(self):
        return self.origin + self._cell_fractions() @ self.lattice

    def _cell_fractions(self):
        return np.array(list(itertools.product(*((0, 1) if p else (0,) for p in self.periodicity))))

    @staticmethod
    def _tolerance(tol):
        if not np.isfinite(tol) or not 0 < tol < .5:
            raise ValueError("tol must be a positive fractional tolerance smaller than 0.5")
        return tol

    def _wrap(self, frac, tol):
        frac = np.array(frac, dtype=float, copy=True).reshape(-1, 3)
        periodic = frac[:, self.periodicity] % 1
        frac[:, self.periodicity] = np.where(np.abs(periodic - 1) < tol, 0, periodic)
        return frac

    def _matching_sites(self, point, source, fractional, indices, tol):
        indices = np.asarray(indices, dtype=int)
        delta = fractional - point
        delta[:, self.periodicity] -= np.round(delta[:, self.periodicity])
        close = np.linalg.norm(delta, axis=-1) <= tol
        same = np.array([self.molecule.atoms[i] == self.molecule.atoms[source] for i in indices], dtype=bool)
        same &= np.isclose(self.molecule.masses[indices], self.molecule.masses[source], atol=1e-8, rtol=0)
        same &= np.isclose(self.occupancies[indices], self.occupancies[source], atol=1e-8, rtol=0)
        return np.flatnonzero(close & same)

    def _new_sites(self, fractional, indices, *, lattice=None, operations=None,
                   symmetry_reduced=False, mappings=()):
        """Rebuild geometry/topology; molecular energy/derivative caches do not transfer."""
        indices = np.asarray(indices, dtype=int)
        lattice = self.lattice if lattice is None else lattice
        # Read stored bonds without triggering molecular bond inference for a CIF.
        bonds = self.molecule._bonds
        if bonds is not None:
            mapped = {}
            for mapping in mappings:
                for i, j, *order in bonds:
                    if i in mapping and j in mapping and mapping[i] != mapping[j]:
                        a, b = sorted((mapping[i], mapping[j]))
                        mapped[a, b] = [a, b, *order]
            bonds = list(mapped.values())
        molecule = Molecule([self.molecule.atoms[i] for i in indices],
                            self.origin + np.asarray(fractional).reshape(-1, 3) @ lattice,
                            masses=self.molecule.masses[indices], bonds=bonds,
                            guess_bonds=self.molecule.guess_bonds, name=self.molecule.name,
                            source_file=self.molecule.source_file,
                            formal_charges=(None if self.molecule.formal_charges is None else
                                            np.asarray(self.molecule.formal_charges)[indices]),
                            display_mode=self.molecule.display_mode,
                            display_settings=self.molecule.display_settings)
        return type(self)(molecule, lattice, origin=self.origin,
                          symmetry_reduced=symmetry_reduced,
                          symmetry_operations=self.symmetry_operations if operations is None else operations,
                          space_group=self.space_group, site_labels=[self.site_labels[i] for i in indices],
                          occupancies=self.occupancies[indices], metadata=self.metadata, periodicity=self.periodicity)

    def copy(self):
        molecule = self.molecule.copy()
        # Molecule.copy shares coordinate storage; detach it for an independent cell.
        molecule.coords = np.array(self.molecule.coords, copy=True)
        return type(self)(molecule, self.lattice, origin=self.origin,
                          symmetry_reduced=self.symmetry_reduced, symmetry_operations=self.symmetry_operations,
                          space_group=self.space_group, site_labels=self.site_labels,
                          occupancies=self.occupancies, metadata=self.metadata, periodicity=self.periodicity)

    def expand_symmetry(self, tol=1e-6):
        """Realize a reduced cell, wrapping sites and removing equivalent duplicates."""
        tol = self._tolerance(tol)
        if not self.symmetry_reduced:
            return self.copy()
        if np.all(self.periodicity):
            _, fractional, sources = expand_cif_sites(self.molecule.atoms, self.fractional_coordinates,
                                                     self.symmetry_operations, tol=tol)
        else:
            images = np.einsum("oij,sj->soi", self.symmetry_operations[:, :3, :3], self.fractional_coordinates)
            fractional = (images + self.symmetry_operations[None, :, :3, 3]).reshape(-1, 3)
            sources = np.repeat(np.arange(len(self.molecule.atoms)), len(self.symmetry_operations))
        fractional = self._wrap(fractional, tol).reshape(-1, 3)
        kept = []
        for k, (point, source) in enumerate(zip(fractional, sources)):
            if not len(self._matching_sites(point, source, fractional[kept], sources[kept], tol)):
                kept.append(k)
        fractional, sources = fractional[kept], sources[kept]
        mappings = []
        for op in self.symmetry_operations:
            images = self.fractional_coordinates @ op[:3, :3].T + op[:3, 3]
            mappings.append({i: int(self._matching_sites(p, i, fractional, sources, tol)[0])
                             for i, p in enumerate(images)})
        return self._new_sites(fractional, sources, mappings=mappings)

    def reduce_symmetry(self, tol=1e-6):
        """Keep one representative per supplied symmetry orbit, preserving species/occupancy.

        The molecule must contain a complete realized cell for these operations.
        This reduces the atom list, not the lattice to a primitive cell.
        """
        tol = self._tolerance(tol)
        if self.symmetry_reduced:
            return self.copy()
        fractional = self.fractional_coordinates
        indices = np.arange(len(fractional))
        covered, kept = set(), []
        for i, point in enumerate(fractional):
            images = self.symmetry_operations[:, :3, :3] @ point + self.symmetry_operations[:, :3, 3]
            orbit = []
            for image in images:
                matches = self._matching_sites(image, i, fractional, indices, tol)
                if not len(matches):
                    raise ValueError("supplied symmetry operations do not preserve the realized cell")
                orbit.extend(matches.tolist())
            if i not in covered:
                kept.append(i)
                covered.update(orbit)
        return self._new_sites(self._wrap(fractional[kept], tol), kept, symmetry_reduced=True,
                               mappings=[{old: new for new, old in enumerate(kept)}])

    def tile(self, repetitions, tol=1e-6):
        """Return a realized n_a x n_b x n_c supercell (or n x n x n for an integer).

        Retain operations compatible with the supercell lattice and add the
        internal cell translations. Supplied bonds are replicated, and Molecule's
        normal bond inference connects copies when guess_bonds is enabled.
        """
        tol = self._tolerance(tol)
        repetitions = np.asarray(repetitions)
        if not np.issubdtype(repetitions.dtype, np.integer):
            raise ValueError("repetitions must be a positive integer or three positive integers")
        if repetitions.ndim == 0:
            repetitions = np.where(self.periodicity, repetitions, 1)
        if (repetitions.shape != (3,) or not np.issubdtype(repetitions.dtype, np.integer)
                or np.any(repetitions <= 0)):
            raise ValueError("repetitions must be a positive integer or three positive integers")
        if np.any(repetitions[~self.periodicity] != 1):
            raise ValueError("nonperiodic directions must have one repetition")
        system = self.expand_symmetry(tol=tol)
        translations = np.array(list(itertools.product(*(range(int(n)) for n in repetitions))))
        frac = self._wrap(system.fractional_coordinates, tol)
        tiled = ((frac[None] + translations[:, None]) / repetitions).reshape(-1, 3)
        n = len(frac)
        sources = np.tile(np.arange(n), len(translations))
        operations = []
        for op in system.symmetry_operations:
            rotation = op[:3, :3] * repetitions[None, :] / repetitions[:, None]
            if np.allclose(rotation, np.round(rotation), atol=1e-8, rtol=0):
                for translation in translations:
                    transformed = np.eye(4)
                    transformed[:3, :3] = np.round(rotation)
                    transformed[:3, 3] = (op[:3, 3] + translation) / repetitions
                    operations.append(transformed)
        mappings = [{i: k * n + i for i in range(n)} for k in range(len(translations))]
        result = system._new_sites(tiled, sources, lattice=system.lattice * repetitions[:, None],
                                    operations=operations, mappings=mappings)
        system._connect_cells(result.molecule, np.repeat(np.arange(len(translations)), n))
        return result

    @classmethod
    def from_cif(cls, file, *, block=0, symmetry_reduced=False, symmetry_operations=None,
                 periodicity=(True, True, True), **molecule_options):
        """Read one CIF block using McUtils (cell lengths in Angstroms, angles in degrees).

        Preserve labels, occupancies, space-group fields, and the source CIF block.
        Non-P1 files without explicit operations require ``symmetry_operations``;
        space-group names alone are not used to guess a setting or origin.
        """
        with _PeriodicCIFParser(file, ignore_comments=True) as parser:
            blocks = parser.parse()
        if not isinstance(block, (int, np.integer)) or not 0 <= block < len(blocks):
            raise ValueError("block must select an existing CIF data block")
        data = {k.lower().replace(".", "_"): v for k, v in blocks[block].items()}
        converter = CIFConverter([data])
        cell, sites = converter.cell_properties, converter.atom_properties
        try:
            parameters = [cell[f"cell_length_{a}"] for a in "abc"] + [cell[f"cell_angle_{a}"] for a in ("alpha", "beta", "gamma")]
            if (not np.all(np.isfinite(parameters)) or np.any(np.asarray(parameters[:3]) <= 0)
                    or np.any((np.asarray(parameters[3:]) <= 0) | (np.asarray(parameters[3:]) >= 180))):
                raise ValueError("invalid CIF cell lengths or angles")
            lattice = cell_matrix(*parameters)
            # cell_matrix clamps imaginary heights; reject impossible cell angles.
            cosines = np.cos(np.deg2rad(parameters[3:]))
            if 1 + 2 * np.prod(cosines) - np.sum(cosines ** 2) <= 0:
                raise ValueError("CIF cell angles do not define a nondegenerate lattice")
            fractional = np.column_stack([np.atleast_1d(sites[f"atom_site_fract_{a}"]) for a in "xyz"])
        except KeyError as exc:
            raise ValueError(f"missing required CIF cell/site field: {exc.args[0]}") from exc
        labels = sites.get("atom_site_label")
        if labels is not None:
            labels = np.atleast_1d(labels).tolist()
        symbols = sites.get("atom_site_type_symbol", labels)
        if symbols is None:
            raise ValueError("CIF needs atom_site_type_symbol or atom_site_label")
        atoms = [element_symbol(s) for s in np.atleast_1d(symbols)]
        space_group = {k: v for k, v in data.items() if k.startswith(("space_group_", "symmetry_"))
                       and k not in ("symmetry_equiv_pos_as_xyz", "space_group_symop_operation_xyz")}
        if symmetry_operations is None:
            ops = data.get("space_group_symop_operation_xyz", data.get("symmetry_equiv_pos_as_xyz"))
            if ops is not None:
                symmetry_operations = ops.transformation
            else:
                number = data.get("space_group_it_number", data.get("symmetry_int_tables_number"))
                name = data.get("space_group_name_h-m_alt", data.get("symmetry_space_group_name_h-m"))
                if ((number not in (None, ".", "?") and float(number) != 1)
                        or (name not in (None, ".", "?") and re.sub(r"\s+", "", name).upper() != "P1")):
                    raise ValueError("non-P1 CIF has no explicit symmetry operations; supply symmetry_operations")
        molecule_options.setdefault("source_file", getattr(file, "name", None) if hasattr(file, "read") else os.fspath(file))
        molecule = Molecule(atoms, fractional @ lattice * UnitsData.convert("Angstroms", "BohrRadius"), **molecule_options)
        system = cls(molecule, lattice, lattice_units="Angstroms", symmetry_reduced=True,
                     symmetry_operations=symmetry_operations, space_group=space_group,
                     site_labels=labels, occupancies=np.atleast_1d(sites.get("atom_site_occupancy", np.ones(len(atoms)))),
                     metadata={"cif": data}, periodicity=periodicity)
        return system if symmetry_reduced else system.expand_symmetry()

    @staticmethod
    def _direction_styles(styles, defaults):
        """Common styles plus a/b/c overrides, or a sequence of three style dicts."""
        if styles is None:
            styles = {}
        if isinstance(styles, dict):
            common = {k: v for k, v in styles.items() if k not in ("a", "b", "c")}
            per_axis = [styles.get(axis, {}) for axis in "abc"]
        else:
            common, per_axis = {}, list(styles)
            if len(per_axis) != 3:
                raise ValueError("direction styles must contain three dictionaries")
        if any(not isinstance(style, dict) for style in per_axis):
            raise ValueError("direction styles must be dictionaries")
        return [copy.deepcopy(defaults[i] | common | per_axis[i]) for i in range(3)]

    def _cell_annotations(self, scale, draw_cell, cell_style, unit_axes, unit_axes_style,
                          unit_axes_label_style, unit_axes_origin, unit_axes_scale):
        primitives = []
        if draw_cell:
            styles = self._direction_styles(cell_style, [dict(color="#697586", thickness=1)] * 3)
            fractions = self._cell_fractions()
            vertices = (self.origin + fractions @ self.lattice) * scale
            for a in range(len(vertices)):
                for b in range(a + 1, len(vertices)):
                    changed = np.flatnonzero(fractions[a] != fractions[b])
                    if len(changed) == 1:
                        style = styles[changed[0]].copy()
                        radius = style.pop("radius", None)
                        if radius is None:
                            primitives.append(Line(vertices[[a, b]], **style))
                        else:
                            if not np.isfinite(radius) or radius <= 0:
                                raise ValueError("cell edge radius must be positive")
                            style.pop("thickness", None)
                            primitives.append(Cylinder(vertices[a], vertices[b], radius, **style))
        axes = np.asarray(unit_axes)
        if axes.ndim == 0 and axes.dtype == np.dtype(bool):
            axes = np.repeat(axes, 3)
        if axes.shape != (3,) or axes.dtype != np.dtype(bool):
            raise ValueError("unit_axes must be a boolean or three booleans")
        axes = axes & self.periodicity
        if np.any(axes):
            if not np.isfinite(unit_axes_scale) or unit_axes_scale <= 0:
                raise ValueError("unit_axes_scale must be positive")
            origin = self.origin if unit_axes_origin is None else np.asarray(unit_axes_origin, dtype=float)
            if origin.shape != (3,) or not np.all(np.isfinite(origin)):
                raise ValueError("unit_axes_origin must be a finite three-vector in Bohr")
            defaults = [dict(color=color, radius=.035) for color in ("#d35454", "#38906b", "#477ec3")]
            styles = self._direction_styles(unit_axes_style, defaults)
            label_styles = self._direction_styles(unit_axes_label_style, [dict(font_size=.35)] * 3)
            for i in np.flatnonzero(axes):
                vector = self.lattice[i] * unit_axes_scale
                primitives.append(Arrow(origin * scale, (origin + vector) * scale, **styles[i]))
                label_style = dict(color=styles[i].get("color", "black")) | label_styles[i]
                primitives.append(Text("abc"[i], (origin + 1.08 * vector) * scale, **label_style))
        return primitives

    def _image_translations(self, periodic_images):
        message = ("periodic_images must have three entries: nonnegative integer/boolean counts "
                   "or [negative, positive] count pairs")
        if periodic_images is None:
            counts = [1] * 3
        elif np.isscalar(periodic_images) or (isinstance(periodic_images, np.ndarray) and periodic_images.ndim == 0):
            counts = [periodic_images] * 3
        else:
            counts = list(periodic_images)
        if len(counts) != 3:
            raise ValueError(message)
        ranges = []
        for count, periodic in zip(counts, self.periodicity):
            pair = np.asarray(count)
            if pair.ndim == 0:
                pair = np.repeat(pair, 2)
            if pair.shape != (2,) or pair.dtype.kind not in "biu" or np.any(pair < 0):
                raise ValueError(message)
            negative, positive = (int(n) for n in pair) if periodic else (0, 0)
            ranges.append(range(-negative, positive + 1))
        shifts = [p for p in itertools.product(*ranges) if any(p)]
        return np.array(shifts, dtype=int).reshape(-1, 3) @ self.lattice

    @staticmethod
    def _image_component_style(original, common, overrides):
        if overrides is False or (original is False and overrides is None):
            return False
        if original is None or isinstance(original, (bool, np.bool_)):
            original = {}
        if not isinstance(original, dict):
            original = dict(enumerate(original))
        styled = copy.deepcopy(original)
        for key, value in styled.items():
            if isinstance(value, dict):
                styled[key] = value | common
        styled.update(common)
        if overrides is not None and overrides is not True:
            if not isinstance(overrides, dict):
                overrides = dict(enumerate(overrides))
            for key, value in overrides.items():
                styled[key] = styled.get(key, {}) | copy.deepcopy(value) if isinstance(value, dict) else value
        return styled

    def _connect_cells(self, molecule, cells):
        """Use the standard Molecule inference; retain supplied intra-cell topology."""
        if not self.molecule.guess_bonds or len(np.unique(cells)) < 2:
            return []
        inferred = molecule.get_guessed_bonds(mode=self.molecule.bond_guessing_mode)
        bridges = {}
        for i, j, *order in inferred:
            if cells[i] != cells[j]:
                a, b = sorted((int(i), int(j)))
                bridges[a, b] = [a, b, *order]
        if molecule._bonds is None:
            molecule.bonds = inferred
        else:
            bonds = {tuple(sorted(b[:2])): list(b) for b in molecule._bonds}
            for key, bond in bridges.items():
                bonds.setdefault(key, bond)
            molecule.bonds = list(bonds.values())
        return list(bridges.values())

    def _image_bond_network(self, translations, bonds):
        """Build the same finite neighborhood used by plotting, with native bonds."""
        shifts = np.concatenate([np.zeros((1, 3)), translations])
        n, count = len(self.molecule.atoms), len(shifts)
        coords = (np.asarray(self.molecule.coords)[None] + shifts[:, None]).reshape(-1, 3)
        if bonds is None:
            bonds = self.molecule.bonds
        elif isinstance(bonds, str):
            if bonds != "recompute":
                raise ValueError("bonds must be a bond list, False, or 'recompute'")
            bonds = self.molecule.get_guessed_bonds()
        repeated = None if bonds is None else [
            [int(i) + k * n, int(j) + k * n, *order]
            for k in range(count) for i, j, *order in bonds]
        charges = self.molecule.formal_charges
        molecule = Molecule(self.molecule.atoms * count, coords, bonds=repeated,
                            masses=np.tile(self.molecule.masses, count),
                            formal_charges=None if charges is None else np.tile(charges, count),
                            guess_bonds=self.molecule.guess_bonds)
        cells = np.repeat(np.arange(count), n)
        return molecule, self._connect_cells(molecule, cells), cells

    @staticmethod
    def _local_material(styles, index, other=None, symbol=None):
        if styles is None or isinstance(styles, (bool, np.bool_)):
            return {}
        if not isinstance(styles, dict):
            styles = dict(enumerate(styles))
        material = {k: v for k, v in styles.items() if isinstance(k, str) and not isinstance(v, dict)}
        if symbol is not None:
            material.update(styles.get(symbol, {}))
        material.update(styles.get(index, {}))
        if other is not None:
            material.update(styles.get((other, index), {}))
            material.update(styles.get((index, other), {}))
        return material

    def _plot_cell_connections(self, network, bridges, cells, central, images, figure, backend, units):
        """Render boundary bonds through Molecule.plot with native split-half styling."""
        n = len(self.molecule.atoms)
        owners = []

        def endpoint_style(pair, symbols, material):
            i, j = pair
            settings = central if cells[i] == 0 else images
            # Native bond primitives invoke this once for each endpoint. Record
            # ownership in that same rendering order for return_objects grouping.
            if i < j:
                owners.append(int(max(cells[i], cells[j])))
            atom_style = self._local_material(settings.get("atom_style"), i % n, symbol=symbols[0])
            material.update({k: atom_style[k] for k in ("color", "glow") if k in atom_style})
            material.update(self._local_material(settings.get("bond_style"), i % n, j % n))
            modifier = material.pop("modifier", None)
            return material if modifier is None else modifier((i % n, j % n), symbols, material)

        options = images.copy()
        options.update(bonds=bridges, atom_style=False, rings=False,
                       bond_style={i: {"modifier": endpoint_style} for i in range(len(network.atoms))})
        # Per-atom radius options still refer to the stored motif's atom indices.
        for key in ("atom_radii", "atom_radius_scaling"):
            value = options.get(key)
            if isinstance(value, dict):
                options[key] = {k: v for k, v in value.items() if isinstance(k, str)} | {
                    k * n + i: v for k in range(len(network.atoms) // n)
                    for i, v in value.items() if isinstance(i, (int, np.integer))}
            elif value is not None and not np.isscalar(value):
                options[key] = np.tile(value, len(network.atoms) // n)
        rendered = network.plot(figure=figure, backend=backend, units=units, **options)
        return rendered, owners

    def plot(self, *, draw_cell=True, cell_style=None, unit_axes=False, unit_axes_style=None,
             unit_axes_label_style=None, unit_axes_origin=None, unit_axes_scale=1,
             show_periodic_images=None, periodic_images=None, periodic_image_styles=None,
             backend="x3d", units="Angstroms", **options):
        """Plot the stored molecule with its periodic boundary and optional neighbors.

        ``cell_style`` and ``unit_axes_style`` accept common styles, a/b/c overrides,
        or three style dictionaries. A cell ``radius`` draws cylindrical edges;
        otherwise edges are lines. Radii use the display units. ``unit_axes=True``
        adds labeled lattice-vector arrows along periodic directions only.

        Images form a full neighborhood outside the central cell: counts [2,0,1]
        mean +/-2 cells along a, none along b, +/-1 along c, including corners.
        An entry [negative,positive] sets the two sides independently: [2,0,[1,0]]
        adds one layer in -c and none in +c. Booleans mean counts 1/0.
        Finite directions never generate images. Supplying
        counts enables images unless ``show_periodic_images=False``. Shared image
        material styles default to transparency=.5; nested atom_style/bond_style
        entries override those defaults. Only the central cell gets axes/annotations.
        With bond inference enabled, copies are connected using the same Molecule
        rules as tile(); native split bonds inherit each endpoint's cell styling.
        With return_objects=True, image objects follow the central geometry's objects.
        """
        if show_periodic_images is not None and not isinstance(show_periodic_images, (bool, np.bool_)):
            raise ValueError("show_periodic_images must be a boolean")
        show_images = periodic_images is not None if show_periodic_images is None else show_periodic_images
        translations = self._image_translations(periodic_images) if show_images else np.empty((0, 3))
        scale = 1 if units is None else UnitsData.convert("BohrRadius", units)
        primitives = self._cell_annotations(scale, draw_cell, cell_style, unit_axes, unit_axes_style,
                                            unit_axes_label_style, unit_axes_origin, unit_axes_scale)
        if primitives or len(translations):
            if backend in ("jsmol", "rdkit", "2d", "rdkit3d", "svg", "svg2D") or options.get("mode") == "jsmol":
                raise ValueError("periodic boundaries, axes, and images require a 3D graphics backend")
        if primitives:
            annotation = options.pop("annotation_function", None)

            def annotate(plotter, index, geometry):
                return primitives + ([] if annotation is None else list(annotation(plotter, index, geometry)))
            options["annotation_function"] = annotate
            if backend in ("matplotlib3D", "plotly3D"):
                bounds = np.concatenate([np.asarray(self.molecule.coords), self.cell_vertices]) * scale
                options.setdefault("plot_range", np.stack([bounds.min(axis=0), bounds.max(axis=0)], axis=-1).tolist())
        if len(translations) and options.get("animate"):
            raise ValueError("periodic image overlays require a static plot")
        # Each image goes through the public molecule renderer, preserving atom/bond
        # styles and return conventions without duplicating lattice/other annotations.
        image_options = options.copy()
        for key in ("figure", "annotation_function", "comparison_styles", "animation_frame_styles"):
            image_options.pop(key, None)
        for key in ("principle_axes", "dipole", "mode_vectors", "draw_coords", "atom_text"):
            image_options[key] = None
        image_options.update(display_atom_numbers=False, animate=False)
        styles = copy.deepcopy(periodic_image_styles or {})
        atom_style = styles.pop("atom_style", None)
        bond_style = styles.pop("bond_style", None)
        common = dict(transparency=.5) | styles
        image_options["atom_style"] = self._image_component_style(options.get("atom_style"), common, atom_style)
        image_options["bond_style"] = self._image_component_style(options.get("bond_style"), common, bond_style)
        if backend in ("matplotlib", "matplotlib3D", "matplotlib3d"):
            image_options["clip_on"] = False
        network, bridges, cells = None, [], None
        if (len(translations) and self.molecule.guess_bonds and image_options["bond_style"] is not False
                and options.get("bonds") is not False and options.get("bond_style") is not False):
            network, bridges, cells = self._image_bond_network(translations, options.get("bonds"))
        if options.get("bonds") is None and self.molecule._bonds is not None and len(self.molecule._bonds) == 0:
            # Molecule's renderer expects False rather than an empty bond list.
            options["bonds"] = image_options["bonds"] = False
        result = self.molecule.plot(backend=backend, units=units, **options)
        figure = result[0] if options.get("return_objects") else result
        for translation in translations:
            image = self.molecule.plot(np.asarray(self.molecule.coords) + translation, figure=figure,
                                       backend=backend, units=units, **copy.deepcopy(image_options))
            if options.get("return_objects"):
                for target, extra in zip(result[1:], image[1:]):
                    target.extend(extra)
        if bridges:
            connections, owners = self._plot_cell_connections(network, bridges, cells, options,
                                                               image_options, figure, backend, units)
            if options.get("return_objects"):
                for owner, bond_objects in zip(owners, connections[2][0]):
                    if result[2][owner] is None:
                        result[2][owner] = []
                    result[2][owner].append(bond_objects)
                if connections[3][0]:
                    if result[3][0] is None:
                        result[3][0] = []
                    result[3][0].extend(connections[3][0])
        return result
