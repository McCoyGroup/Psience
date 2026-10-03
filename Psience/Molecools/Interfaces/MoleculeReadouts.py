"""
Readout interface for `Molecule` (see `McUtils.Jupyter.Readouts`).

Kept out of `Molecule.py`: `Molecule.to_readout(...)` just builds a `MoleculeReadoutInterface`
around the molecule. Sections::

    overview     formula, atom count, charge, spin, mass, energy (+ point group on request)
    source       file path and format
    rdkit        SMILES, InChI, InChIKey, formula, exact mass (needs RDKit)
    structure    3D view of the molecule
    cartesians   atoms, masses and Cartesian coordinates
    internals    bond lengths, angles (and optionally dihedrals)
    bonds        bond list with orders and lengths (off by default)
    normal_modes frequencies, mode animations and displacement vectors (delegated to
                 `NormalModesReadoutInterface`)
"""

import importlib.util
import itertools
import os
import re

import numpy as np

from McUtils.Jupyter import (
    ReadoutAdapter, readout_section, ReadoutSectionUnavailable,
    ReadoutSection, ReadoutFields, ReadoutTable, ReadoutText, ReadoutScene,
    FieldSet, TabularData, Column,
)
import McUtils.Numputils as nput

__all__ = [
    "MoleculeReadoutInterface",
]


def _hill_formula(atoms):
    counts = {}
    for a in atoms:
        sym = re.match(r"[A-Z][a-z]?", str(a))
        sym = sym.group(0) if sym else str(a)
        counts[sym] = counts.get(sym, 0) + 1
    order = (["C", "H"] if "C" in counts else []) + sorted(k for k in counts if not ("C" in counts and k in ("C", "H")))
    return "".join(f"{k}{counts[k] if counts[k] > 1 else ''}" for k in order)


class MoleculeReadoutInterface(ReadoutAdapter):
    """
    **LLM Docstring**

    The readout sections for a `Molecule`; build one with `Molecule.to_readout(...)` or
    `MoleculeReadoutInterface(mol).to_readout(...)`. Coordinates are read in atomic units
    (Bohr, Hartree) and converted to the readout's display units while the readout is built.
    """
    readout_id = "molecule"

    @property
    def mol(self):
        return self.obj

    # ---- labels ------------------------------------------------------------------------- #
    def atom_labels(self):
        return [f"{a}{i + 1}" for i, a in enumerate(self.mol.atoms)]

    def get_readout_title(self):
        name = self.mol.name
        return name if name and name != "Unnamed" else _hill_formula(self.mol.atoms)

    def get_readout_subtitle(self):
        src = getattr(self.mol, "_src", None)
        if src and src.get("file"):
            fmt = f" ({src['mode']})" if src.get("mode") else ""
            return os.path.basename(str(src["file"])) + fmt
        return None

    def get_readout_meta(self):
        return {"atoms": list(self.mol.atoms)}

    # ---- sections ----------------------------------------------------------------------- #
    @readout_section("overview", title="Overview")
    def readout_overview(self, ctx, point_group=False):
        """Formula, size, charge, spin, total mass and (if present) the energy."""
        mol = self.mol
        fields = []
        if mol.name and mol.name != "Unnamed":
            fields.append(ctx.field("name", mol.name, label="Name"))
        fields.append(ctx.field("formula", _hill_formula(mol.atoms), label="Formula"))
        fields.append(ctx.field("num_atoms", int(len(mol.atoms)), label="Atoms", quantity="int"))
        charge = mol.charge
        fields.append(ctx.field("charge", None if charge is None else int(charge), label="Charge", quantity="charge"))
        spin = mol.spin
        if spin is not None:
            fields.append(ctx.field("spin", spin, label="Spin multiplicity"))
        fields.append(ctx.field("mass", float(np.sum(mol.masses)), quantity="mass", unit="AtomicMassUnits",
                                label="Total mass"))
        energy = mol.energy
        if energy is not None:
            fields.append(ctx.field("energy", float(energy), quantity="energy", unit="Hartrees", label="Energy"))
        if point_group:
            try:
                pg = mol.point_group
                fields.append(ctx.field("point_group", str(getattr(pg, "name", pg)), label="Point group"))
            except Exception as e:
                fields.append(ctx.field("point_group", None, label="Point group",
                                        description=f"failed: {type(e).__name__}: {e}"))
        return ReadoutFields(FieldSet(fields, name="overview"))

    def _has_source(self):
        src = getattr(self.mol, "_src", None)
        return (src is not None and src.get("file") is not None), "the molecule was not loaded from a file"

    @readout_section("source", title="Source", available="_has_source")
    def readout_source(self, ctx, absolute=True):
        """The file the molecule was loaded from."""
        src = self.mol._src
        path = str(src["file"])
        full = os.path.abspath(path) if absolute else path
        fields = [
            ctx.field("file", full, label="File", quantity="path"),
            ctx.field("name", os.path.basename(path), label="File name"),
            ctx.field("format", src.get("mode"), label="Format"),
        ]
        if os.path.isfile(path):
            fields.append(ctx.field("size", int(os.path.getsize(path)), label="Size (bytes)", quantity="int"))
        return ReadoutFields(FieldSet(fields, name="source"))

    def _has_rdkit(self):
        if importlib.util.find_spec("rdkit") is None:
            return False, "RDKit is not installed"
        return True, None

    @readout_section("rdkit", title="Identifiers (RDKit)", available="_has_rdkit")
    def readout_rdkit(self, ctx):
        """SMILES, InChI and InChIKey plus RDKit's formula and exact mass."""
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors, Descriptors
        mol = self.mol
        rd = mol.rdmol
        if rd is None:
            raise ReadoutSectionUnavailable("couldn't build an RDKit molecule")
        raw = rd.rdmol
        def safe(f):
            try:
                return f()
            except Exception:
                return None
        fields = [
            ctx.field("smiles", safe(lambda: mol.to_string("smi")), label="SMILES", quantity="identifier"),
            ctx.field("canonical_smiles", safe(lambda: Chem.MolToSmiles(Chem.RemoveHs(raw))),
                      label="Canonical SMILES", quantity="identifier"),
            ctx.field("inchi", safe(lambda: mol.to_string("inchi")), label="InChI", quantity="identifier"),
            ctx.field("inchi_key", safe(lambda: mol.to_string("inchi_key")), label="InChIKey", quantity="identifier"),
            ctx.field("formula", safe(lambda: rdMolDescriptors.CalcMolFormula(raw)), label="Formula"),
            ctx.field("exact_mass", safe(lambda: float(Descriptors.ExactMolWt(raw))), quantity="mass",
                      unit="AtomicMassUnits", label="Exact mass"),
            ctx.field("formal_charge", safe(lambda: int(Chem.GetFormalCharge(raw))), label="Formal charge",
                      quantity="charge"),
            ctx.field("num_bonds", safe(lambda: int(raw.GetNumBonds())), label="Bonds", quantity="int"),
            ctx.field("num_rings", safe(lambda: int(rdMolDescriptors.CalcNumRings(raw))), label="Rings", quantity="int"),
        ]
        return ReadoutFields(FieldSet(fields, name="rdkit"))

    @readout_section("structure", title="Structure")
    def readout_structure(self, ctx, html_backend="x3d", plot_options=None):
        """
        A 3D view: a mesh3D plot (used for PowerPoint `.glb` export and posters) and, for HTML,
        the native X3D plot (``html_backend='x3d'``) or the same mesh (``'mesh'``).
        """
        mol = self.mol
        plot_options = dict(plot_options or {})
        model = lambda: mol.plot(backend="mesh3D", **plot_options)
        html = None
        if html_backend == "x3d":
            def html(view=None):
                opts = dict(plot_options)
                if view is not None:  # share the camera used for the .glb and its poster
                    opts.setdefault("view_settings", view.x3d_viewpoint())
                return mol.plot(backend="x3d", **opts).figure.to_x3d()
        return ReadoutScene(model=model, html=html, id="view", animated=False)

    @readout_section("cartesians", title="Cartesian coordinates")
    def readout_cartesians(self, ctx):
        """Atom labels, masses and Cartesian coordinates."""
        mol = self.mol
        xyz = np.asarray(mol.coords)
        cols = [
            Column("index", np.arange(1, len(mol.atoms) + 1), label="#", quantity="index"),
            Column("atom", np.asarray(mol.atoms, dtype=str), label="Atom"),
            ctx.column("mass", np.asarray(mol.masses, dtype=float), quantity="mass", unit="AtomicMassUnits",
                       label="Mass"),
            ctx.column("x", xyz[:, 0], quantity="length", unit="BohrRadius", label="x"),
            ctx.column("y", xyz[:, 1], quantity="length", unit="BohrRadius", label="y"),
            ctx.column("z", xyz[:, 2], quantity="length", unit="BohrRadius", label="z"),
        ]
        return ReadoutTable(TabularData(cols, name="coordinates"))

    def _internal_specs(self, specs=None, dihedrals=False):
        mol = self.mol
        if specs is None:
            ints = mol.internals
            if isinstance(ints, dict) and ints.get("specs") is not None:
                specs = ints["specs"]
        if specs is None:
            bonds = [tuple(sorted((int(b[0]), int(b[1])))) for b in (mol.bonds or [])]
            nbrs = {}
            for i, j in bonds:
                nbrs.setdefault(i, set()).add(j)
                nbrs.setdefault(j, set()).add(i)
            specs = list(bonds)
            for j in sorted(nbrs):
                for i, k in itertools.combinations(sorted(nbrs[j]), 2):
                    specs.append((i, j, k))
            if dihedrals:
                for j, k in bonds:
                    for i in sorted(nbrs.get(j, ()) - {k}):
                        for l in sorted(nbrs.get(k, ()) - {j, i}):
                            specs.append((i, j, k, l))
        return [tuple(int(x) for x in s) for s in specs]

    def _has_internals(self):
        try:
            return len(self._internal_specs()) > 0, "no internal coordinates or bonds are defined"
        except Exception as e:
            return False, f"{type(e).__name__}: {e}"

    @readout_section("internals", title="Internal coordinates", available="_has_internals")
    def readout_internals(self, ctx, specs=None, dihedrals=False):
        """
        Bond lengths, angles and dihedrals: the molecule's own internal-coordinate specs, or
        (if none are set) coordinates generated from the bond graph. One table per kind, so
        each table has a single unit.
        """
        specs = self._internal_specs(specs, dihedrals=dihedrals)
        coords = np.asarray(self.mol.coords)
        values = np.asarray(nput.internal_coordinate_tensors(coords, specs, order=0)[0]).reshape(-1)
        labels = self.atom_labels()
        kinds = {2: ("bonds", "Bond lengths", "length", "BohrRadius", "r"),
                 3: ("angles", "Bond angles", "angle", "Radians", "a"),
                 4: ("dihedrals", "Dihedral angles", "angle", "Radians", "d")}
        tables = []
        for n, (name, title, quantity, unit, prefix) in kinds.items():
            idx = [i for i, s in enumerate(specs) if len(s) == n]
            if not idx:
                continue
            cols = [
                Column("coordinate", np.array([f"{prefix}({'–'.join(labels[a] for a in specs[i])})" for i in idx]),
                       label="Coordinate"),
                Column("atoms", np.array([" ".join(str(a + 1) for a in specs[i]) for i in idx]), label="Atoms"),
                ctx.column("value", values[idx], quantity=quantity, unit=unit, label="Value"),
            ]
            tables.append(ReadoutTable(TabularData(cols, name=name, title=title)))
        other = [s for s in specs if len(s) not in kinds]
        if other:
            tables.append(ReadoutText(f"{len(other)} coordinate specs of other kinds were skipped", role="note"))
        return tables

    def _has_bonds(self):
        return bool(self.mol.bonds), "no bonds are defined"

    @readout_section("bonds", title="Bonds", default=False, available="_has_bonds")
    def readout_bonds(self, ctx):
        """The bond list with bond orders and lengths."""
        bonds = self.mol.bonds
        labels = self.atom_labels()
        coords = np.asarray(self.mol.coords)
        i = np.array([int(b[0]) for b in bonds])
        j = np.array([int(b[1]) for b in bonds])
        order = np.array([float(b[2]) if len(b) > 2 else 1.0 for b in bonds])
        cols = [
            Column("atom_1", np.array([labels[a] for a in i]), label="Atom 1"),
            Column("atom_2", np.array([labels[a] for a in j]), label="Atom 2"),
            Column("order", order, label="Order", fmt="{:.1f}"),
            ctx.column("length", np.linalg.norm(coords[i] - coords[j], axis=1), quantity="length",
                       unit="BohrRadius", label="Length"),
        ]
        return ReadoutTable(TabularData(cols, name="bonds"))

    def _has_modes(self):
        try:
            nm = self.mol.normal_modes
            modes = nm.modes if nm is not None else None
            return modes is not None, "no normal modes are available"
        except Exception as e:
            return False, f"{type(e).__name__}: {e}"

    @readout_section("normal_modes", title="Normal modes", available="_has_modes", cost="expensive")
    def readout_normal_modes(self, ctx, include=None, exclude=None, **opts):
        """
        Delegates to `NormalModesReadoutInterface` (sections ``frequencies``, ``animations``,
        ``displacements``); nested options such as ``normal_modes={'animations': {'which': [0, 1]}}``
        or dotted includes (``'normal_modes.frequencies'``) are passed through.
        """
        from .ModeReadouts import NormalModesReadoutInterface
        modes = self.mol.normal_modes.modes.basis.to_new_modes()
        sub = NormalModesReadoutInterface(modes, molecule=self.mol)
        return ReadoutSection(*sub.get_readout_sections(include=include, exclude=exclude, ctx=ctx, **opts))
