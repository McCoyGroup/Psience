"""
Unified scan infrastructure

`ScanManager` is the single entry point for both halves of a scan:

  * `.generate(scan_iterator, ...)` -- drains a `scan_iterator` (an iterable
    of `(index, values, atoms, coords)` steps), writes one job file per step
    into `{output_directory}/scan_data_{scan_id}/`, and records a
    `scan_info.json` manifest mapping each step to its file.
  * `.parse(molecular_property_extractor, ...)` -- reads that manifest back,
    rebuilds a `Molecule` from each step's electronic-structure output, runs
    `molecular_property_extractor` on it, and stacks the results into one
    tensor per property key, shaped like the scan grid.

The default implementation targets ORCA (`OrcaJob` for generation, `.out`
files for parsing). To support another package, subclass and override
`default_job_builder` and `default_output_file_generator`:

    class GaussianScanManager(ScanManager):
        job_file_ext = ".gjf"
        output_file_ext = ".log"

        def default_job_builder(self, atoms, coords, charge, **opts):
            return GaussianJob(atoms=atoms, cartesians=..., charge=charge, **opts)

        # default_output_file_generator's extension swap (.gjf -> .log) already
        # works unchanged as long as job_file_ext/output_file_ext are set above

A few free functions are provided to build `scan_iterator`s for the two scan
types worked out previously (a local Cartesian atomic-position scan, and an
internal-coordinate scan), via `structure_scan_iterator`, which zips an
N-dimensional grid of coordinate values with a `structure_generator` into the
`(index, values, atoms, coords)` stream `generate` expects.
"""

import os
import json
import itertools
import warnings
import abc
from collections.abc import Callable # I hate this

import numpy as np

from McUtils.Data import UnitsData
from McUtils.ExternalPrograms import OrcaJob
from ..Molecools import Molecule

__all__ = [
    "ScanManager",
    "iter_scan_grid",
    "structure_scan_iterator",
    "cartesian_position_scan_structures",
    "internal_coordinate_index",
    "internal_coordinate_scan_structures",
]


# region scan-grid / structure-generator helpers


def iter_scan_grid(domains, expand_domains=True):
    """
    Walks every grid point of an N-dimensional scan.

    :param domains: list of `(min, max, npoints)` specs, one per scanned
        coordinate
    :return: yields `(index, values)` pairs, where `index` is the integer
        multi-index of the grid point (e.g. `(2, 5)`) and `values` is the
        corresponding tuple of coordinate values
    """
    if expand_domains:
        axis_values = [np.linspace(*d) for d in domains]
    else:
        axis_values = domains
    shape = tuple(len(v) for v in axis_values)
    dtype = np.dtype(type(axis_values[0]))
    n = len(domains)
    for index in itertools.product(*(range(n) for n in shape)):
        values = np.fromiter(
            (axis_values[k][i] for k, i in enumerate(index)),
            dtype=dtype,
            count=n
        )
        yield index, values

def structure_scan_iterator(domains, structure_generator: Callable[[np.ndarray], Molecule],):
    """
    Bridges the grid/structure-generator world into the `scan_iterator`
    interface `ScanManager.generate` expects, by zipping `iter_scan_grid`
    with a `values -> (atoms, coords)` callable.

    :param domains: list of `(min, max, npoints)` specs, one per scanned
        coordinate
    :param structure_generator: `values -> (atoms, coords)` callable, e.g.
        from `cartesian_position_scan_structures` or
        `internal_coordinate_scan_structures`
    :return: yields `(index, values, atoms, coords)` for every grid point
    """
    for index, values in iter_scan_grid(domains):
        yield index, structure_generator(values)

def molecule_scan_iterator(mol:Molecule, geometries):
    geometries = np.asarray(geometries)
    base_shape = geometries.shape[:-2]
    for index in np.ndindex(*base_shape):
        yield index, mol.modify(coords=geometries[index])

def molecule_displaced_geometries_iterator(
        mol:Molecule,
        domains,
        which,
        **etc
):
    mol.get_scan_coordinates(
        domains,
        which=which,
    )


def axis_embedding_iterator(embedding_axes):
    ...

def cartesian_position_scan_structures(mol, atom_index, embedding_axes, origin, normalize_axes=True):
    """
    Builds a `structure_generator` that displaces a single atom of `mol` over
    a local Cartesian embedding, holding every other atom fixed at its
    reference position.

    :param mol: reference molecule
    :type mol: Psience.Molecools.Molecule
    :param atom_index: 0-based index of the atom to scan
    :param embedding_axes: `(k, 3)` array of local axis directions (`k=2` for
        the usual 2D scan)
    :param origin: `(3,)` Cartesian point (same units as `mol.coords`) that
        anchors the all-zero grid point
    :param normalize_axes: whether to orthonormalize `embedding_axes` first
    :return: `structure_generator(values) -> (atoms, coords)` callable
    """
    origin = np.asarray(origin, dtype=float)
    axes = normalize_embedding_axes(embedding_axes) if normalize_axes else np.asarray(embedding_axes, dtype=float)
    base_coords = np.asarray(mol.coords)
    atoms = mol.atoms

    def structure_generator(values):
        offset = np.tensordot(np.asarray(values, dtype=float), axes, axes=[[0], [0]])
        coords = base_coords.copy()
        coords[atom_index] = origin + offset
        return atoms, coords

    return structure_generator


def internal_coordinate_index(atom_row, component, natoms):
    """
    Builds the flat index (into the `(natoms, 3)` internal-coordinate array,
    mirroring the Cartesian `(natoms, 3)` layout) for a given Z-matrix row and
    coordinate component.

    :param atom_row: 0-based row of the Z-matrix (i.e. atom index in scan order)
    :param component: `0` for the bond length, `1` for the bond angle, `2`
        for the dihedral. Rows 0-2 involve embedding placeholders (see
        `MolecularZMatrixCoordinateSystem.embedding_coords`) rather than
        physical internal coordinates.
    :param natoms: number of atoms in the molecule
    :return: flat index suitable for the `which` argument below
    :rtype: int
    """
    return int(np.ravel_multi_index((atom_row, component), (natoms, 3)))


def internal_coordinate_scan_structures(mol, which, internals=None):
    """
    Builds a `structure_generator` that sets a set of `mol`'s internal
    coordinates to given absolute values and returns the corresponding
    Cartesian structure. Bond lengths are expected in Bohr and angles/
    dihedrals in radians (`mol.internal_coordinates`' native units).

    :param mol: reference molecule
    :type mol: Psience.Molecools.Molecule
    :param which: flat indices (into the `(natoms, 3)` internal-coordinate
        array) of the coordinates to scan; see `internal_coordinate_index`
    :param internals: a Z-matrix spec (`(natoms, 4)` int array) or other spec
        understood by `Molecule.internals`; only used if `mol.internals`
        isn't already set
    :return: `structure_generator(values) -> (atoms, coords)` callable
    """
    if mol.internals is None:
        if internals is None:
            raise ValueError("mol has no internal coordinates set; pass `internals=...`")
        mol.internals = internals

    atoms = mol.atoms

    def structure_generator(values):
        coords = mol.get_displaced_coordinates(
            np.asarray(values, dtype=float),
            which=which,
            use_internals='convert',
            shift=False
        )
        return atoms, np.asarray(coords)

    return structure_generator

# endregion


class ScanManager(metaclass=abc.ABCMeta):
    """
    Generates and parses a Cartesian- or internal-coordinate scan's job files.

    Targets ORCA by default; subclass and override `default_job_builder` /
    `default_output_file_generator` (and `job_file_ext`/`output_file_ext`) to
    support another electronic-structure package.
    """

    job_file_ext: str
    output_file_ext: str

    def __init__(self, output_directory, scan_id=None, job_prefix="scan", index_format="03d"):
        """
        :param output_directory: parent directory in which the scan's storage
            directory lives (or is created)
        :param scan_id: label used to build the storage directory
            `scan_data_{scan_id}`; if `None`, `output_directory` itself is
            used as the storage directory
        :param job_prefix: filename prefix for generated job files
        :param index_format: format spec applied to each integer in a step's
            multi-index when building filenames
        """
        self.output_directory = output_directory
        self.scan_id = scan_id
        self.job_prefix = job_prefix
        self.index_format = index_format

    @property
    def scan_dir(self):
        """Directory jobs are written to / read from."""
        if self.scan_id is None:
            return self.output_directory
        return os.path.join(self.output_directory, f"scan_data_{self.scan_id}")

    @property
    def scan_info_file(self):
        return os.path.join(self.scan_dir, "scan_info.json")

    # region generation

    @abc.abstractmethod
    def default_job_builder(self, mol, level_of_theory='r2scan-3c', MaxCore=3500, **etc):
        ...

    def generate(
            self,
            scan_iterator,
            charge=0,
            job_builder=None,
            coord_labels=None,
            extra_info=None,
            job_kwargs=None
    ):
        """
        Drains `scan_iterator`, writes one job file per step into
        `self.scan_dir`, and records a `scan_info.json` manifest mapping each
        step's multi-index and coordinate values to its job file.

        :param scan_iterator: iterable/generator yielding
            `(index, values, atoms, coords)` for every scan step (see
            `structure_scan_iterator`)
        :param charge: total molecular charge, forwarded to
            `default_job_builder` (ignored if `job_builder` is supplied)
        :param job_builder: `(atoms, coords) -> job` callable; defaults to
            `self.default_job_builder` bound with `charge`/`job_kwargs`
        :param coord_labels: optional human-readable names for the scanned
            coordinates, stored in the manifest
        :param extra_info: extra data folded into `scan_info.json` (e.g.
            `origin`/`embedding_axes`, or the `which`/`zmatrix` used)
        :param job_kwargs: extra keywords forwarded to `default_job_builder`
            (ignored if `job_builder` is supplied)
        :return: path to `self.scan_dir`
        :rtype: str
        """
        if job_builder is None:
            kwargs = job_kwargs or {}
            job_builder = lambda atoms, coords: self.default_job_builder(atoms, coords, charge, **kwargs)

        os.makedirs(self.scan_dir, exist_ok=True)

        steps = []
        for index, values, atoms, coords in scan_iterator:
            job = job_builder(atoms, coords)
            fname = (
                    self.job_prefix + "_"
                    + "_".join(f"{i:{self.index_format}}" for i in index)
                    + self.job_file_ext
            )
            job.write(os.path.join(self.scan_dir, fname))
            steps.append({
                "index": list(index),
                "values": [float(v) for v in values],
                "file": fname
            })

        info = {
            "scan_id": self.scan_id,
            "coord_labels": coord_labels,
            "shape": self._infer_shape(steps),
            "steps": steps
        }
        if extra_info:
            info.update(extra_info)

        with open(self.scan_info_file, "w") as info_file:
            json.dump(info, info_file, indent=2)

        return self.scan_dir

    @staticmethod
    def _infer_shape(steps):
        if not steps:
            return []
        ndim = len(steps[0]["index"])
        return [max(s["index"][d] for s in steps) + 1 for d in range(ndim)]

    # endregion

    # region parsing

    def load_scan_info(self):
        """Loads this scan's `scan_info.json` manifest."""
        with open(self.scan_info_file) as info_file:
            return json.load(info_file)

    def default_output_file_generator(self, input_file):
        """
        Default `output_file_generator`: swaps the input job file's extension
        (`self.job_file_ext`) for the electronic-structure output extension
        (`self.output_file_ext`). Override in a subclass for anything fancier
        (different directories, remote fetches, etc).

        :param input_file: path to the input job file, as recorded in
            `scan_info.json`
        :return: path to the corresponding output file
        :rtype: str
        """
        root, _ = os.path.splitext(input_file)
        return root + self.output_file_ext

    def load_molecules(
            self,
            output_file_generator=None,
            molecule_loader=None,
            scan_info=None,
            skip_missing=True
    ):
        """
        Rebuilds a `Molecule` for every completed step of the scan.

        :param output_file_generator: `input_file_path -> output_file_path`
            callable; defaults to `self.default_output_file_generator`
        :param molecule_loader: `output_file_path -> Molecule` callable;
            defaults to `Molecule.from_file`
        :param scan_info: pre-loaded manifest (loaded from disk if omitted)
        :param skip_missing: if `True`, steps whose output is missing/
            unreadable are skipped with a warning rather than raising
        :return: `dict` mapping each step's multi-index (as a tuple) to its
            `Molecule`
        :rtype: dict
        """

        if output_file_generator is None:
            output_file_generator = self.default_output_file_generator
        if molecule_loader is None:
            molecule_loader = Molecule.from_file
        if scan_info is None:
            scan_info = self.load_scan_info()

        mols = {}
        for step in scan_info["steps"]:
            index = tuple(step["index"])
            input_path = os.path.join(self.scan_dir, step["file"])
            output_path = output_file_generator(input_path)
            try:
                mols[index] = molecule_loader(output_path)
            except (FileNotFoundError, OSError, ValueError) as e:
                if skip_missing:
                    warnings.warn(f"skipping scan step {index} ({output_path}): {e}")
                    continue
                raise

        return mols

    def parse(
            self,
            molecular_property_extractor,
            output_file_generator=None,
            molecule_loader=None,
            scan_info=None,
            skip_missing=True,
            fill_value=np.nan
    ):
        """
        Rebuilds the `Molecule` for every completed scan step, runs
        `molecular_property_extractor` on each, and stacks the results into
        one tensor per property key, shaped like the scan grid.

        :param molecular_property_extractor: `Molecule -> dict[str, np.ndarray]`
            callable; every returned dict must use the same set of keys, and
            a given key's array must have the same shape at every step
        :param output_file_generator: `input_file_path -> output_file_path`
            callable; defaults to `self.default_output_file_generator`
        :param molecule_loader: `output_file_path -> Molecule` callable;
            defaults to `Molecule.from_file`
        :param scan_info: pre-loaded manifest (loaded from disk if omitted)
        :param skip_missing: if `True`, steps whose output is missing/
            unreadable are skipped (leaving `fill_value` in the corresponding
            tensor slots) rather than raising
        :param fill_value: value used for grid points that were never
            populated (e.g. failed or still-running jobs)
        :return: `dict` mapping each property key to an `np.ndarray` of shape
            `scan_info["shape"] + property_shape`
        :rtype: dict
        """
        from Psience.Molecools import Molecule

        if output_file_generator is None:
            output_file_generator = self.default_output_file_generator
        if molecule_loader is None:
            molecule_loader = Molecule.from_file
        if scan_info is None:
            scan_info = self.load_scan_info()

        grid_shape = tuple(scan_info["shape"])

        tensors = {}
        populated = np.zeros(grid_shape, dtype=bool)

        for step in scan_info["steps"]:
            index = tuple(step["index"])
            input_path = os.path.join(self.scan_dir, step["file"])
            output_path = output_file_generator(input_path)

            try:
                mol = molecule_loader(output_path)
                props = molecular_property_extractor(mol)
            except (FileNotFoundError, OSError, ValueError) as e:
                if skip_missing:
                    warnings.warn(f"skipping scan step {index} ({output_path}): {e}")
                    continue
                raise

            for key, val in props.items():
                val = np.asarray(val)
                if key not in tensors:
                    dtype = (
                        val.dtype if np.issubdtype(val.dtype, np.floating)
                        else np.result_type(val.dtype, type(fill_value))
                    )
                    tensors[key] = np.full(grid_shape + val.shape, fill_value, dtype=dtype)
                tensors[key][index] = val

            populated[index] = True

        if not np.all(populated):
            missing = np.argwhere(~populated)
            warnings.warn(f"{len(missing)} / {populated.size} scan points have no extracted properties")

        return tensors

    # endregion

class ORCAScanManager(ScanManager):
    job_file_ext = ".inp"
    output_file_ext = ".out"

    def default_job_builder(self, mol, level_of_theory='r2scan-3c', MaxCore=3500, **etc):
        """
        Wraps a single Cartesian structure as an `OrcaJob`. `coords` is
        assumed to be in Bohr and is converted to Angstroms for ORCA.
        Override in a subclass to target a different package.

        :param atoms: atom symbols
        :param coords: `(natoms, 3)` Cartesian coordinates, in Bohr
        :param charge: total molecular charge
        :param level_of_theory: ORCA level-of-theory keyword
        :param MaxCore: per-core memory in MB
        :param etc: additional `OrcaJob` options
        :return: an un-written job object supporting `.write(path)`
        """
        return OrcaJob(
            "TIGHTSCF",
            level_of_theory=level_of_theory,
            MaxCore=MaxCore,
            atoms=mol.atom,
            cartesians=np.asarray(mol.coords) * UnitsData.convert("BohrRadius", "Angstroms"),
            charge=mol.charge,
            spin=mol.spin,
            **etc
        )