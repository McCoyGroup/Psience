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
import shutil
import warnings

import numpy as np
import itertools

import McUtils.Devutils as dev
import McUtils.Numputils as nput
import McUtils.Iterators as itut
from McUtils.ExternalPrograms import ExternalProgramJob

from ..Molecools import Molecule

__all__ = [
    "ScanManager",
    "shape_scan_iterator",
    "scan_iterator",
    "molecule_scan_geometries_iterator",
    "molecule_displaced_geometries_iterator",
    "molecule_atom_position_scan_iterator",
]


# region scan-grid / structure-generator helpers

def shape_scan_iterator(base_shape, zigzag=False):
    if not zigzag:
        index_iterator = np.ndindex(*base_shape)
    else:
        shape_iterator = [range(s) for s in base_shape]
        index_iterator = itut.zigzag_product(*shape_iterator)
    yield from index_iterator

def scan_iterator(domains, expand_domains=True, index_iterator=None, zigzag=False):
    if expand_domains:
        axis_values = [np.linspace(*d) for d in domains]
    else:
        axis_values = domains
    shape = tuple(len(v) for v in axis_values)
    dtype = np.dtype(type(axis_values[0]))
    n = len(domains)
    if index_iterator is None:
        index_iterator = shape_scan_iterator(shape, zigzag=zigzag)
    for index in index_iterator:
        values = np.fromiter(
            (axis_values[k][i] for k, i in enumerate(index)),
            dtype=dtype,
            count=n
        )
        yield index, values

def molecule_scan_iterator(mol:Molecule, geometries, index_iterator=None, zigzag=False, return_molecules=True):
    geometries = np.asarray(geometries)
    base_shape = geometries.shape[:-2]
    if index_iterator is None:
        index_iterator = shape_scan_iterator(base_shape, zigzag=zigzag)
    if return_molecules:
        for index in index_iterator:
            yield index, mol.modify(coords=geometries[index])
    else:
        for index in index_iterator:
            yield index, geometries[index]

def molecule_scan_geometries_iterator(
        mol:Molecule,
        domains,
        which,
        return_values=False,
        return_molecules=True,
        index_iterator=None,
        zigzag=False,
        coordinate_generator=None,
        **etc
):
    if coordinate_generator is None:
        coordinate_generator = mol.get_scan_coordinates
    geoms = coordinate_generator(
        domains,
        which=which,
        **etc
    )
    iterator = molecule_scan_iterator(mol, geoms,
                                               index_iterator=index_iterator,
                                               zigzag=zigzag,
                                               return_molecules=return_molecules)
    if return_values:
        if index_iterator is not None:
            index_iterator = itertools.tee(index_iterator)
        value_iterator = scan_iterator(
            domains,
            zigzag=zigzag,
            index_iterator=index_iterator
        )
        base_iterator = iterator
        iterator = (
            (idx, v, m)
            for (_, v), (idx, m) in zip(value_iterator, base_iterator)
        )
    return geoms.shape, iterator

def molecule_displaced_geometries_iterator(
        mol: Molecule,
        displacement_positions,
        which,
        return_molecules=True,
        return_values=False,
        index_iterator=None,
        zigzag=False,
        coordinate_generator=None,
        **etc
):
    if coordinate_generator is None:
        coordinate_generator = mol.get_displaced_coordinates
    return molecule_scan_geometries_iterator(
        mol,
        displacement_positions,
        which,
        return_molecules=return_molecules,
        return_values=return_values,
        index_iterator=index_iterator,
        zigzag=zigzag,
        coordinate_generator=coordinate_generator,
        **etc
    )

def molecule_atom_position_scan_iterator(mol, atom_indices, domains, which=None, embedding=None,
                                       return_molecules=True,
                                       **iterator_options
                                       ):
    if nput.is_int(atom_indices):
        atom_indices = [atom_indices]
    if embedding is None:
        _, _, embedding = mol.fragment_embedding(atom_indices, return_axes=True)
    elif nput.is_int(embedding[0]):
        _, _, embedding = mol.fragment_embedding(embedding, return_axes=True)
    else:
        embedding = np.asanyarray(embedding)

    expansion = np.zeros((3, 3 * len(mol.atoms)))
    for i in atom_indices:
        expansion[:, 3 * i:3 * (i + 1)] = embedding.T

    if nput.is_numeric(domains[0]):
        domains = [domains]

    if which is None:
        which = list(range(len(domains)))

    return molecule_scan_geometries_iterator(
        mol,
        domains,
        which,
        coordinate_expansion=[expansion],
        return_molecules=return_molecules,
        **iterator_options
    )


class ScanManager:
    job_file_ext: str
    output_file_ext: str

    def __init__(self, output_directory, scan_id=None, job_prefix="scan", index_format="03d"):
        self.output_directory = output_directory
        self.scan_id = scan_id
        self.job_prefix = job_prefix
        self.index_format = index_format

    scan_data_template = "scan_data_{scan_id}"
    @property
    def scan_dir(self):
        """Directory jobs are written to / read from."""
        if self.scan_id is None:
            return self.output_directory
        return os.path.join(self.output_directory, self.scan_data_template.format(scan_id=self.scan_id))

    info_filename = "scan_info.json"
    @property
    def scan_info_file(self):
        return os.path.join(self.scan_dir, self.info_filename)

    # region generation

    def default_job_builder(self, mol, *, job_type, commands=None, **etc):
        if commands is None:
            commands = []
        elif isinstance(commands, str):
            commands = [commands]
        return mol.setup_job(
            job_type,
            *commands,
            **etc
        )

    job_file_template = "{prefix}_{index}.{ext}"
    def generate(
            self,
            scan_iterator,
            job_builder=None,
            coord_labels=None,
            extra_info=None,
            overwrite=False,
            append=False,
            job_prefix=None,
            job_file_ext=None,
            job_type=None,
            **job_kwargs
    ):
        if job_builder is None:
            job_builder = lambda mol, **kwargs: self.default_job_builder(mol, **kwargs)

        if overwrite and os.path.exists(self.scan_dir):
            shutil.rmtree(self.scan_dir)
        os.makedirs(self.scan_dir, exist_ok=append)

        if job_file_ext is None:
            job_file_ext = getattr(self, 'job_file_ext', None)
            if job_type is not None:
                job_cls, _ = ExternalProgramJob.resolve(job_type)
                job_file_ext = job_cls.extension

        if job_prefix is None:
            job_prefix = self.job_prefix

        if job_type is not None:
            job_kwargs["job_type"] = job_type
        steps = []
        for index, values, mol in scan_iterator:
            job = job_builder(mol, **job_kwargs)
            fname = self.job_file_template.format(
                prefix=job_prefix,
                index="_".join(f"{i:{self.index_format}}" for i in index),
                ext=job_file_ext.strip(".")
            )
            job.write(os.path.join(self.scan_dir, fname))
            steps.append({
                "index": index.tolist() if hasattr(index, "tolist") else index,
                "values": values.tolist() if hasattr(values, "tolist") else index,
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

        dev.write_json(
            self.scan_info_file,
            info,
            indent=2
        )

        return self.scan_dir, {
            'info': info,
            'file': self.scan_info_file
        }, steps

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

# class ORCAScanManager(ScanManager):
#     job_file_ext = ".inp"
#     output_file_ext = ".out"
#
#     def default_job_builder(self, mol, level_of_theory='r2scan-3c', MaxCore=3500, **etc):
#         """
#         Wraps a single Cartesian structure as an `OrcaJob`. `coords` is
#         assumed to be in Bohr and is converted to Angstroms for ORCA.
#         Override in a subclass to target a different package.
#
#         :param atoms: atom symbols
#         :param coords: `(natoms, 3)` Cartesian coordinates, in Bohr
#         :param charge: total molecular charge
#         :param level_of_theory: ORCA level-of-theory keyword
#         :param MaxCore: per-core memory in MB
#         :param etc: additional `OrcaJob` options
#         :return: an un-written job object supporting `.write(path)`
#         """
#         return OrcaJob(
#             "TIGHTSCF",
#             level_of_theory=level_of_theory,
#             MaxCore=MaxCore,
#             atoms=mol.atom,
#             cartesians=np.asarray(mol.coords) * UnitsData.convert("BohrRadius", "Angstroms"),
#             charge=mol.charge,
#             spin=mol.spin,
#             **etc
#         )