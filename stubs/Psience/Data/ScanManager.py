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
__all__ = ['ScanManager', 'shape_scan_iterator', 'scan_iterator', 'molecule_scan_geometries_iterator', 'molecule_displaced_geometries_iterator', 'molecule_atom_position_scan_iterator']

def shape_scan_iterator(base_shape, zigzag=False):
    ...

def scan_iterator(domains, expand_domains=True, index_iterator=None, zigzag=False):
    ...

def molecule_scan_iterator(mol: Molecule, geometries, index_iterator=None, zigzag=False, return_molecules=True):
    ...

def molecule_scan_geometries_iterator(mol: Molecule, domains, which, return_values=False, return_molecules=True, index_iterator=None, zigzag=False, coordinate_generator=None, **etc):
    ...

def molecule_displaced_geometries_iterator(mol: Molecule, displacement_positions, which, return_molecules=True, return_values=False, index_iterator=None, zigzag=False, coordinate_generator=None, **etc):
    ...

def molecule_atom_position_scan_iterator(mol, atom_indices, domains, which=None, embedding=None, return_molecules=True, **iterator_options):
    ...

class ScanManager:
    job_file_ext: str
    output_file_ext: str

    def __init__(self, output_directory, scan_id=None, job_prefix='scan', index_format='03d'):
        ...
    scan_data_template = 'scan_data_{scan_id}'

    @property
    def scan_dir(self):
        """Directory jobs are written to / read from."""
        ...
    info_filename = 'scan_info.json'

    @property
    def scan_info_file(self):
        ...

    def default_job_builder(self, mol, *, job_type, commands=None, **etc):
        ...
    job_file_template = '{prefix}_{index}.{ext}'

    def generate(self, scan_iterator, job_builder=None, coord_labels=None, extra_info=None, overwrite=False, append=False, job_prefix=None, job_file_ext=None, job_type=None, **job_kwargs):
        ...

    @staticmethod
    def _infer_shape(steps):
        ...

    def load_scan_info(self):
        """Loads this scan's `scan_info.json` manifest."""
        ...

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
        ...

    def load_molecules(self, output_file_generator=None, molecule_loader=None, scan_info=None, skip_missing=True):
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
        ...

    def parse(self, molecular_property_extractor, output_file_generator=None, molecule_loader=None, scan_info=None, skip_missing=True, fill_value=np.nan):
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
        ...