# <a id="Psience.Data.ScanManager">Psience.Data.ScanManager</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager.py#L1)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager.py#L1?message=Update%20Docs)]
</div>
    
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

### Members
<div class="container alert alert-secondary bg-light">
  <div class="row">
   <div class="col" markdown="1">
[ScanManager](ScanManager/ScanManager.md)   
</div>
   <div class="col" markdown="1">
[shape_scan_iterator](ScanManager/shape_scan_iterator.md)   
</div>
   <div class="col" markdown="1">
[scan_iterator](ScanManager/scan_iterator.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[molecule_scan_geometries_iterator](ScanManager/molecule_scan_geometries_iterator.md)   
</div>
   <div class="col" markdown="1">
[molecule_displaced_geometries_iterator](ScanManager/molecule_displaced_geometries_iterator.md)   
</div>
   <div class="col" markdown="1">
[molecule_atom_position_scan_iterator](ScanManager/molecule_atom_position_scan_iterator.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>













---


<div markdown="1" class="text-secondary">
<div class="container">
  <div class="row">
   <div class="col" markdown="1">
**Feedback**   
</div>
   <div class="col" markdown="1">
**Examples**   
</div>
   <div class="col" markdown="1">
**Templates**   
</div>
   <div class="col" markdown="1">
**Documentation**   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[Bug](https://github.com/McCoyGroup/Psience/issues/new?title=Documentation%20Improvement%20Needed)/[Request](https://github.com/McCoyGroup/Psience/issues/new?title=Example%20Request)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Data/ScanManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Data/ScanManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Data/ScanManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Data/ScanManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager.py#L1?message=Update%20Docs)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>
</div>