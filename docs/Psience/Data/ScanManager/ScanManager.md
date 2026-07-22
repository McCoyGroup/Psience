## <a id="Psience.Data.ScanManager.ScanManager">ScanManager</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager.py#L196)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager.py#L196?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
scan_data_template: str
info_filename: str
job_file_template: str
```
<a id="Psience.Data.ScanManager.ScanManager.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, output_directory, scan_id=None, job_prefix='scan', index_format='03d'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager.py#L200)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager.py#L200?message=Update%20Docs)]
</div>


<a id="Psience.Data.ScanManager.ScanManager.scan_dir" class="docs-object-method">&nbsp;</a> 
```python
@property
scan_dir(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager/ScanManager.py#L207)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager/ScanManager.py#L207?message=Update%20Docs)]
</div>
Directory jobs are written to / read from.


<a id="Psience.Data.ScanManager.ScanManager.scan_info_file" class="docs-object-method">&nbsp;</a> 
```python
@property
scan_info_file(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager/ScanManager.py#L215)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager/ScanManager.py#L215?message=Update%20Docs)]
</div>


<a id="Psience.Data.ScanManager.ScanManager.default_job_builder" class="docs-object-method">&nbsp;</a> 
```python
default_job_builder(self, mol, *, job_type, commands=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager/ScanManager.py#L221)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager/ScanManager.py#L221?message=Update%20Docs)]
</div>


<a id="Psience.Data.ScanManager.ScanManager.generate" class="docs-object-method">&nbsp;</a> 
```python
generate(self, scan_iterator, job_builder=None, coord_labels=None, extra_info=None, overwrite=False, append=False, job_prefix=None, job_file_ext=None, job_type=None, **job_kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager/ScanManager.py#L233)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager/ScanManager.py#L233?message=Update%20Docs)]
</div>


<a id="Psience.Data.ScanManager.ScanManager.load_scan_info" class="docs-object-method">&nbsp;</a> 
```python
load_scan_info(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager/ScanManager.py#L310)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager/ScanManager.py#L310?message=Update%20Docs)]
</div>
Loads this scan's `scan_info.json` manifest.


<a id="Psience.Data.ScanManager.ScanManager.default_output_file_generator" class="docs-object-method">&nbsp;</a> 
```python
default_output_file_generator(self, input_file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager/ScanManager.py#L315)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager/ScanManager.py#L315?message=Update%20Docs)]
</div>
Default `output_file_generator`: swaps the input job file's extension
(`self.job_file_ext`) for the electronic-structure output extension
(`self.output_file_ext`). Override in a subclass for anything fancier
(different directories, remote fetches, etc).
  - `input_file`: `Any`
    > path to the input job file, as recorded in
    `scan_info.json`
  - `:returns`: `str`
    > path to the corresponding output file


<a id="Psience.Data.ScanManager.ScanManager.load_molecules" class="docs-object-method">&nbsp;</a> 
```python
load_molecules(self, output_file_generator=None, molecule_loader=None, scan_info=None, skip_missing=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager/ScanManager.py#L330)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager/ScanManager.py#L330?message=Update%20Docs)]
</div>
Rebuilds a `Molecule` for every completed step of the scan.
  - `output_file_generator`: `Any`
    > `input_file_path -> output_file_path`
    callable; defaults to `self.default_output_file_generator`
  - `molecule_loader`: `Any`
    > `output_file_path -> Molecule` callable;
    defaults to `Molecule.from_file`
  - `scan_info`: `Any`
    > pre-loaded manifest (loaded from disk if omitted)
  - `skip_missing`: `Any`
    > if `True`, steps whose output is missing/
    unreadable are skipped with a warning rather than raising
  - `:returns`: `dict`
    > `dict` mapping each step's multi-index (as a tuple) to its
    `Molecule`


<a id="Psience.Data.ScanManager.ScanManager.parse" class="docs-object-method">&nbsp;</a> 
```python
parse(self, molecular_property_extractor, output_file_generator=None, molecule_loader=None, scan_info=None, skip_missing=True, fill_value=nan): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/ScanManager/ScanManager.py#L374)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager/ScanManager.py#L374?message=Update%20Docs)]
</div>
Rebuilds the `Molecule` for every completed scan step, runs
`molecular_property_extractor` on each, and stacks the results into
one tensor per property key, shaped like the scan grid.
  - `molecular_property_extractor`: `Any`
    > `Molecule -> dict[str, np.ndarray]`
    callable; every returned dict must use the same set of keys, and
    a given key's array must have the same shape at every step
  - `output_file_generator`: `Any`
    > `input_file_path -> output_file_path`
    callable; defaults to `self.default_output_file_generator`
  - `molecule_loader`: `Any`
    > `output_file_path -> Molecule` callable;
    defaults to `Molecule.from_file`
  - `scan_info`: `Any`
    > pre-loaded manifest (loaded from disk if omitted)
  - `skip_missing`: `Any`
    > if `True`, steps whose output is missing/
    unreadable are skipped (leaving `fill_value` in the corresponding
    tensor slots) rather than raising
  - `fill_value`: `Any`
    > value used for grid points that were never
    populated (e.g. failed or still-running jobs)
  - `:returns`: `dict`
    > `dict` mapping each property key to an `np.ndarray` of shape
    `scan_info["shape"] + property_shape`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Data/ScanManager/ScanManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Data/ScanManager/ScanManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Data/ScanManager/ScanManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Data/ScanManager/ScanManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/ScanManager.py#L196?message=Update%20Docs)   
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