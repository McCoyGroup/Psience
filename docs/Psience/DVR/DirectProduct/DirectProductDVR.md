## <a id="Psience.DVR.DirectProduct.DirectProductDVR">DirectProductDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L17?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.DirectProduct.DirectProductDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, dvrs_1D, zero_threshold=1e-14, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L18)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L18?message=Update%20Docs)]
</div>

  - `dvrs_1D`: `Iterable[AbstractDVR]`
    > a series of 1D DVRs that can provide the inputs we'll product together
  - `base_opts`: `Any`
    >


<a id="Psience.DVR.DirectProduct.DirectProductDVR.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L32?message=Update%20Docs)]
</div>
**LLM Docstring**

Debug string representation showing the class name, the component 1D DVRs, and the potential function.
  - `:returns`: `str`
    > string of the form `ClassName(dvrs, pot=potential_function)`


<a id="Psience.DVR.DirectProduct.DirectProductDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L47)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L47?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the full multi-dimensional DVR grid as the Cartesian product of each component 1D DVR's own grid, handling component DVRs (like `FiniteBasisDVR`) that additionally return a basis transformation alongside their grid points.
  - `domain`: `list[tuple] | None`
    > per-dimension domains to build fresh 1D grids with, instead of each component DVR's own stored grid
  - `divs`: `list[int] | None`
    > per-dimension division counts, paired with `domain`
  - `kwargs`: `dict`
    > extra options, unused
  - `:returns`: `np.ndarray | tuple[np.ndarray, list]`
    > the multi-dimensional grid, or `(grid, transformations)` if any component DVR returned a basis transformation


<a id="Psience.DVR.DirectProduct.DirectProductDVR.grid" class="docs-object-method">&nbsp;</a> 
```python
grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L77)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L77?message=Update%20Docs)]
</div>
**LLM Docstring**

Build (or retrieve) the multi-dimensional DVR grid, falling back to this DVR's stored `domain`/`divs` if not given explicitly, via `get_grid`.
  - `domain`: `list[tuple] | None`
    > per-dimension domains; defaults to `self.domain`
  - `divs`: `list[int] | None`
    > per-dimension division counts; defaults to `self.divs`
  - `kwargs`: `dict`
    > extra options forwarded to `get_grid`
  - `:returns`: `np.ndarray | tuple`
    > the multi-dimensional grid, or `(grid, transformations)`


<a id="Psience.DVR.DirectProduct.DirectProductDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, include_kinetic_coupling=True, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L104)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L104?message=Update%20Docs)]
</div>
**LLM Docstring**

Assemble the full multi-dimensional kinetic-energy operator as a sparse matrix, either as a simple Kronecker sum of the per-dimension 1D kinetic-energy operators (constant-mass case), or, when a (possibly coordinate-dependent) kinetic-coupling tensor `g`/`g_deriv` is supplied, by explicitly building the diagonal G-matrix-weighted kinetic terms, the pseudopotential-like `g_deriv` correction, and the off-diagonal momentum-momentum coupling terms between each pair of dimensions.
  - `grid`: `np.ndarray | tuple`
    > the multi-dimensional grid (or `(grid, transformations)` pair) to build the operator on; each dimension's subgrid is extracted from it
  - `mass`: `float | list[float] | None`
    > the mass (or per-dimension masses) to use if `g` isn't given
  - `hb`: `float | list[float]`
    > the value of hbar (or per-dimension values) to use
  - `g`: `list[list] | None`
    > the (possibly coordinate-dependent) kinetic-coupling tensor, `g[i][j]` giving the coupling between dimensions `i` and `j` as a constant or a callable of the flattened grid
  - `g_deriv`: `list | None`
    > the second-derivative correction terms for the diagonal `g` entries, per dimension
  - `logger`: `Logger | None`
    > logger used to report progress when kinetic coupling is being evaluated
  - `include_kinetic_coupling`: `bool`
    > whether to include the off-diagonal momentum-coupling terms (only relevant when `g` is given and has nonzero off-diagonal entries)
  - `kwargs`: `dict`
    > extra options, unused
  - `:returns`: `scipy.sparse.spmatrix | np.ndarray`
    > the assembled kinetic-energy matrix, as a sparse matrix (or dense array if it ends up sufficiently non-sparse)


<a id="Psience.DVR.DirectProduct.DirectProductDVR.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L307)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct/DirectProductDVR.py#L307?message=Update%20Docs)]
</div>
Computes the N-dimensional kinetic energy
  - `grid`: `Any`
    > 
  - `mass`: `Any`
    > 
  - `hb`: `Any`
    > 
  - `g`: `Any`
    > 
  - `g_deriv`: `Any`
    > 
  - `:returns`: `_`
    >
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DirectProduct/DirectProductDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DirectProduct/DirectProductDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DirectProduct/DirectProductDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DirectProduct/DirectProductDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L17?message=Update%20Docs)   
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