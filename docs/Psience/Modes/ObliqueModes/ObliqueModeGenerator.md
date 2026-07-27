## <a id="Psience.Modes.ObliqueModes.ObliqueModeGenerator">ObliqueModeGenerator</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/ObliqueModes.py#L11)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/ObliqueModes.py#L11?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Modes.ObliqueModes.ObliqueModeGenerator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, f, g, dimensionless=False, sel=None, frequency_scaled=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/ObliqueModes.py#L13)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/ObliqueModes.py#L13?message=Update%20Docs)]
</div>
**LLM Docstring**

Build an "oblique" mode transformation from a force-constant matrix `f` and a G-matrix (or mass vector) `g`: computes ordinary normal modes (removing translation/rotation), then finds the rotation `R` (via the SVD of the mode matrix) that makes the transformation as close as possible to an orthogonal one while still diagonalizing `f`/`g`, optionally after first putting `f`/`g` into a dimensionless form.
  - `f`: `np.ndarray`
    > the force-constant (potential Hessian) matrix
  - `g`: `np.ndarray`
    > the G-matrix, or a vector of atomic masses (broadcast into a diagonal inverse-mass G-matrix)
  - `dimensionless`: `bool`
    > whether to first rescale `f`/`g` into a dimensionless form before computing modes
  - `sel`: `Iterable[int] | None`
    > a subset of coordinate indices to restrict `f`/`g` to
  - `frequency_scaled`: `bool`
    > whether the underlying normal modes should be frequency-scaled (dimensionless)
  - `:returns`: `None`
    > None


<a id="Psience.Modes.ObliqueModes.ObliqueModeGenerator.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, dimensionless=True, sel=None, use_internals=None, frequency_scaled=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L68)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L68?message=Update%20Docs)]
</div>
**LLM Docstring**

Build an `ObliqueModeGenerator` for a molecule, using its internal-coordinate potential derivatives and G-matrix if internal coordinates are available (or requested), otherwise its Cartesian potential derivatives and atomic masses.
  - `mol`: `Molecule`
    > the molecule to build the oblique modes for
  - `dimensionless`: `bool`
    > whether to rescale into a dimensionless form before computing modes
  - `sel`: `Iterable[int] | None`
    > a subset of coordinate indices to restrict to
  - `use_internals`: `bool | None`
    > whether to build the modes in internal coordinates; defaults to whether the molecule has internal coordinates defined
  - `frequency_scaled`: `bool`
    > whether the underlying normal modes should be frequency-scaled
  - `:returns`: `ObliqueModeGenerator`
    > the constructed generator


<a id="Psience.Modes.ObliqueModes.ObliqueModeGenerator.run" class="docs-object-method">&nbsp;</a> 
```python
run(self, scaling_type='normal', remove_frequency_scaling=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/ObliqueModes/ObliqueModeGenerator.py#L109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/ObliqueModes/ObliqueModeGenerator.py#L109?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply the oblique-mode rotation (or its inverse, if `scaling_type='inverse'`) to transform `f` and `g` into the oblique basis, optionally undoing the initial dimensionless rescaling from the constructor.
  - `scaling_type`: `str`
    > `'normal'` to use the forward scaling, or `'inverse'` to swap the scaling/inverse-scaling roles
  - `remove_frequency_scaling`: `bool`
    > whether to undo the dimensionless rescaling applied in the constructor (if any) before returning
  - `:returns`: `tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]`
    > `(f, g, u, ui)` -- the transformed force-constant and G-matrices, and the (possibly rescaled) forward/inverse transformation matrices
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/ObliqueModes/ObliqueModeGenerator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/ObliqueModes/ObliqueModeGenerator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/ObliqueModes/ObliqueModeGenerator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/ObliqueModes/ObliqueModeGenerator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/ObliqueModes.py#L11?message=Update%20Docs)   
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