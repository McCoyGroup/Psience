## <a id="Psience.Molecools.Vibrations.MolecularVibrations">MolecularVibrations</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L20)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L20?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Molecools.Vibrations.MolecularVibrations.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, basis, freqs=None, init=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L22)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L22?message=Update%20Docs)]
</div>
Sets up a vibration for a Molecule object over the CoordinateSystem basis
  - `molecule`: `AbstractMolecule`
    > 
  - `init`: `None | CoordinateSet`
    > 
  - `basis`: `MolecularNormalModes`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.change_mol" class="docs-object-method">&nbsp;</a> 
```python
change_mol(self, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L44)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L44?message=Update%20Docs)]
</div>
**LLM Docstring**

Rebind this set of vibrations to a different molecule, re-expressing the basis for the new molecule while keeping the same frequencies and reference coordinates.
  - `mol`: `AbstractMolecule`
    > the new molecule to associate with the vibrations
  - `:returns`: `MolecularVibrations`
    > a new `MolecularVibrations` built from `mol`, with `basis` changed via `basis.change_mol(mol)` and the same `freqs`/`init` as `self`


<a id="Psience.Molecools.Vibrations.MolecularVibrations.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L62)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L62?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the underlying vibrational basis (a `MolecularNormalModes` object) used to describe the vibrations.
  - `basis`: `MolecularNormalModes`
    > (setter only) the new basis to store
  - `:returns`: `MolecularNormalModes`
    > (getter) the stored basis


<a id="Psience.Molecools.Vibrations.MolecularVibrations.molecule" class="docs-object-method">&nbsp;</a> 
```python
@property
molecule(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L88)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L88?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule associated with these vibrations. Setting it also propagates the new molecule to `self.basis.molecule`.
  - `mol`: `AbstractMolecule`
    > (setter only) the new molecule to associate
  - `:returns`: `AbstractMolecule`
    > (getter) the stored molecule


<a id="Psience.Molecools.Vibrations.MolecularVibrations.freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L116)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L116?message=Update%20Docs)]
</div>
**LLM Docstring**

Frequencies associated with the vibrations. Returns the explicitly stored frequencies if present; otherwise falls back to `self.basis.freqs` when the basis defines that attribute.
  - `:returns`: `np.ndarray | None`
    > the vibrational frequencies, or `None` if neither `self` nor `self.basis` has them


<a id="Psience.Molecools.Vibrations.MolecularVibrations.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L131)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L131?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L145?message=Update%20Docs)]
</div>
**LLM Docstring**

Number of vibrational modes, taken from the number of columns of the basis's mode matrix.
  - `:returns`: `int`
    > the number of vibrational modes


<a id="Psience.Molecools.Vibrations.MolecularVibrations.displace" class="docs-object-method">&nbsp;</a> 
```python
displace(self, displacements=None, amt=0.1, n=1, which=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L156)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L156?message=Update%20Docs)]
</div>
Displaces along the vibrational mode specified by `which`
  - `displacements`: `Any`
    > 
  - `amt`: `Any`
    > 
  - `n`: `Any`
    > 
  - `which`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.visualize" class="docs-object-method">&nbsp;</a> 
```python
visualize(self, step_size=5, steps=(2, 2), which=0, anim_opts=None, mode='fast', **plot_args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L192)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L192?message=Update%20Docs)]
</div>

  - `step_size`: `Any`
    > 
  - `steps`: `Any`
    > 
  - `which`: `Any`
    > 
  - `anim_opts`: `Any`
    > 
  - `mode`: `Any`
    > 
  - `plot_args`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.to_widget" class="docs-object-method">&nbsp;</a> 
```python
to_widget(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L283)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L283?message=Update%20Docs)]
</div>
**LLM Docstring**

Build (and cache) an interactive Jupyter widget for browsing through the vibrational modes: a menu to select which mode (`which`), paired with a live display that calls `self.visualize` for the selected mode.
  - `:returns`: `JHTML.Div | None`
    > the constructed `JHTML.Div` widget, or `None` if a widget was already built and cached in `self._widg`


<a id="Psience.Molecools.Vibrations.MolecularVibrations.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L323)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L323?message=Update%20Docs)]
</div>
Takes a slice of the modes
  - `item`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.embed" class="docs-object-method">&nbsp;</a> 
```python
embed(self, frame): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L357)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L357?message=Update%20Docs)]
</div>

  - `frame`: `MolecularTransformation`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.rescale" class="docs-object-method">&nbsp;</a> 
```python
rescale(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L372)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L372?message=Update%20Docs)]
</div>
Multiplies each mode by some scaling factor
  - `phases`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.rotate" class="docs-object-method">&nbsp;</a> 
```python
rotate(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L386)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L386?message=Update%20Docs)]
</div>
Multiplies each mode by some scaling factor
  - `phases`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L403)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L403?message=Update%20Docs)]
</div>
**LLM Docstring**

Debug string representation showing the class name, basis, and molecule.
  - `:returns`: `str`
    > string of the form `ClassName(basis, molecule)`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Vibrations/MolecularVibrations.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Vibrations/MolecularVibrations.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Vibrations/MolecularVibrations.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Vibrations/MolecularVibrations.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L20?message=Update%20Docs)   
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