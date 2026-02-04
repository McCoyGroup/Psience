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


<a id="Psience.Molecools.Vibrations.MolecularVibrations.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L52)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L52?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.molecule" class="docs-object-method">&nbsp;</a> 
```python
@property
molecule(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L58)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L58?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L66)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L66?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L73)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L73?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L87)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L87?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.displace" class="docs-object-method">&nbsp;</a> 
```python
displace(self, displacements=None, amt=0.1, n=1, which=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L90)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L90?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L126)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L126?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L209)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L209?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L231)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L231?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L265)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L265?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L280)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L280?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L294)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L294?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L311)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularVibrations.py#L311?message=Update%20Docs)]
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