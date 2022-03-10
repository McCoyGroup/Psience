## <a id="Psience.Molecools.Vibrations.MolecularVibrations">MolecularVibrations</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L19)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L19?message=Update%20Docs)]
</div>



<a id="Psience.Molecools.Vibrations.MolecularVibrations.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, basis, freqs=None, init=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L21)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L21?message=Update%20Docs)]
</div>

Sets up a vibration for a Molecule object over the CoordinateSystem basis
- `molecule`: `AbstractMolecule`
    >No description...
- `init`: `None | CoordinateSet`
    >No description...
- `basis`: `MolecularNormalModes`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.molecule" class="docs-object-method">&nbsp;</a> 
```python
@property
molecule(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `CoordinateSet`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L76)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L76?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.displace" class="docs-object-method">&nbsp;</a> 
```python
displace(self, displacements=None, amt=0.1, n=1, which=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L79)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L79?message=Update%20Docs)]
</div>

Displaces along the vibrational mode specified by `which`
- `displacements`: `Any`
    >No description...
- `amt`: `Any`
    >No description...
- `n`: `Any`
    >No description...
- `which`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.visualize" class="docs-object-method">&nbsp;</a> 
```python
visualize(self, step_size=0.1, steps=(5, 5), which=0, anim_opts=None, mode='fast', **plot_args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L111)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L111?message=Update%20Docs)]
</div>


- `step_size`: `Any`
    >No description...
- `steps`: `Any`
    >No description...
- `which`: `Any`
    >No description...
- `anim_opts`: `Any`
    >No description...
- `mode`: `Any`
    >No description...
- `plot_args`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L175)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L175?message=Update%20Docs)]
</div>

Takes a slice of the modes
- `item`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.embed" class="docs-object-method">&nbsp;</a> 
```python
embed(self, frame): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L209)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L209?message=Update%20Docs)]
</div>


- `frame`: `MolecularTransformation`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.rescale" class="docs-object-method">&nbsp;</a> 
```python
rescale(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L224)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L224?message=Update%20Docs)]
</div>

Multiplies each mode by some scaling factor
- `phases`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L239)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L239?message=Update%20Docs)]
</div>



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/Molecools/Vibrations/MolecularVibrations.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/Molecools/Vibrations/MolecularVibrations.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/Molecools/Vibrations/MolecularVibrations.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/Molecools/Vibrations/MolecularVibrations.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L19?message=Update%20Docs)