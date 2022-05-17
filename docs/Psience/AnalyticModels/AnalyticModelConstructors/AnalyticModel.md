## <a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel">AnalyticModel</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L224)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L224?message=Update%20Docs)]
</div>

Provides a symbolic representation of an analytically evaluatable Hamiltonian
which can be used to get derived expressions to evaluate.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
Potential: type
KE: type
```
<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coordinates, potential, dipole=None, values=None, rotation=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L230)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L230?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.normal_modes" class="docs-object-method">&nbsp;</a> 
```python
normal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L262)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L262?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.to_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
to_normal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L278)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L278?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.get_VPT_expansions" class="docs-object-method">&nbsp;</a> 
```python
get_VPT_expansions(self, order=2, evaluate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L291)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L291?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.run_VPT" class="docs-object-method">&nbsp;</a> 
```python
run_VPT(self, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L300)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L300?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, expr): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L323)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L323?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L352)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L352?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.jacobian_inverse" class="docs-object-method">&nbsp;</a> 
```python
jacobian_inverse(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L368)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L368?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.g" class="docs-object-method">&nbsp;</a> 
```python
g(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L407)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L407?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.v" class="docs-object-method">&nbsp;</a> 
```python
v(self, order=2, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L429)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L429?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.vp" class="docs-object-method">&nbsp;</a> 
```python
vp(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L456)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L456?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.mu" class="docs-object-method">&nbsp;</a> 
```python
mu(self, order=1, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L478)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L478?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.morse" class="docs-object-method">&nbsp;</a> 
```python
morse(*args, De=None, a=None, re=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L18)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L18?message=Update%20Docs)]
</div>

Returns a fully symbolic form of a Morse potential
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.harmonic" class="docs-object-method">&nbsp;</a> 
```python
harmonic(*args, k=None, qe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L38)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L38?message=Update%20Docs)]
</div>

Returns a fully symbolic form of a Morse potential
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.linear" class="docs-object-method">&nbsp;</a> 
```python
linear(*args, k=1, xe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L54)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L54?message=Update%20Docs)]
</div>

Returns a fully symbolic form of a linear function
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.sym" class="docs-object-method">&nbsp;</a> 
```python
sym(base, *args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L504)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L504?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.m" class="docs-object-method">&nbsp;</a> 
```python
m(i): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L507)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L507?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.r" class="docs-object-method">&nbsp;</a> 
```python
r(i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L510)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L510?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.a" class="docs-object-method">&nbsp;</a> 
```python
a(i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L513)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L513?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.t" class="docs-object-method">&nbsp;</a> 
```python
t(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L516)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L516?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.y" class="docs-object-method">&nbsp;</a> 
```python
y(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L519)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L519?message=Update%20Docs)]
</div>

<a id="McUtils.McUtils.Data.ConstantsData.UnitsDataHandler.convert" class="docs-object-method">&nbsp;</a> 
```python
convert(unit, target): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L396)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L396?message=Update%20Docs)]
</div>

Converts base unit into target using the scraped NIST data
- `unit`: `Any`
    >No description...
- `target`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.mass" class="docs-object-method">&nbsp;</a> 
```python
mass(atom): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L524)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L524?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticModel.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticModel.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticModel.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticModel.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/AnalyticModelConstructors.py#L224?message=Update%20Docs)