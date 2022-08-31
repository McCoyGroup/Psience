## <a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel">AnalyticModel</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors.py#L264)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L264?message=Update%20Docs)]
</div>

Provides a symbolic representation of an analytically evaluatable Hamiltonian
which can be used to get derived expressions to evaluate.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
SympyExpr: SympyExpr
Potential: AnalyticPotentialConstructor
KE: AnalyticKineticEnergyConstructor
```
<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coordinates, potential, dipole=None, values=None, rotation=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L270)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L270?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L297)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L297?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.constants" class="docs-object-method">&nbsp;</a> 
```python
@property
constants(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L302)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L302?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.normal_modes" class="docs-object-method">&nbsp;</a> 
```python
normal_modes(self, dimensionless=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L306)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L306?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.to_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
to_normal_modes(self, dimensionless=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L327?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.get_VPT_expansions" class="docs-object-method">&nbsp;</a> 
```python
get_VPT_expansions(self, order=2, expansion_order=None, include_potential=None, include_gmatrix=None, include_pseudopotential=None, evaluate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L340)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L340?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.run_VPT" class="docs-object-method">&nbsp;</a> 
```python
run_VPT(self, order=2, states=2, return_analyzer=True, expansion_order=None, include_potential=None, include_gmatrix=None, include_pseudopotential=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L387)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L387?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.wrap_function" class="docs-object-method">&nbsp;</a> 
```python
wrap_function(self, expr, transform_coordinates=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L462)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L462?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.expand_potential" class="docs-object-method">&nbsp;</a> 
```python
expand_potential(self, order, lambdify=True, evaluate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L478)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L478?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.get_DVR_parameters" class="docs-object-method">&nbsp;</a> 
```python
get_DVR_parameters(self, expansion_order=None, lambdify=True, evaluate='constants'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L490)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L490?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.setup_DVR" class="docs-object-method">&nbsp;</a> 
```python
setup_DVR(self, domain=None, divs=None, use_normal_modes=False, expansion_order=None, **params): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L528)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L528?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, expr, mode='all', numericize=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L552)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L552?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L584)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L584?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.jacobian_inverse" class="docs-object-method">&nbsp;</a> 
```python
jacobian_inverse(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L600)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L600?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.g" class="docs-object-method">&nbsp;</a> 
```python
g(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L639)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L639?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.v" class="docs-object-method">&nbsp;</a> 
```python
v(self, order=2, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L661)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L661?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.vp" class="docs-object-method">&nbsp;</a> 
```python
vp(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L688)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L688?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.mu" class="docs-object-method">&nbsp;</a> 
```python
mu(self, order=1, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L710)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L710?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.morse" class="docs-object-method">&nbsp;</a> 
```python
morse(*args, De=None, a=None, re=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L20)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L20?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a Morse potential
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.harmonic" class="docs-object-method">&nbsp;</a> 
```python
harmonic(*args, k=None, qe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L40)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L40?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a Morse potential
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.linear" class="docs-object-method">&nbsp;</a> 
```python
linear(*args, k=1, xe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L56)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L56?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a linear function
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.sym" class="docs-object-method">&nbsp;</a> 
```python
sym(base, *args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L736)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L736?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.m" class="docs-object-method">&nbsp;</a> 
```python
m(i): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L739)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L739?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.r" class="docs-object-method">&nbsp;</a> 
```python
r(i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L742)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L742?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.a" class="docs-object-method">&nbsp;</a> 
```python
a(i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L745)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L745?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.t" class="docs-object-method">&nbsp;</a> 
```python
t(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L748)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L748?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.y" class="docs-object-method">&nbsp;</a> 
```python
y(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L751)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L751?message=Update%20Docs)]
</div>


<a id="McUtils.McUtils.Data.ConstantsData.UnitsDataHandler.convert" class="docs-object-method">&nbsp;</a> 
```python
convert(unit, target): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/McUtils/Data/ConstantsData/UnitsDataHandler.py#L396)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/McUtils/Data/ConstantsData/UnitsDataHandler.py#L396?message=Update%20Docs)]
</div>
Converts base unit into target using the scraped NIST data
  - `unit`: `Any`
    > 
  - `target`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.mass" class="docs-object-method">&nbsp;</a> 
```python
mass(atom): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L756)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L756?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticModel.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticModel.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L264?message=Update%20Docs)   
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