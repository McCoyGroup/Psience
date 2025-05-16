## <a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel">AnalyticModel</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors.py#L595)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L595?message=Update%20Docs)]
</div>

Provides a symbolic representation of an analytically evaluatable Hamiltonian
which can be used to get derived expressions to evaluate.
Supplies methods to automatically run DVR and VPT calculations from the model
specifications as well.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
SympyExpr: SympyExpr
Potential: AnalyticPotentialConstructor
KE: AnalyticKineticEnergyConstructor
NamespaceContext: NamespaceContext
```
<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coordinates, potential, dipole=None, values=None, rotation=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L605)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L605?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.from_potential" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_potential(cls, potential, dipole=None, values=None, rotation=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L632)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L632?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L650)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L650?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.constants" class="docs-object-method">&nbsp;</a> 
```python
@property
constants(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L655)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L655?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.normal_modes" class="docs-object-method">&nbsp;</a> 
```python
normal_modes(self, dimensionless=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L659)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L659?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.to_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
to_normal_modes(self, dimensionless=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L689)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L689?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.get_VPT_expansions" class="docs-object-method">&nbsp;</a> 
```python
get_VPT_expansions(self, order=2, expansion_order=None, include_potential=None, include_gmatrix=None, include_pseudopotential=None, evaluate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L702)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L702?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.run_VPT" class="docs-object-method">&nbsp;</a> 
```python
run_VPT(self, order=2, states=2, return_analyzer=True, expansion_order=None, include_potential=None, include_gmatrix=None, include_pseudopotential=None, atoms=None, coords=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L753)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L753?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.prep_lambda_expr" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_lambda_expr(cls, base_coords, expr, dummify=False, rewrite_trig=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L915)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L915?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.lambdify" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
lambdify(cls, coord_vec, expr, coordinate_transform=None, mode=None, dummify=False, rewrite_trig=True, lambdify_arrays=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L933)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L933?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.wrap_function" class="docs-object-method">&nbsp;</a> 
```python
wrap_function(self, expr, transform_coordinates=True, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L959)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L959?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.expand_potential" class="docs-object-method">&nbsp;</a> 
```python
expand_potential(self, order, lambdify=True, evaluate=True, contract=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L978)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L978?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.get_DVR_parameters" class="docs-object-method">&nbsp;</a> 
```python
get_DVR_parameters(self, expansion_order=None, lambdify=True, evaluate='constants'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L993)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L993?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.setup_DVR" class="docs-object-method">&nbsp;</a> 
```python
setup_DVR(self, domain=None, divs=None, use_normal_modes=False, expansion_order=None, potential_function=None, **params): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1034)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1034?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, expr, mode='all', numericize=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1061)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1061?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.parse_symbol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
parse_symbol(self, sym): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1087)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1087?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, order=0, evaluate=False, lambdify=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1094)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1094?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.jacobian_inverse" class="docs-object-method">&nbsp;</a> 
```python
jacobian_inverse(self, order=0, evaluate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1112)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1112?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.g" class="docs-object-method">&nbsp;</a> 
```python
g(self, order=0, evaluate=False, lambdify=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1152)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1152?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.v" class="docs-object-method">&nbsp;</a> 
```python
v(self, order=2, evaluate=False, lambdify=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1176)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1176?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.vp" class="docs-object-method">&nbsp;</a> 
```python
vp(self, order=0, evaluate=False, lambdify=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1205)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1205?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.mu" class="docs-object-method">&nbsp;</a> 
```python
mu(self, order=1, evaluate=False, lambdify=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1229)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1229?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.morse" class="docs-object-method">&nbsp;</a> 
```python
morse(i, j, De=None, a=None, re=None, eq=None, w=None, wx=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L26)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L26?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a Morse potential
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.harmonic" class="docs-object-method">&nbsp;</a> 
```python
harmonic(*args, k=None, eq=None, qe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L62)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L62?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a Morse potential
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.linear" class="docs-object-method">&nbsp;</a> 
```python
linear(*args, k=1, eq=None, xe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L83)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L83?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a linear function
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.power" class="docs-object-method">&nbsp;</a> 
```python
power(*args, k=1, eq=None, n=None, xe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L101)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L101?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a linear function
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.sin" class="docs-object-method">&nbsp;</a> 
```python
sin(*args, eq=None, qe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L131)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L131?message=Update%20Docs)]
</div>
Returns a fully symbolic form of sin
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.cos" class="docs-object-method">&nbsp;</a> 
```python
cos(*args, eq=None, qe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L117)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L117?message=Update%20Docs)]
</div>
Returns a fully symbolic form of cos
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.sym" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
sym(self, base, *args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1300)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1300?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.m" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
m(self, i): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1303)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1303?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.r" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
r(self, i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1306)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1306?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.a" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
a(self, i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1309)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1309?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.t" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
t(self, i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1312)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1312?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.y" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
y(self, i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1315)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1315?message=Update%20Docs)]
</div>


<a id="McUtils.Data.ConstantsData.UnitsDataHandler.convert" class="docs-object-method">&nbsp;</a> 
```python
convert(unit, target): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/ConstantsData/UnitsDataHandler.py#L515)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/ConstantsData/UnitsDataHandler.py#L515?message=Update%20Docs)]
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
@staticmethod
mass(atom): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1320)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1320?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.molecular_potential" class="docs-object-method">&nbsp;</a> 
```python
molecular_potential(self, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1324)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1324?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.molecular_dipole" class="docs-object-method">&nbsp;</a> 
```python
molecular_dipole(self, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1328)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1328?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticModel.molecular_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
molecular_gmatrix(self, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1332)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticModel.py#L1332?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L595?message=Update%20Docs)   
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