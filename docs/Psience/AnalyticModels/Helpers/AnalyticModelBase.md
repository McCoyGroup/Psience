## <a id="Psience.AnalyticModels.Helpers.AnalyticModelBase">AnalyticModelBase</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers.py#L41)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers.py#L41?message=Update%20Docs)]
</div>

Provides a base class for analytic models







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
sym: SympyShim
numeric_types: tuple
```
<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.take_derivs" class="docs-object-method">&nbsp;</a> 
```python
take_derivs(expr, vars): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L47)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L47?message=Update%20Docs)]
</div>
Takes derivatives of `expr` with respect to `vars` even if `expr` is an array
  - `expr`: `Any`
    > 
  - `vars`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.eval_exprs" class="docs-object-method">&nbsp;</a> 
```python
eval_exprs(expr, subs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L63)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L63?message=Update%20Docs)]
</div>
Evaluates `expr` with the given substitutions
  - `expr`: `Any`
    > 
  - `subs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbol_list" class="docs-object-method">&nbsp;</a> 
```python
symbol_list(names, instance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L98)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L98?message=Update%20Docs)]
</div>
Gets a list of symbols for `names` with a given instance number
  - `names`: `Any`
    > 
  - `instance`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_m" class="docs-object-method">&nbsp;</a> 
```python
symbolic_m(i): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L112)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L112?message=Update%20Docs)]
</div>
Provides a symbolic representation of a mass
  - `i`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbol" class="docs-object-method">&nbsp;</a> 
```python
symbol(base, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L123)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L123?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_r" class="docs-object-method">&nbsp;</a> 
```python
symbolic_r(i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L130)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L130?message=Update%20Docs)]
</div>
Provides a symbolic representation of a bond length
  - `i`: `Any`
    > 
  - `j`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_a" class="docs-object-method">&nbsp;</a> 
```python
symbolic_a(i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L145?message=Update%20Docs)]
</div>
Provides a symbolic representation of a bond angle
  - `i`: `Any`
    > 
  - `j`: `Any`
    > 
  - `k`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_t" class="docs-object-method">&nbsp;</a> 
```python
symbolic_t(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L162)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L162?message=Update%20Docs)]
</div>
Provides a symbolic representation of a dihedral angle
  - `i`: `Any`
    > 
  - `j`: `Any`
    > 
  - `k`: `Any`
    > 
  - `l`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_y" class="docs-object-method">&nbsp;</a> 
```python
symbolic_y(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L181)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L181?message=Update%20Docs)]
</div>
Provides a symbolic representation of a book angle
  - `i`: `Any`
    > 
  - `j`: `Any`
    > 
  - `k`: `Any`
    > 
  - `l`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.var" class="docs-object-method">&nbsp;</a> 
```python
var(*args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L198)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L198?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.reindex_symbol" class="docs-object-method">&nbsp;</a> 
```python
reindex_symbol(symbol, mapping, target_symbols=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L213)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L213?message=Update%20Docs)]
</div>
Changes the indices on symbols using the given mapping
  - `symbol`: `Any`
    > 
  - `mapping`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.lam" class="docs-object-method">&nbsp;</a> 
```python
lam(i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L250)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L250?message=Update%20Docs)]
</div>
Provides the `lambda` expression from Frederick and Woywood
  - `i`: `Any`
    > 
  - `j`: `Any`
    > 
  - `k`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.is_identity" class="docs-object-method">&nbsp;</a> 
```python
is_identity(A): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L269)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L269?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.transpose" class="docs-object-method">&nbsp;</a> 
```python
transpose(A): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L280)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L280?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.dot" class="docs-object-method">&nbsp;</a> 
```python
dot(a, b, axes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L285)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L285?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.contract" class="docs-object-method">&nbsp;</a> 
```python
contract(a, axes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L317)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L317?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.transform_coordinates" class="docs-object-method">&nbsp;</a> 
```python
transform_coordinates(rotation, coord_vec=None, coord_name_fmt='q{id}[{num}]'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/Helpers/AnalyticModelBase.py#L323)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers/AnalyticModelBase.py#L323?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/AnalyticModels/Helpers/AnalyticModelBase.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/AnalyticModels/Helpers/AnalyticModelBase.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/AnalyticModels/Helpers/AnalyticModelBase.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/AnalyticModels/Helpers/AnalyticModelBase.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/Helpers.py#L41?message=Update%20Docs)   
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