## <a id="Psience.AnalyticModels.Helpers.AnalyticModelBase">AnalyticModelBase</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L41)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L41?message=Update%20Docs)]
</div>

Provides a base class for analytic models

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
sym: SympyShim
numeric_types: tuple
```
<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.take_derivs" class="docs-object-method">&nbsp;</a> 
```python
take_derivs(expr, vars): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L47)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L47?message=Update%20Docs)]
</div>

Takes derivatives of `expr` with respect to `vars` even if `expr` is an array
- `expr`: `Any`
    >No description...
- `vars`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.eval_exprs" class="docs-object-method">&nbsp;</a> 
```python
eval_exprs(expr, subs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L63)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L63?message=Update%20Docs)]
</div>

Evaluates `expr` with the given substitutions
- `expr`: `Any`
    >No description...
- `subs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbol_list" class="docs-object-method">&nbsp;</a> 
```python
symbol_list(names, instance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L98)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L98?message=Update%20Docs)]
</div>

Gets a list of symbols for `names` with a given instance number
- `names`: `Any`
    >No description...
- `instance`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_m" class="docs-object-method">&nbsp;</a> 
```python
symbolic_m(i): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L112)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L112?message=Update%20Docs)]
</div>

Provides a symbolic representation of a mass
- `i`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbol" class="docs-object-method">&nbsp;</a> 
```python
symbol(base, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L123)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L123?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_r" class="docs-object-method">&nbsp;</a> 
```python
symbolic_r(i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L130)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L130?message=Update%20Docs)]
</div>

Provides a symbolic representation of a bond length
- `i`: `Any`
    >No description...
- `j`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_a" class="docs-object-method">&nbsp;</a> 
```python
symbolic_a(i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L145)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L145?message=Update%20Docs)]
</div>

Provides a symbolic representation of a bond angle
- `i`: `Any`
    >No description...
- `j`: `Any`
    >No description...
- `k`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_t" class="docs-object-method">&nbsp;</a> 
```python
symbolic_t(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L162)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L162?message=Update%20Docs)]
</div>

Provides a symbolic representation of a dihedral angle
- `i`: `Any`
    >No description...
- `j`: `Any`
    >No description...
- `k`: `Any`
    >No description...
- `l`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_y" class="docs-object-method">&nbsp;</a> 
```python
symbolic_y(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L181)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L181?message=Update%20Docs)]
</div>

Provides a symbolic representation of a book angle
- `i`: `Any`
    >No description...
- `j`: `Any`
    >No description...
- `k`: `Any`
    >No description...
- `l`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.var" class="docs-object-method">&nbsp;</a> 
```python
var(*args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L198)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L198?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.reindex_symbol" class="docs-object-method">&nbsp;</a> 
```python
reindex_symbol(symbol, mapping, target_symbols=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L213)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L213?message=Update%20Docs)]
</div>

Changes the indices on symbols using the given mapping
- `symbol`: `Any`
    >No description...
- `mapping`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.lam" class="docs-object-method">&nbsp;</a> 
```python
lam(i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L250)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L250?message=Update%20Docs)]
</div>

Provides the `lambda` expression from Frederick and Woywood
- `i`: `Any`
    >No description...
- `j`: `Any`
    >No description...
- `k`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.is_identity" class="docs-object-method">&nbsp;</a> 
```python
is_identity(A): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L269)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L269?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.transpose" class="docs-object-method">&nbsp;</a> 
```python
transpose(A): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L280)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L280?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.dot" class="docs-object-method">&nbsp;</a> 
```python
dot(a, b, axes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L285)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L285?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.contract" class="docs-object-method">&nbsp;</a> 
```python
contract(a, axes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L317)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L317?message=Update%20Docs)]
</div>

<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.transform_coordinates" class="docs-object-method">&nbsp;</a> 
```python
transform_coordinates(rotation, coord_vec=None, coord_name_fmt='q{id}[{num}]'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L323)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L323?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/AnalyticModels/Helpers/AnalyticModelBase.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/AnalyticModels/Helpers/AnalyticModelBase.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/AnalyticModels/Helpers/AnalyticModelBase.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/AnalyticModels/Helpers/AnalyticModelBase.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L41?message=Update%20Docs)