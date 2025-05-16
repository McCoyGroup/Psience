## <a id="Psience.AnalyticModels.Helpers.AnalyticModelBase">AnalyticModelBase</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/AnalyticModels/Helpers.py#L44)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L44?message=Update%20Docs)]
</div>

Provides a base class for analytic models







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
sym: SympyShim
```
<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.get_numeric_types" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_numeric_types(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L49)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L49?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.take_derivs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
take_derivs(cls, expr, vars): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L52)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L52?message=Update%20Docs)]
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
@classmethod
eval_exprs(cls, expr, subs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L68)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L68?message=Update%20Docs)]
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
@classmethod
symbol_list(cls, names, instance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L103)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L103?message=Update%20Docs)]
</div>
Gets a list of symbols for `names` with a given instance number
  - `names`: `Any`
    > 
  - `instance`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_x" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
symbolic_x(i): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L117)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L117?message=Update%20Docs)]
</div>
Provides a symbolic representation of a position
  - `i`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_n" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
symbolic_n(i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L128)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L128?message=Update%20Docs)]
</div>
Provides a symbolic representation of a normal to a plane
  - `i`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_m" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
symbolic_m(i): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L139)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L139?message=Update%20Docs)]
</div>
Provides a symbolic representation of a mass
  - `i`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbol" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
symbol(base, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L150)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L150?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.symbolic_r" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
symbolic_r(i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L157)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L157?message=Update%20Docs)]
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
@staticmethod
symbolic_a(i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L172)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L172?message=Update%20Docs)]
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
@staticmethod
symbolic_t(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L189)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L189?message=Update%20Docs)]
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
@staticmethod
symbolic_y(i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L208)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L208?message=Update%20Docs)]
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


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.infer_coord_type" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
infer_coord_type(cls, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L228)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L228?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.var" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
var(cls, *args, coord_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L244)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L244?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.reindex_symbol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
reindex_symbol(cls, symbol, mapping, target_symbols=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L264)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L264?message=Update%20Docs)]
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
@classmethod
lam(cls, i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L301)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L301?message=Update%20Docs)]
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
@classmethod
is_identity(cls, A): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L320)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L320?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.transpose" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
transpose(cls, A): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L331)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L331?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.dot" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
dot(cls, a, b, axes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L336)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L336?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.contract" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
contract(cls, a, axes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L368)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L368?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.Helpers.AnalyticModelBase.transform_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
transform_coordinates(cls, rotation, coord_vec=None, coord_name_fmt='q{id}[{num}]'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L374)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L374?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/Helpers.py#L44?message=Update%20Docs)   
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