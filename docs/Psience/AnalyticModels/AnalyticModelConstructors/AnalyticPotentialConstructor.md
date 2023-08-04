## <a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor">AnalyticPotentialConstructor</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L19?message=Update%20Docs)]
</div>

Provides a set of symbolic potentials for use in models







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.morse" class="docs-object-method">&nbsp;</a> 
```python
morse(i, j, De=None, a=None, re=None, eq=None, w=None, wx=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L25)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L25?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a Morse potential
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.calc_morse" class="docs-object-method">&nbsp;</a> 
```python
calc_morse(De, a, r, re): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L54)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L54?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.harm" class="docs-object-method">&nbsp;</a> 
```python
harm(k, x, x_e): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L58)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L58?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.harmonic" class="docs-object-method">&nbsp;</a> 
```python
harmonic(*args, k=None, eq=None, qe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L61)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L61?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a Morse potential
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.lin" class="docs-object-method">&nbsp;</a> 
```python
lin(k, x, x_e): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L79)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L79?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.linear" class="docs-object-method">&nbsp;</a> 
```python
linear(*args, k=1, eq=None, xe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L82)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L82?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a linear function
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.pow" class="docs-object-method">&nbsp;</a> 
```python
pow(k, x, x_e, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L97)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L97?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.power" class="docs-object-method">&nbsp;</a> 
```python
power(*args, k=1, eq=None, n=None, xe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L100)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L100?message=Update%20Docs)]
</div>
Returns a fully symbolic form of a linear function
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.cos" class="docs-object-method">&nbsp;</a> 
```python
cos(*args, eq=None, qe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L116)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L116?message=Update%20Docs)]
</div>
Returns a fully symbolic form of cos
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.sin" class="docs-object-method">&nbsp;</a> 
```python
sin(*args, eq=None, qe=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L130)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L130?message=Update%20Docs)]
</div>
Returns a fully symbolic form of sin
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticPotentialConstructor.multiwell" class="docs-object-method">&nbsp;</a> 
```python
multiwell(*args, turning_points=None, origin=None, eq=None, minimum=0, depth=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L144)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.py#L144?message=Update%20Docs)]
</div>

  - `args`: `Any`
    > 
  - `turning_points`: `Any`
    > 
  - `depth`: `Any`
    > 
  - `:returns`: `_`
    >
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L19?message=Update%20Docs)   
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