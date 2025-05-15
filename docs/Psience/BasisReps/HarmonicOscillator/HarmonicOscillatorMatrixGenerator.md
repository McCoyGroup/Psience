## <a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator">HarmonicOscillatorMatrixGenerator</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator.py#L327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator.py#L327?message=Update%20Docs)]
</div>

1D evaluator for terms looking like `x`, `p`, `q`, etc.
All of the overall `(-i)^N` info is in the `ProdOp` class that's expected to hold this.
Only maintains phase info & calculates elements.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
state_cache_size: int
default_evaluator_mode: str
```
<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, terms, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L336)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L336?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L354)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L354?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_cache(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L362)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L362?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.set_cache_size" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
set_cache_size(cls, new_size): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L365)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L365?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.load_cached" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
load_cached(cls, terms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L369)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L369?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.selection_rules" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L381)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L381?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L385)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L385?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.state_pair_hash" class="docs-object-method">&nbsp;</a> 
```python
state_pair_hash(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L389)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L389?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.pull_state_groups" class="docs-object-method">&nbsp;</a> 
```python
pull_state_groups(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L400)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L400?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.evaluate_state_terms" class="docs-object-method">&nbsp;</a> 
```python
evaluate_state_terms(self, states, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L414)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L414?message=Update%20Docs)]
</div>
Evaluates terms coming from different state excitations.
Doesn't do any pre-filtering, since that's expected to be in the caller.
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.load_generator" class="docs-object-method">&nbsp;</a> 
```python
load_generator(self, a, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L442)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L442?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.get_paths" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_paths(cls, sizes, change): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L473)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L473?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.get_path_poly" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_path_poly(cls, path, parities=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L499)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L499?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.poly_coeffs" class="docs-object-method">&nbsp;</a> 
```python
poly_coeffs(self, delta, shift=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L565)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L565?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.get_poly_coeffs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_poly_coeffs(cls, terms, delta, shift=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L580)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L580?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.poly_term_generator" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
poly_term_generator(cls, terms, delta, shift=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L587)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L587?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.rho_term_generator" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
rho_term_generator(cls, a, N, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L643)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L643?message=Update%20Docs)]
</div>
Returns a function to be called on a quantum number to get the coefficient associated with exciting that mode by `a` quanta over
`N` steps w/ phase info coming from where the momenta-terms are.
  - `a`: `Any`
    > 
  - `N`: `Any`
    > 
  - `i_phase`: `Any`
    > 
  - `is_complex`: `Any`
    > 
  - `sel`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.rho" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
rho(cls, phases, paths, ni): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L707)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L707?message=Update%20Docs)]
</div>

  - `phases`: `Any`
    > 
  - `paths`: `Any`
    > 
  - `ni`: `Any`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator.py#L327?message=Update%20Docs)   
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