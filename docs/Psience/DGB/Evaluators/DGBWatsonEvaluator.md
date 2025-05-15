## <a id="Psience.DGB.Evaluators.DGBWatsonEvaluator">DGBWatsonEvaluator</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators.py#L1203)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators.py#L1203?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DGB.Evaluators.DGBWatsonEvaluator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, modes, coriolis_inertia_function): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators/DGBWatsonEvaluator.py#L1206)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators/DGBWatsonEvaluator.py#L1206?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBWatsonEvaluator.annoying_coriolis_term" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
annoying_coriolis_term(n, u, m, v, Xc, Dx, Sc, Sp, Gi, Gj, DG): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1217)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1217?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBWatsonEvaluator.annoying_coriolis_momentum_term" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
annoying_coriolis_momentum_term(n, u, m, v, Xc, r, Jp, Dx, Sc, Sp, DG): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1229)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1229?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBWatsonEvaluator.annoying_imaginary_momentum_term" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
annoying_imaginary_momentum_term(n, u, m, v, Xc, r, Jp, Dx, Sc, Sp, DG): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1243)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1243?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBWatsonEvaluator.evaluate_coriolis_contrib" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
evaluate_coriolis_contrib(cls, coriolis_tensors, overlap_data: 'OverlapGaussianData'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1257)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1257?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBWatsonEvaluator.evaluate_watson_term" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
evaluate_watson_term(cls, B_e, overlap_data: 'OverlapGaussianData'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1367)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1367?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBWatsonEvaluator.evaluate_ke" class="docs-object-method">&nbsp;</a> 
```python
evaluate_ke(self, overlap_data: 'OverlapGaussianData', logger=None, include_diagonal_contribution=True, include_coriolis_coupling=True, include_watson_term=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators/DGBWatsonEvaluator.py#L1391)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators/DGBWatsonEvaluator.py#L1391?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Evaluators/DGBWatsonEvaluator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Evaluators/DGBWatsonEvaluator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Evaluators/DGBWatsonEvaluator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Evaluators/DGBWatsonEvaluator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators.py#L1203?message=Update%20Docs)   
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