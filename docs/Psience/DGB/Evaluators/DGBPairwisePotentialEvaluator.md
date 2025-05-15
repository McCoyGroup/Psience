## <a id="Psience.DGB.Evaluators.DGBPairwisePotentialEvaluator">DGBPairwisePotentialEvaluator</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators.py#L1682)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators.py#L1682?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DGB.Evaluators.DGBPairwisePotentialEvaluator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coords, pairwise_potential_functions, quadrature_degree=3, use_with_interpolation='ignored'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1683)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1683?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBPairwisePotentialEvaluator.get_bond_length_deltas" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_bond_length_deltas(cls, natoms, ndim, i, j, full=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1690)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1690?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBPairwisePotentialEvaluator.get_coordinate_bond_length_projection" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_bond_length_projection(self, i, j) -> 'tuple[np.ndarray, np.ndarray]': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1708)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1708?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBPairwisePotentialEvaluator.get_coordinate_change_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_change_transformation(self, coordinate_projection_data) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1711)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1711?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBPairwisePotentialEvaluator.get_bond_length_change_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_bond_length_change_transformation(self, overlap_data: 'OverlapGaussianData', i, j) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1726)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1726?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBPairwisePotentialEvaluator.wrap_distance_function" class="docs-object-method">&nbsp;</a> 
```python
wrap_distance_function(self, i, j, overlap_data: 'OverlapGaussianData', transformations, pairwise_function): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1734)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1734?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBPairwisePotentialEvaluator.evaluate_pairwise_contrib" class="docs-object-method">&nbsp;</a> 
```python
evaluate_pairwise_contrib(self, overlap_data: 'OverlapGaussianData', quadrature_degree=None, expansion_degree=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1791)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators/DGBPairwisePotentialEvaluator.py#L1791?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Evaluators/DGBPairwisePotentialEvaluator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Evaluators/DGBPairwisePotentialEvaluator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Evaluators/DGBPairwisePotentialEvaluator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Evaluators/DGBPairwisePotentialEvaluator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators.py#L1682?message=Update%20Docs)   
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