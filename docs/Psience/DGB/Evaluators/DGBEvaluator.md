## <a id="Psience.DGB.Evaluators.DGBEvaluator">DGBEvaluator</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Evaluators.py#L334)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators.py#L334?message=Update%20Docs)]
</div>

An object that supports evaluating matrix elements in a distributed Gaussian basis.
Provides support for integrating a function via quadrature or as an expansion in a polynomial tensors







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DGB.Evaluators.DGBEvaluator.get_inverse_covariances" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_inverse_covariances(cls, alphas, transformations): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L340)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L340?message=Update%20Docs)]
</div>
Transforms the alphas into proper inverse covariance matrices.
Chosen so that in the case that the transformations, Q, diagonalize S we can write
QT S Q = A
  - `:returns`: `_`
    >


<a id="Psience.DGB.Evaluators.DGBEvaluator.get_covariances" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_covariances(cls, alphas, transformations): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L363)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L363?message=Update%20Docs)]
</div>
Transforms the alphas into proper inverse covariance matrices.
Chosen so that in the case that the transformations, Q, diagonalize S we can write
QT S Q = A
  - `:returns`: `_`
    >


<a id="Psience.DGB.Evaluators.DGBEvaluator.get_momentum_vectors" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_momentum_vectors(cls, phases, transformations): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L386)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L386?message=Update%20Docs)]
</div>
Transforms the momenta so that they're aligned along the Gaussian axes
  - `:returns`: `_`
    >


<a id="Psience.DGB.Evaluators.DGBEvaluator.get_phase_vectors" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_phase_vectors(cls, momenta, transformations): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L406)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L406?message=Update%20Docs)]
</div>
Transforms the momenta so that they're aligned along the Gaussian axes
  - `:returns`: `_`
    >


<a id="Psience.DGB.Evaluators.DGBEvaluator.get_overlap_gaussians" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_overlap_gaussians(cls, centers, alphas, transformations, momenta, *, chunk_size=None, rows_cols=None, logger=None, parallelizer=None) -> 'OverlapGaussianData': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L426)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L426?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBEvaluator.poch" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
poch(cls, n, m): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L445)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L445?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBEvaluator.polyint_1D" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
polyint_1D(cls, centers, alphas, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L460)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L460?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBEvaluator.momentum_coeffient" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
momentum_coeffient(cls, k, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L471)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L471?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBEvaluator.momentum_integral" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
momentum_integral(cls, p, a, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L480)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L480?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBEvaluator.simple_poly_int" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
simple_poly_int(cls, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L490)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L490?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBEvaluator.tensor_expansion_integrate" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
tensor_expansion_integrate(cls, npts, derivs, overlap_data: 'OverlapGaussianData', expansion_type='multicenter', logger=None, reweight=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L493)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L493?message=Update%20Docs)]
</div>
provides an integral from a polynomial expansion with derivs as an expansion in tensors
  - `npts`: `Any`
    > 
  - `derivs`: `Any`
    > 
  - `centers`: `Any`
    > 
  - `alphas`: `Any`
    > 
  - `inds`: `Any`
    > 
  - `rot_data`: `Any`
    > 
  - `expansion_type`: `Any`
    > 
  - `logger`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DGB.Evaluators.DGBEvaluator.quad_weight_eval" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
quad_weight_eval(cls, function, d_chunk, w_chunk, ndim, centers, squa): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L648)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L648?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBEvaluator.quad_nd" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
quad_nd(cls, centers, alphas, function, flatten=False, degree=3, chunk_size=1000000, normalize=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L667)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L667?message=Update%20Docs)]
</div>
N-dimensional quadrature
  - `centers`: `Any`
    > 
  - `alphas`: `Any`
    > 
  - `function`: `Any`
    > 
  - `degree`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DGB.Evaluators.DGBEvaluator.rotated_gaussian_quadrature" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
rotated_gaussian_quadrature(cls, function, alphas, centers, rotations, inverse, momenta, normalize=True, degree=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L784)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L784?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Evaluators.DGBEvaluator.quad_integrate" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
quad_integrate(cls, function, overlap_data: 'OverlapGaussianData', degree=2, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L821)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L821?message=Update%20Docs)]
</div>
Integrate potential over all pairs of Gaussians at once
  - `degree`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DGB.Evaluators.DGBEvaluator.evaluate_overlap" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
evaluate_overlap(cls, overlap_data: 'OverlapGaussianData', logger=None, return_prefactor=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L926)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L926?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Evaluators/DGBEvaluator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Evaluators/DGBEvaluator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Evaluators/DGBEvaluator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Evaluators/DGBEvaluator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Evaluators.py#L334?message=Update%20Docs)   
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