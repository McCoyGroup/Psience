## <a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAGEvaluationPlan">PTTensorCoeffProductDAGEvaluationPlan</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L5569)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L5569?message=Update%20Docs)]
</div>

Streams canonical tensor terms from a DAG without expanding it globally.

Small shared subgraphs are replayed from a bounded, batch-local cache.  A
node that emits more than ``max_cached_terms`` is deliberately not retained;
this prevents a high-order expression from replacing global materialization
with an equally unbounded evaluation cache.

Developer note
--------------
This is an intentionally conservative bridge to the existing evaluator: it
streams symbolic terms in bounded chunks, but each ``mul_along`` still uses
the legacy symbolic product machinery.  Consequently it bounds peak storage
without removing the dominant high-order Cartesian-product work.

A genuinely faster evaluator should compile this DAG to numerical blocks
whose axes describe open/fixed mode indices.  Leaf blocks should substitute
coefficient tensors and evaluate their ``PolyPath`` values once; ``add`` and
``scale`` then operate directly on arrays, while ``permute``, ``shift``, and
``free_up_indices`` transform block metadata/views.  Most importantly,
``mul_along`` should become an indexed relational join/tensor contraction of
its child blocks.  That join must preserve equality/distinctness constraints
between mode indices, energy-denominator metadata, and final square-root
factors; multiplying already-reduced child scalars is not equivalent.

Keep unsupported nodes on this streaming implementation until their direct
kernels have numerical parity tests.  The materialized evaluator must also
remain available as the reference/timing implementation.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAGEvaluationPlan.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, expression, cache=None, max_cached_terms=256): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L5599)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L5599?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAGEvaluationPlan.iter_terms" class="docs-object-method">&nbsp;</a> 
```python
iter_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAGEvaluationPlan.py#L5725)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAGEvaluationPlan.py#L5725?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAGEvaluationPlan.stats" class="docs-object-method">&nbsp;</a> 
```python
stats(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAGEvaluationPlan.py#L5730)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAGEvaluationPlan.py#L5730?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/PTTensorCoeffProductDAGEvaluationPlan.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/PTTensorCoeffProductDAGEvaluationPlan.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/PTTensorCoeffProductDAGEvaluationPlan.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/PTTensorCoeffProductDAGEvaluationPlan.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L5569?message=Update%20Docs)   
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