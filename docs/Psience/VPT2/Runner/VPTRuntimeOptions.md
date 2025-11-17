## <a id="Psience.VPT2.Runner.VPTRuntimeOptions">VPTRuntimeOptions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L920)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L920?message=Update%20Docs)]
</div>

Provides a helper to keep track of the options available
for configuring the way the code runs







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Runner.VPTRuntimeOptions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, operator_chunk_size=None, matrix_element_threshold=None, nondeg_hamiltonian_precision=None, logger=None, verbose=None, checkpoint=None, results=None, parallelizer=None, memory_constrained=None, checkpoint_keys=None, use_cached_representations=None, use_cached_basis=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L939)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L939?message=Update%20Docs)]
</div>

  - `operator_chunk_size`: `int|None default:None`
    > the number of representation matrix elements to calculate in at one time
  - `matrix_element_threshold`: `float|None default:None`
    > the minimum size of matrix element to keep
  - `nondeg_hamiltonian_precision`: `int`
    > the precision with which to print out elements in the degenerate coupling Hamiltonians in the log file
  - `logger`: `str|Logger|bool|None default:None`
    > the `Logger` object to use when logging the status of the calculation (`True` means log normally)
  - `results`: `str|Checkpointer|None default:None`
    > the `Checkpointer` to write corrections out to
  - `parallelizer`: `Parallelizer|None default:None`
    > the `Parallelizer` to use for parallelizing the evaluation of matrix elements
  - `memory_constrained`: `bool|None`
    > whether or not to attempt memory optimizations (`None` means attempt for >20D problems)
  - `checkpoint`: `str|Checkpointer|None default:None`
    > the `Checkpointer` to write Hamiltonians and other bits out to
  - `checkpoint_keys`: `Iterable[str]|None`
    > which keys to save in the checkpoint
  - `use_cached_representations`: `bool`
    > whether other not to use Hamiltonian reps from the checkpoint
  - `use_cached_basis`: `bool`
    > whether other not to use bases from the checkpoint
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTRuntimeOptions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTRuntimeOptions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTRuntimeOptions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTRuntimeOptions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L920?message=Update%20Docs)   
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