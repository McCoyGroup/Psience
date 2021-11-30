## <a id="Psience.VPT2.Runner.VPTRuntimeOptions">VPTRuntimeOptions</a>
Provides a helper to keep track of the options available
for configuring the way the code runs

### Properties and Methods
<a id="Psience.VPT2.Runner.VPTRuntimeOptions.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, operator_chunk_size=None, logger=None, verbose=None, checkpoint=None, results=None, parallelizer=None, memory_constrained=None, checkpoint_keys=None, use_cached_representations=None, use_cached_basis=None): 
```

- `operator_chunk_size`: `int`
    >the number of representation matrix elements to calculate at once
- `logger`: `Logger`
    >the `Logger` object to use when logging the status of the calculation
- `verbose`: `bool`
    >whether or not to be verbose in log output
- `checkpoint`: `str`
    >the checkpoint file or `Checkpointer` object to use
- `parallelizer`: `Parallelizer`
    >the `Parallelizer` object to use when parallelizing pieces of the calculation
- `memory_constrained`: `bool`
    >whether or not to attempt memory optimizations
- `checkpoint_keys`: `Iterable[str]`
    >the keys to write to the checkpoint file
- `use_cached_representations`: `bool`
    >whether or not to try to load representation matrices from the checkpoint
- `use_cached_basis`: `bool`
    >whether or not to try to load the bases to use from the checkpoint

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Runner/VPTRuntimeOptions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Runner/VPTRuntimeOptions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Runner/VPTRuntimeOptions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Runner/VPTRuntimeOptions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py?message=Update%20Docs)