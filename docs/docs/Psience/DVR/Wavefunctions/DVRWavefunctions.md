## <a id="Psience.DVR.Wavefunctions.DVRWavefunctions">DVRWavefunctions</a>


### Properties and Methods
<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, energies=None, wavefunctions=None, wavefunction_class=<class 'Psience.DVR.Wavefunctions.DVRWavefunction'>, results: Psience.DVR.BaseDVR.DVRResults = None, **opts): 
```

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__len__" class="docs-object-method">&nbsp;</a>
```python
__len__(self): 
```

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__iter__" class="docs-object-method">&nbsp;</a>
```python
__iter__(self): 
```

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__getitem__" class="docs-object-method">&nbsp;</a>
```python
__getitem__(self, item): 
```
Provides a single `DVRWavefunction` or slice of `DVRWavefunctions`
- `item`: `Any`
    >No description...
- `:returns`: `DVRWavefunction | DVRWavefunctions`
    >No description...

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.plot" class="docs-object-method">&nbsp;</a>
```python
plot(self, figure=None, graphics_class=None, plot_style=None, scaling=1, shift=0, **opts): 
```
Plots the held wavefunctions
- `figure`: `Any`
    >No description...
- `graphics_class`: `Any`
    >No description...
- `plot_style`: `Any`
    >No description...
- `scaling`: `Any`
    >No description...
- `shift`: `Any`
    >No description...
- `opts`: `Any`
    >No description...
- `:returns`: `Graphics`
    >No description...

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.expectation" class="docs-object-method">&nbsp;</a>
```python
expectation(self, op, other): 
```
Computes the expectation value of operator op over the wavefunction other and self
- `other`: `DVRWavefunctions | np.ndarray`
    >No description...
- `op`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.probability_density" class="docs-object-method">&nbsp;</a>
```python
probability_density(self): 
```
Computes the probability density of the set of wavefunctions
- `:returns`: `_`
    >No description...

### Examples


