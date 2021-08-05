## <a id="Psience.Wavefun.Wavefunctions.Wavefunction">Wavefunction</a>
Represents a single wavefunction object

### Properties and Methods
<a id="Psience.Wavefun.Wavefunctions.Wavefunction.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, energy, data, parent=None, index=None, **opts): 
```

<a id="Psience.Wavefun.Wavefunctions.Wavefunction.plot" class="docs-object-method">&nbsp;</a>
```python
plot(self, figure=None, index=None, **opts): 
```
Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)
- `figure`: `Graphics | Graphics3D`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Wavefun.Wavefunctions.Wavefunction.expectation" class="docs-object-method">&nbsp;</a>
```python
expectation(self, op, other): 
```
Computes the expectation value of operator op over the wavefunction other and self
- `other`: `Wavefunction`
    >No description...
- `op`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Wavefun.Wavefunctions.Wavefunction.probability_density" class="docs-object-method">&nbsp;</a>
```python
probability_density(self): 
```
Computes the probability density of the current wavefunction
- `:returns`: `_`
    >No description...

### Examples


