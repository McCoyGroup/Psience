## <a id="Psience.DVR.Wavefunctions.DVRWavefunctions">DVRWavefunctions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L125)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L125?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
wavefunction_class: DVRWavefunction
```
<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energies=None, wavefunctions=None, grid=None, results: Psience.DVR.BaseDVR.DVRResults = None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L128)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L128?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a collection of DVR wavefunctions sharing a common grid and the `DVRResults` object they were solved from.
  - `energies`: `np.ndarray | None`
    > the energies of each wavefunction in the collection
  - `wavefunctions`: `np.ndarray | None`
    > the matrix of wavefunction values at the grid points, one column per state
  - `grid`: `np.ndarray | None`
    > the DVR grid the wavefunctions are defined on
  - `results`: `DVRResults | None`
    > the `DVRResults` object (grid, kinetic energy, potential energy, etc.) the wavefunctions were solved from
  - `opts`: `dict`
    > extra options forwarded to the base `Wavefunctions.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L151)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L151?message=Update%20Docs)]
</div>
Provides just the state needed to reconstruct this
`DVRWavefunctions` object -- the energies, the (n_grid, n_states)
matrix of wavefunction values, the grid they're defined over, and
any explicit state-space `indices` -- paralleling
`PerturbationTheoryWavefunctions.to_state` /
`PerturbationTheoryCorrections.to_state` in `Psience.VPT2`.

Like `PerturbationTheoryCorrections.to_state` deliberately drops
the (potentially large, and not always cleanly serializable, since
it can carry a raw `potential_function` closure and a `Logger`)
`hamiltonians` payload, this drops `results` (the parent
`DVRResults`/DVR object) rather than trying to serialize it: the
wavefunction data itself is what's needed to reuse the
wavefunctions downstream (e.g. to build a
`ContractedDVRHarmonicRepresentation`).
  - `serializer`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L181)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L181?message=Update%20Docs)]
</div>
Reloads a `DVRWavefunctions` object from the state produced by
`to_state`.
  - `data`: `Any`
    > 
  - `serializer`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L202)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L202?message=Update%20Docs)]
</div>
**LLM Docstring**

Debug string representation showing the class name, the number of wavefunctions, and the DVR object they were solved from.
  - `:returns`: `str`
    > string of the form `ClassName(num=N, DVR=dvr)`


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L247)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L247?message=Update%20Docs)]
</div>
Plots the held wavefunctions
  - `figure`: `Any`
    > 
  - `graphics_class`: `Any`
    > 
  - `plot_style`: `Any`
    > 
  - `scaling`: `Any`
    > 
  - `shift`: `Any`
    > 
  - `opts`: `Any`
    > 
  - `:returns`: `Graphics`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other=None, multiplicative=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L281)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L281?message=Update%20Docs)]
</div>
Computes the expectation value of operator op over the wavefunction other and self
  - `other`: `DVRWavefunctions | np.ndarray`
    > 
  - `op`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.transform_operator" class="docs-object-method">&nbsp;</a> 
```python
transform_operator(self, M): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L310)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L310?message=Update%20Docs)]
</div>
**LLM Docstring**

Transform an operator matrix given in the DVR grid-point basis into the basis of these wavefunctions, by sandwiching it between the wavefunction coefficient matrix and its transpose.
  - `M`: `np.ndarray | scipy.sparse.spmatrix`
    > the operator matrix, in the DVR grid basis (dense or sparse)
  - `:returns`: `np.ndarray`
    > the operator matrix in the wavefunction basis


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.coordinate" class="docs-object-method">&nbsp;</a> 
```python
coordinate(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L325)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L325?message=Update%20Docs)]
</div>
**LLM Docstring**

The position-operator matrix in the wavefunction basis, computed as the expectation value of the DVR grid points themselves.
  - `:returns`: `np.ndarray`
    > the coordinate-operator matrix


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.momentum" class="docs-object-method">&nbsp;</a> 
```python
momentum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L335)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L335?message=Update%20Docs)]
</div>
**LLM Docstring**

The real part of the momentum-operator matrix in the wavefunction basis, computed from the underlying DVR's `real_momentum` operator and transformed into the wavefunction basis.
  - `:returns`: `np.ndarray`
    > the momentum-operator matrix


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.laplacian" class="docs-object-method">&nbsp;</a> 
```python
laplacian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L347)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L347?message=Update%20Docs)]
</div>
**LLM Docstring**

The Laplacian operator matrix in the wavefunction basis, derived from a fresh (unit-mass, zero-potential, uncoupled) kinetic-energy calculation on the underlying DVR and transformed into the wavefunction basis.
  - `:returns`: `np.ndarray`
    > the Laplacian operator matrix


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L360)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L360?message=Update%20Docs)]
</div>
**LLM Docstring**

The kinetic-energy operator matrix in the wavefunction basis, transformed from the stored `DVRResults.kinetic_energy` (in the grid basis).
  - `:returns`: `np.ndarray`
    > the kinetic-energy operator matrix


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.potential_energy" class="docs-object-method">&nbsp;</a> 
```python
potential_energy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L374)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunctions.py#L374?message=Update%20Docs)]
</div>
**LLM Docstring**

The potential-energy operator matrix in the wavefunction basis, transformed from the stored `DVRResults.kinetic_energy` (in the grid basis).
  - `:returns`: `np.ndarray`
    > the (mistakenly computed) kinetic-energy operator matrix
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/Wavefunctions/DVRWavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/Wavefunctions/DVRWavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/Wavefunctions/DVRWavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/Wavefunctions/DVRWavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L125?message=Update%20Docs)   
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