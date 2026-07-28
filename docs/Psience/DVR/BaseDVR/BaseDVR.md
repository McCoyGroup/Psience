## <a id="Psience.DVR.BaseDVR.BaseDVR">BaseDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR.py#L12)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L12?message=Update%20Docs)]
</div>

Provides the abstract interface for creating a
convenient runnable DVR that can be cleanly subclassed to provide
extensions







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.BaseDVR.BaseDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, domain=None, divs=None, potential_function=None, logger=None, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L19?message=Update%20Docs)]
</div>

  - `base_opts`: `Any`
    > base opts to use when running


<a id="Psience.DVR.BaseDVR.BaseDVR.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L47)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L47?message=Update%20Docs)]
</div>
**LLM Docstring**

Debug string representation showing the class name, domain, number of grid points, and potential function.
  - `:returns`: `str`
    > string of the form `ClassName(domain, pts=divs, pot=potential_function)`


<a id="Psience.DVR.BaseDVR.BaseDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L84)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L84?message=Update%20Docs)]
</div>
**LLM Docstring**

Abstract hook for building this DVR's 1D grid over the given domain/division count. Concrete DVR subclasses must implement this.
  - `domain`: `tuple | None`
    > the coordinate domain to build the grid over
  - `divs`: `int | None`
    > the number of grid points
  - `kwargs`: `dict`
    > extra representation-specific options
  - `:returns`: `np.ndarray`
    > never returns on the base class


<a id="Psience.DVR.BaseDVR.BaseDVR.grid" class="docs-object-method">&nbsp;</a> 
```python
grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L102)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L102?message=Update%20Docs)]
</div>
**LLM Docstring**

Build (or retrieve) this DVR's grid, falling back to the stored `domain`/`divs` if not given explicitly, via `get_grid`.
  - `domain`: `tuple | None`
    > the coordinate domain to build the grid over; defaults to `self.domain`
  - `divs`: `int | None`
    > the number of grid points; defaults to `self.divs`
  - `kwargs`: `dict`
    > extra options forwarded to `get_grid`
  - `:returns`: `np.ndarray`
    > the DVR grid


<a id="Psience.DVR.BaseDVR.BaseDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L130)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L130?message=Update%20Docs)]
</div>
**LLM Docstring**

Abstract hook for building this DVR's 1D kinetic-energy operator matrix on the given grid. Concrete DVR subclasses must implement this.
  - `grid`: `np.ndarray | None`
    > the DVR grid to build the operator on
  - `mass`: `float | None`
    > the particle mass
  - `hb`: `float`
    > the value of hbar to use
  - `kwargs`: `dict`
    > extra representation-specific options
  - `:returns`: `np.ndarray`
    > never returns on the base class


<a id="Psience.DVR.BaseDVR.BaseDVR.handle_kinetic_coupling" class="docs-object-method">&nbsp;</a> 
```python
handle_kinetic_coupling(self, grid, ke_1D, g, g_deriv, hb=1, logger=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L150)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L150?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply a (possibly coordinate-dependent) kinetic-coupling correction to a 1D kinetic-energy matrix: multiplies each off-diagonal element by the averaged `g`-value of its two grid points and adds a diagonal correction from `g_deriv`, following the standard variable-mass DVR kinetic-energy formula.
  - `grid`: `np.ndarray`
    > the DVR grid the kinetic-energy matrix is defined on
  - `ke_1D`: `np.ndarray`
    > the base (unit-`g`) 1D kinetic-energy matrix to correct
  - `g`: `callable | float | np.ndarray | None`
    > the kinetic-coupling function/value; if `None`, `ke_1D` is returned unchanged
  - `g_deriv`: `callable | float | np.ndarray | None`
    > the second-derivative correction term for `g`; required if `g` is given
  - `hb`: `float`
    > the value of hbar to use
  - `logger`: `Logger | None`
    > logger for diagnostics
  - `kwargs`: `dict`
    > extra options, unused
  - `:returns`: `np.ndarray`
    > the (possibly `g`-corrected) kinetic-energy matrix


<a id="Psience.DVR.BaseDVR.BaseDVR.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L206)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L206?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the full kinetic-energy operator, computing the base 1D operator (via `get_kinetic_energy`) and then applying any kinetic-coupling correction (via `handle_kinetic_coupling`); when a kinetic-coupling function `g` is supplied, the mass is fixed to `1` since `g` already encodes the effective mass.
  - `grid`: `np.ndarray | None`
    > the DVR grid to build the operator on; defaults to `self.grid()`
  - `mass`: `float | None`
    > the particle mass; required unless `g` is given
  - `hb`: `float`
    > the value of hbar to use
  - `g`: `callable | float | np.ndarray | None`
    > the kinetic-coupling function/value
  - `g_deriv`: `callable | float | np.ndarray | None`
    > the second-derivative correction term for `g`
  - `kwargs`: `dict`
    > extra options forwarded to `get_kinetic_energy`/`handle_kinetic_coupling`
  - `:returns`: `np.ndarray`
    > the kinetic-energy operator matrix


<a id="Psience.DVR.BaseDVR.BaseDVR.real_momentum" class="docs-object-method">&nbsp;</a> 
```python
real_momentum(self, grid=None, mass=None, hb=1, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L243)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L243?message=Update%20Docs)]
</div>
**LLM Docstring**

Abstract hook for the real part of the momentum-operator matrix on the given grid. Not implemented on the base class; concrete DVR subclasses that support it must override this.
  - `grid`: `np.ndarray | None`
    > the DVR grid to build the operator on
  - `mass`: `float | None`
    > the particle mass
  - `hb`: `float`
    > the value of hbar to use
  - `kwargs`: `dict`
    > extra representation-specific options
  - `:returns`: `np.ndarray`
    > never returns on the base class


<a id="Psience.DVR.BaseDVR.BaseDVR.potential_energy" class="docs-object-method">&nbsp;</a> 
```python
potential_energy(self, grid=None, potential_function=None, potential_values=None, potential_grid=None, logger=None, **pars): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L263)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L263?message=Update%20Docs)]
</div>
Calculates the potential energy at the grid points based
on dispatching on the input form of the potential
  - `grid`: `Any`
    > the grid of points built earlier in the DVR
  - `potential_function`: `Any`
    > a function to evaluate the potential energy at the points
  - `potential_values`: `Any`
    > the values of the potential at the DVR points
  - `potential_grid`: `Any`
    > a grid of points and values to be interpolated
  - `pars`: `Any`
    > ignored keyword arguments
  - `:returns`: `_`
    >


<a id="Psience.DVR.BaseDVR.BaseDVR.hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
hamiltonian(self, kinetic_energy=None, potential_energy=None, potential_threshold=None, **pars): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L368)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L368?message=Update%20Docs)]
</div>
Calculates the total Hamiltonian from the kinetic and potential matrices
  - `kinetic_energy`: `Any`
    > 
  - `potential_energy`: `np.ndarray | sp.spmatrix`
    > 
  - `potential_threshold`: `Any`
    > 
  - `pars`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.BaseDVR.BaseDVR.wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
wavefunctions(self, hamiltonian=None, num_wfns=25, nodeless_ground_state=False, diag_mode=None, logger=None, **pars): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L400)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L400?message=Update%20Docs)]
</div>
Calculates the wavefunctions for the given Hamiltonian.
Doesn't support any kind of pruning based on potential values although that might be a good feature
to support explicitly in the future
  - `hamiltonian`: `Any`
    > 
  - `num_wfns`: `Any`
    > 
  - `nodeless_ground_state`: `Any`
    > 
  - `diag_mode`: `Any`
    > 
  - `pars`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.BaseDVR.BaseDVR.run" class="docs-object-method">&nbsp;</a> 
```python
run(self, result='wavefunctions', logger=None, grid=None, potential_energy=None, kinetic_energy=None, hamiltonian=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/BaseDVR.py#L449)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/BaseDVR.py#L449?message=Update%20Docs)]
</div>

  - `:returns`: `DVRResults`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/BaseDVR/BaseDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/BaseDVR/BaseDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/BaseDVR/BaseDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/BaseDVR/BaseDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L12?message=Update%20Docs)   
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