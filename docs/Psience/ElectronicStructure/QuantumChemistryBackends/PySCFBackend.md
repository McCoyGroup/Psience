## <a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend">PySCFBackend</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L168)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L168?message=Update%20Docs)]
</div>

PySCF energies, gradients and (for HF/DFT) Hessians, plus SCF-level charges, dipoles, and
static polarizabilities. Every requested property comes from one converged calculation.
  - `atoms`: `Any`
    > the atom symbols
  - `level_of_theory`: `Any`
    > ``'hf'``/``'rhf'``/``'uhf'``/``'rohf'``, ``'mp2'``, or a DFT functional name
    (``'b3lyp'``, ``'wb97x-d3bj'``, ...; a trailing ``-d3bj``/``-d3zero``/``-d4`` selects a dispersion correction,
    which needs ``pyscf-dispersion``); or a callable ``caller(molecule, **opts)`` returning a finished calculation
  - `basis_set`: `Any`
    > the basis set (``basis`` is accepted as an alias)
  - `charge`: `Any`
    > the total charge
  - `multiplicity`: `Any`
    > the spin multiplicity; PySCF's own ``spin`` (2S) is set to ``multiplicity - 1``.
    Open shells use unrestricted references (``'rohf'`` for restricted open-shell HF).
  - `quiet`: `Any`
    > silence PySCF's output
  - `molecule_options`: `Any`
    > extra options for `gto.Mole` (coordinates are always passed in Bohr)
  - `level_of_theory_options`: `Any`
    > extra options for the calculation (``functional_options``, ``restricted``,
    ``reference`` for MP2, and any SCF attribute such as ``conv_tol`` or ``max_cycle``)
  - `disp`: `Any`
    > a dispersion correction (``'d3bj'``, ``'d3zero'``, ``'d4'``) for DFT
  - `scf_options`: `Any`
    > SCF attributes to set on every calculation (``{'conv_tol': 1e-10, 'max_cycle': 200}``)
  - `grid_level`: `Any`
    > the DFT integration grid level
  - `density_fit`: `Any`
    > `True`, or an auxiliary basis name, to use density fitting
  - `threads`: `Any`
    > the number of OpenMP threads to use for each calculation
  - `max_memory`: `Any`
    > PySCF's memory limit in MB
  - `soscf_fallback`: `Any`
    > if the first-order SCF doesn't converge, continue with a second-order (Newton) SCF
  - `allow_unconverged`: `Any`
    > warn instead of raising when the SCF still hasn't converged
  - `reuse_guess`: `Any`
    > start each SCF from the previous geometry's density (faster along scans and
    optimizations; fresh guesses are used if that fails)
  - `analytic_hessian`: `Any`
    > treat analytic Hessians as available (HF/DFT only)
  - `charge_model`: `Any`
    > ``'mulliken'`` (default) or ``'meta_lowdin'`` population charges






<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
distance_units: str
property_units: dict
analytic_derivative_orders: dict
property_shapes: dict
```
<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, *, level_of_theory, basis_set=None, basis=None, caller=None, charge=None, multiplicity=None, quiet=True, molecule_options=None, level_of_theory_options=None, disp=None, scf_options=None, grid_level=None, density_fit=None, threads=None, max_memory=None, soscf_fallback=True, allow_unconverged=False, reuse_guess=False, analytic_hessian=False, charge_model='mulliken'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L218)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L218?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.handle_specialization" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_specialization(cls, tag): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L282)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L282?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.discover_dependencies" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
discover_dependencies(cls, *, disp=None, level_of_theory=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L289)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L289?message=Update%20Docs)]
</div>
PySCF, plus ``pyscf-dispersion`` if a dispersion correction is requested (via ``disp`` or a ``-d3bj``-style
suffix on the level of theory).


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.caller_dispatch" class="docs-object-method">&nbsp;</a> 
```python
@property
caller_dispatch(self) -> McUtils.Devutils.Options.OptionsMethodDispatch: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L307)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L307?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.get_callers" class="docs-object-method">&nbsp;</a> 
```python
get_callers(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L319)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L319?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.resolve_caller" class="docs-object-method">&nbsp;</a> 
```python
resolve_caller(self, caller, lot): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L329)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L329?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.get_molecule" class="docs-object-method">&nbsp;</a> 
```python
get_molecule(self, coords, charge=None, multiplicity=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L338)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L338?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.run_calculation" class="docs-object-method">&nbsp;</a> 
```python
run_calculation(self, coords, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L442)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L442?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.get_energy_from_calc" class="docs-object-method">&nbsp;</a> 
```python
get_energy_from_calc(self, calc) -> float: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L454)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L454?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.get_gradient_from_calc" class="docs-object-method">&nbsp;</a> 
```python
get_gradient_from_calc(self, calc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L462)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L462?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.get_hessian_from_calc" class="docs-object-method">&nbsp;</a> 
```python
get_hessian_from_calc(self, calc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L467)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L467?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.get_charges_from_calc" class="docs-object-method">&nbsp;</a> 
```python
get_charges_from_calc(self, calc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L483)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L483?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.get_dipole_from_calc" class="docs-object-method">&nbsp;</a> 
```python
get_dipole_from_calc(self, calc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L493)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L493?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.get_polarizability_from_calc" class="docs-object-method">&nbsp;</a> 
```python
get_polarizability_from_calc(self, calc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L497)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L497?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.PySCFBackend.evaluate_properties" class="docs-object-method">&nbsp;</a> 
```python
evaluate_properties(self, coords, orders, charge, multiplicity, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L513)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.py#L513?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/ElectronicStructure/QuantumChemistryBackends/PySCFBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L168?message=Update%20Docs)   
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