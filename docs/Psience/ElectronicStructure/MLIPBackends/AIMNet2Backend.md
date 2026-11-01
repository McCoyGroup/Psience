## <a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend">AIMNet2Backend</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L39)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L39?message=Update%20Docs)]
</div>

AIMNet2 (and AIMNet2-NSE) energies, charges, and fluctuating-charge polarizabilities, with
derivatives by autograd. Borrows structure from `AIMNet2ASE` to call the model directly.







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
polarizability_modes: tuple
```
<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, model='aimnet2', charge=None, device=None, model_dir=None, multiplicity=None, quiet=True, polarizability_mode='coulomb', response_solver='constrained', allow_indefinite=False, stability_tolerance=1e-10): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L64)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L64?message=Update%20Docs)]
</div>

  - `model`: `Any`
    > an AIMNet2 model name (`'aimnet2'`, `'aimnet2nse'`, ...) or a loaded `AIMNet2Calculator`
  - `device`: `Any`
    > the torch device (default: CUDA if available)
  - `model_dir`: `Any`
    > a cache directory for downloaded weights (default `$MODEL_CACHE_DIR`)
  - `polarizability_mode`: `Any`
    > 'coulomb', 'net_coulomb', or experimental 'total_energy'
  - `response_solver`: `Any`
    > 'constrained' (null-space solve) or 'inverse' (full inverse)
  - `allow_indefinite`: `Any`
    > permit unstable curvature for diagnostics; singular kernels still fail
  - `stability_tolerance`: `Any`
    > relative threshold for near-zero constrained eigenvalues


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.handle_specialization" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_specialization(cls, tag): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L89)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L89?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.setup_model" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
setup_model(cls, model, device=None, model_dir=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L173)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L173?message=Update%20Docs)]
</div>
Load an `AIMNet2Calculator` (or return `model` if it already is one).


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.autodiff" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
autodiff(expr, coord, pad_dim, base_shape=None, create_graph=False, retain_graph=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L209)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L209?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.prep_eval" class="docs-object-method">&nbsp;</a> 
```python
prep_eval(self, coords, order, charge=None, multiplicity=None, keep_graph=None, forces=None, hessian=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L269)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L269?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.process_aimnet_derivs" class="docs-object-method">&nbsp;</a> 
```python
process_aimnet_derivs(self, base_shape, coords, data, root_key, order, extra_deriv_coord=None, deriv_key=None, property_shape=(), pad_coord=1, replace_nans=True, output_shape=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L408)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L408?message=Update%20Docs)]
</div>
Collect `data[root_key]` and its coordinate derivatives through `order` (and, with
`extra_deriv_coord=(var, n)`, mixed derivatives through order `n` in `var`) as numpy arrays.


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_energy" class="docs-object-method">&nbsp;</a> 
```python
evaluate_energy(self, coords, order, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L508)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L508?message=Update%20Docs)]
</div>
`[E, dE/dx, ...]` on a flat batch `(nstruct, natoms, 3)` (Hessians are taken one structure at a time)


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_charges" class="docs-object-method">&nbsp;</a> 
```python
evaluate_charges(self, coords, order, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L520)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L520?message=Update%20Docs)]
</div>
`[q, dq/dx, ...]` on a flat batch `(nstruct, natoms, 3)` (derivatives are taken one
structure at a time: the padded per-atom autograd doesn't separate batched molecules)


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_electrostatic_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
evaluate_electrostatic_potential_derivatives(self, coords, coord_order=0, electrostatic_order=2, charge=None, multiplicity=None, polarizability_mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L658)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L658?message=Update%20Docs)]
</div>
Independent-charge derivatives of the selected AIMNet2 energy.

Returns ``[[E, dE/dR], [E_q, dE_q/dR], [E_qq, dE_qq/dR]]`` truncated to
requested orders. Units are eV, elementary charge, and Angstroms. The
total-energy mode intervenes before the final neural pass; the explicit
Coulomb modes differentiate their native modules. No output keys or
Coulomb-module settings are changed. Native charges need not be stationary.


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_charge_hardness" class="docs-object-method">&nbsp;</a> 
```python
evaluate_charge_hardness(self, coords, coord_order=0, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L679)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L679?message=Update%20Docs)]
</div>
``[E_qq, dE_qq/dR]`` in native eV/e^2 and Angstrom units.


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_polarizability" class="docs-object-method">&nbsp;</a> 
```python
evaluate_polarizability(self, coords, coord_order=0, *, polarizability_mode=None, response_solver=None, allow_indefinite=None, stability_tolerance=None, return_diagnostics=False, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L683)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L683?message=Update%20Docs)]
</div>
Conserved-charge polarizability in atomic units, derivatives per Angstrom.

``coulomb`` uses the previous external-module control; ``net_coulomb``
includes embedded short-range corrections; ``total_energy`` includes the
downstream learned energy on a conditional charge surface. The latter is
experimental: native AIMNet2 charges are not variational minimizers.
Singular kernels always fail. Unstable kernels require explicit opt-in.
Diagnostics are available here without changing the standard property API.


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_dipole" class="docs-object-method">&nbsp;</a> 
```python
evaluate_dipole(self, coords, order=0, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L724)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L724?message=Update%20Docs)]
</div>
Native charge dipole in e Bohr, coordinate derivatives per Angstrom.

No charge re-equilibration is implied. Charged-system dipoles use the
supplied coordinate origin. This agrees with AIMNet2DipoleEvaluator.


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_field_response" class="docs-object-method">&nbsp;</a> 
```python
evaluate_field_response(self, coords, coord_order=0, field_order=2, charge=None, multiplicity=None, polarizability_mode=None, response_solver=None, allow_indefinite=None, stability_tolerance=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L740)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L740?message=Update%20Docs)]
</div>
``[native_dipole_expansion, polarizability_expansion][:field_order]``.

Property values are atomic units and nuclear derivative axes are in
Angstroms. Selecting a curvature mode changes alpha, not native charges.


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_properties" class="docs-object-method">&nbsp;</a> 
```python
evaluate_properties(self, coords, orders, charge, multiplicity, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L760)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L760?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L39?message=Update%20Docs)   
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