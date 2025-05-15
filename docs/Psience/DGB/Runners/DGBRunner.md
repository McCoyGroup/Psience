## <a id="Psience.DGB.Runners.DGBRunner">DGBRunner</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Runners.py#L21)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Runners.py#L21?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
default_num_plot_wfns: int
```
<a id="Psience.DGB.Runners.DGBRunner.prep_interpolation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_interpolation(cls, nms, coords, potential_function, symmetrizations=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L26)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L26?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.construct_from_mol_simulation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct_from_mol_simulation(cls, sim, mol, *, potential_function=None, dipole_function=None, use_cartesians=False, use_momenta=False, quadrature_degree=3, expansion_degree=2, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, skip_initial_configurations=True, alphas='virial', modes='normal', transformations='diag', pairwise_potential_functions=None, logger=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L47)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L47?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.construct_from_model_simulation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct_from_model_simulation(cls, sim, model, mol=None, *, use_cartesians=False, use_momenta=False, quadrature_degree=3, expansion_degree=2, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, skip_initial_configurations=True, modes='normal', transformations='diag', **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L163)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L163?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.construct_from_model" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct_from_model(cls, model, trajectories=10, *, sim=None, propagation_time=10, timestep=50, use_cartesians=False, use_momenta=False, pairwise_potential_functions=None, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, track_velocities=True, **aimd_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L246)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L246?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.from_mol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_mol(cls, mol, sim=None, *, potential_function=None, dipole_function=None, trajectories=10, propagation_time=10, timestep=50, use_cartesians=False, use_momenta=False, pairwise_potential_functions=None, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, trajectory_seed=None, total_energy=None, total_energy_scaling=None, sampled_modes=None, initial_energies=None, initial_displacements=None, initial_mode_directions=None, displaced_coords=None, track_velocities=True, **aimd_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L285)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L285?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.run_simple" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
run_simple(cls, system_spec, sim=None, plot_wavefunctions=True, plot_spectrum=True, trajectories=10, propagation_time=10, timestep=50, use_cartesians=False, use_momenta=False, pairwise_potential_functions=None, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, trajectory_seed=None, total_energy=None, total_energy_scaling=None, sampled_modes=None, initial_energies=None, initial_mode_directions=None, initial_displacements=None, displaced_coords=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L408)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L408?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.plot_dgb_potential" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
plot_dgb_potential(cls, dgb, mol, potential, coordinate_sel=None, domain=None, domain_padding=1, potential_cutoff=17000, potential_units='Wavenumbers', potential_min=0, plot_cartesians=None, plot_atoms=True, cmap=None, modes_nearest=False, plot_points=100, levels=24, **plot_styles): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L476)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L476?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.plot_gaussians" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
plot_gaussians(cls, dgb, mol, *, domain=None, domain_padding=1, cmap='RdBu', plot_dir=None, plot_name='gaussian_{i}.pdf', **plot_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L547)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L547?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.plot_wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
plot_wavefunctions(cls, wfns, dgb, mol, which=True, coordinate_sel=None, cartesians=None, plot_dir=None, plot_name='wfn_{i}.pdf', plot_label='{e} cm-1', plot_potential=True, potential_units='Wavenumbers', plot_atoms=None, plot_centers=True, potential_styles=None, scaling=None, **plot_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L579)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L579?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.plot_potential_from_spec" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
plot_potential_from_spec(cls, dgb, mol, spec, plot_centers=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L699)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L699?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.prep_plot_wavefunctions_spec" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_plot_wavefunctions_spec(cls, dgb, spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L725)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L725?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.run_dgb" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
run_dgb(cls, dgb: Psience.DGB.DGB.DGB, mol, plot_centers=True, plot_wavefunctions=True, plot_spectrum=False, pot_cmap='viridis', wfn_cmap='RdBu', wfn_points=100, wfn_contours=12, plot_dir=None, plot_potential=True, pot_points=100, domain=None, domain_padding=1, wavefunction_scaling=None, potential_cutoff=15000, potential_units='Wavenumbers', mode=None, nodeless_ground_state=None, min_singular_value=None, subspace_size=None, plot_similarity=False, similarity_cutoff=None, similarity_chunk_size=None, similar_det_cutoff=None, **plot_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L742)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L742?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.getMorseParameters" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
getMorseParameters(cls, w=None, wx=None, m1=None, m2=None, re=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L872)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L872?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.setupMorseFunction" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
setupMorseFunction(cls, model, i, j, w=None, wx=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L893)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L893?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Runners.DGBRunner.plot_interpolation_error" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
plot_interpolation_error(cls, dgb, pot): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L920)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L920?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Runners/DGBRunner.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Runners/DGBRunner.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Runners/DGBRunner.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Runners/DGBRunner.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Runners.py#L21?message=Update%20Docs)   
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