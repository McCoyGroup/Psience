## <a id="Psience.DGB.Gaussians.DGBGaussians">DGBGaussians</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians.py#L53)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians.py#L53?message=Update%20Docs)]
</div>

A class to set up the actual N-dimensional Gaussians used in a DGB







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
gs_optimization_overlap_cutoff: float
default_energy_cutoff: float
bad_alpha_limit: float
bad_scaling_limit: float
```
<a id="Psience.DGB.Gaussians.DGBGaussians.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coords, alphas, transformations=None, *, momenta=None, poly_coeffs=None, kinetic_options=None, logger=None, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians.py#L58)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians.py#L58?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.overlap_data" class="docs-object-method">&nbsp;</a> 
```python
@property
overlap_data(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L104)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L104?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.get_S" class="docs-object-method">&nbsp;</a> 
```python
get_S(self, return_prefactor=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L116)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L116?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.get_T" class="docs-object-method">&nbsp;</a> 
```python
get_T(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L121)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L121?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.optimize" class="docs-object-method">&nbsp;</a> 
```python
optimize(self, optimizer_options, potential_function=None, logger=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L127)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L127?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.take_gaussian_selection" class="docs-object-method">&nbsp;</a> 
```python
take_gaussian_selection(self, full_good_pos): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L154)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L154?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.construct" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct(cls, coords, alphas, *, potential_expansion=None, potential_function=None, transformations=None, masses=None, atoms=None, modes=None, kinetic_options=None, internals=None, coordinate_selection=None, cartesians=None, gmat_function=None, momenta=None, poly_coeffs=None, logger=None, pairwise_potential_functions=None, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L557)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L557?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_normal_modes(cls, coords, potential_function, masses=None, atoms=None, internals=None, gmat_function=None, reference_structure=None, stationary_point_norm=0.01, project_transrot=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L719)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L719?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.get_reaction_path_transformations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_reaction_path_transformations(cls, coords, potential_function, gmat_function, stationary_point_norm=0.0001, sort_alphas=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L787)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L787?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.get_hessian_diagonalizing_transformations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_hessian_diagonalizing_transformations(cls, coords, potential_function, gmat_function, *, masses=None, project_transrot=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L855)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L855?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.dispatch_get_alphas" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
dispatch_get_alphas(self, alphas, centers, **extra_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L914)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L914?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.get_mass_alphas" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_mass_alphas(cls, centers, *, masses, scaling=10, use_mean=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L964)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L964?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.get_min_distance_alphas" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_min_distance_alphas(cls, masses, centers, scaling=0.25, use_mean=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L974)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L974?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.get_virial_alphas" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_virial_alphas(cls, coords, *, potential_function, gmat_function, transformations, scaling=0.5): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L994)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L994?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.canonicalize_poly_coeffs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
canonicalize_poly_coeffs(cls, coeffs, alphas): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1079)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1079?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.transformations" class="docs-object-method">&nbsp;</a> 
```python
@property
transformations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L1091)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L1091?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.canonicalize_transforms" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
canonicalize_transforms(self, coords, tfs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1097)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1097?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.prefactor" class="docs-object-method">&nbsp;</a> 
```python
@property
prefactor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L1136)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L1136?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.S" class="docs-object-method">&nbsp;</a> 
```python
@property
S(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L1141)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L1141?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.T" class="docs-object-method">&nbsp;</a> 
```python
@property
T(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L1149)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L1149?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.marginalize_out" class="docs-object-method">&nbsp;</a> 
```python
marginalize_out(self, indices, *, bad_alpha_limit=None, bad_scaling_limit=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L1160)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L1160?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.as_cartesians" class="docs-object-method">&nbsp;</a> 
```python
as_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L1218)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L1218?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Gaussians.DGBGaussians.plot_centers" class="docs-object-method">&nbsp;</a> 
```python
plot_centers(self, figure=None, xyz_sel=None, **plot_styles): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Gaussians/DGBGaussians.py#L1246)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians/DGBGaussians.py#L1246?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Gaussians/DGBGaussians.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Gaussians/DGBGaussians.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Gaussians/DGBGaussians.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Gaussians/DGBGaussians.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Gaussians.py#L53?message=Update%20Docs)   
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