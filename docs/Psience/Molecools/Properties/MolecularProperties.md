## <a id="Psience.Molecools.Properties.MolecularProperties">MolecularProperties</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L938)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L938?message=Update%20Docs)]
</div>

An object whose sole purpose in life is to get molecular properties
A property should be implemented in two parts:
1) a classmethod called get_prop_<prop name> that takes raw inputs and uses them to compute a property
2) a classmethod called <prop name> that extracts the property from a passed molecule

All properties should be appropriately vectorized to work on a single configuration or a set of configurations







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Molecools.Properties.MolecularProperties.mass_weighted_coords" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
mass_weighted_coords(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L948)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L948?message=Update%20Docs)]
</div>
Computes the moments of inertia
  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
g_matrix(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L962)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L962?message=Update%20Docs)]
</div>

  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.center_of_mass" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
center_of_mass(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L973)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L973?message=Update%20Docs)]
</div>
Computes the moments of inertia
  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.inertia_tensor" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
inertia_tensor(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L986)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L986?message=Update%20Docs)]
</div>
Computes the inertia tensors for the stored geometries
  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.moments_of_inertia" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
moments_of_inertia(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L999)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L999?message=Update%20Docs)]
</div>
Computes the moments of inertia
  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.principle_axis_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
principle_axis_data(cls, mol, sel=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1012)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1012?message=Update%20Docs)]
</div>
Generates the center of masses and inertial axes
  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.principle_axis_transformation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
principle_axis_transformation(cls, mol, sel=None, inverse=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1029)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1029?message=Update%20Docs)]
</div>
Generates the principle axis transformation for a Molecule
  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.eckart_embedding_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
eckart_embedding_data(cls, mol, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1041)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1041?message=Update%20Docs)]
</div>

  - `mol`: `AbstractMolecule`
    > 
  - `coords`: `Any`
    > 
  - `sel`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.eckart_transformation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
eckart_transformation(cls, mol, ref_mol, sel=None, inverse=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1063)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1063?message=Update%20Docs)]
</div>

  - `ref_mol`: `AbstractMolecule`
    > reference geometry
  - `mol`: `AbstractMolecule`
    > molecules to get Eckart embeddings for
  - `sel`: `Any`
    > coordinate selection to use when doing the Eckart stuff
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.eckart_embedded_coords" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
eckart_embedded_coords(cls, mol, coords, sel=None, in_paf=False, reset_com=True, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1090)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1090?message=Update%20Docs)]
</div>

  - `mol`: `AbstractMolecule`
    > 
  - `coords`: `Any`
    > 
  - `sel`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.coriolis_constants" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
coriolis_constants(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1115)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1115?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.MolecularProperties.translation_rotation_eigenvectors" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
translation_rotation_eigenvectors(cls, mol, sel=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1123)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1123?message=Update%20Docs)]
</div>

  - `mol`: `AbstractMolecule`
    > molecules to get eigenvectors for
  - `sel`: `Any`
    > coordinate selection to use when doing the rotation/translation calculations
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.fragments" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
fragments(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1138)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1138?message=Update%20Docs)]
</div>

  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.fragment_indices" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
fragment_indices(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1174)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1174?message=Update%20Docs)]
</div>

  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.edge_graph" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
edge_graph(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1190)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1190?message=Update%20Docs)]
</div>

  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.guessed_bonds" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
guessed_bonds(cls, mol, tol=1.05, guess_type=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1207)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1207?message=Update%20Docs)]
</div>
Guesses the bonds for the molecule by finding the ones that are less than some percentage of a single bond for that
pair of elements
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.get_prop_chemical_formula" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_chemical_formula(cls, atoms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1219)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1219?message=Update%20Docs)]
</div>

  - `atoms`: `Tuple[str]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.MolecularProperties.chemical_formula" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
chemical_formula(cls, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1229)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1229?message=Update%20Docs)]
</div>

  - `mol`: `AbstractMolecule`
    > 
  - `:returns`: `_`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/MolecularProperties.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/MolecularProperties.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/MolecularProperties.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/MolecularProperties.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L938?message=Update%20Docs)   
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