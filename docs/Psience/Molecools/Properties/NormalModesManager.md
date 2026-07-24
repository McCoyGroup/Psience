## <a id="Psience.Molecools.Properties.NormalModesManager">NormalModesManager</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L2751)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2751?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
recalc_normal_mode_tolerance: float
```
<a id="Psience.Molecools.Properties.NormalModesManager.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, normal_modes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L2752)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2752?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up the normal-modes manager, converting `normal_modes` into a `MolecularVibrations` object from a dict spec, rebinding it to `mol` if it already carries mode data, or leaving it as `None`.
  - `mol`: `AbstractMolecule`
    > the molecule this manager is attached to
  - `normal_modes`: `object | dict | None`
    > the normal modes, given as a `MolecularVibrations`-like object, a dict with `'matrix'`/`'freqs'`/(optionally) `'inverse'`/`'origin'` keys, or `None`
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.Properties.NormalModesManager.from_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_data(cls, mol, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2787)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2787?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `NormalModesManager` from raw normal-mode data of various shapes: an existing `MolecularVibrations` (or `None`), a dict of matrix/frequency data, or generic mode data understood by `NormalModes.prep_modes`.
  - `mol`: `AbstractMolecule`
    > the molecule the manager will be attached to
  - `data`: `MolecularVibrations | dict | object | None`
    > the raw normal-mode data
  - `:returns`: `NormalModesManager`
    > the constructed manager


<a id="Psience.Molecools.Properties.NormalModesManager.set_molecule" class="docs-object-method">&nbsp;</a> 
```python
set_molecule(self, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2823)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2823?message=Update%20Docs)]
</div>
**LLM Docstring**

Rebind this manager (and its stored modes, if any) to a different molecule.
  - `mol`: `AbstractMolecule`
    > the new molecule to associate
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.Properties.NormalModesManager.get_modes" class="docs-object-method">&nbsp;</a> 
```python
get_modes(self, quiet=False, allow_compute=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2837)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2837?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the normal modes, computing them via `get_normal_modes` if not already cached.
  - `quiet`: `bool`
    > if `True`, suppresses errors when the modes can't be found/computed
  - `allow_compute`: `bool`
    > whether force constants may be computed from an energy evaluator if not otherwise available
  - `:returns`: `MolecularVibrations`
    > the normal modes


<a id="Psience.Molecools.Properties.NormalModesManager.modes" class="docs-object-method">&nbsp;</a> 
```python
@property
modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2853)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2853?message=Update%20Docs)]
</div>

  - `:returns`: `MolecularVibrations`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.construct_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
construct_normal_modes(self, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2882)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2882?message=Update%20Docs)]
</div>
**LLM Docstring**

Coerce a variety of mode representations (a dict of matrix/kwargs, a raw `np.ndarray`, an object exposing `modes_by_coords`/`coords_by_modes`, or an existing `MolecularNormalModes`) into a `MolecularVibrations` object bound to this manager's molecule.
  - `modes`: `dict | np.ndarray | MolecularNormalModes | object`
    > the mode data to coerce
  - `:returns`: `MolecularVibrations`
    > the constructed (or passed-through) `MolecularVibrations` object


<a id="Psience.Molecools.Properties.NormalModesManager.load" class="docs-object-method">&nbsp;</a> 
```python
load(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2914)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2914?message=Update%20Docs)]
</div>
**LLM Docstring**

Return the currently cached normal modes without attempting to compute them.
  - `:returns`: `MolecularVibrations | None`
    > the cached modes, or `None` if none have been set


<a id="Psience.Molecools.Properties.NormalModesManager.update" class="docs-object-method">&nbsp;</a> 
```python
update(self, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2924)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2924?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.load_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
load_normal_modes(self, file=None, mode=None, rephase=True, recalculate=False, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3112)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3112?message=Update%20Docs)]
</div>
Loads potential derivatives from a file (or from `source_file` if set)
  - `file`: `Any`
    > 
  - `rephase`: `bool`
    > whether to rephase FChk normal modes or not
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_normal_modes(self, quiet=False, compute_force_constants=True, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3203)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3203?message=Update%20Docs)]
</div>
Loads normal modes from file or calculates
from force constants
  - `kwargs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.get_force_constants" class="docs-object-method">&nbsp;</a> 
```python
get_force_constants(self, compute_force_constants=True, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3231)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3231?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the potential-energy Hessian (force constants), computing it from the molecule's energy evaluator via a numerical order-2 calculation if not already stored on the molecule.
  - `compute_force_constants`: `bool`
    > whether to allow computing the force constants via `self.mol.calculate_energy` if not already available
  - `quiet`: `bool`
    > if `True`, returns `None` instead of raising when force constants aren't available and can't be computed
  - `:returns`: `np.ndarray | None`
    > the force-constant (Hessian) matrix


<a id="Psience.Molecools.Properties.NormalModesManager.get_dipole_derivative_based_rephasing" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_dipole_derivative_based_rephasing(cls, modes, analytic_dipoles, numerical_dipoles, strict=True, allow_swaps=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3258)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3258?message=Update%20Docs)]
</div>
**LLM Docstring**

Derive the sign (and, in limited cases, swap) corrections needed to align a set of normal modes' numerically-computed dipole derivatives with their analytically-computed counterparts, by comparing the inner-product matrix between the two (each normalized) and flagging modes whose overlap magnitude falls below a threshold as needing special handling.
  - `modes`: `MolecularNormalModes`
    > the normal modes whose phase convention is being checked/corrected
  - `analytic_dipoles`: `tuple | None`
    > the analytic dipole derivative expansion (only the first-derivative term is used)
  - `numerical_dipoles`: `tuple | None`
    > the numerical dipole derivative expansion (only the first-derivative term is used)
  - `strict`: `bool`
    > if `True`, raise on shape mismatches between the numerical and analytic derivatives rather than returning `None`
  - `allow_swaps`: `bool`
    > if `True`, attempt to resolve ambiguous (nearly-degenerate) modes by finding candidate index swaps instead of only sign corrections (currently raises with the candidate mapping rather than resolving it automatically)
  - `:returns`: `np.ndarray | None`
    > the rephasing matrix to apply to the modes, or `None` if the necessary dipole data isn't available


<a id="Psience.Molecools.Properties.NormalModesManager.get_fchk_normal_mode_rephasing" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_fchk_normal_mode_rephasing(cls, mol, modes, use_dipoles=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3359)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3359?message=Update%20Docs)]
</div>
Returns the necessary rephasing to make the numerical dipole derivatives
agree with the analytic dipole derivatives as pulled from a Gaussian FChk file
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.get_partial_cubic_based_rephasing" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_partial_cubic_based_rephasing(cls, modes, partial_cubics: numpy.ndarray, degenerate_freq_threshold=1e-06, equivalent_threshold=0.01, nonzero_threshold=4.5e-06, ignore_single_zeros=False, frequency_scale=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3442)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3442?message=Update%20Docs)]
</div>
**LLM Docstring**

Derive a rephasing (sign/permutation) matrix for a set of nearly-degenerate normal modes using consistency of the partial cubic force constants: transforms `partial_cubics` into the mode basis (optionally frequency-scaled), groups modes into degenerate blocks by frequency, resolves pairwise sign/swap ambiguities within each block by comparing triples of cubic force constants that should agree up to the rephasing, and returns the resulting rephasing transformation.
  - `modes`: `MolecularNormalModes`
    > the normal modes to rephase
  - `partial_cubics`: `np.ndarray`
    > the partial cubic force-constant tensor used to disambiguate the phase/ordering of degenerate modes
  - `degenerate_freq_threshold`: `float`
    > maximum frequency difference for two modes to be treated as degenerate
  - `equivalent_threshold`: `float`
    > relative-difference threshold for treating two cubic force-constant magnitudes as equal
  - `nonzero_threshold`: `float`
    > absolute magnitude below which a cubic force constant is treated as zero
  - `ignore_single_zeros`: `bool`
    > accepted parameter but not referenced in the current method body
  - `frequency_scale`: `bool`
    > whether to frequency-scale the cubic tensor before comparing terms
  - `:returns`: `np.ndarray`
    > the rephasing (sign/permutation) matrix to apply to the modes


<a id="Psience.Molecools.Properties.NormalModesManager.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, transf): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3948)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3948?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply a coordinate transformation to a copy of this normal-modes manager by embedding its modes through the transformation.
  - `transf`: `object`
    > the transformation to apply to the modes
  - `:returns`: `NormalModesManager`
    > a new `NormalModesManager` with the transformed modes


<a id="Psience.Molecools.Properties.NormalModesManager.insert_atoms" class="docs-object-method">&nbsp;</a> 
```python
insert_atoms(self, atoms, coords, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3968)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3968?message=Update%20Docs)]
</div>
Handles the insertion of new atoms into the structure
  - `atoms`: `tuple[str]`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `where`: `tuple[int]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.delete_atoms" class="docs-object-method">&nbsp;</a> 
```python
delete_atoms(self, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3991)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3991?message=Update%20Docs)]
</div>
Handles the deletion from the structure
  - `atoms`: `tuple[str]`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `where`: `tuple[int]`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/NormalModesManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/NormalModesManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/NormalModesManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/NormalModesManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2751?message=Update%20Docs)   
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