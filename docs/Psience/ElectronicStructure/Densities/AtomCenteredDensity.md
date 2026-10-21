## <a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity">AtomCenteredDensity</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities.py#L15?message=Update%20Docs)]
</div>

A density written as a sum of atom-centered Gaussian-smeared multipoles,

    rho(r) = sum_i sum_lm p_ilm phi_lm(r - R_i; sigma_i),

in the GTO basis of `graph_longrange`/MACE-POLAR, in which `p_ilm` are the exact multipole
moments of each site (so the `l = 0` coefficients are condensed atomic populations).
Moments use the real Racah-normalized layout `m = -l..l` per `l` (`[q, mu_y, mu_z, mu_x, ...]`).

Everything is stored in atomic units: centers and widths in Bohr, moments in `e Bohr^l`
(or the corresponding units of whatever is being represented, e.g. unpaired electrons for a
spin density). Atoms are the last batch axis.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, multipoles, centers, smearing=None, atoms=None, label=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities.py#L30)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities.py#L30?message=Update%20Docs)]
</div>

  - `multipoles`: `np.ndarray`
    > real moments per site, `(..., natoms, (L+1)^2)` (or `(..., natoms)` for populations)
  - `centers`: `np.ndarray`
    > site positions in Bohr, `(..., natoms, 3)`
  - `smearing`: `float | np.ndarray | None`
    > Gaussian width(s) in Bohr (`None` for point sites: populations only)
  - `atoms`: `list[str] | None`
    > atom symbols
  - `label`: `str | None`
    > a name for the density


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.from_populations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_populations(cls, populations, centers, smearing=None, atoms=None, label=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L60)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L60?message=Update%20Docs)]
</div>

  - `populations`: `np.ndarray`
    > condensed (monopole) values per site, `(..., natoms)`
  - `:returns`: `AtomCenteredDensity`
    >


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.modify" class="docs-object-method">&nbsp;</a> 
```python
modify(self, multipoles=None, label=None, smearing=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L70)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L70?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L79)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L79?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L101)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L101?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.__sub__" class="docs-object-method">&nbsp;</a> 
```python
__sub__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L108)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L108?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.__neg__" class="docs-object-method">&nbsp;</a> 
```python
__neg__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L115)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L115?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, scalar): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L118)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L118?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, scalar): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L118)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L118?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.__truediv__" class="docs-object-method">&nbsp;</a> 
```python
__truediv__(self, scalar): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L123)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L123?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.populations" class="docs-object-method">&nbsp;</a> 
```python
@property
populations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L130)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L130?message=Update%20Docs)]
</div>
the condensed (monopole) values per site, `(..., natoms)`


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.total" class="docs-object-method">&nbsp;</a> 
```python
@property
total(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L135?message=Update%20Docs)]
</div>
the integral of the density


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.dipoles" class="docs-object-method">&nbsp;</a> 
```python
@property
dipoles(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L140)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L140?message=Update%20Docs)]
</div>
site dipoles (Cartesian, `e Bohr`), `(..., natoms, 3)`


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.participation_ratio" class="docs-object-method">&nbsp;</a> 
```python
participation_ratio(self, heavy_atoms_only=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L154)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L154?message=Update%20Docs)]
</div>
The participation ratio `(sum_i |p_i|)^2 / sum_i p_i^2` of the condensed values: the
effective number of sites the density is spread over (1 for a single site, N for
an even spread over N sites). Negative populations count by magnitude.
  - `heavy_atoms_only`: `bool`
    > restrict to non-hydrogen sites
  - `:returns`: `np.ndarray`
    >


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.localization" class="docs-object-method">&nbsp;</a> 
```python
localization(self, heavy_atoms_only=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L167)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L167?message=Update%20Docs)]
</div>
the inverse participation ratio, between `1/N` (delocalized) and 1 (one site)


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.condensed" class="docs-object-method">&nbsp;</a> 
```python
condensed(self, atoms=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L171)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L171?message=Update%20Docs)]
</div>
the condensed values as a list of `(index, atom, value)`


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.expansion" class="docs-object-method">&nbsp;</a> 
```python
@property
expansion(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L183)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L183?message=Update%20Docs)]
</div>
the underlying distributed `MultipoleExpansion` (atomic units)


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.density" class="docs-object-method">&nbsp;</a> 
```python
density(self, points): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L190)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L190?message=Update%20Docs)]
</div>
The density at `points` (Bohr), per Bohr^3. Requires smeared sites.
  - `points`: `np.ndarray`
    > `(N, 3)` or `(..., N, 3)`
  - `:returns`: `np.ndarray`
    >


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.potential" class="docs-object-method">&nbsp;</a> 
```python
potential(self, points, deriv_order=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L202)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L202?message=Update%20Docs)]
</div>
The Coulomb potential of the density at `points` (Bohr), in atomic units, with optional
Cartesian derivatives (as `MultipoleExpansion.evaluate`).


<a id="Psience.ElectronicStructure.Densities.AtomCenteredDensity.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L211)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities/AtomCenteredDensity.py#L211?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/Densities/AtomCenteredDensity.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/Densities/AtomCenteredDensity.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/Densities/AtomCenteredDensity.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/Densities/AtomCenteredDensity.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Densities.py#L15?message=Update%20Docs)   
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