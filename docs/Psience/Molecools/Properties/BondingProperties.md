## <a id="Psience.Molecools.Properties.BondingProperties">BondingProperties</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Properties.py#L755)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Properties.py#L755?message=Update%20Docs)]
</div>

The set of properties that depend only on bonding
and that kind of format







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Molecools.Properties.BondingProperties.get_prop_adjacency_matrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_adjacency_matrix(cls, atoms, bonds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L761)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L761?message=Update%20Docs)]
</div>
Returns the adjacency matrix for the molecule
  - `bonds`: `Iterable[int]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.BondingProperties.get_prop_connectivity" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_connectivity(cls, atoms, bonds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L779)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L779?message=Update%20Docs)]
</div>
Returns the adjacency matrix for the molecule
  - `bonds`: `Iterable[int]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.BondingProperties.get_prop_fragments" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_fragments(cls, atoms, bonds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L797)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L797?message=Update%20Docs)]
</div>
Returns the fragments for the molecule
  - `bonds`: `Iterable[int]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.BondingProperties.get_prop_zmat_ordering" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_zmat_ordering(cls, atoms, bonds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L818)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L818?message=Update%20Docs)]
</div>
Gets a guessed Z-matrix ordering for the molecule with connectivity defined by bonds based on the following:
1. Fragments are separated out
2. The atom with the highest degree of connectivity in each fragment is chosen as the fragment "label"
3. Fragments are ordered by connectivity of the label from high to low
4. Fragment labels reference each other with:
a) the second label on the x-axis
b) the 3rd in the xy-plane
c) all others relative to the first three
5. All other atoms are sorted first by fragment label, then by connection to the fragment label, and then by connectivity
6. Atoms reference each other based on the following:
a) if the atom has one bond:
i)   the atom it is bound to
ii)  the lowest-connectivity atom that one is bound to
iii) the second-lowest-connectivity atom OR the next fragment label
b) if the atom has two bonds:
i)   the highest-connectivity atom it is bound to
ii)  the lowest-connectivity atom it is bound to
iii) the lowest-connectivity atom (i) is bound to
c) if the atom has three bonds:
i)   the highest-connectivity atom it is bound to
ii)  the lowest-connectivity atom it is bound to
iii) the second-highest connectivity atom it is bound to
if any of these atoms do not exist, the fragment labels will be used in their place
  - `bonds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.BondingProperties.get_prop_guessed_bonds" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_guessed_bonds(cls, mol, tol=1.05, guess_type=True, covalent_radius_scaling=1.1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L857)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L857?message=Update%20Docs)]
</div>
Guesses the bonds for the molecule by finding the ones that are less than some percentage of a single bond for that
pair of elements
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/BondingProperties.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/BondingProperties.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/BondingProperties.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/BondingProperties.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Properties.py#L755?message=Update%20Docs)   
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