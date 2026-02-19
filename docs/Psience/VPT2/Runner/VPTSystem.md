## <a id="Psience.VPT2.Runner.VPTSystem">VPTSystem</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L43?message=Update%20Docs)]
</div>

Provides a little helper for setting up the input
system for a VPT job







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Runner.VPTSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, internals=None, dummy_atoms=None, modes=None, local_modes=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, potential_derivatives=None, potential_function=None, order=2, dipole_derivatives=None, eckart_embed=False, copy_mol=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L76)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L76?message=Update%20Docs)]
</div>

  - `mol`: `str | list | Molecule`
    > the molecule or system specification to use (doesn't really even need to be a molecule)
  - `internals`: `list | dict`
    > the Z-matrix for the internal coordinates optionally with a specification of a `conversion` and `inverse`
To supply a conversion function, provide a `dict` like so
```python
{
    'zmatrix': [[atom1, bond1, angle1, dihed1], [atom2, bond2, angle2, dihed2], ...] or None,
    'conversion': 'a function to convert from Z-matrix coordinates to desired coordinates',
    'inverse': 'the inverse conversion'
}
```
  - `modes`: `MolecularVibrations|dict`
    > the normal modes to use if not already supplied by the Molecule
  - `potential_derivatives`: `Iterable[np.ndarray]`
    > the derivatives of the potential to use for expansions
  - `dipole_derivatives`: `Iterable[np.ndarray]`
    > the set of dipole derivatives to use for expansions


<a id="Psience.VPT2.Runner.VPTSystem.prep_local_modes" class="docs-object-method">&nbsp;</a> 
```python
prep_local_modes(self, dRdX, dXdR=None, sort_freqs=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTSystem.py#L224)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTSystem.py#L224?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTSystem.nmodes" class="docs-object-method">&nbsp;</a> 
```python
@property
nmodes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTSystem.py#L254)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTSystem.py#L254?message=Update%20Docs)]
</div>
Provides the number of modes in the system
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTSystem.get_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_potential_derivatives(self, potential_function, order=2, **fd_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTSystem.py#L275)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTSystem.py#L275?message=Update%20Docs)]
</div>
Computes potential derivatives for the given function through finite difference
  - `potential_function`: `Any`
    > 
  - `order`: `Any`
    > 
  - `fd_opts`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTSystem.from_harmonic_scan" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_harmonic_scan(cls, scan_array): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L296)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L296?message=Update%20Docs)]
</div>
 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Details-03c80a" markdown="1"> Details</a> <a class="float-right" data-toggle="collapse" href="#Details-03c80a"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Details-03c80a" markdown="1">
 When using functions of internal (Z-matrix/polyspherical) coordinates, a sample form of the conversion function is
```python
def conv(r, t, f, **kwargs):
    '''
    Takes the bond lengths (`r`), angles `(t)` and dihedrals `(f)`,
    and returns new coordinates that are functions of these coordinates
    '''
    ... # convert the coordinates
    return np.array([r, t, f])
```
and then the inverse function will take the output of `conv` and return the original Z-matrix/polyspherical coordinates.
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L43?message=Update%20Docs)   
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