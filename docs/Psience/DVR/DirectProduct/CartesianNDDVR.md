## <a id="Psience.DVR.DirectProduct.CartesianNDDVR">CartesianNDDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct.py#L269)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct.py#L269?message=Update%20Docs)]
</div>

Provides an ND-DVR over different domains







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.DirectProduct.CartesianNDDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, domains, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct/CartesianNDDVR.py#L273)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct/CartesianNDDVR.py#L273?message=Update%20Docs)]
</div>
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-540bdf" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-540bdf"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-540bdf" markdown="1">
 - [energies_2D](#energies_2D)
- [energies_3D](#energies_3D)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-63be66" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-63be66"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-63be66" markdown="1">
 
Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.

All tests are wrapped in a test class
```python
class DVRTests(TestCase):
    def ho(self, grid, k=1):
        return k/2*np.power(grid, 2)
    def ho_2D(self, grid, k1=1, k2=1):
        return k1/2*np.power(grid[:, 0], 2) + k2/2*np.power(grid[:, 1], 2)
    def ho_3D(self, grid, k1=1, k2=1, k3=1):
        return k1/2*np.power(grid[:, 0], 2) + k2/2*np.power(grid[:, 1], 2) + k3/2*np.power(grid[:, 2], 2)
    def cos2D(self, grid):
        return np.cos(grid[..., 0]) * np.cos(grid[..., 1])
    def cos3D(self, grid):
        return np.cos(grid[..., 0]) * np.cos(grid[..., 1]) * np.cos(grid[..., 2])
    def cos_sin_pot(self, grid):
        return UnitsData.convert("Wavenumbers", "Hartrees")* 2500 / 8 * ((2 + np.cos(grid[..., :, 0])) * (2 + np.sin(grid[..., :, 1])) - 1)
```

 </div>
</div>

#### <a name="energies_2D">energies_2D</a>
```python
    def test_energies_2D(self):
        dvr_2D = CartesianNDDVR(((-5, 5, 25), (-5, 5, 25)))
        res = dvr_2D.run(potential_function=self.ho_2D, mass=1)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
```

#### <a name="energies_3D">energies_3D</a>
```python
    def test_energies_3D(self):
        dvr_3D = CartesianNDDVR(((-5, 5, 25), (-5, 5, 25), (-5, 5, 25)))
        res = dvr_3D.run(potential_function=self.ho_3D, mass=1)
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
```

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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DirectProduct/CartesianNDDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DirectProduct/CartesianNDDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DirectProduct/CartesianNDDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DirectProduct/CartesianNDDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct.py#L269?message=Update%20Docs)   
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