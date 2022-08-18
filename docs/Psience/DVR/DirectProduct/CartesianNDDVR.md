## <a id="Psience.DVR.DirectProduct.CartesianNDDVR">CartesianNDDVR</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L269)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L269?message=Update%20Docs)]
</div>

Provides an ND-DVR over different domains



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.DVR.DirectProduct.CartesianNDDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, domains, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L273)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L273?message=Update%20Docs)]
</div>

 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [energies_2D](#energies_2D)
- [energies_3D](#energies_3D)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
#### <a class="collapse-link" data-toggle="collapse" href="#test-setup">Setup</a> <a class="float-right" data-toggle="collapse" href="#test-setup"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="test-setup" markdown="1">

Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.
```python
from Peeves.TestUtils import *
from unittest import TestCase
from McUtils.Data import UnitsData, PotentialData
from McUtils.Zachary import Interpolator
import McUtils.Plots as plt
from Psience.DVR import *
from Psience.Molecools import Molecule
import numpy as np
```

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

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DirectProduct/CartesianNDDVR.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DirectProduct/CartesianNDDVR.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DirectProduct/CartesianNDDVR.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DirectProduct/CartesianNDDVR.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L269?message=Update%20Docs)