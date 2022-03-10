## <a id="Psience.DVR.ColbertMiller.CartesianDVR">CartesianDVR</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/ColbertMiller.py#L17)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/ColbertMiller.py#L17?message=Update%20Docs)]
</div>

Provides the Colbert Miller DVR on the Cartesian [-inf, inf] range

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.DVR.ColbertMiller.CartesianDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/ColbertMiller.py#L22)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/ColbertMiller.py#L22?message=Update%20Docs)]
</div>

Provides the Colbert-Miller DVR grid for the [-inf, inf] range
- `domain`: `Any`
    >No description...
- `divs`: `Any`
    >No description...
- `flavor`: `Any`
    >No description...
- `kw`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.DVR.ColbertMiller.CartesianDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/ColbertMiller.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/ColbertMiller.py#L39?message=Update%20Docs)]
</div>

<a id="Psience.DVR.ColbertMiller.CartesianDVR.real_momentum" class="docs-object-method">&nbsp;</a> 
```python
real_momentum(self, grid=None, mass=None, hb=1, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/ColbertMiller.py#L64)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/ColbertMiller.py#L64?message=Update%20Docs)]
</div>

Provides the real part of the momentum for the [0, 2pi] range
- `grid`: `Any`
    >No description...
- `hb`: `Any`
    >No description...
- `kw`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [1D](#1D)
- [energies_1D](#energies_1D)
- [MoleculeDVR](#MoleculeDVR)

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

#### <a name="1D">1D</a>
```python
    def test_1D(self):
        dvr_1D = CartesianDVR(domain=(-5, 5), divs=250)
        pot = dvr_1D.run(potential_function=self.ho, result='potential_energy')
        self.assertIsInstance(pot.potential_energy, np.ndarray)
```
#### <a name="energies_1D">energies_1D</a>
```python
    def test_energies_1D(self):
        dvr_1D = CartesianDVR()
        res = dvr_1D.run(potential_function=self.ho, domain=(-5, 5), divs=250, mass=1)
        # print(e[:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions.energies, np.ndarray)
        self.assertTrue(np.allclose(res.wavefunctions.energies[:5].tolist(), [1/2, 3/2, 5/2, 7/2, 9/2]))
```
#### <a name="MoleculeDVR">MoleculeDVR</a>
```python
    def test_MoleculeDVR(self):

        scan_coords = Molecule.from_file(TestManager.test_data("water_HOH_scan.log"))
        scan_coords.zmatrix = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        g = scan_coords.g_matrix

        a_vals = scan_coords.bond_angle(1, 0, 2)

        g_func = Interpolator(a_vals, g[:, 2, 2])
        g_deriv = g_func.derivative(2)

        # Get potential
        pot = scan_coords.potential_surface.unshare(coordinates=((1, 0, 2),))
        min_pos = np.argmin(pot.base.interp_data[1])
        g_eq = g[min_pos][2, 2]

        carts = CartesianDVR(domain=(np.min(a_vals), np.max(a_vals)), divs=251,
                             mass=1/g_eq,
                             potential_function=pot,
                             nodeless_ground_state=True
                             )
        res_const = carts.run()

        carts = CartesianDVR(domain=(np.min(a_vals), np.max(a_vals)),
                             divs=251,
                             g=g_func,
                             g_deriv=g_deriv,
                             potential_function=pot,
                             nodeless_ground_state=True
                             )
        res = carts.run()

        print(
            (res.wavefunctions.energies[1] - res.wavefunctions.energies[0])*UnitsData.convert("Hartrees", "Wavenumbers"),
            (res_const.wavefunctions.energies[1] - res_const.wavefunctions.energies[0])*UnitsData.convert("Hartrees", "Wavenumbers")
        )

        grid = plt.GraphicsGrid(nrows=1, ncols=3,
                                subimage_size=(400, 400),
                                spacings=[70, 0],
                                padding=[[50, 0], [50, 50]],
                                figure_label='Water HOH DVR'
                                )

        res.plot_potential(figure=grid[0, 0], zero_shift=True, plot_units='wavenumbers'); grid[0, 0].plot_label = 'HOH Potential'
        res.wavefunctions[(0, 3, 7),].plot(figure=grid[0, 1]); grid[0, 1].plot_label = 'HOH Wavefunctions'
        res_const.wavefunctions[(0, 3, 7),].plot(figure=grid[0, 2]); grid[0, 2].plot_label = 'Constant G'
        # wf_ploot = res.wavefunctions[(0, 3, 7),]
        # wf_ploot.plot(figure=grid[0, 2]); grid[0, 2].plot_label ='HOH G-Matrix'

        grid.show()
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/ColbertMiller/CartesianDVR.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/ColbertMiller/CartesianDVR.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/ColbertMiller/CartesianDVR.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/ColbertMiller/CartesianDVR.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/ColbertMiller.py#L17?message=Update%20Docs)