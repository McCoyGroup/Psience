## <a id="Psience.DVR.ColbertMiller.CartesianDVR">CartesianDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller.py#L17?message=Update%20Docs)]
</div>

Provides the Colbert Miller DVR on the Cartesian [-inf, inf] range







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.ColbertMiller.CartesianDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller/CartesianDVR.py#L22)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller/CartesianDVR.py#L22?message=Update%20Docs)]
</div>
Provides the Colbert-Miller DVR grid for the [-inf, inf] range
  - `domain`: `Any`
    > 
  - `divs`: `Any`
    > 
  - `flavor`: `Any`
    > 
  - `kw`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.ColbertMiller.CartesianDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller/CartesianDVR.py#L39)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller/CartesianDVR.py#L39?message=Update%20Docs)]
</div>


<a id="Psience.DVR.ColbertMiller.CartesianDVR.real_momentum" class="docs-object-method">&nbsp;</a> 
```python
real_momentum(self, grid=None, mass=None, hb=1, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller/CartesianDVR.py#L64)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller/CartesianDVR.py#L64?message=Update%20Docs)]
</div>
Provides the real part of the momentum for the [0, 2pi] range
  - `grid`: `Any`
    > 
  - `hb`: `Any`
    > 
  - `kw`: `Any`
    > 
  - `:returns`: `_`
    >
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-0f273d" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-0f273d"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-0f273d" markdown="1">
 - [1D](#1D)
- [energies_1D](#energies_1D)
- [MoleculeDVR](#MoleculeDVR)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-f96716" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-f96716"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-f96716" markdown="1">
 
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/ColbertMiller/CartesianDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/ColbertMiller/CartesianDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/ColbertMiller/CartesianDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/ColbertMiller/CartesianDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller.py#L17?message=Update%20Docs)   
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