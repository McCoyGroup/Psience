# <a id="Psience.DVR.DVR.DVR">DVR</a>
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DVR.py#L126)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DVR.py#L126?message=Update%20Docs)]
</div>

```python
DVR(domain=None, divs=None, classes=None, potential_function=None, g=None, g_deriv=None, scf=False, potential_optimize=False, **base_opts): 
```
Constructs a DVR object
  - `domain`: `Any`
    > 
  - `divs`: `Any`
    > 
  - `classes`: `Any`
    > 
  - `potential_function`: `Any`
    > 
  - `g`: `Any`
    > 
  - `g_deriv`: `Any`
    > 
  - `base_opts`: `Any`
    > 
  - `:returns`: `_`
    > 



## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-5a01e2" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-5a01e2"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-5a01e2" markdown="1">
 - [1D](#1D)
- [energies_1D](#energies_1D)
- [energies_2D](#energies_2D)
- [energies_3D](#energies_3D)
- [RingDVR1D](#RingDVR1D)
- [RingDVR1DExplicitMass](#RingDVR1DExplicitMass)
- [RingDVR2DExplicitMass](#RingDVR2DExplicitMass)
- [RingDVR1DCosMass](#RingDVR1DCosMass)
- [Ring3D](#Ring3D)
- [Ring3DCosMass3D](#Ring3DCosMass3D)
- [Ring2DDifferentMass](#Ring2DDifferentMass)
- [MBPolDVR](#MBPolDVR)
- [MoleculeDVR](#MoleculeDVR)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-98dd31" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-98dd31"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-98dd31" markdown="1">
 
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
    def setupMBPolModel(self):
        ...
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

#### <a name="RingDVR1D">RingDVR1D</a>
```python
    def test_RingDVR1D(self):
        dvr_1D = RingDVR()
        npts = 5
        n = (5 - 1) // 2
        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=npts,
                         mass=1,
                         result='grid'
                         )
        self.assertTrue(np.allclose(
            res.grid,
            (2 * np.pi) * np.arange(1, npts + 1) / npts
        ))

        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=5,
                         mass=1,
                         result='kinetic_energy'
                         )
        self.assertTrue(np.allclose(res.kinetic_energy,
                                    [
                                        [1.0, -0.5854101966249685, 8.541019662496847e-2, 8.54101966249685e-2,
                                         -0.5854101966249681],
                                        [-0.5854101966249685, 1.0, -0.5854101966249685, 8.541019662496847e-2,
                                         8.54101966249685e-2],
                                        [8.541019662496847e-2, -0.5854101966249685, 1.0, -0.5854101966249685,
                                         8.541019662496847e-2],
                                        [8.54101966249685e-2, 8.541019662496847e-2, -0.5854101966249685, 1.0,
                                         -0.5854101966249685],
                                        [-0.5854101966249681, 8.54101966249685e-2, 8.541019662496847e-2,
                                         -0.5854101966249685, 1.0]
                                    ]
                                    ))

        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=5,
                         mass=1,
                         result='potential_energy'
                         )
        self.assertTrue(np.allclose(np.diag(res.potential_energy), np.sin(res.grid)))
```

#### <a name="RingDVR1DExplicitMass">RingDVR1DExplicitMass</a>
```python
    def test_RingDVR1DExplicitMass(self):

        dvr_1D = RingDVR()
        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         mass=1/(2*0.000197),
                         divs=251,
                         flavor='[0,2pi]'
                         )

        print(
            UnitsData.convert("Hartrees", "Wavenumbers")*(
                    res.wavefunctions[:5].energies[1:] -
                    res.wavefunctions[:5].energies[0]
            )
            )
```

#### <a name="RingDVR2DExplicitMass">RingDVR2DExplicitMass</a>
```python
    def test_RingDVR2DExplicitMass(self):

        dvr_2D = RingNDDVR((25, 25))
        res = dvr_2D.run(potential_function=self.cos_sin_pot,
                         domain=((0, 2 * np.pi),)*2,
                         mass=[1/(2*0.000197), 1/(2*0.000197)],
                         divs=(25, 25),
                         flavor='[0,2pi]',
                         diag_mode='dense'
                         )

        print(
            UnitsData.convert("Hartrees", "Wavenumbers") * (
                    res.wavefunctions[:5].energies[1:] -
                    res.wavefunctions[:5].energies[0]
            )
        )
        # self.assertTrue(np.allclose(
        #     res.wavefunctions[:5].energies.tolist(), [-0.536281, 0.341958, 0.854909, 2.05781, 2.08047],
        #     atol=.03  # different eigensolvers?
        # ))
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
```

#### <a name="RingDVR1DCosMass">RingDVR1DCosMass</a>
```python
    def test_RingDVR1DCosMass(self):
        dvr_1D = RingDVR()
        res = dvr_1D.run(potential_function=np.sin,
                         g=np.cos,
                         g_deriv=lambda g:-np.cos(g),
                         domain=(0, 2*np.pi),
                         divs=251
                         )
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
```

#### <a name="Ring3D">Ring3D</a>
```python
    def test_Ring3D(self):
        dvr_3D = RingNDDVR((15,) * 3)
        res = dvr_3D.run(mass=1, potential_function=self.cos3D)
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
```

#### <a name="Ring3DCosMass3D">Ring3DCosMass3D</a>
```python
    def test_Ring3DCosMass3D(self):
        dvr_3D = RingNDDVR((15,) * 3)

        res_basic = dvr_3D.run(potential_function=self.cos3D,
                               mass=1,
                               domain=((0, 2 * np.pi),) * 3,
                               divs=(15,) * 3,
                               flavor='[0,2pi]'
                               )

        g_el = lambda vals: np.full(len(vals), 1/2)
        gd_el = lambda vals: np.zeros(len(vals))

        res = dvr_3D.run(potential_function=self.cos3D,
                         domain=((0, 2 * np.pi),) * 3,
                         divs=(15,) * 3,
                         g=[
                             [g_el, 0, 0],
                             [0, g_el, 0],
                             [0, 0, g_el]
                         ],
                         g_deriv=[gd_el, gd_el, gd_el],
                         flavor='[0,2pi]',
                         num_wavefunctions=2
                         )

        self.assertTrue(np.allclose(res.wavefunctions.energies, res_basic.wavefunctions.energies))

        g_el = lambda vals: (2 + np.cos(vals[..., 0])) # plausible in size????
        gd_el = lambda vals: -np.cos(vals)/10
        res = dvr_3D.run(potential_function=self.cos3D,
                         domain=((0, 2*np.pi),) * 3,
                         divs=(15,) * 3,
                         g=[
                             [g_el, 0, 0],
                             [0, g_el, 0],
                             [0, 0, g_el]
                         ],
                         g_deriv=[gd_el, gd_el, gd_el],
                         flavor='[0,2pi]',
                         num_wavefunctions=2
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
```

#### <a name="Ring2DDifferentMass">Ring2DDifferentMass</a>
```python
    def test_Ring2DDifferentMass(self):

        dvr_2D = RingNDDVR((45, 45))


        g_tt = lambda vals: np.full(len(vals), 2); gd_tt = lambda vals: np.zeros(len(vals))
        g_HH = lambda vals: np.full(len(vals), 3); gd_HH = lambda vals: np.zeros(len(vals))

        # g_tt = lambda vals: (2 + np.cos(vals[..., 0])); gd_tt = lambda vals: -np.cos(vals[..., 0])
        # g_HH = lambda vals: (2 + np.cos(2*vals[..., 1])); gd_HH = lambda vals: -np.sin(vals[..., 1])

        zero_pot = lambda v: np.zeros(len(v))
        res = dvr_2D.run(potential_function=zero_pot,
                         g=[
                             [g_tt, 0],
                             [0, g_HH]
                         ],
                         g_deriv=[gd_tt, gd_HH],
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
        f1 = res.wavefunctions.energies
        ke1 = res.kinetic_energy

        # g_tH = lambda vals: np.cos(vals[..., 0])*np.cos(2*vals[..., 1]) # test
        # g_tH = lambda vals: np.zeros(len(vals)) # zeros
        g_tH = lambda vals: np.full(len(vals), 1) # constant
        res = dvr_2D.run(potential_function=zero_pot,
                         g=[
                             [g_tt, g_tH],
                             [g_tH, g_HH]
                         ],
                         g_deriv=[gd_tt, gd_HH],
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
        f2 = res.wavefunctions.energies
        ke2 = res.kinetic_energy

        # raise Exception(f1, f2)
        # raise Exception(np.round(f1), np.round(f2))

        self.assertLess(np.max(np.abs(f2-f1)), 5)

        # ke1 = ke1.toarray()
        # ke_coupling = ke2 - ke1
        #
        # import McUtils.Plots as plt
        #
        # plt.ArrayPlot(ke1, colorbar=True)
        # plt.ArrayPlot(ke2, colorbar=True)
        # plt.ArrayPlot(ke_coupling, colorbar=True)
        # ke_coupling_2 = np.zeros_like(ke_coupling)
        # wat = np.nonzero(ke1)
        # ke_coupling_2[wat] = ke_coupling[wat]
        #
        # raise Exception(np.max(ke_coupling_2), np.min(ke_coupling_2))
        # plt.ArrayPlot(ke_coupling_2, colorbar=True).show()

        # print(ke2 - ke1)

        dvr_3D = RingNDDVR((15, 15, 15))
        g_tH = lambda vals: (2 + np.cos(vals[..., 0])*np.cos(2*vals[..., 1]))
        res = dvr_3D.run(potential_function=self.cos3D,
                         g=[
                             [g_tt, g_tH,    0],
                             [g_tH, g_HH,    0],
                             [0,       0, g_HH]
                         ],
                         g_deriv=[gd_tt, gd_HH, gd_HH],
                         flavor='[0,2pi]',
                         diag_mode='dense'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
```

#### <a name="MBPolDVR">MBPolDVR</a>
```python
    def test_MBPolDVR(self):
        loader = ModuleLoader(TestManager.current_manager().test_data_dir)
        mbpol = loader.load("LegacyMBPol").potential
        mol = Molecule.from_file(TestManager.test_data("water_freq.fchk")) # we won't bother to reoptimize


        mol.zmatrix = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]

        disps_r = np.linspace(-.2, .5, 25) / UnitsData.bohr_to_angstroms
        # # disps_r = np.zeros_like(disps_r)
        # dr_coords0 = mol.get_displaced_coordinates(
        #     disps_r[:, np.newaxis],
        #     [0],
        #     internals='convert',
        #     strip_embedding=True
        # )
        # dr_coords1 = mol.get_displaced_coordinates(
        #     disps_r[:, np.newaxis],
        #     [1],
        #     internals='convert',
        #     strip_embedding=True
        # )

        # pot0 = mbpol(dr_coords0)[0]
        # pot1 = mbpol(dr_coords1)[0]
        # p1 = plt.Plot(disps_r, pot0 * 219475.6)
        # plt.Plot(disps_r, pot1 * 219475.6).show()

        disps_grid = np.array(np.meshgrid(disps_r, disps_r))
        dr_coords = mol.get_displaced_coordinates(
            np.moveaxis(disps_grid, 0, -1).reshape(-1, disps_grid.shape[0]),
            [0, 1],
            use_internals='convert',
            strip_embedding=True
        )

        # raise Exception(
        #     np.linalg.norm(dr_coords[:, 0] - dr_coords[:, 1], axis=1),
        #     np.linalg.norm(dr_coords[:, 0] - dr_coords[:, 2], axis=1)
        # )

        pot = mbpol(dr_coords)[0]
        plt.ContourPlot(*disps_grid, pot.reshape(disps_grid[0].shape) * 219475.6).show()
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DVR/DVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DVR/DVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DVR/DVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DVR/DVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DVR.py#L126?message=Update%20Docs)   
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