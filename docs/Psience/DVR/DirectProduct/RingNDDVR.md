## <a id="Psience.DVR.DirectProduct.RingNDDVR">RingNDDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct.py#L284)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct.py#L284?message=Update%20Docs)]
</div>

Provides an ND-DVR for products of periodic (0, 2Pi) ranges







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.DirectProduct.RingNDDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, divs, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct/RingNDDVR.py#L289)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct/RingNDDVR.py#L289?message=Update%20Docs)]
</div>
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-9d5bbd" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-9d5bbd"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-9d5bbd" markdown="1">
 - [RingDVR2DExplicitMass](#RingDVR2DExplicitMass)
- [Ring3D](#Ring3D)
- [Ring3DCosMass3D](#Ring3DCosMass3D)
- [Ring2DDifferentMass](#Ring2DDifferentMass)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-e55a69" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-e55a69"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-e55a69" markdown="1">
 
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DirectProduct/RingNDDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DirectProduct/RingNDDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DirectProduct/RingNDDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DirectProduct/RingNDDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct.py#L284?message=Update%20Docs)   
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