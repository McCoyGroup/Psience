# <a id="Psience.Data">Psience.Data</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/__init__.py#L1)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/__init__.py#L1?message=Update%20Docs)]
</div>
    
Provides core Data-related types and data structures.
Intended to be high-level data as opposed to the lower-level stuff in `McUtils.Data`.
That means including stuff like dipole and potential energy surfaces that know how to compute their own properties.
We also have expressions for G-matrix elements from Frederick and Woywood to use with `sympy`.

### Members
<div class="container alert alert-secondary bg-light">
  <div class="row">
   <div class="col" markdown="1">
[DipoleSurface](Data/Surfaces/DipoleSurface.md)   
</div>
   <div class="col" markdown="1">
[PotentialSurface](Data/Surfaces/PotentialSurface.md)   
</div>
   <div class="col" markdown="1">
[KEData](Data/KEData/KEData.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[KEDataHandler](Data/KEData/KEDataHandler.md)   
</div>
   <div class="col" markdown="1">
[PotentialRegistryAPI](Data/PotentialRegistry/PotentialRegistryAPI.md)   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>





## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-272360" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-272360"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-272360" markdown="1">
 - [FChkFileDipoleSurface](#FChkFileDipoleSurface)
- [LogFileDipoleSurface](#LogFileDipoleSurface)
- [LogFilePotentialSurface](#LogFilePotentialSurface)
- [PotentialRegistry](#PotentialRegistry)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-c33c96" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-c33c96"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-c33c96" markdown="1">
 
Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.

All tests are wrapped in a test class
```python
class DataTests(TestCase):
    maxDiff = None
```

 </div>
</div>

#### <a name="FChkFileDipoleSurface">FChkFileDipoleSurface</a>
```python
    def test_FChkFileDipoleSurface(self):
        fchk = TestManager.test_data("HOD_freq.fchk")
        surf = DipoleSurface.from_fchk_file(fchk)
        surf_center = surf.surfs[0].base.data['center']
        self.assertIsInstance(surf_center, np.ndarray)
        self.assertTrue(
            np.allclose(surf(surf_center) - np.array([s.base.data['ref'] for s in surf.surfs]), 0.)
        )
        self.assertEquals(surf([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).shape, (1, 3))
        self.assertEquals(surf([
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        ]).shape, (2, 3))
```

#### <a name="LogFileDipoleSurface">LogFileDipoleSurface</a>
```python
    def test_LogFileDipoleSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: cartesian_to_zmatrix(
            x, ordering=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        ).coords[:, 0, 0]
        surf = DipoleSurface.from_log_file(log, conv)
        dips = surf(np.arange(.5, 2, .1))
        self.assertEquals(dips.shape, ((2-.5)/.1, 3))
```

#### <a name="LogFilePotentialSurface">LogFilePotentialSurface</a>
```python
    def test_LogFilePotentialSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: np.linalg.norm(x[:, 0] - x[:, 1], axis=1)
        surf = PotentialSurface.from_log_file(log, conv)
        pots = surf(np.arange(.5, 2, .1))
        self.assertEquals(pots.shape, ((2-.5)/.1,))
```

#### <a name="PotentialRegistry">PotentialRegistry</a>
```python
    def test_PotentialRegistry(self):

        from Psience.Molecools import Molecule

        h2co_mod = PotentialRegistryAPI().get_potential('H2COPot')
        def cart_pot(coords, order=None):
            return h2co_mod.Potential.get_pot(coords)
        def internal_pot(coords, order=None):
            # internals = np.moveaxis(
            #     np.array([rOC, rCH1, rCH2, aOCH1, aOCH2, dOCHH]),
            #     0, -1
            # )
            coords = coords[..., (0, 1, 3, 2, 4, 5)]
            vals = h2co_mod.InternalsPotential.get_pot(coords, threading_mode='serial')
            # if coords.ndim > 1:
            #     vv = vals.reshape(-1)
            #     cc = coords.reshape((-1, 6))
            # else:
            #     vv = [vals]
            #     cc = [coords]
            # for c, v in zip(cc, vv):
            #     print(c, "==>", v)
            return vals

        ochh_base = Molecule.from_file(
            TestManager.test_data('OCHH_freq.fchk'),
            energy_evaluator={
                'potential_function':cart_pot,
                'permutation':[2, 3, 1, 0],
                "distance_units": "Angstroms",
                "energy_units": "Wavenumbers"
            }
        )

        cart_eng = ochh_base.calculate_energy()

        #
        # opt_ochh = ochh_base.optimize(max_displacement=.1, max_iterations=50)
        # opt_plot = opt_ochh.plot(backend='x3d')
        # base_plot = ochh_base.plot(backend='x3d', highlight_atoms=[0, 1, 2, 3], figure=opt_plot)
        # base_plot.show()
        # return
        # print(ochh_base.calculate_energy(), opt_ochh.calculate_energy())
        # print(opt_ochh.coords - ochh_base.coords)

        ochh_base = Molecule.from_file(
            TestManager.test_data('OCHH_freq.fchk'),
            energy_evaluator={
                'potential_function': internal_pot,
                # 'permutation': [2, 3, 1, 0],
                "distance_units": "Angstroms",
                "energy_units": "Wavenumbers",
                "strip_embedding": True,
                # "flatten_internals": True
            },
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  1,  0, -1],
                [3,  1,  0,  2],
            ]
        )

        opt_ochh = ochh_base.optimize(
            # method='quasi-newton'
            method='conjugate-gradient'
            # method='gradient-descent'
            # , max_iterations=100
            , stencil=3
            # , logger=True
            # , max_displacement=.01
            , prevent_oscillations=3
            , restart_interval=15
            # , mesh_spacing=1e-2
        )
        # print("...")
        b1 = ochh_base.calculate_energy()
        b2 = opt_ochh.calculate_energy()
        print(b1, b2)
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Data.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Data.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Data.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Data.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/__init__.py#L1?message=Update%20Docs)   
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