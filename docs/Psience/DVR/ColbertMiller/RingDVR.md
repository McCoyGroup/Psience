## <a id="Psience.DVR.ColbertMiller.RingDVR">RingDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller.py#L92)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller.py#L92?message=Update%20Docs)]
</div>

Provides a DVR for working on the (0, 2Pi) range with periodicity from Colbert and Miller







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.ColbertMiller.RingDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, domain=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller/RingDVR.py#L97)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller/RingDVR.py#L97?message=Update%20Docs)]
</div>


<a id="Psience.DVR.ColbertMiller.RingDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller/RingDVR.py#L102)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller/RingDVR.py#L102?message=Update%20Docs)]
</div>
Provides the Colbert-Miller 1D grid for the [0, 2Pi] range
  - `domain`: `Any`
    > 
  - `divs`: `Any`
    > 
  - `kw`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.ColbertMiller.RingDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=1, hb=1, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller/RingDVR.py#L120)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller/RingDVR.py#L120?message=Update%20Docs)]
</div>
Colbert-Miller kinetic energy for the [0, 2pi] range
  - `grid`: `Any`
    > 
  - `mass`: `Any`
    > 
  - `hb`: `Any`
    > 
  - `kw`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.ColbertMiller.RingDVR.real_momentum" class="docs-object-method">&nbsp;</a> 
```python
real_momentum(self, grid=None, hb=1, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/ColbertMiller/RingDVR.py#L154)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller/RingDVR.py#L154?message=Update%20Docs)]
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
## <a class="collapse-link" data-toggle="collapse" href="#Tests-d65b73" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-d65b73"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-d65b73" markdown="1">
 - [RingDVR1D](#RingDVR1D)
- [RingDVR1DExplicitMass](#RingDVR1DExplicitMass)
- [RingDVR1DCosMass](#RingDVR1DCosMass)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-d48d86" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-d48d86"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-d48d86" markdown="1">
 
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/ColbertMiller/RingDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/ColbertMiller/RingDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/ColbertMiller/RingDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/ColbertMiller/RingDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/ColbertMiller.py#L92?message=Update%20Docs)   
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