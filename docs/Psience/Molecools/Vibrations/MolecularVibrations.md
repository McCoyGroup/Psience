## <a id="Psience.Molecools.Vibrations.MolecularVibrations">MolecularVibrations</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations.py#L19?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Molecools.Vibrations.MolecularVibrations.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, basis, freqs=None, init=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L21)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L21?message=Update%20Docs)]
</div>
Sets up a vibration for a Molecule object over the CoordinateSystem basis
  - `molecule`: `AbstractMolecule`
    > 
  - `init`: `None | CoordinateSet`
    > 
  - `basis`: `MolecularNormalModes`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L43?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.molecule" class="docs-object-method">&nbsp;</a> 
```python
@property
molecule(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L49)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L49?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L57)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L57?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L64)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L64?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L78)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L78?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.displace" class="docs-object-method">&nbsp;</a> 
```python
displace(self, displacements=None, amt=0.1, n=1, which=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L81)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L81?message=Update%20Docs)]
</div>
Displaces along the vibrational mode specified by `which`
  - `displacements`: `Any`
    > 
  - `amt`: `Any`
    > 
  - `n`: `Any`
    > 
  - `which`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.visualize" class="docs-object-method">&nbsp;</a> 
```python
visualize(self, step_size=5, steps=(2, 2), which=0, anim_opts=None, mode='fast', **plot_args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L113)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L113?message=Update%20Docs)]
</div>

  - `step_size`: `Any`
    > 
  - `steps`: `Any`
    > 
  - `which`: `Any`
    > 
  - `anim_opts`: `Any`
    > 
  - `mode`: `Any`
    > 
  - `plot_args`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.to_widget" class="docs-object-method">&nbsp;</a> 
```python
to_widget(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L196)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L196?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularVibrations.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L206)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L206?message=Update%20Docs)]
</div>
Takes a slice of the modes
  - `item`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.embed" class="docs-object-method">&nbsp;</a> 
```python
embed(self, frame): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L240)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L240?message=Update%20Docs)]
</div>

  - `frame`: `MolecularTransformation`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.rescale" class="docs-object-method">&nbsp;</a> 
```python
rescale(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L255)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L255?message=Update%20Docs)]
</div>
Multiplies each mode by some scaling factor
  - `phases`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.rotate" class="docs-object-method">&nbsp;</a> 
```python
rotate(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L269)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L269?message=Update%20Docs)]
</div>
Multiplies each mode by some scaling factor
  - `phases`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularVibrations.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularVibrations.py#L284)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularVibrations.py#L284?message=Update%20Docs)]
</div>
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-3cc617" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-3cc617"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-3cc617" markdown="1">
 - [VisualizeNormalModes](#VisualizeNormalModes)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-9ee9b8" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-9ee9b8"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-9ee9b8" markdown="1">
 
Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.

All tests are wrapped in a test class
```python
class MolecoolsTests(TestCase):
    def setUp(self):
        self.test_log_water = TestManager.test_data("water_OH_scan.log")
        self.test_log_freq = TestManager.test_data("water_freq.log")
        self.test_HOD = TestManager.test_data("HOD_freq.fchk")
        self.test_fchk = TestManager.test_data("water_freq.fchk")
        self.test_log_h2 = TestManager.test_data("outer_H2_scan_new.log")
```

 </div>
</div>

#### <a name="VisualizeNormalModes">VisualizeNormalModes</a>
```python
    def test_VisualizeNormalModes(self):

        from Psience.Molecools.Vibrations import MolecularVibrations, MolecularNormalModes
        from McUtils.Plots import GraphicsGrid, Graphics3D

        m = Molecule.from_file(self.test_fchk, bonds = [[0, 1, 1], [0, 2, 1]])

        with GaussianFChkReader(self.test_fchk) as reader:
            parse = reader.parse(("VibrationalModes", "VibrationalData"))
        modes = parse["VibrationalModes"].T

        test_freqs = parse["VibrationalData"]["Frequencies"]

        nms = m.normal_modes
        realvibs = MolecularVibrations(m, basis=MolecularNormalModes(m, modes, freqs=test_freqs))

        realvibs.visualize(mode='jupyter') # get no bugs

        plot_vibrations = False
        if plot_vibrations:
            nmodes = 1
            mode_start = 0
            g = GraphicsGrid(nrows=2, ncols=nmodes,
                             graphics_class=Graphics3D,
                             plot_range = [[-2, 2], [-2, 2], [-2, 2]],
                             fig_kw = dict(figsize = (17, 5)),
                             tighten = True
                             )

            for i in range(nmodes):
                nms.visualize(step_size=.1, figure = g[0, i], which=mode_start + i,
                              anim_opts= dict(interval = 10)
                              )

            for i in range(nmodes):
                realvibs.visualize(step_size=.1, figure = g[1, i], which= mode_start+i,
                                   anim_opts= dict(interval = 10)
                                   )

            g.show()

        self.assertEquals(
            tuple(np.round(UnitsData.convert("Hartrees", "Wavenumbers")*nms.modes.freqs, 4)),
            tuple(np.round(test_freqs, 4))
        )
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Vibrations/MolecularVibrations.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Vibrations/MolecularVibrations.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Vibrations/MolecularVibrations.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Vibrations/MolecularVibrations.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations.py#L19?message=Update%20Docs)   
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