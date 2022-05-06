## <a id="Psience.Molecools.Vibrations.MolecularVibrations">MolecularVibrations</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L19)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L19?message=Update%20Docs)]
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, basis, freqs=None, init=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L21)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L21?message=Update%20Docs)]
</div>

Sets up a vibration for a Molecule object over the CoordinateSystem basis
- `molecule`: `AbstractMolecule`
    >No description...
- `init`: `None | CoordinateSet`
    >No description...
- `basis`: `MolecularNormalModes`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.molecule" class="docs-object-method">&nbsp;</a> 
```python
@property
molecule(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `CoordinateSet`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L78)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L78?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.displace" class="docs-object-method">&nbsp;</a> 
```python
displace(self, displacements=None, amt=0.1, n=1, which=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L81)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L81?message=Update%20Docs)]
</div>

Displaces along the vibrational mode specified by `which`
- `displacements`: `Any`
    >No description...
- `amt`: `Any`
    >No description...
- `n`: `Any`
    >No description...
- `which`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.visualize" class="docs-object-method">&nbsp;</a> 
```python
visualize(self, step_size=5, steps=(2, 2), which=0, anim_opts=None, mode='fast', **plot_args): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L113)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L113?message=Update%20Docs)]
</div>


- `step_size`: `Any`
    >No description...
- `steps`: `Any`
    >No description...
- `which`: `Any`
    >No description...
- `anim_opts`: `Any`
    >No description...
- `mode`: `Any`
    >No description...
- `plot_args`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.to_widget" class="docs-object-method">&nbsp;</a> 
```python
to_widget(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L196)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L196?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L206)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L206?message=Update%20Docs)]
</div>

Takes a slice of the modes
- `item`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.embed" class="docs-object-method">&nbsp;</a> 
```python
embed(self, frame): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L240)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L240?message=Update%20Docs)]
</div>


- `frame`: `MolecularTransformation`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.rescale" class="docs-object-method">&nbsp;</a> 
```python
rescale(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L255)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L255?message=Update%20Docs)]
</div>

Multiplies each mode by some scaling factor
- `phases`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L270)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L270?message=Update%20Docs)]
</div>

 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [VisualizeNormalModes](#VisualizeNormalModes)

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
from Peeves import BlockProfiler
from Psience.Molecools import Molecule, MolecularNormalModes
from Psience.Data import DipoleSurface
from McUtils.GaussianInterface import GaussianFChkReader, GaussianLogReader
from McUtils.Plots import *
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Data import UnitsData
import numpy as np
import McUtils.Numputils as nput
```

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

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Vibrations/MolecularVibrations.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Vibrations/MolecularVibrations.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Vibrations/MolecularVibrations.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Vibrations/MolecularVibrations.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L19?message=Update%20Docs)