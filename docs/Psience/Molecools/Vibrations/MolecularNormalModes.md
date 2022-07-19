## <a id="Psience.Molecools.Vibrations.MolecularNormalModes">MolecularNormalModes</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L286)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L286?message=Update%20Docs)]
</div>

A Coordinerds CoordinateSystem object that manages all of the data needed to
work with normal mode coordinates + some convenience functions for generating and whatnot

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
name: str
```
<a id="Psience.Molecools.Vibrations.MolecularNormalModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, coeffs, name=None, freqs=None, internal=False, origin=None, basis=None, inverse=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L292)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L292?message=Update%20Docs)]
</div>


- `molecule`: `AbstractMolecule`
    >No description...
- `coeffs`: `Any`
    >No description...
- `name`: `Any`
    >No description...
- `freqs`: `Any`
    >No description...
- `internal`: `Any`
    >No description...
- `origin`: `Any`
    >No description...
- `basis`: `Any`
    >No description...
- `inverse`: `Any`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.molecule" class="docs-object-method">&nbsp;</a> 
```python
@property
molecule(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.to_internals" class="docs-object-method">&nbsp;</a> 
```python
to_internals(self, intcrds=None, dYdR=None, dRdY=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L357)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L357?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.origin" class="docs-object-method">&nbsp;</a> 
```python
@property
origin(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.embed" class="docs-object-method">&nbsp;</a> 
```python
embed(self, frame): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L397)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L397?message=Update%20Docs)]
</div>


- `frame`: `MolecularTransformation`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.insert" class="docs-object-method">&nbsp;</a> 
```python
insert(self, val, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L461)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L461?message=Update%20Docs)]
</div>

Inserts values into the appropriate positions in the mode matrix
- `val`: `Any`
    >No description...
- `where`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.from_force_constants" class="docs-object-method">&nbsp;</a> 
```python
from_force_constants(molecule, fcs, atoms=None, masses=None, mass_units='AtomicMassUnits', inverse_mass_matrix=False, remove_transrot=True, normalize=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L514)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L514?message=Update%20Docs)]
</div>

Generates normal modes from the specified force constants
- `molecule`: `AbstractMolecule`
    >No description...
- `fcs`: `np.ndarray`
    >force constants array
- `atoms`: `Iterable[str]`
    >atom list
- `masses`: `Iterable[float]`
    >mass list
- `mass_units`: `str`
    >units for the masses...not clear if this is useful or a distraction
- `inverse_mass_matrix`: `bool`
    >whether or not we have G or G^-1 (default: `False`)
- `remove_transrot`: `bool`
    >whether or not to remove the translations and rotations (default: `True`)
- `normalize`: `bool`
    >whether or not to normalize the modes (default: `True`)
- `opts`: `Any`
    >No description...
- `:returns`: `MolecularNormalModes`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L592)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L592?message=Update%20Docs)]
</div>

Takes a slice of the modes
- `item`: `Any`
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

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Vibrations/MolecularNormalModes.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Vibrations/MolecularNormalModes.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Vibrations/MolecularNormalModes.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L286?message=Update%20Docs)