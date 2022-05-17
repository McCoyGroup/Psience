# <a id="Psience.AnalyticModels">Psience.AnalyticModels</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/tree/master/Psience/AnalyticModels)]
</div>
    


<div class="container alert alert-secondary bg-light">
  <div class="row">
   <div class="col" markdown="1">
[AnalyticPotentialConstructor](AnalyticModels/AnalyticModelConstructors/AnalyticPotentialConstructor.md)   
</div>
   <div class="col" markdown="1">
[AnalyticKineticEnergyConstructor](AnalyticModels/AnalyticModelConstructors/AnalyticKineticEnergyConstructor.md)   
</div>
   <div class="col" markdown="1">
[AnalyticModel](AnalyticModels/AnalyticModelConstructors/AnalyticModel.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[SymbolicCaller](AnalyticModels/Helpers/SymbolicCaller.md)   
</div>
   <div class="col" markdown="1">
[AnalyticModelBase](AnalyticModels/Helpers/AnalyticModelBase.md)   
</div>
</div>
</div>




<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [GmatrixElements](#GmatrixElements)

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
from Psience.AnalyticModels import *
from McUtils.Plots import *
from unittest import TestCase
import sys, h5py, math, numpy as np
```

All tests are wrapped in a test class
```python
class AnalyticModelsTests(TestCase):
    maxDiff = None
    def check_expr(self, test, real, raise_error=True):
        if not isinstance(test, int):
            test = test.expand().simplify().expand()
        if not isinstance(real, int):
            real = real.expand().simplify().expand()
        diff = (test - real)
        if not isinstance(diff, int):
            diff = diff.expand().simplify().expand()
        msg = "\nTest: {}\nReal: {}\nDiff: {}".format(test, real, diff)
        if not raise_error:
            return msg
        else:
            self.assertEquals(diff, 0, msg=msg)
```

 </div>
</div>

#### <a name="GmatrixElements">GmatrixElements</a>
```python
    def test_GmatrixElements(self):
        import sympy as sym
        class SymbolicCaller:
            def __init__(self, sym):
                self.sym = sym
            def __getitem__(self, item):
                if isinstance(item, int):
                    return self.sym(item)
                else:
                    return self.sym(*item)
        m = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_m)
        r = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_r)
        a = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_a)
        t = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_t)
        y = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_y)
        sin = sym.sin; cos = sym.cos; cot = sym.cot; tan = sym.tan
        L = SymbolicCaller(AnalyticKineticEnergyConstructor.lam)

        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 2]),
            1 / m[1] + 1 / m[2]
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 3]),
            cos(a[2,1,3])/m[1]
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [3, 4]),
            0
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 2, 3]),
            -sin(a[1,2,3])/(m[2]*r[2,3])
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 3, 4]),
            sin(a[2,1,3])*cos(t[2,1,3,4])/(m[1]*r[1,3])
        )
        self.check_expr(
            # This requires some non-automatic simplifications...
            AnalyticKineticEnergyConstructor.g([1, 2], [3, 1, 4]),
            -sin(a[3, 1, 4])/m[1]*(
                              1/r[1, 4]*cos(a[2, 1, 3])*sin(a[3, 1, 4])
                              + L[3, 1, 4]*sin(a[2, 1, 3])*cos(y[3, 1, 2, 4])
                          )
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 2, 3, 4]),
            -sin(a[1, 2, 3])*sin(t[1, 2, 3, 4])*cot(a[2, 3, 4])/(m[2]*r[2,3])
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [3, 2, 1, 4]),
            0
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 3, 4, 5]),
            -sin(a[2, 1, 3])*sin(t[2, 1, 3, 4])/(m[1]*r[1,3]*sin(a[1,3,4]))
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [4, 3, 1, 5]),
            -sin(a[2, 1, 3]) / m[1] * (
                    cot(a[1, 3, 4]) / r[1, 3] * sin(t[4, 3, 1, 2])
                     + L[5, 1, 3] * sin(t[4, 3, 1, 5] - t[4, 3, 1, 2])
            )
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2, 3], [1, 2, 3]),
            1/(m[1]*r[1,2]**2) + 1/(m[3]*r[2,3]**2)
            + 1/m[2]*(1/r[1,2]**2 + 1/r[2,3]**2 - 2*cos(a[1,2,3])/(r[1,2]*r[2,3]))
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2, 3], [2, 1, 4]),
            -1 / r[1, 2] * cos(t[3, 2, 1, 4]) * (
                L[2, 1, 4] / m[1] * sin(a[2, 1, 4]) + L[1, 2, 3] / m[2] * sin(a[1, 2, 3])
            )
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2, 3], [1, 2, 4]),
            1/m[2]*sin(a[1, 2, 3])*sin(a[1, 2, 4])*(
                L[1, 2, 3]*L[1, 2, 4]*cos(y[1, 2, 3, 4]) + 1/(r[2, 3]*r[2, 4])
            ) + 1/(m[1] * r[1, 2]**2)*cos(y[1, 2, 3, 4])
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2, 3], [1, 2, 3, 4]),
            1/r[2, 3]*sin(t[1, 2, 3, 4])*(
                L[3, 2, 1]/m[2]*sin(a[1, 2, 3])*cot(a[2, 3, 4]) - L[4, 3, 2]/m[3]
            )
        )
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/master/ci/examples/Psience/AnalyticModels.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/examples/Psience/AnalyticModels.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/Psience/AnalyticModels.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/AnalyticModels.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/AnalyticModels/__init__.py?message=Update%20Docs)