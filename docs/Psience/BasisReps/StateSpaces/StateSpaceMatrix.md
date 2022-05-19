## <a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix">StateSpaceMatrix</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4058)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4058?message=Update%20Docs)]
</div>

A `SparseArray` that holds onto a `BasisStateSpace` that keeps track of the
total set of states involved.
By default is assumed real-symmetric. This can be relaxed in the future.

TODO: The idea is good, but calculating what is "in" the array and what is "out"
        every single time this is applied could be slow...
      We'll need to test to see how slow

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, initial_basis, initial_vals=None, column_space=None, symmetric=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4069)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4069?message=Update%20Docs)]
</div>


- `initial_basis`: `BasisStateSpace | RepresentationBasis`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

Returns the basis for the matrix rep
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.brakets" class="docs-object-method">&nbsp;</a> 
```python
@property
brakets(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

Returns the BraKetSpace for the held indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.identity_from_space" class="docs-object-method">&nbsp;</a> 
```python
identity_from_space(space, column_space=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4155)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4155?message=Update%20Docs)]
</div>

Returns a StateSpaceMatrix where the diagonal is filled with 1s
- `space`: `Any`
    >No description...
- `column_space`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.extend_basis" class="docs-object-method">&nbsp;</a> 
```python
extend_basis(self, states, extend_columns=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4193)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4193?message=Update%20Docs)]
</div>

Extends the held state space and resizes the held array if need be
- `states`: `BasisStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.compute_values" class="docs-object-method">&nbsp;</a> 
```python
compute_values(self, func, brakets): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4276)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4276?message=Update%20Docs)]
</div>

Computes new values into the held `SparseArray` based on the function and brakets provided
        and returns the entire array of values
- `func`: `Any`
    >A function that can take a braket spec and compute values
- `brakets`: `Any`
    >A set of brakets to compute values for
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.dot" class="docs-object-method">&nbsp;</a> 
```python
dot(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4328)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4328?message=Update%20Docs)]
</div>

Performs a dot product between the held SparseArray and another
        StateSpaceMatrix
- `other`: `StateSpaceMatrix`
    >other matrix
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4364)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4364?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__setitem__" class="docs-object-method">&nbsp;</a> 
```python
__setitem__(self, item, vals): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4375)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4375?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L4386)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4386?message=Update%20Docs)]
</div>

 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [BasisRepMatrixOps](#BasisRepMatrixOps)
- [ImprovedRepresentations](#ImprovedRepresentations)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
#### <a class="collapse-link" data-toggle="collapse" href="#test-setup">Setup</a> <a class="float-right" data-toggle="collapse" href="#test-setup"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="test-setup" markdown="1">

Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.
```python
from Peeves import Timer, BlockProfiler
from McUtils.Scaffolding import *
import McUtils.Plots as plt
from McUtils.Combinatorics import CompleteSymmetricGroupSpace
from Peeves.TestUtils import *
from unittest import TestCase
from Psience.BasisReps import *
import sys, os, numpy as np
```

All tests are wrapped in a test class
```python
class BasisSetTests(TestCase):
    def get_states(self, n_quanta, n_modes, max_quanta=None):
        return [np.flip(x) for x in BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(n_modes),
            range(n_quanta)
        ).excitations]
```

 </div>
</div>

#### <a name="BasisRepMatrixOps">BasisRepMatrixOps</a>
```python
    def test_BasisRepMatrixOps(self):

        n = 15 # totally meaningless these days
        m = 4
        basis = HarmonicOscillatorProductBasis((n,) * m)

        mat = StateSpaceMatrix(basis)

        self.assertEquals(mat.array.shape[0], 0)

        states = BasisStateSpace.from_quanta(basis, range(10))
        brakets = BraKetSpace(states, states)
        vals = np.ones(len(brakets))
        mat_2 = StateSpaceMatrix(brakets, vals)

        def wat(state_space):
            return np.ones(len(state_space))
        sub_brakets = BasisStateSpace.from_quanta(basis, range(4)).get_representation_brakets()
        mat2_vals = mat_2.compute_values(wat, sub_brakets)

        self.assertEquals(mat2_vals.tolist(), mat_2[sub_brakets].tolist())
```
#### <a name="ImprovedRepresentations">ImprovedRepresentations</a>
```python
    def test_ImprovedRepresentations(self):
        n = 15  # totally meaningless these days
        m = 4
        basis = HarmonicOscillatorProductBasis((n,) * m)

        x_rep = basis.representation('x', coeffs=np.ones((m,)))
        initial_space = basis.get_state_space(range(4))
        initial_states = StateSpaceMatrix.identity_from_space(initial_space)

        x = x_rep.apply(initial_states)
        x2 = x_rep.apply(x)
        # plt.ArrayPlot(x.array.asarray(), aspect_ratio=x.shape[0]/x.shape[1], image_size=300)
        # plt.ArrayPlot(x2.array.asarray(), aspect_ratio=x2.shape[0]/x2.shape[1], image_size=300).show()

        x2_rep = basis.representation('x', 'x', coeffs=np.ones((m, m)))
        x22 = x2_rep.apply(initial_states)

        plt.ArrayPlot(x2.array.asarray(), aspect_ratio=x2.shape[0]/x2.shape[1], image_size=200)
        plt.ArrayPlot(x22.array.asarray(), aspect_ratio=x22.shape[0]/x22.shape[1], image_size=200).show()

        self.assertTrue(np.allclose(x2.array.asarray(), x22.array.asarray()))
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L4058?message=Update%20Docs)