## <a id="Psience.BasisReps.Representations.Representation">Representation</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L25)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L25?message=Update%20Docs)]
</div>

A `Representation` provides a simple interface to build matrix representations of operators expressed
in high-dimensional spaces.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.Representations.Representation.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, compute, basis, name=None, logger=None, selection_rules=None, selection_rule_steps=None, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L32)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L32?message=Update%20Docs)]
</div>


- `logger`: `None | Logger`
    >logger for printing out debug info
- `basis`: `RepresentationBasis`
    >the basis quanta used in the representations
- `compute`: `callable | Operator`
    >the function that turns indices into values

<a id="Psience.BasisReps.Representations.Representation.parallelizer" class="docs-object-method">&nbsp;</a> 
```python
@property
parallelizer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.compute" class="docs-object-method">&nbsp;</a> 
```python
compute(self, inds, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L71)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L71?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.compute_cached" class="docs-object-method">&nbsp;</a> 
```python
compute_cached(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L77)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L77?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.chunk_size" class="docs-object-method">&nbsp;</a> 
```python
@property
chunk_size(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L90)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L90?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.diag" class="docs-object-method">&nbsp;</a> 
```python
@property
diag(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.ndims" class="docs-object-method">&nbsp;</a> 
```python
@property
ndims(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.dim_inds" class="docs-object-method">&nbsp;</a> 
```python
@property
dim_inds(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.get_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_brakets(self, states, check_orthogonality=True, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L139)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L139?message=Update%20Docs)]
</div>

Computes term elements based on getting a BraKetSpace.
Can directly pass element specs through, since the shape management shit
is managed by the BraKetSpace
- `:returns`: `_`
    >
- `states`: `BraKetSpace | Tuple[np.ndarray, np.ndarray]`
    >

<a id="Psience.BasisReps.Representations.Representation.get_element" class="docs-object-method">&nbsp;</a> 
```python
get_element(self, n, m): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L160)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L160?message=Update%20Docs)]
</div>

Computes term elements.
Determines first whether it needs to pull single elements or blocks of them.
- `:returns`: `_`
    >
- `m`: `Any`
    >
- `n`: `Any`
    >

<a id="Psience.BasisReps.Representations.Representation.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L250)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L250?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L267)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L267?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L279)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L279?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L292)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L292?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L325)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L325?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.selection_rules" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >

<a id="Psience.BasisReps.Representations.Representation.selection_rule_steps" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rule_steps(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >

<a id="Psience.BasisReps.Representations.Representation.is_diagonal" class="docs-object-method">&nbsp;</a> 
```python
@property
is_diagonal(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.is_zero" class="docs-object-method">&nbsp;</a> 
```python
@property
is_zero(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.skipped_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
skipped_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >

<a id="Psience.BasisReps.Representations.Representation.get_transformed_space" class="docs-object-method">&nbsp;</a> 
```python
get_transformed_space(self, space, parallelizer=None, logger=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L398)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L398?message=Update%20Docs)]
</div>

Returns the state space obtained by using the
held operator to transform `space`
- `:returns`: `SelectionRuleStateSpace`
    >c
o
n
n
e
c
t
e
d
 
s
t
a
t
e
 
s
p
a
c
e
s
- `space`: `Any`
    >

<a id="Psience.BasisReps.Representations.Representation.apply" class="docs-object-method">&nbsp;</a> 
```python
apply(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L428)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L428?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.get_representation_matrix" class="docs-object-method">&nbsp;</a> 
```python
get_representation_matrix(self, coupled_space, total_space, filter_space=None, diagonal=False, logger=None, zero_element_warning=True, clear_sparse_caches=True, clear_operator_caches=True, assume_symmetric=True, remove_duplicates=True, memory_constrained=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L589)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L589?message=Update%20Docs)]
</div>

Actively constructs a perturbation theory Hamiltonian representation
- `:returns`: `_`
    >
- `cs`: `Any`
    >
- `h`: `Any`
    >

<a id="Psience.BasisReps.Representations.Representation.get_diagonal_representation" class="docs-object-method">&nbsp;</a> 
```python
get_diagonal_representation(self, coupled_space, total_space, logger=None, zero_element_warning=True, clear_sparse_caches=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Representations.py#L826)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L826?message=Update%20Docs)]
</div>

Actively constructs a perturbation theory Hamiltonian representation
- `:returns`: `_`
    >
- `cs`: `Any`
    >
- `h`: `Any`
    >

 </div>
</div>





<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [HOBasis1DX](#HOBasis1DX)
- [HOBasis1DXX](#HOBasis1DXX)
- [HOBasis1DPXP](#HOBasis1DPXP)
- [HOBasis1DPP](#HOBasis1DPP)
- [HOBasis1DXXX](#HOBasis1DXXX)
- [HOBasis1DPPXX](#HOBasis1DPPXX)

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

#### <a name="HOBasis1DX">HOBasis1DX</a>
```python
    def test_HOBasis1DX(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        term = ['x']
        iphase = (-1) ** (term.count("p") // 2)
        rep1 = basis.representation(*term)
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator(*term), basis)
        xx = rep1[:, :].todense()
        x2 = iphase * rep2[:, :].todense()

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)
```
#### <a name="HOBasis1DXX">HOBasis1DXX</a>
```python
    def test_HOBasis1DXX(self):

        n = 7
        basis = HarmonicOscillatorBasis(n)

        rep1 = basis.representation('x', 'x')
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator('x', 'x'), basis)
        xx = rep1[:, :].todense()
        x2 = rep2[:, :].todense()
        targ =np.zeros((n, n))
        targ[np.arange(n), np.arange(n)] = np.arange(n) + 1/2
        targ[np.arange(n-2), np.arange(2, n)] = np.sqrt(np.arange(1, n-1)*(np.arange(1, n-1)+1)/4)
        targ[np.arange(2, n), np.arange(n-2)] = np.sqrt(np.arange(1, n-1)*(np.arange(1, n-1)+1)/4)

        # raise Exception([
        #     targ**2,
        #     xx**2
        #     ])

        self.assertLess(np.average(np.abs(x2 - targ)), 1e-14)
        self.assertLess(np.average(np.abs(xx - targ)), 1e-14)
```
#### <a name="HOBasis1DPXP">HOBasis1DPXP</a>
```python
    def test_HOBasis1DPXP(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        term = ['p', 'x', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        rep1 = basis.representation(*term)
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator(*term), basis)
        xx = rep1[:, :].todense()
        x2 = iphase * rep2[:, :].todense()

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)
```
#### <a name="HOBasis1DPP">HOBasis1DPP</a>
```python
    def test_HOBasis1DPP(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        rep1 = basis.representation('p', 'p')
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator('p', 'p'), basis)
        xx = rep1[:, :].todense()
        x2 = -rep2[:, :].todense()
        # targ = np.zeros((n, n))
        # targ[np.arange(n), np.arange(n)] = np.arange(n) + 1 / 2
        # targ[np.arange(n - 2), np.arange(2, n)] = np.sqrt(np.arange(1, n - 1) * (np.arange(1, n - 1) + 1) / 4)
        # targ[np.arange(2, n), np.arange(n - 2)] = np.sqrt(np.arange(1, n - 1) * (np.arange(1, n - 1) + 1) / 4)

        # raise Exception([
        #     targ**2,
        #     xx**2
        #     ])

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)
```
#### <a name="HOBasis1DXXX">HOBasis1DXXX</a>
```python
    def test_HOBasis1DXXX(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)
        rep1 = basis.representation(*term)
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator(*term), basis)
        xx = rep1[:, :].todense()
        x2 = iphase * rep2[:, :].todense()

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)
```
#### <a name="HOBasis1DPPXX">HOBasis1DPPXX</a>
```python
    def test_HOBasis1DPPXX(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        term = ['p', 'p', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)
        rep1 = basis.representation(*term)
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator(*term), basis)
        xx = rep1[:, :].todense()
        x2 = iphase * rep2[:, :].todense()
        # targ = np.zeros((n, n))
        # targ[np.arange(n), np.arange(n)] = np.arange(n) + 1 / 2
        # targ[np.arange(n - 2), np.arange(2, n)] = np.sqrt(np.arange(1, n - 1) * (np.arange(1, n - 1) + 1) / 4)
        # targ[np.arange(2, n), np.arange(n - 2)] = np.sqrt(np.arange(1, n - 1) * (np.arange(1, n - 1) + 1) / 4)

        # raise Exception([
        #     targ**2,
        #     xx**2
        #     ])

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Representations/Representation.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Representations/Representation.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Representations/Representation.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Representations/Representation.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Representations.py#L25?message=Update%20Docs)