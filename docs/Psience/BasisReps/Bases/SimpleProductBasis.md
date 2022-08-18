## <a id="Psience.BasisReps.Bases.SimpleProductBasis">SimpleProductBasis</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L310)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L310?message=Update%20Docs)]
</div>

Defines a direct product basis from a 1D basis.
Mixed product bases aren't currently supported, but this provides
at least a sample for how that kind of things could be
generated.



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
array_indexer_cutoff: int
```
<a id="Psience.BasisReps.Bases.SimpleProductBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis_type, n_quanta, indexer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L320)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L320?message=Update%20Docs)]
</div>


- `indexer`: `BaseStateIndexer`
    >an object that can turn state specs into indices and indices into state specs
- `n_quanta`: `Iterable[int]`
    >the number of quanta for the representations
- `basis_type`: `type`
    >the type of basis to do a product over

<a id="Psience.BasisReps.Bases.SimpleProductBasis.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L348)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L348?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Bases.SimpleProductBasis.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L354)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L354?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Bases.SimpleProductBasis.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L?message=Update%20Docs)]
</div>

Provides the number of dimensions of the basis
- `:returns`: `_`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.dimensions" class="docs-object-method">&nbsp;</a> 
```python
@property
dimensions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L?message=Update%20Docs)]
</div>

Provides the dimensions of the basis
- `:returns`: `_`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L380)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L380?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >
- `other`: `SimpleProductBasis`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.quanta" class="docs-object-method">&nbsp;</a> 
```python
@property
quanta(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L?message=Update%20Docs)]
</div>

Provides the quanta in each dimension of the basis
- `:returns`: `_`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.selection_rules_mapping" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules_mapping(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Bases.SimpleProductBasis.ravel_state_inds" class="docs-object-method">&nbsp;</a> 
```python
ravel_state_inds(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L410)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L410?message=Update%20Docs)]
</div>

Converts state indices from an array of quanta to an array of indices
- `:returns`: `tuple[int]`
    >a
r
r
a
y
 
o
f
 
s
t
a
t
e
 
i
n
d
i
c
e
s
 
i
n
 
t
h
e
 
b
a
s
i
s
- `idx`: `Iterable[Iterable[int]]`
    >indices

<a id="Psience.BasisReps.Bases.SimpleProductBasis.unravel_state_inds" class="docs-object-method">&nbsp;</a> 
```python
unravel_state_inds(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L422)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L422?message=Update%20Docs)]
</div>

Converts state indices from an array of ints to an array of quanta
- `:returns`: `tuple[tuple[int]]`
    >a
r
r
a
y
 
o
f
 
s
t
a
t
e
 
t
u
p
l
e
s
 
i
n
 
t
h
e
 
b
a
s
i
s
- `idx`: `Iterable[int]`
    >indices

<a id="Psience.BasisReps.Bases.SimpleProductBasis.get_function" class="docs-object-method">&nbsp;</a> 
```python
get_function(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L436)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L436?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Bases.SimpleProductBasis.operator" class="docs-object-method">&nbsp;</a> 
```python
operator(self, *terms, coeffs=None, axes=None, parallelizer=None, logger=None, chunk_size=None, **operator_settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L440)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L440?message=Update%20Docs)]
</div>

Builds an operator based on supplied terms, remapping names where possible.
If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.
- `:returns`: `_`
    >
- `axes`: `Any`
    >
- `coeffs`: `Any`
    >
- `terms`: `Any`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.representation" class="docs-object-method">&nbsp;</a> 
```python
representation(self, *terms, coeffs=None, name=None, axes=None, logger=None, parallelizer=None, chunk_size=None, memory_constrained=False, **operator_settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L484)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L484?message=Update%20Docs)]
</div>

Provides a representation of a product operator specified by _terms_.
If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.
- `:returns`: `_`
    >
- `terms`: `Any`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.x" class="docs-object-method">&nbsp;</a> 
```python
x(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L506)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L506?message=Update%20Docs)]
</div>

Returns the representation of x in the multi-dimensional basis with every term evaluated up to n quanta
Whether this is what we want or not is still TBD
- `:returns`: `_`
    >
- `n`: `Any`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.p" class="docs-object-method">&nbsp;</a> 
```python
p(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L516)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L516?message=Update%20Docs)]
</div>

Returns the representation of p in the multi-dimensional basis with every term evaluated up to n quanta
Whether this is what we want or not is still TBD
- `:returns`: `_`
    >
- `n`: `Any`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, dims): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L527)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L527?message=Update%20Docs)]
</div>

Casts down to lower dimensional space
- `:returns`: `_`
    >
- `dims`: `Any`
    >

<a id="Psience.BasisReps.Bases.SimpleProductBasis.get_state_space" class="docs-object-method">&nbsp;</a> 
```python
get_state_space(self, quanta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Bases.py#L539)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L539?message=Update%20Docs)]
</div>

 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [HOBasis2DXX](#HOBasis2DXX)
- [HOBasis2DPP](#HOBasis2DPP)
- [HOBasis3DPXP](#HOBasis3DPXP)
- [HOBasis3DXXX](#HOBasis3DXXX)
- [HOBasis3DXXX2D](#HOBasis3DXXX2D)
- [HOBasis3DXXX2DContracted](#HOBasis3DXXX2DContracted)
- [HOBasis4DPXXP](#HOBasis4DPXXP)

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

#### <a name="HOBasis2DXX">HOBasis2DXX</a>
```python
    def test_HOBasis2DXX(self):

        n = 7
        m = 10
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)
        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        states = (
            (0, 0, 0, 0, 0),
            (0, 1, 2, 3, 4)
        )
        # with Timer("New style"):
            # with BlockProfiler("New Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals1 = xxpp1[states]
        self.assertEquals(vals1.shape, (m,) * n_terms + (len(states[0]),))

        # with Timer("Old style"):
            # with BlockProfiler("Old Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals2 = xxpp2[states]
        self.assertEquals(vals2.shape, (m,) * n_terms + (len(states[0]),))

        wat = np.roll(np.arange(n_terms + 1), 1)
        # print(wat)
        v1 = vals1.toarray().transpose(wat)
        v2 = iphase * vals2.toarray().transpose(wat)

        # print([np.max(v) for v in v1])
        # print([np.max(v) for v in v2])
        # print([np.max(v) for v in np.abs(v1 - v2)])
        # print(v1[0], v2[0])
        # print(v1[1], v2[1])
        # print(vals1.toarray()[:, :, -1] - vals2.toarray()))

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)
```
#### <a name="HOBasis2DPP">HOBasis2DPP</a>
```python
    def test_HOBasis2DPP(self):
        from Peeves import Timer, BlockProfiler

        n = 10
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)

        g1 = np.array(
            [[-1.81146079e-04, 3.97836803e-05],
             [3.97836803e-05, 2.63572358e-05]])
        xxpp1 = 2 * oppo.representation(*term, coeffs=g1, axes=[[0, 1], [1, 0]])
        xxpp1 = xxpp1 + xxpp1
        xxpp2 = 2 * oppo2.representation(*term, coeffs=g1, axes=[[0, 1], [1, 0]])
        xxpp2 = xxpp2 + xxpp2

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        quant_states = BasisStateSpace(
            oppo,
            self.get_states(9, m, max_quanta=10)
        )
        brakets = quant_states.get_representation_brakets()

        # with Timer("New style"):
            # with BlockProfiler("New Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals1 = xxpp1[brakets]

        # with Timer("Old style"):
            # with BlockProfiler("Old Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals2 = xxpp2[brakets]

        v1 = vals1
        v2 = iphase * vals2

        # n = len(quant_states)
        # plt.ArrayPlot(v1.reshape((n, n)))
        # plt.ArrayPlot(v2.reshape((n, n))).show()

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)
```
#### <a name="HOBasis3DPXP">HOBasis3DPXP</a>
```python
    def test_HOBasis3DPXP(self):
        from Peeves import Timer, BlockProfiler

        n = 10
        m = 3
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'x', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)

        g1 = np.array([
            [[-1.81146079e-04,  3.97836803e-05,  2.91649691e-05],
             [ 3.97836803e-05,  2.63572358e-05,  2.37597837e-04],
             [ 2.91649691e-05,  2.37597837e-04, -3.38457268e-05]],

            [[-4.36589189e-04,  2.79004059e-05, -1.50059967e-05],
             [ 2.79004059e-05, -1.44188965e-06,  3.49657651e-06],
             [-1.50059967e-05,  3.49657652e-06,  3.11501367e-06]],

            [[-8.10821036e-04,  6.31615150e-06,  5.50255712e-05],
             [ 6.31615151e-06,  4.05569426e-06,  3.51303496e-08],
             [ 5.50255712e-05,  3.51303696e-08, -3.80070492e-06]]])
        xxpp1 = oppo.representation(*term, coeffs=g1, axes=[[0, 1, 2], [1, 0, 2]])
        # xxpp1 = xxpp1 + xxpp1
        xxpp2 = oppo2.representation(*term, coeffs=g1,  axes=[[0, 1, 2], [1, 0, 2]])
        # xxpp2 = xxpp2 + xxpp2

        quant_states = BasisStateSpace(
            oppo,
            self.get_states(3, 3, max_quanta=10)
        )
        inds = quant_states.get_representation_brakets()

        #
        # # raise Exception(inds.bras.indices)
        #
        # quant_states = BasisStateSpace(
        #     oppo,
        #     self.get_states(3, 3, max_quanta=10)
        # )
        # new_stuff = quant_states.apply_selection_rules([[-1, 1]])
        # inds2 = new_stuff.get_representation_brakets()
        #
        # plt.ArrayPlot(inds2.adjacency_matrix().toarray()).show()


        # inds = BasisStateSpace(
        #     oppo,
        #     (
        #         [0, 0, 0],
        #         [1, 0, 0],
        #     )
        # ).get_representation_brakets()


        # with Timer("New style"):
        vals1 = xxpp1[inds]
        # with Timer("Old style"):
        vals2 = xxpp2[inds]

        v1 = vals1
        v2 = iphase * vals2

        # n = len(quant_states)
        # plt.ArrayPlot(v1.reshape((n, n)))
        # plt.ArrayPlot(v2.reshape((n, n)))
        # plt.ArrayPlot(v1.reshape((n, n)) - v1.reshape((n, n)).T,
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5))
        # plt.ArrayPlot(v2.reshape((n, n)) - v2.reshape((n, n)).T,
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5))
        # plt.ArrayPlot((v1 - v2).reshape((n, n)),
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5)).show()

        self.assertTrue(
            np.allclose(
                v1[:15],
                [0.00000000e+00, -2.86578374e-04, 0.00000000e+00, 3.29150701e-06,
                 -1.53766049e-04, 0.00000000e+00, -1.59263719e-06, 0.00000000e+00,
                 -5.52442364e-06, 1.24871307e-06, -6.66923918e-05, 0.00000000e+00,
                 -3.81027078e-05, 0.00000000e+00, -1.61862393e-04],
                atol=1.0e-5
            ))

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)
```
#### <a name="HOBasis3DXXX">HOBasis3DXXX</a>
```python
    def test_HOBasis3DXXX(self):

        n = 15
        m = 5
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        quant_states = self.get_states(4, m)
        states = oppo.ravel_state_inds(quant_states)
        # print(quant_states)
        import itertools as ip
        wat = np.array(list(ip.product(states, states))).T
        # with Timer("New style"):
        vals1 = xxpp1[wat[0], wat[1]]

        # with Timer("Old style"):
        vals2 = xxpp2[wat[0], wat[1]]

        v1 = vals1.asarray()
        v2 = iphase * vals2.asarray()

        # with JSONCheckpointer(os.path.expanduser("~/Desktop/test_terms.json")) as chk:
        #     chk['XXX_3D_new'] = v1
        #     chk['XXX_3D_old'] = v2

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)
```
#### <a name="HOBasis3DXXX2D">HOBasis3DXXX2D</a>
```python
    def test_HOBasis3DXXX2D(self):

        n = 15
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        states = BasisStateSpace.from_quanta(oppo, range(10))
        brakets = states.get_representation_brakets()
        vals1 = xxpp1[brakets]
        vals2 = xxpp2[brakets]

        v1 = vals1.asarray()
        v2 = iphase * vals2.asarray()

        # with JSONCheckpointer(os.path.expanduser("~/Desktop/test_terms.json")) as chk:
        #     chk['XXX_exc'] = states.excitations
        #     chk['XXX_3D_new'] = v1
        #     chk['XXX_3D_old'] = v2

        self.assertLess(np.max(np.abs(v1 - v2)), 2.0e-14)
```
#### <a name="HOBasis3DXXX2DContracted">HOBasis3DXXX2DContracted</a>
```python
    def test_HOBasis3DXXX2DContracted(self):
        n = 15
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term, coeffs=np.ones((m, m, m)))
        xxpp2 = oppo2.representation(*term, coeffs=np.ones((m, m, m)))


        states = BasisStateSpace.from_quanta(oppo, range(10))
        brakets = states.get_representation_brakets()
        vals1 = xxpp1[brakets]
        vals2 = xxpp2[brakets]

        v1 = vals1
        v2 = iphase * vals2

        # with JSONCheckpointer(os.path.expanduser("~/Desktop/test_terms.json")) as chk:
        #     chk['XXX_exc'] = states.excitations
        #     chk['XXX_3D_new'] = v1
        #     chk['XXX_3D_old'] = v2

        self.assertLess(np.max(np.abs(v1 - v2)), 2.0e-14)
```
#### <a name="HOBasis4DPXXP">HOBasis4DPXXP</a>
```python
    def test_HOBasis4DPXXP(self):
        from Peeves import Timer, BlockProfiler

        n = 15
        m = 5
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'x', 'x', 'p']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        quant_states = self.get_states(4, m)
        states = oppo.ravel_state_inds(quant_states)
        # print(quant_states)
        import itertools as ip
        wat = np.array(list(ip.product(states, states))).T
        # with Timer("New style"):
        vals1 = xxpp1[wat[0], wat[1]]

        # with Timer("Old style"):
        vals2 = xxpp2[wat[0], wat[1]]

        v1 = vals1.toarray()
        v2 = iphase * vals2.toarray()

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Bases/SimpleProductBasis.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Bases/SimpleProductBasis.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Bases/SimpleProductBasis.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Bases/SimpleProductBasis.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Bases.py#L310?message=Update%20Docs)