## <a id="Psience.BasisReps.Representations.Representation">Representation</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations.py#L25)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations.py#L25?message=Update%20Docs)]
</div>

A `Representation` provides a simple interface to build matrix representations of operators expressed
in high-dimensional spaces.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.Representations.Representation.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, compute, basis, name=None, logger=None, selection_rules=None, selection_rule_steps=None, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L32?message=Update%20Docs)]
</div>

  - `compute`: `callable | Operator`
    > the function that turns indices into values
  - `basis`: `RepresentationBasis`
    > the basis quanta used in the representations
  - `logger`: `None | Logger`
    > logger for printing out debug info


<a id="Psience.BasisReps.Representations.Representation.parallelizer" class="docs-object-method">&nbsp;</a> 
```python
@property
parallelizer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L64)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L64?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.compute" class="docs-object-method">&nbsp;</a> 
```python
compute(self, inds, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L71)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L71?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.compute_cached" class="docs-object-method">&nbsp;</a> 
```python
compute_cached(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L77)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L77?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.chunk_size" class="docs-object-method">&nbsp;</a> 
```python
@property
chunk_size(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L83)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L83?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L90)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L90?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.diag" class="docs-object-method">&nbsp;</a> 
```python
@property
diag(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L105)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L105?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.ndims" class="docs-object-method">&nbsp;</a> 
```python
@property
ndims(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L108)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L108?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.dim_inds" class="docs-object-method">&nbsp;</a> 
```python
@property
dim_inds(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L111)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L111?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.get_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_brakets(self, states, check_orthogonality=True, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L139)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L139?message=Update%20Docs)]
</div>
Computes term elements based on getting a BraKetSpace.
Can directly pass element specs through, since the shape management shit
is managed by the BraKetSpace
  - `states`: `BraKetSpace | Tuple[np.ndarray, np.ndarray]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.Representation.get_element" class="docs-object-method">&nbsp;</a> 
```python
get_element(self, n, m): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L160)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L160?message=Update%20Docs)]
</div>
Computes term elements.
Determines first whether it needs to pull single elements or blocks of them.
  - `n`: `Any`
    > 
  - `m`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.Representation.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L250)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L250?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L267)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L267?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L279)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L279?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L292)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L292?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L325)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L325?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.selection_rules" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L337)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L337?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.Representation.selection_rule_steps" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rule_steps(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L358)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L358?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.Representation.is_diagonal" class="docs-object-method">&nbsp;</a> 
```python
@property
is_diagonal(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L373)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L373?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.is_zero" class="docs-object-method">&nbsp;</a> 
```python
@property
is_zero(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L379)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L379?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.skipped_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
skipped_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L386)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L386?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.Representation.get_transformed_space" class="docs-object-method">&nbsp;</a> 
```python
get_transformed_space(self, space, parallelizer=None, logger=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L398)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L398?message=Update%20Docs)]
</div>
Returns the state space obtained by using the
held operator to transform `space`
  - `space`: `Any`
    > 
  - `:returns`: `SelectionRuleStateSpace`
    > c
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


<a id="Psience.BasisReps.Representations.Representation.apply" class="docs-object-method">&nbsp;</a> 
```python
apply(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L428)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L428?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.Representation.get_representation_matrix" class="docs-object-method">&nbsp;</a> 
```python
get_representation_matrix(self, coupled_space, total_space, filter_space=None, diagonal=False, logger=None, zero_element_warning=True, clear_sparse_caches=True, clear_operator_caches=True, assume_symmetric=True, remove_duplicates=True, memory_constrained=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L589)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L589?message=Update%20Docs)]
</div>
Actively constructs a perturbation theory Hamiltonian representation
  - `h`: `Any`
    > 
  - `cs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.Representation.get_diagonal_representation" class="docs-object-method">&nbsp;</a> 
```python
get_diagonal_representation(self, coupled_space, total_space, logger=None, zero_element_warning=True, clear_sparse_caches=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/Representation.py#L826)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/Representation.py#L826?message=Update%20Docs)]
</div>
Actively constructs a perturbation theory Hamiltonian representation
  - `h`: `Any`
    > 
  - `cs`: `Any`
    > 
  - `:returns`: `_`
    >
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-a79f87" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-a79f87"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-a79f87" markdown="1">
 - [HOBasis1DX](#HOBasis1DX)
- [HOBasis1DXX](#HOBasis1DXX)
- [HOBasis1DPXP](#HOBasis1DPXP)
- [HOBasis1DPP](#HOBasis1DPP)
- [HOBasis1DXXX](#HOBasis1DXXX)
- [HOBasis1DPPXX](#HOBasis1DPPXX)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-e068b5" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-e068b5"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-e068b5" markdown="1">
 
Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.

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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Representations/Representation.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Representations/Representation.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Representations/Representation.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Representations/Representation.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations.py#L25?message=Update%20Docs)   
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