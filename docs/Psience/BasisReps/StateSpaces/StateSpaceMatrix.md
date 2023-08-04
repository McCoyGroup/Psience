## <a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix">StateSpaceMatrix</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces.py#L4397)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L4397?message=Update%20Docs)]
</div>

A `SparseArray` that holds onto a `BasisStateSpace` that keeps track of the
total set of states involved.
By default is assumed real-symmetric. This can be relaxed in the future.

TODO: The idea is good, but calculating what is "in" the array and what is "out"
every single time this is applied could be slow...
We'll need to test to see how slow







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, initial_basis, initial_vals=None, column_space=None, symmetric=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4408)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4408?message=Update%20Docs)]
</div>

  - `initial_basis`: `BasisStateSpace | RepresentationBasis`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4467)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4467?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4473)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4473?message=Update%20Docs)]
</div>
Returns the basis for the matrix rep
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.brakets" class="docs-object-method">&nbsp;</a> 
```python
@property
brakets(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4483)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4483?message=Update%20Docs)]
</div>
Returns the BraKetSpace for the held indices
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.identity_from_space" class="docs-object-method">&nbsp;</a> 
```python
identity_from_space(space, column_space=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4494)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4494?message=Update%20Docs)]
</div>
Returns a StateSpaceMatrix where the diagonal is filled with 1s
  - `space`: `Any`
    > 
  - `column_space`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.extend_basis" class="docs-object-method">&nbsp;</a> 
```python
extend_basis(self, states, extend_columns=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4532)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4532?message=Update%20Docs)]
</div>
Extends the held state space and resizes the held array if need be
  - `states`: `BasisStateSpace`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.compute_values" class="docs-object-method">&nbsp;</a> 
```python
compute_values(self, func, brakets): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4615)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4615?message=Update%20Docs)]
</div>
Computes new values into the held `SparseArray` based on the function and brakets provided
and returns the entire array of values
  - `func`: `Any`
    > A function that can take a braket spec and compute values
  - `brakets`: `Any`
    > A set of brakets to compute values for
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.dot" class="docs-object-method">&nbsp;</a> 
```python
dot(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4667)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4667?message=Update%20Docs)]
</div>
Performs a dot product between the held SparseArray and another
StateSpaceMatrix
  - `other`: `StateSpaceMatrix`
    > other matrix
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4703)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4703?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__setitem__" class="docs-object-method">&nbsp;</a> 
```python
__setitem__(self, item, vals): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4714)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4714?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4725)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4725?message=Update%20Docs)]
</div>
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-27546c" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-27546c"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-27546c" markdown="1">
 - [BasisRepMatrixOps](#BasisRepMatrixOps)
- [ImprovedRepresentations](#ImprovedRepresentations)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-a7fe4c" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-a7fe4c"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-a7fe4c" markdown="1">
 
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L4397?message=Update%20Docs)   
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