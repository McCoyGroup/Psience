## <a id="Psience.BasisReps.StateSpaces.BraKetSpace">BraKetSpace</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces.py#L3503)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L3503?message=Update%20Docs)]
</div>

Represents a set of pairs of states that can be fed into a `Representation` or `Operator`
to efficiently tell it what terms it need to calculate.
This basically just implements a bunch of stuff for generating a Graph defining
the connections between states.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
aggressive_caching_enabled: bool
preindex_trie_enabled: bool
OrthogoIndexerTrie: OrthogoIndexerTrie
CachingOrthogonalIndexCalculator: CachingOrthogonalIndexCalculator
OrthogonalIndexSparseCalculator: OrthogonalIndexSparseCalculator
OrthogonalIndexCalculator: OrthogonalIndexCalculator
use_change_indices: bool
```
<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, bra_space, ket_space, changes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3510)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3510?message=Update%20Docs)]
</div>

  - `bra_space`: `BasisStateSpace`
    > 
  - `ket_space`: `BasisStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.state_pairs" class="docs-object-method">&nbsp;</a> 
```python
@property
state_pairs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3533)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3533?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.from_indices" class="docs-object-method">&nbsp;</a> 
```python
from_indices(inds, basis=None, quanta=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3548)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3548?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3580)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3580?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3588)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3588?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.remove_duplicates" class="docs-object-method">&nbsp;</a> 
```python
remove_duplicates(self, assume_symmetric=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3592)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3592?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_space_diffs" class="docs-object-method">&nbsp;</a> 
```python
load_space_diffs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3623)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3623?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_non_orthog" class="docs-object-method">&nbsp;</a> 
```python
load_non_orthog(self, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, shared_memory_manager=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3629)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3629?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.share" class="docs-object-method">&nbsp;</a> 
```python
share(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3965)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3965?message=Update%20Docs)]
</div>
Creates a shared memory version of the `BraKetSpace`
  - `shared_memory_manager`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.unshare" class="docs-object-method">&nbsp;</a> 
```python
unshare(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3986)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3986?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3991)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3991?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.free" class="docs-object-method">&nbsp;</a> 
```python
free(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3996)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3996?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_non_orthog" class="docs-object-method">&nbsp;</a> 
```python
get_non_orthog(self, inds, assume_unique=False, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, shared_memory_manager=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4032)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4032?message=Update%20Docs)]
</div>
Returns whether the states are non-orthogonal under the set of indices.
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rules_from1d" class="docs-object-method">&nbsp;</a> 
```python
get_sel_rules_from1d(self, inds, rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4062)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4062?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rule_filter" class="docs-object-method">&nbsp;</a> 
```python
get_sel_rule_filter(self, rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4084)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4084?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4109?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4126)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4126?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_non_orthogonality" class="docs-object-method">&nbsp;</a> 
```python
apply_non_orthogonality(self, inds, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, assume_unique=False, use_change_indices=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4138)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4138?message=Update%20Docs)]
</div>
Takes the bra-ket pairs that are non-orthogonal under the indices `inds`
  - `inds`: `Any`
    > 
  - `assume_unique`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_rules_along" class="docs-object-method">&nbsp;</a> 
```python
apply_sel_rules_along(self, rules, inds, permute=True, dim=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4254)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4254?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_sums" class="docs-object-method">&nbsp;</a> 
```python
apply_sel_sums(self, rules, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4305)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4305?message=Update%20Docs)]
</div>
We reckon it's fast enough to just determine if the number
of quanta in the bra is compatible with the number of
quanta in the ket...
  - `rules`: `Any`
    > 
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_sel_rules(self, rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4327?message=Update%20Docs)]
</div>
Applies selections rules
  - `rules`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.adjacency_matrix" class="docs-object-method">&nbsp;</a> 
```python
adjacency_matrix(self, total_space=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4339)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4339?message=Update%20Docs)]
</div>
Generates the (sparse) unweighted adjacency matrix for the bras & kets
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.split" class="docs-object-method">&nbsp;</a> 
```python
split(self, chunksize): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4368)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4368?message=Update%20Docs)]
</div>
splits the brakets into blocks of at max chunksize
  - `chunksize`: `int`
    > 
  - `:returns`: `Iterable[BraKetSpace]`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.concatenate" class="docs-object-method">&nbsp;</a> 
```python
concatenate(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4384)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4384?message=Update%20Docs)]
</div>
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-08542a" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-08542a"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-08542a" markdown="1">
 - [HarmHam](#HarmHam)
- [BasisRepMatrixOps](#BasisRepMatrixOps)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-49f86d" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-49f86d"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-49f86d" markdown="1">
 
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

#### <a name="HarmHam">HarmHam</a>
```python
    def test_HarmHam(self):

        n = 10
        m = 3
        basis = HarmonicOscillatorProductBasis((n,) * m)
        G, V = [
            np.array([[6.47886479e-03, 5.17641431e-12, -1.12922679e-12],
                      [5.17641431e-12, 1.28034398e-02, -3.15629792e-12],
                      [-1.12922679e-12, -3.15629792e-12, 1.76505371e-02]]),
            np.array([[6.47886478e-03, -8.45595180e-13, -1.01327126e-11],
                      [-8.45595549e-13, 1.28034398e-02, -4.72136245e-12],
                      [-1.01327124e-11, -4.72136255e-12, 1.76505372e-02]])]

        mommy = (1 / 2) * basis.representation('p', 'p', coeffs=G)
        possy = (1 / 2) * basis.representation('x', 'x', coeffs=V)
        H0 = ( mommy + possy )

        states = BasisStateSpace(basis, self.get_states(2, 3, max_quanta=10), mode='excitations')

        diag_inds = BraKetSpace(states, states)

        # raise Exception(diag_inds.state_pairs)

        diags = H0[diag_inds]

        self.assertEquals(np.average(diags), 0.036932841734999985)
```

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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/BraKetSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/BraKetSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/BraKetSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L3503?message=Update%20Docs)   
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