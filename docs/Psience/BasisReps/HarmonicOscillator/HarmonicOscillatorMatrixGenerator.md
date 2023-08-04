## <a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator">HarmonicOscillatorMatrixGenerator</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator.py#L328)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator.py#L328?message=Update%20Docs)]
</div>

1D evaluator for terms looking like `x`, `p`, `q`, etc.
All of the overall `(-i)^N` info is in the `ProdOp` class that's expected to hold this.
Only maintains phase info & calculates elements.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
state_cache_size: int
default_evaluator_mode: str
```
<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, terms, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L337)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L337?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L355)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L355?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L363)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L363?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.set_cache_size" class="docs-object-method">&nbsp;</a> 
```python
set_cache_size(new_size): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L366)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L366?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.load_cached" class="docs-object-method">&nbsp;</a> 
```python
load_cached(terms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L370)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L370?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.selection_rules" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L382)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L382?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L386)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L386?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.state_pair_hash" class="docs-object-method">&nbsp;</a> 
```python
state_pair_hash(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L390)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L390?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.pull_state_groups" class="docs-object-method">&nbsp;</a> 
```python
pull_state_groups(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L401)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L401?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.evaluate_state_terms" class="docs-object-method">&nbsp;</a> 
```python
evaluate_state_terms(self, states, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L415)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L415?message=Update%20Docs)]
</div>
Evaluates terms coming from different state excitations.
Doesn't do any pre-filtering, since that's expected to be in the caller.
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.load_generator" class="docs-object-method">&nbsp;</a> 
```python
load_generator(self, a, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L443)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L443?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.get_paths" class="docs-object-method">&nbsp;</a> 
```python
get_paths(sizes, change): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L474)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L474?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.get_path_poly" class="docs-object-method">&nbsp;</a> 
```python
get_path_poly(path, parities=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L500)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L500?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.poly_coeffs" class="docs-object-method">&nbsp;</a> 
```python
poly_coeffs(self, delta, shift=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L566)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L566?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.get_poly_coeffs" class="docs-object-method">&nbsp;</a> 
```python
get_poly_coeffs(terms, delta, shift=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L581)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L581?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.poly_term_generator" class="docs-object-method">&nbsp;</a> 
```python
poly_term_generator(terms, delta, shift=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L588)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L588?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.rho_term_generator" class="docs-object-method">&nbsp;</a> 
```python
rho_term_generator(a, N, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L644)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L644?message=Update%20Docs)]
</div>
Returns a function to be called on a quantum number to get the coefficient associated with exciting that mode by `a` quanta over
`N` steps w/ phase info coming from where the momenta-terms are.
  - `a`: `Any`
    > 
  - `N`: `Any`
    > 
  - `i_phase`: `Any`
    > 
  - `is_complex`: `Any`
    > 
  - `sel`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorMatrixGenerator.rho" class="docs-object-method">&nbsp;</a> 
```python
rho(phases, paths, ni): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L708)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.py#L708?message=Update%20Docs)]
</div>

  - `phases`: `Any`
    > 
  - `paths`: `Any`
    > 
  - `ni`: `Any`
    > 
  - `:returns`: `_`
    >
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-6b4a30" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-6b4a30"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-6b4a30" markdown="1">
 - [HOElements](#HOElements)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-b039a0" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-b039a0"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-b039a0" markdown="1">
 
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

#### <a name="HOElements">HOElements</a>
```python
    def test_HOElements(self):

        # mat_gen = HarmonicOscillatorMatrixGenerator(['x', 'x', 'x', 'x', 'p', 'p'])
        # raise Exception(
        #     mat_gen._get_poly_coeffs(0)
        # )
        terms = ['p', 'x', 'x', 'p', 'x', 'x']
        mat_gen_old = HarmonicOscillatorMatrixGenerator(terms, mode='rho')
        mat_gen_new = HarmonicOscillatorMatrixGenerator(terms, mode='poly')
        # raise Exception(mat_gen_new(np.array([
        #     [0, 2, 4],
        #     [0, 2, 4]
        # ])))
        n = 5
        rows, cols = np.triu_indices(n)
        states = (
            np.concatenate([rows, cols]),
            np.concatenate([cols, rows])
        )
        old_vals = mat_gen_old.evaluate_state_terms(states)
        new_vals = mat_gen_new.evaluate_state_terms(states)

        old_mat = np.zeros((n, n))
        new_mat = np.zeros((n, n))

        old_mat[np.concatenate([rows, cols]), np.concatenate([cols, rows])] = old_vals
        new_mat[np.concatenate([rows, cols]), np.concatenate([cols, rows])] = new_vals
        print("==="*10)
        print(np.round(old_mat, 6))
        print(np.round(new_mat, 6))
        raise Exception(...)

        self.assertTrue(
            np.allclose(
                mat_gen_old.evaluate_state_terms(states),
                mat_gen_new.evaluate_state_terms(states)
            )
        )

        mat_gen_old = HarmonicOscillatorMatrixGenerator(['p', 'x', 'x', 'p'], mode='rho')
        mat_gen_new = HarmonicOscillatorMatrixGenerator(['p', 'x', 'x', 'p'], mode='poly')
        states = np.array([
            [0, 0, 0, 0],
            [0, 2, 4, 6]
            ])
        raise Exception(
            mat_gen_old.evaluate_state_terms(states),
            mat_gen_new.evaluate_state_terms(states)
        )

        # mat_gen = HarmonicOscillatorMatrixGenerator(['p', 'p', 'x', 'x'])
        # mat_gen = HarmonicOscillatorMatrixGenerator(['p', 'p', 'p', 'p'])
        d = 4
        n = np.arange(4)
        cf = mat_gen._get_poly_coeffs(d)
        if not isinstance(cf, np.ndarray) and cf == 0:
            raise Exception(0)
        else:
            raise Exception(
                cf,
                np.dot(
                    cf,
                    np.power(n[np.newaxis, :], np.arange(len(cf))[:, np.newaxis])
                ) * np.sqrt(np.prod([n+i for i in range(1, abs(d)+1)], axis=0))
            )

        mat_gen = HarmonicOscillatorMatrixGenerator(['p', 'x', 'x', 'p'])
        raise Exception(
            mat_gen._get_poly_coeffs(1)
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorMatrixGenerator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator.py#L328?message=Update%20Docs)   
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