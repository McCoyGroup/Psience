## <a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian">PerturbationTheoryHamiltonian</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian.py#L27)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian.py#L27?message=Update%20Docs)]
</div>

Represents the main Hamiltonian used in the perturbation theory calculation.
Uses a harmonic oscillator basis for representing H0, H1, H2 etc.
The PT process is split into a PerturbationTheorySolver and a PerturbationTheoryHamiltonian
where the Hamiltonian just implements the split of the perturbation and the Solver manages the equations.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
TermGetter: TermGetter
CoriolisTermGetter: CoriolisTermGetter
```
<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule=None, n_quanta=None, modes=None, mode_selection=None, potential_derivatives=None, include_potential=True, include_gmatrix=True, include_coriolis_coupling=True, include_pseudopotential=True, include_only_mode_couplings=None, potential_terms=None, allow_higher_potential_terms=False, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, selection_rules=None, operator_chunk_size=None, operator_coefficient_threshold=None, matrix_element_threshold=None, logger=None, checkpoint=None, results=None, parallelizer=None, **expansion_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L35)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L35?message=Update%20Docs)]
</div>

  - `molecule`: `Molecule`
    > the molecule on which we're doing perturbation theory
  - `n_quanta`: `int | None`
    > the numbers of quanta to use when representing the entire state space
  - `modes`: `None | MolecularNormalModes`
    > the set of modes to use as the basis
  - `mode_selection`: `None | Iterable[int]`
    > the subset of modes to use when doing expansions
  - `include_coriolis_coupling`: `bool`
    > whether to add coriolis coupling if not in internals
  - `parallelizer`: `Parallelizer`
    > parallelism manager
  - `logger`: `str | Logger`
    > log file or logger to write to
  - `checkpoint`: `str | Checkpointer`
    > checkpoint file or checkpointer to store intermediate results


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.from_fchk" class="docs-object-method">&nbsp;</a> 
```python
from_fchk(file, internals=None, mode_selection=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L185)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L185?message=Update%20Docs)]
</div>

  - `file`: `str`
    > fchk file to load from
  - `internals`: `Iterable[Iterable[int]]`
    > internal coordinate specification as a Z-matrix ordering
  - `n_quanta`: `int | Iterable[int]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.H0" class="docs-object-method">&nbsp;</a> 
```python
@property
H0(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L251)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L251?message=Update%20Docs)]
</div>
Provides the representation for H0 in this basis


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_perturbations" class="docs-object-method">&nbsp;</a> 
```python
get_perturbations(self, expansion_orders, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L424)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L424?message=Update%20Docs)]
</div>
Gets the `Representation` objects for the perturbations up through second order
  - `order`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_Nielsen_xmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_Nielsen_xmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L589)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L589?message=Update%20Docs)]
</div>
Provides Nielsen's X-Matrix when working in Cartesian coordinates
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_Nielsen_energies" class="docs-object-method">&nbsp;</a> 
```python
get_Nielsen_energies(self, states, x_mat=None, return_split=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L615)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L615?message=Update%20Docs)]
</div>

  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_coupled_space" class="docs-object-method">&nbsp;</a> 
```python
get_coupled_space(self, states, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L648)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L648?message=Update%20Docs)]
</div>
Returns the set of states that couple the given states up to the given order at each level of perturbation (beyond zero order).
We keep track of how each individual state in states is transformed, as we only need to compute elements within those
blocks, allowing for relatively dramatic speed-ups.
  - `state`: `BasisStateSpace`
    > the states of interest
  - `order`: `int`
    > the order of perturbation theory we're doing
  - `freqs`: `Iterable[float]`
    > the zero-order frequencies in each vibrational mode being coupled
  - `freq_threshold`: `None | float`
    > the threshold for the maximum frequency difference between states to be considered
  - `:returns`: `tuple[BasisMultiStateSpace]`
    > t
h
e
 
s
e
t
s
 
o
f
 
c
o
u
p
l
e
d
 
s
t
a
t
e
s


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_representations" class="docs-object-method">&nbsp;</a> 
```python
get_representations(self, states, coupled_states=None, degeneracies=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L694)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L694?message=Update%20Docs)]
</div>
Returns the representations of the perturbations in the basis of coupled states
  - `coupled_states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_input_state_spaces" class="docs-object-method">&nbsp;</a> 
```python
get_input_state_spaces(self, states, coupled_states=None, degeneracies=None, order=2, deg_extra_order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L729)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L729?message=Update%20Docs)]
</div>
Converts the input state specs into proper `BasisStateSpace` specs that
will directly feed into the code
  - `states`: `Any`
    > 
  - `coupled_states`: `Any`
    > 
  - `degeneracies`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_solver" class="docs-object-method">&nbsp;</a> 
```python
get_solver(self, states, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L847)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L847?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunctions(self, states, initial_states=None, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, results=None, degenerate_transformation_layout=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L893)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L893?message=Update%20Docs)]
</div>
Gets a set of `PerturbationTheoryWavefunctions` from the perturbations defined by the Hamiltonian
  - `states`: `BasisStateSpace | Iterable[int] | Iterable[Iterable[int]]`
    > the states to get the index for, given either as indices or as a numbers of quanta
  - `coupled_states`: `BasisStateSpace | Iterable[int] | Iterable[Iterable[int]]`
    > the list of states to explicitly allow to couple in
  - `degeneracies`: `(Iterable[int], Iterable[int])  | (Iterable[Iterable[int]], Iterable[Iterable[int]])`
    > the pairs of states to be treated via degenerate perturbation theory
  - `:returns`: `PerturbationTheoryWavefunctions`
    > g
e
n
e
r
a
t
e
d
 
w
a
v
e
 
f
u
n
c
t
i
o
n
s


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_action_expansion" class="docs-object-method">&nbsp;</a> 
```python
get_action_expansion(self, coupled_states=None, degeneracies=None, allow_sakurai_degs=False, allow_post_PT_calc=True, modify_degenerate_perturbations=False, intermediate_normalization=False, ignore_odd_order_energies=True, zero_element_warning=True, state_space_iterations=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L1060)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L1060?message=Update%20Docs)]
</div>
Gets the expansion of the energies in terms of Miller's "classical actions" by
doing just enough PT to invert the matrix
  - `order`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_breakdown" class="docs-object-method">&nbsp;</a> 
```python
get_breakdown(self, states, coupled_states=None, degeneracies=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L1105)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.py#L1105?message=Update%20Docs)]
</div>
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Hamiltonian.py#L27?message=Update%20Docs)   
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