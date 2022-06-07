## <a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian">PerturbationTheoryHamiltonian</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L27)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L27?message=Update%20Docs)]
</div>

Represents the main Hamiltonian used in the perturbation theory calculation.
Uses a harmonic oscillator basis for representing H0, H1, H2 etc.
The PT process is split into a PerturbationTheorySolver and a PerturbationTheoryHamiltonian
where the Hamiltonian just implements the split of the perturbation and the Solver manages the equations.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule=None, n_quanta=None, modes=None, mode_selection=None, potential_derivatives=None, include_coriolis_coupling=True, include_pseudopotential=True, potential_terms=None, allow_higher_potential_terms=False, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, selection_rules=None, operator_chunk_size=None, logger=None, checkpoint=None, results=None, parallelizer=None, **expansion_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L35)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L35?message=Update%20Docs)]
</div>


- `molecule`: `Molecule`
    >the molecule on which we're doing perturbation theory
- `n_quanta`: `int | None`
    >the numbers of quanta to use when representing the entire state space
- `modes`: `None | MolecularNormalModes`
    >the set of modes to use as the basis
- `mode_selection`: `None | Iterable[int]`
    >the subset of modes to use when doing expansions
- `include_coriolis_coupling`: `bool`
    >whether to add coriolis coupling if not in internals
- `parallelizer`: `Parallelizer`
    >parallelism manager
- `logger`: `str | Logger`
    >log file or logger to write to
- `checkpoint`: `str | Checkpointer`
    >checkpoint file or checkpointer to store intermediate results

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.from_fchk" class="docs-object-method">&nbsp;</a> 
```python
from_fchk(file, internals=None, mode_selection=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L168)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L168?message=Update%20Docs)]
</div>


- `file`: `str`
    >fchk file to load from
- `internals`: `Iterable[Iterable[int]]`
    >internal coordinate specification as a Z-matrix ordering
- `n_quanta`: `int | Iterable[int]`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.H0" class="docs-object-method">&nbsp;</a> 
```python
@property
H0(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L?message=Update%20Docs)]
</div>

Provides the representation for H0 in this basis

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.H1" class="docs-object-method">&nbsp;</a> 
```python
@property
H1(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L?message=Update%20Docs)]
</div>

Provides the representation for H1 in this basis

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.H2" class="docs-object-method">&nbsp;</a> 
```python
@property
H2(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L?message=Update%20Docs)]
</div>

Provides the representation for H2 in this basis

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_perturbations" class="docs-object-method">&nbsp;</a> 
```python
get_perturbations(self, expansion_orders): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L460)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L460?message=Update%20Docs)]
</div>

Gets the `Representation` objects for the perturbations up through second order
- `order`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_Nielsen_xmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_Nielsen_xmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L623)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L623?message=Update%20Docs)]
</div>

Provides Nielsen's X-Matrix when working in Cartesian coordinates
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_Nielsen_energies" class="docs-object-method">&nbsp;</a> 
```python
get_Nielsen_energies(self, states, x_mat=None, return_split=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L649)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L649?message=Update%20Docs)]
</div>


- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_coupled_space" class="docs-object-method">&nbsp;</a> 
```python
get_coupled_space(self, states, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L682)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L682?message=Update%20Docs)]
</div>

Returns the set of states that couple the given states up to the given order at each level of perturbation (beyond zero order).
        We keep track of how each individual state in states is transformed, as we only need to compute elements within those
        blocks, allowing for relatively dramatic speed-ups.
- `state`: `BasisStateSpace`
    >the states of interest
- `order`: `int`
    >the order of perturbation theory we're doing
- `freqs`: `Iterable[float]`
    >the zero-order frequencies in each vibrational mode being coupled
- `freq_threshold`: `None | float`
    >the threshold for the maximum frequency difference between states to be considered
- `:returns`: `tuple[BasisMultiStateSpace]`
    >the sets of coupled states

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_representations" class="docs-object-method">&nbsp;</a> 
```python
get_representations(self, states, coupled_states=None, degeneracies=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L728)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L728?message=Update%20Docs)]
</div>

Returns the representations of the perturbations in the basis of coupled states
- `coupled_states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_input_state_spaces" class="docs-object-method">&nbsp;</a> 
```python
get_input_state_spaces(self, states, coupled_states=None, degeneracies=None, order=2, deg_extra_order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L763)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L763?message=Update%20Docs)]
</div>

Converts the input state specs into proper `BasisStateSpace` specs that
        will directly feed into the code
- `states`: `Any`
    >No description...
- `coupled_states`: `Any`
    >No description...
- `degeneracies`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_solver" class="docs-object-method">&nbsp;</a> 
```python
get_solver(self, states, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L881)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L881?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunctions(self, states, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, results=None, degenerate_transformation_layout=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L924)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L924?message=Update%20Docs)]
</div>

Gets a set of `PerturbationTheoryWavefunctions` from the perturbations defined by the Hamiltonian
- `states`: `BasisStateSpace | Iterable[int] | Iterable[Iterable[int]]`
    >the states to get the index for, given either as indices or as a numbers of quanta
- `coupled_states`: `BasisStateSpace | Iterable[int] | Iterable[Iterable[int]]`
    >the list of states to explicitly allow to couple in
- `degeneracies`: `(Iterable[int], Iterable[int])  | (Iterable[Iterable[int]], Iterable[Iterable[int]])`
    >the pairs of states to be treated via degenerate perturbation theory
- `:returns`: `PerturbationTheoryWavefunctions`
    >generated wave functions

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_action_expansion" class="docs-object-method">&nbsp;</a> 
```python
get_action_expansion(self, coupled_states=None, degeneracies=None, allow_sakurai_degs=False, allow_post_PT_calc=True, modify_degenerate_perturbations=False, intermediate_normalization=False, ignore_odd_order_energies=True, zero_element_warning=True, state_space_iterations=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L1088)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L1088?message=Update%20Docs)]
</div>

Gets the expansion of the energies in terms of Miller's "classical actions" by
        doing just enough PT to invert the matrix
- `order`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_breakdown" class="docs-object-method">&nbsp;</a> 
```python
get_breakdown(self, states, coupled_states=None, degeneracies=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Hamiltonian.py#L1133)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L1133?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Hamiltonian.py#L27?message=Update%20Docs)