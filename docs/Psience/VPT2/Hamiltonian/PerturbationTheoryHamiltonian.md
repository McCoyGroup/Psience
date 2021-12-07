## <a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian">PerturbationTheoryHamiltonian</a>
Represents the main Hamiltonian used in the perturbation theory calculation.
Uses a harmonic oscillator basis for representing H0, H1, H2 etc.
The PT process is split into a PerturbationTheorySolver and a PerturbationTheoryHamiltonian
where the Hamiltonian just implements the split of the perturbation and the Solver manages the equations.

### Properties and Methods
<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule=None, n_quanta=None, modes=None, mode_selection=None, potential_derivatives=None, include_coriolis_coupling=True, include_pseudopotential=True, potential_terms=None, allow_higher_potential_terms=False, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, selection_rules=None, operator_chunk_size=None, logger=None, checkpoint=None, parallelizer=None, **expansion_options): 
```

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
Provides the representation for H0 in this basis

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.H1" class="docs-object-method">&nbsp;</a>
```python
@property
H1(self): 
```
Provides the representation for H1 in this basis

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.H2" class="docs-object-method">&nbsp;</a>
```python
@property
H2(self): 
```
Provides the representation for H2 in this basis

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_perturbations" class="docs-object-method">&nbsp;</a>
```python
get_perturbations(self, order): 
```

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.perturbations" class="docs-object-method">&nbsp;</a>
```python
@property
perturbations(self): 
```

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_Nielsen_xmatrix" class="docs-object-method">&nbsp;</a>
```python
get_Nielsen_xmatrix(self): 
```

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_Nielsen_energies" class="docs-object-method">&nbsp;</a>
```python
get_Nielsen_energies(self, states, x_mat=None, return_split=False): 
```

- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_coupled_space" class="docs-object-method">&nbsp;</a>
```python
get_coupled_space(self, states, order): 
```
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
Returns the representations of the perturbations in the basis of coupled states
- `coupled_states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_input_state_spaces" class="docs-object-method">&nbsp;</a>
```python
get_input_state_spaces(self, states, coupled_states=None, degeneracies=None, order=2, deg_extra_order=2): 
```
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

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryHamiltonian.get_wavefunctions" class="docs-object-method">&nbsp;</a>
```python
get_wavefunctions(self, states, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, results=None, degenerate_transformation_layout=None, **opts): 
```
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

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Hamiltonian.py?message=Update%20Docs)