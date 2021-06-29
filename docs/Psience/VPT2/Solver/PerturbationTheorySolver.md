## <a id="Psience.VPT2.Solver.PerturbationTheorySolver">PerturbationTheorySolver</a>
A solver that applies perturbation theory
given a series of corrections and population of states.
Supports degenerate and non-degenerate PT.

### Properties and Methods
```python
PastIndexableTuple: type
use_cached_representations: bool
use_cached_basis: bool
StateSpaceWrapper: type
ProjectionOperatorWrapper: type
ProjectedOperator: type
```
<a id="Psience.VPT2.Solver.PerturbationTheorySolver.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, perturbations, states, coupled_states=None, order=2, total_space=None, flat_total_space=None, state_space_iterations=None, state_space_terms=None, allow_sakurai_degs=False, allow_post_PT_calc=True, modify_degenerate_perturbations=False, gaussian_resonance_handling=False, ignore_odd_order_energies=False, intermediate_normalization=False, zero_element_warning=True, degenerate_states=None, memory_constrained=False, logger=None, verbose=False, parallelizer=None, checkpointer=None): 
```

- `perturbations`: `Iterable[Representation]`
    >No description...
- `states`: `BasisStateSpace`
    >No description...
- `coupled_states`: `BasisMultiStateSpace`
    >No description...
- `order`: `Any`
    >No description...
- `degenerate_states`: `Any`
    >No description...
- `degeneracy_mode`: `Any`
    >No description...
- `logger`: `Any`
    >No description...
- `parallelizer`: `Any`
    >No description...
- `checkpointer`: `Any`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.coupled_states" class="docs-object-method">&nbsp;</a>
```python
@property
coupled_states(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.total_space_dim" class="docs-object-method">&nbsp;</a>
```python
@property
total_space_dim(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.flat_total_space" class="docs-object-method">&nbsp;</a>
```python
@property
flat_total_space(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.total_state_space" class="docs-object-method">&nbsp;</a>
```python
@property
total_state_space(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.representations" class="docs-object-method">&nbsp;</a>
```python
@property
representations(self): 
```

- `:returns`: `Iterable[SparseArray]`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.degenerate_spaces" class="docs-object-method">&nbsp;</a>
```python
@property
degenerate_spaces(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.zero_order_energies" class="docs-object-method">&nbsp;</a>
```python
@property
zero_order_energies(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT" class="docs-object-method">&nbsp;</a>
```python
apply_VPT(self): 
```
Applies perturbation theory to the held basis of states using the
        built representations and degenerate state spaces
- `:returns`: `PerturbationTheoryCorrections`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_VPT_representations" class="docs-object-method">&nbsp;</a>
```python
get_VPT_representations(self): 
```
Gets the sparse representations of the passed perturbation inside the basis of coupled states.
- `:returns`: `Iterable[SparseArray]`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.load_state_spaces" class="docs-object-method">&nbsp;</a>
```python
load_state_spaces(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.load_coupled_spaces" class="docs-object-method">&nbsp;</a>
```python
load_coupled_spaces(self): 
```
Determines which states need to be coupled at which levels of correction
        to handle the PT
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_coupled_space" class="docs-object-method">&nbsp;</a>
```python
get_coupled_space(self, input_state_space, degenerate_space, use_second_deg, allow_PT_degs=True, wavefunction_terms=None, spaces=None): 
```
Applies the VPT equations semi-symbolically, dispatching based on how many
        degeneracies we need to handle
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_nondeg_coupled_space" class="docs-object-method">&nbsp;</a>
```python
get_nondeg_coupled_space(self, input_state_space, degenerate_space=None, spaces=None, wavefunction_terms=None): 
```
Applies the non-degenerate equations in semi-symbolic form to determine
        which states needs to be calculated.
        This will always be the initial input to a calculation and then
        certain additional states may be calculated on the fly if they are needed to handle
        truly degenerate stuff.
        The only difference there will be to add on
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_deg_coupled_space" class="docs-object-method">&nbsp;</a>
```python
get_deg_coupled_space(self, degenerate_space, spaces=None): 
```
Applies the degenerate equations in semi-symbolic form to determine
        which states needs to be calculated at which orders of correction
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_second_deg_coupled_space" class="docs-object-method">&nbsp;</a>
```python
get_second_deg_coupled_space(self, degenerate_space, spaces=None): 
```
Does the dirty work of doing the VPT iterative equations.
        Needs to be adapted to include the two types of degeneracies that can
        be introduced in Sakurai's approach.
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_corrections" class="docs-object-method">&nbsp;</a>
```python
get_corrections(self, non_zero_cutoff=1e-14): 
```
Applies the perturbation theory equations to obtain
        corrections to the wave functions and energies
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_equations" class="docs-object-method">&nbsp;</a>
```python
apply_VPT_equations(self, state_index, degenerate_space_indices, degenerate_energies, zero_order_state, degenerate_subspace, degenerate_subsubspace, perturbations=None, allow_PT_degs=None, ignore_odd_orders=None, intermediate_normalization=None, non_zero_cutoff=1e-14): 
```
Applies VPT equations, dispatching based on how many
        degeneracies we need to handle
- `state_index`: `int`
    >the index of the primary state being treated using the PT
- `degenerate_space_indices`: `np.ndarray[int]`
    >the indices corresponding to degeneracies with the primary state in the zero-order picture
- `degenerate_energies`: `Iterable[float | None]`
    >the first and (possibly) second order correction to the energies
- `zero_order_states`: `np.ndarray[float]`
    >the vector for the proper zero-order state corresponding ot state_index
- `degenerate_subsubspace`: `tuple[np.ndarray[float], np.ndarray[int]]`
    >the set of vectors for the zero-order states in the secondary degenerate subspace
- `non_zero_cutoff`: `float`
    >cutoff for when a term can be called zero for performance reasons
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_nondeg_equations" class="docs-object-method">&nbsp;</a>
```python
apply_VPT_nondeg_equations(self, state_index, deg_group, perturbations=None, non_zero_cutoff=1e-14, check_overlap=True, intermediate_normalization=False, ignore_odd_orders=False): 
```
Does the dirty work of doing the VPT iterative equations.
        Needs to be adapted to include the two types of degeneracies that can
        be introduced in Sakurai's approach.
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_deg_equations" class="docs-object-method">&nbsp;</a>
```python
apply_VPT_deg_equations(self, state_index, degenerate_space_indices, degenerate_energy, zero_order_state, degenerate_subspace, non_zero_cutoff=1e-14, check_overlap=True): 
```
Does the dirty work of doing the VPT iterative equations.
        Needs to be adapted to include the two types of degeneracies that can
        be introduced in Sakurai's approach.
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_second_deg_equations" class="docs-object-method">&nbsp;</a>
```python
apply_VPT_second_deg_equations(self, state_index, degenerate_space_indices, degenerate_energies, zero_order_state, degenerate_subspace, degenerate_subsubspace, non_zero_cutoff=1e-14, check_overlap=True): 
```
Does the dirty work of doing the VPT iterative equations.
        Needs to be adapted to include the two types of degeneracies that can
        be introduced in Sakurai's approach.
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_post_PT_variational_calc" class="docs-object-method">&nbsp;</a>
```python
apply_post_PT_variational_calc(self, degenerate_states, corrs): 
```
Applies degenerate perturbation theory by building a representation
        for the degenerate terms in the Hamiltonian.
        This is then diagonalized, allowing the degenerate states to be expressed
        in the basis of non-degenerate states
- `H`: `Iterable[SparseArray]`
    >No description...
- `corrs`: `PerturbationTheoryCorrections`
    >the standard PerturbationTheory Corrections object that comes out of the application of non-deg PT
- `degenerate_states`: `Any`
    >population of degenerate states
- `logger`: `Logger`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.drop_deg_pert_els" class="docs-object-method">&nbsp;</a>
```python
drop_deg_pert_els(self, perts, deg_groups): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_transformed_Hamiltonians" class="docs-object-method">&nbsp;</a>
```python
get_transformed_Hamiltonians(self, corrs, deg_group=None): 
```

<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_degenerate_rotation" class="docs-object-method">&nbsp;</a>
```python
get_degenerate_rotation(self, deg_group, corrs): 
```

### Examples


