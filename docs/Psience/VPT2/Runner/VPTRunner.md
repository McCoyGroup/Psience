## <a id="Psience.VPT2.Runner.VPTRunner">VPTRunner</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1283)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1283?message=Update%20Docs)]
</div>

A helper class to make it easier to run jobs by making the inputs/options
clear and making it easier to customize run options







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
helpers: AnneInputHelpers
```
<a id="Psience.VPT2.Runner.VPTRunner.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, system, states, initial_states=None, hamiltonian_options=None, solver_options=None, runtime_options=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1289)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1289?message=Update%20Docs)]
</div>

  - `system`: `VPTSystem`
    > the system to run perturbation theory on
  - `hamiltonian_options`: `VPTHamiltonianOptions`
    > options to configure the Hamiltonian
  - `solver_options`: `VPTSolverOptions`
    > options to configure the way the perturbation theory is applied
  - `runtime_options`: `VPTRuntimeOptions`
    > options to configure the way the code runs


<a id="Psience.VPT2.Runner.VPTRunner.get_Hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
get_Hamiltonian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTRunner.py#L1334)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTRunner.py#L1334?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `PerturbationTheoryHamiltonian` for this runner's system, combining the configured Hamiltonian and runtime options.
  - `:returns`: `PerturbationTheoryHamiltonian`
    > the constructed Hamiltonian


<a id="Psience.VPT2.Runner.VPTRunner.hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
@property
hamiltonian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTRunner.py#L1349)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTRunner.py#L1349?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) `PerturbationTheoryHamiltonian` for this runner, built lazily via `get_Hamiltonian` the first time it's needed.
  - `:returns`: `PerturbationTheoryHamiltonian`
    > the Hamiltonian


<a id="Psience.VPT2.Runner.VPTRunner.get_wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunctions(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTRunner.py#L1363)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTRunner.py#L1363?message=Update%20Docs)]
</div>
**LLM Docstring**

Run the full VPT calculation and return the resulting wavefunctions, combining the solver options, degenerate-state/degeneracy-handler settings, and runtime options configured on this runner (with any explicitly passed `opts` taking precedence).
  - `opts`: `dict`
    > extra options overriding the runner's configured solver/runtime options
  - `:returns`: `PerturbationTheoryWavefunctions`
    > the computed perturbation-theory wavefunctions


<a id="Psience.VPT2.Runner.VPTRunner.get_solver" class="docs-object-method">&nbsp;</a> 
```python
get_solver(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTRunner.py#L1391)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTRunner.py#L1391?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `PerturbationTheorySolver` for this runner's target states, without running the full wavefunction calculation, combining the configured solver and runtime options.
  - `:returns`: `PerturbationTheorySolver`
    > the constructed solver


<a id="Psience.VPT2.Runner.VPTRunner.print_output_tables" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
print_output_tables(cls, wfns=None, file=None, print_intensities=True, print_energies=True, print_energy_corrections=True, print_transition_moments=True, operators=None, logger=None, sep_char='=', sep_len=100): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1410)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1410?message=Update%20Docs)]
</div>
Prints a bunch of formatted output data from a PT run
  - `wfns`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTRunner.print_tables" class="docs-object-method">&nbsp;</a> 
```python
print_tables(self, wfns=None, file=None, print_intensities=True, print_energy_corrections=True, print_transition_moments=True, operators=None, logger=None, sep_char='=', sep_len=100): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTRunner.py#L1543)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTRunner.py#L1543?message=Update%20Docs)]
</div>
Prints a bunch of formatted output data from a PT run
  - `wfns`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTRunner.get_Nielsen_energies" class="docs-object-method">&nbsp;</a> 
```python
get_Nielsen_energies(self, return_split=False, return_X=False, **potential_params): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTRunner.py#L1582)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTRunner.py#L1582?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute harmonic and anharmonic (Nielsen-formula) vibrational energies for the target states, via the Hamiltonian's own `get_Nielsen_energies`.
  - `return_split`: `bool`
    > whether to also return the anharmonicity split out separately
  - `return_X`: `bool`
    > whether to also return the underlying X anharmonicity-constant matrix
  - `potential_params`: `dict`
    > extra options forwarded to `PerturbationTheoryHamiltonian.get_Nielsen_energies`
  - `:returns`: `tuple`
    > `(harmonic, total)` energies, or `(harmonic, total, x_matrix)` if `return_split`/`return_X` is set


<a id="Psience.VPT2.Runner.VPTRunner.print_Nielsen_frequencies" class="docs-object-method">&nbsp;</a> 
```python
print_Nielsen_frequencies(self, logger=None, state_formatting='vector', **potential_params): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTRunner.py#L1609)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTRunner.py#L1609?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Nielsen-formula harmonic and anharmonic transition frequencies (relative to the ground state) and print them as a formatted state/frequency table.
  - `logger`: `Logger | bool | None`
    > `None` to print directly, `True` to use the Hamiltonian's own logger, or an explicit logger to print into a block
  - `state_formatting`: `str`
    > `'vector'` to display states as raw excitation vectors, or another format understood by `VPTStateMaker.parse_state`
  - `potential_params`: `dict`
    > extra options forwarded to `get_Nielsen_energies`
  - `:returns`: `None`
    > None


<a id="Psience.VPT2.Runner.VPTRunner.construct" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct(cls, system, states, target_property=None, extended_space_target_property=None, basis_filters=None, initial_states=None, corrected_fundamental_frequencies=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1737)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1737?message=Update%20Docs)]
</div>
**LLM Docstring**

Top-level constructor that assembles a fully configured `VPTRunner` (and, if target states are given, its resolved `VPTStateSpace`) from a molecule/system specification, target states, and a large set of options split automatically across `VPTSystem`, `VPTStateSpace`, `VPTHamiltonianOptions`, `VPTRuntimeOptions`, and `VPTSolverOptions`.
  - `system`: `str | list | Molecule | VPTSystem`
    > the molecule or system specification to run VPT on
  - `states`: `VPTStateSpace | int | Iterable`
    > the target states, given as a states object, quanta cutoff, or explicit list
  - `target_property`: `str | None`
    > the property (e.g. `'intensities'`) used to build a default state-space filter if `state_space_filters` isn't explicitly given
  - `extended_space_target_property`: `object | None`
    > accepted for interface consistency but not directly used in this method's body
  - `basis_filters`: `object | None`
    > extra post-transformation filters merged into the generated state-space filter
  - `initial_states`: `VPTStateSpace | int | Iterable | None`
    > the initial states used both for filter construction and as the runner's own initial states
  - `corrected_fundamental_frequencies`: `np.ndarray | None`
    > replacement fundamental frequencies to use for the underlying molecule's normal modes
  - `opts`: `dict`
    > the full set of system/state-space/Hamiltonian/runtime/solver options, validated against the union of all their `__props__`
  - `:returns`: `VPTRunner | tuple[VPTRunner, VPTStateSpace]`
    > the constructed runner, or `(runner, states)` if target states were given


<a id="Psience.VPT2.Runner.VPTRunner.run_simple" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
run_simple(cls, system, states, target_property=None, corrected_fundamental_frequencies=None, calculate_intensities=True, plot_spectrum=False, operators=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1882)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1882?message=Update%20Docs)]
</div>
The standard runner for VPT.
Makes a runner using the `construct` method and then calls that
runner's `print_tables` method after printing out run info.
  - `system`: `list|str|Molecule`
    > the system spec, either as a `Molecule`, molecule spec (atoms, coords, opts) or a file to construct a `Molecule`
  - `states`: `int|list`
    > the states to get corrections for either an `int` (up to that many quanta) or an explicit state list
  - `target_property`: `str`
    > the target property to get corrections for (one of 'frequencies', 'intensities', 'wavefunctions')
  - `corrected_fundamental_frequencies`: `Iterable[float]|None`
    > a set of fundamental frequencies to use to get new zero-order energies
  - `calculate_intensities`: `bool default:True`
    > whether or not to calculate energies
  - `opts`: `Any`
    > options that work for a `VPTSystem`, `VPTStateSpace`, `VPTRuntimeOptions`, `VPTSolverOptions`, or `VPTHamiltonianOptions` object which will be filtered automatically
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTRunner.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTRunner.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTRunner.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTRunner.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1283?message=Update%20Docs)   
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