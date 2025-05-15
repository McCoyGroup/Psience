## <a id="Psience.VPT2.Runner.VPTRunner">VPTRunner</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner.py#L1145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner.py#L1145?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTRunner.py#L1151)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTRunner.py#L1151?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTRunner.py#L1196)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTRunner.py#L1196?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTRunner.hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
@property
hamiltonian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTRunner.py#L1203)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTRunner.py#L1203?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTRunner.get_wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunctions(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTRunner.py#L1209)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTRunner.py#L1209?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTRunner.get_solver" class="docs-object-method">&nbsp;</a> 
```python
get_solver(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTRunner.py#L1227)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTRunner.py#L1227?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTRunner.print_output_tables" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
print_output_tables(cls, wfns=None, file=None, print_intensities=True, print_energies=True, print_energy_corrections=True, print_transition_moments=True, operators=None, logger=None, sep_char='=', sep_len=100): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1238)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1238?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTRunner.py#L1315)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTRunner.py#L1315?message=Update%20Docs)]
</div>
Prints a bunch of formatted output data from a PT run
  - `wfns`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTRunner.construct" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct(cls, system, states, target_property=None, extended_space_target_property=None, basis_filters=None, initial_states=None, corrected_fundamental_frequencies=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1412)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1412?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTRunner.run_simple" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
run_simple(cls, system, states, target_property=None, corrected_fundamental_frequencies=None, calculate_intensities=True, plot_spectrum=False, operators=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1532)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1532?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner.py#L1145?message=Update%20Docs)   
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