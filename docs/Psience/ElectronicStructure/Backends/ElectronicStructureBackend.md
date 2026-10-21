## <a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend">ElectronicStructureBackend</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends.py#L111)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends.py#L111?message=Update%20Docs)]
</div>

One electronic-structure model evaluated on Cartesian geometries of a fixed set of atoms.

`evaluate(coords, order, properties)` returns, for each requested property, the list
`[value, d value/dx, d^2 value/dx^2, ...]` up to `order`, where `x` are the flattened
Cartesians (so the `k`-th term has shape `base_shape + (3 * natoms,) * k + property_shape`).
Coordinates are in `distance_units`; property values are in `property_units[prop]`.
Derivatives beyond `analytic_derivative_orders[prop]` are not available from a backend
(the `Molecools` evaluators fill those in by finite differences).

The total `charge` and spin `multiplicity` set at construction can be overridden per call.

Subclasses implement `evaluate_properties` on a flat batch of structures.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
distance_units: str
property_units: dict
analytic_derivative_orders: dict
property_shapes: dict
default_charge: int
default_multiplicity: NoneType
backend_registry: dict
```
<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, charge=None, multiplicity=None, quiet=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends.py#L135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends.py#L135?message=Update%20Docs)]
</div>

  - `atoms`: `Sequence[str]`
    > atom symbols
  - `charge`: `int | None`
    > total charge (default `default_charge`)
  - `multiplicity`: `int | None`
    > spin multiplicity (default `default_multiplicity`)
  - `quiet`: `bool`
    > silence the underlying program's output


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.register" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
register(cls, name, backend=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L154)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L154?message=Update%20Docs)]
</div>
Register a backend class under `name` (usable as a decorator).


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.get_backends" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_backends(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L166)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L166?message=Update%20Docs)]
</div>
The registered backends, by name.
  - `:returns`: `dict[str, type]`
    >


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.handle_specialization" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_specialization(cls, tag): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L175)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L175?message=Update%20Docs)]
</div>
Options implied by a `name:tag` specialization (e.g. `'mace:small'`).


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.resolve_backend" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
resolve_backend(cls, spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L182)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L182?message=Update%20Docs)]
</div>
Resolve `'name'` or `'name:tag'` to `(backend_class, options)`.
  - `:returns`: `(type, dict)`
    >


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.from_spec" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_spec(cls, spec, atoms, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L200)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L200?message=Update%20Docs)]
</div>
Build a backend from a name/`'name:tag'` (or return `spec` if it already is one).
  - `:returns`: `ElectronicStructureBackend`
    >


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.discover_dependencies" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
discover_dependencies(cls, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L212)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L212?message=Update%20Docs)]
</div>
The packages/programs this backend needs, as a `McUtils.ExternalPrograms.JobDependencies`.


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.get_init_options" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_init_options(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L220)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L220?message=Update%20Docs)]
</div>
The names this backend's constructor accepts and whether it also takes arbitrary
keywords (which are then passed on to the underlying model).
  - `:returns`: `(set[str], bool)`
    >


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.natoms" class="docs-object-method">&nbsp;</a> 
```python
@property
natoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L247)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L247?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.supported_properties" class="docs-object-method">&nbsp;</a> 
```python
@property
supported_properties(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L251)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L251?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.get_analytic_derivative_order" class="docs-object-method">&nbsp;</a> 
```python
get_analytic_derivative_order(self, prop): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L255)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L255?message=Update%20Docs)]
</div>
The highest coordinate-derivative order available analytically for `prop`.


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.get_property_units" class="docs-object-method">&nbsp;</a> 
```python
get_property_units(self, prop): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L264)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L264?message=Update%20Docs)]
</div>
The native units of `prop` (a unit name, a tuple of unit names, or `None` for model-native/atomic).


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.get_property_shape" class="docs-object-method">&nbsp;</a> 
```python
get_property_shape(self, prop): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L270)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L270?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.resolve_state" class="docs-object-method">&nbsp;</a> 
```python
resolve_state(self, charge=None, multiplicity=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L276)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L276?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.check_electron_count" class="docs-object-method">&nbsp;</a> 
```python
check_electron_count(self, charge, multiplicity): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L282)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L282?message=Update%20Docs)]
</div>
Raise if `charge` and `multiplicity` are inconsistent with the electron count.


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.prep_coords" class="docs-object-method">&nbsp;</a> 
```python
prep_coords(self, coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L297)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L297?message=Update%20Docs)]
</div>
Flatten `(..., natoms, 3)` or `(..., 3 natoms)` coordinates to `(nstruct, natoms, 3)`.
  - `:returns`: `_`
    > `(base_shape, flat_coords)`


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.canonicalize_request" class="docs-object-method">&nbsp;</a> 
```python
canonicalize_request(self, order, properties): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L313)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L313?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.evaluate_properties" class="docs-object-method">&nbsp;</a> 
```python
evaluate_properties(self, coords, orders, charge, multiplicity, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L332)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L332?message=Update%20Docs)]
</div>
Evaluate properties on a flat batch of structures.
  - `coords`: `np.ndarray`
    > `(nstruct, natoms, 3)`, in `distance_units`
  - `orders`: `dict[str, int]`
    > the derivative order wanted for each property
  - `charge`: `Any`
    > the total charge
  - `multiplicity`: `Any`
    > the spin multiplicity
  - `:returns`: `dict[str, list[np.ndarray]]`
    > for each property, `[value, d1, ...]` with a leading `nstruct` axis and
    derivative axes of size `3 * natoms`


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, coords, order=0, properties=None, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L349)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L349?message=Update%20Docs)]
</div>
Evaluate properties and their analytic Cartesian derivatives.
  - `coords`: `np.ndarray`
    > `(..., natoms, 3)` or `(..., 3 natoms)`, in `distance_units`
  - `order`: `int | dict[str, int]`
    > the derivative order (for every property), or `{property: order}`
  - `properties`: `str | Iterable[str] | None`
    > property names (default: the backend's first supported property)
  - `charge`: `Any`
    > total charge override for this call
  - `multiplicity`: `Any`
    > spin multiplicity override for this call
  - `opts`: `Any`
    > backend-specific evaluation options
  - `:returns`: `dict[str, list[np.ndarray]]`
    > `{property: [value, d1, ...]}`, with the leading batch shape of `coords`


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.evaluate_property" class="docs-object-method">&nbsp;</a> 
```python
evaluate_property(self, coords, prop, order=0, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L374)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L374?message=Update%20Docs)]
</div>
`evaluate` for a single property.
  - `:returns`: `list[np.ndarray]`
    >


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.evaluate_field_response" class="docs-object-method">&nbsp;</a> 
```python
evaluate_field_response(self, coords, coord_order=0, field_order=2, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L382)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L382?message=Update%20Docs)]
</div>
The response to a uniform electric field, `[dipole_expansion, polarizability_expansion][:field_order]`,
each a list over coordinate-derivative orders through `coord_order`. Backends with their own
response models (charge equilibration) override this.
  - `:returns`: `list[list[np.ndarray]]`
    >


<a id="Psience.ElectronicStructure.Backends.ElectronicStructureBackend.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L396)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.py#L396?message=Update%20Docs)]
</div>
 </div>
</div>




## Examples


















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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/Backends/ElectronicStructureBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/Backends.py#L111?message=Update%20Docs)   
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