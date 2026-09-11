## <a id="Psience.Molecools.Builder.MoleculeBuilder">MoleculeBuilder</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Builder.py#L36)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Builder.py#L36?message=Update%20Docs)]
</div>

Namespace of classmethods for constructing molecules from functional groups or fragments.

Not instantiated: all functionality is exposed through classmethods, so no constructor is
needed. The scaffold molecule (formerly `self`) is passed in explicitly, and the molecule
class used to build new objects is passed as `molecule_type` (defaulting to the standard
`Molecule`).







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Molecools.Builder.MoleculeBuilder.fragment_embedding" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
fragment_embedding(cls, mol, fragment_indices, ref=None, return_axes=False, order=1, view_inds=(1, 2), excluded=None, use_moments=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L59)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L59?message=Update%20Docs)]
</div>
Compute a local coordinate frame (origin, offset vector, and an up-vector or full axis
set) anchored at a fragment of `mol`, used as the reference frame for attaching or
orienting substituents; falls back to center-of-mass/principal-axis reference points
(encoded as indices `-1`/`-2`/`-3`) when the fragment doesn't have enough atoms of its
own to define a frame.
  - `mol`: `Molecule`
    > the scaffold molecule the fragment lives in
  - `fragment_indices`: `int | Iterable[int]`
    > the atom index (or indices) defining the fragment to embed
  - `ref`: `Iterable[int] | None`
    > reference atom(s) (outside the fragment) used to anchor the origin/frame;
    computed from the local neighborhood if not given
  - `return_axes`: `bool`
    > whether to return a full 3x3 axis frame instead of just an up-vector
  - `view_inds`: `tuple[int, int]`
    > which two fragment-atom positions define the "view" direction used to
    build the axis frame
  - `use_moments`: `bool`
    > whether to derive the up-vector from the fragment's moments of
    inertia rather than from its first three atom positions
  - `:returns`: `tuple`
    > `(origin, offset, up_or_axes)` -- the reference origin point, the offset from
    origin to the fragment's first atom, and either an up-vector or (if `return_axes`) a
    full axis frame


<a id="Psience.Molecools.Builder.MoleculeBuilder.resolve_stereo_hydrogen" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
resolve_stereo_hydrogen(cls, atoms, coords, bonds, stereo_pos, stereos, ref_exclude=()): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L218)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L218?message=Update%20Docs)]
</div>
atoms, coords, bonds : the geometry the stereocenter currently lives in
                        (either `mol` or an unattached fragment -- both
                        are static at this point, so their existing
                        torsions are ground truth)
stereo_pos            : local index of the atom being functionalized
                         (i.e. one of the two atoms of the double bond)
stereos                : {(i, j): 'cis'/'trans'}
ref_exclude             : local indices to exclude when hunting for the
                          retained reference substituent on the far atom
                          (pass [target_pos] etc. as needed)

Returns the local index of the hydrogen to use as the attachment site
(i.e. the one to pass as `target_fragment` or `group_site`), or None if
there's no ambiguity to resolve (caller should fall back to its default
H-picking logic).


<a id="Psience.Molecools.Builder.MoleculeBuilder.attach_functional_group" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
attach_functional_group(cls, mol, target_fragment, atoms, new_coords, bonds='recompute', ref=None, masses=None, distance='auto', angle=0, dihedral='auto', bond_sites=None, dihedral_search_steps=36, dihedral_distance_metric=None, embedding='auto', bond_order=None, use_absolue_posititions=False, group_site=None, molecule_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L356)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L356?message=Update%20Docs)]
</div>
Build a copy of `mol` with a new group of atoms (`atoms`/`new_coords`) attached at
`target_fragment`, positioning and orienting the new group using the fragment's local
reference frame (bond distance/angle/dihedral, or an explicit embedding), and splicing
the corresponding bonds into the result; supports designating a `group_site` atom within
the new group as its attachment point, in which case the method recurses after
re-deriving the embedding/bonds relative to that site.

**Multiple bond sites.** More than one bond can be formed at once by passing a list of
binding sites: `target_fragment` becomes a list (one target fragment per site), `ref`
and `bond_order` may be given as matching per-site lists, and `group_site` becomes a
list of the group atoms to replace at each corresponding site. As in the single-site
`group_site` path, each `group_site` atom is a *placeholder that is removed* -- typically
a hydrogen -- and the heavy atom it was bonded to becomes the real attachment atom that
bonds to the scaffold. This lets a bidentate (or polydentate) ligand be attached to an
existing scaffold at several points at once, closing a ring or forming a linker. The
group is oriented so its attachment atoms line up with the scaffold's per-site bond
directions, with the overall frame anchored on the reference atoms of the *first* binding
site by default (`embedding='auto'`); an explicit `embedding` overrides this.
  - `mol`: `Molecule`
    > the scaffold molecule the group is attached to
  - `target_fragment`: `int | Iterable[int] | Iterable[Iterable[int]]`
    > the atom(s) of `mol` the new group attaches to/replaces; a list
    of target fragments (one per site) selects multi-site mode
  - `atoms`: `Iterable[str]`
    > the element symbols of the atoms in the new group
  - `new_coords`: `np.ndarray`
    > the (local) coordinates of the new group's atoms
  - `bonds`: `str | list | None`
    > bonds within the new group; `'recompute'` to guess them fresh, `None` to
    reuse `mol.bonds` remapped, or an explicit bond list
  - `ref`: `Iterable[int] | Iterable[Iterable[int]] | None`
    > reference atom(s) used to anchor the attachment frame; computed automatically
    if not given. In multi-site mode this may be a per-site list of reference-atom lists
  - `masses`: `np.ndarray | None`
    > masses for the new group's atoms; looked up from `atoms` if not given
  - `distance`: `str | float | None`
    > the bond distance to place the new group at; `'auto'` to look it up from
    `BondData`, or `None`/a number
  - `angle`: `float`
    > rotation angle (about the up-vector) to apply to the new group
  - `dihedral`: `float | str`
    > rotation angle (about the offset axis) to apply to the new group;
    `'auto'` to instead scan `dihedral_search_steps` evenly-spaced angles and keep
    whichever maximizes `dihedral_distance_metric` between the new group and the rest of
    the scaffold
  - `dihedral_search_steps`: `int`
    > number of evenly-spaced angles (over 360°) to try when
    `dihedral='auto'`
  - `dihedral_distance_metric`: `Callable | None`
    > `(frag_coords, other_coords) -> float` scoring function
    used when `dihedral='auto'`; higher is better. Defaults to the average pairwise
    distance between the new group's atoms and the surviving scaffold atoms
  - `embedding`: `str | tuple | np.ndarray | None`
    > the reference orientation for the new group; `'auto'` to derive it from
    moments of inertia (single site) or from the per-site bond geometry anchored on the
    first site (multi-site), or an explicit `(origin, axes)`/axes specification
  - `bond_order`: `float | Iterable[float] | None`
    > the bond order connecting the new group to the target fragment;
    defaults to `1` (or inferred when `group_site` is used). In multi-site mode may be a
    per-site list
  - `use_absolue_posititions`: `bool`
    > whether `new_coords` should be used as absolute
    coordinates rather than being repositioned relative to the fragment frame
  - `group_site`: `int | Iterable[int] | None`
    > index (within `atoms`/`new_coords`) of the placeholder atom that
    marks the attachment point; the placeholder is removed and its heavy neighbor becomes
    the atom that bonds to the scaffold. A single int uses the single-site path; a list
    of ints selects multi-site mode (each placeholder removed and its neighbor bonded at
    the corresponding site), e.g. two hydrogens replaced when a bidentate ligand closes a
    ring or forms a linker
  - `molecule_type`: `type | None`
    > the molecule class used to build any intermediate molecules;
    defaults to the standard `Molecule`
  - `:returns`: `Molecule`
    > the molecule with the new group attached


<a id="Psience.Molecools.Builder.MoleculeBuilder.from_fragments" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_fragments(cls, scaffold, *replacements, active_sites=None, chiralities=None, stereos=None, bond_orders=None, atom_replacements=None, cache=None, add_implicit_hydrogens='full', remove_sites=False, recompute_properties=True, reorder_from_atom_map=False, molecule_type=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L676)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L676?message=Update%20Docs)]
</div>
Construct a molecule from a `scaffold` plus a set of replacement fragments, driving the
connectivity via a templated-SMILES join (`build_templated_smiles`) and the 3D geometry
via repeated `attach_functional_group` calls, then reordering the atoms to match the
final SMILES numbering.
  - `scaffold`: `str | Molecule | dict`
    > the base fragment: a SMILES string, a `Molecule`, or a dict with any of
    `molecule`/`smiles`/`atoms`/`coords`/`bonds`
  - `replacements`: `str | Molecule | dict`
    > the replacement fragments, each in the same accepted forms as
    `scaffold`
  - `active_sites`: `Any`
    > forwarded to `build_templated_smiles`
  - `chiralities`: `Any`
    > forwarded to `build_templated_smiles`
  - `stereos`: `Any`
    > `{(i, j): 'cis'/'trans'}` stereochemistry constraints
  - `bond_orders`: `Any`
    > forwarded to `build_templated_smiles`
  - `atom_replacements`: `Any`
    > forwarded to `build_templated_smiles`
  - `cache`: `Any`
    > SMILES-parsing cache forwarded through
  - `add_implicit_hydrogens`: `Any`
    > implicit-hydrogen handling forwarded through
  - `remove_sites`: `Any`
    > forwarded to `build_templated_smiles`
  - `recompute_properties`: `bool`
    > whether to rebuild the final molecule cleanly from the joined
    SMILES (carrying over the assembled coordinates), or just reorder the assembled
    molecule in place
  - `molecule_type`: `type | None`
    > the molecule class to build with; defaults to the standard
    `Molecule`
  - `opts`: `dict`
    > extra options forwarded to the molecule constructor / `from_string`
  - `:returns`: `Molecule`
    > the assembled molecule
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Builder/MoleculeBuilder.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Builder/MoleculeBuilder.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Builder/MoleculeBuilder.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Builder/MoleculeBuilder.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Builder.py#L36?message=Update%20Docs)   
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