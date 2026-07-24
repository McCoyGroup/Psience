"""
Provides a simple Molecule class that we can adapt as we need
Uses AtomData to get properties and whatnot
"""
from __future__ import annotations

import io
import os, numpy as np
import tempfile
import typing
import collections

from McUtils.Data import AtomData, UnitsData, BondData
from McUtils.Coordinerds import (
    CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D,
    PrimitiveCoordinatePicker, RedundantCoordinateGenerator
)
import McUtils.Numputils as nput
import McUtils.Iterators as itut
import McUtils.Devutils as dev
import McUtils.Coordinerds as coordops
import McUtils.Zachary as zach
from McUtils.Zachary import Mesh
from McUtils.Graphs import EdgeGraph
import McUtils.Plots as plt
from McUtils.ExternalPrograms import RDMolecule, ExternalProgramJob
from McUtils.Scaffolding import Logger
import McUtils.Symmetry as symm

from .MoleculeInterface import *

from .CoordinateSystems import MolecularEmbedding, ModeEmbedding
from .Evaluator import (
    MolecularEvaluator, EnergyEvaluator, DipoleEvaluator, ChargeEvaluator,
    ReducedDimensionalPotentialHandler, DipolePolarizabilityEvaluator
)
from .Hamiltonian import MolecularHamiltonian
from .Properties import *
from .Serializers import MoleculePropertyCache

__all__ = [
    "Molecule",
    "MolecoolException"
]

__reload_hook__ = ["..Modes", ".MoleculeInterface", '.CoordinateSystems', '.Hamiltonian', '.Evaluator', '.Properties']

from .Transformations import MolecularTransformation


class Molecule(AbstractMolecule):
    """
    General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all
    """

    def __init__(self,
                 atoms,
                 coords,
                 bonds=None,
                 masses=None,
                 name=None,
                 internals=None,
                 rdmol=None,
                 dipole_surface=None,
                 dipole_derivatives=None,
                 potential_surface=None,
                 potential_derivatives=None,
                 normal_modes=None,
                 source_file=None,
                 guess_bonds=True,
                 charge=None,
                 formal_charges=None,
                 spin=None,
                 display_mode=None,
                 display_settings=None,
                 energy=None,
                 energy_evaluator=None,
                 dipole_evaluator=None,
                 charge_evaluator=None,
                 polarizability_evaluator=None,
                 polarizability_derivatives=None,
                 checkpoint_file=None,
                 **metadata
                 ):
        """
        :param atoms: atoms specified by name, either full name or short
        :type atoms: Iterable[str]
        :param coords: coordinates for the molecule, assumed to be in Bohr by default
        :type coords: np.ndarray | Iterable[Iterable[float]]
        :param bonds: bond specification for the molecule
        :type bonds: Iterable[Iterable[int]] | None
        :param obmol: OpenBabel molecule for doing conversions
        :type obmol:
        :param charge: Net charge on the molecule
        :type charge: int | None
        :param name: Name for the molecule
        :type name: str | None
        :param name: The internal coordinate specification for the molecule
        :type name: np.ndarray[int] | None
        :param dipole_surface: The dipole surface for the system
        :type dipole_surface: DipoleSurface | None
        :param dipole_derivatives: Derivatives of the dipole surface
        :type dipole_derivatives: Iterable[np.ndarray] | None
        :param potential_surface: The potential surface for the system
        :type potential_surface: PotentialSurface | None
        :param potential_derivatives: Derivatives of the potential surface
        :type potential_derivatives: Iterable[np.ndarray] | None
        :param guess_bonds: Whether or not to guess the bonding arrangement when that would be used
        :type guess_bonds: bool
        :param source_file: The data file the molecule was loaded from
        :type source_file: str
        :param kw: Other bound parameters that might be useful
        :type kw:
        """
        # convert "atoms" into list of atom data
        if isinstance(atoms, np.ndarray):
            atoms = atoms.tolist()
        self._ats = [
            AtomData[atom]
                if isinstance(atom, (int, np.integer, str, np.str_)) else atom
            for atom in atoms
        ]
        self._mass = np.array([a["Mass"] for a in self._ats]) if masses is None else np.asanyarray(masses)
        coords = CoordinateSet(coords, CartesianCoordinates3D)


        # properties to be returned
        self._bonds = bonds
        self._edge_graph = None
        self.guess_bonds = guess_bonds

        if charge is not None:
            metadata['charge'] = charge
        if spin is not None:
            metadata['spin'] = spin
        if formal_charges is not None:
            metadata['formal_charges'] = formal_charges
        self._meta = metadata
        self._fragment_indices = None

        self._rdmol = rdmol
        # a little messy
        self._embedding = MolecularEmbedding(self.atomic_masses, coords, None)
        internals = self.canonicalize_internals(internals, self.atoms, coords, bonds, masses=self._mass)
        self._embedding = MolecularEmbedding(self.atomic_masses, coords, internals)
        self._mode_embedding = None

        self._name = name


        self._src = None
        self.source_file = source_file

        self._rdmol = rdmol
        self._dips = DipoleSurfaceManager(self,
                                          surface=dipole_surface,
                                          derivatives=dipole_derivatives,
                                          polarizability_derivatives=polarizability_derivatives
                                          )
        self._pes = PotentialSurfaceManager(self,
                                            surface=potential_surface,
                                            derivatives=potential_derivatives
                                            )

        self._normal_modes = NormalModesManager(self, normal_modes=normal_modes)

        self._evaluator = None
        self._hamiltonian = None

        if display_mode is None:
            display_mode = self.default_display_mode
        self.display_mode = display_mode
        self.display_settings = display_settings
        self._energy = energy
        self.energy_evaluator = energy_evaluator
        self.dipole_evaluator = dipole_evaluator
        self.charge_evaluator = charge_evaluator
        self.polarizability_evaluator = polarizability_evaluator

        # if checkpoint_file is True:
        #     checkpoint_file = self.molecule_hash() # TODO: define this based on DB work
        self.checkpoint = MoleculePropertyCache(self, checkpoint_file)

    def modify(self,
               atoms=dev.default,
               coords=dev.default,
               *,
               internals=dev.default,
               masses=dev.default,
               bonds=dev.default,
               guess_bonds=dev.default,
               energy=dev.default,
               energy_evaluator=dev.default,
               dipole_evaluator=dev.default,
               charge_evaluator=dev.default,
               polarizability_evaluator=dev.default,
               charge=dev.default,
               spin=dev.default,
               rdmol=dev.default,
               display_mode=dev.default,
               display_settings=dev.default,
               normal_modes=dev.default,
               dipole_surface=dev.default,
               potential_surface=dev.default,
               dipole_derivatives=dev.default,
               potential_derivatives=dev.default,
               polarizability_derivatives=dev.default,
               meta=dev.default,
               source_file=dev.default
               ):
        """
        **LLM Docstring**

        Build a new `Molecule` that is a copy of this one with the given fields overridden, treating any argument left at its `dev.default` sentinel as "keep the current value" (with some fields, like masses/energy/surfaces, only carried over automatically when the arguments they depend on -- e.g. `atoms`, `coords`, `energy_evaluator` -- are also left unspecified).

        :param atoms: replacement atoms, or `dev.default` to keep the current ones
        :type atoms: Iterable[str] | object
        :param coords: replacement coordinates, or `dev.default` to keep the current ones
        :type coords: np.ndarray | object
        :param internals: replacement internal-coordinate specification
        :type internals: object
        :param masses: replacement masses
        :type masses: np.ndarray | object
        :param bonds: replacement bonds
        :type bonds: object
        :param guess_bonds: replacement bond-guessing flag
        :type guess_bonds: bool | object
        :param energy: replacement cached energy value
        :type energy: float | object
        :param energy_evaluator: replacement energy evaluator
        :type energy_evaluator: object
        :param dipole_evaluator: replacement dipole evaluator
        :type dipole_evaluator: object
        :param charge_evaluator: replacement charge evaluator
        :type charge_evaluator: object
        :param polarizability_evaluator: replacement polarizability evaluator
        :type polarizability_evaluator: object
        :param charge: replacement net charge
        :type charge: int | object
        :param spin: replacement spin
        :type spin: object
        :param rdmol: replacement RDKit molecule
        :type rdmol: object
        :param display_mode: replacement display mode
        :type display_mode: str | object
        :param display_settings: replacement display settings
        :type display_settings: dict | object
        :param normal_modes: replacement normal modes
        :type normal_modes: object
        :param dipole_surface: replacement dipole surface
        :type dipole_surface: object
        :param potential_surface: replacement potential surface
        :type potential_surface: object
        :param dipole_derivatives: replacement dipole derivatives
        :type dipole_derivatives: object
        :param potential_derivatives: replacement potential derivatives
        :type potential_derivatives: object
        :param polarizability_derivatives: replacement polarizability derivatives
        :type polarizability_derivatives: object
        :param meta: replacement/merged metadata
        :type meta: dict | object
        :param source_file: replacement source file path
        :type source_file: str | object
        :return: the new, modified `Molecule`
        :rtype: Molecule
        """
        return type(self)(
            self.atoms if dev.is_default(atoms) else atoms,
            self.coords if dev.is_default(coords) else coords,
            masses=(
                self.masses if
                    dev.is_default(masses, allow_None=False) and dev.is_default(atoms)
                else None if
                    dev.is_default(masses)
                else masses
            ),
            bonds=self._bonds if dev.is_default(bonds, allow_None=False) else bonds,
            guess_bonds=self.guess_bonds if dev.is_default(guess_bonds) else guess_bonds,
            energy_evaluator=self.energy_evaluator if dev.is_default(energy_evaluator, allow_None=False) else energy_evaluator,
            dipole_evaluator=self.dipole_evaluator if dev.is_default(dipole_evaluator, allow_None=False) else dipole_evaluator,
            charge_evaluator=self.charge_evaluator if dev.is_default(charge_evaluator, allow_None=False) else charge_evaluator,
            polarizability_evaluator=self.polarizability_evaluator if dev.is_default(polarizability_evaluator, allow_None=False) else polarizability_evaluator,
            charge=self.charge if dev.is_default(charge) else charge,
            spin=self.spin if dev.is_default(spin) else spin,
            internals=self.internals if dev.is_default(internals, allow_None=False) else internals,
            normal_modes=self.normal_modes if dev.is_default(normal_modes, allow_None=False) else normal_modes,
            energy=self._energy if (
                    dev.is_default(energy, allow_None=False)
                    and dev.is_default(coords)
                    and dev.is_default(energy_evaluator, allow_None=False)
            ) else energy,
            dipole_surface=self.dipole_surface if (
                    dev.is_default(dipole_surface, allow_None=False)
                    and dev.is_default(dipole_derivatives, allow_None=False)
                    and dev.is_default(coords)
                    and dev.is_default(energy_evaluator, allow_None=False)
            ) else (
                None
                    if dev.is_default(dipole_surface) else
                dipole_surface
            ),
            dipole_derivatives=None if dev.is_default(dipole_derivatives) else dipole_derivatives,
            polarizability_derivatives=None if dev.is_default(polarizability_derivatives) else polarizability_derivatives,
            potential_surface=self.potential_surface if (
                    dev.is_default(potential_surface, allow_None=False)
                    and dev.is_default(potential_derivatives, allow_None=False)
                    and dev.is_default(coords)
                    and dev.is_default(energy_evaluator, allow_None=False)
            ) else (
                None
                    if dev.is_default(potential_surface) else
                potential_surface
            ),
            rdmol=(
                self._rdmol
                   if (
                        dev.is_default(rdmol)
                        and dev.is_default(atoms)
                        and dev.is_default(coords)
                        and dev.is_default(bonds)
                    ) else
                None
            ),
            display_mode=(
                self.display_mode
                    if dev.is_default(display_mode, allow_None=False) else
                display_mode
            ),
            display_settings=(
                self.display_settings
                    if dev.is_default(display_settings, allow_None=False) else
                display_settings
            ),
            potential_derivatives=(
                None
                    if dev.is_default(potential_derivatives) else
                potential_derivatives
            ),
            meta=(
                     self._meta
                        if dev.is_default(meta) else
                     {}
                        if meta is None else
                     dev.merge_dicts(self._meta, meta)
            ),
            source_file=self._src if dev.is_default(source_file) else source_file
        )

    def __del__(self):
        """
        **LLM Docstring**

        Clean up the molecule's coordinate embedding (deregistering any converters it registered) when the object is garbage collected.

        :return: None
        :rtype: None
        """
        if hasattr(self, '_embedding'):
            self._embedding.cleanup()

    def to_state(self, serializer=None):
        """
        **LLM Docstring**

        Serialize this molecule's essential data (atoms, masses, coordinates, bonds, internal-coordinate spec, evaluators, potential/dipole derivatives, charge, spin) into a plain dict, stripping out non-serializable embedding-specific converter options first.

        :param serializer: accepted for interface consistency but not used in this method's body
        :type serializer: object | None
        :return: the serialized state dict
        :rtype: dict
        """
        internals = self.internals
        if internals is not None:
            internals = internals.copy()
            converter_options = internals.get('converter_options')
            if converter_options is not None:
                converter_options = converter_options.copy()
                for k in [
                    'embedding_coords',
                    'jacobian_prep',
                    'origin',
                    'axes'
                ]:
                    converter_options.pop(k, None)
                internals['converter_options'] = converter_options

            # if internals.get('zmatrix') is not None: # stores a jacobian prepper we don't need to serialize
            #     internals = {
            #         'zmatrix':internals['zmatrix']
            #     }
            # elif internals.get('specs') is not None:
            #     if 'redundant_transformation' in internals:
            #
            #         internals = {
            #             'primitives':internals['specs']
            #         }
            #     else:
            #         internals = {
            #             'specs':internals['specs']
            #         }
            #     converter_options = internals.get('converter_options')

        data = {
            'atoms':self.atoms,
            'masses':self.masses,
            'coords':self.coords,
            'bonds':self.bonds,
            'internals':internals,
            'energy_evaluator':self.energy_evaluator,
            'dipole_evaluator':self.dipole_evaluator,
            'charge_evaluator':self.charge_evaluator,
            'polarizability_evaluator':self.polarizability_evaluator,
            'potential_derivatives':self.potential_derivatives,
            'dipole_derivatives':self.dipole_derivatives,
            'charge':self.charge,
            'spin':self.spin
        }
        return data
    @classmethod
    def from_state(cls, data, serializer=None):
        """
        **LLM Docstring**

        Reconstruct a `Molecule` from a previously serialized state dict, by passing its entries directly as constructor keyword arguments.

        :param data: the serialized state, as produced by `to_state`
        :type data: dict
        :param serializer: accepted for interface consistency but not used in this method's body
        :type serializer: object | None
        :return: the reconstructed molecule
        :rtype: Molecule
        """
        return cls(**data)

    def cached_eval(self,
                    key, generator,
                    *,
                    condition=None,
                    args=(),
                    kwargs=None):
        """
        **LLM Docstring**

        Evaluate (and cache, via this molecule's on-disk/in-memory checkpoint) a value under `key`, delegating to `self.checkpoint.cached_eval`.

        :param key: the cache key to look up or populate
        :type key: str
        :param generator: callable used to compute the value when it is not already cached
        :type generator: callable
        :param condition: optional predicate controlling whether the cached value should be recomputed
        :type condition: callable | None
        :param args: positional arguments passed to `generator` if it is called
        :type args: tuple
        :param kwargs: keyword arguments passed to `generator` if it is called
        :type kwargs: dict | None
        :return: the cached or newly computed value
        :rtype: object
        """
        return self.checkpoint.cached_eval(
            key,
            generator,
            condition=condition,
            args=args,
            kwargs=kwargs
        )

    @classmethod
    def _generate_auto_spec(cls, atoms, bonds, base_coords=None, **opts):
        """
        **LLM Docstring**

        Automatically pick a set of primitive internal coordinates (bond stretches, angles, dihedrals) from the bonding graph, via `PrimitiveCoordinatePicker`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param bonds: the bonds to build coordinates from
        :type bonds: Iterable[Iterable[int]]
        :param base_coords: coordinates to seed/prioritize the picker with
        :type base_coords: Iterable | None
        :param opts: extra options forwarded to `PrimitiveCoordinatePicker`
        :type opts: dict
        :return: the picked coordinate specs
        :rtype: list
        """
        return PrimitiveCoordinatePicker(
            atoms,
            [b[:2] for b in bonds],
            base_coords=base_coords,
            **opts
        ).coords

    @classmethod
    def _generate_stretch_spec(cls, atoms, bonds, **opts):
        """
        **LLM Docstring**

        Build a "natural"-coordinate specification consisting only of the bond-stretch coordinates implied by `bonds`.

        :param atoms: the atom labels (unused directly, kept for interface consistency with `_generate_auto_spec`)
        :type atoms: Iterable[str]
        :param bonds: the bonds to build stretch coordinates from
        :type bonds: Iterable[Iterable[int]]
        :param opts: extra options, unused
        :type opts: dict
        :return: the list of bond-stretch coordinate specs
        :rtype: list
        """
        return sum(coordops.get_stretch_coordinate_system([tuple(s[:2]) for s in bonds]), [])

    @classmethod
    def _auto_auto_spec(cls, spec_generator, atoms, coords, bonds, redundant=False, base_coordinates=None,
                        masses=None,
                        untransformed_coordinates=None,
                        prune_coordinates=True,
                        pruning_options=None,
                        formal_charges=None,
                        **opts):
        """
        **LLM Docstring**

        Shared driver behind `_auto_spec`/`_stretch_spec`: guesses bonds if not given, optionally sets up a redundant-coordinate specification (folding in any `untransformed_coordinates`), generates the primitive coordinate specs via `spec_generator`, and (optionally) prunes them down to a well-conditioned, non-redundant subset via `RedundantCoordinateGenerator.prune_coordinate_specs`.

        :param spec_generator: the coordinate-generating function to use (`_generate_auto_spec` or `_generate_stretch_spec`)
        :type spec_generator: callable
        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use; guessed via RDKit if `None`
        :type bonds: Iterable[Iterable[int]] | None
        :param redundant: whether to build a redundant coordinate specification
        :type redundant: bool
        :param base_coordinates: seed coordinates to prioritize/include
        :type base_coordinates: Iterable | None
        :param masses: atomic masses, used for pruning; computed from `atoms` if not given
        :type masses: np.ndarray | None
        :param untransformed_coordinates: coordinates that should remain untransformed under the redundant transformation
        :type untransformed_coordinates: Iterable | None
        :param prune_coordinates: whether to prune the generated coordinate specs down to a well-conditioned subset
        :type prune_coordinates: bool
        :param pruning_options: extra options forwarded to `RedundantCoordinateGenerator.prune_coordinate_specs`
        :type pruning_options: dict | None
        :param formal_charges: formal charges used when guessing bonds
        :type formal_charges: Iterable | None
        :param opts: extra options forwarded to `spec_generator`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict (with `'specs'`, and `'redundant'`/`'untransformed_coordinates'` if applicable)
        :rtype: dict
        """
        base_coords = base_coordinates
        if bonds is None:
            bonds = RDMolecule.from_coords(
                                           atoms,
                                           coords * UnitsData.convert("BohrRadius", "Angstroms"),
                                           bonds,
                                           formal_charges=formal_charges,
                                           guess_bonds=True
                                           ).bonds
        if redundant and untransformed_coordinates is None:
            untransformed_coordinates = base_coords
            base_coords = None
        if redundant and untransformed_coordinates is not None:
            untransformed_coordinates = PrimitiveCoordinatePicker.prep_unique_coords(untransformed_coordinates)
            base_coords = untransformed_coordinates + ([] if base_coords is None else list(base_coords))
        if base_coords is not None:
            base_coords = PrimitiveCoordinatePicker.prep_unique_coords(base_coords)
        specs = spec_generator(atoms, bonds, base_coords=base_coords, **opts)
        if prune_coordinates:
            if pruning_options is None:
                pruning_options = {}
            expansion = nput.internal_coordinate_tensors(coords, specs, order=1)[1:]
            if masses is None:
                ats = [AtomData[atom] if isinstance(atom, (int, np.integer, str)) else atom for atom in atoms]
                masses = np.array([a["Mass"] for a in ats])
            prune_pos = RedundantCoordinateGenerator.prune_coordinate_specs(
                expansion,
                masses=masses,
                untransformed_coordinates=np.arange(len(base_coords)) if base_coords is not None else None,
                **pruning_options
            )
            specs = [specs[i] for i in prune_pos]
        spec = {'specs':specs}
        if redundant:
            spec['redundant'] = True
            if untransformed_coordinates is not None:
                spec['untransformed_coordinates'] = np.arange(len(untransformed_coordinates))
        return spec

    @classmethod
    def _auto_spec(cls, atoms, coords, bonds, **opts):
        """
        **LLM Docstring**

        Build an automatically-chosen internal-coordinate specification from the bonding graph, via `_auto_auto_spec` with `_generate_auto_spec`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use
        :type bonds: Iterable[Iterable[int]] | None
        :param opts: extra options forwarded to `_auto_auto_spec`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict
        :rtype: dict
        """
        return cls._auto_auto_spec(cls._generate_auto_spec, atoms, coords, bonds, **opts)
    @classmethod
    def _stretch_spec(cls, atoms, coords, bonds, **opts):
        """
        **LLM Docstring**

        Build a "natural"/stretch-only internal-coordinate specification from the bonding graph, via `_auto_auto_spec` with `_generate_stretch_spec`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use
        :type bonds: Iterable[Iterable[int]] | None
        :param opts: extra options forwarded to `_auto_auto_spec`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict
        :rtype: dict
        """
        return cls._auto_auto_spec(cls._generate_stretch_spec, atoms, coords, bonds, **opts)

    # @classmethod
    def canonicalize_internals(self, spec, atoms, coords, bonds, relocalize=True, masses=None):
        """
        **LLM Docstring**

        Normalize the many accepted forms of an internal-coordinate specification (the strings `'auto'`/`'zmatrix'`, a dict with `'primitives'`/`'specs'`/`'zmatrix'` keys where `'specs'` may itself be `'auto'`/`'natural'`, a bare Z-matrix-like array, or a bare list of primitive specs) down into the canonical dict form expected by `MolecularEmbedding`, recursively re-dispatching as needed.

        :param spec: the internal-coordinate specification to canonicalize
        :type spec: str | dict | Iterable | None
        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use when auto-generating coordinates
        :type bonds: Iterable[Iterable[int]] | None
        :param relocalize: whether redundant coordinates should be relocalized by default
        :type relocalize: bool
        :param masses: atomic masses, forwarded to the auto-generation routines
        :type masses: np.ndarray | None
        :return: the canonicalized specification
        :rtype: dict | None
        :raises ValueError: if `spec` is a string that isn't recognized (`'auto'`/`'zmatrix'`/`'natural'`)
        """
        if isinstance(spec, str) and spec.lower() == 'auto':
            spec = {
                'primitives': 'auto'
            }
        elif dev.str_is(spec, 'zmatrix'):
            spec = self.get_bond_zmatrix()

        if isinstance(spec, str):
            # if spec.lower() == 'auto':
            #     spec = cls._auto_spec(atoms, coords, bonds)
            # else:
            raise ValueError(f"can't understand internal spec '{spec}'")
        elif isinstance(spec, dict):
            if 'zmatrix' in spec: return spec
            prims = spec.pop('primitives', None)
            if prims is not None:
                spec = spec.copy()
                spec['specs'] = prims
                spec['redundant'] = True
            subspec = spec.get('specs', '')
            if isinstance(subspec, str):
                if subspec.lower() == 'auto':
                    opts = spec.copy()
                    del opts['specs']
                    if 'relocalize' in opts:
                        relocalize = spec.get('relocalize', relocalize)
                        del opts['relocalize']
                    spec = self._auto_spec(atoms, coords, bonds, masses=masses, **opts)
                elif subspec.lower() == 'natural':
                    opts = spec.copy()
                    del opts['specs']
                    if 'relocalize' in opts:
                        relocalize = spec.get('relocalize', relocalize)
                        del opts['relocalize']
                    spec = self._stretch_spec(atoms, coords, bonds, masses=masses, **opts)
                else:
                    raise ValueError(f"can't understand internal spec '{spec}'")
            else:
                untransformed_coordinates = spec.get('untransformed_coordinates')
                if untransformed_coordinates is not None:
                    if not nput.is_int(untransformed_coordinates[0]):
                        prims = spec.get('specs')
                        untransformed_coordinates = [
                            prims.index(u)
                            for u in untransformed_coordinates
                        ]
                    spec['untransformed_coordinates'] = untransformed_coordinates
            if spec.get('redundant'):
                spec['relocalize'] = spec.get('relocalize', relocalize)
        elif not isinstance(spec, dict) and spec is not None:
            if all(not isinstance(x, dict) and len(x) == 4 for x in spec):
                spec = {'zmatrix': spec}
            else:
                spec = {'primitives':spec}
            spec = self.canonicalize_internals(spec, atoms, coords, bonds, relocalize=relocalize, masses=masses)
        return spec
    def prep_internal_spec(self, spec, relocalize=True, masses=None):
        """
        **LLM Docstring**

        Canonicalize an internal-coordinate specification against this molecule's own atoms, coordinates, bonds, and masses, via `canonicalize_internals`.

        :param spec: the internal-coordinate specification to canonicalize
        :type spec: object
        :param relocalize: whether redundant coordinates should be relocalized by default
        :type relocalize: bool
        :param masses: atomic masses to use instead of this molecule's own
        :type masses: np.ndarray | None
        :return: the canonicalized specification
        :rtype: dict | None
        """
        return self.canonicalize_internals(
            spec,
            self.atoms,
            self.coords,
            self._bonds,
            relocalize=relocalize,
            masses=masses
        )

    @property
    def embedding(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's `MolecularEmbedding`. Setting it also resets the cached evaluator and Hamiltonian, since both depend on the embedding.

        :param e: (setter only) the new embedding
        :type e: MolecularEmbedding
        :return: (getter) the current embedding
        :rtype: MolecularEmbedding
        """
        return self._embedding
    @embedding.setter
    def embedding(self, e):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's `MolecularEmbedding`. Setting it also resets the cached evaluator and Hamiltonian, since both depend on the embedding.

        :param e: (setter only) the new embedding
        :type e: MolecularEmbedding
        :return: (getter) the current embedding
        :rtype: MolecularEmbedding
        """
        self._embedding = e
        self.evaluator = None
        self.hamiltonian = None
    def get_evaluator(self, embedding=None, normal_modes=dev.default):
        """
        **LLM Docstring**

        Build a `MolecularEvaluator` for this molecule's embedding (or an alternate one), using either the given `normal_modes` or this molecule's own.

        :param embedding: an alternate embedding to build the evaluator for; defaults to `self.embedding`
        :type embedding: MolecularEmbedding | None
        :param normal_modes: normal modes to use instead of `self._normal_modes`
        :type normal_modes: object
        :return: the constructed evaluator
        :rtype: MolecularEvaluator
        """
        if embedding is None: embedding = self.embedding
        if dev.is_default(normal_modes, allow_None=False): normal_modes = self._normal_modes
        return MolecularEvaluator(embedding, normal_modes)
    @property
    def evaluator(self):
        """
        **LLM Docstring**

        Property getter/setter for the `MolecularEvaluator` used for energy/derivative calculations. The getter lazily builds one via `get_evaluator` the first time it's needed.

        :param e: (setter only) the new evaluator
        :type e: MolecularEvaluator
        :return: (getter) the cached (or newly built) evaluator
        :rtype: MolecularEvaluator
        """
        if self._evaluator is None:
            self._evaluator = self.get_evaluator()
        return self._evaluator
    @evaluator.setter
    def evaluator(self, e):
        """
        **LLM Docstring**

        Property getter/setter for the `MolecularEvaluator` used for energy/derivative calculations. The getter lazily builds one via `get_evaluator` the first time it's needed.

        :param e: (setter only) the new evaluator
        :type e: MolecularEvaluator
        :return: (getter) the cached (or newly built) evaluator
        :rtype: MolecularEvaluator
        """
        self._evaluator = e
    #region Base Coords
    @property
    def coords(self):
        """
        **LLM Docstring**

        Property getter/setter for the Cartesian coordinates, delegating to `self.embedding.coords`.

        :param coords: (setter only) the new Cartesian coordinates
        :type coords: np.ndarray
        :return: (getter) the Cartesian coordinates
        :rtype: CoordinateSet
        """
        return self.embedding.coords
    @coords.setter
    def coords(self, coords):
        """
        **LLM Docstring**

        Property getter/setter for the Cartesian coordinates, delegating to `self.embedding.coords`.

        :param coords: (setter only) the new Cartesian coordinates
        :type coords: np.ndarray
        :return: (getter) the Cartesian coordinates
        :rtype: CoordinateSet
        """
        self.embedding.coords = coords
    @property
    def masses(self):
        """
        **LLM Docstring**

        Property getter/setter for the atomic masses. The setter also updates the embedding's masses (via `atomic_masses`, i.e. in atomic units).

        :param masses: (setter only) the new masses
        :type masses: np.ndarray
        :return: (getter) the atomic masses
        :rtype: np.ndarray
        """
        return self._mass
    @masses.setter
    def masses(self, masses):
        """
        **LLM Docstring**

        Property getter/setter for the atomic masses. The setter also updates the embedding's masses (via `atomic_masses`, i.e. in atomic units).

        :param masses: (setter only) the new masses
        :type masses: np.ndarray
        :return: (getter) the atomic masses
        :rtype: np.ndarray
        """
        self._mass = masses
        self.embedding.masses = self.atomic_masses
    @property
    def internals(self):
        """
        **LLM Docstring**

        Getter for the raw (canonicalized) internal-coordinate specification, delegating to `self.embedding.internals`. (A setter with the same name separately rebuilds the embedding from a new specification.)

        :return: the internal-coordinate specification, or `None` if none is set
        :rtype: dict | None
        """
        return self.embedding.internals
    @property
    def charge(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's net charge, stored in its metadata (getter defaults to `0` if unset).

        :param c: (setter only) the new net charge
        :type c: int
        :return: (getter) the net charge
        :rtype: int
        """
        return self._meta.get('charge', 0)
    @charge.setter
    def charge(self, c):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's net charge, stored in its metadata (getter defaults to `0` if unset).

        :param c: (setter only) the new net charge
        :type c: int
        :return: (getter) the net charge
        :rtype: int
        """
        self._meta['charge'] = c
    @property
    def spin(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's spin, stored in its metadata.

        :param c: (setter only) the new spin value
        :type c: object
        :return: (getter) the spin, or `None` if unset
        :rtype: object | None
        """
        return self._meta.get('spin')
    @spin.setter
    def spin(self, c):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's spin, stored in its metadata.

        :param c: (setter only) the new spin value
        :type c: object
        :return: (getter) the spin, or `None` if unset
        :rtype: object | None
        """
        self._meta['spin'] = c
    @property
    def charges(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's per-atom partial charges, stored in its metadata.

        :param c: (setter only) the new per-atom charges
        :type c: np.ndarray
        :return: (getter) the per-atom charges, or `None` if unset
        :rtype: np.ndarray | None
        """
        return self._meta.get('charges', None)
    @charges.setter
    def charges(self, c):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's per-atom partial charges, stored in its metadata.

        :param c: (setter only) the new per-atom charges
        :type c: np.ndarray
        :return: (getter) the per-atom charges, or `None` if unset
        :rtype: np.ndarray | None
        """
        self._meta['charges'] = c
    @property
    def formal_charges(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's per-atom formal charges, stored in its metadata.

        :param c: (setter only) the new per-atom formal charges
        :type c: np.ndarray
        :return: (getter) the per-atom formal charges, or `None` if unset
        :rtype: np.ndarray | None
        """
        return self._meta.get('formal_charges', None)
    @formal_charges.setter
    def formal_charges(self, c):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's per-atom formal charges, stored in its metadata.

        :param c: (setter only) the new per-atom formal charges
        :type c: np.ndarray
        :return: (getter) the per-atom formal charges, or `None` if unset
        :rtype: np.ndarray | None
        """
        self._meta['formal_charges'] = c

    def get_charge_evaluator(self, evaluator=None, **opts):
        """
        **LLM Docstring**

        Resolve (and, if needed, instantiate from this molecule) a charge-evaluator object, defaulting to `self.charge_evaluator` if none is given explicitly.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to resolve; defaults to `self.charge_evaluator`
        :type evaluator: object | None
        :param opts: extra options forwarded to the evaluator's `from_mol` constructor, if applicable
        :type opts: dict
        :return: the resolved charge-evaluator instance
        :rtype: object
        :raises ValueError: if the evaluator type can't be resolved
        """
        if evaluator is None:
            evaluator = self.charge_evaluator
        eval_type, new_opts = ChargeEvaluator.resolve_evaluator(evaluator)
        if eval_type is None:
            raise ValueError(f"can't resolve charge evaluator type for {evaluator}")
        if hasattr(eval_type, 'from_mol') and isinstance(eval_type, type):
            return eval_type.from_mol(self, **dict(opts, **new_opts))
        else:
            return eval_type
    def calculate_charges(self, evaluator=None, order=None, **opts):
        """
        **LLM Docstring**

        Compute partial-charge values (and, optionally, their derivatives) using the resolved charge evaluator at this molecule's current geometry.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use; defaults to `self.charge_evaluator`
        :type evaluator: object | None
        :param order: the highest derivative order to compute; if `None`, only the charges themselves are returned
        :type order: int | None
        :param opts: extra options forwarded to `get_charge_evaluator`
        :type opts: dict
        :return: the charges (if `order` is `None`) or the full charge/derivative expansion
        :rtype: np.ndarray | list[np.ndarray]
        """
        evaluator = self.get_charge_evaluator(evaluator, **opts)
        smol = order is None
        if smol: order = 0
        expansion = evaluator.evaluate(
            self.coords * UnitsData.convert("BohrRadius", evaluator.distance_units),
            order=order
        )
        if smol: expansion = expansion[0]
        return expansion
    @internals.setter
    def internals(self, spec):
        """
        **LLM Docstring**

        Getter for the raw (canonicalized) internal-coordinate specification, delegating to `self.embedding.internals`. (A setter with the same name separately rebuilds the embedding from a new specification.)

        :return: the internal-coordinate specification, or `None` if none is set
        :rtype: dict | None
        """
        self.embedding = MolecularEmbedding(
            self.atomic_masses,
            self.coords,
            self.canonicalize_internals(spec, self.atoms, self.coords, self._bonds, masses=self._mass)
        )
    @property
    def internal_coordinates(self):
        """
        **LLM Docstring**

        Getter for the internal coordinates at the molecule's current geometry, delegating to `self.embedding.internal_coordinates`.

        :return: the internal coordinates, or `None` if none are defined
        :rtype: CoordinateSet | None
        """
        return self.embedding.internal_coordinates
    @property
    def redundant_internal_transformation(self):
        """
        **LLM Docstring**

        Getter for the redundant-to-non-redundant internal-coordinate transformation, delegating to `self.embedding.redundant_internal_transformation`.

        :return: the redundant transformation, or `None` if not applicable
        :rtype: np.ndarray | None
        """
        return self.embedding.redundant_internal_transformation

    @classmethod
    def _check_label(cls, label,
                     allowed_coordinate_types=None,
                     excluded_coordinate_types=None,
                     allowed_ring_types=None,
                     excluded_ring_types=None,
                     allowed_group_types=None,
                     excluded_group_types=None,
                     ):
        """
        **LLM Docstring**

        Test whether a coordinate label passes a set of allow/exclude filters on its atom types, ring membership, and functional-group membership.

        :param label: the coordinate label to test (exposing `.atoms`, `.ring`, `.group` attributes)
        :type label: object
        :param allowed_coordinate_types: if given, `label.atoms` must be among these to pass
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: if given, `label.atoms` must not be among these to pass
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: if given, `label.ring` must be among these to pass
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: if given, `label.ring` must not be among these to pass
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: if given, `label.group` must be among these to pass
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: if given, `label.group` must not be among these to pass
        :type excluded_group_types: Iterable | None
        :return: whether the label passes every specified filter
        :rtype: bool
        """
        if allowed_coordinate_types is not None:
            if label.atoms not in allowed_coordinate_types: return False
        if excluded_coordinate_types is not None:
            if label.atoms in excluded_coordinate_types: return False
        if allowed_ring_types is not None:
            if label.ring not in allowed_ring_types: return False
        if excluded_ring_types is not None:
            if label.ring in excluded_ring_types: return False
        if allowed_group_types is not None:
            if label.group not in allowed_group_types: return False
        if excluded_group_types is not None:
            if label.group in excluded_group_types: return False
        return True

    @classmethod
    def get_coordinate_filer(cls,
                             allowed_coordinate_types=None,
                             excluded_coordinate_types=None,
                             allowed_ring_types=None,
                             excluded_ring_types=None,
                             allowed_group_types=None,
                             excluded_group_types=None
                             ):
        """
        **LLM Docstring**

        Build a filter function (closing over the given allow/exclude criteria) that, given a dict of coordinate-to-label mappings, returns only the entries whose label passes `_check_label`.

        :param allowed_coordinate_types: forwarded to `_check_label`
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: forwarded to `_check_label`
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: forwarded to `_check_label`
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: forwarded to `_check_label`
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: forwarded to `_check_label`
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: forwarded to `_check_label`
        :type excluded_group_types: Iterable | None
        :return: the constructed coordinate-filtering function
        :rtype: callable
        """
        def coordinate_filter(coords):
            """
            **LLM Docstring**

            Filter a dict of coordinate-to-label mappings down to just the entries whose label satisfies the enclosing allow/exclude criteria (via `_check_label`).

            :param coords: the coordinate-to-label mapping to filter
            :type coords: dict
            :return: the filtered mapping
            :rtype: dict
            """
            return {
                c: l
                for c, l in coords.items()
                if cls._check_label(l,
                                    allowed_coordinate_types=allowed_coordinate_types,
                                    excluded_coordinate_types=excluded_coordinate_types,
                                    allowed_ring_types=allowed_ring_types,
                                    excluded_ring_types=excluded_ring_types,
                                    allowed_group_types=allowed_group_types,
                                    excluded_group_types=excluded_group_types
                                    )
            }

        return coordinate_filter

    default_coordinate_pruning = 'graph'
    def get_bond_graph_internals(self,
                                 include_stretches=True,
                                 include_bends=True,
                                 include_dihedrals=True,
                                 include_fragments=True,
                                 pruning=None,
                                 fragment=None,
                                 base_internals=None,
                                 use_distance_matrix=True,
                                 concatenate=True
                                 ):
        """
        **LLM Docstring**

        Build a set of internal coordinates (bond stretches, bends, dihedrals, and/or inter-fragment coordinates) directly from the bonding graph, optionally restricted to a single fragment (recursively, with the result permuted back into the full atom indexing) and/or pruned down to a well-conditioned subset.

        :param include_stretches: whether to include bond-stretch coordinates
        :type include_stretches: bool
        :param include_bends: whether to include bond-angle coordinates
        :type include_bends: bool
        :param include_dihedrals: whether to include dihedral-angle coordinates
        :type include_dihedrals: bool
        :param include_fragments: whether to include coordinates connecting separate molecular fragments
        :type include_fragments: bool
        :param pruning: whether/how to prune the resulting coordinates (`True` for the default method, or an explicit method spec), forwarded to `prune_internals`
        :type pruning: bool | str | dict | None
        :param fragment: restrict to a single fragment, given as a fragment index or an explicit list of atom indices
        :type fragment: int | Iterable[int] | None
        :param base_internals: accepted and forwarded when recursing on a fragment, but not otherwise used directly in this method's own body
        :type base_internals: object | None
        :param use_distance_matrix: whether to precompute a distance matrix for the fragment-coordinate generation
        :type use_distance_matrix: bool
        :param concatenate: whether to concatenate the different coordinate categories (stretches/bends/dihedrals/fragments) into a single list, or return them as separate groups
        :type concatenate: bool
        :return: the generated internal coordinates, as a single concatenated list or a list of category groups depending on `concatenate`
        :rtype: list
        :raises ValueError: if `pruning` is requested while `concatenate` is `False`
        """
        if fragment is not None:
            if nput.is_int(fragment):
                fragment = self.fragment_indices[fragment]
            base_ints = self.take_submolecule(fragment).get_bond_graph_internals(
                include_stretches=include_stretches,
                include_bends=include_bends,
                include_dihedrals=include_dihedrals,
                include_fragments=include_fragments,
                base_internals=base_internals,
                pruning=pruning,
                concatenate=concatenate
            )
            if concatenate:
                return coordops.permute_internals(base_ints, fragment)
            else:
                return [
                    coordops.permute_internals(b, fragment)
                    for b in base_ints
                ]
        else:
            st, bo, di = coordops.get_stretch_coordinate_system(
                [tuple(b[:2]) for b in self.bonds],
                include_bends=include_bends,
                include_dihedrals=include_dihedrals
            )
            bits = []
            if include_fragments:
                if use_distance_matrix:
                    dm = nput.distance_matrix(self.coords)
                else:
                    dm = None
                frag_bits = coordops.get_fragment_coordinate_system(
                    self.edge_graph,
                    masses=self.masses,
                    distance_matrix=dm
                )
                bits.append(frag_bits)
            if include_stretches:
                bits.append(st)
            if include_bends:
                bits.append(bo)
            if include_dihedrals:
                bits.append(di)

            if concatenate:
                internals = bits[0]
                for b in bits[1:]:
                    internals = internals + b

                if pruning:
                    if pruning is True:
                        pruning = self.default_coordinate_pruning
                    internals = self.prune_internals(internals, method=pruning)
            else:
                if pruning:
                    raise ValueError("can't prune without concatenating")
                internals = bits

            return internals
    def prune_internals(self, coords, method='b_matrix', check_rigidity=True):
        """
        **LLM Docstring**

        Reduce a set of internal coordinates down to a non-redundant, well-conditioned subset, defaulting to a B-matrix-rank-based method (building the necessary translation/rotation-projected B-matrix generator and a sensible `max_coords` cap) if no custom method is supplied.

        :param coords: the internal-coordinate specs to prune
        :type coords: list
        :param method: the pruning method: a method-name string, or a dict of method options (with a `'method'` key defaulting to `'b_matrix'`)
        :type method: str | dict
        :param check_rigidity: whether to check that the pruned coordinate set spans a rigid (non-redundant) representation
        :type check_rigidity: bool
        :return: the pruned coordinate specs
        :rtype: list
        """
        if isinstance(method, str):
            method = {'method':method}
        if hasattr(method, 'items'):
            meth = method.get('method')
            if meth is None:
                method = method.copy()
                method['method'] = 'b_matrix'
                meth = 'b_matrix'
            if dev.str_is(meth, 'b_matrix'):
                if 'b_matrix' not in method:
                    g12 = self.get_gmatrix(power=1 / 2)
                    proj = nput.translation_rotation_projector(self.coords, self.atomic_masses, mass_weighted=True)
                    def b_gen(pos, crds):
                        """
                        **LLM Docstring**

                        Compute the (translation/rotation-projected, mass-weighted) B-matrix for a candidate set of internal coordinates at this molecule's current geometry, used by the default `'b_matrix'` pruning method to assess rank/conditioning.

                        :param pos: the coordinate index/indices under consideration (unused directly in the body, but part of the callback signature expected by the pruning routine)
                        :type pos: object
                        :param crds: the candidate coordinate specs to build the B-matrix for
                        :type crds: list
                        :return: the projected, mass-weighted B-matrix
                        :rtype: np.ndarray
                        """
                        return proj @ g12 @ nput.internal_coordinate_tensors(self.coords, crds, order=1)[1]
                    method = method.copy()
                    method['b_matrix'] = b_gen
                if 'max_coords' not in method:
                    method = method.copy()
                    method['max_coords'] =  min([3 * len(self.atoms) - 6, len(coords)])
        return coordops.prune_internal_coordinates(
            coords,
            method=method,
            check_rigidity=check_rigidity,
        )
    def get_labeled_internals(self,
                              coordinate_filter=None,
                              allowed_coordinate_types=None,
                              excluded_coordinate_types=None,
                              allowed_ring_types=None,
                              excluded_ring_types=None,
                              allowed_group_types=None,
                              excluded_group_types=None,
                              include_stretches=True,
                              include_bends=True,
                              include_dihedrals=True,
                              include_fragments=True,
                              coordinate_sorting=None,
                              pruning=False
                              ):
        """
        **LLM Docstring**

        Build the internal coordinates from the bonding graph (via `get_bond_graph_internals`) and label each one by its atom types/ring/functional-group membership (via `edge_graph.get_label_types` and `coordops.get_coordinate_label`), then filter and sort them.

        :param coordinate_filter: an explicit filter function to apply instead of building one from the allow/exclude arguments
        :type coordinate_filter: callable | None
        :param allowed_coordinate_types: forwarded to `get_coordinate_filer` if `coordinate_filter` is not given
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: forwarded to `get_coordinate_filer`
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: forwarded to `get_coordinate_filer`
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: forwarded to `get_coordinate_filer`
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: forwarded to `get_coordinate_filer`
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: forwarded to `get_coordinate_filer`
        :type excluded_group_types: Iterable | None
        :param include_stretches: whether to include bond-stretch coordinates
        :type include_stretches: bool
        :param include_bends: whether to include bond-angle coordinates
        :type include_bends: bool
        :param include_dihedrals: whether to include dihedral-angle coordinates
        :type include_dihedrals: bool
        :param include_fragments: whether to include inter-fragment coordinates
        :type include_fragments: bool
        :param coordinate_sorting: a custom sorting function to apply to the labeled coordinates instead of the default `coordops.sort_internal_coordinates`; pass a falsy value to skip sorting
        :type coordinate_sorting: callable | bool | None
        :param pruning: whether/how to prune the coordinates, forwarded to `get_bond_graph_internals`
        :type pruning: bool | str | dict
        :return: a mapping from coordinate spec to its label, filtered and sorted
        :rtype: dict
        """
        internals = self.get_bond_graph_internals(
            include_stretches=include_stretches,
            include_bends=include_bends,
            include_dihedrals=include_dihedrals,
            include_fragments=include_fragments,
            pruning=pruning
        )

        labels = self.edge_graph.get_label_types()
        internals = {
            (c if isinstance(c, tuple) else coordops.InternalCoordinateType.resolve(c)): coordops.get_coordinate_label(
                c,
                labels
            )
            for c in internals
        }

        if coordinate_filter is None:
            coordinate_filter = self.get_coordinate_filer(
                allowed_coordinate_types=allowed_coordinate_types,
                excluded_coordinate_types=excluded_coordinate_types,
                allowed_ring_types=allowed_ring_types,
                excluded_ring_types=excluded_ring_types,
                allowed_group_types=allowed_group_types,
                excluded_group_types=excluded_group_types,
            )

        if coordinate_filter:
            internals = coordinate_filter(internals)

        if coordinate_sorting is None:
            coordinate_sorting = coordops.sort_internal_coordinates

        if coordinate_sorting:
            internals = coordinate_sorting(internals)

        return internals

    def get_mode_labels(self,
                        internals=None,
                        modes=None,
                        use_redundants=True,
                        expansions=None,
                        return_modes=False,
                        **internals_opts
                        ):
        """
        **LLM Docstring**

        Assign human-readable labels (e.g. "C-H stretch") to a set of normal modes by projecting them onto labeled internal coordinates, handling both redundant and non-redundant internal-coordinate expansions and both Cartesian- and internal-coordinate-basis modes.

        :param internals: the labeled internal coordinates to project onto; computed via `get_labeled_internals` if not given
        :type internals: dict | None
        :param modes: the normal modes to label; computed via `get_normal_modes` if not given
        :type modes: object | None
        :param use_redundants: whether to build a redundant-coordinate expansion (with relocalization) for the projection, rather than using the internal coordinates directly
        :type use_redundants: bool
        :param expansions: precomputed `(expansions, inverse_expansion)` internal-coordinate Jacobian data to reuse instead of recomputing it
        :type expansions: tuple | None
        :param return_modes: whether to also return the internal-coordinate-basis mode matrix alongside the labels
        :type return_modes: bool
        :param internals_opts: extra options forwarded to `get_labeled_internals` if `internals` is not given
        :type internals_opts: dict
        :return: the mode labels, or `(internal_modes, labels)` if `return_modes` is set
        :rtype: list | tuple
        """
        if modes is None:
            modes = self.get_normal_modes()
        modes = modes.remove_mass_weighting()

        if internals is None:
            internals = self.get_labeled_internals(**internals_opts)

        if modes.is_cartesian:
            if expansions is not None:
                expansions, inv_expansion = expansions
            else:
                expansions = inv_expansion = None

            if use_redundants:
                redundant_tf, expansions = coordops.RedundantCoordinateGenerator(
                    internals,
                    masses=self.atomic_masses,
                    relocalize=True
                ).compute_redundant_expansions(self.coords,
                                               expansions=expansions
                                               )

                redund_labs = coordops.get_mode_labels(
                    internals,
                    redundant_tf,
                    norm_cutoff=.3
                )

                inv_expansion = nput.inverse_internal_coordinate_tensors(
                    expansions,
                    coords=self.coords,
                    masses=self.atomic_masses,
                    order=1,
                    remove_translation_rotation=True
                )

            else:
                redund_labs = internals
                if expansions is None:
                    expansions, inv_expansion = nput.internal_coordinate_tensors(
                        self.coords,
                        internals,
                        order=1,
                        masses=self.atomic_masses,
                        return_inverse=True
                    )
                    expansions = expansions[1:]

            g = expansions[0].T @ self.get_gmatrix() @ expansions[0]
            g12 = nput.fractional_power(g, 1 / 2)
            internal_modes = g12 @ inv_expansion[0] @ modes.modes_by_coords
        else:
            redund_labs = internals
            internal_modes = modes

        labels = coordops.get_mode_labels(
            redund_labs,
            internal_modes,
            norm_cutoff=.8
        )
        if return_modes:
            return internal_modes, labels
        else:
            return labels

    @property
    def mode_embedding(self):
        """
        **LLM Docstring**

        The (cached) `ModeEmbedding` combining this molecule's coordinate embedding with its normal modes, built lazily the first time it's needed.

        :return: the mode embedding
        :rtype: ModeEmbedding
        """
        if self._mode_embedding is None:
            self._mode_embedding = ModeEmbedding(self.embedding, self.normal_modes)
        return self._mode_embedding
    def get_internals(self, coords=None, *, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch internal coordinates, either the molecule's own cached ones or those for an alternate set of Cartesian `coords`, via `self.embedding.get_internals`.

        :param coords: alternate Cartesian coordinates to convert instead of using the cached internal coordinates
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the internal coordinates, or `None` if none are defined
        :rtype: CoordinateSet | None
        """
        return self.embedding.get_internals(coords=coords, strip_embedding=strip_embedding)

    def get_cartesians_by_internals(self, order=None, coords=None, *, strip_embedding=False, **kw):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-internals Jacobian expansion, via `self.embedding.get_cartesians_by_internals`.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param coords: alternate coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param kw: extra options forwarded to the embedding
        :type kw: dict
        :return: the Cartesians-by-internals Jacobian tensors
        :rtype: list[np.ndarray]
        """
        return self.embedding.get_cartesians_by_internals(coords=coords, order=order, strip_embedding=strip_embedding, **kw)

    def get_internals_by_cartesians(self, order=None, *, coords=None, strip_embedding=False, **kw):
        """
        **LLM Docstring**

        Fetch the internals-by-Cartesians Jacobian expansion, via `self.embedding.get_internals_by_cartesians`.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param coords: alternate coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param kw: extra options forwarded to the embedding
        :type kw: dict
        :return: the internals-by-Cartesians Jacobian tensors
        :rtype: list[np.ndarray]
        """
        return self.embedding.get_internals_by_cartesians(order=order, coords=coords, strip_embedding=strip_embedding, **kw)

    def get_cartesians_by_modes(self, order=None, **kw):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-normal-modes Jacobian expansion, via `self.mode_embedding.get_cartesians_by_internals`.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param kw: extra options forwarded to the mode embedding
        :type kw: dict
        :return: the Cartesians-by-modes Jacobian tensors
        :rtype: list[np.ndarray]
        """
        return self.mode_embedding.get_cartesians_by_internals(order=order, **kw)

    def get_modes_by_cartesians(self, order=None, **kw):
        """
        **LLM Docstring**

        Fetch the normal-modes-by-Cartesians Jacobian expansion, via `self.mode_embedding.get_internals_by_cartesians`.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param kw: extra options forwarded to the mode embedding
        :type kw: dict
        :return: the modes-by-Cartesians Jacobian tensors
        :rtype: list[np.ndarray]
        """
        return self.mode_embedding.get_internals_by_cartesians(order=order, **kw)

    #endregion

    #region Properties
    @property
    def dipole_surface(self):
        """
        :return:
        :rtype: DipoleSurfaceManager
        """
        return self._dips
    @dipole_surface.setter
    def dipole_surface(self, val):
        """
        **LLM Docstring**

        Setter for the `DipoleSurfaceManager`, requiring the new value to already be one. (The getter simply returns `self._dips`.)

        :param val: the new dipole-surface manager
        :type val: DipoleSurfaceManager
        :return: None
        :rtype: None
        :raises TypeError: if `val` isn't a `DipoleSurfaceManager`
        """
        if not isinstance(val, DipoleSurfaceManager):
            raise TypeError("`dipole_surface` must be {}".format(
                DipoleSurfaceManager.__name__
            ))
        self._dips = val
    @property
    def dipole_derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the dipole derivative tensors. The getter delegates to `self.dipole_surface.get_derivatives(quiet=True)`; the setter assigns to `self.dipole_surface.derivatives`.

        :param derivs: (setter only) the new dipole derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the dipole derivative tensors, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        return self.dipole_surface.get_derivatives(quiet=True)
    @dipole_derivatives.setter
    def dipole_derivatives(self, derivs):
        """
        **LLM Docstring**

        Property getter/setter for the dipole derivative tensors. The getter delegates to `self.dipole_surface.get_derivatives(quiet=True)`; the setter assigns to `self.dipole_surface.derivatives`.

        :param derivs: (setter only) the new dipole derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the dipole derivative tensors, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        self.dipole_surface.derivatives = derivs
    def get_cartesian_dipole_derivatives(self, order=None, evaluator=None, include_constant_term=False):
        """
        **LLM Docstring**

        Fetch the dipole derivatives in Cartesian coordinates, computing them via `calculate_dipole` (and caching the result on the molecule, if the same evaluator is configured as the default) if not already available to the requested order.

        :param order: the highest derivative order needed; if `None`, whatever is available is returned
        :type order: int | None
        :param evaluator: an explicit dipole evaluator to use instead of `self.dipole_evaluator`
        :type evaluator: object | None
        :param include_constant_term: whether to include the zeroth-order (reference dipole) term in the result
        :type include_constant_term: bool
        :return: the dipole derivative tensors (from first order, or zeroth if `include_constant_term`), or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        dipole_derivatives = self.dipole_derivatives
        # if dipole_derivatives is None and self.energy_evaluator is not None:
        #     if order is None: order = 1
        #     dipole_derivatives = self.calculate_dipole(order=order)
        if dipole_derivatives is None or (order is not None and len(dipole_derivatives) < order):
            if evaluator is None: evaluator = self.dipole_evaluator
            if evaluator is not None:
                if order is None: order = 1
                opts = dict(evaluator=evaluator, order=order)
                if dev.str_is(evaluator, 'expansion'):
                    opts['use_modes'] = False
                dipole_derivatives = self.calculate_dipole(**opts)
                if (
                        isinstance(evaluator, str)
                        and isinstance(self.dipole_evaluator, str)
                        and evaluator == self.dipole_evaluator
                ) or evaluator is self.dipole_evaluator:
                    self.dipole_derivatives = dipole_derivatives
        o = 0 if include_constant_term else 1
        if dipole_derivatives is None:
            return dipole_derivatives
        if order is None:
            return dipole_derivatives[o:]
        else:
            return dipole_derivatives[o:order+1]
    def get_internal_dipole_derivatives(self, order=None, reembed=True, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the dipole derivatives re-expressed in internal coordinates, by re-expanding the Cartesian dipole derivatives through the Cartesians-by-internals Jacobian.

        :param order: the highest derivative order needed
        :type order: int | None
        :param reembed: whether to use the Eckart-reembedded Cartesians-by-internals Jacobian
        :type reembed: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the Jacobian
        :type strip_embedding: bool
        :return: the internal-coordinate dipole derivative tensors
        :rtype: list[np.ndarray]
        """
        derivs = self.get_cartesian_dipole_derivatives(order=order)
        if order is None:
            order = len(derivs)
        return nput.tensor_reexpand(
            self.get_cartesians_by_internals(order, reembed=reembed, strip_embedding=strip_embedding),
            derivs,
            order
        )
    def get_cartesian_polarizability_derivatives(self, order=None, evaluator=None, include_constant_term=False):
        """
        **LLM Docstring**

        Fetch the dipole-polarizability derivatives in Cartesian coordinates, computing them via `calculate_dipole_polarizability` (and caching the result, if the same evaluator is configured as the default) if not already available to the requested order, or re-expanding them through the normal modes if they were stored in a mode basis smaller than the full Cartesian space.

        :param order: the highest derivative order needed; if `None`, whatever is available is returned
        :type order: int | None
        :param evaluator: an explicit polarizability evaluator to use instead of `self.polarizability_evaluator`
        :type evaluator: object | None
        :param include_constant_term: accepted for interface consistency with `get_cartesian_dipole_derivatives` but not used in this method's body
        :type include_constant_term: bool
        :return: the polarizability derivative tensors
        :rtype: list[np.ndarray]
        """
        derivs = self.polarizability_derivatives
        # if dipole_derivatives is None and self.energy_evaluator is not None:
        #     if order is None: order = 1
        #     dipole_derivatives = self.calculate_dipole(order=order)
        if derivs is None or (order is not None and len(derivs) < order):
            if evaluator is None: evaluator = self.polarizability_evaluator
            if evaluator is not None:
                if order is None: order = 1
                opts = dict(evaluator=evaluator, order=order)
                # if dev.str_is(evaluator, 'expansion'):
                #     opts['use_modes'] = False
                derivs = self.calculate_dipole_polarizability(**opts)[1]
                if (
                        isinstance(evaluator, str)
                        and isinstance(self.polarizability_evaluator, str)
                        and evaluator == self.polarizability_evaluator
                ) or evaluator is self.polarizability_evaluator:
                    self.polarizability_derivatives = derivs
        elif len(derivs) > 1 and derivs[1].shape[0] < (3*len(self._ats)):
            derivs = [derivs[0]] + nput.tensor_reexpand(
                [self.get_normal_modes(use_internals=False, project_transrot=False).modes_by_coords],
                derivs[1:],
                order=len(derivs) - 1
            )

        return derivs

    def get_hamiltonian(self,
                        embedding=None,
                        potential_derivatives=None,
                        modes=None,
                        dipole_derivatives=None,
                        **etc
                        ):
        """
        **LLM Docstring**

        Build a `MolecularHamiltonian` for this molecule, defaulting to its own embedding, potential derivatives, normal modes, and dipole derivatives wherever not explicitly overridden.

        :param embedding: an alternate coordinate embedding to use
        :type embedding: MolecularEmbedding | None
        :param potential_derivatives: alternate potential-energy derivative tensors to use
        :type potential_derivatives: list[np.ndarray] | None
        :param modes: alternate normal modes to use
        :type modes: object | None
        :param dipole_derivatives: alternate dipole derivative tensors to use
        :type dipole_derivatives: list[np.ndarray] | None
        :param etc: extra options forwarded to the `MolecularHamiltonian` constructor
        :type etc: dict
        :return: the constructed Hamiltonian
        :rtype: MolecularHamiltonian
        """
        if embedding is None:
            embedding = self.embedding
        if potential_derivatives is None:
            potential_derivatives = self.potential_derivatives
        if modes is None:
            modes = self.normal_modes
        if dipole_derivatives is None:
            dipole_derivatives = self.dipole_derivatives
        return MolecularHamiltonian(embedding,
                                    potential_derivatives=potential_derivatives,
                                    modes=modes,
                                    dipole_derivatives=dipole_derivatives,
                                    **etc
                                    )
    @property
    def hamiltonian(self):
        """
        **LLM Docstring**

        Property getter/setter for the `MolecularHamiltonian`. The getter lazily builds one via `get_hamiltonian` the first time it's needed.

        :param e: (setter only) the new Hamiltonian
        :type e: MolecularHamiltonian
        :return: (getter) the cached (or newly built) Hamiltonian
        :rtype: MolecularHamiltonian
        """
        if self._hamiltonian is None:
            self._hamiltonian = self.get_hamiltonian()
        return self._hamiltonian
    @hamiltonian.setter
    def hamiltonian(self, e):
        """
        **LLM Docstring**

        Property getter/setter for the `MolecularHamiltonian`. The getter lazily builds one via `get_hamiltonian` the first time it's needed.

        :param e: (setter only) the new Hamiltonian
        :type e: MolecularHamiltonian
        :return: (getter) the cached (or newly built) Hamiltonian
        :rtype: MolecularHamiltonian
        """
        self._hamiltonian = e

    @property
    def potential_surface(self):
        """
        :return:
        :rtype: PotentialSurfaceManager
        """
        return self._pes
    @potential_surface.setter
    def potential_surface(self, val):
        """
        **LLM Docstring**

        Setter for the `PotentialSurfaceManager`, requiring the new value to already be one. (The getter simply returns `self._pes`.)

        :param val: the new potential-surface manager
        :type val: PotentialSurfaceManager
        :return: None
        :rtype: None
        :raises TypeError: if `val` isn't a `PotentialSurfaceManager`
        """
        if not isinstance(val, PotentialSurfaceManager):
            raise TypeError("`potential_surface` must be {}".format(
                PotentialSurfaceManager.__name__
            ))
        self._pes = val
    @property
    def potential_derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the potential-energy derivative tensors. The getter fetches them from `self.potential_surface`, normalizing missing/placeholder entries to `0` and trimming any trailing zero-padded (unset) higher-order terms; the setter assigns to `self.potential_surface.derivatives`.

        :param derivs: (setter only) the new potential derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the potential derivative tensors with trailing unset orders trimmed, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        base_derivs = self.potential_surface.get_derivs(quiet=True)
        if base_derivs is None:
            return None
        base_derivs = [
            (
                (0 if nput.is_numeric(p) else np.asanyarray(p))
                if p is not None else 0
            )
            for p in base_derivs
        ]
        n = len(base_derivs)
        for i in range(n): # remove 0 padding
            if isinstance(base_derivs[n - (i+1)], np.ndarray):
                n = n - i
                break
        else:
            n = 0
        return base_derivs[:n]
    @potential_derivatives.setter
    def potential_derivatives(self, derivs):
        """
        **LLM Docstring**

        Property getter/setter for the potential-energy derivative tensors. The getter fetches them from `self.potential_surface`, normalizing missing/placeholder entries to `0` and trimming any trailing zero-padded (unset) higher-order terms; the setter assigns to `self.potential_surface.derivatives`.

        :param derivs: (setter only) the new potential derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the potential derivative tensors with trailing unset orders trimmed, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        self.potential_surface.derivatives = derivs

    def get_cartesian_potential_derivatives(self, order=None, evaluator=None, use_cached=True):
        """
        **LLM Docstring**

        Fetch the potential-energy derivatives in Cartesian coordinates, computing them via `calculate_energy` (and caching the result on the molecule, if the same evaluator is configured as the default) if not already available to the requested order.

        :param order: the highest derivative order needed; defaults to `2` if a fresh calculation is required
        :type order: int | None
        :param evaluator: an explicit energy evaluator to use instead of `self.energy_evaluator`
        :type evaluator: object | None
        :param use_cached: accepted for interface consistency but not used in this method's body
        :type use_cached: bool
        :return: the potential derivative tensors, truncated to `order` if given, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        potential_derivatives = self.potential_derivatives
        if potential_derivatives is None or (order is not None and len(potential_derivatives) < order):
            if evaluator is None: evaluator = self.energy_evaluator
            if evaluator is not None:
                if order is None: order = 2
                potential_derivatives = self.calculate_energy(evaluator=evaluator, order=order)[1:]
                if (
                        isinstance(evaluator, str)
                        and isinstance(self.energy_evaluator, str)
                        and evaluator == self.energy_evaluator
                ) or evaluator is self.energy_evaluator:
                    self.potential_derivatives = potential_derivatives
        if order is None or potential_derivatives is None:
            return potential_derivatives
        else:
            return potential_derivatives[:order]
    def get_internal_potential_derivatives(self, order=None, reembed=True, strip_embedding=True, zero_gradient=False):
        """
        **LLM Docstring**

        Fetch the potential-energy derivatives re-expressed in internal coordinates, by re-expanding the Cartesian potential derivatives through the Cartesians-by-internals Jacobian.

        :param order: the highest derivative order needed
        :type order: int | None
        :param reembed: whether to use the Eckart-reembedded Cartesians-by-internals Jacobian
        :type reembed: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the Jacobian
        :type strip_embedding: bool
        :param zero_gradient: whether to zero out the first-order (gradient) term before re-expanding
        :type zero_gradient: bool
        :return: the internal-coordinate potential derivative tensors, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        derivs = self.get_cartesian_potential_derivatives(order)
        if derivs is None: return None
        if zero_gradient:
            derivs = [0] + list(derivs[1:])
        if order is None:
            order = len(derivs)
        return nput.tensor_reexpand(
            self.get_cartesians_by_internals(order, reembed=reembed, strip_embedding=strip_embedding),
            derivs,
            order
        )

    @property
    def normal_modes(self):
        """
        :return:
        :rtype: NormalModesManager
        """
        return self._normal_modes
    @normal_modes.setter
    def normal_modes(self, val):
        """
        **LLM Docstring**

        Setter for the normal modes: wraps the given value in a `NormalModesManager` (via `NormalModesManager.from_data`) unless it already is one. (A separate getter, not shown here, returns the current modes.)

        :param val: the new normal modes, in any form accepted by `NormalModesManager.from_data`
        :type val: object
        :return: None
        :rtype: None
        """
        if not isinstance(val, NormalModesManager):
            val = NormalModesManager.from_data(self, val)
            # raise TypeError("`normal_modes` must be {}".format(
            #     NormalModesManager.__name__
            # ))
        self._normal_modes = val
    def get_normal_modes(self, masses=None,
                         potential_derivatives=None,
                         use_internals=None,
                         project_transrot=True,
                         **opts):
        """
        **LLM Docstring**

        Compute this molecule's normal modes (via `NormalModes.from_molecule`), optionally using alternate masses/potential derivatives and controlling whether internal coordinates and translation/rotation projection are used.

        :param masses: masses to use instead of this molecule's own
        :type masses: np.ndarray | None
        :param potential_derivatives: potential derivatives to use instead of this molecule's own
        :type potential_derivatives: list[np.ndarray] | None
        :param use_internals: whether to compute the modes in internal coordinates rather than Cartesians
        :type use_internals: bool | None
        :param project_transrot: whether to project out translational/rotational degrees of freedom
        :type project_transrot: bool
        :param opts: extra options forwarded to `NormalModes.from_molecule`
        :type opts: dict
        :return: the computed normal modes
        :rtype: NormalModes
        """
        from ..Modes import NormalModes
        return NormalModes.from_molecule(self,
                                         masses=masses,
                                         potential_derivatives=potential_derivatives,
                                         use_internals=use_internals,
                                         project_transrot=project_transrot,
                                         **opts
                                         )
    def get_reaction_path_modes(self, masses=None,
                                potential_derivatives=None,
                                **opts):
        """
        **LLM Docstring**

        Compute reaction-path-following normal modes for this molecule, via `ReactionPathModes.from_molecule`.

        :param masses: masses to use instead of this molecule's own
        :type masses: np.ndarray | None
        :param potential_derivatives: potential derivatives to use instead of this molecule's own
        :type potential_derivatives: list[np.ndarray] | None
        :param opts: extra options forwarded to `ReactionPathModes.from_molecule`
        :type opts: dict
        :return: the computed reaction-path modes
        :rtype: ReactionPathModes
        """
        from ..Modes import ReactionPathModes
        return ReactionPathModes.from_molecule(self,
                                               masses=masses,
                                               potential_derivatives=potential_derivatives,
                                               **opts)
    @property
    def metadata(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's metadata dict. The setter requires the new value to already be a `dict`.

        :param val: (setter only) the new metadata dict
        :type val: dict
        :return: (getter) the metadata dict
        :rtype: dict
        :raises TypeError: if the setter is given something that isn't a `dict`
        """
        return self._meta
    @metadata.setter
    def metadata(self, val):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's metadata dict. The setter requires the new value to already be a `dict`.

        :param val: (setter only) the new metadata dict
        :type val: dict
        :return: (getter) the metadata dict
        :rtype: dict
        :raises TypeError: if the setter is given something that isn't a `dict`
        """
        if not isinstance(val, dict):
            raise TypeError("metadata must be {}".format(
                dict.__name__
            ))
        self._meta = val


    def get_harmonic_spectrum(self, **opts):
        """
        **LLM Docstring**

        Build a harmonic IR spectrum for this molecule, via `HarmonicSpectrum.from_mol`.

        :param opts: extra options forwarded to `HarmonicSpectrum.from_mol`
        :type opts: dict
        :return: the constructed harmonic spectrum
        :rtype: HarmonicSpectrum
        """
        from ..Spectra import HarmonicSpectrum
        return HarmonicSpectrum.from_mol(self, **opts)

    def get_harmonic_raman_spectrum(self, **opts):
        """
        **LLM Docstring**

        Build a harmonic Raman spectrum for this molecule, via `HarmonicSpectrum.raman_from_mol`.

        :param opts: extra options forwarded to `HarmonicSpectrum.raman_from_mol`
        :type opts: dict
        :return: the constructed harmonic Raman spectrum
        :rtype: HarmonicSpectrum
        """
        from ..Spectra import HarmonicSpectrum
        return HarmonicSpectrum.raman_from_mol(self, **opts)
    #endregion

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation using the molecule's name (falling back to the source-file basename, then an SMILES string from any attached RDKit molecule, then the chemical formula) and its coordinate shape.

        :return: string of the form `ClassName('name', shape=coord_shape)`
        :rtype: str
        """
        name = self._name
        if name is None:
            src = self.source_file
            if src is not None:
                name = os.path.basename(src)
        if name is None:
            rdmol = self._rdmol
            if rdmol is not None:
                name = rdmol.to_smiles(remove_hydrogens=True)
        if name is None:
            name = self.formula

        return "{cls}('{name}', shape={shape})".format(
            cls=type(self).__name__,
            name=name,
            # formula=self.formula,
            shape=self.coords.shape,
            # coord_sys=self.coords.system.name if isinstance(self.coords, CoordinateSet) else 'undefined'
        )

    @property
    def num_atoms(self):
        """
        **LLM Docstring**

        The number of atoms in the molecule.

        :return: the atom count
        :rtype: int
        """
        return len(self._ats)
    @property
    def atom_positions(self):
        """
        **LLM Docstring**

        A mapping from element symbol to the list of atom indices having that symbol.

        :return: the symbol-to-indices mapping
        :rtype: dict
        """
        pos_map = {}
        for i,a in enumerate(self._ats):
            if a["Symbol"] in pos_map:
                pos_map[a["Symbol"]].append(i)
            else:
                pos_map[a["Symbol"]] = [i]
        return pos_map
    @property
    def dummy_positions(self):
        """
        **LLM Docstring**

        The indices of any dummy (`'X'`) atoms in the molecule.

        :return: the dummy-atom indices, or an empty list if there are none
        :rtype: list[int]
        """
        ats = self.atom_positions
        return ats['X'] if 'X' in ats else []
    @property
    def atoms(self):
        """
        **LLM Docstring**

        The element symbols of the molecule's atoms.

        :return: the atom symbols
        :rtype: tuple[str]
        """
        return tuple(a["Symbol"] for a in self._ats)
    # @property
    # def masses(self):
    #     if self._mass is None:
    #         return np.array([a["Mass"] for a in self._ats])
    #     else:
    #         return self._mass
    def _atomic_masses(self):
        """
        **LLM Docstring**

        The atomic masses converted to atomic units (electron masses) if they appear to be given in amu (heuristically, if the smallest mass is below 100).

        :return: the masses in atomic units
        :rtype: np.ndarray
        """
        m = self.masses
        if np.min(m) < 100:
            m = m*UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        return m
    @property
    def atomic_masses(self):
        """
        **LLM Docstring**

        The atomic masses in atomic units (electron masses), via `_atomic_masses`.

        :return: the masses in atomic units
        :rtype: np.ndarray
        """
        return self._atomic_masses()
    @property
    def bonds(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's bonds. The getter lazily guesses the bonds (via `get_guessed_bonds`) if none are set and `self.guess_bonds` is enabled.

        :param b: (setter only) the new bonds
        :type b: list[tuple] | None
        :return: (getter) the bonds, or `None` if unset and bond-guessing is disabled
        :rtype: list[tuple] | None
        """
        if self._bonds is None and self.guess_bonds:
            self._bonds = self.get_guessed_bonds()
        return self._bonds
    @bonds.setter
    def bonds(self, b):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's bonds. The getter lazily guesses the bonds (via `get_guessed_bonds`) if none are set and `self.guess_bonds` is enabled.

        :param b: (setter only) the new bonds
        :type b: list[tuple] | None
        :return: (getter) the bonds, or `None` if unset and bond-guessing is disabled
        :rtype: list[tuple] | None
        """
        self._bonds = b
    def break_bonds(self, bonds, use_rdkit=False, **rdopts):
        """
        **LLM Docstring**

        Build a copy of this molecule with the specified bonds removed, either by delegating to the attached RDKit molecule's `break_bonds` or by filtering the bond list directly.

        :param bonds: the bonds to remove, each an atom-index pair
        :type bonds: Iterable[tuple]
        :param use_rdkit: whether to perform the break via the RDKit molecule instead of filtering `self.bonds` directly
        :type use_rdkit: bool
        :param rdopts: extra options forwarded to the RDKit `break_bonds` call
        :type rdopts: dict
        :return: the new molecule with the given bonds removed
        :rtype: Molecule
        """
        if use_rdkit:
            return self.from_rdmol(self.rdmol.break_bonds(bonds, **rdopts))
        else:
            bond_sets = [{b[0], b[1]} for b in bonds]
            return self.modify(
                bonds=[b for b in self.bonds if
                       all(b[0] not in bs or b[1] not in bs for bs in bond_sets)]
            )
    @property
    def formula(self):
        """
        **LLM Docstring**

        The molecule's chemical formula, via `self.prop('chemical_formula')`.

        :return: the chemical formula
        :rtype: str
        """
        return self.prop('chemical_formula')
    @property
    def multiconfig(self):
        """
        **LLM Docstring**

        Whether this molecule holds multiple geometry configurations at once, delegating to `self.coords.multiconfig`.

        :return: whether multiple configurations are stored
        :rtype: bool
        """
        return self.coords.multiconfig
    @property
    def name(self):
        """
        **LLM Docstring**

        The molecule's name, falling back to `"Unnamed"` if none was set.

        :return: the molecule's name
        :rtype: str
        """
        if self._name is None:
            return "Unnamed"
        else:
            return self._name
    @property
    def source_file(self):
        """
        **LLM Docstring**

        Setter for the source file: if given a plain path string, infers the source `mode` from its extension and wraps both into a `{'file':..., 'mode':...}` dict; an already-structured dict is stored as-is.

        :param src: the new source file, as a path string or a `{'file', 'mode'}` dict
        :type src: str | dict | None
        :return: None
        :rtype: None
        """
        if self._src is not None:
            return self._src['file']
    @source_file.setter
    def source_file(self, src):
        """
        **LLM Docstring**

        Setter for the source file: if given a plain path string, infers the source `mode` from its extension and wraps both into a `{'file':..., 'mode':...}` dict; an already-structured dict is stored as-is.

        :param src: the new source file, as a path string or a `{'file', 'mode'}` dict
        :type src: str | dict | None
        :return: None
        :rtype: None
        """
        if isinstance(src, str):
            path, ext = os.path.splitext(src)
            ext = ext.lower()
            mode = ext.strip(".")

            src = {
                'file':src,
                'mode':mode
            }
        self._src = src

    @property
    def source_mode(self):
        """
        **LLM Docstring**

        The inferred/stored source-file mode (e.g. `'fchk'`, `'log'`), if a source file is set.

        :return: the source mode, or `None` if no source file is set or no mode was recorded
        :rtype: str | None
        """
        if self._src is not None:
            return self._src.get('mode')

    @property
    def shape(self):
        """
        **LLM Docstring**

        The shape of the molecule's Cartesian coordinates.

        :return: the coordinate array's shape
        :rtype: tuple[int]
        """
        return self.coords.shape
    def __len__(self):
        """
        **LLM Docstring**

        The number of geometry configurations held by this molecule: the leading dimension of its coordinates if `multiconfig`, otherwise `1`.

        :return: the number of configurations
        :rtype: int
        """
        if self.multiconfig:
            return self.coords.shape[0]
        else:
            return 1

    def copy(self):
        """
        **LLM Docstring**

        Make a copy of this molecule by calling `modify()` with no overrides.

        :return: the copied molecule
        :rtype: Molecule
        """
        return self.modify()

        # import copy
        # # mostly just use the default and don't be fancy
        # new = copy.copy(self)
        # # but we also need to do some stuff where we store objects that
        # # reference the molecule
        # new.normal_modes = new.normal_modes.copy()
        # new.normal_modes.set_molecule(new)
        # new.potential_surface = new.potential_surface.copy()
        # new.potential_surface.set_molecule(new)
        # new.dipole_surface = new.dipole_surface.copy()
        # new.dipole_surface.set_molecule(new)
        # # new._rdmol = new.rdmol.copy()
        # # new.rdmol.set_molecule(new)
        # return new

    def take_submolecule(self, pos):
        """
        **LLM Docstring**

        Build a new molecule containing only the atoms at the given positions, remapping bonds to the new (sub)indexing and dropping any bonds that reference atoms outside the subset.

        :param pos: the atom indices to keep, in the desired order for the submolecule
        :type pos: Iterable[int]
        :return: the constructed submolecule
        :rtype: Molecule
        """
        ats = self.atoms
        atoms = [ats[i] for i in pos]
        masses = self.masses[pos,]
        coords = self.coords[..., pos, :]
        pos = {p:i for i,p in enumerate(pos)}
        if self.bonds is not None:
            bonds = [
                (
                    [pos[b[0]], pos[b[1]]]
                        if len(b) == 2 else
                    [pos[b[0]], pos[b[1]], b[2]]
                ) for b in self.bonds
                if b[0] in pos and b[1] in pos
            ]
        else:
            bonds = None
        return type(self)(
            atoms,
            coords,
            masses=masses,
            bonds=bonds
        )

    def prop(self, name, *args, **kwargs):
        """
        **LLM Docstring**

        Compute a named derived molecular property by dispatching to the corresponding function on `MolecularProperties`.

        :param name: the property name, matching an attribute of `MolecularProperties`
        :type name: str
        :param args: positional arguments forwarded to the property function
        :type args: tuple
        :param kwargs: keyword arguments forwarded to the property function
        :type kwargs: dict
        :return: the computed property value
        :rtype: object
        :raises MolecularPropertyError: if `name` doesn't match a known property
        """
        from .Properties import MolecularProperties, MolecularPropertyError
        if hasattr(MolecularProperties, name):
            return getattr(MolecularProperties, name)(self, *args, **kwargs)
        else:
            raise MolecularPropertyError("{}.{}: property '{}' unknown".format(
                type(self).__name__,
                'prop',
                name
            ))

    bond_guessing_mode = 'rdkit'
    def get_guessed_bonds(self, mode=None, **opts):
        """
        **LLM Docstring**

        Guess the bonding arrangement for this molecule, either via RDKit (from the Cartesian geometry) or via `MolecularProperties.guessed_bonds`, depending on `mode`.

        :param mode: the bond-guessing strategy to use (`'rdkit'` or another mode understood by `MolecularProperties.guessed_bonds`); defaults to `self.bond_guessing_mode`
        :type mode: str | None
        :param opts: extra options forwarded to the underlying bond-guessing routine
        :type opts: dict
        :return: the guessed bonds
        :rtype: list[tuple]
        """
        if mode is None:
            mode = self.bond_guessing_mode
        if mode == 'rdkit':
            return RDMolecule.from_coords(
                self.atoms,
                self.coords * UnitsData.convert("BohrRadius", "Angstroms"),
                None,
                guess_bonds=True,
                formal_charges=self.formal_charges,
                **opts
            ).bonds
        else:
            return MolecularProperties.guessed_bonds(self, **opts)
        # self._bonds = self.prop("guessed_bonds", tol=1.05, guess_type=True)

    @property
    def edge_graph(self) -> EdgeGraph:
        """
        **LLM Docstring**

        The (cached) `EdgeGraph` representation of the molecule's bonding structure, built lazily via `MolecularProperties.edge_graph`.

        :return: the edge graph
        :rtype: EdgeGraph
        """
        if self._edge_graph is None:
            self._edge_graph = MolecularProperties.edge_graph(self)
        return self._edge_graph

    def find_path(self, atom1, atom2):
        """
        **LLM Docstring**

        Find a path between two atoms through the bonding graph.

        :param atom1: the starting atom index
        :type atom1: int
        :param atom2: the ending atom index
        :type atom2: int
        :return: the path between the two atoms
        :rtype: list[int]
        """
        return self.edge_graph.get_path(atom1, atom2)

    def find_substructure(self, pattern):
        """
        **LLM Docstring**

        Find matches of a substructure pattern within the molecule, delegating to the attached RDKit molecule.

        :param pattern: the substructure pattern to search for (e.g. a SMARTS string)
        :type pattern: str
        :return: the matching substructures
        :rtype: object
        """
        return self.rdmol.find_substructure(pattern)

    def apply_smarts(self, pattern):
        """
        **LLM Docstring**

        Apply a SMARTS reaction/transformation pattern to the molecule via RDKit, returning each resulting product as a new `Molecule`.

        :param pattern: the SMARTS pattern to apply
        :type pattern: str
        :return: the resulting molecules
        :rtype: list[Molecule]
        """
        new_mols = self.rdmol.apply_smarts(pattern)
        return [
            self.from_rdmol(m)
            for m in new_mols
        ]

    def neighborhood(self, loc, size=1):
        """
        **LLM Docstring**

        Find the atoms within a given graph-distance of a location in the bonding graph.

        :param loc: the atom index to center the neighborhood on
        :type loc: int
        :param size: the neighborhood radius (in bond-graph steps)
        :type size: int
        :return: the neighboring atom indices
        :rtype: tuple[int]
        """
        return tuple(l for l in self.edge_graph.neighbor_iterator(loc, num=size))

    def remove_hydrogens(self, positions=None, max=None, *, hydrogen_types=None):
        """
        **LLM Docstring**

        Build a copy of this molecule with hydrogen-type atoms (`H`/`D`/`T` by default) removed, either all of them or only those neighboring specified positions (optionally capped at `max` per position).

        :param positions: atom(s) whose neighboring hydrogens should be removed; if `None`, every hydrogen-type atom in the molecule is removed
        :type positions: int | Iterable[int] | None
        :param max: maximum number of hydrogens to remove per position in `positions`
        :type max: int | None
        :param hydrogen_types: the set of element symbols treated as hydrogen isotopes; defaults to `{'H', 'D', 'T'}`
        :type hydrogen_types: set[str] | None
        :return: the molecule with the selected hydrogens removed
        :rtype: Molecule
        """
        if hydrogen_types is None: hydrogen_types = {'H', 'D', 'T'}
        if positions is None:
            dropped = [i for i,a in enumerate(self.atoms) if a in hydrogen_types]
        else:
            dropped = set()
            if nput.is_int(positions): positions = [positions]
            for p in positions:
                subhs = [
                    i for i in self.neighborhood(p) if self.atoms[i] in hydrogen_types
                ]
                if max is not None:
                    subhs = subhs[:max]
                dropped.update(subhs)
            dropped = list(dropped)
        return self.take_submolecule(np.setdiff1d(np.arange(len(self.atoms)), dropped))

    def fragment_embedding(self, fragment_indices,
                           ref=None,
                           return_axes=False,
                           view_inds=(1, 2),
                           use_moments=False):
        """
        **LLM Docstring**

        Compute a local coordinate frame (origin, offset vector, and an up-vector or full axis set) anchored at a fragment of the molecule, used as the reference frame for attaching or orienting substituents; falls back to center-of-mass/principal-axis reference points (encoded as indices `-1`/`-2`/`-3`) when the fragment doesn't have enough atoms of its own to define a frame.

        :param fragment_indices: the atom index (or indices) defining the fragment to embed
        :type fragment_indices: int | Iterable[int]
        :param ref: reference atom(s) (outside the fragment) used to anchor the origin/frame; computed from the local neighborhood if not given
        :type ref: Iterable[int] | None
        :param return_axes: whether to return a full 3x3 axis frame instead of just an up-vector
        :type return_axes: bool
        :param view_inds: which two fragment-atom positions define the "view" direction used to build the axis frame
        :type view_inds: tuple[int, int]
        :param use_moments: whether to derive the up-vector from the fragment's moments of inertia rather than from its first three atom positions
        :type use_moments: bool
        :return: `(origin, offset, up_or_axes)` -- the reference origin point, the offset from origin to the fragment's first atom, and either an up-vector or (if `return_axes`) a full axis frame
        :rtype: tuple
        """
        if nput.is_int(fragment_indices):
            fragment_indices = [fragment_indices]
        if ref is None:
            ref = self.neighborhood(fragment_indices[0], size=2)
        if len(fragment_indices) < 3:
            if len(fragment_indices) + len(ref) < 3:
                extra = self.neighborhood(fragment_indices[0], size=2)
                ref = np.concatenate([ref, [r for r in extra if r not in ref]])
            ref = [r for r in ref if r not in fragment_indices]

            fragment_indices = np.concatenate([fragment_indices, ref])
            if len(fragment_indices) < 3: # these refer to COM and principle axis positions
                fragment_indices = np.concatenate([fragment_indices, [-1, -2, -3]])

        ref = np.asanyarray(ref)[0]
        r_coords = self.coords[ref,]
        if ref == -1: r_coords = self.center_of_mass
        if ref == -2: r_coords = self.center_of_mass + self.inertial_axes[:, 2]
        if ref == -3: r_coords = self.center_of_mass + self.inertial_axes[:, 0]
        if use_moments:
            f_coords = self.coords[fragment_indices,]
            _, axes = nput.moments_of_inertia(f_coords, np.asanyarray(self.masses)[fragment_indices,])
            up = axes[:, 2]
        else:
            # need origin, displacement vector, and up-vector
            fragment_indices = np.asanyarray(fragment_indices)[:3]
            f_coords = self.coords[fragment_indices,]
            com_f = ref == -1
            if np.any(com_f):
                f_coords[com_f] = self.center_of_mass
            paxc_f = ref == -2
            if np.any(paxc_f):
                f_coords[paxc_f] = self.center_of_mass + self.inertial_axes[:, 2]
            paxa_f = ref == -3
            if np.any(paxa_f):
                f_coords[paxa_f] = self.center_of_mass + self.inertial_axes[:, 0]

            up = nput.pts_normals(*f_coords, normalize=True)
            if return_axes:
                axes = nput.view_matrix(
                    up_vector=up,
                    view_vector=f_coords[view_inds[0]] - f_coords[view_inds[1]]
                )
            else:
                axes = np.eye(3)


        origin = r_coords
        offset = f_coords[0] - r_coords

        if return_axes:
            return origin, offset, axes
        else:
            return origin, offset, up

    def attach_functional_group(self,
                                target_fragment,
                                atoms,
                                new_coords,
                                bonds='recompute',
                                ref=None,
                                masses=None,
                                distance='auto',
                                angle=0,
                                dihedral=0,
                                embedding='auto',
                                bond_order=None,
                                use_absolue_posititions=False,
                                group_site=None
                                ) -> 'typing.Self':
        """
        **LLM Docstring**

        Build a copy of this molecule with a new group of atoms (`atoms`/`new_coords`) attached at `target_fragment`, positioning and orienting the new group using the fragment's local reference frame (bond distance/angle/dihedral, or an explicit embedding), and splicing the corresponding bonds into the result; supports designating a `group_site` atom within the new group as its attachment point, in which case the method recurses after re-deriving the embedding/bonds relative to that site.

        :param target_fragment: the atom(s) of this molecule the new group attaches to/replaces
        :type target_fragment: int | Iterable[int]
        :param atoms: the element symbols of the atoms in the new group
        :type atoms: Iterable[str]
        :param new_coords: the (local) coordinates of the new group's atoms
        :type new_coords: np.ndarray
        :param bonds: bonds within the new group; `'recompute'` to guess them fresh, `None` to reuse `self.bonds` remapped, or an explicit bond list
        :type bonds: str | list | None
        :param ref: reference atom(s) used to anchor the attachment frame; computed automatically if not given
        :type ref: Iterable[int] | None
        :param masses: masses for the new group's atoms; looked up from `atoms` if not given
        :type masses: np.ndarray | None
        :param distance: the bond distance to place the new group at; `'auto'` to look it up from `BondData`, or `None`/a number
        :type distance: str | float | None
        :param angle: rotation angle (about the up-vector) to apply to the new group
        :type angle: float
        :param dihedral: rotation angle (about the offset axis) to apply to the new group
        :type dihedral: float
        :param embedding: the reference orientation for the new group; `'auto'` to derive it from moments of inertia, or an explicit `(origin, axes)`/axes specification
        :type embedding: str | tuple | np.ndarray | None
        :param bond_order: the bond order connecting the new group to the target fragment; defaults to `1` (or inferred when `group_site` is used)
        :type bond_order: float | None
        :param use_absolue_posititions: whether `new_coords` should be used as absolute coordinates rather than being repositioned relative to the fragment frame
        :type use_absolue_posititions: bool
        :param group_site: index (within `atoms`/`new_coords`) of the atom that should serve as the attachment point; if given, the method recurses with the group re-anchored at this site and that atom excluded from the final group
        :type group_site: int | None
        :return: the molecule with the new group attached
        :rtype: Molecule
        """
        if group_site is not None:
            if dev.str_is(embedding, 'auto'):
                if dev.str_is(bonds, 'recompute') or bonds is not None:
                    if dev.str_is(bonds, 'recompute'): bonds = None
                    mol = Molecule(atoms, new_coords, bonds=bonds)
                    if bond_order is None:
                        bond_order = next(
                            (b[2] for b in mol.bonds if b[1] in (0, group_site) and b[0] in (0, group_site)),
                            None
                        )
                    embedding = mol.fragment_embedding(
                        0,
                        ref=[group_site]
                    )
                    embedding = nput.view_matrix(up_vector=embedding[-1])
                elif bonds == None:
                    _, embedding = nput.moments_of_inertia(new_coords, masses=masses)
                # new_coords = (new_coords - new_coords[(0,)]) @ embedding
                # origin, offset, up = self.fragment_embedding(target_fragment, ref=ref)
                # embedding = nput.view_matrix(up_vector=up)
                # new_coords = new_coords @ embedding
                # u = new_coords[0] - new_coords[group_site]
                # new_coords = new_coords @ nput.rotation_matrix(u, offset)
                embedding = new_coords[group_site], embedding
            rem = np.setdiff1d(np.arange(len(atoms)), [group_site])
            if bonds is not None and not dev.str_is(bonds, 'recompute'):
                remapping = {i: n for n, i in enumerate(rem)}
                bonds = [
                    [remapping[b[0]], remapping[b[1]], 1 if len(b) == 2 else b[2]]
                    for b in bonds
                    if b[0] in remapping and b[1] in remapping
                ]
            if masses is not None:
                masses = np.asanyarray(masses)[rem,]
            return self.attach_functional_group(
                target_fragment,
                [atoms[i] for i in rem],
                new_coords[rem,],
                bonds=bonds,
                ref=ref,
                masses=masses,
                distance=distance,
                angle=angle,
                dihedral=dihedral,
                embedding=embedding,
                bond_order=bond_order,
                use_absolue_posititions=use_absolue_posititions,
                group_site=None
            )

        if bond_order is None:
            bond_order = 1
        if dev.str_is(distance, 'auto'):
            if ref is None:
                ref = self.neighborhood(target_fragment[0], size=2)
            if ref[0] > -1:
                ref_type = self.atoms[ref[0]]
                distance = BondData[(ref_type, atoms[0], bond_order)] * UnitsData.convert("Angstroms", "BohrRadius")
            else:
                distance = None

        if masses is None:
            masses = np.array([AtomData[a, "Mass"] for a in atoms])

        if not use_absolue_posititions:
            origin, offset, up = self.fragment_embedding(target_fragment, ref=ref)
            if distance is not None:
                offset = nput.vec_normalize(offset) * distance
            if len(new_coords) == 1:
                new_coords = (origin + offset)[np.newaxis]
            else:
                if dev.str_is(embedding, 'auto'):
                    _, embedding = nput.moments_of_inertia(new_coords, masses=masses)
                shift = new_coords[0]
                new_coords = new_coords - shift[np.newaxis]
                if embedding is not None:
                    if len(embedding) == 2:
                        cent, embedding = embedding
                        cent = cent - shift
                    else:
                        cent = nput.center_of_mass(new_coords, masses=masses)
                    u = new_coords[0] - cent
                    inv = nput.view_matrix(
                        up_vector=embedding[:, 2],
                        view_vector=u
                    )
                    rot = nput.view_matrix(
                        up_vector=up,
                        view_vector=offset
                    )
                    new_coords = new_coords @ (inv @ rot.T)
                    if angle != 0:
                        new_coords = new_coords @ nput.rotation_matrix(up, angle)
                    if dihedral != 0:
                        # rotate about offset axis, assumed to be perp to up
                        new_coords = new_coords @ nput.rotation_matrix(offset, dihedral)
                new_coords = (origin + offset)[np.newaxis] + new_coords

        rem = np.setdiff1d(np.arange(len(self.atoms)), target_fragment)
        remapping = {i:n for n,i in enumerate(rem)}
        if dev.str_is(bonds, 'recompute'):
            total_bonds = None
        elif bonds is None:
            total_bonds = [
                [remapping[b[0]], remapping[b[1]], 1 if len(b) == 2 else b[2]]
                for b in self.bonds
                if b[0] in remapping and b[1] in remapping
            ]
        else:
            n = len(rem)
            total_bonds = [
                [remapping[b[0]], remapping[b[1]], 1 if len(b) == 2 else b[2]]
                for b in self.bonds
                if b[0] in remapping and b[1] in remapping
            ] + [
                [remapping[ref[0]], n, bond_order]
            ] + [
                [b[0]+n, b[1]+n, 1 if len(b) == 2 else b[2]]
                for b in bonds
            ]

        return self.modify(
            atoms=tuple(self.atoms[i] for i in rem) + tuple(atoms),
            coords=np.concatenate([
                self.coords[rem],
                new_coords
            ]),
            bonds=total_bonds,
            masses=np.concatenate([[self.masses[i] for i in rem], masses])
        )

    def find_heavy_atom_backbone(self, root=None):
        """
        **LLM Docstring**

        Find the longest chain of atoms in the bonding graph (the heavy-atom backbone), via `edge_graph.find_longest_chain`.

        :param root: an atom index to force as one end of the chain
        :type root: int | None
        :return: the backbone atom indices, in chain order
        :rtype: list[int]
        """
        return self.edge_graph.find_longest_chain(root=root)

    def find_backbone_segments(self, root=None, initial_backbone=None):
        """
        **LLM Docstring**

        Split the bonding graph into backbone-connected segments, via `edge_graph.segment_by_chains`.

        :param root: an atom index to anchor the segmentation at
        :type root: int | None
        :param initial_backbone: an initial backbone chain to build the segmentation around
        :type initial_backbone: Iterable[int] | None
        :return: the resulting chain segments
        :rtype: list
        """
        return self.edge_graph.segment_by_chains(root=root, backbone=initial_backbone)

    def get_backbone_zmatrix(self, root=None,
                             segments=None,
                             return_remainder=False,
                             return_segments=False,
                             required_coordinates=None,
                             isolated_coordinates=None,
                             root_coordinates=None,
                             initial_backbone=None,
                             validate=True
                             ):
        """
        **LLM Docstring**

        Build a Z-matrix for a (typically single-fragment) molecule by first segmenting its bonding graph into backbone chains (via `find_backbone_segments`, validating there are no duplicate atoms across segments), constructing the base Z-matrix graph from the bonds and segments, filling in any bonds missing from the initial graph, and (if requested) enforcing required/isolated/root coordinate constraints.

        :param root: an atom index to anchor the backbone segmentation at
        :type root: int | None
        :param segments: precomputed backbone segments to use instead of calling `find_backbone_segments`
        :type segments: list | None
        :param return_remainder: whether to also return the bonds that had to be added beyond the base backbone graph
        :type return_remainder: bool
        :param return_segments: whether to also return the backbone segments used
        :type return_segments: bool
        :param required_coordinates: internal coordinates that must appear in the resulting Z-matrix
        :type required_coordinates: Iterable | None
        :param isolated_coordinates: coordinates that must be built in isolation from others
        :type isolated_coordinates: Iterable | None
        :param root_coordinates: coordinates to anchor at the root of the Z-matrix
        :type root_coordinates: Iterable | None
        :param initial_backbone: an initial backbone chain to seed the segmentation with
        :type initial_backbone: Iterable[int] | None
        :param validate: whether to validate the Z-matrix construction at each step (duplicate atoms, valid additions)
        :type validate: bool
        :return: the Z-matrix, or a tuple additionally including the segments and/or remainder bonds depending on the `return_*` flags
        :rtype: np.ndarray | tuple
        :raises ValueError: if `validate` is set and duplicate atoms are found across backbone segments
        """
        if segments is None:
            segments = self.find_backbone_segments(root=root, initial_backbone=initial_backbone)
            if validate:
                flat_frags = list(itut.flatten(segments))
                frag_counts = itut.counts(flat_frags)
                bad_frags = {k: v for k, v in frag_counts.items() if v > 1}
                if len(bad_frags) > 0:
                    raise ValueError(f"diplicate atoms {list(bad_frags.keys())} encountered in {segments}")

        bond_list = [b[:2] for b in self.bonds]
        base_graph = coordops.bond_graph_zmatrix(
            bond_list,
            segments,
            validate_additions=validate,
            required_coordinates=required_coordinates,
            isolated_coordinates=isolated_coordinates,
            root_coordinates=root_coordinates,
            enforce_requirements=False
        )
        zmat, new_bonds = coordops.add_missing_zmatrix_bonds(
            base_graph,
            bond_list,
            validate_additions=validate
        )
        if (
                required_coordinates is not None
                or isolated_coordinates is not None
                or root_coordinates is not None
        ):
            zmat = coordops.enforce_required_zmatrix_coordinates(zmat,
                                                                 required_coordinates,
                                                                 isolated_coordinates=isolated_coordinates,
                                                                 root_coordinates=root_coordinates,
                                                                 validate=validate)

        if return_segments or return_remainder:
            res = (zmat,)
            if return_segments:
                res = res + (segments,)
            if return_remainder:
                res = res + (new_bonds,)

            return res
        else:
            return zmat

    def get_canonical_zmatrix(self, ordering=None, validate=True):
        """
        **LLM Docstring**

        Build a canonical Z-matrix ordering for the molecule from a canonical fragmentation of the bonding graph.

        :param ordering: the atom ordering to use as the basis for canonicalization; defaults to the natural `0..N` ordering
        :type ordering: np.ndarray | None
        :param validate: whether to validate each Z-matrix addition
        :type validate: bool
        :return: the canonical Z-matrix
        :rtype: np.ndarray
        """
        if ordering is None: ordering = np.arange(len(self.atoms))
        frags = self.edge_graph.get_canonical_fragments(ordering)
        return coordops.canonical_fragment_zmatrix(frags, validate_additions=validate)

    @classmethod
    def _filter_coordinates_by_fragments(cls, inds, frags, required_coordinates):
        """
        **LLM Docstring**

        Split a list of required internal-coordinate specs into those fully contained within a single fragment (reindexed to that fragment's local atom numbering) versus those that span multiple fragments and must be merged in afterward.

        :param inds: the atom-index lists defining each fragment
        :type inds: list[Iterable[int]]
        :param frags: the corresponding submolecules for each fragment (unused directly in the body but kept for interface/length consistency)
        :type frags: list[Molecule]
        :param required_coordinates: the coordinate specs to classify, each a tuple of atom indices in the full-molecule numbering
        :type required_coordinates: Iterable[tuple] | None
        :return: `(merge_coordinates, fragment_requireds)` -- the coordinates that couldn't be assigned to a single fragment, and a per-fragment list of the (locally reindexed) coordinates assigned to each fragment
        :rtype: tuple[list, list]
        """
        merge_coordinates = []
        if required_coordinates is not None:
            fragment_requireds = [[] for _ in range(len(frags))]
            for c in required_coordinates:
                for i, f in enumerate(inds):
                    ff = list(f)
                    sub = []
                    for j in c:
                        try:
                            x = ff.index(j)
                        except ValueError:
                            break
                        else:
                            sub.append(x)
                    if len(sub) == len(c):
                        fragment_requireds[i].append(tuple(sub))
                        break
                else:
                    merge_coordinates.append(c)
        else:
            fragment_requireds = [None] * len(inds)
        return merge_coordinates, fragment_requireds
    def get_bond_zmatrix(self,
                         fragments=None,
                         segments=None,
                         root=None,
                         required_coordinates=None,
                         isolated_coordinates=None,
                         root_coordinates=None,
                         attachment_points=None,
                         check_attachment_points=True,
                         validate=True,
                         for_fragment=None,
                         fragment_ordering=None,
                         connect_fragments=True,
                         initial_backbone=None
                         ):
        """
        **LLM Docstring**

        Build a full Z-matrix for the molecule from its bonding graph, handling the single-fragment case via `get_backbone_zmatrix` directly and the multi-fragment case by building a per-fragment Z-matrix for each fragment (optionally reordering/rooting/filtering required coordinates per fragment) and then splicing them together into one connected Z-matrix (via `coordops.complex_zmatrix`) unless `connect_fragments` is `False`; can also be restricted to build the Z-matrix for just one fragment (`for_fragment`), in which case it recurses on the corresponding submolecule and reindexes the result back to the full atom numbering.

        :param fragments: explicit fragment atom-index groups to use instead of `self.fragment_indices`
        :type fragments: list[Iterable[int]] | None
        :param segments: precomputed backbone segments for the single-fragment case
        :type segments: list | None
        :param root: root atom(s) to anchor the Z-matrix construction at (per fragment, in the multi-fragment case)
        :type root: int | Iterable | None
        :param required_coordinates: internal coordinates that must appear in the resulting Z-matrix
        :type required_coordinates: Iterable | None
        :param isolated_coordinates: coordinates that must be built in isolation from others
        :type isolated_coordinates: Iterable | None
        :param root_coordinates: coordinates to anchor at the root of the Z-matrix
        :type root_coordinates: Iterable | None
        :param attachment_points: explicit inter-fragment attachment points to use when connecting fragments
        :type attachment_points: dict | Iterable | None
        :param check_attachment_points: whether to validate the attachment points used when connecting fragments
        :type check_attachment_points: bool
        :param validate: whether to validate each Z-matrix addition
        :type validate: bool
        :param for_fragment: restrict the Z-matrix construction to just this fragment (an index into `self.fragment_indices`, or an explicit list of atom indices), returning the result reindexed to the full molecule
        :type for_fragment: int | Iterable[int] | None
        :param fragment_ordering: explicit ordering to apply to the fragments before connecting them
        :type fragment_ordering: Iterable[int] | None
        :param connect_fragments: whether to splice the per-fragment Z-matrices into one connected Z-matrix, or return them as a list of separate Z-matrices
        :type connect_fragments: bool
        :param initial_backbone: an initial backbone chain to seed the segmentation of each fragment with
        :type initial_backbone: Iterable[int] | None
        :return: the (connected) Z-matrix, or a list of per-fragment Z-matrices if `connect_fragments` is `False`
        :rtype: np.ndarray | list[np.ndarray]
        """
        if for_fragment is not None:
            if nput.is_int(for_fragment):
                for_fragment = self.fragment_indices[for_fragment]
            if attachment_points is not None:
                frag_map = dict(zip(for_fragment, np.arange(len(for_fragment))))
                if hasattr(attachment_points, 'items'):
                    attachment_points = {
                        frag_map[i]: (frag_map[k] if nput.is_numeric(k) else tuple(frag_map[kk] for kk in k))
                        for i,k in attachment_points.items()
                    }
                else:
                    attachment_points = [
                        frag_map[i] for i in attachment_points
                    ]

            if initial_backbone is not None:
                ff = list(for_fragment)
                initial_backbone = [ff.index(i) for i in initial_backbone]
            base_ints = self.take_submolecule(for_fragment).get_bond_zmatrix(
                fragments=fragments, segments=segments, root=root,
                attachment_points=attachment_points,
                check_attachment_points=check_attachment_points,
                fragment_ordering=fragment_ordering,
                required_coordinates=required_coordinates,
                isolated_coordinates=isolated_coordinates,
                root_coordinates=root_coordinates,
                initial_backbone=initial_backbone,
                validate=validate
            )
            zm = coordops.reindex_zmatrix(base_ints, for_fragment)
            return np.asarray(zm)
        else:
            no_frag = fragments is None
            if no_frag:
                fragments = self.fragment_indices

            if len(fragments) == 1:
                if segments is not None and len(segments) == 1:
                    segments = segments[0]
                zm = self.get_backbone_zmatrix(
                    root=root, segments=segments,
                    required_coordinates=required_coordinates,
                    isolated_coordinates=isolated_coordinates,
                    root_coordinates=root_coordinates,
                    validate=validate,
                    initial_backbone=initial_backbone
                )
                zm = np.asarray(zm)
                if connect_fragments:
                    return zm
                else:
                    return [zm]
            else:
                inds = fragments
                if no_frag:
                    if fragment_ordering is None:
                        fragment_ordering = np.argsort([-len(x) for x in inds])
                    inds = [inds[i] for i in fragment_ordering]
                if root is not None:
                    if nput.is_numeric(root):
                        inds = list(sorted(inds, key=lambda x:root not in x))
                    else:
                        inds = list(
                            sorted(inds,
                                   key=lambda x:sum(i if r is not None and r in x else len(inds) for i,r in enumerate(root))
                                   )
                        )

                sort_attch = isinstance(attachment_points, dict)
                if sort_attch:
                    check_attachment_points = False
                    inds, attachment_points = coordops.sort_complex_attachment_points(
                        inds,
                        attachment_points
                    )

                frags = [self.take_submolecule(ix) for ix in inds]
                if root is None and sort_attch:
                    root = [ix[0] for ix in inds]
                if root is None:
                    root = [root]


                if initial_backbone is not None:
                    initial_backbones = []
                    for frag in inds:
                        ff = list(frag)
                        sub = []
                        for i in initial_backbone:
                            try:
                                x = ff.index(i)
                            except ValueError:
                                continue
                            else:
                                sub.append(x)
                        if len(sub) > 0:
                            initial_backbones.append(sub)
                        else:
                            initial_backbones.append(None)
                else:
                    initial_backbones = [None] * len(inds)

                root = list(root) + [None] * (len(inds) - len(root))
                merge_coordinates, fragment_requireds = self._filter_coordinates_by_fragments(inds, frags, required_coordinates)
                merge_isolated, fragment_isolated = self._filter_coordinates_by_fragments(inds, frags, isolated_coordinates)
                merge_root, fragment_root = self._filter_coordinates_by_fragments(inds, frags, root_coordinates)
                if len(merge_coordinates) == 0:
                    merge_coordinates = None
                if len(merge_isolated) == 0:
                    merge_isolated = None
                if len(merge_root) == 0:
                    merge_root = None
                zmats = [
                    f.get_backbone_zmatrix(root=r, initial_backbone=bb,
                                           required_coordinates=rq,
                                           isolated_coordinates=iso,
                                           root_coordinates=rot
                                           )
                    for r,f,bb,rq,iso,rot in zip(root, frags, initial_backbones,
                                                 fragment_requireds, fragment_isolated, fragment_root)
                ]

                if connect_fragments:

                    # inds = [inds[i] for i in ordering]
                    # zmats = [zmats[i] for i in ordering]

                    dm = nput.distance_matrix(self.coords)
                    h_pos = [i for i,a in enumerate(self.atoms) if a in {'H', 'D'}]
                    dm[:, h_pos] = 1e8
                    dm[h_pos, :] = 1e8

                    zm = coordops.complex_zmatrix(
                        [b[:2] for b in self.bonds],
                        inds,
                        zmats,
                        distance_matrix=dm,
                        attachment_points=attachment_points,
                        check_attachment_points=check_attachment_points,
                        required_coordinates=required_coordinates,
                        isolated_coordinates=isolated_coordinates,
                        root_coordinates=root_coordinates,
                        validate_additions=validate
                    )
                    return np.asarray(zm)
                else:
                    shift_mats = []
                    offsets = 0
                    for zmat in zmats:
                        shift_mats.append([
                            [(z + offsets) if z >= 0 else z for z in zm]
                            for zm in zmat
                        ])
                        offsets += len(zmat)
                    return [np.asarray(zm) for zm in shift_mats]

    @property
    def fragment_indices(self):
        """
        **LLM Docstring**

        The (cached) grouping of atom indices into connected molecular fragments, computed lazily via `MolecularProperties.fragment_indices`.

        :return: the list of per-fragment atom-index groups
        :rtype: list[np.ndarray]
        """
        if self._fragment_indices is None:
            self._fragment_indices = MolecularProperties.fragment_indices(self)
        return self._fragment_indices

    @property
    def fragments(self):
        """
        **LLM Docstring**

        The molecule split into its connected fragments (as separate `Molecule` objects), via `MolecularProperties.fragments`.

        :return: the list of fragment molecules
        :rtype: list[Molecule]
        """
        return MolecularProperties.fragments(self)
    #region Coordinate Embeddings

    @property
    def mass_weighted_coords(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        return MolecularProperties.mass_weighted_coords(self)
        # return self.prop('mass_weighted_coords')
    @property
    def center_of_mass(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        return MolecularProperties.center_of_mass(self)
    @property
    def inertia_tensor(self):
        """
        :return:
        :rtype: (np.ndarray, np.ndarray)
        """
        return MolecularProperties.inertia_tensor(self)
    @property
    def inertial_eigensystem(self):
        """
        :return:
        :rtype: (np.ndarray, np.ndarray)
        """
        return MolecularProperties.moments_of_inertia(self)
    @property
    def moments_of_inertia(self):
        """
        :return:
        :rtype: np.ndarray
        """
        return MolecularProperties.moments_of_inertia(self)[0]
    @property
    def inertial_axes(self):
        """
        :return:
        :rtype: np.ndarray
        """
        return MolecularProperties.moments_of_inertia(self)[1]

    @property
    def translation_rotation_modes(self):
        """
        :return:
        :rtype: np.ndarray
        """
        return MolecularProperties.translation_rotation_eigenvectors(self)

    def get_translation_rotation_projector(self, mass_weighted=False):
        """
        **LLM Docstring**

        Build the projector matrix that removes translational and rotational displacement components from a Cartesian displacement.

        :param mass_weighted: whether the projector should act on mass-weighted coordinates
        :type mass_weighted: bool
        :return: the translation/rotation projector
        :rtype: np.ndarray
        """
        return nput.translation_rotation_projector(
            self.coords,
            self.atomic_masses,
            mass_weighted=mass_weighted,
            return_modes=False
        )
        # L_tr = self.translation_rotation_modes[1]
        # A = np.eye(L_tr.shape[0]) - (L_tr @ L_tr.T)
        # if not mass_weighted:
        #     M = np.diag(np.repeat(1 / np.sqrt(self.masses), 3))
        #     A = M @ A @ M
        # return A

    def get_translation_rotation_invariant_transformation(self,
                                           mass_weighted=False,
                                           strip_embedding=True
                                           ):
        """
        **LLM Docstring**

        Build the transformation (and its inverse) that projects out the translational and rotational degrees of freedom from this molecule's Cartesian coordinates.

        :param mass_weighted: whether the transformation should act on mass-weighted coordinates
        :type mass_weighted: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the translation/rotation-invariant transformation and its inverse
        :rtype: tuple
        """
        return nput.translation_rotation_invariant_transformation(
            self.coords,
            self.atomic_masses,
            mass_weighted=mass_weighted,
            strip_embedding=strip_embedding
        )


    # @property
    # def zmatrix(self):
    #     """
    #     :return:
    #     :rtype:
    #     """
    #     return self._zmat
    # @zmatrix.setter
    # def zmatrix(self, zmatrix):
    #     """
    #     :return:
    #     :rtype:
    #     """
    #     #TODO: add some validation
    #     zmatrix = np.asanyarray(zmatrix).astype(int)
    #     if zmatrix.shape[1] != 4:
    #         raise ValueError("can't understand Z-matrix {}".format(zmatrix))
    #     self._zmat = zmatrix

    #region Evaluation

    def _load_energy(self):
        """
        **LLM Docstring**

        Load the molecule's total energy from its source file, currently only supported for Gaussian FChk files.

        :return: the total energy, or `None` if there is no source file or it isn't a supported format
        :rtype: float | None
        """
        if self._src is None:
            return None
        elif (
                dev.is_dict_like(self._src)
                and self._src['mode'] == '.fchk'
        ) or (
                os.path.splitext(self.source_file)[1] == '.fchk'
        ):
            from McUtils.GaussianInterface import GaussianFChkReader
            with GaussianFChkReader(self.source_file) as gr:
                parse = gr.parse(
                    ['Total Energy']
                )
            return parse['Total Energy']
        else:
            return None

    @property
    def energy(self):
        """
        **LLM Docstring**

        The molecule's total energy: loaded from its source file if possible, else computed via `calculate_energy` if an energy evaluator is available, and cached thereafter.

        :return: the total energy, or `None` if neither a source-file energy nor an evaluator is available
        :rtype: float | None
        """
        if self._energy is None:
            self._energy = self._load_energy()
            if self._energy is None and self.energy_evaluator is not None:
                self._energy = self.calculate_energy()
        return self._energy

    default_energy_evalutor = 'rdkit'
    def get_energy_evaluator(self, evaluator=None, **opts):
        """
        **LLM Docstring**

        Resolve (and, if needed, instantiate from this molecule) an energy-evaluator object, defaulting to `self.energy_evaluator` (then `self.default_energy_evalutor`) if none is given explicitly.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to resolve
        :type evaluator: object | None
        :param opts: extra options forwarded to the evaluator's `from_mol` constructor, if applicable
        :type opts: dict
        :return: the resolved energy-evaluator instance
        :rtype: object
        :raises ValueError: if the evaluator type can't be resolved
        """
        if evaluator is None:
            evaluator = self.energy_evaluator
        if evaluator is None:
            evaluator = self.default_energy_evalutor
        eval_type, new_opts = EnergyEvaluator.resolve_evaluator(evaluator)
        if eval_type is None:
            raise ValueError(f"can't resolve energy evaluator type for {evaluator}")
        if hasattr(eval_type, 'from_mol') and isinstance(eval_type, type):
            return eval_type.from_mol(self, **dict(opts, **new_opts))
        else:
            return eval_type

    def get_energy_function(self, evaluator=None, *, order=None, **opts):
        """
        **LLM Docstring**

        Build a standalone callable that evaluates this molecule's energy (and derivatives) at arbitrary coordinates using the resolved energy evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the derivative order the returned function should evaluate to by default
        :type order: int | None
        :param opts: extra options forwarded to `get_energy_evaluator`
        :type opts: dict
        :return: a function `(coords, order=order) -> energy_or_expansion`
        :rtype: callable
        """
        evaluator = self.get_energy_evaluator(evaluator, **opts)
        def evaluate_energy(coords, order=order):
            """
            **LLM Docstring**

            Evaluate the energy (and, if `order` is given, its derivatives) at `coords` using the evaluator from the enclosing scope, defaulting to this molecule's own coordinates if none are given.

            :param coords: the coordinates to evaluate at, in Bohr; defaults to `self.coords`
            :type coords: np.ndarray | None
            :param order: the highest derivative order to evaluate; if `None`, only the energy itself is returned
            :type order: int | None
            :return: the energy, or the full energy/derivative expansion
            :rtype: float | list
            """
            smol = order is None
            if smol: order = 0
            if coords is None:
                coords = self.coords
            expansion = evaluator.evaluate(
                np.asanyarray(coords) * UnitsData.convert("BohrRadius", evaluator.distance_units),
                order=order
            )
            if smol: expansion = expansion[0]
            return expansion
        return evaluate_energy

    def calculate_energy(self, coords=None, *, evaluator=None, order=None, **opts):
        """
        **LLM Docstring**

        Compute the energy (and, optionally, its derivatives) at the given (or this molecule's own) coordinates using the resolved energy evaluator.

        :param coords: the coordinates to evaluate at, in Bohr; defaults to `self.coords`
        :type coords: np.ndarray | None
        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the highest derivative order to compute; if `None`, only the energy itself is returned
        :type order: int | None
        :param opts: extra options forwarded to `get_energy_evaluator`
        :type opts: dict
        :return: the energy, or the full energy/derivative expansion
        :rtype: float | list
        """
        evaluator = self.get_energy_evaluator(evaluator, **opts)
        smol = order is None
        if smol: order = 0
        if coords is None:
            coords = self.coords
        expansion = evaluator.evaluate(
            np.asanyarray(coords) * UnitsData.convert("BohrRadius", evaluator.distance_units),
            order=order
        )
        if smol: expansion = expansion[0]
        return expansion
    def partial_force_field(self,
                            coords=None, modes=None, *,
                            evaluator=None,
                            order=4,
                            mesh_spacing=1,
                            analytic_derivative_order=None,
                            **opts):
        """
        **LLM Docstring**

        Compute a partial (mode-selected) force-field expansion of the potential, evaluated in the given normal modes via the resolved energy evaluator's `partial_force_field` method.

        :param coords: the coordinates to evaluate at; defaults to `self.coords`
        :type coords: np.ndarray | None
        :param modes: the normal modes to expand in; computed via `get_normal_modes` if not given
        :type modes: object | None
        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the highest derivative order to compute
        :type order: int
        :param mesh_spacing: finite-difference step size to use for the underlying force-field evaluation
        :type mesh_spacing: float
        :param analytic_derivative_order: order up to which derivatives should be computed analytically rather than numerically
        :type analytic_derivative_order: int | None
        :param opts: extra options forwarded to `get_energy_evaluator`
        :type opts: dict
        :return: the partial force-field expansion
        :rtype: object
        """
        if modes is None:
            modes = self.get_normal_modes()
        evaluator = self.get_energy_evaluator(evaluator, **opts)
        if coords is None:
            coords = self.coords
        expansion = evaluator.partial_force_field(
            np.asanyarray(coords) * UnitsData.convert("BohrRadius", evaluator.distance_units),
            order=order,
            modes=modes,
            mesh_spacing=mesh_spacing,
            analytic_derivative_order=analytic_derivative_order
        )
        return expansion

    def optimize(self,
                 evaluator=None,
                 *,
                 method=None,
                 tol=None,
                 max_iterations=None,
                 logger=None,
                 reembed=True,
                 **opts):
        """
        **LLM Docstring**

        Geometry-optimize the molecule using the resolved energy evaluator, splitting `opts` into optimizer-specific and evaluator-specific options, re-embedding the optimized (and any trajectory) coordinates back onto the original geometry via an Eckart embedding, and returning a modified copy of the molecule at the optimized geometry.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param method: the optimization method/algorithm to use
        :type method: str | None
        :param tol: convergence tolerance for the optimizer
        :type tol: float | None
        :param max_iterations: maximum number of optimizer iterations
        :type max_iterations: int | None
        :param logger: logger to report optimization progress/warnings to
        :type logger: Logger | str | None
        :param reembed: whether to Eckart-reembed the optimized geometry (and trajectory) onto the original one
        :type reembed: bool
        :param opts: extra options, split between optimizer options and evaluator-construction options
        :type opts: dict
        :return: the optimized molecule, or `(optimized_molecule, trajectory)` if the optimizer returned a trajectory
        :rtype: Molecule | tuple[Molecule, list]
        """
        opts = dev.OptionsSet(opts)
        base_opts = EnergyEvaluator.get_optimizer_options() + ('force_field_type',)
        optimizer_opts = opts.filter(None, props=base_opts)
        eval_opts = opts.exclude(None, props=base_opts)
        evaluator = self.get_energy_evaluator(evaluator, **eval_opts)
        conv = UnitsData.convert("BohrRadius", evaluator.distance_units)
        opt_params = dict(
            method=method,
            tol=tol,
            max_iterations=max_iterations,
            logger=logger
        )
        opt, opt_coords, settings = evaluator.optimize(
            self.coords * conv,
            **{k:v for k,v in opt_params.items() if v is not None},
            **optimizer_opts
        )
        if not opt:
            if logger is not None:
                logger = Logger.lookup(logger)
                logger.log_print('WARNING: failed to optimize molecule')
        settings['optimized'] = opt
        no_traj = isinstance(opt_coords, np.ndarray)
        if no_traj: opt_coords = (opt_coords, [])
        opt_coords, traj = opt_coords
        opt_coords = opt_coords / conv
        if reembed:
            opt_coords = nput.eckart_embedding(
                self.coords, opt_coords, masses=self.masses
            ).coordinates
        _ = []
        for t in traj:
            t = nput.eckart_embedding(
                self.coords, t, masses=self.masses
            ).coordinates
            _.append(t)
        traj = _
        if no_traj:
            return self.modify(coords=opt_coords, meta=settings)
        else:
            return self.modify(coords=opt_coords, meta=settings), traj

    def relaxed_scan(self,
                     scan_values,
                     scan_coordinates,
                     evaluator=None,
                     *,
                     method=None,
                     tol=None,
                     max_iterations=None,
                     logger=None,
                     reembed=False,
                     **opts):
        """
        **LLM Docstring**

        Perform a relaxed (constrained) potential-energy scan over the given coordinates/values using the resolved energy evaluator, optionally Eckart-reembedding the resulting trajectory onto the original geometry.

        :param scan_values: the coordinate value(s) to scan over
        :type scan_values: object
        :param scan_coordinates: the coordinate(s) to hold fixed/scan along
        :type scan_coordinates: object
        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param method: the optimization method/algorithm to use for the constrained relaxations
        :type method: str | None
        :param tol: convergence tolerance for the optimizer
        :type tol: float | None
        :param max_iterations: maximum number of optimizer iterations
        :type max_iterations: int | None
        :param logger: logger to report progress/warnings to
        :type logger: Logger | str | None
        :param reembed: whether to Eckart-reembed the resulting trajectory onto the original geometry
        :type reembed: bool
        :param opts: extra options, split between optimizer/scan options and evaluator-construction options
        :type opts: dict
        :return: `(opts, traj, meta)` -- the resolved options, the resulting geometry trajectory, and scan metadata
        :rtype: tuple
        """
        #TODO: minimized duplication
        opts = dev.OptionsSet(opts)
        base_opts = EnergyEvaluator.get_relaxed_scan_options() + ('force_field_type',)
        optimizer_opts = opts.filter(None, props=base_opts)
        eval_opts = opts.exclude(None, props=base_opts)
        evaluator = self.get_energy_evaluator(evaluator, **eval_opts)
        conv = UnitsData.convert("BohrRadius", evaluator.distance_units)
        opt_params = dict(
            method=method,
            tol=tol,
            max_iterations=max_iterations,
            logger=logger
        )
        # if self.internals is None:
        #     if nput.is_numeric(scan_values[0]):
        #         s, e, n = scan_values
        #         scan_values = (conv * s, conv * e, n)
        #     else:
        #         scan_values = [
        #             (conv * s, conv * e, n)
        #             for s,e,n in scan_values
        #         ]
        opts, traj, meta = evaluator.relaxed_scan(
            self.coords * conv,
            scan_values,
            scan_coordinates,
            **{k: v for k, v in opt_params.items() if v is not None},
            **optimizer_opts
        )
        traj = traj / conv
        if reembed:
            traj = nput.eckart_embedding(
                self.coords, traj, masses=self.masses
            ).coordinates
        return opts, traj, meta

    def get_dipole_evaluator(self, evaluator=None, **opts):
        """
        **LLM Docstring**

        Resolve (and, if needed, instantiate from this molecule) a dipole-evaluator object, defaulting to `self.dipole_evaluator` if none is given explicitly.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to resolve
        :type evaluator: object | None
        :param opts: extra options forwarded to the evaluator's `from_mol` constructor, if applicable
        :type opts: dict
        :return: the resolved dipole-evaluator instance
        :rtype: object
        :raises ValueError: if the evaluator type can't be resolved
        """
        if evaluator is None:
            evaluator = self.dipole_evaluator
        eval_type, new_opts = DipoleEvaluator.resolve_evaluator(evaluator)
        if eval_type is None:
            raise ValueError(f"can't resolve dipole evaluator type for {evaluator}")
        if hasattr(eval_type, 'from_mol') and isinstance(eval_type, type):
            return eval_type.from_mol(self, **dict(opts, **new_opts))
        else:
            return eval_type

    def get_dipole_function(self, evaluator=None, *, order=None, **opts):
        """
        **LLM Docstring**

        Build a standalone callable that evaluates this molecule's dipole (and derivatives) at arbitrary coordinates using the resolved dipole evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the derivative order the returned function should evaluate to by default
        :type order: int | None
        :param opts: extra options forwarded to `get_dipole_evaluator`
        :type opts: dict
        :return: a function `(coords, order=order) -> dipole_or_expansion`
        :rtype: callable
        """
        evaluator = self.get_dipole_evaluator(evaluator, **opts)
        def evaluate_dipole(coords, order=order):
            """
            **LLM Docstring**

            Evaluate the dipole (and, if `order` is given, its derivatives) at `coords` using the evaluator from the enclosing scope, defaulting to this molecule's own coordinates if none are given.

            :param coords: the coordinates to evaluate at, in Bohr; defaults to `self.coords`
            :type coords: np.ndarray | None
            :param order: the highest derivative order to evaluate; if `None`, only the dipole itself is returned
            :type order: int | None
            :return: the dipole, or the full dipole/derivative expansion
            :rtype: np.ndarray | list
            """
            smol = order is None
            if smol: order = 0
            if coords is None:
                coords = self.coords
            expansion = evaluator.evaluate(
                np.asanyarray(coords)
                    * UnitsData.convert("BohrRadius", evaluator.distance_units),
                order=order
            )
            if smol: expansion = expansion[0]
            return expansion
        return evaluate_dipole

    def calculate_dipole(self, evaluator=None, order=None, **opts):
        """
        **LLM Docstring**

        Compute the dipole moment (and, optionally, its derivatives) at this molecule's current coordinates using the resolved dipole evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the highest derivative order to compute; if `None`, only the dipole itself is returned
        :type order: int | None
        :param opts: extra options forwarded to `get_dipole_evaluator`
        :type opts: dict
        :return: the dipole, or the full dipole/derivative expansion
        :rtype: np.ndarray | list
        """
        evaluator = self.get_dipole_evaluator(evaluator, **opts)
        smol = order is None
        if smol: order = 0
        expansion = evaluator.evaluate(
            self.coords * UnitsData.convert("BohrRadius", evaluator.distance_units),
            order=order
        )
        if smol: expansion = expansion[0]
        return expansion

    def get_polarizability_evaluator(self, evaluator=None, **opts):
        """
        **LLM Docstring**

        Resolve (and, if needed, instantiate from this molecule) a dipole-polarizability-evaluator object, defaulting to `self.polarizability_evaluator` if none is given explicitly.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to resolve
        :type evaluator: object | None
        :param opts: extra options forwarded to the evaluator's `from_mol` constructor, if applicable
        :type opts: dict
        :return: the resolved polarizability-evaluator instance
        :rtype: object
        :raises ValueError: if the evaluator type can't be resolved
        """
        if evaluator is None:
            evaluator = self.polarizability_evaluator
        eval_type, new_opts = DipolePolarizabilityEvaluator.resolve_evaluator(evaluator)
        if eval_type is None:
            raise ValueError(f"can't resolve polarizability evaluator type for {evaluator}")
        if hasattr(eval_type, 'from_mol') and isinstance(eval_type, type):
            return eval_type.from_mol(self, **dict(opts, **new_opts))
        else:
            return eval_type
    def get_dipole_polarizability_function(self, evaluator=None, *, order=None, **opts):
        """
        **LLM Docstring**

        Build a standalone callable that evaluates this molecule's dipole polarizability (and derivatives) at arbitrary coordinates using the resolved polarizability evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the derivative order the returned function should evaluate to by default
        :type order: int | None
        :param opts: extra options forwarded to `get_polarizability_evaluator`
        :type opts: dict
        :return: a function `(coords, order=order) -> polarizability_or_expansion`
        :rtype: callable
        """
        evaluator = self.get_polarizability_evaluator(evaluator, **opts)
        def evaluate_polarizability(coords, order=order):
            """
            **LLM Docstring**

            Evaluate the dipole polarizability (and, if `order` is given, its derivatives) at `coords` using the evaluator from the enclosing scope, defaulting to this molecule's own coordinates if none are given.

            :param coords: the coordinates to evaluate at, in Bohr; defaults to `self.coords`
            :type coords: np.ndarray | None
            :param order: the highest derivative order to evaluate; if `None`, only the polarizability itself is returned
            :type order: int | None
            :return: the polarizability, or the full polarizability/derivative expansion
            :rtype: np.ndarray | list
            """
            smol = order is None
            if smol: order = 0
            if coords is None:
                coords = self.coords
            expansion = evaluator.evaluate(
                np.asanyarray(coords)
                    * UnitsData.convert("BohrRadius", evaluator.distance_units),
                order=order
            )
            if smol: expansion = expansion[0]
            return expansion
        return evaluate_polarizability
    def calculate_dipole_polarizability(self, evaluator=None, order=None, **opts):
        """
        **LLM Docstring**

        Compute the dipole polarizability (and, optionally, its derivatives) at this molecule's current coordinates using the resolved polarizability evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the highest derivative order to compute; if `None`, only the zeroth-order term of each returned quantity is kept
        :type order: int | None
        :param opts: extra options forwarded to `get_polarizability_evaluator`
        :type opts: dict
        :return: the polarizability expansion
        :rtype: list
        """
        evaluator = self.get_polarizability_evaluator(evaluator, **opts)
        smol = order is None
        if smol: order = 0
        expansion = evaluator.evaluate(
            self.coords * UnitsData.convert("BohrRadius", evaluator.distance_units),
            order=order
        )
        if smol: expansion = [e[0] for e in expansion]
        return expansion
    @property
    def polarizability_derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the dipole-polarizability derivative tensors, delegating to `self.dipole_surface.polarizability_derivatives`.

        :param derivs: (setter only) the new polarizability derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the polarizability derivative tensors
        :rtype: list[np.ndarray]
        """
        return self.dipole_surface.polarizability_derivatives
    @polarizability_derivatives.setter
    def polarizability_derivatives(self, derivs):
        """
        **LLM Docstring**

        Property getter/setter for the dipole-polarizability derivative tensors, delegating to `self.dipole_surface.polarizability_derivatives`.

        :param derivs: (setter only) the new polarizability derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the polarizability derivative tensors
        :rtype: list[np.ndarray]
        """
        self.dipole_surface.polarizability_derivatives = derivs

    def get_reduced_potential_generator(self):
        """
        **LLM Docstring**

        Build a `ReducedDimensionalPotentialHandler` for generating reduced-dimensionality potential slices of this molecule.

        :return: the constructed handler
        :rtype: ReducedDimensionalPotentialHandler
        """
        return ReducedDimensionalPotentialHandler(self)
    def get_1d_potentials(self,
                          spec,
                          evaluator=None,
                          energy_expansion=None,
                          potential_params=None,
                          **opts):
        """
        **LLM Docstring**

        Generate 1D potential-energy slices along the specified coordinate(s), via a `ReducedDimensionalPotentialHandler`.

        :param spec: the coordinate specification(s) to slice along
        :type spec: object
        :param evaluator: an explicit energy evaluator to use for the potential evaluations
        :type evaluator: object | None
        :param energy_expansion: a precomputed energy expansion to use instead of evaluating fresh
        :type energy_expansion: list[np.ndarray] | None
        :param potential_params: extra parameters controlling the potential generation
        :type potential_params: dict | None
        :param opts: extra options forwarded to `ReducedDimensionalPotentialHandler.get_1d_potentials`
        :type opts: dict
        :return: the generated 1D potentials
        :rtype: object
        """
        pot_gen = self.get_reduced_potential_generator()
        return pot_gen.get_1d_potentials(
            spec,
            evaluator=evaluator,
            energy_expansion=energy_expansion,
            potential_params=potential_params,
            **opts
        )

    def evaluate(self,
                 func,
                 use_internals=None,
                 order=None,
                 strip_embedding=False
                 ):
        """
        **LLM Docstring**

        Evaluate an arbitrary function of the molecule's coordinates (and its derivatives) via the molecule's `MolecularEvaluator`.

        :param func: the function to evaluate
        :type func: callable
        :param use_internals: whether to evaluate/differentiate in internal coordinates rather than Cartesians
        :type use_internals: bool | None
        :param order: the highest derivative order to compute
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the evaluated function value(s)/derivatives
        :rtype: object
        """
        return self.evaluator.evaluate(
            func,
            use_internals=use_internals,
            order=order,
            strip_embedding=strip_embedding
        )
    def evaluate_at(self,
                    func,
                    coords,
                    use_internals=None,
                    order=None,
                    strip_embedding=False
                    ):
        """
        **LLM Docstring**

        Evaluate an arbitrary function of the molecule's coordinates (and its derivatives) at a specific set of coordinates, via the molecule's `MolecularEvaluator`.

        :param func: the function to evaluate
        :type func: callable
        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param use_internals: whether to evaluate/differentiate in internal coordinates rather than Cartesians
        :type use_internals: bool | None
        :param order: the highest derivative order to compute
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the evaluated function value(s)/derivatives
        :rtype: object
        """
        return self.evaluator.evaluate_at(
            func,
            coords,
            use_internals=use_internals,
            order=order,
            strip_embedding=strip_embedding
        )

    def get_displaced_coordinates(self, displacements,
                                  which=None, sel=None, axes=None,
                                  use_internals=False,
                                  coordinate_expansion=None,
                                  strip_embedding=False,
                                  shift=True
                                  ):
        """
        **LLM Docstring**

        Build displaced copies of the molecule's coordinates along specified atoms/axes, via the molecule's `MolecularEvaluator`.

        :param displacements: the displacement values to apply
        :type displacements: np.ndarray
        :param which: which atoms/coordinates to displace
        :type which: object | None
        :param sel: a selection of atoms to restrict the displacement to
        :type sel: Iterable[int] | None
        :param axes: which axes to displace along
        :type axes: object | None
        :param use_internals: whether the displacements are given in internal coordinates rather than Cartesians
        :type use_internals: bool
        :param coordinate_expansion: a coordinate-transformation expansion to apply to the displacements before applying them
        :type coordinate_expansion: list[np.ndarray] | None
        :param strip_embedding: whether the coordinate expansion has had its embedding coordinates stripped
        :type strip_embedding: bool
        :param shift: whether the displacements are relative shifts (`True`) or absolute target values
        :type shift: bool
        :return: the displaced coordinates
        :rtype: np.ndarray
        """
        return self.evaluator.get_displaced_coordinates(
            displacements,
            which=which, sel=sel, axes=axes,
            use_internals=use_internals,
            strip_embedding=strip_embedding,
            coordinate_expansion=coordinate_expansion,
            shift=shift
        )

    def get_scan_coordinates(self,
                             domains,
                             internals=False,
                             modes=None,
                             order=None,
                             which=None, sel=None, axes=None,
                             shift=True,
                             coordinate_expansion=None,
                             strip_embedding=False,
                             return_displacements=False
                             ):
        """
        **LLM Docstring**

        Build a grid of displaced coordinates spanning the given coordinate domains, optionally re-expressed through a normal-mode (or other) coordinate expansion, via the molecule's `MolecularEvaluator`.

        :param domains: the coordinate ranges to scan over
        :type domains: object
        :param internals: whether the scan coordinates are internal coordinates rather than Cartesians
        :type internals: bool
        :param modes: normal modes to scan along instead of raw coordinates; `True` to use this molecule's own modes (with translation/rotation projection disabled)
        :type modes: object | bool | None
        :param order: the Jacobian order to use when building the mode-based coordinate expansion
        :type order: int | None
        :param which: which atoms/coordinates to scan
        :type which: object | None
        :param sel: a selection of atoms to restrict the scan to
        :type sel: Iterable[int] | None
        :param axes: which axes to scan along
        :type axes: object | None
        :param shift: whether the scan values are relative shifts or absolute target values
        :type shift: bool
        :param coordinate_expansion: an explicit coordinate-transformation expansion to combine with the mode-based one (if `modes` is given) or use directly
        :type coordinate_expansion: list[np.ndarray] | None
        :param strip_embedding: whether the coordinate expansion has had its embedding coordinates stripped
        :type strip_embedding: bool
        :param return_displacements: whether to also return the raw displacement values used
        :type return_displacements: bool
        :return: the scan coordinates, or additionally the displacements if `return_displacements` is set
        :rtype: np.ndarray | tuple
        """
        from ..Modes import NormalModes

        if modes is not None:
            if modes is True:
                modes = self.get_normal_modes(project_transrot=False, use_internals=False)
            modes = NormalModes.prep_modes(modes).remove_mass_weighting()
            # if order is None:
            #     base_expansion = [modes.coords_by_modes]
            # else:
            #     order = 2
            base_expansion = ModeEmbedding(self.embedding, modes).get_cartesians_by_internals(order=order)
            if coordinate_expansion is not None:
                coordinate_expansion = nput.tensor_reexpand(coordinate_expansion, base_expansion)
            else:
                coordinate_expansion = base_expansion
        return self.evaluator.get_scan_coordinates(
            domains,
            internals=internals,
            which=which, sel=sel, axes=axes,
            shift=shift, coordinate_expansion=coordinate_expansion,
            strip_embedding=strip_embedding,
            return_displacements=return_displacements
        )

    def get_nearest_displacement_atoms(self,
                                       points,
                                       sel=None, axes=None, weighting_function=None,
                                       return_distances=False
                                       ):
        """
        **LLM Docstring**

        Find the atoms nearest to a set of query points (optionally restricted to a selection/axes, with a custom distance weighting), via the molecule's `MolecularEvaluator`.

        :param points: the query points to find nearest atoms for
        :type points: np.ndarray
        :param sel: a selection of atoms to restrict the search to
        :type sel: Iterable[int] | None
        :param axes: which coordinate axes to compute distances over
        :type axes: object | None
        :param weighting_function: a custom function to weight distances by
        :type weighting_function: callable | None
        :param return_distances: whether to also return the computed distances
        :type return_distances: bool
        :return: the nearest atom indices, or additionally the distances if `return_distances` is set
        :rtype: np.ndarray | tuple
        """
        return self.evaluator.get_nearest_displacement_atoms(
            points,
            sel=sel, axes=axes, weighting_function=weighting_function,
            return_distances=return_distances
        )
    def get_nearest_displacement_coordinates(self,
                                             points,
                                             sel=None, axes=None, weighting_function=None,
                                             modes_nearest=False,
                                             return_distances=False
                                             ):
        """
        **LLM Docstring**

        Find the displaced coordinates nearest to a set of query points, optionally in normal-mode space, via the molecule's `MolecularEvaluator`.

        :param points: the query points to find nearest coordinates for
        :type points: np.ndarray
        :param sel: a selection of atoms to restrict the search to
        :type sel: Iterable[int] | None
        :param axes: which coordinate axes to compute distances over
        :type axes: object | None
        :param weighting_function: a custom function to weight distances by
        :type weighting_function: callable | None
        :param modes_nearest: whether to find the nearest point in normal-mode space rather than Cartesian/internal space
        :type modes_nearest: bool
        :param return_distances: whether to also return the computed distances
        :type return_distances: bool
        :return: the nearest coordinates, or additionally the distances if `return_distances` is set
        :rtype: np.ndarray | tuple
        """
        return self.evaluator.get_nearest_displacement_coordinates(
            points,
            sel=sel, axes=axes, weighting_function=weighting_function,
            modes_nearest=modes_nearest,
            return_distances=return_distances
        )
    def get_nearest_scan_coordinates(self, domains, sel=None, axes=None):
        """
        **LLM Docstring**

        Find the scan-grid coordinates nearest to a set of query domains, via the molecule's `MolecularEvaluator`.

        :param domains: the coordinate domains to build the nearest scan grid for
        :type domains: object
        :param sel: a selection of atoms to restrict the search to
        :type sel: Iterable[int] | None
        :param axes: which coordinate axes to compute distances over
        :type axes: object | None
        :return: the nearest scan coordinates
        :rtype: np.ndarray
        """
        return self.evaluator.get_nearest_scan_coordinates(domains, sel=sel, axes=axes)

    @classmethod
    def _get_atomic_radius(cls, atom_data, radius_type=None):
        """
        **LLM Docstring**

        Resolve the display radius for an atom-data record: the icon radius (falling back to the van der Waals radius if too small), or an explicitly named radius field.

        :param atom_data: the atom-data record to get a radius for
        :type atom_data: dict
        :param radius_type: the specific radius field to use; if `None`, uses the icon radius with a van der Waals fallback
        :type radius_type: str | None
        :return: the resolved atomic radius
        :rtype: float
        """
        if radius_type is None:
            rad = atom_data["IconRadius"]
            if rad < .8:
                rad = atom_data["VanDerWaalsRadius"]
        else:
            rad = atom_data[radius_type]
        return rad
    def plot_molecule_function(self,
                               function,
                               *,
                               axes,
                               sel=None,
                               embed=False,
                               modes_nearest=False,
                               domain=None,
                               domain_padding=1,
                               plot_points=500,
                               weighting_function=None,
                               mask_function=None,
                               mask_value=0,
                               plot_atoms=False,
                               atom_colors=None,
                               atom_radii=None,
                               plotter=None,
                               epilog=None,
                               **plot_options
                               ):
        """
        **LLM Docstring**

        Plot an arbitrary scalar function of the molecule's geometry over a 1D or 2D grid of displaced coordinates, evaluating the function at the grid points nearest to each displacement, optionally masking out-of-range values and overlaying the atoms' positions.

        :param function: the function to evaluate, taking a batch of coordinates and returning scalar values
        :type function: callable
        :param axes: which atom/coordinate axis (or axes, up to 2) to build the plot grid over
        :type axes: int | Iterable[int]
        :param sel: a selection of atoms to restrict the displacement search to
        :type sel: Iterable[int] | None
        :param embed: whether to Eckart-embed the evaluation points before calling `function`
        :type embed: bool
        :param modes_nearest: whether to find the nearest displacement in normal-mode space
        :type modes_nearest: bool
        :param domain: explicit `(min, max)` bounds per axis; computed from the geometry's bounding box (plus `domain_padding`) if not given
        :type domain: np.ndarray | None
        :param domain_padding: padding to add around the automatically computed domain
        :type domain_padding: float | Iterable[float] | None
        :param plot_points: number of grid points per axis
        :type plot_points: int | Iterable[int]
        :param weighting_function: custom distance-weighting function used when finding nearest displacements
        :type weighting_function: callable | None
        :param mask_function: a function `(values, eval_points, dists) -> mask` used to blank out certain grid values
        :type mask_function: callable | None
        :param mask_value: the value to substitute where `mask_function` returns `True`
        :type mask_value: float
        :param plot_atoms: whether to overlay disks at the atoms' projected positions
        :type plot_atoms: bool
        :param atom_colors: per-atom overlay colors; defaults to each atom's icon color
        :type atom_colors: list | None
        :param atom_radii: per-atom overlay radii; defaults to each atom's resolved radius
        :type atom_radii: list | None
        :param plotter: the plotting class/function to use; defaults to `plt.Plot` (1D) or `plt.TriContourPlot` (2D)
        :type plotter: type | callable | None
        :param epilog: extra plot elements to draw on top of the function plot
        :type epilog: list | None
        :param plot_options: extra options forwarded to the plotting call
        :type plot_options: dict
        :return: the resulting plot
        :rtype: object
        :raises NotImplementedError: if the molecule holds more than one geometry configuration
        :raises ValueError: if more than 2 axes are requested
        """

        if self.coords.ndim > 2:
            raise NotImplementedError("function plotting only supported for one structure at a time")

        axes = np.asanyarray(axes)
        if axes.ndim == 0:
            axes = np.array([axes[:]])
        if len(axes) > 2:
            raise ValueError("can only plot up to 2 axes at a time")

        if domain is None:
            domain = Mesh(self.coords[:, axes]).bounding_box
            if domain_padding is not None:
                if isinstance(domain_padding, (int, float, np.integer, np.floating)):
                    domain_padding = [-domain_padding, domain_padding]
                domain_padding = np.asanyarray(domain_padding)
                if domain_padding.ndim == 1:
                    domain_padding = domain_padding[np.newaxis, :]
                domain = domain + domain_padding

        if isinstance(plot_points, (int, np.integer)):
            plot_points = [plot_points] * len(domain)

        grids = []
        for dom, pts in zip(domain, plot_points):
            grids.append(np.linspace(*dom, pts))
        grids = np.array(np.meshgrid(*grids, indexing='xy'))
        grid_points = np.moveaxis(grids, 0, -1).reshape(-1, len(domain))  # vector of points

        eval_points, dists = self.evaluator.get_nearest_displacement_coordinates(
            grid_points,
            axes=axes,
            sel=sel,
            weighting_function=weighting_function,
            modes_nearest=modes_nearest,
            return_distances=True
        )

        if embed:
            eval_points = self.embed_coords(eval_points)

        values = function(eval_points)
        if mask_function is not None:
            mask = mask_function(values, eval_points, dists)
            values[mask] = mask_value

        if plotter is None:
            if len(grids) == 1:
                plotter = plt.Plot
            else:
                plotter = plt.TriContourPlot

        if plot_atoms:
            if epilog is None:
                epilog = []
            ref = self.coords[:, axes]
            atoms = self._ats
            if sel is not None:
                ref = ref[sel, :]
                atoms = [atoms[i] for i in sel]
            if atom_colors is None:
                atom_colors = [None] * len(atoms)
            atom_colors = [at['IconColor'] if c is None else c for c,at in zip(atom_colors, atoms)]
            if atom_radii is None:
                atom_radii = [None] * len(atoms)
            atom_radii = [
                self._get_atomic_radius(at)
                    if c is None else c for c,at in
                zip(atom_radii, atoms)
            ]
            epilog = list(epilog) + [
                plt.Disk(
                    crd,
                    r,
                    color=c
                )
                for crd,r,c in zip(ref, atom_radii, atom_colors)
            ]

        if values.ndim > 1 and isinstance(plotter, type) and issubclass(plotter, (plt.Plot, plt.Plot2D)):
            values = np.moveaxis(values, 0, -1)[..., np.newaxis]
            def plotter(subvals,
                        _baseclass=plotter,
                        _subgrids=tuple(np.moveaxis(grid_points, -1, 0)),
                        method=None,
                        **opts
                        ):
                """
                **LLM Docstring**

                Wrap a batch of function values (for one component of a multi-component/tensor-valued function) back into a call to the original plotting class over the grid coordinates from the enclosing scope.

                :param subvals: the values for this component, to be flattened and plotted
                :type subvals: np.ndarray
                :param _baseclass: the original plotting class to delegate to (captured from the enclosing scope)
                :type _baseclass: type
                :param _subgrids: the per-axis grid coordinate arrays (captured from the enclosing scope)
                :type _subgrids: tuple[np.ndarray]
                :param method: accepted but not used in this method's body
                :type method: object | None
                :param opts: extra options forwarded to the base plotting class
                :type opts: dict
                :return: the constructed plot for this component
                :rtype: object
                """
                return _baseclass(
                    *_subgrids,
                    subvals.reshape(-1),
                    **opts
                )
            return plt.TensorPlot(values, plot_class=plotter, epilog=epilog, **plot_options)
        else:
            return plotter(*np.moveaxis(grid_points, -1, 0), values, epilog=epilog, **plot_options)

    def get_model(self, potential_specs, dipole=None):
        """
        **LLM Docstring**

        Build an analytic `MolecularModel` (a symbolic potential/dipole surface expressed in terms of this molecule's Z-matrix internal coordinates) from a specification of per-coordinate potential (and optional dipole) function contributions.

        :param potential_specs: either a raw potential expression, or a dict mapping internal-coordinate index (or index tuples, for coupled terms) to a function-type specification (a dict with a `'function_type'` key or a single `{function_type: params}` entry) describing that coordinate's contribution
        :type potential_specs: dict | object
        :param dipole: an optional 3-element (x/y/z) sequence of dipole-component specifications, each following the same per-coordinate specification format as `potential_specs`
        :type dipole: Iterable | None
        :return: the constructed analytic molecular model
        :rtype: MolecularModel
        :raises ValueError: if the molecule has no internal coordinates, if its internal-coordinate spec isn't a plain Z-matrix, or if `dipole` isn't given exactly 3 components
        """
        from ..AnalyticModels import MolecularModel

        if self.internals is None:
            raise ValueError("need internal coordinates to generate analytic model")

        internals = self.internals
        if internals['zmatrix'] is None or internals['conversion'] is not None:
            raise ValueError("only plain Z-matrix coordinate specs currently supported")

        eq_vals = self.internal_coordinates
        embedding_coords = [0, 1, 2, 4, 5, 8]
        good_coords = np.setdiff1d(np.arange(3 * len(self.atoms)), embedding_coords)
        eq_vals = eq_vals.flatten()[good_coords]

        zmat = np.array(internals['zmatrix'])
        def canonicalize_coord(coord_index):
            """
            **LLM Docstring**

            Resolve which atoms (and hence coordinate type -- stretch/bend/dihedral) a given flat Z-matrix coordinate index corresponds to, and build the corresponding `MolecularModel` coordinate object.

            :param coord_index: the flat internal-coordinate index (into the stripped, embedding-coordinate-free Z-matrix coordinate vector)
            :type coord_index: int
            :return: `(coord, atoms)` -- the constructed `MolecularModel` coordinate object and the atom indices defining it
            :rtype: tuple
            :raises ValueError: if the resolved atom count doesn't correspond to a stretch, bend, or dihedral
            """
            if coord_index == 0:
                atoms = zmat[1, :2]
            elif coord_index == 1:
                atoms = zmat[2, :2]
            elif coord_index == 2:
                atoms = zmat[2, :3]
            else:
                row = 3 + (coord_index - 3) // 3
                num = 2 + (coord_index - 3) % 3
                atoms = zmat[row, :num]

            if len(atoms) == 2:
                coord_type = 'r'
            elif len(atoms) == 3:
                coord_type = 'a'
            elif len(atoms) == 4:
                coord_type = 't'
            else:
                raise ValueError("bad coord type")
            coord = getattr(MolecularModel, coord_type)(*atoms)

            return coord, atoms

        def canonicalize_spec(coord_index, spec_dict):
            """
            **LLM Docstring**

            Build the analytic potential/dipole contribution function for a single internal coordinate from its specification dict, filling in the coordinate's equilibrium value as the default `'eq'` parameter if not given, and applying any `'scaling'` factor.

            :param coord_index: the internal-coordinate index this specification applies to
            :type coord_index: int
            :param spec_dict: the specification, either containing an explicit `'function_type'` key or (if it has exactly one entry) a single `{function_type: params}` mapping
            :type spec_dict: dict
            :return: the constructed (and scaled) analytic function contribution
            :rtype: object
            """
            if 'function_type' not in spec_dict and len(spec_dict) == 1:
                function_type, spec_dict = next(iter(spec_dict.items()))
            else:
                function_type = spec_dict['function_type']
            atoms = canonicalize_coord(coord_index)[1]

            scaling = spec_dict.get('scaling', 1)
            params = spec_dict.copy()
            if 'scaling' in params:
                del params['scaling']
            if 'function_type' in params:
                del params['function_type']

            if 'eq' not in params:
                params['eq'] = eq_vals[coord_index]

            return scaling * getattr(MolecularModel, function_type)(*atoms, **params)

        coord_indices = set()

        if isinstance(potential_specs, dict):
            potential_contribs = []
            for idx, spec in potential_specs.items():
                if nput.is_int(idx):
                    idx = [idx]
                    spec = [spec]
                coord_indices.update(idx)
                fns = [canonicalize_spec(i, s) for i,s in zip(idx, spec)]
                f = fns[0]
                for fn in fns[1:]:
                    f = f * fn
                potential_contribs.append(f)
            pot = sum(potential_contribs)
        else:
            pot = potential_specs

        if dipole is not None:
            if len(dipole) != 3:
                raise ValueError("need xyz for dipole contribs")

            dip = []
            for d in dipole:
                if not isinstance(d, dict):
                    dip.append(d)
                    continue

                contribs = []
                for idx, spec in d.items():
                    if nput.is_int(idx):
                        idx = [idx]
                        spec = [spec]
                    coord_indices.update(idx)
                    fns = [canonicalize_spec(i, s) for i, s in zip(idx, spec)]
                    f = fns[0]
                    for fn in fns[1:]:
                        f = f * fn
                    contribs.append(f)
                dip.append(sum(contribs))
        else:
            dip = None

        atoms = set()
        coords = []
        vals = {}
        for idx in sorted(list(coord_indices)):
            c, a = canonicalize_coord(idx)
            atoms.update(a)
            coords.append(c)
            vals[c] = eq_vals[idx]
        masses = self._atomic_masses()

        vals.update({
                MolecularModel.m(i): masses[i]
                for i in sorted(list(atoms))
            })

        return MolecularModel(
            self,
            coords,
            pot,
            dipole=dip,
            values=vals
        )
    def get_point_group(self, coords=None, masses=None, *, sel=None, verbose=False, return_identifier=False, **tols):
        """
        **LLM Docstring**

        Identify the molecular point group at the given (or this molecule's own) geometry, via `McUtils.Symmetry.PointGroupIdentifier`.

        :param coords: the coordinates to identify the point group for; defaults to `self.coords`
        :type coords: np.ndarray | None
        :param masses: the masses to use; defaults to `self.atomic_masses`
        :type masses: np.ndarray | None
        :param sel: a selection of atoms to restrict the identification to
        :type sel: Iterable[int] | None
        :param verbose: whether the identifier should print diagnostic output
        :type verbose: bool
        :param return_identifier: whether to also return the underlying `PointGroupIdentifier` object
        :type return_identifier: bool
        :param tols: extra tolerance options forwarded to `PointGroupIdentifier`
        :type tols: dict
        :return: the identified point group, or `(identifier, point_group)` if `return_identifier` is set
        :rtype: object | tuple
        """
        if coords is None:
            coords = self.coords
        else:
            coords = np.asanyarray(coords)
        if masses is None:
            masses = self.atomic_masses
        else:
            masses = np.asanyarray(masses)
        if sel is not None:
            coords = coords[sel,]
            masses = masses[sel,]
        pg_id = symm.PointGroupIdentifier(coords, masses, verbose=verbose, **tols)
        _, pg = pg_id.identify_point_group()
        if return_identifier:
            return pg_id, pg
        return pg

    @property
    def point_group(self) -> symm.PointGroup:
        """
        **LLM Docstring**

        Property getter/setter for the molecule's (cached) point group. The getter computes it lazily via `get_point_group` the first time it's needed.

        :param point_group: (setter only) the new point group to cache
        :type point_group: symm.PointGroup
        :return: (getter) the cached (or newly computed) point group
        :rtype: symm.PointGroup
        """
        if self._pg is None:
            self._pg = self.get_point_group()
        return self._pg
    @point_group.setter
    def point_group(self, point_group):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's (cached) point group. The getter computes it lazily via `get_point_group` the first time it's needed.

        :param point_group: (setter only) the new point group to cache
        :type point_group: symm.PointGroup
        :return: (getter) the cached (or newly computed) point group
        :rtype: symm.PointGroup
        """
        self._pg = point_group

    # def get_standard_orientation(self,
    #                              coords=None, masses=None,
    #                              *,
    #                              point_group=None,
    #                              origin=None,
    #                              moments_of_inertia=None,
    #                              axes=None
    #                              ):
    #     default = coords is None and masses is None
    #     if coords is None:
    #         coords = self.coords
    #     if masses is None:
    #         masses = self.masses
    #     if origin is None:
    #         if default:
    #             origin = self.center_of_mass
    #         else:
    #             origin = nput.center_of_mass(coords, masses)
    #     moms = moments_of_inertia
    #     if axes is None:
    #         if default:
    #             moms, axes = self.inertial_eigensystem
    #         else:
    #             moms, axes = nput.moments_of_inertia(coords, masses)
    #
    #     point_group:symm.PointGroup # help out pycharm...
    #     if point_group is None:
    #         if default:
    #             point_group = self.point_group
    #         else:
    #             point_group = self.get_point_group(coords=coords, masses=masses)
    #
    #     if moms is None:
    #         if default:
    #             moms = self.moments_of_inertia
    #         else:
    #             moms, _ = nput.moments_of_inertia(coords, masses)
    #     rotor_type, planar = symm.identify_rotor_type(moms)
    #
    #     def _handle_spherical_tops(point_group, coords, rotor_type, axes, masses):
    #         if rotor_type == symm.RotorTypes.Prolate: # a axis unique
    #             axes = axes[:, (1, 2, 0)]
    #
    #         # find Gaussian's "circular sets"
    #         test_coords = coords @ axes
    #         axis_dists = np.linalg.norm(test_coords[:, 2], axis=1)
    #         (_, groups), _ = nput.group_by(np.arange(len(coords)), np.round(masses))
    #         groups = [
    #             g
    #             for gg in groups
    #             for g in nput.group_by(gg, np.round(test_coords[gg, 3], 1))[0][1]
    #         ]
    #         groups = [
    #             g
    #             for gg in groups
    #             for g in nput.group_by(gg, np.round(axis_dists[gg,], 1))[0][1]
    #         ]
    #         key_group = min(groups,
    #                         lambda g: (
    #                             abs(test_coords[g[0], 3]),
    #                             -test_coords[g[0], 3],
    #                             axis_dists[g[0]],
    #                             masses[g[0]]
    #                         ))
    #         key_atom = min(key_group)
    #
    #
    #     if point_group:
    #         if (
    #                 point_group.group.value == "Cv"
    #                 and point_group.n == 2
    #         ):
    #             if not planar:
    #                 raise NotImplementedError("non-planar C2v structures not handled")
    #         elif planar or (
    #                 point_group.group.value == "Dh"
    #                 and point_group.n == 2
    #         ):
    #             ...
    #         elif (
    #                 point_group.group.value == "Ci"
    #                 or (
    #                         point_group.group.value == "C"
    #                         and point_group.n == 1
    #                 )
    #         ):
    #             return coords - origin[np.newaxis], MolecularTransformation(-origin)
    #
    #     return MolecularTransformation(-origin, axes)

    def get_point_group_embedded_coordinates(self,
                                             pg=None, sel=None,
                                             return_point_group=False,
                                             return_identifier=False,
                                             **tols):
        """
        **LLM Docstring**

        Compute this molecule's coordinates re-embedded into the standard reference frame of its molecular point group (translating to the center of mass, aligning to the principal axes, then aligning to the point group's own axis convention).

        :param pg: an explicit point group to embed into; identified automatically via `get_point_group` if not given
        :type pg: symm.PointGroup | None
        :param sel: a selection of atoms to restrict the point-group identification to
        :type sel: Iterable[int] | None
        :param return_point_group: whether to also return the (embedded) point group
        :type return_point_group: bool
        :param return_identifier: whether to also return the underlying `PointGroupIdentifier`
        :type return_identifier: bool
        :param tols: extra tolerance options forwarded to `PointGroupIdentifier`
        :type tols: dict
        :return: the point-group-embedded coordinates, or additionally the point group and/or identifier depending on the `return_*` flags
        :rtype: np.ndarray | tuple
        """
        coords = self.coords
        masses = self.atomic_masses
        if sel is not None:
            coords = coords[sel,]
            masses = masses[sel,]

        com = nput.center_of_mass(coords, masses)
        coords = coords - com[np.newaxis]
        _, axes = nput.moments_of_inertia(coords, masses)
        coords = coords @ axes

        if no_pg := pg is None:
            pg_id, pg = self.get_point_group(sel=sel,  **tols, return_identifier=True)
        else:
            # make sure this stays in sync with above
            pg_id = symm.PointGroupIdentifier(coords, masses, **tols)
            pg = pg_id.embed_point_group(pg)

        emb_coords = (self.coords - com[np.newaxis]) @ pg.axes.T @ pg.base_axes
        if return_point_group or return_identifier:
            res = (emb_coords,)
            if return_point_group:
                res = res + (pg,)
            if return_identifier:
                res = res + (pg_id,)
        else:
            res = emb_coords

        return res

    def symmetrize(self, pg=None,
                   return_identifier=False,
                   tol=1e-1,
                   sel=None,
                   return_coordinates=None,
                   return_point_group=False,
                   **tols):
        """
        **LLM Docstring**

        Symmetrize the molecule's geometry to conform exactly to a (identified or given) point group, adding/relabeling atoms as needed if the symmetrization changes the atom count, and matching the symmetrized geometry back onto the original atom ordering when the count is unchanged.

        :param pg: an explicit point group to symmetrize to; identified automatically via `get_point_group` if not given
        :type pg: symm.PointGroup | None
        :param return_identifier: whether to also return the underlying `PointGroupIdentifier`
        :type return_identifier: bool
        :param tol: tolerance used for point-group identification and symmetrization
        :type tol: float
        :param sel: a selection of atoms to restrict the symmetrization to
        :type sel: Iterable[int] | None
        :param return_coordinates: whether to return raw coordinates instead of a modified `Molecule` (when the atom count is unchanged)
        :type return_coordinates: bool | None
        :param return_point_group: whether to also return the point group used
        :type return_point_group: bool
        :param tols: extra tolerance options forwarded to `PointGroupIdentifier`
        :type tols: dict
        :return: a tuple starting with the symmetrized molecule/coordinates (or `None` if the atom count changed and no direct match was possible), followed by `(new_atoms, new_coords)`, and optionally the identifier and/or point group
        :rtype: tuple
        """

        coords = self.coords
        masses = self.atomic_masses
        labels = self.atoms
        if sel is not None:
            coords = coords[sel,]
            masses = masses[sel,]
            labels = [labels[s] for s in sel]

        com = nput.center_of_mass(coords, masses)
        coords = coords - com[np.newaxis]
        _, axes = nput.moments_of_inertia(coords, masses)
        coords = coords @ axes

        if no_pg := pg is None:
            pg_id, pg = self.get_point_group(sel=sel, tol=tol, **tols, return_identifier=True)
        else:
            # make sure this stays in sync with above
            pg_id = symm.PointGroupIdentifier(coords, masses, tol=tol, **tols)
            pg = pg_id.embed_point_group(pg)

        new_coords, new_atoms = symm.symmetrize_structure(
            coords,
            pg,
            labels=labels,
            groups=pg_id.groups,
            tol=tol
        )

        new_coords = new_coords @ axes.T + com[np.newaxis]
        if len(coords) == len(new_coords): # nothing added, find best match to old coords
            perm = nput.find_coordinate_matching_permutation(new_coords, coords)
            all_coords = new_coords[perm,]
            if sel is not None:
                all_coords = self.coords.copy()
                all_coords[sel,] = new_coords
            if return_coordinates:
                res = (all_coords,)
            else:
                res = (self.modify(coords=all_coords),)
        else:
            res = (None,)

        res = res + ((new_atoms, new_coords),)
        if return_identifier:
            res = res + (pg_id,)
        if return_point_group:
            res = res + (pg,)
        return res

    def get_symmetrized_internals(self,
                                  point_group=None,
                                  *,
                                  internals=None,
                                  extra_internals=None,
                                  masses=None,
                                  return_expansions=False,
                                  atom_selection=None,
                                  as_characters=True,
                                  normalize=None,
                                  drop_empty_modes=None,
                                  perms=None,
                                  return_base_expansion=False,
                                  return_point_group=False,
                                  reduce_redundant_coordinates=None,
                                  ops=None,
                                  permutation_tol=1e-2,
                                  **opts
                                  ):
        """
        **LLM Docstring**

        Build symmetry-adapted linear combinations of internal coordinates for the molecule under a (identified or given) point group, via `McUtils.Symmetry.symmetrize_internals`.

        :param point_group: an explicit point group to symmetrize under; identified automatically via `get_point_group` (or on a submolecule, if `atom_selection` is given) if not given
        :type point_group: symm.PointGroup | None
        :param internals: the internal coordinates to symmetrize; extracted from `self.internals` (Z-matrix or generic spec) if not given
        :type internals: Iterable | None
        :param extra_internals: extra internal coordinates to prepend to the extracted set (deduplicated)
        :type extra_internals: Iterable | None
        :param masses: masses to use; defaults to `self.atomic_masses`
        :type masses: np.ndarray | None
        :param return_expansions: whether to also return the coordinate-derivative expansions used
        :type return_expansions: bool
        :param atom_selection: restrict the point-group identification to a submolecule
        :type atom_selection: Iterable[int] | None
        :param as_characters: whether to express the symmetrized coordinates in terms of irrep characters
        :type as_characters: bool
        :param normalize: whether to normalize the symmetrized coordinate combinations
        :type normalize: bool | None
        :param drop_empty_modes: whether to drop symmetry combinations that come out identically zero
        :type drop_empty_modes: bool | None
        :param perms: explicit atom permutations to use for the symmetry operations
        :type perms: Iterable | None
        :param return_base_expansion: whether to also return the base (pre-symmetrization) coordinate expansion
        :type return_base_expansion: bool
        :param return_point_group: whether to also return the point group used
        :type return_point_group: bool
        :param reduce_redundant_coordinates: whether to reduce a redundant coordinate set down to a non-redundant symmetrized basis
        :type reduce_redundant_coordinates: bool | None
        :param ops: explicit symmetry operations to use instead of the full point-group operation set
        :type ops: Iterable | None
        :param permutation_tol: tolerance used when matching atom permutations to symmetry operations
        :type permutation_tol: float
        :param opts: extra options forwarded to point-group identification/`symmetrize_internals`
        :type opts: dict
        :return: the symmetrized internal coordinates, plus any additionally requested return values, prefixed with the point group if `return_point_group` is set
        :rtype: object | tuple
        :raises ValueError: if no internal coordinates are available or given
        """
        if internals is None:
            if self.internals is None:
                raise ValueError("can't get symmetrized internals without an initial internal set")

            if self.internals['zmatrix'] is None:
                internals = self.internals['spec']
            else:
                internals = coordops.extract_zmatrix_internals(self.internals['zmatrix'])

            if extra_internals is not None:
                internals = list(itut.delete_duplicates(list(extra_internals) + list(internals)))

        if point_group is None:
            if atom_selection is not None:
                point_group = self.take_submolecule(atom_selection).get_point_group(**opts)
            else:
                point_group = self.get_point_group(**opts)
        else:
            # make sure this stays in sync with above
            pg_id = symm.PointGroupIdentifier(self.coords, self.atomic_masses, **opts)
            point_group = pg_id.embed_point_group(point_group)

        if masses is None:
            masses = self.atomic_masses

        res = symm.symmetrize_internals(
            point_group, internals,
            self.coords,
            masses=masses,
            return_expansions=return_expansions,
            atom_selection=atom_selection,
            as_characters=as_characters,
            normalize=normalize,
            drop_empty_modes=drop_empty_modes,
            reduce_redundant_coordinates=reduce_redundant_coordinates,
            perms=perms,
            return_base_expansion=return_base_expansion,
            ops=ops,
            permutation_tol=permutation_tol,
            **opts
        )

        if return_point_group:
            res = (point_group,) + res

        return res

    def get_surface(self,
                    radius_type='VanDerWaalsRadius',
                    *,
                    surface_type=None,
                    radius_units="Angstroms",
                    samples=100,
                    radius_scaling=1,
                    **etc):
        """
        **LLM Docstring**

        Build a surface representation of the molecule (by default, a union of per-atom spheres) using per-atom radii of the requested type.

        :param radius_type: the `AtomData` radius field to use for each atom
        :type radius_type: str
        :param surface_type: the surface class to construct; defaults to `zach.SphereUnionSurface`
        :type surface_type: type | None
        :param radius_units: the units the radii are given in before converting to Bohr
        :type radius_units: str
        :param samples: number of samples used to represent the surface
        :type samples: int
        :param radius_scaling: uniform scale factor applied to all radii
        :type radius_scaling: float
        :param etc: extra options forwarded to the surface constructor
        :type etc: dict
        :return: the constructed surface
        :rtype: object
        """
        if surface_type is None:
            surface_type = zach.SphereUnionSurface

        radii = np.array([
            AtomData[a, radius_type] * UnitsData.convert(radius_units, "BohrRadius") * radius_scaling
            for a in self.atoms
        ])

        return surface_type(
            self.coords,
            radii,
            samples=samples,
            **etc
        )

    def get_surface_mesh(self,
                         radius_type='VanDerWaalsRadius',
                         *,
                         surface_type=None,
                         radius_units="Angstroms",
                         samples=50,
                         expansion=.01,
                         mesh_options=None,
                         **etc
                         ):
        """
        **LLM Docstring**

        Build a triangulated mesh of the molecule's surface, via `get_surface` followed by `generate_mesh`.

        :param radius_type: the `AtomData` radius field to use for each atom
        :type radius_type: str
        :param surface_type: the surface class to construct; defaults to `zach.SphereUnionSurface`
        :type surface_type: type | None
        :param radius_units: the units the radii are given in before converting to Bohr
        :type radius_units: str
        :param samples: number of samples used to represent the surface
        :type samples: int
        :param expansion: expansion factor applied to the surface before meshing
        :type expansion: float
        :param mesh_options: extra options forwarded to `generate_mesh`
        :type mesh_options: dict | None
        :param etc: extra options forwarded to `get_surface`
        :type etc: dict
        :return: the generated surface mesh
        :rtype: object
        """

        if mesh_options is None:
            mesh_options = {}

        surf:zach.SphereUnionSurface = self.get_surface(
            radius_type=radius_type,
            surface_type=surface_type,
            radius_units=radius_units,
            samples=samples,
            expansion=expansion,
            **etc
        )

        return surf.generate_mesh(**mesh_options)

    def setup_AIMD(self,
                   potential_function=None,
                   timestep=.5,
                   seed=None,
                   total_energy=None,
                   total_energy_scaling=None,
                   trajectories=1,
                   sampled_modes=None,
                   initial_energies=None,
                   initial_displacements=None,
                   initial_mode_directions=None,
                   displaced_coords=None,
                   track_kinetic_energy=False,
                   track_velocities=False,
                   modes=None,
                   **etc
                   ):
        """
        **LLM Docstring**

        Set up an `AIMDSimulator` (ab initio / classical molecular dynamics) for this molecule, either propagating from an explicit set of displaced initial positions, or sampling initial normal-mode velocities/energies (from explicit mode directions, explicit per-trajectory energies, or randomly-sampled directions distributing a target total energy) when starting from the equilibrium geometry.

        :param potential_function: the energy/gradient function to use; built via `get_energy_function` if not given
        :type potential_function: callable | None
        :param timestep: the simulation timestep
        :type timestep: float
        :param seed: random seed for sampling initial mode directions
        :type seed: int | None
        :param total_energy: the total vibrational energy to distribute across sampled trajectories/modes
        :type total_energy: float | None
        :param total_energy_scaling: scale factor applied to the default total energy (sum of mode frequencies) if `total_energy` isn't given
        :type total_energy_scaling: float | None
        :param trajectories: number of trajectories to sample when randomly generating mode directions
        :type trajectories: int
        :param sampled_modes: which normal modes to sample energy into; defaults to all of them
        :type sampled_modes: Iterable[int] | None
        :param initial_energies: explicit per-trajectory, per-mode energies to use instead of sampling them
        :type initial_energies: np.ndarray | None
        :param initial_displacements: explicit initial Cartesian displacements to start the trajectories from, bypassing the normal-mode energy-sampling path entirely
        :type initial_displacements: np.ndarray | None
        :param initial_mode_directions: explicit per-trajectory mode-direction vectors to convert into initial energies
        :type initial_mode_directions: np.ndarray | None
        :param displaced_coords: which coordinates `initial_displacements` applies to
        :type displaced_coords: object | None
        :param track_kinetic_energy: whether the simulator should track kinetic energy over the trajectory
        :type track_kinetic_energy: bool
        :param track_velocities: whether the simulator should track velocities over the trajectory
        :type track_velocities: bool
        :param modes: normal modes to use for the energy-to-velocity conversion; computed via `get_normal_modes` if not given
        :type modes: object | None
        :param etc: extra options forwarded to the `AIMDSimulator` constructor
        :type etc: dict
        :return: the constructed simulator
        :rtype: AIMDSimulator
        :raises ValueError: if both `initial_energies` and `initial_mode_directions` are given
        """
        from ..AIMD import AIMDSimulator

        if potential_function is None:
            potential_function = self.get_energy_function()

        if initial_displacements is not None:
            init_pos = self.get_displaced_coordinates(
                initial_displacements,
                which=displaced_coords,
                internals='reembed'
            )
            sim = AIMDSimulator(
                self.masses,
                init_pos,
                lambda c: -potential_function(c, order=1)[1].reshape(c.shape),
                timestep=timestep,
                track_kinetic_energy=track_kinetic_energy,
                track_velocities=track_velocities
            )
        else:
            new = self.modify(potential_derivatives=potential_function(self.coords, order=2)[1:])

            if initial_mode_directions is not None:
                if initial_energies is not None:
                    raise ValueError("definitions for both `initial_energies` and `initial_mode_directions`")
                freqs = new.normal_modes.modes.freqs
                if sampled_modes is None:
                    sampled_modes = list(range(freqs.shape[0]))

                initial_energies = np.zeros((len(initial_mode_directions), freqs.shape[0]))
                initial_mode_directions = np.asanyarray(initial_mode_directions)
                subdirs = initial_mode_directions[:, sampled_modes] * freqs[sampled_modes,][np.newaxis]
                initial_energies[:, sampled_modes] = subdirs

                if total_energy is not None:
                    dirs = initial_energies
                    dirs = dirs / np.sum(np.abs(dirs), axis=1)[:, np.newaxis]  # weights in each dimension
                    energies = dirs * freqs[np.newaxis, :]
                    initial_energies = total_energy * energies / np.sum(np.abs(energies), axis=1)[:, np.newaxis]

            if initial_energies is None:
                if modes is None:
                    modes = new.get_normal_modes()
                freqs = np.abs(modes.freqs)

                if total_energy is None:
                    if total_energy_scaling is None:
                        total_energy_scaling = 1/2
                    total_energy = np.sum(freqs) * total_energy_scaling

                if seed is not None:
                    np.random.seed(seed)
                if sampled_modes is None:
                    sampled_modes = list(range(freqs.shape[0]))
                subdirs = np.random.normal(0, 1, size=(trajectories, len(sampled_modes)))

                dirs = np.zeros((trajectories, freqs.shape[0]))
                dirs[:, sampled_modes] = subdirs
                dirs = dirs / np.linalg.norm(dirs, axis=1)[:, np.newaxis] # random unbiased directions

                dirs = dirs / np.sum(np.abs(dirs), axis=1)[:, np.newaxis] # weights in each dimension
                energies = dirs * freqs[np.newaxis, :]
                initial_energies = total_energy * energies / np.sum(np.abs(energies), axis=1)[:, np.newaxis]

            nms = new.get_normal_modes(use_internals=False, mass_weighted=False)
            sim = AIMDSimulator(
                self.atomic_masses,
                [self.coords] * len(initial_energies),
                lambda c: -potential_function(c, order=1)[1].reshape(c.shape),
                velocities=AIMDSimulator.mode_energies_to_velocities(
                    nms.coords_by_modes,
                    self.atomic_masses,
                    initial_energies,
                    inverse=nms.modes_by_coords
                ),
                timestep=timestep,
                track_kinetic_energy=track_kinetic_energy,
                track_velocities=track_velocities,
                **etc
            )

        return sim

    def setup_VPT(self,
                  *,
                  states=2,
                  order=2,
                  use_internals=None,
                  potential_derivatives=None,
                  energy_evaluator=None,
                  dipole_derivatives=None,
                  dipole_evaluator=None,
                  runner='matrix',
                  use_reaction_path=False,
                  modes=None,
                  projected_modes=None,
                  mode_transformation=None,
                  potential_terms=None,
                  dipole_terms=None,
                  **opts
                  ):
        """
        **LLM Docstring**

        Set up a vibrational perturbation theory (VPT) calculation for this molecule, resolving the potential/dipole derivative expansions (computing them via the configured evaluators if not already available), the normal modes (optionally localized against reaction-path-projected modes), and dispatching to the appropriate VPT runner (`VPTRunner` for the default matrix-based approach, or `AnalyticVPTRunner`).

        :param states: the target vibrational states specification (e.g. max quanta) to compute
        :type states: int | object
        :param order: the perturbation-theory order to compute the potential/dipole expansions to
        :type order: int
        :param use_internals: whether to run VPT in internal coordinates (constructing from `self.modify()`) rather than bare atoms/coordinates
        :type use_internals: bool | None
        :param potential_derivatives: explicit potential-energy derivative tensors to use instead of computing them
        :type potential_derivatives: list[np.ndarray] | None
        :param energy_evaluator: an explicit energy evaluator to use when computing the potential derivatives
        :type energy_evaluator: object | None
        :param dipole_derivatives: explicit dipole derivative tensors to use instead of computing them
        :type dipole_derivatives: list[np.ndarray] | None
        :param dipole_evaluator: an explicit dipole evaluator to use when computing the dipole derivatives
        :type dipole_evaluator: object | None
        :param runner: which VPT runner to use: `'matrix'` for `VPTRunner`, or an explicit runner class/object
        :type runner: str | object
        :param use_reaction_path: whether to project the normal modes against reaction-path modes
        :type use_reaction_path: bool
        :param modes: explicit normal modes to use instead of computing them
        :type modes: object | None
        :param projected_modes: explicit modes to localize the normal modes against
        :type projected_modes: object | None
        :param mode_transformation: an explicit mode-coordinate transformation to use
        :type mode_transformation: object | None
        :param potential_terms: precomputed potential expansion terms, bypassing derivative computation entirely
        :type potential_terms: object | None
        :param dipole_terms: precomputed dipole expansion terms, bypassing derivative computation entirely
        :type dipole_terms: object | None
        :param opts: extra options forwarded to the VPT runner's `construct` method
        :type opts: dict
        :return: the constructed VPT runner/calculation
        :rtype: object
        """
        from ..VPT2 import VPTRunner, AnalyticVPTRunner
        if dev.str_is(runner, 'matrix'):
            runner = VPTRunner
        elif not hasattr(runner, 'construct'):
            runner = AnalyticVPTRunner

        og_pot_der = potential_derivatives
        if potential_terms is None:
            if potential_derivatives is None:
                potential_derivatives = self.get_cartesian_potential_derivatives(
                    evaluator=energy_evaluator,
                    order=order+2
                )

        if modes is None:
            modes = self.get_normal_modes(
                potential_derivatives=og_pot_der,
                use_internals=False, project_transrot=False
            )
            if use_reaction_path:
                rpnms, rpnm_status = self.get_reaction_path_modes(
                    potential_derivatives=og_pot_der,
                    return_status=True
                )
                if rpnms.status:
                    projected_modes = rpnms[1:]

        if projected_modes is not None:
            #TODO: integrate this into the VPT infrastructure directly...
            from ..Modes import NormalModes

            if modes.mass_weighted:
                projected_modes = modes.make_mass_weighted()
            modes = modes.localize(target_modes=projected_modes.modes_by_coords)
            # modes = NormalModes(
            #     loc.basis,
            #     loc.matrix,
            #     inverse=loc.inverse,
            #     masses=loc.masses,
            #     freqs=loc.freqs,
            #     mass_weighted=loc.mass_weighted
            # )[list(range(1, len(loc.freqs)))]

        if dipole_terms is None:
            if dipole_derivatives is None:
                dipole_derivatives = self.get_cartesian_dipole_derivatives(
                    evaluator=dipole_evaluator,
                    order=order + 1,
                    include_constant_term=True
                )

        if use_internals or use_internals is None:
            return runner.construct(self.modify(),
                                    states,
                                    order=order,
                                    potential_derivatives=potential_derivatives,
                                    dipole_derivatives=dipole_derivatives,
                                    modes=modes,
                                    mode_transformation=mode_transformation,
                                    potential_terms=potential_terms,
                                    dipole_terms=dipole_terms,
                                    **opts
                                    )
        else:
            return runner.construct(
                [self.atoms, self.coords],
                states,
                potential_derivatives=potential_derivatives,
                dipole_derivatives=dipole_derivatives,
                modes=modes,
                mode_transformation=mode_transformation,
                order=order,
                potential_terms=potential_terms,
                dipole_terms=dipole_terms,
                **opts
            )

    def setup_job(self, job_type, *args, **kwargs):
        """
        **LLM Docstring**

        Set up an external-program computational job (e.g. a Gaussian/ORCA job) for this molecule, via `ExternalProgramJob.resolve` and `from_mol`.

        :param job_type: the job type to resolve (a name or job-type object)
        :type job_type: str | object
        :param args: positional arguments forwarded to the resolved job type's `from_mol`
        :type args: tuple
        :param kwargs: keyword arguments forwarded to the resolved job type's `from_mol`, merged over its default options
        :type kwargs: dict
        :return: the constructed job
        :rtype: object
        """
        job_type, opts = ExternalProgramJob.resolve(job_type)
        return job_type.from_mol(
            self,
            *args,
            **(opts | kwargs)
        )

    def get_gmatrix(self,
                    masses=None, coords=None, use_internals=None, power=None,
                    **internals_opts
                    ):
        """
        **LLM Docstring**

        Compute the G-matrix (inverse effective-mass metric tensor) for this molecule, either the trivial diagonal inverse-mass form (in plain Cartesians) or the projected `B G0 B^T`-style form built from the internals-by-Cartesians Jacobian (when using internal coordinates), optionally raised to a fractional power.

        :param masses: masses to use instead of this molecule's own
        :type masses: np.ndarray | None
        :param coords: alternate coordinates to compute the G-matrix at
        :type coords: np.ndarray | None
        :param use_internals: whether to compute the G-matrix in internal coordinates; defaults to `True` if the molecule has internal coordinates defined
        :type use_internals: bool | None
        :param power: an optional power to raise the resulting G-matrix to
        :type power: float | None
        :param internals_opts: extra options forwarded to `get_internals_by_cartesians` when using internal coordinates
        :type internals_opts: dict
        :return: the G-matrix
        :rtype: np.ndarray
        """
        if use_internals is None:
            use_internals = self.internals is not None

        if not use_internals:
            if masses is None:
                masses = self._atomic_masses()
            else:
                masses = np.asanyarray(masses)
            mass_spec = np.broadcast_to(masses[:, np.newaxis], (len(masses), 3)).flatten()
            if power is not None:
                mass_spec = np.power(mass_spec, power)
            g = np.diag(1 / mass_spec)
            if coords is not None:
                g = np.expand_dims(g, list(range(coords.ndim-2)))
                return np.broadcast_to(g, coords.shape[:-2] + g.shape[-2:])
            else:
                return g
            # raise ValueError("need internal coordinates to calculate the G-matrix")
        else:
            if masses is None:
                masses = self.atomic_masses

            bT = np.tensordot(
                self.get_internals_by_cartesians(1, coords=coords, strip_embedding=True,
                                                 **internals_opts
                                                 )[0],
                np.diag(np.repeat(1 / np.sqrt(masses), 3)),
                axes=[-2, 0]
            )
            # g = self.hamiltonian.gmatrix_expansion(0, masses=masses, coords=coords, modes=None)[0]
            g = bT @ np.moveaxis(bT, -1, -2)
            if power is not None:
                g = nput.fractional_power(g, power)
            return g
    @property
    def g_matrix(self):
        """
        Returns the molecular g-matrix for the system
        :return:
        :rtype:
        """
        return self.get_gmatrix()

    @property
    def coriolis_constants(self):
        """
        Returns the molecular g-matrix for the system
        :return:
        :rtype:
        """
        return self.prop('coriolis_constants')

    def principle_axis_frame(self, sel=None, inverse=False):
        """
        Gets the principle axis frame(s) for the molecule
        :param mol:
        :type mol:
        :param sel: selection of atoms to use when getting the Eckart frame
        :type sel:
        :param inverse: whether to return the inverse of the rotations or not
        :type inverse: bool
        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """
        return self.prop('principle_axis_transformation', sel=sel, inverse=inverse)

    @property
    def principle_axis_data(self):
        """
        Gets the principle axis embedded coords and embedding parameters for the molecule

        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """
        return self.prop('principle_axis_data')

    def permute_atoms(self, perm):
        """
        **LLM Docstring**

        Build a copy of this molecule with its atoms reordered according to a permutation, remapping bonds to the new indexing accordingly.

        :param perm: the new atom ordering, given as the sequence of old indices in their new order
        :type perm: Iterable[int]
        :return: the permuted molecule
        :rtype: Molecule
        """
        inv_perm = np.argsort(perm)
        #TODO: handle applying affine transformation to permute
        return self.modify(
            atoms=[self._ats[i] for i in perm],
            masses=[self._mass[i] for i in perm],
            coords=self.coords[perm, :],
            bonds=[
                [inv_perm[b[0]], inv_perm[b[1]]] + list(b[2:])
                for b in self._bonds
            ] if self._bonds is not None else None
        )

    def apply_affine_transformation(self, transformation, load_properties=False, embed_properties=True):
        """
        **LLM Docstring**

        Apply an affine (rotation/translation) transformation to the molecule's coordinates, optionally propagating the transformation onto its normal modes, potential surface, and dipole surface as well.

        :param transformation: the affine transformation to apply, either a raw transformation matrix or an object exposing an `apply` method (e.g. a `MolecularTransformation`)
        :type transformation: np.ndarray | object
        :param load_properties: whether to force-load the normal modes/potential/dipole derivatives before transforming (so they get carried over even if not already computed); `None` to load only if cheaply available
        :type load_properties: bool | None
        :param embed_properties: whether to propagate the transformation onto already-loaded (or newly loaded) normal modes/potential/dipole surfaces
        :type embed_properties: bool
        :return: the transformed molecule (with `source_file` cleared)
        :rtype: Molecule
        """
        from .Transformations import MolecularTransformation

        if not hasattr(transformation, 'apply'): # To support affine transformations
            transformation = MolecularTransformation(transformation)
        if load_properties:
            pot_derivs = self.potential_derivatives
            dip_derivs = self.dipole_derivatives
            modes = self.normal_modes.modes.basis
        new = transformation.apply(self)
        if embed_properties:
            if load_properties is None:
                _ = self.normal_modes.get_normal_modes(quiet=True, compute_force_constants=False)
                if _ is None:
                    load_properties = False

            if load_properties or new._normal_modes._modes is not None:
                new.normal_modes = new.normal_modes.apply_transformation(transformation)

            if load_properties is None:
                try:
                    _ = self.potential_derivatives
                except ValueError:
                    load_properties = False

            if (
                    load_properties
                    or new._pes._surf is not None
                    or new._pes._derivs is not None
            ):
                new.potential_surface = new.potential_surface.apply_transformation(transformation)

            if load_properties is None:
                try:
                    _ = self.dipole_derivatives
                except ValueError:
                    load_properties = False
            if load_properties or new._dips._surf is not None or new._dips._derivs is not None:
                new.dipole_surface = new.dipole_surface.apply_transformation(transformation)
        new.source_file = None  # for safety

        return new

    def apply_rotation(self, rotation_matrix, shift_com=None, load_properties=None, embed_properties=True):
        """
        **LLM Docstring**

        Apply a rotation (optionally shifted to rotate about the center of mass) to the molecule, via `apply_affine_transformation`.

        :param rotation_matrix: the rotation matrix (or full affine matrix) to apply
        :type rotation_matrix: np.ndarray
        :param shift_com: whether to first shift so the rotation is performed about the center of mass (only applied if `rotation_matrix` isn't already a 4x4 affine matrix)
        :type shift_com: bool | None
        :param load_properties: forwarded to `apply_affine_transformation`
        :type load_properties: bool | None
        :param embed_properties: forwarded to `apply_affine_transformation`
        :type embed_properties: bool
        :return: the rotated molecule
        :rtype: Molecule
        """

        if shift_com and rotation_matrix.shape[-1] != 4:
            com = self.center_of_mass
            rotation_matrix = nput.affine_matrix(rotation_matrix, -com)

        return self.apply_affine_transformation(rotation_matrix,
                                                load_properties=load_properties,
                                                embed_properties=embed_properties)

    def eckart_frame(self,
                     mol,
                     sel=None, inverse=False, planar_ref_tolerance=None,
                     proper_rotation=False,
                     reset_com=True
                     ):
        """
        Gets the Eckart frame(s) for the molecule

        :param mol:
        :type mol:
        :param sel: selection of atoms to use when getting the Eckart frame
        :type sel:
        :param inverse: whether to return the inverse of the rotations or not
        :type inverse: bool
        :return:
        :rtype: MolecularTransformation
        """
        return MolecularProperties.eckart_transformation(
            self,
            mol,
            sel=sel, inverse=inverse,
            planar_ref_tolerance=planar_ref_tolerance,
            proper_rotation=proper_rotation,
            reset_com=reset_com
        )

    def embed_coords(self, crds, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False):
        """
        Embeds coords in the Eckart frame using `self` as a reference

        :param crds:
        :type crds:
        :return:
        :rtype:
        """

        return self.prop('eckart_embedded_coords', crds, sel=sel, in_paf=in_paf,
                         planar_ref_tolerance=planar_ref_tolerance,
                         proper_rotation=proper_rotation
                         )
    def get_embedding_data(self, crds, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False):
        """
        Gets the necessary data to embed crds in the Eckart frame using `self` as a reference

        :param crds:
        :type crds:
        :return:
        :rtype:
        """
        return self.prop('eckart_embedding_data', crds, sel=sel, in_paf=in_paf,
                         planar_ref_tolerance=planar_ref_tolerance,
                         proper_rotation=proper_rotation)
    def get_embedded_molecule(self,
                              ref=None,
                              sel=None,
                              planar_ref_tolerance=None,
                              proper_rotation=False,
                              embed_properties=True,
                              load_properties=None,
                              reset_com=True
                              ):
        """
        Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
        a principle-axis embedded Molecule
        :return:
        :rtype: Molecule
        """

        if ref is None:
            frame = self.principle_axis_frame(sel=sel, inverse=False)
        else:
            frame = self.eckart_frame(ref, sel=sel, planar_ref_tolerance=planar_ref_tolerance, inverse=False,
                                      proper_rotation=proper_rotation, reset_com=reset_com
                                      )
        return self.apply_rotation(frame, load_properties=load_properties, embed_properties=embed_properties)

    def get_rmsd(self, other:'typing.Self | np.ndarray', sel=None,
                 embed=True,
                 embedding_sel=None,
                 mass_weighted=False):
        """
        **LLM Docstring**

        Compute the (optionally mass-weighted) root-mean-square displacement between this molecule's geometry and another geometry/molecule, optionally Eckart-embedding both onto a common reference frame first and restricting the comparison to a subset of atoms.

        :param other: the geometry (or molecule) to compare against
        :type other: Molecule | np.ndarray
        :param sel: a selection of atoms to restrict the RMSD computation to
        :type sel: Iterable[int] | None
        :param embed: whether to Eckart-embed both geometries onto a common frame before comparing
        :type embed: bool
        :param embedding_sel: a selection of atoms to use for the embedding fit, if different from `sel`
        :type embedding_sel: Iterable[int] | None
        :param mass_weighted: whether to mass-weight the displacement before computing the norm
        :type mass_weighted: bool
        :return: the RMSD value(s)
        :rtype: float | np.ndarray
        """
        if embed:
            if embedding_sel is None: embedding_sel = sel
            self = self.get_embedded_molecule(embed_properties=False, sel=embedding_sel)
            if isinstance(other, Molecule):
                other = other.get_embedded_molecule(ref=self, sel=embedding_sel, embed_properties=False).coords
            else:
                other = self.embed_coords(other, sel=embedding_sel)
        elif isinstance(other, Molecule):
            other = other.coords

        ref = self.coords
        other = np.asanyarray(other)
        base_shape = other.shape[:-2]
        ref = np.expand_dims(ref, list(range(other.ndim-2)))
        if mass_weighted:
            mass_scaling = np.asanyarray(self.masses) / np.sum(self.masses)
            mass_scaling = np.expand_dims(mass_scaling[:, np.newaxis], list(range(other.ndim-2)))
            ref = ref * mass_scaling
            other = other * mass_scaling

        if sel is not None:
            ref = ref[..., sel, :]
            other = other[..., sel, :]

        ref = ref.reshape((-1, np.prod(ref.shape[-2:], dtype=int)))
        other = other.reshape((-1, np.prod(other.shape[-2:], dtype=int)))

        return np.linalg.norm(ref - other, axis=-1).reshape(base_shape)

    def align_molecule(self, other:'typing.Self',
                       reindex_bonds=True,
                       permute_atoms=True,
                       align_structures=True,
                       sel=None,
                       embed_properties=True,
                       load_properties=False
                       ):
        """
        Aligns `other` with `self` by first finding the reindexing of the bonds of `other` that
        lead to the best graph overlap with `self`, then determining which atoms can be permuted based on their graph
        structures, then determining which permutation of equivalent atoms leads to the best agreement between the structures,
        and then finally finding the Eckart/min-RMSD transformation after this transformation has been applied

        :param other:
        :param reindex_bonds:
        :return:
        """
        from .Transformations import MolecularTransformation

        if len(itut.dict_diff(itut.counts(other.atoms), itut.counts(self.atoms))) > 0:
            raise ValueError(f"{self} and {other} have different atoms and can't be aligned")

        if reindex_bonds:
            perm = other.edge_graph.get_reindexing(self.edge_graph)
            other = other.permute_atoms(perm)
        else:
            perm = np.arange(len(self.atoms))

        if permute_atoms:
            all_perms = permute_atoms == 'all'
            permutable_atoms = []
            bond_set = other.edge_graph.map
            rem = list(range(0, len(other.atoms)))
            for i,a in enumerate(other.atoms):
                if i not in rem: continue
                group = [i]
                for j in rem:
                    if j == i: continue
                    if a == other.atoms[j]:
                        if all_perms or (bond_set[i] - {j} == bond_set[j] - {i}): # same bonds
                            group.append(j)
                rem = np.setdiff1d(rem, group)
                if len(group) > 1:
                    permutable_atoms.append(group)
        else:
            permutable_atoms = None

        if align_structures:
            embedding_data = nput.eckart_embedding(
                self.coords,
                other.coords,
                masses=self.atomic_masses,
                sel=sel,
                permutable_groups=permutable_atoms
            )

            rot, new_coords, ref_stuff, coord_stuff = embedding_data
            ref, ref_com, ref_rot = ref_stuff
            crd, com, crd_rot = coord_stuff
            transf = MolecularTransformation(ref_com, ref_rot, rot, crd_rot.T, -com)

            other = other.apply_affine_transformation(transf,
                                                      embed_properties=embed_properties,
                                                      load_properties=load_properties)

        if permute_atoms:
            perm_2 = nput.eckart_permutation(
                self.coords,
                other.coords,
                masses=self.atomic_masses,
                sel=sel,
                permutable_groups=permutable_atoms
            )

            # perm = perm[perm_2]
            other = other.permute_atoms(perm_2)

        # elif align_structures:
        #     other = other.get_embedded_molecule(self,
        #                                         sel=sel,
        #                                         embed_properties=embed_properties,
        #                                         load_properties=load_properties)

        return other
    #endregion

    #region Input Formats
    @property
    def rdmol(self):
        """
        **LLM Docstring**

        The (cached) attached RDKit molecule representation, built lazily via `RDMolecule.from_mol` the first time it's needed (returning `None` if RDKit isn't available), and kept in sync with the current coordinates on subsequent access.

        :return: the RDKit molecule, or `None` if RDKit isn't available
        :rtype: object | None
        """
        if self._rdmol is None:
            from McUtils.ExternalPrograms import RDMolecule

            try:
                self._rdmol = RDMolecule.from_mol(self, coord_unit="BohrRadius")
            except ImportError:
                ...
        else:
            self._rdmol.coords = self.coords * UnitsData.convert("BohrRadius", "Angstroms")
        return self._rdmol

    def to_ase(self, **kwargs):
        """
        **LLM Docstring**

        Convert this molecule to an ASE (Atomic Simulation Environment) molecule object, via `ASEMolecule.from_mol`.

        :param kwargs: extra options forwarded to `ASEMolecule.from_mol`
        :type kwargs: dict
        :return: the constructed ASE molecule
        :rtype: object
        """
        from McUtils.ExternalPrograms import ASEMolecule

        return ASEMolecule.from_mol(self, coord_unit='BohrRadius', **kwargs)
    @classmethod
    def from_ase(cls, ase_mol, **kwargs):
        """
        **LLM Docstring**

        Build a `Molecule` from an ASE molecule object.

        :param ase_mol: the ASE molecule to convert
        :type ase_mol: object
        :param kwargs: extra keyword arguments merged over the ASE molecule's metadata and forwarded to the constructor
        :type kwargs: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        return cls(
            ase_mol.atoms,
            ase_mol.coords * UnitsData.convert("Angstroms", "BohrRadius"),
            **(ase_mol.meta | kwargs)
        )


    @classmethod
    def from_zmat(cls, zmat, internals=None, axes=None, origin=None, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a Z-matrix specification (either a Z-matrix string or an explicit `(atoms, (ordering, coords))` tuple), converting the internal Z-matrix coordinates to Cartesians.

        :param zmat: the Z-matrix specification, as a string or an `(atoms, (ordering, coords))` tuple
        :type zmat: str | tuple
        :param internals: the internal-coordinate specification to store on the resulting molecule; defaults to the Z-matrix's own atom ordering
        :type internals: object | None
        :param axes: reference axes to use when converting to Cartesians
        :type axes: np.ndarray | None
        :param origin: reference origin to use when converting to Cartesians
        :type origin: np.ndarray | None
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        if isinstance(zmat, str):
            (atoms, ordering, coords) = coordops.parse_zmatrix_string(zmat)
            zmat = (atoms, (ordering, coords))
        (atoms, (ordering, coords)) = zmat
        coords = CoordinateSet(coords[1:], ZMatrixCoordinates).convert(CartesianCoordinates3D,
                                                                       ordering=ordering,
                                                                       axes=axes,
                                                                       origin=origin)
        if internals is None: internals = ordering
        return cls(atoms, coords, internals=internals, **opts)
    @classmethod
    def from_openbabel(cls, mol, **opts):
        """

        :param mol:
        :type mol: pybel.mol
        :return:
        :rtype:
        """
        return cls(
            mol.atoms,
            mol.coords * UnitsData.convert("Angstroms", "BohrRadius"),
            bonds=mol.bonds,
            **opts
            # **dict(
            #     rdmol.meta,
            #     **opts
            # )
        )
    def get_obmol(self, **opts):
        """
        **LLM Docstring**

        Convert this molecule to an OpenBabel molecule object, via `OBMolecule.from_mol`.

        :param opts: extra options forwarded to `OBMolecule.from_mol`
        :type opts: dict
        :return: the constructed OpenBabel molecule
        :rtype: object
        """
        from McUtils.ExternalPrograms import OBMolecule

        return OBMolecule.from_mol(self, **opts)
    @classmethod
    def _from_log_file(cls, file, num=None, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from the final Cartesian geometry parsed out of a Gaussian `.log` file.

        :param file: path to the Gaussian log file
        :type file: str
        :param num: which structure(s) to parse from the file (forwarded to the log reader)
        :type num: int | None
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.GaussianInterface import GaussianLogReader
        with GaussianLogReader(file) as gr:
            parse = gr.parse('CartesianCoordinates', num=num)
        spec, coords = parse['CartesianCoordinates']
        ang2bohr = UnitsData.convert("Angstroms", "AtomicUnitOfLength")
        return cls(
            spec[-1][:, 1],
            ang2bohr*coords[-1],
            **opts
        )
    @classmethod
    def _from_fchk_file(cls, file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a Gaussian FChk file, using its atomic numbers and integer atomic weights to build isotope-labeled atom symbols and its real atomic weights as masses.

        :param file: path to the Gaussian FChk file
        :type file: str
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.GaussianInterface import GaussianFChkReader
        with GaussianFChkReader(file) as gr:
            parse = gr.parse(
                ['Coordinates', 'AtomicNumbers', 'Integer atomic weights', "Real atomic weights"]
            )
        nums = parse["AtomicNumbers"]
        wts = parse['Integer atomic weights']
        masses = parse["Real atomic weights"]

        mol = cls(
            [AtomData[a]["Symbol"] + str(b) for a, b in zip(nums, wts)],
            parse["Coordinates"],
            masses=masses,
            **opts
        )
        return mol

    @classmethod
    def _from_orca_file(cls, file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from the final Cartesian geometry (in atomic units) parsed out of an ORCA log file.

        :param file: path to the ORCA log file
        :type file: str
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import OrcaLogReader
        with OrcaLogReader(file) as gr:
            parse = gr.parse(['CartesianAUCoordinates'])['CartesianAUCoordinates']

        mol = cls(
            parse.atoms[-1],
            parse.coords[-1],
            masses=parse.masses[-1],
            **opts
        )
        return mol

    @classmethod
    def _from_molpro_file(cls, file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from the final Cartesian geometry parsed out of a MOLPRO log file.

        :param file: path to the MOLPRO log file
        :type file: str
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import MOLPROLogReader
        with MOLPROLogReader(file) as gr:
            parse = gr.parse(['CartesianCoordinates'])['CartesianCoordinates']

        mol = cls(
            parse.atoms[-1],
            parse.coords[-1],
            # masses=parse.masses[-1],
            **opts
        )
        return mol

    @classmethod
    def _from_hess_file(cls, file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from the atoms/coordinates/masses parsed out of an ORCA `.hess` file.

        :param file: path to the ORCA `.hess` file
        :type file: str
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import OrcaHessReader
        with OrcaHessReader(file) as gr:
            parse = gr.parse(["atoms"])['atoms']

        mol = cls(
            parse.atoms,
            parse.coords, # * UnitsData.convert("Angstroms", "BohrRadius"),
            masses=parse.mass,
            **opts
        )
        return mol
    @classmethod
    def from_rdmol(cls, rdmol, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from an RDKit molecule (or a raw RDKit `Mol`/owning-mol object, which is first wrapped in an `RDMolecule`), carrying over its atoms, coordinates, bonds, charges, and metadata.

        :param rdmol: the RDKit molecule (or wrappable RDKit object) to convert
        :type rdmol: object
        :param opts: extra options merged over the RDKit molecule's metadata and forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        if hasattr(rdmol, 'GetOwningMol'):
            rdmol = RDMolecule(rdmol, **opts)
        elif hasattr(rdmol, 'GetAtoms'):
            rdmol = RDMolecule.from_base_mol(rdmol, **opts)

        return cls(
            rdmol.atoms,
            rdmol.coords * UnitsData.convert("Angstroms", "BohrRadius"),
            bonds=rdmol.bonds,
            charge=rdmol.charge,
            formal_charges=rdmol.formal_charges,
            rdmol=rdmol,
            **dict(
                rdmol.meta,
                **opts
            )
        )

    @classmethod
    def _from_smiles(cls, smi,
                     add_implicit_hydrogens=True,
                     num_confs=None,
                     conf_id=None,
                     take_min=None,
                     optimize=False,
                     sanitize=False,
                     parse_name=True,
                     allow_cxsmiles=True,
                     strict_cxsmiles=True,
                     remove_hydrogens=False,
                     replacements=None,
                     parser_options=None,
                     confgen_opts=None,
                     coords=None,
                     **opts):
        """
        **LLM Docstring**

        Build one or more `Molecule` objects from a SMILES string, generating 3D conformer(s) via RDKit and converting each resulting RDKit molecule via `from_rdmol`.

        :param smi: the SMILES string to parse
        :type smi: str
        :param add_implicit_hydrogens: whether to add implicit hydrogens
        :type add_implicit_hydrogens: bool
        :param num_confs: number of conformers to generate
        :type num_confs: int | None
        :param conf_id: which conformer id to select
        :type conf_id: int | None
        :param take_min: whether to keep only the lowest-energy conformer
        :type take_min: bool | None
        :param optimize: whether to geometry-optimize the generated conformer(s)
        :type optimize: bool
        :param sanitize: whether to sanitize the parsed RDKit molecule
        :type sanitize: bool
        :param parse_name: whether to parse a name out of the SMILES string
        :type parse_name: bool
        :param allow_cxsmiles: whether to allow CXSMILES extensions
        :type allow_cxsmiles: bool
        :param strict_cxsmiles: whether to strictly enforce CXSMILES syntax
        :type strict_cxsmiles: bool
        :param remove_hydrogens: whether to remove explicit hydrogens after parsing
        :type remove_hydrogens: bool
        :param replacements: SMILES replacement/abbreviation dict
        :type replacements: dict | None
        :param parser_options: extra low-level parser options forwarded to `RDMolecule.from_smiles`
        :type parser_options: dict | None
        :param confgen_opts: extra conformer-generation options
        :type confgen_opts: dict | None
        :param coords: explicit coordinates to use instead of generating conformers
        :type coords: np.ndarray | None
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: a single molecule, or a list of molecules if multiple conformers were generated
        :rtype: Molecule | list[Molecule]
        """
        from McUtils.ExternalPrograms import RDMolecule

        if parser_options is None:
            parser_options = {}

        confs = RDMolecule.from_smiles(smi,
                                       add_implicit_hydrogens=add_implicit_hydrogens,
                                       num_confs=num_confs,
                                       conf_id=conf_id,
                                       take_min=take_min,
                                       optimize=optimize,
                                       sanitize=sanitize,
                                       parse_name=parse_name,
                                       allow_cxsmiles=allow_cxsmiles,
                                       strict_cxsmiles=strict_cxsmiles,
                                       remove_hydrogens=remove_hydrogens,
                                       replacements=replacements,
                                       confgen_opts=confgen_opts,
                                       coords=coords,
                                       **parser_options
                                       )

        if isinstance(confs, list):
            return [cls.from_rdmol(s, **opts) for s in confs]
        else:
            return cls.from_rdmol(confs, **opts)
    @classmethod
    def _from_name(cls, name, api_key=None, add_implicit_hydrogens=True, method='pubchem', **opts):
        """
        **LLM Docstring**

        Build one or more `Molecule` objects by looking up a compound name via a chemical database API (PubChem by default, or ChemSpider) to get its SMILES string(s), then parsing those via `_from_smiles`.

        :param name: the compound name to look up
        :type name: str
        :param api_key: API key for the chosen lookup service (used for ChemSpider)
        :type api_key: str | None
        :param add_implicit_hydrogens: whether to add implicit hydrogens when parsing the resulting SMILES
        :type add_implicit_hydrogens: bool
        :param method: which lookup service to use (`'pubchem'` or `'chemspider'`)
        :type method: str
        :param opts: extra options forwarded to `_from_smiles`
        :type opts: dict
        :return: a single molecule, or a list of molecules if multiple compounds matched
        :rtype: Molecule | list[Molecule]
        :raises ValueError: if the name doesn't resolve to any known compound
        """
        from McUtils.ExternalPrograms import ChemSpiderAPI, PubChemAPI

        if dev.str_is(method, 'chemspider', ignore_case=True):
            smiles = ChemSpiderAPI(api_key).get_compounds_by_name(name, fields="SMILES")
        else:
            smiles = PubChemAPI().get_compounds_by_name(name, fields="SMILES")
        if len(smiles) == 1:
            return cls._from_smiles(smiles[0]['smiles'], add_implicit_hydrogens=add_implicit_hydrogens, **opts)
        elif len(smiles) > 0:
            return [
                cls._from_smiles(s['smiles'], add_implicit_hydrogens=add_implicit_hydrogens, **opts)
                for s in smiles
            ]
        else:
            raise ValueError(f"{name} didn't resolve to any known compounds with {method}")
    @classmethod
    def _from_sdf(cls, sdf, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from an SDF string/file, via RDKit.

        :param sdf: the SDF content to parse
        :type sdf: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_sdf(sdf), **opts)
    @classmethod
    def _from_molblock(cls, sdf, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a MOL block string, via RDKit.

        :param sdf: the MOL block content to parse
        :type sdf: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_molblock(sdf), **opts)
    @classmethod
    def _from_mol2(cls, mol, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a MOL2-format string/file, via RDKit.

        :param mol: the MOL2 content to parse
        :type mol: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_mol2(mol), **opts)
    @classmethod
    def _from_pdb(cls, pdb, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a PDB-format string/file, via RDKit.

        :param pdb: the PDB content to parse
        :type pdb: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_pdb(pdb), **opts)
    @classmethod
    def _from_mrv(cls, mrv, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a ChemDraw/Marvin (MRV) format string/file, via RDKit.

        :param mrv: the MRV content to parse
        :type mrv: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_mrv(mrv), **opts)
    @classmethod
    def _from_fasta(cls, fasta, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a FASTA-format sequence string, via RDKit.

        :param fasta: the FASTA content to parse
        :type fasta: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_fasta(fasta), **opts)
    @classmethod
    def _from_helm(cls, helm, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a HELM-format string, via RDKit.

        :param helm: the HELM content to parse
        :type helm: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_helm(helm), **opts)
    @classmethod
    def _from_cdxml(cls, cdxml, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a ChemDraw XML (CDXML) string/file, via RDKit.

        :param cdxml: the CDXML content to parse
        :type cdxml: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_cdxml(cdxml), **opts)

    @classmethod
    def _from_xyz(cls, xyz, units=None, max_blocks=None, **opts):
        """
        **LLM Docstring**

        Build one or more `Molecule` objects by parsing one or more geometry blocks from XYZ-format text (or a stream/file-like object), attaching any per-block comment line as molecule metadata and converting coordinates to Bohr.

        :param xyz: the XYZ content (string, stream, or file-like object) to parse
        :type xyz: str | object
        :param units: the units the coordinates are given in; converted to Bohr (`"BohrRadius"`). If `None`, coordinates are assumed to already be in Bohr
        :type units: str | None
        :param max_blocks: maximum number of XYZ blocks (geometries) to parse; `None` parses exactly one block and returns a single `Molecule`, a negative value parses all blocks, and a positive value caps the count (returning a list)
        :type max_blocks: int | None
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: a single molecule (if `max_blocks` is `None`) or a list of molecules (one per parsed block)
        :rtype: Molecule | list[Molecule]
        """
        from McUtils.Devutils import StreamInterface
        from McUtils.Parsers import XYZParser

        single = max_blocks is None
        if single:
            max_blocks = 1
        if max_blocks < 0:
            max_blocks = None

        with StreamInterface(xyz, file_backed=True, mode='r+') as stream:
            with XYZParser(stream) as parser:
                blocks = parser.parse(max_blocks=max_blocks)

        if units is not None:
            units = UnitsData.convert(units, "BohrRadius")
        else:
            units = 1

        structs = []
        for comment, atoms, coords in blocks:
            comment = comment.strip()
            if len(comment.strip()) > 0:
                mol_opts = {'comment':comment} | opts
            else:
                mol_opts = opts
            structs.append(
                cls(
                    atoms,
                    coords * units,
                    **mol_opts
                )
            )

        if single: structs = structs[0]
        return structs

    @classmethod
    def _from_xyz_file(cls, xyz_file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` (or list of molecules) from an XYZ file, via `_from_xyz`.

        :param xyz_file: path to (or stream of) the XYZ file
        :type xyz_file: str | object
        :param opts: extra options forwarded to `_from_xyz`
        :type opts: dict
        :return: a single molecule or a list of molecules, depending on `opts['max_blocks']`
        :rtype: Molecule | list[Molecule]
        """
        return cls._from_xyz(xyz_file, **opts)
        # with open(xyz_file) as xyz:
        #     return cls._from_xyz(xyz.read(), **opts)


    @classmethod
    def _from_gspec(cls, gspec:str, charge=None, spin=None, report=None, units='Angstroms', format_options=None, **etc):
        """
        **LLM Docstring**

        Build a `Molecule` from a Gaussian "molecule specification" string (as embedded in a Gaussian job report), substituting any report-variable placeholders, normalizing formatting quirks, and dispatching to `from_string` (as a `'gspec'`-format string) or, if the string looks like a full Gaussian archive entry, to `_from_gspec_file`.

        :param gspec: the Gaussian molecule-spec string to parse
        :type gspec: str
        :param charge: the net charge; defaults to the value parsed from the spec's header if not given
        :type charge: int | None
        :param spin: the spin multiplicity; defaults to the value parsed from the spec's header if not given
        :type spin: int | None
        :param report: a parsed Gaussian job-report dict whose entries are substituted into the spec text
        :type report: dict | None
        :param units: the units the coordinates are given in
        :type units: str
        :param format_options: per-format parsing options merged with the default Z-matrix axis convention
        :type format_options: dict | None
        :param etc: extra options forwarded to `from_string`/`_from_gspec_file`
        :type etc: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        # parses from Gaussian report molecule specs
        if gspec.startswith("1\\"):
            return cls._from_gspec_file(
                io.StringIO(gspec),
                charge=charge,
                spin=spin,
                **etc
            )

        gspec = gspec.replace(",", " ")

        header, gspec = gspec.strip().split("\n", 1)
        gspec:str # idk why pycharm needs this hint...

        c, m = [int(x) for x in header.split()]
        if report is not None:
            for k,v in report.items():
                if hasattr(v, 'value'):
                    gspec = gspec.replace(
                        " " + k + " ",
                        " " + str(v.value) + " "
                    )
                    gspec = gspec.replace(
                        " " + k + "\n",
                        " " + str(v.value) + "\n"
                    )
                elif nput.is_numeric(v):
                    gspec = gspec.replace(" " + k + " ", " " + str(v) + " ")
                    gspec = gspec.replace(" " + k + "\n", " " + str(v) + "\n")
        gspec = gspec.replace(" 0\n", "\n")
        if gspec.endswith(" 0"):
            gspec = gspec[:-2]

        if format_options is None:
            format_options = {}
        format_options = dev.merge_dicts(
            {
                'zmat': {'axes': np.eye(3)[(2, 0), :]}
            },
            format_options
        )

        return cls.from_string(
            gspec,
            charge=c if charge is None else charge,
            spin=m if spin is None else spin,
            units=units,
            format_options=format_options,
            **etc
        )

    @classmethod
    def _from_gspec_file(cls,
                         logfile,
                         charge=None, spin=None,
                         potential_derivatives=None,
                         dipole_derivatives=None,
                         target_job='Freq',
                         # use_standard_orientation_coords=True,
                         **etc):
        """
        **LLM Docstring**

        Build a `Molecule` from a full Gaussian log/archive file, extracting the job report matching `target_job` (or the last report if none matches), and, for frequency jobs, the potential-energy and dipole derivative tensors (with the appropriate atomic-mass-unit conversions applied) to attach to the constructed molecule.

        :param logfile: path to (or stream of) the Gaussian log file
        :type logfile: str | object
        :param charge: the net charge to use instead of the report's own value
        :type charge: int | None
        :param spin: the spin multiplicity to use instead of the report's own value
        :type spin: int | None
        :param potential_derivatives: explicit potential derivative tensors to use instead of parsing them from the report
        :type potential_derivatives: list[np.ndarray] | None
        :param dipole_derivatives: explicit dipole derivative tensors to use instead of parsing them from the report
        :type dipole_derivatives: list[np.ndarray] | None
        :param target_job: which Gaussian job type's report to use (e.g. `'Freq'`)
        :type target_job: str
        :param etc: extra options forwarded to `from_string`
        :type etc: dict
        :return: the constructed molecule
        :rtype: Molecule
        :raises ValueError: if no job report is found in the file
        """
        # parses from Gaussian report molecule specs
        from McUtils.ExternalPrograms import (
            GaussianLogReader,
            FchkForceConstants, FchkForceDerivatives,
            FchkDipoleDerivatives, FchkDipoleHigherDerivatives
        )

        specs = ['Reports']
        # if use_standard_orientation_coords:
        #     specs.append('StandardCartesianCoordinates')
        with GaussianLogReader(logfile) as parser:
            parse = parser.parse(specs)

        reports = parse['Reports']
        if len(reports) == 0:
            raise ValueError(f"no job report found in file {logfile}")

        for report in reports:
            if report['job'] == target_job:
                break
        else:
            report = reports[-1]

        if report['job'] in {'Freq'}:
            if potential_derivatives is None:
                potential_derivatives = list(report['PotentialDeriv'])
                if len(potential_derivatives) > 1:
                    potential_derivatives[1] = FchkForceConstants(potential_derivatives[1]).array
                if len(potential_derivatives) > 2:
                    higher = FchkForceDerivatives(potential_derivatives[2])
                    amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
                    potential_derivatives = potential_derivatives[:2] + [
                        higher.third_deriv_array / np.sqrt(amu_conv),
                        higher.fourth_deriv_array.asarray() / amu_conv
                    ]

            if dipole_derivatives is None:
                dipole_derivatives = [
                    FchkDipoleDerivatives(report['DipoleDeriv']).array
                ]

        base_mol = cls.from_string(
            report['molecule'],
            'gspec',
            charge=charge,
            spin=spin,
            potential_derivatives=potential_derivatives,
            dipole_derivatives=dipole_derivatives,
            report=report,
            **etc
        )

        # if use_standard_orientation_coords:
        #     coords = parse['StandardCartesianCoordinates']
        #     if len(coords) > 0:
        #         base_mol = base_mol.modify(
        #             atoms=[a for a in base_mol.atoms if a != "X"],
        #             coords=coords[1][-1],
        #             potential_derivatives=base_mol.potential_derivatives,
        #             dipole_derivatives=base_mol.dipole_derivatives
        #         )

        return base_mol


    @classmethod
    def from_name(cls, name, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` by looking up a compound name, via `from_string` with format `'name'`.

        :param name: the compound name to look up
        :type name: str
        :param opts: extra options forwarded to `from_string`
        :type opts: dict
        :return: the constructed molecule (or list of molecules, if multiple matches are found)
        :rtype: Molecule | list[Molecule]
        """
        return cls.from_string(name, 'name', **opts)

    _atom_strs = None
    @classmethod
    def get_atom_strings(cls):
        """
        **LLM Docstring**

        The (cached) set of up-to-2-character atomic element symbols known to `AtomData`, used for heuristically detecting SMILES/Z-matrix strings.

        :return: the set of element symbol prefixes
        :rtype: set[str]
        """
        if cls._atom_strs is None:
            cls._atom_strs = {d["Symbol"][:2] for d in AtomData.data.values()}
        return cls._atom_strs
    _smi_punct=(
        'c', 'n', 'o', 's', '*', '[', ']', '(', ')', '+',
        '.', '-', '=', '#', '@', '$', ':', '/', '\\', '0','1','2','3','4','5','6','7','8','9')
    @classmethod
    def _check_smi(cls, string, atom_types, other_syms=None):
        """
        **LLM Docstring**

        Heuristically check whether a string could be a SMILES string by stripping out all known atom-type substrings and SMILES punctuation characters and checking that nothing is left over.

        :param string: the string to check
        :type string: str
        :param atom_types: the atom-symbol substrings to strip out
        :type atom_types: Iterable[str]
        :param other_syms: the punctuation/bond-symbol characters to strip out; defaults to `cls._smi_punct`
        :type other_syms: Iterable[str] | None
        :return: whether the string is fully consumed by atom types and punctuation (i.e. looks like a SMILES string)
        :rtype: bool
        """
        for s in atom_types:
            string = string.replace(s, '')
        if other_syms is None:
            other_syms = cls._smi_punct
        for s in other_syms:
            string = string.replace(s, '')
        return len(string.strip()) == 0
    @classmethod
    def _infer_str_format(cls, string:str, allow_names=False, **opts):
        """
        **LLM Docstring**

        Heuristically guess which structural-file format a raw string represents (SMILES, a compound name, MOL/SDF, Z-matrix, Gaussian molecule-spec, or XYZ), by inspecting its line count/structure and content.

        :param string: the string to classify
        :type string: str
        :param allow_names: whether a single, non-SMILES-looking word/line should be inferred as a compound name
        :type allow_names: bool
        :param opts: accepted for interface consistency but not used in this method's body
        :type opts: dict
        :return: the inferred format key (e.g. `'smi'`, `'name'`, `'mol'`, `'zmat'`, `'gspec'`, `'xyz'`)
        :rtype: str
        :raises ValueError: if no format can be inferred from the string
        """
        from McUtils.Parsers import Number, Word

        lines = string.strip().split('\n', 3)
        at_strs = cls.get_atom_strings()
        if len(lines) == 1:
            if len(string.strip().split()) == 1 and cls._check_smi(string, at_strs):
                return 'smi'
            elif allow_names:
                return 'name'
        elif 'V2000' in string or 'V3000' in string:
            return 'mol'
        elif (
            len(lines[0].split()) == 1
            and all(l.split()[0] in at_strs for l in lines[:3])
        ):
            return 'zmat'
        elif (
                len(lines[0].split()) == 2
            and lines[0].strip().isdigit() and lines[0].strip()[-1].isdigit()
        ):
            return 'gspec'
        else:
            try:
                int(lines[0].strip())
            except (TypeError, ValueError):
                pass
            else:
                return 'xyz'

            if all(
                    len(l.split()) == 4
                    and all(Number.fullmatch(t) for t in l.split()[1:])
                for l in lines[:3]
            ):
                return 'xyz'

        raise ValueError(f"can't infer molecule spec from '''{string}'''")

    @classmethod
    def get_string_format_dispatchers(cls):
        """
        **LLM Docstring**

        The mapping from string-format key to the constructor method that parses that format, used by `from_string`.

        :return: the format-to-constructor mapping
        :rtype: dict
        """
        return {
            "smi": cls._from_smiles,
            "name": cls._from_name,
            "mol": cls._from_molblock,
            "sdf": cls._from_sdf,
            "xyz": cls._from_xyz,
            "zmat": cls.from_zmat,
            "gspec": cls._from_gspec,
            "mol2": cls._from_mol2,
            "pdb": cls._from_pdb,
            "fasta": cls._from_fasta,
            "helm": cls._from_helm,
            "cdxml": cls._from_cdxml
        }
    @classmethod
    def from_string(cls, string, fmt=None, allow_names=False, format_options=None, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a raw string in any supported structural format, inferring the format automatically if not given, dispatching to the matching in-memory parser if one exists, or otherwise writing the string to a temporary file and dispatching through `from_file`.

        :param string: the structural string to parse
        :type string: str
        :param fmt: the format to parse as; inferred via `_infer_str_format` if not given
        :type fmt: str | None
        :param allow_names: whether to allow inferring the format as a compound name
        :type allow_names: bool
        :param format_options: per-format parsing options, keyed by format, merged under `opts`
        :type format_options: dict | None
        :param opts: extra options forwarded to the format-specific parser
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        if fmt is None:
            fmt = cls._infer_str_format(string, allow_names=allow_names)
        format_dispatcher = cls.get_string_format_dispatchers()
        if format_options is not None:
            opts = collections.ChainMap(opts, format_options.get(fmt, {}))

        if fmt in format_dispatcher:
            return format_dispatcher[fmt](string, **opts)
        else:
            with tempfile.NamedTemporaryFile("w+") as tf:
                tf.write(string)
                new = cls.from_file(tf.name, mode=fmt, **opts)
            new.source_file = None # don't want to keep this lying around...
        return new
    @classmethod
    def _from_ob_import(cls, file, fmt=None, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a file via OpenBabel's format-import machinery.

        :param file: path to the file to import
        :type file: str
        :param fmt: the OpenBabel format code to use; auto-detected if not given
        :type fmt: str | None
        :param opts: extra options forwarded to `from_openbabel`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        from McUtils.ExternalPrograms import OBMolecule

        return cls.from_openbabel(
            OBMolecule.from_file(file, fmt=fmt),
            **opts
        )

    @classmethod
    def get_file_format_dispatchers(cls):
        """
        **LLM Docstring**

        The mapping from file-format key (typically a file extension) to the constructor method that parses that format, used by `from_file`.

        :return: the format-to-constructor mapping
        :rtype: dict
        """
        return {
            "log": cls._from_log_file,
            "gspec": cls._from_gspec_file,
            "fchk": cls._from_fchk_file,
            "orca": cls._from_orca_file,
            "molpro": cls._from_molpro_file,
            "hess": cls._from_hess_file,
            "smi": cls._from_smiles,
            "mol": cls._from_molblock,
            "sdf": cls._from_sdf,
            "xyz": cls._from_xyz_file
        }
    @classmethod
    def from_file(cls, file, mode=None, format_options=None, use_ob_fallback=False, **opts) -> Molecule:
        """In general we'll delegate to pybel except for like Fchk and Log files

        :param file:
        :type file:
        :return:
        :rtype:
        """
        import os

        opts['source_file'] = file
        format_dispatcher = cls.get_file_format_dispatchers()

        if mode == None:
            path, ext = os.path.splitext(file)
            ext = ext.lower()
            mode = ext.strip(".")

        if format_options is not None:
            opts = collections.ChainMap(opts, format_options.get(mode, {}))

        opts['source_file'] = {'file':file, 'mode':mode}
        if mode in format_dispatcher:
            loader = format_dispatcher[mode]
            return loader(file, **opts)
        elif use_ob_fallback:
            return cls._from_ob_import(file, fmt=mode, **opts)
        else:
            raise ValueError(f"couldn't parse file with format {mode}")

            # try:
            #     pybel = OpenBabelInterface().pybel
            # except ImportError:
            #     pybel = None
            # if pybel is None:
            #     raise IOError("{} doesn't support file type {} without OpenBabel installed.".format(cls.__name__, mode))
            # else:
            #     mol = next(pybel.readfile(mode, file))
            #     return cls.from_openbabel(mol)

    @classmethod
    def _to_smiles(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a SMILES string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_smiles`
        :type opts: dict
        :return: the SMILES string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        rdmol = mol.rdmol
        if rdmol is not None:
            return rdmol.to_smiles(**opts)
        else:
            raise ValueError(f"couldn't get `rdmol` for {mol}")

    @classmethod
    def _to_molblock(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a MOL block string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_molblock`
        :type opts: dict
        :return: the MOL block string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        rdmol = mol.rdmol
        if rdmol is not None:
            return rdmol.to_molblock(**opts)
        else:
            raise ValueError(f"couldn't get `rdmol` for {mol}")

    @classmethod
    def _to_sdf(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to an SDF string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_sdf`
        :type opts: dict
        :return: the SDF string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        rdmol = mol.rdmol
        if rdmol is not None:
            return rdmol.to_sdf(**opts)
        else:
            raise ValueError(f"couldn't get `rdmol` for {mol}")

    @classmethod
    def _to_pdb(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a PDB string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_pdb`
        :type opts: dict
        :return: the PDB string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        rdmol = mol.rdmol
        if rdmol is not None:
            return rdmol.to_pdb(**opts)
        else:
            raise ValueError(f"couldn't get `rdmol` for {mol}")

    @classmethod
    def _to_cml(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a CML string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_cml`
        :type opts: dict
        :return: the CML string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        rdmol = mol.rdmol
        if rdmol is not None:
            return rdmol.to_cml(**opts)
        else:
            raise ValueError(f"couldn't get `rdmol` for {mol}")

    @classmethod
    def _to_mrv(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a ChemDraw/Marvin (MRV) string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_mrv`
        :type opts: dict
        :return: the MRV string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        rdmol = mol.rdmol
        if rdmol is not None:
            return rdmol.to_mrv(**opts)
        else:
            raise ValueError(f"couldn't get `rdmol` for {mol}")

    @classmethod
    def _to_xyz_string(cls, mol, comment=None, units=None, num_prec=8):
        """
        **LLM Docstring**

        Format a molecule's geometry as an XYZ-format string, using a comment line (defaulting to the molecule's `repr`) and computing a column width wide enough for the largest coordinate value.

        :param mol: the molecule to export
        :type mol: Molecule
        :param comment: the comment line to use; defaults to `repr(mol)`
        :type comment: str | None
        :param units: the units to convert the coordinates into before formatting; if `None`, coordinates are left in Bohr
        :type units: str | None
        :param num_prec: number of digits after the decimal point for each coordinate
        :type num_prec: int
        :return: the formatted XYZ string
        :rtype: str
        """
        ats = mol.atoms
        crds = mol.coords
        if units is not None:
            crds = crds * UnitsData.convert("BohrRadius", units)
        num_ats = len(ats)
        x_width = 1 + np.ceil(np.max(np.log10(np.abs(crds.flatten()) + 1e-6)))
        total_width = int(2 + x_width + num_prec)
        return "\n".join([
            str(num_ats),
            repr(mol) if comment is None else comment
        ] + [
            f"{at:<3} {c[0]:>{total_width}.{num_prec}f} {c[1]:>{total_width}.{num_prec}f} {c[2]:>{total_width}.{num_prec}f}"
            for at, c in zip(ats, crds)
        ])

    @classmethod
    def _to_zmat_string(cls, mol, units='Angstroms', float_fmt="{:8.4f}",
                        variables=None, variable_modifications=None, **etc):
        """
        **LLM Docstring**

        Format a molecule's internal (Z-matrix) coordinates as a Z-matrix string, via `McUtils.Coordinerds.format_zmatrix_string`.

        :param mol: the molecule to export
        :type mol: Molecule
        :param units: the units to format distances in
        :type units: str
        :param float_fmt: the format string used for numeric values
        :type float_fmt: str
        :param variables: named variables to substitute for coordinate values
        :type variables: dict | None
        :param variable_modifications: modifications to apply to specific variable values
        :type variable_modifications: dict | None
        :param etc: extra options forwarded to `format_zmatrix_string`
        :type etc: dict
        :return: the formatted Z-matrix string
        :rtype: str
        :raises ValueError: if the molecule has no internal coordinates, or if they aren't a Z-matrix system
        """
        from McUtils.Coordinerds import format_zmatrix_string
        ics = mol.internal_coordinates
        if ics is None:
            raise ValueError("can't write ZMatrix without internal coordinates")
        if 'ZMatrix' not in mol.internal_coordinates.system.name:
            raise ValueError(f"{mol.internal_coordinates.system} isn't a Z-matrix coordinate system")

        return format_zmatrix_string(
            mol.atoms,
            mol.internal_coordinates,
            ordering=mol.internals.get('zmatrix'),
            units=units,
            float_fmt=float_fmt,
            variables=variables,
            variable_modifications=variable_modifications,
            **etc
        )

    @classmethod
    def get_string_export_dispatchers(cls):
        """
        **LLM Docstring**

        The mapping from string-export format key to the exporter method that produces that format, used by `to_string`.

        :return: the format-to-exporter mapping
        :rtype: dict
        """
        return {
            "smi": cls._to_smiles,
            "mol": cls._to_molblock,
            "sdf": cls._to_sdf,
            "pdb": cls._to_pdb,
            "cml": cls._to_cml,
            "mrv": cls._to_mrv,
            "xyz": cls._to_xyz_string,
            "zmat": cls._to_zmat_string,
        }
    def to_string(self, fmt, **opts):
        """
        **LLM Docstring**

        Export this molecule to a string in the given format, dispatching to an in-memory string exporter if one exists, otherwise round-tripping through a temporary file via `to_file`, or falling back to OpenBabel's string-export machinery.

        :param fmt: the export format key
        :type fmt: str
        :param opts: extra options forwarded to the format-specific exporter
        :type opts: dict
        :return: the exported string
        :rtype: str
        """
        format_dispatcher = self.get_string_export_dispatchers()
        file_format_dispatcher = self.get_file_export_dispatchers()

        if fmt in format_dispatcher:
            exporter = format_dispatcher[fmt]
            return exporter(self, **opts)
        elif fmt in file_format_dispatcher:
            import tempfile as tf
            with tf.NamedTemporaryFile(mode='w+') as file:
                name = file.name
            name = self.to_file(name, fmt, **opts)
            try:
                with open(name) as file:
                    return file.read()
            finally:
                try:
                    os.remove(name)
                except:
                    ...
        else:
            obmol = self.get_obmol()
            return obmol.to_string(fmt)



    @classmethod
    def get_file_export_dispatchers(cls):
        """
        **LLM Docstring**

        The mapping from file-export format key to the exporter method that writes that format to disk, used by `to_file`. Currently empty (no dedicated file exporters beyond the string-based ones and OpenBabel's fallback).

        :return: the (currently empty) format-to-exporter mapping
        :rtype: dict
        """
        return {
        }
    def to_file(self, file, mode=None, use_ob_fallback=False, **opts):
        """
        :param file:
        :type file:
        :return:
        :rtype:
        """
        import os

        format_dispatcher = self.get_file_export_dispatchers()
        string_format_dispatcher = self.get_string_export_dispatchers()

        if mode == None:
            path, ext = os.path.splitext(file)
            ext = ext.lower()
            mode = ext.strip(".")

        if mode in format_dispatcher:
            exporter = format_dispatcher[mode]
            return exporter(self, file, **opts)
        elif mode in string_format_dispatcher:
            exporter = string_format_dispatcher[mode]
            data = exporter(self, **opts)
            return dev.write_file(file, data)
        elif use_ob_fallback:
            obmol = self.get_obmol()
            return obmol.to_file(file, mode)
        else:
            raise ValueError(f"couldn't convert to file with format {mode}")

    @classmethod
    def _infer_spec_format(cls, spec, **opts):
        """
        **LLM Docstring**

        Heuristically classify a raw molecule-construction `spec` (an RDKit-like object, a file path, a raw string, a dict of constructor kwargs, or an `(atoms, coords[, opts])` tuple) so `construct` can dispatch it to the right constructor.

        :param spec: the specification to classify
        :type spec: object
        :param opts: accepted for interface consistency but not used in this method's body
        :type opts: dict
        :return: `(fmt, subopts)` where `fmt` is one of `'rdmol'`, `'file'`, `'str'`, `'dict'`, or an `(atoms, coords)` pair, and `subopts` is any extra options bundled with the spec
        :rtype: tuple
        :raises ValueError: if `spec` doesn't match any recognized shape
        """
        if all(hasattr(spec, k) for k in ['atoms', 'coords', 'bonds', 'meta']):
            return 'rdmol', {}
        elif isinstance(spec, str):
            if os.path.isfile(spec) or any(spec.endswith('.'+fmt) for fmt in cls.get_file_format_dispatchers()):
                return 'file', {}
            else:
                return 'str', {}
        elif dev.is_dict_like(spec):
            return 'dict'

        atoms = coords = None
        opts = {}
        if len(spec) == 2:
            atoms, coords = spec
        elif len(spec) == 3:
            atoms, coords, opts = spec

        if atoms is None:
            raise ValueError(f"don't know how to build a molecule from {spec}")

        return (atoms, coords), opts

    @classmethod
    def construct(cls, spec, fmt=None, **opts):
        """
        **LLM Docstring**

        Universal `Molecule` constructor: builds a molecule from essentially any reasonable input (an existing `Molecule` to copy/modify, an RDKit/ASE object, a file path, a raw structural string, a dict of constructor kwargs, or an `(atoms, coords)`/Z-matrix pair), inferring the format automatically if not given.

        :param spec: the specification to build the molecule from
        :type spec: object
        :param fmt: an explicit format to use instead of inferring one; can also be an `(atoms, coords)` pair for direct construction
        :type fmt: str | tuple | None
        :param opts: extra options forwarded to the resolved constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        :raises NotImplementedError: if `spec`/`fmt` don't match any supported construction path
        """
        if isinstance(spec, Molecule):
            return spec.modify(**opts)

        if fmt is None:
            fmt, subopts = cls._infer_spec_format(spec, **opts)
        else:
            subopts = {}
        opts = collections.ChainMap(opts, subopts)
        if fmt == 'rdmol':
            return cls.from_rdmol(spec, **opts)
        elif fmt == 'ase':
            return cls.from_ase(spec, **opts)
        elif fmt == 'file':
            return cls.from_file(spec, **opts)
        elif fmt == 'str':
            return cls.from_string(spec, **opts)
        elif fmt == 'dict':
            return cls(**dict(spec, **opts))
        elif isinstance(fmt, str):
            if isinstance(spec, str):
                if os.path.isfile(spec):
                    return cls.from_file(spec, fmt, **opts)
                else:
                    return cls.from_string(spec, fmt, **opts)
            else:
                raise NotImplementedError(f"constructing a mol from {spec} with fmt {fmt} not supported")
        else:
            atoms, coords = fmt
            if isinstance(coords, tuple) and len(coords) == 2:
                return cls.from_zmat(fmt, **dict(subopts, **opts))
            else:
                return cls(atoms, coords, **dict(subopts, **opts))
    #endregion

    #region Visualization
    def plot(self,
             *geometries,
             figure=None,
             return_objects=False,
             bonds=None,
             bond_radius=None,
             atom_radius_scaling=None,
             atom_style=None,
             atom_radii=None,
             atom_text=None,
             display_atom_numbers=False,
             radius_type=None,
             bond_style=None,
             reconcile_bonds=True,
             capped_bonds=None,
             reflectiveness=None,
             vector_style=None,
             highlight_atoms=None,
             highlight_bonds=None,
             highlight_rings=None,
             highlight_styles=None,
             comparison_styles=None,
             animation_frame_styles=None,
             mode_vectors=None,
             mode_vector_origins=None,
             mode_vector_origin_mode='set',
             mode_vector_display_cutoff=1e-2,
             principle_axes=None,
             principle_axes_origin=None,
             principle_axes_origin_mode='set',
             principle_axes_style=None,
             dipole=None,
             dipole_origin=None,
             dipole_origin_mode='set',
             render_multiple_bonds=None,
             render_fractional_bonds=None,
             fractional_bond_offset=None,
             bond_center_radius_offset=None,
             draw_coords=None,
             draw_coords_style=None,
             up_vector=None,
             multiple_bond_spacing=None,
             mode=None,#'quality',
             backend=None,
             include_save_buttons=None,
             objects=False,
             graphics_class=None,
             cylinder_class=None,
             cylinder_options=None,
             sphere_class=None,
             sphere_options=None,
             arrow_class=None,
             arrow_options=None,
             line_class=None,
             line_options=None,
             disk_class=None,
             disk_options=None,
             animate=None,
             recording_options=None,
             animation_options=None,
             jsmol_load_script=None,
             include_jsmol_script_interface=None,
             dynamic_loading=None,
             units="Angstroms",
             label_style=None,
             theme='default',
             theme_function=None,
             plot_range_padding='auto',
             annotation_function=None,
             **plot_ops
             ):
        """
        Dispatches to the appropriate `MoleculePlotter` for the resolved backend/mode.

        The full keyword surface now lives on `MoleculePlotter.plot_molecule`; calls,
        defaults, and return conventions are unchanged.
        """
        from .Visualizations import MoleculePlotter
        return MoleculePlotter.plot_molecule(self, *geometries,
                                             figure=figure,
                                             mode=mode,
                                             backend=backend,
                                             theme=theme,
                                             return_objects=return_objects,
                                             units=units,
                                             bonds=bonds,
                                             reconcile_bonds=reconcile_bonds,
                                             bond_radius=bond_radius,
                                             atom_radius_scaling=atom_radius_scaling,
                                             atom_style=atom_style,
                                             atom_radii=atom_radii,
                                             atom_text=atom_text,
                                             display_atom_numbers=display_atom_numbers,
                                             radius_type=radius_type,
                                             bond_style=bond_style,
                                             capped_bonds=capped_bonds,
                                             reflectiveness=reflectiveness,
                                             vector_style=vector_style,
                                             highlight_atoms=highlight_atoms,
                                             highlight_bonds=highlight_bonds,
                                             highlight_rings=highlight_rings,
                                             highlight_styles=highlight_styles,
                                             comparison_styles=comparison_styles,
                                             animation_frame_styles=animation_frame_styles,
                                             mode_vectors=mode_vectors,
                                             mode_vector_origins=mode_vector_origins,
                                             mode_vector_origin_mode=mode_vector_origin_mode,
                                             mode_vector_display_cutoff=mode_vector_display_cutoff,
                                             principle_axes=principle_axes,
                                             principle_axes_origin=principle_axes_origin,
                                             principle_axes_origin_mode=principle_axes_origin_mode,
                                             principle_axes_style=principle_axes_style,
                                             dipole=dipole,
                                             dipole_origin=dipole_origin,
                                             dipole_origin_mode=dipole_origin_mode,
                                             render_multiple_bonds=render_multiple_bonds,
                                             render_fractional_bonds=render_fractional_bonds,
                                             fractional_bond_offset=fractional_bond_offset,
                                             bond_center_radius_offset=bond_center_radius_offset,
                                             draw_coords=draw_coords,
                                             draw_coords_style=draw_coords_style,
                                             up_vector=up_vector,
                                             multiple_bond_spacing=multiple_bond_spacing,
                                             include_save_buttons=include_save_buttons,
                                             objects=objects,
                                             graphics_class=graphics_class,
                                             cylinder_class=cylinder_class,
                                             cylinder_options=cylinder_options,
                                             sphere_class=sphere_class,
                                             sphere_options=sphere_options,
                                             arrow_class=arrow_class,
                                             arrow_options=arrow_options,
                                             line_class=line_class,
                                             line_options=line_options,
                                             disk_class=disk_class,
                                             disk_options=disk_options,
                                             animate=animate,
                                             recording_options=recording_options,
                                             animation_options=animation_options,
                                             jsmol_load_script=jsmol_load_script,
                                             include_jsmol_script_interface=include_jsmol_script_interface,
                                             dynamic_loading=dynamic_loading,
                                             label_style=label_style,
                                             theme_function=theme_function,
                                             plot_range_padding=plot_range_padding,
                                             annotation_function=annotation_function,
                                             **plot_ops)

    def get_animation_geoms(self, which, extent=.35, steps=8, strip_embedding=True, units=None,
                            coordinate_expansion=None):
        """
        **LLM Docstring**

        Build a back-and-forth looping sequence of displaced geometries for animating a single coordinate (or an arbitrary coordinate-expansion direction), via `get_scan_coordinates`.

        :param which: the coordinate index to animate (if `coordinate_expansion` is not given directly), or an explicit displacement-direction array/list of arrays (in which case this becomes the coordinate index within that expansion, defaulting to `0`)
        :type which: int | np.ndarray | list[np.ndarray]
        :param extent: the maximum displacement magnitude in each direction
        :type extent: float
        :param steps: number of steps from the equilibrium geometry out to `extent`
        :type steps: int
        :param strip_embedding: whether to strip the fixed embedding coordinates from the default Cartesians-by-internals expansion
        :type strip_embedding: bool
        :param units: units to convert the resulting geometries into; left in Bohr if `None`
        :type units: str | None
        :param coordinate_expansion: an explicit coordinate-transformation expansion to displace along, instead of the default internal-coordinate Jacobian
        :type coordinate_expansion: list[np.ndarray] | np.ndarray | None
        :return: the looping sequence of displaced geometries (out to `extent` and back)
        :rtype: np.ndarray
        """
        if nput.is_int(which):
            if coordinate_expansion is None:
                coordinate_expansion = self.get_cartesians_by_internals(2, strip_embedding=strip_embedding)
        else:
            if isinstance(which, np.ndarray):
                if which.ndim == 1:
                    which = which[np.newaxis]
                which = [which]
            coordinate_expansion = which
            which = 0
        geoms = self.get_scan_coordinates(
            [[-extent, extent, steps]],
            which=[which],
            coordinate_expansion=coordinate_expansion
        )
        geoms = np.concatenate([geoms, np.flip(geoms, axis=0)[1:]], axis=0)
        if units is not None:
            geoms = geoms * UnitsData.convert("BohrRadius", units)
        return geoms
    def animate_coordinate(self, which, extent=.5, steps=3, return_objects=False, strip_embedding=True,
                           units="Angstroms",
                           backend=None,
                           mode=None,
                           jsmol_load_script=None,
                           coordinate_expansion=None,
                           **plot_opts
                           ):
        """
        **LLM Docstring**

        Build an animation of a displaced coordinate and return it as a plottable/displayable object, either as a JSMol vibration animation or as a sequence of rendered frames via `plot`.

        :param which: the coordinate index (or explicit direction) to animate, forwarded to `get_animation_geoms`/`format_animation_file`
        :type which: int | np.ndarray | list[np.ndarray]
        :param extent: the maximum displacement magnitude
        :type extent: float
        :param steps: number of steps out to `extent`
        :type steps: int
        :param return_objects: whether to return the constructed graphics objects alongside the figure (non-JSMol path only)
        :type return_objects: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the default coordinate expansion
        :type strip_embedding: bool
        :param units: units to display the geometries in
        :type units: str
        :param backend: the rendering backend to use; defaults to `self.display_mode`
        :type backend: str | None
        :param mode: the display mode to use; defaults to `backend`
        :type mode: str | None
        :param jsmol_load_script: extra JSMol load-script text, used only in the JSMol path
        :type jsmol_load_script: str | list[str] | None
        :param coordinate_expansion: an explicit coordinate-transformation expansion to animate along
        :type coordinate_expansion: list[np.ndarray] | np.ndarray | None
        :param plot_opts: extra options forwarded to `plot`
        :type plot_opts: dict
        :return: the resulting animation figure/widget
        :rtype: object
        """
        if backend is None:
            backend = self.display_mode
        if mode is None:
            mode = backend
        if mode == 'jsmol':
            disps = self.format_animation_file(which,
                                               format="jmol",
                                               extent=extent, steps=steps, strip_embedding=strip_embedding,
                                               units="Angstroms" if units is None else units,
                                               coordinate_expansion=coordinate_expansion
                                               )
            return self.plot(xyz=disps,
                             mode='jsmol',
                             vibrate=True,
                             jsmol_load_script=jsmol_load_script,
                             **plot_opts
                             )
        else:
            geoms = self.get_animation_geoms(which, extent=extent, steps=steps, strip_embedding=strip_embedding, units=units,
                                             coordinate_expansion=coordinate_expansion)
            return self.plot(geoms, return_objects=return_objects, units=None, backend=backend, **plot_opts)

    def animate_mode(self,
                     which,
                     extent=.5,
                     steps=3,
                     modes=None,
                     coordinate_expansion=None,
                     order=None,
                     normalize=True,
                     mass_weight=False,
                     mass_scale=True,
                     frequency_scale=False,
                     **opts
                     ):
        """
        **LLM Docstring**

        Build an animation of a normal mode's displacement, converting the mode into a coordinate-expansion direction (with optional normalization, mass-weighting/scaling, and frequency scaling of the displacement extent) and delegating to `animate_coordinate`.

        :param which: the mode index to animate
        :type which: int
        :param extent: the base displacement extent, before any mass/frequency scaling
        :type extent: float
        :param steps: number of steps out to `extent`
        :type steps: int
        :param modes: the normal modes to animate; defaults to this molecule's own (converted to a fresh mode basis)
        :type modes: object | None
        :param coordinate_expansion: an additional coordinate-transformation expansion to combine with the mode-displacement direction
        :type coordinate_expansion: list[np.ndarray] | None
        :param order: the Jacobian order to use when building the mode-based coordinate expansion (only used if given)
        :type order: int | None
        :param normalize: whether to normalize the mode-displacement direction before use
        :type normalize: bool
        :param mass_weight: whether to keep the modes mass-weighted rather than removing the mass-weighting
        :type mass_weight: bool
        :param mass_scale: whether to scale the displacement extent by the mode's effective mass (only applied when not mass-weighted)
        :type mass_scale: bool
        :param frequency_scale: whether to scale the displacement extent relative to the mode's frequency
        :type frequency_scale: bool
        :param opts: extra options forwarded to `animate_coordinate`
        :type opts: dict
        :return: the resulting animation figure/widget
        :rtype: object
        """
        from ..Modes import NormalModes

        if modes is None:
            modes = self.normal_modes.modes.basis.to_new_modes()
        modes = NormalModes.prep_modes(modes)
        base_expansion = [modes.coords_by_modes]
        if normalize:
            base_expansion = [nput.vec_normalize(base_expansion[0])]
        if mass_weight:
            modes = modes.make_mass_weighted()
        else:
            modes = modes.remove_mass_weighting()
            if mass_scale:
                extent /= np.sqrt(
                    np.dot(np.repeat(self.atomic_masses, 3), base_expansion[0][which]**2)
                        * UnitsData.convert("ElectronMass", "AtomicMassUnits")
                )
        if frequency_scale:
            freqs = np.abs(modes.freqs)
            extent *= np.sqrt(np.min(freqs) / freqs[which])
        if order is not None:
            base_expansion = ModeEmbedding(self.embedding, modes).get_cartesians_by_internals(order=order)
        if coordinate_expansion is not None:
            base_expansion = nput.tensor_reexpand(coordinate_expansion, base_expansion)
        return self.animate_coordinate(which,
                                       extent=extent, steps=steps,
                                       coordinate_expansion=base_expansion,
                                       **opts
                                       )


    def _format_xyz(self, which, nat, atoms, geom, float_format='10.3f'):
        """
        **LLM Docstring**

        Format a single geometry frame as an XYZ-format block (atom count, comment line naming the frame, then one line per atom).

        :param which: the frame index/label used in the comment line
        :type which: int
        :param nat: the number of atoms
        :type nat: int
        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param geom: the Cartesian coordinates for this frame
        :type geom: np.ndarray
        :param float_format: the format spec used for each coordinate value
        :type float_format: str
        :return: the formatted XYZ block
        :rtype: str
        """
        xyz_elems = len(geom[0])
        template = "{atom} " + " ".join(f"{{xyz[{i}]:{float_format}}}" for i in range(xyz_elems))
        body = "\n".join(template.format(atom=atom, xyz=xyz) for atom, xyz in zip(atoms, geom))
        return f"{nat}\nstruct {which}\n{body}"

    def format_structs(self,
                      geoms,
                      format='xyz'
                    ):
        """
        **LLM Docstring**

        Format a batch of geometries as a concatenated multi-frame string in the given format.

        :param geoms: the geometries to format, reshaped to `(-1, natoms, 3)`
        :type geoms: np.ndarray
        :param format: the output format; currently only `'xyz'` is supported
        :type format: str
        :return: the concatenated multi-frame string
        :rtype: str
        :raises NotImplementedError: if `format` isn't `'xyz'`
        """
        geom_data = np.reshape(np.asanyarray(geoms), (-1, len(self.atoms), 3))
        atoms = self.atoms
        if format == 'xyz':
            nat = len(atoms)
            return "\n".join(
                self._format_xyz(which + 1, nat, atoms, geom)
                for which, geom in enumerate(geom_data)
            )
        else:
            raise NotImplementedError(format)
    def _format_jmol_displacements_files(self, coords, expansion, units='Angstroms'):
        """
        **LLM Docstring**

        Format a set of per-atom displacement vectors as JSMol-style XYZ-with-vibration blocks (base coordinates plus a zero placeholder column and the displacement vector per atom), one block per displacement direction supplied.

        :param coords: the base Cartesian coordinates
        :type coords: np.ndarray
        :param expansion: a length-1 list containing the displacement-direction tensor, reshaped to `(-1, natoms, 3)`
        :type expansion: list[np.ndarray]
        :param units: units to convert the coordinates/displacements into
        :type units: str
        :return: the list of formatted displacement blocks, one per direction
        :rtype: list[str]
        """
        atoms = self.atoms
        nat = len(atoms)
        disps = expansion[0].reshape(-1, nat, 3)
        if units is not None:
            coords = coords * UnitsData.convert("BohrRadius", units)
            disps = disps / UnitsData.convert("BohrRadius", units)
        # raise Exception(coords.shape, disp[0].shape)
        return [
            self._format_xyz(
                i + 1,
                nat,
                atoms,
                np.concatenate([coords, np.zeros((nat, 1)), disp], axis=1)
            )
            for i,disp in enumerate(disps)
        ]

    def format_animation_file(self, which, format='xyz', extent=.35, steps=8, strip_embedding=True, units='Angstroms',
                              coordinate_expansion=None
                              ):
        """
        **LLM Docstring**

        Build a formatted animation string/block for a displaced coordinate, either as a single JSMol vibration block (`'jmol'` format) or as a looping multi-frame structure string via `get_animation_geoms`/`format_structs`.

        :param which: the coordinate index (or explicit direction) to animate
        :type which: int | np.ndarray | list[np.ndarray]
        :param format: the output format (`'jmol'` for a JSMol vibration block, otherwise forwarded to `format_structs`)
        :type format: str
        :param extent: the maximum displacement magnitude
        :type extent: float
        :param steps: number of steps out to `extent` (non-`'jmol'` formats only)
        :type steps: int
        :param strip_embedding: whether to strip the fixed embedding coordinates from the default coordinate expansion
        :type strip_embedding: bool
        :param units: units to format the geometries in
        :type units: str
        :param coordinate_expansion: an explicit coordinate-transformation expansion to animate along
        :type coordinate_expansion: list[np.ndarray] | None
        :return: the formatted animation string
        :rtype: str
        """
        if format == 'jmol':
            if coordinate_expansion is None:
                coordinate_expansion = self.get_cartesians_by_internals(1, strip_embedding=strip_embedding)
            coordinate_expansion = coordinate_expansion[0]
            return self._format_jmol_displacements_files(
                self.coords,
                [coordinate_expansion[(which,), :] * extent],
                units
            )[0]
        else:
            geoms = self.get_animation_geoms(which, extent=extent, steps=steps, strip_embedding=strip_embedding, units=units,
                                             coordinate_expansion=coordinate_expansion)
            return self.format_structs(geoms, format)


    def to_widget(self):
        """
        **LLM Docstring**

        Build a displayable widget for this molecule, either a JSMol applet or an X3D scene, depending on `self.display_mode`.

        :return: the constructed widget/scene object
        :rtype: object
        """
        display_opts = self.display_settings
        if display_opts is None: display_opts = {}
        if self.display_mode == 'jsmol':
            obj = self.plot(mode='jsmol', return_objects=False, **display_opts)
        else:
            obj = self.plot(backend='x3d', return_objects=False, **display_opts).figure.to_x3d()
        return obj

    default_display_mode = 'jsmol'
    def _ipython_display_(self):
        """
        **LLM Docstring**

        IPython/Jupyter display hook: builds the display widget via `to_widget` and displays it.

        :return: the result of displaying the widget
        :rtype: object
        """
        return self.to_widget()._ipython_display_()
        # display_opts = self.display_settings
        # if display_opts is None: display_opts = {}
        # obj = self.plot(mode=self.display_mode, return_objects=False, **display_opts)
        # return obj._ipython_display_()
        # return self.jupyter_viz()._ipython_display_()
    #endregion

    #region External Program Properties
    def _get_ob_attr(self, item):
        """
        **LLM Docstring**

        Look up an attribute on the attached `pybel` (OpenBabel) molecule, if one is set.

        :param item: the attribute name to look up
        :type item: str
        :return: the attribute's value
        :rtype: object
        :raises AttributeError: if no pybel molecule is attached
        """
        if self._mol is None:
            raise AttributeError("No pybel molecule")
        else:
            return getattr(self._mol, item)
    #endregion

class MolecoolException(Exception):
    pass
