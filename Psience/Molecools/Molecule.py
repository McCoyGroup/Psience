"""
Provides a simple Molecule class that we can adapt as we need
Most standard functionality should be served by OpenBabel
Uses AtomData to get properties and whatnot
"""
import io
import itertools
import math
import os, numpy as np
import tempfile
import typing

from McUtils.Data import AtomData, UnitsData
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
import McUtils.Plots as plt
from McUtils.ExternalPrograms import RDMolecule
from McUtils.Scaffolding import Logger
import McUtils.Symmetry as symm

from .MoleculeInterface import *

from .CoordinateSystems import MolecularEmbedding, ModeEmbedding
from .Evaluator import (
    MolecularEvaluator, EnergyEvaluator, DipoleEvaluator, ChargeEvaluator,
    ReducedDimensionalPotentialHandler
)
from .Hamiltonian import MolecularHamiltonian
from .Properties import *
from .Serializers import MoleculePropertyCache

__all__ = [
    "Molecule",
    "MolecoolException"
]

__reload_hook__ = ["..Modes", ".MoleculeInterface", '.CoordinateSystems', '.Hamiltonian', '.Evaluator', '.Properties']

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
                 spin=None,
                 display_mode=None,
                 energy=None,
                 energy_evaluator=None,
                 dipole_evaluator=None,
                 charge_evaluator=None,
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
        self._ats = [AtomData[atom] if isinstance(atom, (int, np.integer, str)) else atom for atom in atoms]
        self._mass = np.array([a["Mass"] for a in self._ats]) if masses is None else np.asanyarray(masses)
        coords = CoordinateSet(coords, CartesianCoordinates3D)

        internals = self.canonicalize_internals(internals, self.atoms, coords, bonds, masses=self._mass)
        self.embedding = MolecularEmbedding(self.atomic_masses, coords, internals)
        self._mode_embedding = None

        self._name = name

        # properties to be returned

        self._bonds = bonds
        self.guess_bonds = guess_bonds

        self._src = None
        self.source_file = source_file

        self._rdmol = rdmol
        self._dips = DipoleSurfaceManager(self,
                                          surface=dipole_surface,
                                          derivatives=dipole_derivatives
                                          )
        self._pes = PotentialSurfaceManager(self,
                                            surface=potential_surface,
                                            derivatives=potential_derivatives
                                            )

        self._normal_modes = NormalModesManager(self, normal_modes=normal_modes)

        self.evaluator = MolecularEvaluator(self.embedding, self._normal_modes)
        self.hamiltonian = MolecularHamiltonian(self.embedding,
                                                potential_manager=self._pes,
                                                modes_manager=self._normal_modes,
                                                dipole_manager=self._dips,
                                                )

        metadata['charge'] = charge
        metadata['spin'] = spin
        self._meta = metadata

        self.display_mode = display_mode
        self._energy = energy
        self.energy_evaluator = energy_evaluator
        self.dipole_evaluator = dipole_evaluator
        self.charge_evaluator = charge_evaluator

        if checkpoint_file is True:
            checkpoint_file = self.molecule_hash() # TODO: define this based on DB work
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
               display_mode=dev.default,
               charge=dev.default,
               spin=dev.default,
               normal_modes=dev.default,
               dipole_surface=dev.default,
               potential_surface=dev.default,
               dipole_derivatives=dev.default,
               potential_derivatives=dev.default,
               meta=dev.default
               ):
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
            bonds=self._bonds if dev.is_default(bonds) else bonds,
            guess_bonds=self.guess_bonds if dev.is_default(guess_bonds) else guess_bonds,
            energy_evaluator=self.energy_evaluator if dev.is_default(energy_evaluator, allow_None=False) else energy_evaluator,
            dipole_evaluator=self.dipole_evaluator if dev.is_default(dipole_evaluator, allow_None=False) else dipole_evaluator,
            charge_evaluator=self.charge_evaluator if dev.is_default(charge_evaluator, allow_None=False) else charge_evaluator,
            display_mode=self.display_mode if dev.is_default(display_mode) else display_mode,
            charge=self.charge if dev.is_default(charge) else charge,
            spin=self.spin if dev.is_default(spin) else spin,
            internals=self.internals if dev.is_default(internals, allow_None=False) else internals,
            normal_modes=self.normal_modes if dev.is_default(normal_modes, allow_None=False) else normal_modes,
            energy=self.energy if (
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
            )
        )

    def to_state(self, serializer=None):
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
            'potential_derivatives':self.potential_derivatives,
            'dipole_derivatives':self.dipole_derivatives,
            'charge':self.charge,
            'spin':self.spin
        }
        return data
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls(**data)

    def cached_eval(self,
                    key, generator,
                    *,
                    condition=None,
                    args=(),
                    kwargs=None):
        return self.checkpoint.cached_eval(
            key,
            generator,
            condition=condition,
            args=args,
            kwargs=kwargs
        )

    @classmethod
    def _generate_auto_spec(cls, atoms, bonds, base_coords=None, **opts):
        return PrimitiveCoordinatePicker(
            atoms,
            [b[:2] for b in bonds],
            base_coords=base_coords,
            **opts
        ).coords

    @classmethod
    def _generate_stretch_spec(cls, atoms, bonds, **opts):
        return sum(coordops.get_stretch_coordinate_system([tuple(s[:2]) for s in bonds]), [])

    @classmethod
    def _auto_auto_spec(cls, spec_generator, atoms, coords, bonds, redundant=False, base_coordinates=None,
                   masses=None,
                   untransformed_coordinates=None,
                   prune_coordinates=True,
                   pruning_options=None,
                   **opts):
        base_coords = base_coordinates
        if bonds is None:
            bonds = RDMolecule.from_coords(
                                           atoms,
                                           coords * UnitsData.convert("BohrRadius", "Angstroms"),
                                           bonds,
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
        return cls._auto_auto_spec(cls._generate_auto_spec, atoms, coords, bonds, **opts)
    @classmethod
    def _stretch_spec(cls, atoms, coords, bonds, **opts):
        return cls._auto_auto_spec(cls._generate_stretch_spec, atoms, coords, bonds, **opts)

    @classmethod
    def canonicalize_internals(cls, spec, atoms, coords, bonds, relocalize=True, masses=None):
        if isinstance(spec, str) and spec.lower() == 'auto':
            spec = {
                'primitives': 'auto'
            }

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
                    spec = cls._auto_spec(atoms, coords, bonds, masses=masses, **opts)
                elif subspec.lower() == 'natural':
                    opts = spec.copy()
                    del opts['specs']
                    if 'relocalize' in opts:
                        relocalize = spec.get('relocalize', relocalize)
                        del opts['relocalize']
                    spec = cls._stretch_spec(atoms, coords, bonds, masses=masses, **opts)
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
        return spec
    def prep_internal_spec(self, spec, relocalize=True, masses=None):
        return self.canonicalize_internals(
            spec,
            self.atoms,
            self.coords,
            self.bonds,
            relocalize=relocalize,
            masses=masses
        )

    #region Base Coords
    @property
    def coords(self):
        return self.embedding.coords
    @coords.setter
    def coords(self, coords):
        self.embedding.coords = coords
    @property
    def masses(self):
        return self._mass
    @masses.setter
    def masses(self, masses):
        self._mass = masses
        self.embedding.masses = self.atomic_masses
    @property
    def internals(self):
        return self.embedding.internals
    @property
    def charge(self):
        return self._meta.get('charge', 0)
    @charge.setter
    def charge(self, c):
        self._meta['charge'] = c
    @property
    def spin(self):
        return self._meta.get('spin')
    @spin.setter
    def spin(self, c):
        self._meta['spin'] = c
    @property
    def charges(self):
        return self._meta.get('charges', None)
    @charges.setter
    def charges(self, c):
        self._meta['charges'] = c

    def get_charge_evaluator(self, evaluator=None, **opts):
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
        self.embedding = MolecularEmbedding(
            self.masses,
            self.coords,
            spec
        )
    @property
    def internal_coordinates(self):
        return self.embedding.internal_coordinates
    @property
    def redundant_internal_transformation(self):
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
        def coordinate_filter(coords):
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

    def get_bond_graph_internals(self,
                                 include_stretches=True,
                                 include_bends=True,
                                 include_dihedrals=True,
                                 pruning=None
                                 ):
        st, bo, di = coordops.get_stretch_coordinate_system(
            [tuple(b[:2]) for b in self.bonds]
        )
        bits = []
        if include_stretches:
            bits.append(st)
        if include_bends:
            bits.append(bo)
        if include_dihedrals:
            bits.append(di)
        internals = bits[0]
        for b in bits[1:]:
            internals = internals + b

        if pruning:
            if pruning is True:
                g12 = self.get_gmatrix(power=1 / 2)
                def b_gen(pos, crds):
                    return g12 @ nput.internal_coordinate_tensors(self.coords, crds, order=1)[1]
                pruning = {'method':'b_matrix', 'b_matrix':b_gen, 'max_coords':3*len(self.atoms) - 6}
            internals = coordops.prune_internal_coordinates(
                internals,
                method=pruning
            )

        return internals
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
                              coordinate_sorting=None,
                              pruning=False
                              ):
        internals = self.get_bond_graph_internals(
            include_stretches=include_stretches,
            include_bends=include_bends,
            include_dihedrals=include_dihedrals,
            pruning=pruning
        )

        labels = self.edge_graph.get_label_types()
        internals = {
            c: coordops.get_coordinate_label(
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
        if self._mode_embedding is None:
            self._mode_embedding = ModeEmbedding(self.embedding, self.normal_modes)
        return self._mode_embedding
    def get_internals(self, strip_embedding=True):
        return self.embedding.get_internals(strip_embedding=strip_embedding)

    def get_cartesians_by_internals(self, order=None, strip_embedding=False, **kw):
        return self.embedding.get_cartesians_by_internals(order=order, strip_embedding=strip_embedding, **kw)

    def get_internals_by_cartesians(self, order=None, strip_embedding=False, **kw):
        return self.embedding.get_internals_by_cartesians(order=order, strip_embedding=strip_embedding, **kw)

    def get_cartesians_by_modes(self, order=None, **kw):
        return self.mode_embedding.get_cartesians_by_internals(order=order, **kw)

    def get_modes_by_cartesians(self, order=None, **kw):
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
        if not isinstance(val, DipoleSurfaceManager):
            raise TypeError("`dipole_surface` must be {}".format(
                DipoleSurfaceManager.__name__
            ))
        self._dips = val
    @property
    def dipole_derivatives(self):
        return self.dipole_surface.derivatives
    @dipole_derivatives.setter
    def dipole_derivatives(self, derivs):
        self.dipole_surface.derivatives = derivs
    def get_cartesian_dipole_derivatives(self, order=None, evaluator=None, include_constant_term=False):
        dipole_derivatives = self.dipole_derivatives
        # if dipole_derivatives is None and self.energy_evaluator is not None:
        #     if order is None: order = 1
        #     dipole_derivatives = self.calculate_dipole(order=order)
        if dipole_derivatives is None or (order is not None and len(dipole_derivatives) < order):
            if evaluator is None: evaluator = self.dipole_evaluator
            if evaluator is not None:
                if order is None: order = 1
                dipole_derivatives = self.calculate_dipole(evaluator=evaluator, order=order)
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
        derivs = self.get_cartesian_dipole_derivatives(order=order)
        if order is None:
            order = len(derivs)
        return nput.tensor_reexpand(
            self.get_cartesians_by_internals(order, reembed=reembed, strip_embedding=strip_embedding),
            derivs,
            order
        )

    def get_hamiltonian(self,
                        embedding=None,
                        potential_derivatives=None,
                        modes=None,
                        dipole_derivatives=None,
                        **etc
                        ):
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
    def potential_surface(self):
        """
        :return:
        :rtype: PotentialSurfaceManager
        """
        return self._pes
    @potential_surface.setter
    def potential_surface(self, val):
        if not isinstance(val, PotentialSurfaceManager):
            raise TypeError("`potential_surface` must be {}".format(
                PotentialSurfaceManager.__name__
            ))
        self._pes = val
    @property
    def potential_derivatives(self):
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
        self.potential_surface.derivatives = derivs

    def get_cartesian_potential_derivatives(self, order=None, evaluator=None, use_cached=True):
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
        from ..Modes import NormalModes
        return NormalModes.from_molecule(self,
                                         masses=masses,
                                         potential_derivatives=potential_derivatives,
                                         use_internals=use_internals,
                                         project_transrot=project_transrot,
                                         **opts
                                         )
    def get_reaction_path_modes(self, masses=None, potential_derivatives=None, **opts):
        from ..Modes import ReactionPathModes
        return ReactionPathModes.from_molecule(self, masses=masses, potential_derivatives=potential_derivatives, **opts)
    @property
    def metadata(self):
        return self._meta
    @metadata.setter
    def metadata(self, val):
        if not isinstance(val, dict):
            raise TypeError("metadata must be {}".format(
                dict.__name__
            ))
        self._meta = val


    def get_harmonic_spectrum(self, **opts):
        from ..Spectra import HarmonicSpectrum
        return HarmonicSpectrum.from_mol(self, **opts)
    #endregion

    def __repr__(self):
        name = self._name
        if name is None:
            src = self.source_file
            if src is not None:
                name = os.path.basename(src)
        if name is None:
            rdmol = self.rdmol
            if rdmol is not None:
                name = rdmol.to_smiles()
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
        return len(self._ats)
    @property
    def atom_positions(self):
        """
        A mapping of atom types to positions

        :param spec:
        :type spec:
        :return:
        :rtype:
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
        ats = self.atom_positions
        return ats['X'] if 'X' in ats else []
    @property
    def atoms(self):
        return tuple(a["Symbol"] for a in self._ats)
    # @property
    # def masses(self):
    #     if self._mass is None:
    #         return np.array([a["Mass"] for a in self._ats])
    #     else:
    #         return self._mass
    def _atomic_masses(self):
        m = self.masses
        if np.min(m) < 100:
            m = m*UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        return m
    @property
    def atomic_masses(self):
        return self._atomic_masses()
    @property
    def bonds(self):
        if self._bonds is None and self.guess_bonds:
            self._bonds = self.get_guessed_bonds()
        return self._bonds
    @bonds.setter
    def bonds(self, b):
        self._bonds = b
    @property
    def formula(self):
        return self.prop('chemical_formula')
    @property
    def multiconfig(self):
        return self.coords.multiconfig
    @property
    def name(self):
        if self._name is None:
            return "Unnamed"
        else:
            return self._name
    @property
    def source_file(self):
        if self._src is not None:
            return self._src['file']
    @source_file.setter
    def source_file(self, src):
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
        if self._src is not None:
            return self._src.get('mode')

    @property
    def shape(self):
        return self.coords.shape
    def __len__(self):
        if self.multiconfig:
            return self.coords.shape[0]
        else:
            return 1

    def copy(self):
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
        if mode is None:
            mode = self.bond_guessing_mode
        if mode == 'rdkit':
            return RDMolecule.from_coords(
                self.atoms,
                self.coords * UnitsData.convert("BohrRadius", "Angstroms"),
                None,
                guess_bonds=True,
                **opts
            ).bonds
        else:
            return MolecularProperties.guessed_bonds(self, **opts)
        # self._bonds = self.prop("guessed_bonds", tol=1.05, guess_type=True)

    @property
    def edge_graph(self):
        return MolecularProperties.edge_graph(self)

    def find_heavy_atom_backbone(self, root=None):
        return self.edge_graph.find_longest_chain(root=root)

    def find_backbone_segments(self, root=None):
        return self.edge_graph.segment_by_chains(root=root)

    def get_backbone_zmatrix(self, root=None, segments=None, return_remainder=False, return_segments=False):
        if segments is None:
            segments = self.find_backbone_segments(root=root)

        bond_list = [b[:2] for b in self.bonds]
        base_graph = coordops.bond_graph_zmatrix(
            bond_list,
            segments
        )
        zmat, new_bonds = coordops.add_missing_zmatrix_bonds(
            base_graph,
            bond_list
        )

        if return_segments or return_remainder:
            res = (zmat,)
            if return_segments:
                res = res + (segments,)
            if return_remainder:
                res = res + (new_bonds,)

            return res
        else:
            return zmat

    def get_bond_zmatrix(self, fragments=None, segments=None, root=None,
                         attachment_points=None,
                         check_attachment_points=True
                         ):
        no_frag = fragments is None
        if no_frag:
            fragments = self.fragment_indices

        if len(fragments) == 1:
            if segments is not None and len(segments) == 1:
                segments = segments[0]
            return self.get_backbone_zmatrix(root=root, segments=segments)
        else:
            inds = fragments
            if no_frag:
                ordering = np.argsort([-len(x) for x in inds])
                inds = [inds[i] for i in ordering]
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

            root = list(root) + [None] * (len(inds) - len(root))
            zmats = [
                f.get_backbone_zmatrix(root=r)
                for r,f in zip(root, frags)
            ]

            # inds = [inds[i] for i in ordering]
            # zmats = [zmats[i] for i in ordering]

            dm = nput.distance_matrix(self.coords)
            h_pos = [i for i,a in enumerate(self.atoms) if a in {'H', 'D'}]
            dm[:, h_pos] = 1e8
            dm[h_pos, :] = 1e8

            return coordops.complex_zmatrix(
                [b[:2] for b in self.bonds],
                inds,
                zmats,
                distance_matrix=dm,
                attachment_points=attachment_points,
                check_attachment_points=check_attachment_points
            )

    @property
    def fragment_indices(self):
        return MolecularProperties.fragment_indices(self)

    @property
    def fragments(self):
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
        if self._energy is None:
            self._energy = self._load_energy()
            if self._energy is None and self.energy_evaluator is not None:
                self._energy = self.calculate_energy()
        return self._energy

    default_energy_evalutor = 'rdkit'
    def get_energy_evaluator(self, evaluator=None, **opts):
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
        evaluator = self.get_energy_evaluator(evaluator, **opts)
        def evaluate_energy(coords, order=order):
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
        opt_coords = opt_coords / conv
        if reembed:
            opt_coords = nput.eckart_embedding(
                self.coords, opt_coords, masses=self.masses
            ).coordinates
        return self.modify(coords=opt_coords, meta=settings)

    def get_dipole_evaluator(self, evaluator=None, **opts):
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
        evaluator = self.get_dipole_evaluator(evaluator, **opts)
        def evaluate_dipole(coords, order=order):
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
        evaluator = self.get_dipole_evaluator(evaluator, **opts)
        smol = order is None
        if smol: order = 0
        expansion = evaluator.evaluate(
            self.coords * UnitsData.convert("BohrRadius", evaluator.distance_units),
            order=order
        )
        if smol: expansion = expansion[0]
        return expansion

    def get_reduced_potential_generator(self):
        return ReducedDimensionalPotentialHandler(self)
    def get_1d_potentials(self,
                          spec,
                          evaluator=None,
                          energy_expansion=None,
                          potential_params=None,
                          **opts):
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
        return self.evaluator.get_nearest_displacement_coordinates(
            points,
            sel=sel, axes=axes, weighting_function=weighting_function,
            modes_nearest=modes_nearest,
            return_distances=return_distances
        )
    def get_nearest_scan_coordinates(self, domains, sel=None, axes=None):
        return self.evaluator.get_nearest_scan_coordinates(domains, sel=sel, axes=axes)

    @classmethod
    def _get_atomic_radius(cls, atom_data, radius_type=None):
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
                return _baseclass(
                    *_subgrids,
                    subvals.reshape(-1),
                    **opts
                )
            return plt.TensorPlot(values, plot_class=plotter, epilog=epilog, **plot_options)
        else:
            return plotter(*np.moveaxis(grid_points, -1, 0), values, epilog=epilog, **plot_options)

    def get_model(self, potential_specs, dipole=None):
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
                if isinstance(idx, int):
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
                    if isinstance(idx, int):
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
    def get_point_group(self, *, sel=None, verbose=False, return_identifier=False, **tols):
        coords = self.coords
        masses = self.atomic_masses
        if sel is not None:
            coords = coords[sel,]
            masses = masses[sel,]
        pg_id = symm.PointGroupIdentifier(coords, masses, verbose=verbose, **tols)
        _, pg = pg_id.identify_point_group()
        if return_identifier:
            return pg_id, pg
        # print(pg)
        return pg

    def get_point_group_embedded_coordinates(self,
                                             pg=None, sel=None,
                                             return_point_group=False,
                                             return_identifier=False,
                                             **tols):
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
                    samples=50,
                    radius_scaling=1,
                    **etc):
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
                   track_velocities=False
                   ):
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
                subdirs = np.asanyarray(initial_mode_directions) * freqs[sampled_modes,][np.newaxis]
                initial_energies[:, sampled_modes] = subdirs

            if initial_energies is None:
                freqs = new.normal_modes.modes.freqs

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
                track_velocities=track_velocities
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
                  **opts
                  ):
        from ..VPT2 import VPTRunner, AnalyticVPTRunner
        if dev.str_is(runner, 'matrix'):
            runner = VPTRunner
        elif not hasattr(runner, 'construct'):
            runner = AnalyticVPTRunner

        og_pot_der = potential_derivatives
        if potential_derivatives is None:
            potential_derivatives = self.get_cartesian_potential_derivatives(
                evaluator=energy_evaluator,
                order=order+2
            )
        if dipole_derivatives is None:
            dipole_derivatives = self.get_cartesian_dipole_derivatives(
                evaluator=dipole_evaluator,
                order=order+1,
                include_constant_term=True
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

        if use_internals or use_internals is None:
            return runner.construct(self.modify(),
                                    states,
                                    order=order,
                                    potential_derivatives=potential_derivatives,
                                    dipole_derivatives=dipole_derivatives,
                                    modes=modes,
                                    mode_transformation=mode_transformation,
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
                **opts
            )

    def get_gmatrix(self,
                    masses=None, coords=None, use_internals=None, power=None,
                    **internals_opts
                    ):
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
                try:
                    _ = self.normal_modes.modes
                except ValueError:
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

        if shift_com and rotation_matrix.shape[-1] != 4:
            com = self.center_of_mass
            rotation_matrix = nput.affine_matrix(rotation_matrix, -com)

        return self.apply_affine_transformation(rotation_matrix,
                                                load_properties=load_properties,
                                                embed_properties=embed_properties)

    def eckart_frame(self, mol, sel=None, inverse=False, planar_ref_tolerance=None, proper_rotation=False):
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
            proper_rotation=proper_rotation
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
                              load_properties=None
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
                                      proper_rotation=proper_rotation
                                      )
        return self.apply_rotation(frame, load_properties=load_properties, embed_properties=embed_properties)

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
            # print(self.atoms)
            # print(other.atoms)
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
        from McUtils.ExternalPrograms import RDMolecule

        try:
            return RDMolecule.from_mol(self, coord_unit="BohrRadius")
        except ImportError:
            return None

    @classmethod
    def from_zmat(cls, zmat, internals=None, **opts):
        if isinstance(zmat, str):
            (atoms, ordering, coords) = coordops.parse_zmatrix_string(zmat)
            zmat = (atoms, (ordering, coords))
        (atoms, (ordering, coords)) = zmat
        coords = CoordinateSet(coords[1:], ZMatrixCoordinates).convert(CartesianCoordinates3D, ordering=ordering)
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
        from McUtils.ExternalPrograms import OBMolecule

        return OBMolecule.from_mol(self, **opts)
    @classmethod
    def _from_log_file(cls, file, num=None, **opts):
        from McUtils.GaussianInterface import GaussianLogReader
        with GaussianLogReader(file) as gr:
            parse = gr.parse('CartesianCoordinates', num=num)
        spec, coords = parse['CartesianCoordinates']
        ang2bohr = UnitsData.convert("Angstroms", "AtomicUnitOfLength")
        return cls(
            spec[:, 1],
            ang2bohr*coords[-1],
            **opts
        )
    @classmethod
    def _from_fchk_file(cls, file, **opts):
        from McUtils.GaussianInterface import GaussianFChkReader
        with GaussianFChkReader(file) as gr:
            parse = gr.parse(
                ['Coordinates', 'AtomicNumbers', 'Integer atomic weights', "Real atomic weights"]
            )
        nums = parse["AtomicNumbers"]
        wts = parse['Integer atomic weights']
        masses = parse["Real atomic weights"]

        # print(nums, wts)
        mol = cls(
            [AtomData[a]["Symbol"] + str(b) for a, b in zip(nums, wts)],
            parse["Coordinates"],
            masses=masses,
            **opts
        )
        return mol

    @classmethod
    def _from_orca_file(cls, file, **opts):
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
    def _from_hess_file(cls, file, **opts):
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
        if hasattr(rdmol, 'GetOwningMol'):
            rdmol = RDMolecule(rdmol, **opts)
        elif hasattr(rdmol, 'GetAtoms'):
            rdmol = RDMolecule.from_base_mol(rdmol, **opts)
        return cls(
            rdmol.atoms,
            rdmol.coords * UnitsData.convert("Angstroms", "BohrRadius"),
            bonds=rdmol.bonds,
            **dict(
                rdmol.meta,
                **opts
            )
        )

    @classmethod
    def _from_smiles(cls, smi,
                     add_implicit_hydrogens=True,
                     num_confs=1,
                     optimize=False,
                     sanitize=False,
                     parse_name=True,
                     allow_cxsmiles=True,
                     strict_cxsmiles=True,
                     remove_hydrogens=False,
                     replacements=None,
                     parser_options=None,
                     **opts):
        from McUtils.ExternalPrograms import RDMolecule

        if parser_options is None:
            parser_options = {}

        return cls.from_rdmol(
            RDMolecule.from_smiles(smi,
                                   add_implicit_hydrogens=add_implicit_hydrogens,
                                   num_confs=num_confs,
                                   optimize=optimize,
                                   sanitize=sanitize,
                                   parse_name=parse_name,
                                   allow_cxsmiles=allow_cxsmiles,
                                   strict_cxsmiles=strict_cxsmiles,
                                   remove_hydrogens=remove_hydrogens,
                                   replacements=replacements,
                                   **parser_options
                                   ),
            **opts
        )
    @classmethod
    def _from_name(cls, name, api_key=None, add_implicit_hydrogens=True, **opts):
        from McUtils.ExternalPrograms import ChemSpiderAPI

        smiles = ChemSpiderAPI(api_key).get_compounds_by_name(name, fields="SMILES")
        if len(smiles) == 1:
            return cls._from_smiles(smiles[0]['smiles'], add_implicit_hydrogens=add_implicit_hydrogens, **opts)
        elif len(smiles) > 0:
            return [
                cls._from_smiles(s['smiles'], add_implicit_hydrogens=add_implicit_hydrogens, **opts)
                for s in smiles
            ]
        else:
            raise ValueError(f"{name} didn't resolve to any known compounds with ChemSpider")
    @classmethod
    def _from_sdf(cls, sdf, **opts):
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_sdf(sdf), **opts)
    @classmethod
    def _from_molblock(cls, sdf, **opts):
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_molblock(sdf), **opts)

    @classmethod
    def _from_xyz(cls, xyz, units=None, **opts):
        from McUtils.Parsers import Number, Word
        if xyz[0].isdigit():
            xyz = "\n".join(xyz.splitlines()[2:]) # should confirm that regular `xyz.split("\n", 2)[-1]` would work
        atoms = Word.findall(xyz)
        coords = np.array(Number.findall(xyz)).astype(float).reshape(-1, 3)
        if units is not None:
            coords *= UnitsData.convert(units, "BohrRadius")
        return cls(atoms, coords, **opts)
        # return cls.from_rdmol(RDMolecule.from_molblock(sdf), **opts)

    @classmethod
    def _from_xyz_file(cls, xyz_file, **opts):
        with open(xyz_file) as xyz:
            return cls._from_xyz(xyz.read(), **opts)


    @classmethod
    def _from_gspec(cls, gspec:str, charge=None, spin=None, report=None, units='Angstroms', **etc):
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

        return cls.from_string(
            gspec,
            charge=c if charge is None else charge,
            spin=m if spin is None else spin,
            units=units,
            **etc
        )

    @classmethod
    def _from_gspec_file(cls,
                         logfile,
                         charge=None, spin=None,
                         potential_derivatives=None,
                         dipole_derivatives=None,
                         **etc):
        # parses from Gaussian report molecule specs
        from McUtils.ExternalPrograms import (
            GaussianLogReader,
            FchkForceConstants, FchkForceDerivatives,
            FchkDipoleDerivatives, FchkDipoleHigherDerivatives
        )

        with GaussianLogReader(logfile) as parser:
            reports = parser.parse('Reports')

        if len(reports) == 0:
            raise ValueError(f"no job report found in file {logfile}")

        reports = reports['Reports']
        if len(reports) == 0:
            raise ValueError(f"no job report found in file {logfile}")

        report = reports[-1]
        if report['job'] == 'Freq':
            if potential_derivatives is None:
                potential_derivatives = list(report['PotentialDeriv'])
                if len(potential_derivatives) > 1:
                    potential_derivatives[1] = FchkForceConstants(potential_derivatives[1]).array
                if len(potential_derivatives) > 2:
                    higher = FchkForceDerivatives(potential_derivatives[2])
                    potential_derivatives = potential_derivatives[:2] + [
                        higher.third_deriv_array,
                        higher.fourth_deriv_array
                    ]
            if dipole_derivatives is None:
                dipole_derivatives = [FchkDipoleDerivatives(report['DipoleDeriv']).array]

        return cls.from_string(
            report['molecule'],
            'gspec',
            charge=charge,
            spin=spin,
            potential_derivatives=potential_derivatives,
            dipole_derivatives=dipole_derivatives,
            report=report,
            **etc
        )

    @classmethod
    def from_name(cls, name, **opts):
        return cls.from_string(name, 'name', **opts)

    _atom_strs = None
    @classmethod
    def get_atom_strings(cls):
        if cls._atom_strs is None:
            cls._atom_strs = {d["Symbol"][:2] for d in AtomData.data.values()}
        return cls._atom_strs
    _smi_punct=(
        'c', 'n', 'o', '*', '[', ']', '(', ')', '+',
        '.', '-', '=', '#', '$' ':', '/', '\\', '0','1','2','3','4','5','6','7','8','9')
    @classmethod
    def _check_smi(cls, string, atom_types, other_syms=None):
        for s in atom_types:
            string = string.replace(s, '')
        if other_syms is None:
            other_syms = cls._smi_punct
        for s in other_syms:
            string = string.replace(s, '')
        return len(string.strip()) == 0
    @classmethod
    def _infer_str_format(cls, string:str, **opts):
        from McUtils.Parsers import Number, Word

        lines = string.strip().split('\n', 3)
        at_strs = cls.get_atom_strings()
        if len(lines) == 1:
            if len(string.split()) == 1 and cls._check_smi(string, at_strs):
                return 'smi'
            else:
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
        return {
            "smi": cls._from_smiles,
            "name": cls._from_name,
            "mol": cls._from_molblock,
            "sdf": cls._from_sdf,
            "xyz": cls._from_xyz,
            "zmat": cls.from_zmat,
            "gspec": cls._from_gspec
        }
    @classmethod
    def from_string(cls, string, fmt=None, **opts):
        if fmt is None:
            fmt = cls._infer_str_format(string)
        format_dispatcher = cls.get_string_format_dispatchers()

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
        from McUtils.ExternalPrograms import OBMolecule

        return cls.from_openbabel(
            OBMolecule.from_file(file, fmt=fmt),
            **opts
        )

    @classmethod
    def get_file_format_dispatchers(cls):
        return {
            "log": cls._from_log_file,
            "gspec": cls._from_gspec_file,
            "fchk": cls._from_fchk_file,
            "orca": cls._from_orca_file,
            "hess": cls._from_hess_file,
            "smi": cls._from_smiles,
            "mol": cls._from_molblock,
            "sdf": cls._from_sdf,
            "xyz": cls._from_xyz_file
        }
    @classmethod
    def from_file(cls, file, mode=None, **opts):
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

        opts['source_file'] = {'file':file, 'mode':mode}
        if mode in format_dispatcher:
            loader = format_dispatcher[mode]
            return loader(file, **opts)
        else:
            return cls._from_ob_import(file, fmt=mode, **opts)

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
        rdmol = mol.rdmol
        if rdmol is not None:
            return rdmol.to_smiles(**opts)
        else:
            raise ValueError(f"couldn't get `rdmol` for {mol}")

    @classmethod
    def _to_xyz_string(cls, mol, comment=None, units=None, num_prec=8):
        ats = mol.atoms
        crds = mol.coords
        if units is not None:
            crds = crds * UnitsData.convert("BohrRadius", units)
        num_ats = len(ats)
        x_width = 1 + np.ceil(np.max(np.log10(np.abs(crds.flatten()))))
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
        return {
            "smi": cls._to_smiles,
            # "mol": cls._to_molblock,
            # "sdf": cls._to_sdf,
            "xyz": cls._to_xyz_string,
            "zmat": cls._to_zmat_string,
        }
    def to_string(self, fmt, **opts):
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
        return {
        }
    def to_file(self, file, mode=None, **opts):
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
        else:
            obmol = self.get_obmol()
            return obmol.to_file(file, mode)

    @classmethod
    def _infer_spec_format(cls, spec, **opts):
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
        if isinstance(spec, Molecule):
            return spec.modify(**opts)

        if fmt is None:
            fmt, subopts = cls._infer_spec_format(spec, **opts)
        else:
            subopts = {}
        if fmt == 'rdmol':
            return cls.from_rdmol(spec, **opts)
        elif fmt == 'file':
            return cls.from_file(spec, **opts)
        elif fmt == 'str':
            return cls.from_string(spec, **opts)
        elif fmt == 'dict':
            return cls(**dict(spec, **opts))
        else:
            atoms, coords = fmt
            if isinstance(coords, tuple) and len(coords) == 2:
                return cls.from_zmat(fmt, **dict(subopts, **opts))
            else:
                return cls(atoms, coords, **dict(subopts, **opts))
    #endregion

    #region Visualization
    highlight_styles = {
        "glow":"green",
        "color":"white"
    }
    vector_style = {
        'color': 'black',
        'radius': .1
    }
    principle_axes_style = [
        {'color':'green'},
        {'color':'red'},
        {'color':'blue'}
    ]
    def plot(self,
             *geometries,
             figure=None,
             return_objects=False,
             bond_radius=.1,
             atom_radius_scaling=.25,
             atom_style=None,
             atom_radii=None,
             radius_type=None,
             bond_style=None,
             capped_bonds=False,
             reflectiveness=None,
             vector_style=None,
             highlight_atoms=None,
             highlight_bonds=None,
             highlight_rings=None,
             highlight_styles=None,
             mode_vectors=None,
             mode_vector_origins=None,
             mode_vector_origin_mode='set',
             mode_vector_display_cutoff=1e-2,
             principle_axes=None,
             principle_axes_origin=None,
             principle_axes_style=None,
             dipole=None,
             dipole_origin=None,
             dipole_origin_mode='set',
             render_multiple_bonds=True,
             up_vector=None,
             multiple_bond_spacing=None,
             mode=None,#'quality',
             backend=None,
             include_save_buttons=None,
             objects=False,
             graphics_class=None,
             cylinder_class=None,
             sphere_class=None,
             arrow_class=None,
             animate=None,
             animation_options=None,
             jsmol_load_script=None,
             units="Angstroms",
             **plot_ops
             ):

        from McUtils.Plots import Graphics3D, Sphere, Cylinder, Line, Disk, Arrow

        if backend is None:
            backend = self.display_mode
        if backend is None:
            backend = self.default_display_mode
        if mode is None:
            mode = backend

        if len(geometries) > 0 or any(
            opt is not None
            for opt in [
                highlight_atoms,
                highlight_bonds,
                highlight_rings,
                highlight_styles,
                mode_vectors,
                dipole
            ]
        ):
            if mode in {'jupyter', 'jsmol'}:
                mode = 'x3d'

        if mode == 'jupyter':
            return self.jupyter_viz()
        elif mode == 'jsmol':
            return self.jsmol_viz(script=jsmol_load_script)

        if backend in {'jupyter', 'jsmol'}:
            backend = 'x3d'

        graphics_keys = Graphics3D.known_keys | Graphics3D.opt_keys | Graphics3D.figure_keys
        graphics_opts = {k:plot_ops[k] for k in plot_ops.keys() & graphics_keys}
        plot_ops = {k:plot_ops[k] for k in plot_ops.keys() - graphics_keys}

        if graphics_class is None:
            graphics_class = Graphics3D
        if cylinder_class is None:
            cylinder_class = Line if mode == 'fast' else Cylinder
        if sphere_class is None:
            sphere_class = Disk if mode == 'fast' else Sphere
        if arrow_class is None:
            arrow_class = Arrow

        if len(geometries) == 0:
            geometries = self.coords
        elif len(geometries) == 1:
            geometries = CoordinateSet(np.asanyarray(geometries[0]), self.coords.system)
        else:
            geometries = CoordinateSet([np.asanyarray(g) for g in geometries], self.coords.system)

        if units is not None:
            geometries = geometries * UnitsData.convert("BohrRadius", units)

        if geometries.ndim == 2:
            geometries = geometries[np.newaxis]
        if animate is None:
            animate = geometries.shape[0] > 1

        if dipole is not None:
            dipole = np.asanyarray(dipole)
            if dipole.ndim == 1:
                dipole = np.broadcast_to(dipole[np.newaxis], (len(geometries),) + dipole.shape)

            if units is not None:
                dipole = dipole * UnitsData.convert("BohrRadius", units)

        if principle_axes is True:
            _, principle_axes = nput.moments_of_inertia(geometries, self.atomic_masses)
        elif principle_axes is False:
            principle_axes = None
        if principle_axes is not None:
            if isinstance(principle_axes_style, dict):
                principle_axes_style = [principle_axes_style] * 3
            principle_axes_style = [
                dict(self.principle_axes_style[i], **principle_axes_style[i])
                for i in range(3)
            ]
            if principle_axes_origin is not None:
                if principle_axes_origin.ndim == 1:
                    principle_axes_origin = np.broadcast_to(
                        principle_axes_origin[np.newaxis],
                        (len(geometries),) + principle_axes_origin.shape
                    )

                if units is not None:
                    principle_axes_origin = principle_axes_origin * UnitsData.convert("BohrRadius", units)

        if dipole_origin is not None:
            dipole_origin = np.asanyarray(dipole_origin)
            if dipole_origin.ndim == 1:
                dipole_origin = np.broadcast_to(
                    dipole_origin[np.newaxis],
                    (len(geometries),) + dipole_origin.shape
                )

            if units is not None:
                dipole_origin = dipole_origin * UnitsData.convert("BohrRadius", units)

        if mode_vectors is not None:
            mode_vectors = np.asanyarray(mode_vectors)
            if mode_vectors.ndim == 1:
                mode_vectors = np.reshape(mode_vectors, (-1, 3))
            if mode_vectors.ndim == 2:
                mode_vectors = np.broadcast_to(mode_vectors[np.newaxis], (len(geometries),) + mode_vectors.shape)

            if units is not None:
                mode_vectors = mode_vectors * UnitsData.convert("BohrRadius", units)

        if mode_vector_origins is not None:
            mode_vector_origins = np.asanyarray(mode_vector_origins)
            if mode_vector_origins.ndim == 1:
                mode_vector_origins = np.reshape(mode_vector_origins, (-1, 3))
            if mode_vector_origins.ndim == 2:
                mode_vector_origins = np.broadcast_to(mode_vector_origins[np.newaxis],
                                                      (len(geometries),) + mode_vector_origins.shape)

            if units is not None:
                mode_vector_origins = mode_vector_origins * UnitsData.convert("BohrRadius", units)

        geometries = geometries.convert(CartesianCoordinates3D)

        if figure is None:
            figure = graphics_class(backend=backend, **graphics_opts)

        colors = [ at["IconColor"] for at in self._ats ]
        if atom_radii is None:
            atom_radii = [None] * len(self._ats)
        elif dev.is_dict_like(atom_radii):
            atom_radii = [atom_radii.get(a["ElementSymbol"]) for a in self._ats]
        atom_radii = [
            self._get_atomic_radius(at, radius_type)
                if c is None else c for c, at in
            zip(atom_radii, self._ats)
        ]
        radii = [ atom_radius_scaling * r for r in atom_radii ]

        bonds = [None] * len(geometries)
        atoms = [None] * len(geometries)
        arrows = [None] * len(geometries)

        if vector_style is None:
            vector_style = {}
        vector_style = dict(self.vector_style, **vector_style)
        if atom_style is None or atom_style is True:
            atom_style = {}
        elif atom_style is False:
            ...
        elif not isinstance(atom_style, dict):
            atom_style = {i:a for i,a in enumerate(atom_style)}

        if atom_style is not False:
            base_atom_style = {}
            if reflectiveness is not None and backend == 'x3d':
                base_atom_style.update({
                    'specularity': 'white',
                    'shininess': 100 * np.clip(1.1 - reflectiveness, 0, 1)
                })
            _atom_style = {i:{} for i in range(len(self._ats))}
            for k,v in atom_style.items():
                if isinstance(k, str):
                    base_atom_style[k] = v
                else:
                    _atom_style[k] = v
            for k,v in _atom_style.items():
                _atom_style[k] = dict(base_atom_style, **v)
            atom_style = _atom_style

        if highlight_styles is None:
            highlight_styles = self.highlight_styles

        if highlight_atoms is True:
            highlight_atoms = list(range(len(self.atoms)))
        if highlight_rings is not None:
            if highlight_atoms is None:
                highlight_atoms = []
            if highlight_bonds is None:
                highlight_bonds = []
            highlight_atoms = list(highlight_atoms)
            highlight_bonds = list(highlight_bonds)
            for r in highlight_rings:
                highlight_atoms.extend(r)
                highlight_bonds.extend(zip(
                    r, r[1:] + r[:1]
                ))

        if highlight_atoms is not None:
            for k in highlight_atoms:
                if k not in atom_style: atom_style[k] = {}
                atom_style[k].update(highlight_styles)

        if bond_style is None or bond_style is True:
            bond_style = {}
        elif bond_style is False:
            ...
        elif not isinstance(bond_style, dict):
            bond_style = {i:a for i,a in enumerate(bond_style)}

        if bond_style is not False:
            base_bond_style = {}
            if capped_bonds:
                base_bond_style['capped'] = True
            if reflectiveness is not None and backend == 'x3d':
                base_bond_style.update({
                    'specularity': 'white',
                    'shininess': 100 * np.clip(1.1 - reflectiveness, 0, 1)
                })
            _bond_style = {(i,j):{} for i,j in itertools.combinations(range(len(self._ats)), 2)}
            for k,v in bond_style.items():
                if isinstance(k, str):
                    base_bond_style[k] = v
                else:
                    _bond_style[k] = v
            for k,v in _bond_style.items():
                _bond_style[k] = dict(base_bond_style, **v)
            bond_style = _bond_style
            if highlight_bonds is None and highlight_atoms is not None:
                highlight_bonds = highlight_atoms
            if highlight_bonds is not None:
                for k in highlight_bonds:
                    if not nput.is_numeric(k): k = tuple(k)
                    if k not in bond_style: bond_style[k] = {}
                    bond_style[k].update(highlight_styles)

        for i, geom in enumerate(geometries):
            bond_list = self.bonds
            if bond_style is not False and bond_list is not None:
                bonds[i] = [None] * len(bond_list)
                for j, b in enumerate(bond_list):
                    atom1 = b[0]
                    atom2 = b[1]

                    c1 = colors[atom1]
                    c2 = colors[atom2]
                    base_bstyle = dict(
                        bond_style.get((atom2, atom1), {}),
                        **bond_style.get((atom1, atom2), {})
                    )
                    b_sty_1 = dict(
                        bond_style.get(atom1, {}),
                        **base_bstyle
                    )
                    if b_sty_1.get('color') is None:
                        b_sty_1['color'] = c1
                    b_sty_2 = dict(
                        bond_style.get(atom2, {}),
                        **base_bstyle
                    )
                    if b_sty_2.get('color') is None:
                        b_sty_2['color'] = c2

                    p1 = geom[atom1]
                    p2 = geom[atom2]
                    disp_vector = p2 - p1
                    midpoint = disp_vector/2 + p1

                    if not render_multiple_bonds or len(b) == 2 or b[2] < 2 or b[2] > 3:
                        bond_point_list = [
                            [p1, p2, midpoint]
                        ]
                    elif b[2] == 2:
                        if up_vector is None:
                            up_vector = [0, 0, 1]
                        if multiple_bond_spacing is None:
                            multiple_bond_spacing = bond_radius * 1.1

                        axis = nput.vec_normalize(
                            np.cross(disp_vector, up_vector)
                        )

                        p11 = p1 - axis * multiple_bond_spacing
                        p21 = p2 - axis * multiple_bond_spacing
                        mp1 = midpoint - axis * multiple_bond_spacing


                        p12 = p1 + axis * multiple_bond_spacing
                        p22 = p2 + axis * multiple_bond_spacing
                        mp2 = midpoint + axis * multiple_bond_spacing
                        bond_point_list = [
                            [p11, p21, mp1],
                            [p12, p22, mp2]
                        ]


                    bond_objs = []
                    for pp1, pp2, mp in bond_point_list:
                        cc1 = cylinder_class(
                            pp1,
                            mp,
                            bond_radius,
                            **plot_ops,
                            **b_sty_1
                        )
                        cc2 = cylinder_class(
                            mp,
                            pp2,
                            bond_radius,
                            **plot_ops,
                            **b_sty_2
                        )

                        if objects:
                            bond_objs.extend(( cc1, cc2 ))
                        else:
                            cyl_1 = cc1.plot(figure)
                            cyl_2 = cc2.plot(figure)
                            if isinstance(cyl_1, (list, tuple)):
                                cyl_1 = cyl_1[0]
                            if isinstance(cyl_2, (list, tuple)):
                                cyl_2 = cyl_2[0]
                            bond_objs.extend(( cyl_1, cyl_2 ))

                    bonds[i][j] = bond_objs

            if atom_style is not False:
                atoms[i] = [None] * len(geom)
                for j, stuff in enumerate(zip(colors, radii, geom)):
                    color, radius, coord = stuff
                    a_sty = atom_style.get(j, {})
                    if a_sty.get('color') is None:
                        a_sty['color'] = color

                    sphere = sphere_class(coord, radius, **plot_ops, **a_sty)
                    if objects:
                        atoms[i][j] = sphere
                    else:
                        plops = sphere.plot(figure)
                        if isinstance(plops, tuple):
                            atoms[i][j] = plops[0]
                        else:
                            atoms[i][j] = plops

            if dipole is not None:
                arrows[i] = []
                dip = dipole[i]
                if np.linalg.norm(dip) > mode_vector_display_cutoff:
                    if dipole_origin is None or dipole_origin_mode == 'shift':
                        com = np.tensordot(self.masses, geom, axes=[0, 0]) / np.sum(self.masses)
                        if dipole_origin is not None:
                            com = com + dipole_origin[i]
                    else:
                        com = dipole_origin[i]
                    dipole_arrow = arrow_class(
                        com,
                        com + dip,
                        **vector_style
                    )
                    if objects:
                        arrows[i].append(dipole_arrow)
                    else:
                        plops = dipole_arrow.plot(figure)
                        if isinstance(plops, tuple):
                            arrows[i].append(plops[0])
                        else:
                            arrows[i].append(plops)

            if principle_axes is not None:
                arrows[i] = []
                pax:np.ndarray = principle_axes[i]
                if principle_axes_origin is None or principle_axes_origin == 'shift':
                    com = np.tensordot(self.masses, geom, axes=[0, 0]) / np.sum(self.masses)
                    if principle_axes_origin is not None:
                        com = com + principle_axes_origin[i]
                else:
                    com = dipole_origin[i]
                for ax, sty in zip(pax.T, principle_axes_style):
                    pax_arrow = arrow_class(
                        com,
                        com + ax,
                        **dict(vector_style, **sty)
                    )
                    if objects:
                        arrows[i].append(pax_arrow)
                    else:
                        plops = pax_arrow.plot(figure)
                        if isinstance(plops, tuple):
                            arrows[i].append(plops[0])
                        else:
                            arrows[i].append(plops)

            if mode_vectors is not None:
                if arrows[i] is None:
                    arrows[i] = []
                for j,v in enumerate(mode_vectors[i]):
                    if np.linalg.norm(v) > mode_vector_display_cutoff:
                        if mode_vector_origins is None or mode_vector_origin_mode == 'shift':
                            com = geom[j]
                            if mode_vector_origins is not None:
                                com = com + mode_vector_origins[i][j]
                        else:
                            com = mode_vector_origins[i][j]
                        mode_arrow = arrow_class(
                            com,
                            com + v,
                            **vector_style
                        )
                        if objects:
                            arrows[i].append(mode_arrow)
                        else:
                            plops = mode_arrow.plot(figure)
                            if isinstance(plops, tuple):
                                arrows[i].append(plops[0])
                            else:
                                arrows[i].append(plops)

        if animate:
            if animation_options is None: animation_options = {}
            figure = figure.animate_frames(
                [
                    (a if a is not None else [])
                    + (sum([list(b) for b in bl], []) if bl is not None else [])
                    + (ar if ar is not None else [])
                    for a,bl,ar in zip(atoms, bonds, arrows)
                ],
                **animation_options
            )

        if include_save_buttons is not None:
            plt.X3D.include_export_button = include_save_buttons
            plt.X3D.include_record_button = include_save_buttons

        if return_objects:
            return figure, atoms, bonds, arrows
        else:
            return figure

    def get_animation_geoms(self, which, extent=.35, steps=8, strip_embedding=True, units=None,
                            coordinate_expansion=None
                            ):
        if isinstance(which, int):
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
        geoms = np.concatenate([geoms, np.flip(geoms, axis=0)], axis=0)
        if units is not None:
            geoms = geoms * UnitsData.convert("BohrRadius", units)
        return geoms
    def animate_coordinate(self, which, extent=.5, steps=8, return_objects=False, strip_embedding=True,
                           units="Angstroms",
                           backend=None,
                           mode=None,
                           jsmol_load_script=None,
                           coordinate_expansion=None,
                           **plot_opts
                           ):
        if backend is None:
            backend = self.display_mode
        if backend is None:
            backend = self.default_display_mode
        if mode is None:
            mode = backend
        if mode == 'jsmol':
            disps = self.format_animation_file(which, format="jmol",
                                               extent=extent, steps=steps, strip_embedding=strip_embedding,
                                               units="Angstroms" if units is None else units,
                                               coordinate_expansion=coordinate_expansion
                                               )
            return self.jsmol_viz(disps, vibrate=True, script=jsmol_load_script)
        else:
            geoms = self.get_animation_geoms(which, extent=extent, steps=steps, strip_embedding=strip_embedding, units=units,
                                             coordinate_expansion=coordinate_expansion)
            return self.plot(geoms, return_objects=return_objects, units=None, backend=backend, **plot_opts)

    def animate_mode(self,
                     which,
                     extent=.5,
                     steps=8,
                     modes=None,
                     coordinate_expansion=None,
                     order=None,
                     normalize=True,
                     mass_weight=False,
                     mass_scale=True,
                     frequency_scale=False,
                     **opts
                     ):
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
        xyz_elems = len(geom[0])
        template = "{atom} " + " ".join(f"{{xyz[{i}]:{float_format}}}" for i in range(xyz_elems))
        body = "\n".join(template.format(atom=atom, xyz=xyz) for atom, xyz in zip(atoms, geom))
        return f"{nat}\nstruct {which}\n{body}"

    def format_structs(self,
                      geoms,
                      format='xyz'
                    ):
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

    def jsmol_viz(self, xyz=None, animate=False, vibrate=False, script=None):
        from McUtils.Jupyter import JSMol
        if xyz is None:
            xyz = self._format_xyz(0,
                                   len(self.atoms),
                                   self.atoms,
                                   self.coords * UnitsData.convert("BohrRadius", "Angstroms")
                                   )
        return JSMol.Applet(xyz, animate=animate, vibrate=vibrate, load_script=script)

    def jupyter_viz(self):
        from McUtils.Jupyter import MoleculeGraphics

        return MoleculeGraphics(self.atoms,
                                np.ndarray.view(self.coords.convert(CartesianCoordinates3D)),
                                bonds=self.bonds
                                )

    def to_widget(self):
        if self.display_mode == 'jsmol':
            obj = self.plot(mode='jsmol', return_objects=False)
        else:
            obj = self.plot(backend='x3d', return_objects=False).figure.to_x3d()
        return obj

    default_display_mode = 'jsmol'
    def _ipython_display_(self):
        obj = self.plot(mode=self.display_mode, return_objects=False)
        return obj._ipython_display_()
        # return self.jupyter_viz()._ipython_display_()
    #endregion

    #region External Program Properties
    def _get_ob_attr(self, item):
        if self._mol is None:
            raise AttributeError("No pybel molecule")
        else:
            return getattr(self._mol, item)
    #endregion

class MolecoolException(Exception):
    pass
