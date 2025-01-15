"""
Provides a simple Molecule class that we can adapt as we need
Most standard functionality should be served by OpenBabel
Uses AtomData to get properties and whatnot
"""

import os, numpy as np
import tempfile

from McUtils.Data import AtomData, UnitsData
from McUtils.Coordinerds import CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D
import McUtils.Numputils as nput
from McUtils.Zachary import Mesh
import McUtils.Plots as plt
from McUtils.ExternalPrograms import RDMolecule

from ..Modes import PrimitiveCoordinatePicker, RedundantCoordinateGenerator

from .MoleculeInterface import *

from .CoordinateSystems import MolecularEmbedding
from .Evaluator import MolecularEvaluator
from .Hamiltonian import MolecularHamiltonian
from .Properties import *

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
                 **metadata
                 ):
        """
        :param atoms: atoms specified by name, either full name or short
        :type atoms: Iterable[str]
        :param coords: coordinates for the molecule, assumed to be in Bohr by default
        :type coords: np.ndarray
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
        self._mass = np.array([a["Mass"] for a in self._ats]) if masses is None else masses
        coords = CoordinateSet(coords, CartesianCoordinates3D)

        internals = self.canonicalize_internals(internals, self.atoms, coords, bonds, masses=self._mass)
        self.embedding = MolecularEmbedding(self.atomic_masses, coords, internals)

        self._name = name

        # properties to be returned

        self._bonds = bonds
        self.guess_bonds = guess_bonds

        self._src = source_file

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
        self._meta = metadata

    @classmethod
    def _auto_spec(cls, atoms, coords, bonds, redundant=False, base_coords=None,
                   masses=None,
                   nonredundant_coordinates=None,
                   prune_coordinates=True,
                   pruning_options=None,
                   **opts):
        if bonds is None:
            bonds = RDMolecule.from_coords(
                                           atoms,
                                           coords * UnitsData.convert("BohrRadius", "Angstroms"),
                                           bonds,
                                           guess_bonds=True
                                           ).bonds
        if redundant and nonredundant_coordinates is None:
            nonredundant_coordinates = base_coords
            base_coords = None
        if redundant and nonredundant_coordinates is not None:
            nonredundant_coordinates = list(nonredundant_coordinates)
            base_coords = nonredundant_coordinates + ([] if base_coords is None else list(base_coords))
        specs = PrimitiveCoordinatePicker(
                atoms,
                [b[:2] for b in bonds],
                base_coords=base_coords,
                **opts
            ).coords
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

            if nonredundant_coordinates is not None:
                spec['untransformed_coordinates'] = np.arange(len(nonredundant_coordinates))
        return spec
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
            prims = spec.get('primitives')
            if prims is not None:
                spec = spec.copy()
                spec['specs'] = prims
                spec['redundant'] = True
                del spec['primitives']
            subspec = spec.get('specs', '')
            if isinstance(subspec, str):
                if subspec.lower() == 'auto':
                    opts = spec.copy()
                    del opts['specs']
                    if 'relocalize' in opts:
                        relocalize = spec.get('relocalize', relocalize)
                        del opts['relocalize']
                    spec = cls._auto_spec(atoms, coords, bonds, masses=masses, **opts)
                else:
                    raise ValueError(f"can't understand internal spec '{spec}'")
            if spec.get('redundant'):
                spec['relocalize'] = spec.get('relocalize', relocalize)
        return spec

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

    def get_internals(self, strip_embedding=True):
        return self.embedding.get_internals(strip_embedding=strip_embedding)

    def get_cartesians_by_internals(self, order=None, strip_embedding=False, **kw):
        return self.embedding.get_cartesians_by_internals(order=order, strip_embedding=strip_embedding, **kw)

    def get_internals_by_cartesians(self, order=None, strip_embedding=False, **kw):
        return self.embedding.get_internals_by_cartesians(order=order, strip_embedding=strip_embedding, **kw)

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
        return self.potential_surface.derivatives
    @potential_derivatives.setter
    def potential_derivatives(self, derivs):
        self.potential_surface.derivatives = derivs

    def get_internal_potential_derivatives(self, order=None):
        raise NotImplementedError("use Numputils tensor reexpand instead")
        if self.potential_derivatives is None:
            raise ValueError("no potential derivatives")
        if order is None:
            order = len(self.potential_derivatives)
        derivs = self.potential_surface.apply_transformation(
            self.embedding.get_cartesians_by_internals(order)
        ).derivatives
        ecs = self.embedding.embedding_coords
        comp_coords = np.setdiff1d(np.arange(derivs[0].shape[0]), ecs)
        _ = []
        for d in derivs:
            for i in range(d.ndim):
                d = np.take(d, comp_coords, axis=i)
            _.append(d)
        derivs = _
        return derivs

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
            raise TypeError("`normal_modes` must be {}".format(
                NormalModesManager.__name__
            ))
        self._normal_modes = val
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
    #endregion

    def __repr__(self):
        return "{cls}('{name}', formula='{formula}', shape={shape}, coord_sys={coord_sys})".format(
            cls=type(self).__name__,
            name=self.name,
            formula=self.formula,
            shape=self.coords.shape,
            coord_sys=self.coords.system.name if isinstance(self.coords, CoordinateSet) else 'undefined'
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
        if min(m) < 100:
            m = m*UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        return m
    @property
    def atomic_masses(self):
        return self._atomic_masses()
    @property
    def bonds(self):
        if self._bonds is None and self.guess_bonds:
            self._bonds = RDMolecule.from_coords(
                                           self.atoms,
                                           self.coords * UnitsData.convert("BohrRadius", "Angstroms"),
                                           None,
                                           guess_bonds=True
                                           ).bonds
            # self._bonds = self.prop("guessed_bonds", tol=1.05, guess_type=True)
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
        return self._src
    @source_file.setter
    def source_file(self, src):
        self._src = src

    @property
    def shape(self):
        return self.coords.shape
    def __len__(self):
        if self.multiconfig:
            return self.coords.shape[0]
        else:
            return 1

    def copy(self):
        import copy
        # mostly just use the default and don't be fancy
        new = copy.copy(self)
        # but we also need to do some stuff where we store objects that
        # reference the molecule
        new.normal_modes = new.normal_modes.copy()
        new.normal_modes.set_molecule(new)
        new.potential_surface = new.potential_surface.copy()
        new.potential_surface.set_molecule(new)
        new.dipole_surface = new.dipole_surface.copy()
        new.dipole_surface.set_molecule(new)
        # new._rdmol = new.rdmol.copy()
        # new.rdmol.set_molecule(new)
        return new

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

    #region Coordinate Embeddings

    @property
    def mass_weighted_coords(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        return self.prop('mass_weighted_coords')
    @property
    def center_of_mass(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        return self.prop('center_of_mass')
    @property
    def inertia_tensor(self):
        """
        :return:
        :rtype: (np.ndarray, np.ndarray)
        """
        return self.prop('inertia_tensor')
    @property
    def inertial_eigensystem(self):
        """
        :return:
        :rtype: (np.ndarray, np.ndarray)
        """
        return self.prop('moments_of_inertia')
    @property
    def moments_of_inertia(self):
        """
        :return:
        :rtype: np.ndarray
        """
        return self.prop('moments_of_inertia')[0]
    @property
    def inertial_axes(self):
        """
        :return:
        :rtype: np.ndarray
        """
        return self.prop('moments_of_inertia')[1]

    @property
    def translation_rotation_modes(self):
        """
        :return:
        :rtype: np.ndarray
        """
        return self.prop('translation_rotation_eigenvectors')

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
    def evaluate(self,
                 func,
                 use_internals=None,
                 deriv_order=None,
                 strip_embedding=False
                 ):
        return self.evaluator.evaluate(
            func,
            use_internals=use_internals,
            deriv_order=deriv_order,
            strip_embedding=strip_embedding
        )
    def evaluate_at(self,
                    func,
                    coords,
                    use_internals=None,
                    deriv_order=None,
                    strip_embedding=False
                    ):
        return self.evaluator.evaluate_at(
            func,
            coords,
            use_internals=use_internals,
            deriv_order=deriv_order,
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
                             which=None, sel=None, axes=None,
                             shift=True,
                             coordinate_expansion=None,
                             strip_embedding=False
                             ):
        return self.evaluator.get_scan_coordinates(
            domains,
            internals=internals,
            which=which, sel=sel, axes=axes,
            shift=shift, coordinate_expansion=coordinate_expansion,
            strip_embedding=strip_embedding
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
            atom_radii = [at['IconRadius'] if c is None else c for c,at in zip(atom_radii, atoms)]
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

    def setup_AIMD(self,
                   potential_function,
                   timestep=.5,
                   seed=None,
                   total_energy=None,
                   trajectories=1,
                   sampled_modes=None,
                   initial_energies=None,
                   initial_displacements=None,
                   displaced_coords=None,
                   track_kinetic_energy=False,
                   track_velocities=False
                   ):
        from ..AIMD import AIMDSimulator

        if initial_displacements is not None:
            init_pos = self.get_displaced_coordinates(
                initial_displacements,
                which=displaced_coords,
                internals='reembed'
            )
            sim = AIMDSimulator(
                self.masses,
                init_pos,
                lambda c: -potential_function(c, deriv_order=1)[1].reshape(c.shape),
                timestep=timestep,
                track_kinetic_energy=track_kinetic_energy,
                track_velocities=track_velocities
            )
        else:
            self.potential_derivatives = potential_function(self.coords, deriv_order=2)[1:]

            if total_energy is not None:
                if seed is not None:
                    np.random.seed(seed)
                freqs = self.normal_modes.modes.freqs
                if sampled_modes is None:
                    sampled_modes = list(range(freqs.shape[0]))
                subdirs = np.random.normal(0, 1, size=(trajectories, len(sampled_modes)))

                dirs = np.zeros((trajectories, freqs.shape[0]))
                dirs[:, sampled_modes] = subdirs
                dirs = dirs / np.linalg.norm(dirs, axis=1)[:, np.newaxis] # random unbiased directions

                dirs = dirs / np.sum(np.abs(dirs), axis=1)[:, np.newaxis] # weights in each dimension
                energies = dirs * freqs[np.newaxis, :]
                initial_energies = total_energy * energies / np.sum(np.abs(energies), axis=1)[:, np.newaxis]

            nms = self.normal_modes.modes.basis
            sim = AIMDSimulator(
                self.atomic_masses,
                [self.coords] * len(initial_energies),
                lambda c: -potential_function(c, deriv_order=1)[1].reshape(c.shape),
                velocities=AIMDSimulator.mode_energies_to_velocities(
                    nms.inverse.T,
                    self.atomic_masses,
                    initial_energies,
                    inverse=nms.matrix.T
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
                  **opts
                  ):
        from ..VPT2 import VPTRunner

        if use_internals or use_internals is None:
            return VPTRunner.construct(self, states, order=order, **opts)
        else:
            return VPTRunner.construct(
                [self.atoms, self.coords],
                states,
                potential_derivatives=self.potential_derivatives,
                modes=self.normal_modes.modes.basis,
                order=order,
                **opts
            )

    @property
    def g_matrix(self):
        """
        Returns the molecular g-matrix for the system
        :return:
        :rtype:
        """
        if self.internal_coordinates is None:
            atma = self._atomic_masses()
            mass_spec = np.broadcast_to(atma[:, np.newaxis], (len(atma), 3)).flatten()
            return np.diag(1 / mass_spec)
            # raise ValueError("need internal coordinates to calculate the G-matrix")
        return self.hamiltonian.gmatrix_expansion(0, modes=None)
        return self.prop('g_matrix')

    @property
    def coriolis_constants(self):
        """
        Returns the molecular g-matrix for the system
        :return:
        :rtype:
        """
        return self.prop('coriolis_constants')

    def bond_length(self, i, j):
        """
        Returns the bond length of the coordinates

        :param i:
        :type i:
        :param j:
        :type j:
        :return:
        :rtype:
        """
        return nput.pts_norms(self.coords[..., i, :], self.coords[..., j, :])
    def bond_angle(self, i, j, k):
        """
        Returns the bond angle of the specified coordinates

        :param i:
        :type i:
        :param j:
        :type j:
        :return:
        :rtype:
        """
        return nput.pts_angles(self.coords[..., i, :], self.coords[..., j, :], self.coords[..., k, :])[0]
    def dihedral(self, i, j, k, l):
        """
        Returns the dihedral angle of the specified coordinates

        :param i:
        :type i:
        :param j:
        :type j:
        :return:
        :rtype:
        """
        return nput.pts_dihedrals(
            self.coords[..., i, :], self.coords[..., j, :],
            self.coords[..., k, :], self.coords[..., l, :]
        )[0]

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

    def eckart_frame(self, mol, sel=None, inverse=False, planar_ref_tolerance=None):
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
        return self.prop('eckart_transformation', mol, sel=sel, inverse=inverse, planar_ref_tolerance=planar_ref_tolerance)

    def embed_coords(self, crds, sel=None, in_paf=False, planar_ref_tolerance=None):
        """
        Embeds coords in the Eckart frame using `self` as a reference

        :param crds:
        :type crds:
        :return:
        :rtype:
        """

        return self.prop('eckart_embedded_coords', crds, sel=sel, in_paf=in_paf, planar_ref_tolerance=planar_ref_tolerance)
    def get_embedding_data(self, crds, sel=None, in_paf=False, planar_ref_tolerance=None):
        """
        Gets the necessary data to embed crds in the Eckart frame using `self` as a reference

        :param crds:
        :type crds:
        :return:
        :rtype:
        """
        return self.prop('eckart_embedding_data', crds, sel=sel, in_paf=in_paf, planar_ref_tolerance=planar_ref_tolerance)
    def get_embedded_molecule(self,
                              ref=None,
                              sel=None, planar_ref_tolerance=None,
                              embed_properties=True,
                              load_properties=True
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
            frame = self.eckart_frame(ref, sel=sel, planar_ref_tolerance=planar_ref_tolerance, inverse=False)
        # self.normal_modes.modes
        new = frame.apply(self)
        if embed_properties:
            # inv_frame = frame
            # if ref is None:
            #     inv_frame = self.principle_axis_frame(inverse=True)
            # else:
            #     inv_frame = self.eckart_frame(ref, inverse=True)
            if load_properties or new._normal_modes._modes is not None:
                new.normal_modes = new.normal_modes.apply_transformation(frame)
            if load_properties or new._pes._surf is not None or new._pes._derivs is not None:
                new.potential_surface = new.potential_surface.apply_transformation(frame)
            if load_properties or new._dips._surf is not None or new._dips._derivs is not None:
                new.dipole_surface = new.dipole_surface.apply_transformation(frame)
        new.source_file = None # for safety

        return new

    #endregion

    #region Input Formats
    @property
    def rdmol(self):
        from McUtils.ExternalPrograms import RDMolecule

        return RDMolecule.from_mol(self, coord_unit="BohrRadius")

    @classmethod
    def from_zmat(cls, zmat, **opts):
        """Little z-matrix importer

        :param zmat:
        :type zmat: str | tuple
        :return:
        :rtype: Molecule
        """

        if isinstance(zmat, str):
            from McUtils.Parsers import StringParser, ZMatPattern
            zmcs = StringParser(ZMatPattern).parse(zmat)
        else:
            zmcs = zmat
        coords = CoordinateSet([zmcs[1]], ZMatrixCoordinates).convert(CartesianCoordinates3D)
        return cls(zmcs[0], coords, zmatrix=zmat[:, (1, 3, 5)], **opts)
    @classmethod
    def from_pybel(cls, mol, **opts):
        """

        :param mol:
        :type mol: pybel.mol
        :return:
        :rtype:
        """
        from McUtils.ExternalPrograms import OpenBabelInterface

        ob = OpenBabelInterface().openbabel
        bonds = list(ob.OBMolBondIter(mol.OBMol))
        atoms = list(mol.atoms)

        opts = dict({'bonds':bonds, 'obmol':mol}, **opts)
        return cls(atoms, [a.coords for a in atoms], **opts)
    @classmethod
    def _from_log_file(cls, file, num=None, **opts):
        from McUtils.GaussianInterface import GaussianLogReader
        with GaussianLogReader(file) as gr:
            parse = gr.parse('StandardCartesianCoordinates', num=num)
        spec, coords = parse['StandardCartesianCoordinates']
        ang2bohr = UnitsData.convert("Angstroms", "AtomicUnitOfLength")
        return cls(
            [int(a[1]) for a in spec],
            CoordinateSet(ang2bohr*np.array(coords), CartesianCoordinates3D),
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
    def from_rdmol(cls, rdmol, **opts):
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
    def _from_smiles(cls, smi, add_implicit_hydrogens=False, **opts):
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_smiles(smi, add_implicit_hydrogens=add_implicit_hydrogens), **opts)
    @classmethod
    def _from_sdf(cls, sdf, **opts):
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_sdf(sdf), **opts)
    @classmethod
    def _from_molblock(cls, sdf, **opts):
        from McUtils.ExternalPrograms import RDMolecule

        return cls.from_rdmol(RDMolecule.from_molblock(sdf), **opts)

    @classmethod
    def from_string(cls, string, fmt, **opts):
        format_dispatcher = {
            "smi": cls._from_smiles,
            "mol": cls._from_molblock,
            "sdf": cls._from_sdf
        }

        if fmt in format_dispatcher:
            return format_dispatcher[fmt](string, **opts)
        else:
            with tempfile.NamedTemporaryFile("w+") as tf:
                tf.write(string)
                new = cls.from_file(tf.name, mode=fmt, **opts)
            new.source_file = None # don't want to keep this lying around...
        return new

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
        format_dispatcher = {
            "log": cls._from_log_file,
            "fchk": cls._from_fchk_file,
            "smi": cls._from_smiles,
            "mol":cls._from_molblock,
            "sdf":cls._from_sdf
        }

        if mode == None:
            path, ext = os.path.splitext(file)
            ext = ext.lower()
            mode = ext.strip(".")

        if mode in format_dispatcher:
            loader = format_dispatcher[mode]
            return loader(file, **opts)
        else:
            from McUtils.ExternalPrograms import OpenBabelInterface
            try:
                pybel = OpenBabelInterface().pybel
            except ImportError:
                pybel = None
            if pybel is None:
                raise IOError("{} doesn't support file type {} without OpenBabel installed.".format(cls.__name__, mode))
            else:
                mol = next(pybel.readfile(mode, file))
                return cls.from_pybel(mol)

    @classmethod
    def _infer_spec_format(cls, spec):
        if isinstance(spec, str):
            return 'file'

        fmt = None
        try:
            atoms = spec[0]
            if all(isinstance(a, str) for a in atoms):
                return 'standard'
        except:
            pass

        if fmt is None:
            raise ValueError("don't know how to build a molecule from spec {}".format(spec))

    @classmethod
    def from_spec(cls, spec):
        fmt = cls._infer_spec_format(spec)
        if fmt == 'file':
            return cls.from_file(spec)
        elif fmt == 'standard':
            atoms = spec[0]
            coords = spec[1]
            if len(spec) == 2:
                opts = {}
            elif len(spec) == 3:
                opts = spec[2]
            else:
                raise ValueError("too many arguments in {} to build {}".format(spec, cls.__name__))
            return cls(atoms, coords, **opts)
        elif fmt == 'zmat':
            if isinstance(spec[0], str):
                opts = {}
                zmat = spec
            else:
                if len(spec) == 1:
                    opts = {}
                    zmat = spec[0]
                elif len(spec) == 2:
                    opts = spec[1]
                    zmat = spec[0]
                else:
                    raise ValueError("too many arguments in {} to build {}".format(spec, cls.__name__))
            return cls.from_zmat(zmat, **opts)
        else:
            raise NotImplementedError("don't have {} loading from format {} for spec {}".format(cls.__name__, fmt, spec))
    #endregion

    #region Visualization
    highlight_styles = {
        "glow":"green",
        "color":"white"
    }
    def plot(self,
             *geometries,
             figure=None,
             return_objects=True,
             bond_radius=.1,
             atom_radius_scaling=.25,
             atom_style=None,
             bond_style=None,
             highlight_atoms=None,
             highlight_bonds=None,
             highlight_rings=None,
             highlight_styles=None,
             mode='quality',
             backend='x3d',
             objects=False,
             graphics_class=None,
             cylinder_class=None,
             sphere_class=None,
             animate=None,
             animation_options=None,
             jsmol_load_script=None,
             **plot_ops
             ):

        if mode == 'jupyter':
            return self.jupyter_viz()
        elif mode == 'jsmol':
            return self.jsmol_viz(script=jsmol_load_script)

        from McUtils.Plots import Graphics3D, Sphere, Cylinder, Line, Disk

        graphics_keys = Graphics3D.known_keys | Graphics3D.opt_keys | Graphics3D.figure_keys
        graphics_opts = {k:plot_ops[k] for k in plot_ops.keys() & graphics_keys}
        plot_ops = {k:plot_ops[k] for k in plot_ops.keys() - graphics_keys}

        if graphics_class is None:
            graphics_class = Graphics3D
        if cylinder_class is None:
            cylinder_class = Line if mode == 'fast' else Cylinder
        if sphere_class is None:
            sphere_class = Disk if mode == 'fast' else Sphere

        if len(geometries) == 0:
            geometries = self.coords
        elif len(geometries) == 1:
            geometries = CoordinateSet(np.asanyarray(geometries[0]), self.coords.system)
        else:
            geometries = CoordinateSet([np.asanyarray(g) for g in geometries], self.coords.system)

        if geometries.ndim == 2:
            geometries = geometries[np.newaxis]
        if animate is None:
            animate = geometries.shape[0] > 1

        geometries = geometries.convert(CartesianCoordinates3D)

        if figure is None:
            figure = graphics_class(backend=backend, **graphics_opts)

        colors = [ at["IconColor"] for at in self._ats ]
        radii = [ atom_radius_scaling * at["IconRadius"] for at in self._ats ]

        bonds = [None] * len(geometries)
        atoms = [None] * len(geometries)

        if atom_style is None:
            atom_style = {}
        elif not isinstance(atom_style, dict):
            atom_style = {i:a for i,a in enumerate(atom_style)}
        base_atom_style = {}
        _atom_style = {}
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

        if bond_style is None:
            bond_style = {}
        elif not isinstance(bond_style, dict):
            bond_style = {i:a for i,a in enumerate(bond_style)}
        base_bond_style = {}
        _bond_style = {}
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
                    # i'm not supporting double or triple bonds for not because they're just too much of a pain...
                    # in Mathematica I have some clever code for finding the midpoint between the vectors to draw to
                    # I'm not gonna do that for now because it's not worth it for the quick-and-dirty stuff I want to do
                    p1 = geom[atom1]
                    p2 = geom[atom2]
                    midpoint = (p2 - p1)/2 + p1
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
                    cc1 = cylinder_class(
                        p1,
                        midpoint,
                        bond_radius,
                        **plot_ops,
                        **b_sty_1
                    )

                    b_sty_2 = dict(
                        bond_style.get(atom2, {}),
                        **base_bstyle
                    )
                    if b_sty_2.get('color') is None:
                        b_sty_2['color'] = c2
                    cc2 = cylinder_class(
                        midpoint,
                        p2,
                        bond_radius,
                        **plot_ops,
                        **b_sty_2
                    )
                    if objects:
                        bonds[i][j] = (( cc1, cc2 ))
                    else:
                        cyl_1 = cc1.plot(figure)
                        cyl_2 = cc2.plot(figure)
                        if isinstance(cyl_1, (list, tuple)):
                            cyl_1 = cyl_1[0]
                        if isinstance(cyl_2, (list, tuple)):
                            cyl_2 = cyl_2[0]
                        bonds[i][j] = (( cyl_1, cyl_2 ))

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

        if animate:
            if animation_options is None: animation_options = {}
            figure = figure.animate_frames(
                [
                    a + sum([list(b) for b in bl], [])
                    for a,bl in zip(atoms, bonds)
                ],
                **animation_options
            )

        if return_objects:
            return figure, atoms, bonds
        else:
            return figure

    def get_animation_geoms(self, which, extent=.35, steps=8, strip_embedding=True, units=None):
        if isinstance(which, int):
            coord_expansion = self.get_cartesians_by_internals(2, strip_embedding=strip_embedding)
        else:
            if isinstance(which, np.ndarray):
                if which.ndim == 1:
                    which = which[np.newaxis]
                which = [which]
            coord_expansion = which
            which = 0
        geoms = self.get_scan_coordinates(
            [[-extent, extent, steps]],
            which=[which],
            coordinate_expansion=coord_expansion
        )
        geoms = np.concatenate([geoms, np.flip(geoms, axis=0)], axis=0)
        if units is not None:
            geoms = geoms * UnitsData.convert("BohrRadius", units)
        return geoms
    def animate_coordinate(self, which, extent=.35, steps=8, return_objects=False, strip_embedding=True,
                           units=None,
                           mode='x3d',
                           jsmol_load_script=None,
                           **plot_opts
                           ):
        if mode == 'jsmol':
            disps = self.format_animation_file(which, format="jmol",
                                               extent=extent, steps=steps, strip_embedding=strip_embedding,
                                               units="Angstroms" if units is None else units)
            return self.jsmol_viz(disps, vibrate=True, script=jsmol_load_script)
        else:
            geoms = self.get_animation_geoms(which, extent=extent, steps=steps, strip_embedding=strip_embedding, units=units)
            return self.plot(geoms, return_objects=return_objects, **plot_opts)

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

    def format_animation_file(self, which, format='xyz', extent=.35, steps=8, strip_embedding=True, units='Angstroms'):
        if format == 'jmol':
            coord_expansion = self.get_cartesians_by_internals(1, strip_embedding=strip_embedding)[0]
            return self._format_jmol_displacements_files(
                self.coords,
                [coord_expansion[(which,), :] * extent],
                units
            )[0]
        else:
            geoms = self.get_animation_geoms(which, extent=extent, steps=steps, strip_embedding=strip_embedding, units=units)
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
        return self.jupyter_viz().to_widget()

    display_jsmol = True
    def _ipython_display_(self):
        if self.display_jsmol:
            obj = self.plot(mode='jsmol')
        else:
            obj = self.plot(backend='x3d', return_objects=False)
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
