"""
Provides a simple Molecule class that we can adapt as we need
Most standard functionality should be served by OpenBabel
Uses AtomData to get properties and whatnot
"""

import os, numpy as np

from McUtils.Data import AtomData, UnitsData
from McUtils.Coordinerds import CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D, CompositeCoordinateSystem
from McUtils.Parallelizers import Parallelizer
from McUtils.Zachary import TensorDerivativeConverter
import McUtils.Numputils as nput

from .MoleculeInterface import *
from .CoordinateSystems import (
    MolecularCartesianCoordinateSystem, MolecularZMatrixCoordinateSystem,
    MolecularZMatrixToCartesianConverter, MolecularCartesianToZMatrixConverter,
    MolecularCartesianToRegularCartesianConverter, RegularCartesianToMolecularCartesianConverter,
    MolecularZMatrixToRegularZMatrixConverter, RegularZMatrixToMolecularZMatrixConverter
)
from .Properties import *

__all__ = [
    "Molecule",
    "MolecoolException"
]

__reload_hook__ = [".MoleculeInterface", '.Properties']

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
                 obmol=None,
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
        self._mass = masses

        coords = CoordinateSet(coords, CartesianCoordinates3D)

        self._name = name

        # properties to be returned
        self._coords = coords
        self._sys = MolecularCartesianCoordinateSystem(self)
        self._coords = CoordinateSet(self._coords, self._sys)
        MolecularCartesianToRegularCartesianConverter(self.coords.system).register()
        RegularCartesianToMolecularCartesianConverter(self.coords.system).register()
        if isinstance(internals, CoordinateSet):
            self._int_spec = None
            self._ints = internals
        else:
            self._int_spec = self.canonicalize_internal_coordinate_spec(internals)
            self._ints = None
        self._jacobians = {
            'internal':[],
            'cartesian':[]
        }

        self._bonds = bonds

        self._src = source_file

        self.ext_mol = OpenBabelMolManager(self, obmol)
        self._dips = DipoleSurfaceManager(self,
                                                   surface=dipole_surface,
                                                   derivatives=dipole_derivatives
                                                   )
        self._pes = PotentialSurfaceManager(self,
                                                   surface=potential_surface,
                                                   derivatives=potential_derivatives
                                                   )

        self._normal_modes = NormalModesManager(self, normal_modes=normal_modes)

        metadata['charge'] = charge
        self._meta = metadata

        self.guess_bonds=guess_bonds

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
    @property
    def masses(self):
        if self._mass is None:
            return np.array([a["Mass"] for a in self._ats])
        else:
            return self._mass
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
            self._bonds = self.prop("guessed_bonds", tol=1.05, guess_type=True)
        return self._bonds
    @property
    def coords(self):
        return self._coords
    @coords.setter
    def coords(self, new):
        if not isinstance(new, CoordinateSet):
            new = CoordinateSet(new, self.sys)
        self._ints = None
        self._coords = new
    @property
    def sys(self):
        return self._coords.system
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

    #region Structure Modification
    def insert_atoms(self, atoms, coords, where, handle_properties=True):
        new = self.copy()

        #awkwardly these need to come first...?
        if handle_properties:
            new.normal_modes = new.normal_modes.insert_atoms(atoms, coords, where)
            new.dipole_surface = new.dipole_surface.insert_atoms(atoms, coords, where)
            new.potential_surface = new.potential_surface.insert_atoms(atoms, coords, where)

        new._coords = np.insert(new.coords, where, coords,
                                axis=1 if self.multiconfig else 0
                                )
        new._ats = np.insert(np.array(new._ats, dtype=object),
                             where,
                             [AtomData[atom] if isinstance(atom, (int, np.integer, str)) else atom for atom in atoms]
                             )
        if new._mass is not None:
            new._mass = np.insert(np.array(new._mass),
                             where,
                             [AtomData[atom]["Mass"] if isinstance(atom, (int, np.integer, str)) else atom["Mass"] for atom in atoms]
                             )

        new.coords.system = MolecularCartesianCoordinateSystem(new)

        return new

    def delete_atoms(self, where, handle_properties=True):
        new = self.copy()
        new._coords = np.delete(new.coords, where,
                                axis=1 if self.multiconfig else 0
                                )
        new._ats = np.delete(np.array(new._ats, dtype=object),
                             where
                             )
        new.coords.system = MolecularCartesianCoordinateSystem(new)
        if handle_properties:
            new.dipole_surface = new.dipole_surface.delete_atoms(where)
            new.pes = new.pes.delete_atoms(where)
            new.normal_modes = new.normal_modes.delete_atoms(where)

        return new

    def take_submolecule(self, spec):
        """
        Takes a 'slice' of a molecule if working with Cartesian coords.
        If not, need to do some corner case handling for that.

        :param spec:
        :type spec:
        :return:
        :rtype:
        """
        new_coords = self.coords[spec]
        new_shape = new_coords.shape
        cur_shape = self.coords.shape
        # if we're no longer working with Cartesians, then we say "Abort!"
        if new_shape[-1] != 3:
            return new_coords
        elif new_shape[-2] != cur_shape[-2]:
            # we have a different number of atoms now...
            raise IndexError("I haven't implemented slicing for molecules that changes the # of atoms")
        else:
            new = self.copy()
            new._coords = new_coords
            return new
    #endregion

    @property
    def shape(self):
        return self.coords.shape
    def __len__(self):
        if self.multiconfig:
            return self.coords.shape[0]
        else:
            return 1
    def __iter__(self):
        if self.multiconfig:
            for i in range(len(self)):
                yield self[i]
        else:
            yield self
    def __getitem__(self, item):
        return self.take_submolecule(item)

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
        new.ext_mol = new.ext_mol.copy()
        new.ext_mol.set_molecule(new)
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
    @classmethod
    def canonicalize_internal_coordinate_spec(cls, spec):
        if spec is not None:
            if hasattr(spec, 'items'):
                try:
                    zmatrix = spec['zmatrix']
                except KeyError:
                    zmatrix = None
                else:
                    zmatrix = np.asanyarray(zmatrix).astype(int)
                    if zmatrix.shape[1] != 4:
                        raise ValueError("can't understand Z-matrix {}".format(zmatrix))
                spec['zmatrix'] = zmatrix
                try:
                    conversion = spec['conversion']
                except KeyError:
                    conversion = None
                spec['conversion'] = cls._wrap_conv(conversion)
                try:
                    inverse = spec['inverse']
                except KeyError:
                    inverse = None
                spec['inverse'] = cls._wrap_conv(inverse)
                try:
                    converter_options = spec['converter_options']
                except KeyError:
                    converter_options = {}
                else:
                    if converter_options is None:
                        converter_options = {}
                if 'embedding_coords' not in converter_options:
                    if spec['zmatrix'] is not None:
                        converter_options['embedding_coords'] = MolecularZMatrixCoordinateSystem.embedding_coords
                if 'jacobian_prep' not in converter_options:
                    if spec['zmatrix'] is not None:
                        converter_options['jacobian_prep'] = ZMatrixCoordinates.jacobian_prep_coordinates
                spec['converter_options'] = converter_options
            elif callable(spec):
                zmatrix = None
                conversion = cls._wrap_conv(spec)
                spec = {'zmatrix':zmatrix, 'conversion':conversion, 'inverse':None, 'converter_options':{}}
            else:
                conversion = None
                zmatrix = np.asanyarray(spec).astype(int)
                if zmatrix.shape[1] != 4:
                    raise ValueError("can't understand Z-matrix {}".format(zmatrix))
                spec = {'zmatrix':zmatrix, 'conversion':conversion, 'inverse':None, 'converter_options':{}}
        return spec
    @staticmethod
    def _wrap_conv(f):
        if f is None:
            return f
        def wrapped(*args, **kwargs):
            vals = f(*args, **kwargs)
            if not isinstance(vals, np.ndarray):
                vals, opts = vals
            else:
                opts = {}
            opts = dict(kwargs, **opts)
            return vals, opts
        return wrapped

    @property
    def internals(self):
        if self._int_spec is not None:
            return self._int_spec
    @internals.setter
    def internals(self, internals):
        self._int_spec = self.canonicalize_internal_coordinate_spec(internals)
        self._ints = None
    @property
    def zmatrix(self):
        if self._int_spec is not None:
            return self._int_spec['zmatrix']
    @zmatrix.setter
    def zmatrix(self, zmat):
        if zmat is not None:
            zmat = np.asanyarray(zmat).astype(int)
            if zmat.shape[1] != 4:
                raise ValueError("can't understand Z-matrix {}".format(zmat))
        if self._int_spec is None:
            self._int_spec = self.canonicalize_internal_coordinate_spec(zmat)
        else:
            self._int_spec['zmatrix'] = zmat
        self._ints = None
    @property
    def internal_coordinates(self):
        if self._ints is None and (
                self._int_spec is not None
                and (self._int_spec['zmatrix'] is not None or self._int_spec['conversion'] is not None)
        ):
            coords = self.coords
            if self._int_spec['zmatrix'] is not None:
                zms = MolecularZMatrixCoordinateSystem(self, ordering=self._int_spec['zmatrix'])
                MolecularCartesianToZMatrixConverter(self.coords.system, zms).register()
                MolecularZMatrixToCartesianConverter(zms, self.coords.system).register()
                MolecularZMatrixToRegularZMatrixConverter(zms).register()
                RegularZMatrixToMolecularZMatrixConverter(zms).register()
                coords = self.coords.convert(zms)
            if self._int_spec['conversion'] is not None:
                conv = CompositeCoordinateSystem.register(
                    coords.system,
                    self._int_spec['conversion'],
                    inverse_conversion=self._int_spec['inverse'],
                    **self._int_spec['converter_options']
                )
                coords = coords.convert(conv)
            # print(zms)
            # print(zms, self.coords, self.coords.system.converter(zms))
            self._ints = coords
        return self._ints
    @internal_coordinates.setter
    def internal_coordinates(self, ics):
        if not isinstance(ics, CoordinateSet):
            raise ValueError("{} must be a {} to be valid internal coordinates".format(
                ics, CoordinateSet.__name__
            ))
        self._ints = ics

    def _get_int_jacobs(self,
                       jacs,
                       strip_dummies=False,
                       stencil=None, mesh_spacing=1.0e-3,
                       all_numerical=True, reembed=False,
                       planar_ref_tolerance=None,
                       parallelizer=None
                       ):
        """
        Gets the specified dX/dRs

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        intcds = self.internal_coordinates
        ccoords = self.coords
        carts = ccoords.system
        internals = intcds.system

        if isinstance(jacs, int):
            jacs = list(range(1, jacs + 1))

        exist_jacs = self._jacobians['internal']
        max_jac = max(jacs)
        need_jacs = [x+1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
        if len(need_jacs) > 0:
            stencil = (max(need_jacs) + 2 + (1+max(need_jacs))%2) if stencil is None else stencil
            # odd behaves better
            with Parallelizer.lookup(parallelizer) as par:
                new_jacs = [
                    x.squeeze() if isinstance(x, np.ndarray) else x
                    for x in intcds.jacobian(carts, need_jacs,
                                             # odd behaves better
                                             mesh_spacing=mesh_spacing,
                                             stencil=stencil,
                                             all_numerical=all_numerical,
                                             converter_options=dict(
                                                 reembed=reembed,
                                                 planar_ref_tolerance=planar_ref_tolerance,
                                                 strip_dummies=strip_dummies
                                             ),
                                             parallelizer=par
                                             )
                ]
                # np.set_printoptions
                # with np.printoptions(linewidth=1e8, threshold=1e8, floatmode='fixed', precision=10):
                #     raise Exception(str(np.round(new_jacs[0].reshape(9, 9)[(3, 6, 7), :], 12)))
            for j,v in zip(need_jacs, new_jacs):
                for d in range(j-len(exist_jacs)):
                    exist_jacs.append(None)
                exist_jacs[j-1] = v

        return [exist_jacs[j-1] for j in jacs]

    def _get_cart_jacobs(self, jacs,
                         strip_dummies=False,
                         stencil=None, mesh_spacing=1.0e-3,
                         all_numerical=True,
                         parallelizer=None
                         ):
        """
        Gets the specified dR/dXs

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        intcds = self.internal_coordinates
        ccoords = self.coords
        carts = ccoords.system
        internals = intcds.system

        if isinstance(jacs, int):
            jacs = list(range(1, jacs + 1))

        exist_jacs = self._jacobians['cartesian']
        max_jac = max(jacs)
        need_jacs = [x+1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
        if len(need_jacs) > 0:
            stencil = (max(need_jacs) + 2 + (1+max(need_jacs))%2) if stencil is None else stencil
            # odd behaves better
            with Parallelizer.lookup(parallelizer) as par:
                new_jacs = [
                    x.squeeze() if isinstance(x, np.ndarray) else x
                    for x in ccoords.jacobian(internals, need_jacs,
                                                          mesh_spacing=mesh_spacing,
                                                          stencil=stencil,
                                                          all_numerical=all_numerical,
                                                          converter_options=dict(strip_dummies=strip_dummies),
                                                          parallelizer=par
                                                          )
                ]

                for j, v in zip(need_jacs, new_jacs):
                    for d in range(j - len(exist_jacs)):
                        exist_jacs.append(None)
                    exist_jacs[j - 1] = v

        return [exist_jacs[j-1] for j in jacs]

    def _get_embedding_coords(self):
        try:
            embedding = self.internal_coordinates.system.embedding_coords
        except AttributeError:
            try:
                embedding = self.internal_coordinates.system.converter_options['embedding_coords']
            except KeyError:
                embedding = None
        return embedding

    def get_cartesians_by_internals(self, order=None, strip_embedding=False):
        base = self._get_int_jacobs(order) if order is not None else self._jacobians['internals']
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'CartesiansByInternals',
                    len(base),
                    order
                ))
            base = base[:order]

        _ = []
        sh = self.coords.shape[:-2]
        nc = 3 * len(self.atoms)
        for i, b in enumerate(base):
            b = b.reshape(sh + (nc,) * (i + 2))
            _.append(b)
        base = _

        embedding_coords = self._get_embedding_coords() if strip_embedding else None
        if embedding_coords is not None and strip_embedding:
            good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
            base = [t[np.ix_(*((good_coords,) * (t.ndim - 1)))] for t in base]
        return base
    def get_internals_by_cartesians(self, order=None, strip_embedding=False):
        base = self._get_cart_jacobs(order) if order is not None else self._jacobians['cartesian']
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'InternalsByCartesians',
                    len(base),
                    order
                ))
            base = base[:order]

        _ = []
        sh = self.coords.shape[:-2]
        nc = 3*len(self.atoms)
        for i,b in enumerate(base):
            b = b.reshape(sh+(nc,)*(i+2))
            _.append(b)
        base = _

        if strip_embedding:
            embedding_coords = self._get_embedding_coords()
            if embedding_coords is not None:
                good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
                base = [t[..., good_coords] for t in base]
        return base

    def evaluate(self,
                 func,
                 internals=None,
                 deriv_order=None,
                 strip_embedding=False
                 ):
        if internals is None:
            internals = self.internals is not None
        if internals:
            coords = self.internal_coordinates
            if strip_embedding:
                embedding_coords = self._get_embedding_coords()
                if embedding_coords is not None:
                    good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
                    coords = coords.reshape(coords.shape[:-2] + (3 * len(self.masses),))
                    coords = coords[..., good_coords]
            if deriv_order is None:
                return func(coords).view(np.ndarray)

            terms = func(coords, deriv_order=deriv_order)
            if strip_embedding:
                embedding_coords = self._get_embedding_coords()
                if embedding_coords is not None:
                    good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)

                    const = terms[0]
                    terms = terms[1:]
                    new = []
                    ncs = 3 * len(self.masses)
                    if self.coords.ndim > 2:
                        npts = len(coords)
                    else:
                        npts = 1
                    for n, ders in enumerate(terms):
                        dt = np.zeros((npts,) + (ncs,)*(n+1))
                        idx_pos = (...,) + np.ix_(*[good_coords]*(n+1))
                        dt[idx_pos] = ders.view(np.ndarray)
                        if self.coords.ndim == 2:
                            dt = dt[0]
                        new.append(dt)
                    terms = [const] + new

            const = terms[0]
            jacs = self.get_internals_by_cartesians(deriv_order)

            terms = TensorDerivativeConverter(
                jacs,
                terms[1:],
                jacobians_name='dXdR',
                values_name='f'
            ).convert()  # , check_arrays=True)

            return [const.view(np.ndarray)] + [t.view(np.ndarray) for t in terms]
        else:
            if deriv_order is None:
                return func(self.coords).view(np.ndarray)
            else:
                return [x.view(np.ndarray) for x in func(self.coords, deriv_order=deriv_order)]

    def evaluate_at(self,
                    func,
                    coords,
                    internals=None,
                    deriv_order=None,
                    strip_embedding=False
                    ):
        return type(self)(
            self.atoms,
            coords,
            internals=self.internals
        ).evaluate(func, internals=internals, deriv_order=deriv_order, strip_embedding=strip_embedding)

    def get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None, internals=False, shift=True):
        displacements = np.asanyarray(displacements)

        if which is not None:
            which = tuple(
                np.ravel_multi_index(idx, (3, len(self._ats)))
                    if not isinstance(idx, (int, np.integer)) else
                idx
                for idx in which
            )

        if internals and self.internals is None:
            raise ValueError("can't displace in internals without internal coordinate spec")
        base_coords = self.coords if not internals else self.internal_coordinates

        if which is not None:
            if displacements.shape[-1] != len(which):  # displacements provided in atom coordinates
                displacements = displacements.reshape(
                    displacements.shape[:-2] +
                    (np.prod(displacements.shape[-2:], dtype=int),)
                )
            if displacements.ndim > 1:
                for _ in range(displacements.ndim - 1):
                    base_coords = np.expand_dims(base_coords, 0)
                base_coords = np.broadcast_to(base_coords, displacements.shape[:-1] + base_coords.shape[-2:])
            base_coords = base_coords.copy()
            flat_coords = base_coords.reshape(displacements.shape[:-1] + (-1,))

            if shift:
                flat_coords[..., which] += displacements
            else:
                flat_coords[..., which] = displacements
            base_coords = flat_coords.reshape(base_coords.shape)

        elif sel is not None or axes is not None:
            if displacements.ndim > 2:
                for _ in range(displacements.ndim - 2):
                    base_coords = np.expand_dims(base_coords, 0)
                base_coords = np.broadcast_to(base_coords, displacements.shape[:-2] + base_coords.shape[-2:])
            base_coords = base_coords.copy()

            if sel is None:
                sel = np.arange(len(self.masses))
            if axes is None:
                axes = np.arange(3)

            sel = np.asanyarray(sel)[:, np.newaxis]
            axes = np.asanyarray(axes)[np.newaxis, :]

            if shift:
                base_coords[..., sel, axes] += displacements
            else:
                base_coords[..., sel, axes] = displacements

        else:
            if displacements.shape[-2:] != base_coords.shape:
                raise ValueError("displacements with shape {} passed but coordinates have shape {}".format(
                    displacements.shape,
                    base_coords.shape
                ))
            base_coords = displacements

        if internals:
            # track the embedding info...
            base_coords = self.internal_coordinates.system(base_coords, **self.internal_coordinates.converter_options)
            if isinstance(internals, str):
                if internals == 'convert':
                    base_coords = base_coords.convert(self.coords.system)
                elif internals == 'reembed':
                    base_coords = self.embed_coords(base_coords.convert(self.coords.system))
        else:
            base_coords = self.coords.system(base_coords)
        return base_coords

    def get_scan_coordinates(self,
                             domains,
                             internals=False,
                             which=None, sel=None, axes=None,
                             shift=True
                             ):

        displacement_mesh = np.moveaxis(
            np.array(
                np.meshgrid(*[np.linspace(*d) for d in domains], indexing='ij')
            ),
            0, -1
        )
        return self.get_displaced_coordinates(displacement_mesh, shift=shift,
                                              internals=internals, which=which, sel=sel, axes=axes)

    def get_nearest_displacement_coordinates(self, points, sel=None, axes=None, weighting_function=None):
        pts = np.asanyarray(points)
        smol = pts.ndim == 1
        if smol: pts = pts[np.newaxis]
        base_shape = pts.shape[:-1]
        pts = pts.reshape(-1, pts.shape[-1])

        if axes is None:
            axes = [0, 1, 2]
        axes = np.asanyarray(axes)

        if sel is None:
            sel = np.arange(len(self.masses))
        sel = np.asanyarray(sel)

        ref = self.coords
        masses = self.masses
        dists = np.linalg.norm(
            pts[:, np.newaxis, :] - ref[np.newaxis, sel[:, np.newaxis], axes[np.newaxis, :]],
            axis=-1
        )
        if weighting_function is None:
            weighting_function = np.sqrt
        dists = dists * weighting_function(masses)[np.newaxis, sel]
        nearest = np.argmin(dists, axis=-1)
        atom_idx = sel[nearest]
        coords = np.broadcast_to(ref[np.newaxis], (len(pts),) + ref.shape).copy()
        coords[np.arange(len(nearest))[:, np.newaxis], atom_idx[:, np.newaxis], axes[np.newaxis, :]] = pts

        coords = coords.reshape(base_shape + coords.shape[-2:])
        if smol:
            coords = coords[0]

        return coords

    def get_nearest_scan_coordinates(self, domains, sel=None, axes=None):
        displacement_mesh = np.moveaxis(
            np.array(
                np.meshgrid(*[np.linspace(*d) for d in domains], indexing='ij')
            ),
            0, -1
        )
        return self.get_nearest_displacement_coordinates(displacement_mesh, axes=axes, sel=sel)

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

    @property
    def g_matrix(self):
        """
        Returns the molecular g-matrix for the system
        :return:
        :rtype:
        """
        if self.internal_coordinates is None:
            raise ValueError("need internal coordinates to calculate the G-matrix")
        return self.prop('g_matrix')

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
        :rtype: tuple[np.ndarray, tuple[np.ndarray], tuple[np.ndarray]]
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
            "fchk": cls._from_fchk_file
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
    def plot(self,
             *geometries,
             figure=None,
             bond_radius=.1,
             atom_radius_scaling=.25,
             atom_style=None,
             bond_style=None,
             mode='fast',
             objects=False,
             **plot_ops
             ):

        if mode == 'jupyter':
            return self.jupyter_viz()

        from McUtils.Plots import Graphics3D, Sphere, Cylinder, Line, Disk

        if len(geometries) == 0:
            geometries = CoordinateSet([self._coords], self._coords.system)
        elif len(geometries) == 1 and isinstance(geometries[0], CoordinateSet):
            geometries = geometries[0]
        else:
            geometries = CoordinateSet(geometries, self._coords.system)

        # from McUtils.Coordinerds import CoordinateSystemConverters
        # raise Exception(">>>>>", [(k.__name__, b.__name__) for k, b in CoordinateSystemConverters.converters])

        geometries = geometries.convert(CartesianCoordinates3D)

        if figure is None:
            #backend = "VTK" -- once the VTK backend is more fleshed out we can use it...
            figure = Graphics3D(**plot_ops)

        colors = [ at["IconColor"] for at in self._ats ]
        radii = [ atom_radius_scaling * at["IconRadius"] for at in self._ats ]

        bonds = [None] * len(geometries)
        atoms = [None] * len(geometries)

        if atom_style is None:
            atom_style = {}
        if bond_style is None:
            bond_style = {}

        c_class = Line if mode == 'fast' else Cylinder
        s_class = Disk if mode == 'fast' else Sphere
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

                    cc1 = c_class(
                        p1,
                        midpoint,
                        bond_radius,
                        color = c1,
                        **bond_style
                    )
                    cc2 = c_class(
                        midpoint,
                        p2,
                        bond_radius,
                        color = c2,
                        **bond_style
                    )
                    if objects:
                        bonds[i][j] = (( cc1, cc2 ))
                    else:
                        bonds[i][j] = (( cc1.plot(figure)[0], cc2.plot(figure)[0] ))

            if atom_style is not False:
                atoms[i] = [None] * len(geom)
                for j, stuff in enumerate(zip(colors, radii, geom)):
                    color, radius, coord = stuff
                    if 'color' not in atom_style:
                        a_sty = atom_style.copy()
                        a_sty['color'] = color
                    else:
                        a_sty = atom_style
                    sphere = s_class(coord, radius, **a_sty)
                    if objects:
                        atoms[i][j] = sphere
                    else:
                        plops = sphere.plot(figure)
                        if isinstance(plops, tuple):
                            atoms[i][j] = plops[0]
                        else:
                            atoms[i][j] = plops

        return figure, atoms, bonds

    def jupyter_viz(self):
        from McUtils.Jupyter import MoleculeGraphics

        return MoleculeGraphics(self.atoms,
                                np.ndarray.view(self.coords.convert(CartesianCoordinates3D)),
                                bonds=self.bonds
                                )
    def to_widget(self):
        return self.jupyter_viz().to_widget()
    def _ipython_display_(self):
        return self.jupyter_viz()._ipython_display_()
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
