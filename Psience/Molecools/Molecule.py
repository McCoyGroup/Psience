"""
Provides a simple Molecule class that we can adapt as we need
Most standard functionality should be served by OpenBabel
Uses AtomData to get properties and whatnot
"""

import os, numpy as np

from McUtils.Data import AtomData, UnitsData
from McUtils.Coordinerds import CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D
import McUtils.Numputils as nput

from .MoleculeInterface import *
from .CoordinateSystems import MolecularCartesianCoordinateSystem, MolecularZMatrixCoordinateSystem
from .Properties import *
from .Transformations import *

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
                 zmatrix=None,
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
        :param name: The internal coordinate Z-matrix specification for the molecule
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
        if zmatrix is not None:
            zmatrix = np.asanyarray(zmatrix).astype(int)
            if zmatrix.shape[1] != 4:
                raise ValueError("can't understand Z-matrix {}".format(zmatrix))
        self._zmat = zmatrix
        self._ints = None

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

    @property
    def zmatrix(self):
        """
        :return:
        :rtype:
        """
        return self._zmat
    @zmatrix.setter
    def zmatrix(self, zmatrix):
        """
        :return:
        :rtype:
        """
        #TODO: add some validation
        zmatrix = np.asanyarray(zmatrix).astype(int)
        if zmatrix.shape[1] != 4:
            raise ValueError("can't understand Z-matrix {}".format(zmatrix))
        self._zmat = zmatrix
    @property
    def internal_coordinates(self):
        if self._ints is None and self._zmat is not None:
            zms = MolecularZMatrixCoordinateSystem(self, ordering=self._zmat)
            # print(zms)
            # print(zms, self.coords, self.coords.system.converter(zms))
            self._ints = self.coords.convert(zms)
        return self._ints
    @property
    def g_matrix(self):
        """
        Returns the molecular g-matrix for the system
        :return:
        :rtype:
        """
        if self.internal_coordinates is None:
            raise ValueError("need an internal coordinate Z-matrix to calculate the G-matrix")
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

    def embed_coords(self, crds, sel=None, planar_ref_tolerance=None):
        """
        Embeds coords in the Eckart frame using `self` as a reference
        :param crds:
        :type crds:
        :return:
        :rtype:
        """

        return self.prop('eckart_embedded_coords', crds, sel=sel, planar_ref_tolerance=planar_ref_tolerance)
    def get_embedding_data(self, crds, sel=None):
        """
        Gets the necessary data to embed crds in the Eckart frame using `self` as a reference
        :param crds:
        :type crds:
        :return:
        :rtype: tuple[np.ndarray, tuple[np.ndarray], tuple[np.ndarray]]
        """
        return self.prop('eckart_embedding_data', crds, sel=sel)
    def get_embedded_molecule(self,
                              ref=None,
                              embed_properties=True
                              ):
        """
        Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
        a principle-axis embedded Molecule
        :return:
        :rtype: Molecule
        """

        if ref is None:
            frame = self.principle_axis_frame(inverse=False)
        else:
            frame = self.eckart_frame(ref, inverse=False)
        # self.normal_modes.modes
        new = frame.apply(self)
        if embed_properties:
            # inv_frame = frame
            # if ref is None:
            #     inv_frame = self.principle_axis_frame(inverse=True)
            # else:
            #     inv_frame = self.eckart_frame(ref, inverse=True)
            new.normal_modes = new.normal_modes.apply_transformation(frame)
            new.potential_surface = new.potential_surface.apply_transformation(frame)
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
             figure = None,
             bond_radius = .1, atom_radius_scaling = .25,
             atom_style = None,
             bond_style = None,
             mode = 'fast',
             objects = False,
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
