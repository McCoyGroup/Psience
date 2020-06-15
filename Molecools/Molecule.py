"""
Provides a simple Molecule class that we can adapt as we need
Most standard functionality should be served by OpenBabel
Uses AtomData to get properties and whatnot
"""

import os
from McUtils.Data import AtomData, UnitsData
from McUtils.Coordinerds import CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D

__all__ = [
    "Molecule",
    "MolecoolException"
]

class Molecule:
    """
    General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all
    """

    PYBEL_SUPPORTED = None
    OC_SUPPORTED = None
    def __init__(self,
                 atoms,
                 coords,
                 bonds=None,
                 mol=None,
                 charge=None,
                 name=None,
                 force_constants=None,
                 dipole_surface=None,
                 potential_surface=None,
                 **kw
                 ):
        import numpy as np
        # convert "atoms" into list of atom data
        self._ats = [ AtomData[atom] if isinstance(atom, (int, np.integer, str)) else atom for atom in atoms ]

        # turn coordinates into a proper CoordinateSet object if necessary
        if not isinstance(coords, CoordinateSet):
            coords = CoordinateSet(coords)
        if not coords.system is CartesianCoordinates3D:
            sys = CartesianCoordinates3D
            coords = coords.convert(CartesianCoordinates3D)
        else:
            sys = coords.system

        # properties to be returned
        self._coords = coords
        self._sys = sys
        self._bonds = bonds
        self._charge = charge
        self._name = name
        self._mol = mol
        self._kw = kw
        self._fcs = None
        self._normal_modes = None

    def __repr__(self):
        return "{cls}('{name}', formula='{formula}', shape={shape}, coord_sys={coord_sys})".format(
            cls=type(self).__name__,
            name=self.name,
            formula=self.formula,
            shape=self.coords.shape,
            coord_sys=self.coords.system
        )

    @property
    def atoms(self):
        return tuple(a["Symbol"] for a in self._ats)
    @property
    def bonds(self):
        if self._bonds is None:
            self._bonds = self.prop("guessed_bonds", tol=1.05, guess_type=True)
        return self._bonds
    @property
    def coords(self):
        return self._coords
    @property
    def formula(self):
        return self.prop('chemical_formula')
    @property
    def multiconfig(self):
        return self.coords.multiconfig
    @property
    def force_constants(self):
        if self._fcs is None:
            self._fcs = self.load_force_constants()
        return self._fcs

    @property
    def dipole_surface(self):
        if self._fcs is None:
            self._fcs = self.load_force_constants()
        return self._fcs

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
        return cls(zmcs[0], coords, **opts)

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

    def copy(self):
        import copy
        # just use the default and don't be fancy
        return copy.copy(self)

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

    def load_force_constants(self, file=None):
        """
        Loads force constants from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """

        path, ext = os.path.splitext(file)
        ext = ext.lower()

        if ext == ".fchk":
            from McUtils.GaussianInterface import GaussianFChkReader
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(
                    ['ForceConstants']#, 'ForceDerivatives']
                )
            return parse["ForceConstants"].array
        elif ext == ".log":
            raise NotImplementedError("{}: support for loading force constants from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
            from McUtils.GaussianInterface import GaussianLogReader
        else:
            raise NotImplementedError("{}: support for loading force constants from {} files not there yet".format(
                type(self).__name__,
                ext
            ))

    #TODO: I should put pybel support into McUtils so that it can be used outside the context of a Molecule object

    @classmethod
    def from_pybel(cls, mol, **opts):
        """

        :param mol:
        :type mol: pybel.mol
        :return:
        :rtype:
        """
        import openbabel.openbabel as ob
        bonds = list(ob.OBMolBondIter(mol.OBMol))
        atoms = list(mol.atoms)

        opts = dict({'bonds':bonds, 'mol':mol}, **opts)
        return cls(atoms, [a.coords for a in atoms], **opts)

    @classmethod
    def _from_log_file(cls, file, **opts):
        from McUtils.GaussianInterface import GaussianLogReader
        with GaussianLogReader(file) as gr:
            parse = gr.parse('StandardCartesianCoordinates', num=1)
        spec, coords = parse['StandardCartesianCoordinates']
        return cls(
            [int(a[1]) for a in spec],
            coords[0],
            **opts
        )
    @classmethod
    def _from_fchk_file(cls, file, **opts):
        from McUtils.GaussianInterface import GaussianFChkReader
        with GaussianFChkReader(file) as gr:
            parse = gr.parse(
                ['Coordinates', 'AtomicNumbers', 'Integer atomic weights']
            )
        nums = parse["AtomicNumbers"]
        wts = parse['Integer atomic weights']
        # print(nums, wts)
        mol = cls(
            [AtomData[a]["Symbol"] + str(b) for a, b in zip(nums, wts)],
            UnitsData.convert("BohrRadius", "Angstroms") * parse["Coordinates"],
            **opts
        )
        # mol.force_constants = parse["ForceConstants"].array
        return mol

    @classmethod
    def from_file(cls, file, **opts):
        """In general we'll delegate to pybel except for like Fchk and Log files

        :param file:
        :type file:
        :return:
        :rtype:
        """
        import os

        path, ext = os.path.splitext(file)
        ext = ext.lower()
        opts['source_file'] = file

        format_dispatcher = {
            ".log": cls._from_log_file,
            ".fchk": cls._from_fchk_file
        }
        if ext in format_dispatcher:
            loader = format_dispatcher[ext]
            return loader(file, **opts)

        elif cls._pybel_installed():
            import openbabel.pybel as pybel
            mol = next(pybel.readfile(ext, file))
            return cls.from_pybel(mol)

        else:
            raise IOError("{} doesn't support file type {} without OpenBabel installed.".format(cls.__name__, ext))

    @classmethod
    def _pybel_installed(cls):
        if cls.PYBEL_SUPPORTED is None:
            try:
                import openbabel.pybel
            except ImportError:
                cls.PYBEL_SUPPORTED = False
            else:
                cls.PYBEL_SUPPORTED = True

        return cls.PYBEL_SUPPORTED

    @classmethod
    def _oc_installed(cls):
        if cls.OC_SUPPORTED is None:
            try:
                import openchemistry.io
            except ImportError:
                cls.OC_SUPPORTED = False
            else:
                cls.OC_SUPPORTED = True

        return cls.OC_SUPPORTED

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
        from McUtils.Plots import Graphics3D, Sphere, Cylinder, Line, Disk

        if len(geometries) == 0:
            geometries = CoordinateSet([self._coords], self._coords.system)
        elif len(geometries) == 1 and isinstance(geometries[0], CoordinateSet):
            geometries = geometries[0]
        else:
            geometries = CoordinateSet(geometries, self._coords.system)

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
                    sphere = s_class(coord, radius, color = color, **atom_style)
                    if objects:
                        atoms[i][j] = sphere
                    else:
                        plops = sphere.plot(figure)
                        if isinstance(plops, tuple):
                            atoms[i][j] = plops[0]
                        else:
                            atoms[i][j] = plops

        return figure, atoms, bonds

    def get_normal_modes(self, **kwargs):
        from .Vibrations import NormalModeCoordinates, VibrationalModes

        try:
            fcs = self.force_constants
        except AttributeError:
            raise MolecoolException("{} needs '{}' bound to calculate {}".format(
                type(self).__name__,
                'force_constants',
                'normal_modes'
            ))

        return VibrationalModes(self, NormalModeCoordinates.from_force_constants(fcs, self.atoms, **kwargs))
    @property
    def normal_modes(self):
        if self._normal_modes is None:
            self._normal_modes = self.get_normal_modes()
        return self._normal_modes

    def _get_ob_attr(self, item):
        if self._mol is not None:
            raise AttributeError("No pybel molecule")
        else:
            return getattr(self._mol, item)

    def __getattr__(self, item):
        try:
            res = self._kw[item]
        except KeyError:
            try:
                res = self._get_ob_attr(item)
            except AttributeError:
                raise AttributeError("{} has no attribute '{}'".format(type(self). __name__, item))
        return res

class MolecoolException(Exception):
    pass
