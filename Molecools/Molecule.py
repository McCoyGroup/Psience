"""
Provides a simple Molecule class that we can adapt as we need
Most standard functionality should be served by OpenBabel
Uses AtomData to get properties and whatnot
"""

import os, numpy as np
from McUtils.Data import AtomData, UnitsData
from McUtils.Coordinerds import CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D
from .CoordinateSystems import MolecularCartesianCoordinateSystem, MolecularZMatrixCoordinateSystem

__all__ = [
    "Molecule",
    "MolecoolException"
]

class Molecule:
    """
    General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all
    """

    def __init__(self,
                 atoms,
                 coords,
                 bonds=None,
                 obmol=None,
                 charge=None,
                 name=None,
                 zmatrix=None,
                 dipole_surface=None,
                 potential_surface=None,
                 potential_derivatives=None,
                 source_file=None,
                 guess_bonds=True,
                 **kw
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

        coords = CoordinateSet(coords, CartesianCoordinates3D)

        # properties to be returned
        self._coords = coords
        self._sys = MolecularCartesianCoordinateSystem(self)
        self._coords = CoordinateSet(self._coords, self._sys)
        self._bonds = bonds
        self._charge = charge
        self._name = name
        self._mol = obmol
        self._kw = kw
        self.source_file = source_file
        if potential_derivatives is None or len(potential_derivatives) == 1:
            force_constants = None
        else:
            force_constants = potential_derivatives[1]
        self._fcs = force_constants
        self._pds = potential_derivatives
        self._dipoles = dipole_surface
        self._pes = potential_surface
        self._normal_modes = None
        self._zmat = zmatrix
        self._ints = None
        self.guess_bonds=guess_bonds

    def __repr__(self):
        return "{cls}('{name}', formula='{formula}', shape={shape}, coord_sys={coord_sys})".format(
            cls=type(self).__name__,
            name=self.name,
            formula=self.formula,
            shape=self.coords.shape,
            coord_sys=self.coords.system.name if isinstance(self.coords, CoordinateSet) else 'undefined'
        )

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

    @property
    def num_atoms(self):
        return len(self._ats)
    @property
    def atoms(self):
        return tuple(a["Symbol"] for a in self._ats)
    @property
    def masses(self):
        return np.array([a["Mass"] for a in self._ats])
    @property
    def bonds(self):
        if self._bonds is None and self.guess_bonds:
            self._bonds = self.prop("guessed_bonds", tol=1.05, guess_type=True)
        return self._bonds
    @property
    def coords(self):
        return self._coords
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
    def force_constants(self):
        if self._fcs is None:
            self._fcs = self.load_force_constants()
        return self._fcs
    @property
    def potential_derivatives(self):
        if self._pds is None:
            self._pds = self.load_potential_derivatives()
        return self._pds
    @property
    def potential_surface(self):
        if self._pes is None:
            self._pes = self.load_potential_surface()
        return self._pes
    @property
    def dipole_surface(self):
        if self._dipoles is None:
            self._dipoles = self.load_dipole_surface()
        return self._dipoles

    @property
    def center_of_mass(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        return self.prop('center_of_mass')
    @property
    def inertial_axes(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        return self.prop('moments_of_inertia')[1]
    @property
    def internal_coordinates(self):
        if self._ints is None and self._zmat is not None:
            self._ints = self.coords.convert(MolecularZMatrixCoordinateSystem(self, ordering=self._zmat))
        return self._ints
    @property
    def normal_modes(self):
        """

        :return:
        :rtype: VibrationalModes
        """
        if self._normal_modes is None:
            self._normal_modes = self.get_normal_modes()
        return self._normal_modes
    @normal_modes.setter
    def normal_modes(self, modes):
        """

        :return:
        :rtype: VibrationalModes
        """
        from .Vibrations import MolecularVibrations
        if not isinstance(modes, MolecularVibrations):
            raise TypeError("{}.{}: '{}' is expected to be a MolecularVibrations object".format(
                type(self).__name__,
                'normal_modes',
                modes
            ))
        self._normal_modes = modes

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
        if self._normal_modes is not None:
            self._normal_modes.molecule = self
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

    def load_force_constants(self, file=None):
        """
        Loads force constants from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """

        if file is None:
            file = self.source_file

        path, ext = os.path.splitext(file)
        ext = ext.lower()

        if ext == ".fchk":
            from McUtils.GaussianInterface import GaussianFChkReader
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(['ForceConstants'])
            return parse["ForceConstants"].array
        elif ext == ".log":
            raise NotImplementedError("{}: support for loading force constants from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
        else:
            raise NotImplementedError("{}: support for loading force constants from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
    def load_potential_derivatives(self, file=None):
        """
        Loads potential derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """

        if file is None:
            file = self.source_file
        path, ext = os.path.splitext(file)
        ext = ext.lower()

        if ext == ".fchk":
            from McUtils.GaussianInterface import GaussianFChkReader
            keys= ['Gradient', 'ForceConstants', 'ForceDerivatives']
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(keys)

            return tuple(parse[k] for k in keys)
        elif ext == ".log":
            raise NotImplementedError("{}: support for loading force constants from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
        else:
            raise NotImplementedError("{}: support for loading force constants from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
    def load_normal_modes(self, file=None):
        """
        Loads potential derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """
        from .Vibrations import MolecularNormalModes

        if file is None:
            file = self.source_file
        path, ext = os.path.splitext(file)
        ext = ext.lower()

        if ext == ".fchk":
            from McUtils.GaussianInterface import GaussianFChkReader
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(
                    ['Real atomic weights', 'VibrationalModes', 'ForceConstants', 'VibrationalData']
                )

            modes = parse["VibrationalModes"]
            freqs = parse["VibrationalData"]["Frequencies"] * UnitsData.convert("Wavenumbers", "Hartrees")
            masses = parse['Real atomic weights'] * UnitsData.convert("AtomicMassUnits", "ElectronMass")
            fcs = self.force_constants

            internal_F = np.dot(np.dot(modes, fcs), modes.T)
            raw = np.sqrt(np.diag(internal_F))
            reweight = freqs / raw
            modes = modes * reweight[:, np.newaxis]

            # mass_conv = np.sqrt(np.broadcast_to(masses[:, np.newaxis], (len(masses), 3)).flatten())
            # modes = modes * mass_conv[np.newaxis, :]

            # add in translations and rotations
            # tr_freqs, tr_vecs = self.prop("translation_rotation_eigenvectors")
            # freqs = np.concatenate([tr_freqs, freqs])
            # modes = np.concatenate([tr_vecs, modes.T], axis=1)
            modes = modes.T

            return MolecularNormalModes(self, modes, inverse=modes.T, freqs=freqs)
        elif ext == ".log":
            raise NotImplementedError("{}: support for loading normal modes from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
        else:
            raise NotImplementedError("{}: support for loading normal modes from {} files not there yet".format(
                type(self).__name__,
                ext
            ))

    def load_potential_surface(self):
        raise NotImplemented
    def load_dipole_surface(self):
        raise NotImplemented

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
    def eckart_frame(self, mol, sel=None, inverse=False):
        """
        Gets the Eckart frame(s) for the molecule
        :param mol:
        :type mol:
        :param sel: selection of atoms to use when getting the Eckart frame
        :type sel:
        :param inverse: whether to return the inverse of the rotations or not
        :type inverse: bool
        :return:
        :rtype:
        """
        return self.prop('eckart_transformation', mol, sel=sel, inverse=inverse)

    def get_embedded_molecule(self, ref=None):
        """
        Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
        a principle-axis embedded Molecule
        :return:
        :rtype: Molecule
        """

        if ref is None:
            frame = self.principle_axis_frame(inverse=True)
        else:
            frame = self.eckart_frame(ref, inverse=True)
        new = frame.apply(self)
        modes = new.normal_modes
        if modes is not None:
            modes = modes.embed(frame)
            new.normal_modes = modes
        pot_d = new.potential_derivatives
        if pot_d is not None:
            derivs = [None]*len(pot_d)
            for i, d in enumerate(pot_d):
                dim = d.ndim
                derivs[i] = frame.apply(d)
            new.potential_derivatives = derivs
        return new

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

        opts = dict({'bonds':bonds, 'mol':mol}, **opts)
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
                ['Coordinates', 'AtomicNumbers', 'Integer atomic weights']
            )
        nums = parse["AtomicNumbers"]
        wts = parse['Integer atomic weights']
        # print(nums, wts)
        mol = cls(
            [AtomData[a]["Symbol"] + str(b) for a, b in zip(nums, wts)],
            parse["Coordinates"],
            **opts
        )
        return mol

    @classmethod
    def from_file(cls, file, mode = None, **opts):
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

    def get_normal_modes(self, **kwargs):
        from .Vibrations import MolecularNormalModes, MolecularVibrations

        if self.source_file is not None:
            vibs = MolecularVibrations(self, self.load_normal_modes())
        else:
            try:
                fcs = self.force_constants
            except AttributeError:
                raise MolecoolException("{} needs '{}' bound to calculate {}".format(
                    type(self).__name__,
                    'force_constants',
                    'normal_modes'
                ))
            vibs = MolecularVibrations(self, MolecularNormalModes.from_force_constants(self, fcs, self.atoms, **kwargs))

        return vibs

    def _get_ob_attr(self, item):
        if self._mol is None:
            raise AttributeError("No pybel molecule")
        else:
            return getattr(self._mol, item)

    # def __getattr__(self, item):
    #     if '_kw' in self.__dict__:
    #         needs_raise = False
    #         try:
    #             res = self._kw[item]
    #         except KeyError:
    #             if self._mol is not None:
    #                 try:
    #                     res = self._get_ob_attr(item)
    #                 except AttributeError:
    #                     res = None
    #                     needs_raise = True
    #             else:
    #                 needs_raise = True
    #         if needs_raise:
    #             raise AttributeError("{} has no attribute '{}'".format(type(self).__name__, item))
    #         return res
    #     else:
    #         raise AttributeError("{} has no attribute '{}'".format(type(self).__name__, item))

class MolecoolException(Exception):
    pass
