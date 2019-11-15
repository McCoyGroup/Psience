"""
Provides a simple Molecule class that we can adapt as we need
Most standard functionality should be served by OpenBabel
Uses AtomData to get properties and whatnot
"""

from McUtils.Data import AtomData
from ..Coordinerds import CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D

__all__ = [
    "Molecule"
]

class Molecule:
    """General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all

    """
    # TODO:
    #   We'll need a) a set of atoms b) a coordinate set
    #    there might be some point at which connectivity would be helpful so I guess we can include that too
    #   The coordinate set should also allow for multiconfiguration systems I think
    #    that way we can store many copies of a molecule at once

    PYBEL_SUPPORTED = None
    OC_SUPPORTED = None
    def __init__(self, atoms, coords,
                 bonds = None,
                 mol = None,
                 **kw
                 ):
        import numpy as np
        self._ats = [ AtomData[atom] if isinstance(atom, (int, np.integer, str)) else atom for atom in atoms ]
        if not isinstance(coords, CoordinateSet):
            coords = CoordinateSet(coords)
        self._coords = coords
        self._bonds = bonds
        self._mol = mol
        self._kw = kw
        self._normal_modes = None

    @property
    def atoms(self):
        return tuple(a["Symbol"] for a in self._ats)

    @classmethod
    def from_zmat(cls, zmat, **opts):
        """Little z-matrix importer

        :param zmat:
        :type zmat: str | tuple
        :return:
        :rtype: Molecule
        """

        if isinstance(zmat, str):
            from McUtils.Parsers.RegexPatterns import pull_zmat
            zmcs = pull_zmat(zmat)
        else:
            zmcs = zmat

        coords = CoordinateSet([ zmcs[1] ], ZMatrixCoordinates).convert(CartesianCoordinates3D)

        return cls(zmcs[0], coords, **opts)

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
    def from_file(cls, file, **opts):
        """In general we'll delegate to pybel except for like Fchk and Log files

        :param file:
        :type file:
        :return:
        :rtype:
        """
        import os
        path, ext = os.path.splitext(file)

        if ext == '.log':
            from McUtils.GaussianInterface import GaussianLogReader
            with GaussianLogReader(file) as gr:
                parse = gr.parse('CartesianCoordinates', num=1)
            spec, coords = parse['CartesianCoordinates']
            return cls(
                [ int(a[1]) for a in spec ],
                coords[0],
                **opts
            )
        elif ext == '.fchk':
            from McUtils.GaussianInterface import GaussianFChkReader
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(['Coordinates', 'AtomicNumbers', "ForceConstants", "ForceDerivatives"])
            mol = cls(
                parse["AtomicNumbers"],
                parse["Coordinates"],
                **opts
            )
            mol.force_constants = parse["ForceConstants"].array
            return mol

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
             mode = 'normal',
             objects = False,
             **plot_ops
             ):
        from McUtils.Plots import Graphics3D, Sphere, Cylinder, Line, Disk

        if len(geometries) == 0:
            geometries = (self._coords, )

        if figure is None:
            figure = Graphics3D(backend="VTK", **plot_ops)

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
            if bond_style is not False and self._bonds is not None:

                bonds[i] = [None] * len(self._bonds)
                for j, b in enumerate(self._bonds):
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
                        atoms[i][j] = sphere.plot(figure)[0]

        return figure, atoms, bonds

    @property
    def normal_modes(self):
        if self._normal_modes is None:
            from .Vibrations import NormalModeCoordiantes, VibrationalModes

            try:
                fcs = self.force_constants
            except AttributeError:
                raise MolecoolException("{} needs '{}' bound to calculate {}".format(
                    type(self).__name__,
                    'force_constants',
                    'normal_modes'
                ))

            self._normal_modes = VibrationalModes(self, NormalModeCoordiantes.from_force_constants(fcs, self.atoms))
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
