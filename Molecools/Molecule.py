"""
Provides a simple Molecule class that we can adapt as we need
Most standard functionality should be served by OpenBabel
Uses AtomData to get properties and whatnot
"""

from McUtils.Data import AtomData

class Molecule:
    """General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all

    """
    # TODO:
    #   We'll need a) a set of atoms b) a coordinate set
    #    there might be some point at which connectivity would be helpful so I guess we can include that too
    #   The coordinate set should also allow for multiconfiguration systems I think
    #    that way we can store many copies of a molecule at once

    PYBEL_SUPPORTED = None
    def __init__(self, atoms, coords, bonds = None, mol = None, **kw):
        import numpy as np
        self._ats = [ AtomData[atom] if isinstance(atom, (int, np.integer, str)) else atom for atom in atoms ]
        self._coords = coords
        self._bonds = bonds
        self._mol = mol
        self._kw = kw

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
        from ..Coordinerds import CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D

        if isinstance(zmat, str):
            from McUtils.Parsers.ParserUtils import pull_zmat
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
                parse = gr.parse(['Coordinates', 'AtomicNumbers'])
            return cls(
                parse["AtomicNumbers"],
                parse["Coordinates"],
                **opts
            )
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

    def plot(self,
             figure = None,
             bond_radius = .1, atom_radius_scaling = .25,
             atom_style = None,
             bond_style = None,
             mode = 'fast'
             ):
        from McUtils.Plots import Graphics3D, Sphere, Cylinder, Line, Disk

        if figure is None:
            figure = Graphics3D()

        colors = [ at["IconColor"] for at in self._ats ]
        radii = [ atom_radius_scaling * at["IconRadius"] for at in self._ats ]

        bonds = []
        atoms = []

        if atom_style is None:
            atom_style = {}
        if bond_style is None:
            bond_style = {}

        if bond_style is not False and self._bonds is not None:
            c_class = Line if mode == 'fast' else Cylinder
            for b in self._bonds:
                atom1 = b[0]
                atom2 = b[1]
                # i'm not supporting double or triple bonds for not because they're just too much of a pain...
                # in Mathematica I have some clever code for finding the midpoint between the vectors to draw to
                # I'm not gonna do that for now because it's not worth it for the quick-and-dirty stuff I want to do
                p1 = self._coords[atom1]
                p2 = self._coords[atom2]
                midpoint = (p2 - p1)/2 + p1
                c1 = colors[atom1]
                c2 = colors[atom2]


                c = c_class(
                    p1,
                    midpoint,
                    bond_radius,
                    color = c1,
                    **bond_style
                )
                bonds.append(c.plot(figure))
                c = c_class(
                    midpoint,
                    p2,
                    bond_radius,
                    color = c2,
                    **bond_style
                )
                bonds.append(c.plot(figure))

        if atom_style is not False:
            s_class = Disk if mode == 'fast' else Sphere
            for color, radius, coord in zip(colors, radii, self._coords):
                sphere = s_class(coord, radius, color = color, **atom_style)
                atoms.append(sphere.plot(figure))

        return figure, atoms, bonds


