from collections import OrderedDict
from McUtils.Data import AtomData

class Molecule:
    """General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all

    """
    # TODO:
    #   We'll need a) a set of atoms b) a coordinate set
    #    there might be some point at which connectivity would be helpful so I guess we can include that too
    #   The coordinate set should also allow for multiconfiguration systems I think.
    #    that way we can store the
    def __init__(self, atoms, coords, bonds = None):

        self._ats = [ AtomData[atom] for atom in atoms ]
        self._coords = coords
        self._bonds = bonds

    @property
    def atoms(self):
        return tuple(a["Symbol"] for a in self._ats)

    def from_zmat(self, zmat):
        pass
    def from_file(self, file):
        pass

    def plot(self,
             figure = None,
             bond_radius = .1, atom_radius_scaling = .25,
             atom_style = None,
             bond_style = None,
             mode = 'real'
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

        if bond_style is not False:
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


