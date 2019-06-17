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
        self._ats = OrderedDict( (( at["Symbol"], at) for at in self._ats) )
        self._coords = coords
        self._bonds = bonds

    def from_zmat(self, zmat):
        pass
    def from_file(self, file):
        pass

