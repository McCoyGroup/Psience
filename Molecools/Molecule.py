
class Molecule:
    """General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all

    """
    # TODO:
    #   We'll need a) a set of atoms b) a coordinate set
    #    there might be some point at which connectivity would be helpful so I guess we can include that too
    #   The coordinate set should also allow for multiconfiguration systems I think.
    #    that way we can store the
    def __init__(self, atoms, coords):
        self._ats = self._canonicalize_atoms(atoms)
    def from_zmat(self, zmat):
        pass
    def from_file(self, file):
        pass

    class Atom: # atom subclass just for providing a convenient wrapper
        def __init__(self, name, **props):
            self._kw = props
            self._name = name
class Modes:
    """A prettied up version of a Coordinerds CoordinateSet object
    Most common case of course will be working with NormalModes
    """
