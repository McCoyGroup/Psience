from ..Coordinerds.CoordinateSystems import CoordinateSystem, CartesianCoordinates3D
from McUtils.Data import AtomData
import numpy as np

class NormalModes(CoordinateSystem):
    """A prettied up version of a Coordinerds CoordinateSystem object
    Has function for generating, though, too
    """
    def __init__(self, coeffs, basis = CartesianCoordinates3D, name = None, freqs = None):
        super().__init__(
            matrix = coeffs,
            name = name,
            basis = basis,
            dimension = basis.dimension
        )
        self.freqs = freqs

    @classmethod
    def from_force_constants(cls, fcs, atoms = None, **opts):
        if atoms is None:
            # intended for the internal coordinate case where it's hard to define masses in general
            weighted_fcs = fcs
        else:
            masses = np.array([ AtomData[a, "Mass"] if isinstance(a, str) else a for a in atoms ])
            masses = np.broadcast_to(masses, (len(masses), 3)).T.flatten()
            g = 1/np.sqrt(masses)
            g = np.outer(g, g)
            weighted_fcs = g * fcs

        freqs, modes = np.linalg.eigh(weighted_fcs)
        cls(modes.T, freqs = freqs, **opts)
