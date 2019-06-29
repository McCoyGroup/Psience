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
        return cls(modes.T, freqs = freqs, **opts)

#TODO: make it possible to just extract certain G-matrix elements without computing the whole thing
class GMatrix:
    """Represents Wilson's G Matrix between two coordinate systems"""
    def __init__(self, system1, system2, masses, **fd_opts):
        self.sys1 = system1
        self.sys2 = system2
        self.masses = masses,
        self.opts = fd_opts
        self._jacobian = None

    @property
    def array(self):
        """Returns the numpy array form of the G-matrix

        :return:
        :rtype: np.ndarray
        """
        return self.asarray()
    @property
    def jacobian(self):
        if self._jacobian is None:
            self._jacobian = self.sys1.jacobian(self.sys2, **self.opts)
        return self._jacobian
    def asarry(self, **opts):
        if len(opts) == 0:
            jacobian = self.jacobian
        else:
            opts = dict(self.opts, **opts)
            jacobian = self.sys1.jacobian(self.sys2, **opts)
        jj = np.matmul(jacobian, jacobian.T)
        if self.masses is not None:
            jj = self.masses * jj # mass weight by the rows
        return jj


