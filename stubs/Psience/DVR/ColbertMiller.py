"""
Provides implementations of the Colbert-Miller DVR types defined in
http://xbeams.chem.yale.edu/~batista/v572/ColbertMiller.pdf
"""
import numpy as np
from .BaseDVR import BaseDVR
__all__ = ['CartesianDVR', 'RingDVR', 'RadialDVR', 'PolarDVR']

class CartesianDVR(BaseDVR):
    """
    Provides the Colbert Miller DVR on the Cartesian [-inf, inf] range
    """

    def get_grid(self, domain=None, divs=None, **kw):
        """
        Provides the Colbert-Miller DVR grid for the [-inf, inf] range
        :param domain:
        :type domain:
        :param divs:
        :type divs:
        :param flavor:
        :type flavor:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """
        ...

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        """
        **LLM Docstring**

        Build the Colbert-Miller kinetic-energy matrix for the free-particle (Cartesian, `[-inf, inf]`) DVR grid, with diagonal elements `coeff * pi^2/3` and off-diagonal band elements `coeff * (-1)^(j-i) * 2/(j-i)^2`.

        :param grid: the DVR grid points
        :type grid: np.ndarray
        :param mass: the particle mass
        :type mass: float
        :param hb: the value of hbar to use
        :type hb: float
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the kinetic-energy matrix
        :rtype: np.ndarray
        """
        ...

    def real_momentum(self, grid=None, mass=None, hb=1, **kwargs):
        """
        Provides the real part of the momentum for the [-inf, inf] range
        :param grid:
        :type grid:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """
        ...

class RingDVR(BaseDVR):
    """
    Provides a DVR for working on the (0, 2Pi) range with periodicity from Colbert and Miller
    """

    def __init__(self, domain=None, **opts):
        """
        **LLM Docstring**

        Set up a Colbert-Miller DVR on the periodic `(0, 2*pi)` ring domain, defaulting the domain to `(0, 2*pi)` if not given.

        :param domain: the coordinate domain; defaults to `(0, 2*pi)`
        :type domain: tuple | None
        :param opts: extra options forwarded to the base `BaseDVR.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def get_grid(self, domain=None, divs=None, **kw):
        """
        Provides the Colbert-Miller 1D grid for the [0, 2Pi] range
        :param domain:
        :type domain:
        :param divs:
        :type divs:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """
        ...

    def get_kinetic_energy(self, grid=None, mass=1, hb=1, **kw):
        """
        Colbert-Miller kinetic energy for the [0, 2pi] range
        :param grid:
        :type grid:
        :param mass:
        :type mass:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """
        ...

    def real_momentum(self, grid=None, hb=1, **kw):
        """
        Provides the real part of the momentum for the [0, 2pi] range
        :param grid:
        :type grid:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """
        ...

class PolarDVR(BaseDVR):
    """
    Provides a DVR for working on the (0, pi) range from Colbert and Miller
    """

    def __init__(self, domain=None, **opts):
        """
        **LLM Docstring**

        Set up a Colbert-Miller DVR on the polar/angular `(0, pi)` domain, defaulting the domain to `(0, pi)` if not given.

        :param domain: the coordinate domain; defaults to `(0, pi)`
        :type domain: tuple | None
        :param opts: extra options forwarded to the base `BaseDVR.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def get_grid(self, domain=(0, np.pi), divs=None, **kwargs):
        """
        Provides the grid appropriate for the Colbert-Miller (0, Pi) range

        :param domain:
        :type domain:
        :param divs:
        :type divs:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """
        ...

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        """
        Colbert-Miller kinetic energy for the [0, pi] range
        :param grid:
        :type grid:
        :param mass:
        :type mass:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """
        ...

class RadialDVR(BaseDVR):
    """
    Provides a DVR for working on the (0, inf) range from Colbert and Miller
    """

    def get_grid(self, domain=(0, np.pi), divs=None, **kwargs):
        """
        Provides the grid appropriate for the Colbert-Miller (0, Pi) range

        :param domain:
        :type domain:
        :param divs:
        :type divs:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """
        ...

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        """
        Colbert-Miller kinetic energy for the (0, inf) range

        :param grid:
        :type grid:
        :param mass:
        :type mass:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """
        ...