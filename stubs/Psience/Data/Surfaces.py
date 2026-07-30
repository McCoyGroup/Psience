"""
Provides concrete tools for dealing with two of the most useful types of surfaces we have
"""
import numpy as np
from collections import namedtuple
from McUtils.GaussianInterface import GaussianLogReader, GaussianFChkReader
from McUtils.Zachary import Surface, MultiSurface, InterpolatedSurface, TaylorSeriesSurface
import McUtils.Numputils as nput
__all__ = ['DipoleSurface', 'PotentialSurface']

class DipoleSurface(MultiSurface):
    """
    Provides a unified interface to working with dipole surfaces.
    Currently basically no fancier than a regular surface (although with convenient loading functions), but dipole-specific
    stuff could come
    """

    def __init__(self, mu_x, mu_y, mu_z):
        """

        :param mu_x: X-component of dipole moment
        :type mu_x: Surface
        :param mu_y: Y-component of dipole moment
        :type mu_y: Surface
        :param mu_z: Z-component of dipole moment
        :type mu_z: Surface
        """
        ...

    @property
    def center(self):
        ...

    @property
    def ref(self):
        ...

    @property
    def expansion_tensors(self):
        ...

    @staticmethod
    def get_log_values(log_file, keys=('StandardCartesianCoordinates', 'DipoleMoments')):
        ...

    @classmethod
    def from_log_file(cls, log_file, coord_transf, keys=('StandardCartesianCoordinates', 'DipoleMoments'), tol=0.001, **opts):
        """
        Loads dipoles from a Gaussian log file and builds a dipole surface by interpolating.
        Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
        to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
        Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions

        :param log_file: a Gaussian log file to pull from
        :type log_file: str
        :return:
        :rtype:
        """
        ...

    @classmethod
    def from_fchk_file(cls, fchk_file, **opts):
        """
        Loads dipoles from a Gaussian formatted checkpoint file and builds a dipole surface via a linear approximation

        :param fchk_file: a Gaussian fchk file to pull from
        :type log_file: str
        :return:
        :rtype:
        """
        ...

    @classmethod
    def from_derivatives(cls, expansion, center=None, **opts):
        ...

    @classmethod
    def from_mol(cls, mol, expansion=None, center=None, transforms=None, use_internals=True, **opts):
        ...

    def __call__(self, gridpoints, **opts):
        """
        Explicitly overrides the Surface-level evaluation because we know the Taylor surface needs us to flatten our gridpoints

        :param gridpoints:
        :type gridpoints:
        :param opts:
        :type opts:
        :return:
        :rtype:
        """
        ...

class PotentialSurface(Surface):
    """
    A potential surface structure to go along with the DipoleSurface.
    Provides convenient access to potential data + a unified interface to things like energy minimization
    """

    @staticmethod
    def get_log_values(log_file, keys=('StandardCartesianCoordinates', 'ScanEnergies')):
        ...

    @classmethod
    def from_log_file(cls, log_file, coord_transf, keys=('StandardCartesianCoordinates', 'ScanEnergies'), tol=0.001, **opts):
        """
        Loads dipoles from a Gaussian log file and builds a potential surface by interpolating.
        Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
        to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
        Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions.

        :param log_file: a Gaussian log file to pull from
        :type log_file: str
        :return:
        :rtype:
        """
        ...

    @classmethod
    def from_fchk_file(cls, fchk_file, ref=None, **opts):
        """
        Loads potential from a Gaussian formatted checkpoint file and builds a potential surface via a quartic approximation

        :param fchk_file: a Gaussian fchk file to pull from
        :type log_file: str
        :return:
        :rtype:
        """
        ...

    @classmethod
    def from_mol(cls, mol, expansion=None, center=None, transforms=None, transformed_derivatives=False, use_internals=True, **opts):
        ...

    @classmethod
    def from_derivatives(cls, expansion, center=None, ref=None, transforms=None, transformed_derivatives=False, **opts):
        ...

    def __call__(self, gridpoints, **opts):
        """
        Explicitly overrides the Surface-level evaluation because we know the Taylor surface needs us to flatten our gridpoints

        :param gridpoints:
        :type gridpoints:
        :param opts:
        :type opts:
        :return:
        :rtype:
        """
        ...