"""
Defines the fundamental interface we expect all molecules to support
without any concrete implementation details so that we can separate property
calculation from explicit molecule instantiation
"""

import abc, typing, numpy as np
from McUtils.Coordinerds import CoordinateSet

__all__ = ["AbstractMolecule"]

class PotentialSurfaceType(typing.Protocol):
    @property
    def surface(self):
        raise NotImplementedError('abstract interface')
    @property
    def derivatives(self):
        raise NotImplementedError('abstract interface')
    def load(self):
        raise NotImplementedError('abstract interface')
    def update(self, val):
        raise NotImplementedError('abstract interface')
    def insert_atoms(self, atoms, coords, where):
        raise NotImplementedError('abstract interface')
    def delete_atoms(self, where):
        raise NotImplementedError('abstract interface')
    def apply_transformation(self, frame):
        raise NotImplementedError('abstract interface')

class DipoleSurfaceType(typing.Protocol):
    @property
    def surface(self):
        raise NotImplementedError('abstract interface')
    @property
    def derivatives(self):
        raise NotImplementedError('abstract interface')
    def load(self):
        raise NotImplementedError('abstract interface')
    def update(self, val):
        raise NotImplementedError('abstract interface')
    def insert_atoms(self, atoms, coords, where):
        raise NotImplementedError('abstract interface')
    def delete_atoms(self, where):
        raise NotImplementedError('abstract interface')
    def apply_transformation(self, frame):
        raise NotImplementedError('abstract interface')

class NormalModesType(typing.Protocol):
    @property
    def modes(self):
        raise NotImplementedError('abstract interface')
    def load(self):
        raise NotImplementedError('abstract interface')
    def update(self, val):
        raise NotImplementedError('abstract interface')
    def insert_atoms(self, atoms, coords, where):
        raise NotImplementedError('abstract interface')
    def delete_atoms(self, where):
        raise NotImplementedError('abstract interface')
    def apply_transformation(self, frame):
        raise NotImplementedError('abstract interface')

class AbstractMolecule(metaclass=abc.ABCMeta):
    """
    A molecule base class to provide type constraints
    """

    @property
    @abc.abstractmethod
    def name(self) -> typing.AnyStr:
        """
        Provides the name of the molecule
        """
        raise NotImplementedError('abstract class')

    @property
    @abc.abstractmethod
    def metadata(self) -> typing.Dict:
        """
        Provides the metadta the molecule
        """
        raise NotImplementedError('abstract class')

    @property
    @abc.abstractmethod
    def atoms(self) -> typing.Tuple[typing.AnyStr]:
        """
        Provides the atom names for the molecule
        """
        raise NotImplementedError('abstract class')

    @property
    @abc.abstractmethod
    def masses(self) -> np.ndarray:
        """
        Provides the atomic masses for the molecule
        """
        raise NotImplementedError('abstract class')

    @property
    @abc.abstractmethod
    def coords(self) -> CoordinateSet:
        """
        Provides the atomic coordinates for the molecule
        """
        raise NotImplementedError('abstract class')
    @coords.setter
    def coords(self, new):
        """
        Sets the atomic coordinates for the molecule
        """
        raise NotImplementedError('abstract class')
    @property
    @abc.abstractmethod
    def internal_coordinates(self) -> CoordinateSet:
        """
        Provides the internal coordinates for the molecule
        """
        raise NotImplementedError('abstract class')
    @property
    @abc.abstractmethod
    def center_of_mass(self) -> np.ndarray:
        """
        Provides the center of mass for the molecule
        """
        raise NotImplementedError('abstract class')
    @property
    @abc.abstractmethod
    def inertial_axes(self) -> np.ndarray:
        """
        Provides the principal axes for the molecule
        """
        raise NotImplementedError('abstract class')


    @property
    @abc.abstractmethod
    def bonds(self) -> typing.Iterable[typing.Tuple]:
        """
        Provides the bonds for the molecule
        """
        raise NotImplementedError('abstract class')

    @property
    @abc.abstractmethod
    def potential_surface(self) -> PotentialSurfaceType:
        """
        Provides the potential surface for the molecule
        """
        raise NotImplementedError('abstract class')

    @property
    @abc.abstractmethod
    def dipole_surface(self) -> DipoleSurfaceType:
        """
        Provides the dipole surface for the molecule
        """
        raise NotImplementedError('abstract class')
    @property
    @abc.abstractmethod
    def normal_modes(self) -> NormalModesType:
        """
        Provides the normal modes for the molecule
        """
        raise NotImplementedError('abstract class')

    @property
    @abc.abstractmethod
    def source_file(self) -> str:
        """
        Provides the source file for the molecule
        """
        raise NotImplementedError('abstract class')

    @abc.abstractmethod
    def copy(self) -> 'AbstractMolecule':
        """
        Copies the molecule
        """
        raise NotImplementedError('abstract class')
