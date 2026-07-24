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
        """
        **LLM Docstring**

        Protocol stub for a property that provides access to the underlying potential energy surface data. Not implemented here; concrete types satisfying this protocol must provide their own `surface` property.

        :return: the potential surface
        :rtype: object
        """
        raise NotImplementedError('abstract interface')
    @property
    def derivatives(self):
        """
        **LLM Docstring**

        Protocol stub for a property that provides the derivatives of the potential surface (e.g. gradients, force constants). Not implemented here.

        :return: the potential surface derivatives
        :rtype: object
        """
        raise NotImplementedError('abstract interface')
    def load(self):
        """
        **LLM Docstring**

        Protocol stub for loading the potential surface data. Not implemented here; concrete types must define how the surface is loaded (e.g. from a file or calculation).

        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def update(self, val):
        """
        **LLM Docstring**

        Protocol stub for updating the stored potential surface with a new value. Not implemented here.

        :param val: the new value to update the surface with
        :type val: object
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def insert_atoms(self, atoms, coords, where):
        """
        **LLM Docstring**

        Protocol stub describing how the potential surface should be updated when atoms are inserted into the underlying molecule. Not implemented here.

        :param atoms: the atoms being inserted
        :type atoms: Iterable
        :param coords: coordinates of the inserted atoms
        :type coords: np.ndarray
        :param where: index/indices at which the atoms are inserted
        :type where: int or Iterable[int]
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def delete_atoms(self, where):
        """
        **LLM Docstring**

        Protocol stub describing how the potential surface should be updated when atoms are removed from the underlying molecule. Not implemented here.

        :param where: index/indices of the atoms being removed
        :type where: int or Iterable[int]
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def apply_transformation(self, frame):
        """
        **LLM Docstring**

        Protocol stub describing how the potential surface should respond to a coordinate frame transformation being applied to the molecule. Not implemented here.

        :param frame: the transformation/frame to apply
        :type frame: object
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')

class DipoleSurfaceType(typing.Protocol):
    @property
    def surface(self):
        """
        **LLM Docstring**

        Protocol stub for a property that provides access to the underlying dipole surface data. Not implemented here; concrete types satisfying this protocol must provide their own `surface` property.

        :return: the dipole surface
        :rtype: object
        """
        raise NotImplementedError('abstract interface')
    @property
    def derivatives(self):
        """
        **LLM Docstring**

        Protocol stub for a property that provides the derivatives of the dipole surface. Not implemented here.

        :return: the dipole surface derivatives
        :rtype: object
        """
        raise NotImplementedError('abstract interface')
    def load(self):
        """
        **LLM Docstring**

        Protocol stub for loading the dipole surface data. Not implemented here.

        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def update(self, val):
        """
        **LLM Docstring**

        Protocol stub for updating the stored dipole surface with a new value. Not implemented here.

        :param val: the new value to update the surface with
        :type val: object
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def insert_atoms(self, atoms, coords, where):
        """
        **LLM Docstring**

        Protocol stub describing how the dipole surface should be updated when atoms are inserted into the underlying molecule. Not implemented here.

        :param atoms: the atoms being inserted
        :type atoms: Iterable
        :param coords: coordinates of the inserted atoms
        :type coords: np.ndarray
        :param where: index/indices at which the atoms are inserted
        :type where: int or Iterable[int]
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def delete_atoms(self, where):
        """
        **LLM Docstring**

        Protocol stub describing how the dipole surface should be updated when atoms are removed from the underlying molecule. Not implemented here.

        :param where: index/indices of the atoms being removed
        :type where: int or Iterable[int]
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def apply_transformation(self, frame):
        """
        **LLM Docstring**

        Protocol stub describing how the dipole surface should respond to a coordinate frame transformation being applied to the molecule. Not implemented here.

        :param frame: the transformation/frame to apply
        :type frame: object
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')

class NormalModesType(typing.Protocol):
    @property
    def modes(self):
        """
        **LLM Docstring**

        Protocol stub for a property that provides access to the underlying normal-mode data. Not implemented here; concrete types satisfying this protocol must provide their own `modes` property.

        :return: the normal modes
        :rtype: object
        """
        raise NotImplementedError('abstract interface')
    def load(self):
        """
        **LLM Docstring**

        Protocol stub for loading the normal-mode data. Not implemented here.

        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def update(self, val):
        """
        **LLM Docstring**

        Protocol stub for updating the stored normal modes with a new value. Not implemented here.

        :param val: the new value to update the normal modes with
        :type val: object
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def insert_atoms(self, atoms, coords, where):
        """
        **LLM Docstring**

        Protocol stub describing how the normal modes should be updated when atoms are inserted into the underlying molecule. Not implemented here.

        :param atoms: the atoms being inserted
        :type atoms: Iterable
        :param coords: coordinates of the inserted atoms
        :type coords: np.ndarray
        :param where: index/indices at which the atoms are inserted
        :type where: int or Iterable[int]
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def delete_atoms(self, where):
        """
        **LLM Docstring**

        Protocol stub describing how the normal modes should be updated when atoms are removed from the underlying molecule. Not implemented here.

        :param where: index/indices of the atoms being removed
        :type where: int or Iterable[int]
        :return: None
        :rtype: None
        """
        raise NotImplementedError('abstract interface')
    def apply_transformation(self, frame):
        """
        **LLM Docstring**

        Protocol stub describing how the normal modes should respond to a coordinate frame transformation being applied to the molecule. Not implemented here.

        :param frame: the transformation/frame to apply
        :type frame: object
        :return: None
        :rtype: None
        """
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
        Provides the masses for the molecule
        """
        raise NotImplementedError('abstract class')

    @property
    @abc.abstractmethod
    def atomic_masses(self) -> np.ndarray:
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
    def copy(self) -> 'typing.Self':
        """
        Copies the molecule
        """
        raise NotImplementedError('abstract class')

    @abc.abstractmethod
    def take_submolecule(self, pos) -> 'typing.Self':
        """
        Takes a submolecule
        """
        raise NotImplementedError('abstract class')
