"""
Defines a MolecularTransformation class that uses Coordinerds to describe a physical transformation of a molecule
Then it layers some common transformations on top of that
"""
from .MoleculeInterface import AbstractMolecule
from McUtils.Coordinerds import CoordinateSet
from McUtils.Numputils import GeometricTransformation

__all__ = ["MolecularTransformation"]
__reload_hook__ = [".MoleculeInterface"]

class MolecularTransformation(GeometricTransformation):
    def apply(self, mol, shift=True):
        """

        :param mol:
        :type mol: AbstractMolecule | np.ndarray | CoordinateTransform
        :return:
        :rtype:
        """
        if isinstance(mol, AbstractMolecule):
            new = mol.copy()
            new_coords = super().apply(new.coords)
            if isinstance(mol.coords, CoordinateSet):
                new.coords = CoordinateSet(new_coords, mol.coords.system)
            else:
                new.coords = CoordinateSet(new_coords)
        elif isinstance(mol, GeometricTransformation):
            new = type(self)(self, mol)#, shift=shift)
        else:
            new = super().apply(mol, shift=shift)
        return new
    def __call__(self, mol, shift=True):
        """
        **LLM Docstring**

        Alias for `apply`. Lets a `MolecularTransformation` instance be used directly as a callable to transform `mol`.

        :param mol: the object to transform, either an `AbstractMolecule`, a raw coordinate array, or another `GeometricTransformation`
        :type mol: AbstractMolecule | np.ndarray | CoordinateTransform
        :param shift: whether to apply the translation part of the transformation (passed through to `apply`)
        :type shift: bool
        :return: the transformed object, same handling as `apply`
        :rtype: AbstractMolecule | np.ndarray | GeometricTransformation
        """
        return self.apply(mol, shift=shift)
