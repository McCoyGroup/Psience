"""
Defines a MolecularTransformation class that uses Coordinerds to describe a physical transformation of a molecule
Then it layers some common transformations on top of that
"""
from .MoleculeInterface import AbstractMolecule
from McUtils.Coordinerds import CoordinateSet, CoordinateTransform

__all__ = ["MolecularTransformation"]
__reload_hook__ = [".MoleculeInterface"]

class MolecularTransformation(CoordinateTransform):
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
        elif isinstance(mol, CoordinateTransform):
            new = type(self)(self, mol)#, shift=shift)
        else:
            new = super().apply(mol, shift=shift)
        return new
    def __call__(self, mol, shift=True):
        return self.apply(mol, shift=shift)
