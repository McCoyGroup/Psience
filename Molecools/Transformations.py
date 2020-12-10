"""
Defines a MolecularTransformation class that uses Coordinerds to describe a physical transformation of a molecule
Then it layers some common transformations on top of that
"""
from .Molecule import Molecule
from McUtils.Coordinerds import CoordinateSet, CoordinateTransform

class MolecularTransformation(CoordinateTransform):
    def apply(self, mol, shift=True):
        """

        :param mol:
        :type mol: Molecule | np.ndarray | CoordinateTransform
        :return:
        :rtype:
        """
        if isinstance(mol, Molecule):
            new_coords = super().apply(mol.coords)
            new = mol.copy()
            if isinstance(mol.coords, CoordinateSet):
                new._coords = CoordinateSet(new_coords, mol.coords.system)
            else:
                new._coords = CoordinateSet(new_coords)
        elif isinstance(mol, CoordinateTransform):
            new = super().__call__(mol, shift=shift)
        else:
            new = super().apply(mol, shift=shift)
        return new
    def __call__(self, mol, shift=True):
        return self.apply(mol, shift=shift)
