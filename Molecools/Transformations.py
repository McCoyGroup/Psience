"""
Defines a MolecularTransformation class that uses Coordinerds to describe a physical transformation of a molecule
Then it layers some common transformations on top of that
"""
from .Molecule import Molecule
from McUtils.Coordinerds.CoordinateTransformations import CoordinateTransform

class MolecularTransformation(CoordinateTransform):
    def apply(self, mol):
        """

        :param mol:
        :type mol: Molecule | np.ndarray
        :return:
        :rtype:
        """
        if isinstance(mol, Molecule):
            new_coords = super().apply(mol.coords)
            new = mol.copy()
            new._coords = new_coords
        else:
            new = super().apply(mol)
        return new
