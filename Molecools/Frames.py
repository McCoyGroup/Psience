"""
Defines useful extended internal coordinate frames
"""

__all__ = [
    "MolecularZMatrixCoordinateSystem",
    "MolecularCartesianCoordinateSystem"
]

import numpy as np

from McUtils.Coordinerds import (
    ZMatrixCoordinateSystem, CartesianCoordinateSystem, CoordinateSystemConverter,
    ZMatrixCoordinates, CartesianCoordinates3D, CoordinateSet
)

class MolecularZMatrixCoordinateSystem(ZMatrixCoordinateSystem):
    """
    Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding
    """
    name = "MolecularZMatrixCoordinates"
    def __init__(self, molecule, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: Molecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        from .Molecule import Molecule
        molecule = molecule #type: Molecule

        com = molecule.center_of_mass
        axes = molecule.inertial_axes
        if converter_options is None:
            converter_options = opts
            opts = {}
        converter_options['origin'] = com
        converter_options['axes'] = axes
        converter_options['molecule'] = molecule
        nats = len(molecule.atoms)
        super().__init__(converter_options=converter_options, dimension=(nats, 3), coordinate_shape=(nats, 3), opts=opts)
    @property
    def origin(self):
        return self.converter_options['origin']
    @property
    def axes(self):
        return self.converter_options['axes']

class MolecularCartesianCoordinateSystem(CartesianCoordinateSystem):
    """
    Mirrors the standard Cartesian coordinate system in _almost_ all regards, but forces an embedding
    """
    def __init__(self, molecule, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: Molecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """

        from .Molecule import Molecule
        molecule = molecule #type: Molecule
        com = molecule.center_of_mass
        axes = molecule.inertial_axes
        if converter_options is None:
            converter_options = opts
            opts = {}
        converter_options['origin'] = com
        converter_options['axes'] = axes
        converter_options['molecule'] = molecule
        nats = len(molecule.atoms)
        super().__init__(converter_options=converter_options, dimension=(nats, 3), opts=opts)

class MolecularCartesianToMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """
    types = (MolecularCartesianCoordinateSystem, MolecularZMatrixCoordinateSystem)
    def convert(self, coords, molecule=None, origins=None, axes=None, ordering=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        n_coords = len(coords)
        n_atoms = len(molecule.atoms)
        if n_coords != n_atoms:
            raise ValueError('Embedding unclear when num_coords ({}) < num_atoms ({})'.format(
                n_coords,
                n_atoms
            ))

        carts, opts = ZMatrixCoordinates(coords).convert(CartesianCoordinates3D, ordering=ordering, **kwargs)
        opts['origins'] = origins
        opts['axes'] = axes
        carts = carts[3:]
        return carts, opts

class MolecularZMatrixToCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """
    types = (MolecularZMatrixCoordinateSystem, MolecularCartesianCoordinateSystem)
    def convert(self, coords, molecule=None, origins=None, axes=None, ordering=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        n_coords = len(coords)
        n_atoms = len(molecule.atoms)
        if n_coords != n_atoms:
            # we add three dummy atoms at the origin and along the axes before doing the conversion
            coords = np.concatenate([origins, axes, coords])
            if ordering is not None:
                ordering = np.ndarray(ordering, dtype=int)
                ordering[0, 1] = 0; ordering[0, 2] = -1; ordering[1, 2] = 0
                ordering = ordering + 2

        carts, opts = ZMatrixCoordinates(coords).convert(CartesianCoordinates3D, ordering=ordering, **kwargs)
        opts['origins'] = origins
        opts['axes'] = axes
        carts = carts[3:]
        return carts, opts



