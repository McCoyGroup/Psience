"""
Defines useful extended internal coordinate frames
"""

__all__ = [
    "MolecularZMatrixCoordinateSystem",
    "MolecularCartesianCoordinateSystem"
]

import numpy as np
import McUtils.Numputils as nput

from McUtils.Coordinerds import (
    ZMatrixCoordinateSystem, CartesianCoordinateSystem, CoordinateSystemConverter,
    ZMatrixCoordinates, CartesianCoordinates3D, CoordinateSet
)

class MolecularZMatrixCoordinateSystem(ZMatrixCoordinateSystem):
    """
    Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding
    """
    name = "MolecularZMatrix"
    def __init__(self, molecule, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: Molecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        # from .Molecule import Molecule
        # molecule = molecule #type: Molecule

        com = molecule.center_of_mass
        axes = molecule.inertial_axes
        if converter_options is None:
            converter_options = opts
            opts = {}
        if 'ordering' in converter_options:
            ordering = np.array(converter_options['ordering'], dtype=int)
            ordering[0, 1] = -3; ordering[0, 2] = -1; ordering[0, 3] = -2
            ordering[1, 2] = -1; ordering[1, 3] = -2
            ordering[2, 3] = -2
            converter_options['ordering'] = ordering
            first = ordering[0, 0]
        else:
            first = 0

        first_pos = molecule.coords[first]
        fp_norm = np.linalg.norm(first_pos)
        if fp_norm > 1.0e-10:  # not chilling at the origin...
            first_pos = first_pos / fp_norm
            if np.abs(np.dot(first_pos, axes[0])) == 1:
                axes = axes[1:]
            elif np.abs(np.dot(first_pos, axes[1])) == 1:
                axes = axes[(0, 2)]
        else:
            axes = axes[:2]

        # np.random.seed(100)
        # com = np.random.rand(3)
        # ax_1 = np.random.rand(3)
        # ax_1 = ax_1 / np.linalg.norm(ax_1)
        # ax_2 = np.random.rand(3)
        # ax_2 = ax_2 / np.linalg.norm(ax_2)
        # ax_2 = np.cross(ax_1, ax_2)
        # ax_2 = ax_2 / np.linalg.norm(ax_2)
        # axes = np.array([ax_1, ax_2])

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['molecule'] = molecule
        nats = len(molecule.atoms)
        super().__init__(converter_options=converter_options, dimension=(nats, 3), coordinate_shape=(nats, 3), opts=opts)
    @property
    def origins(self):
        return self.converter_options['origins']
    @property
    def axes(self):
        return self.converter_options['axes']

class MolecularCartesianCoordinateSystem(CartesianCoordinateSystem):
    """
    Mirrors the standard Cartesian coordinate system in _almost_ all regards, but forces an embedding
    """
    name= "MolecularCartesians"
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
        if 'ordering' in converter_options:
            ordering = np.array(converter_options['ordering'], dtype=int)
            ordering[0, 1] = -3; ordering[0, 2] = -2; ordering[0, 3] = -1
            ordering[1, 2] = -1; ordering[1, 3] = -2
            ordering[2, 3] = -2
            converter_options['ordering'] = ordering
            first = ordering[0, 0]
        else:
            first = 0

        # make sure first atom isn't along one of the axes...
        first_pos = molecule.coords[first]
        fp_norm = np.linalg.norm(first_pos)
        if fp_norm > 1.0e-10: # not chilling at the origin...
            first_pos = first_pos / fp_norm
            if np.abs(np.dot(first_pos, axes[0])) > .95:
                axes = axes[1:]
            elif np.abs(np.dot(first_pos, axes[1])) > .95:
                axes = axes[(0, 2)]
        else:
            axes = axes[:2]

        # np.random.seed(100)
        # com = np.random.rand(3)
        # ax_1 = np.random.rand(3)
        # ax_1 = ax_1 / np.linalg.norm(ax_1)
        # ax_2 = np.random.rand(3)
        # ax_2 = ax_2 / np.linalg.norm(ax_2)
        # ax_2 = np.cross(ax_1, ax_2)
        # ax_2 = ax_2 / np.linalg.norm(ax_2)
        # axes = np.array([ax_1, ax_2])

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['molecule'] = molecule
        nats = len(molecule.atoms)
        super().__init__(converter_options=converter_options, dimension=(nats, 3), opts=opts)

class MolecularCartesianToZMatrixConverter(CoordinateSystemConverter):
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

        # we add three dummy atoms at the origins and along the axes before doing the conversion
        if origins.ndim == 1:
            origins = origins[np.newaxis]
        coords = np.concatenate([origins, origins+axes, coords], axis=0)
        # print(coords)
        if ordering is not None:
            ordering = np.array(ordering, dtype=int)
            ordering[0, 1] = -3; ordering[0, 2] = -1; ordering[0, 3] = -2
            ordering[1, 2] = -2; ordering[1, 3] = -1
            ordering[2, 3] = -1
            ordering = ordering + 3
            ordering = np.concatenate([
                [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]],
                ordering
            ])

            # raise Exception(ordering)

        res = CoordinateSet(coords, CartesianCoordinates3D).convert(ZMatrixCoordinates,
                                                                    ordering=ordering,
                                                                    **kwargs
                                                                    )

        if isinstance(res, tuple):
            zmcs, opts = res
        else:
            zmcs = res
            opts=res.converter_options


        # opts['origins'] = origins
        # opts['axes'] = axes
        opts['ordering'] = opts['ordering'][3:] - 3
        zmcs = zmcs[2:]
        if 'derivs' in opts:
            derivs = opts['derivs']
            reshaped_derivs = [None] * len(derivs)
            for i, v in enumerate(derivs):
                # drop all terms relating to the embedding of the embedding
                sel_1 = (
                        (slice(3, None, None), slice(None, None, None)) * (i + 1)
                        + (slice(2, None, None), slice(None, None, None))
                )
                reshaped_derivs[i] = v[sel_1]

                # print("> watwat", v.shape, reshaped_derivs[i].shape)


            opts['derivs'] = reshaped_derivs
        return zmcs, opts
    def convert_many(self, coords, molecule=None, origins=None, axes=None, ordering=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """

        n_sys = coords.shape[0]
        n_coords = coords.shape[1]
        n_atoms = len(molecule.atoms)

        # we add three dummy atoms at the origins and along the axes before doing the conversion
        if origins.ndim == 1:
            origins = np.broadcast_to(origins[np.newaxis, np.newaxis], (n_sys, 1, 3))
        if axes.ndim == 2:
            axes = np.broadcast_to(axes[np.newaxis], (n_sys, 2, 3))
        coords = np.concatenate([origins, origins+axes, coords], axis=1)
        if ordering is not None:
            ordering = np.array(ordering, dtype=int)
            ordering[0, 1] = -3; ordering[0, 2] = -2; ordering[0, 3] = -1
            ordering[1, 2] = -2; ordering[1, 3] = -1
            ordering[2, 3] = -1
            ordering = ordering + 3
            ordering = np.concatenate([ [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], ordering])
            # print("...?", ordering)

        res = CoordinateSet(coords, CartesianCoordinates3D).convert(ZMatrixCoordinates,
                                                                    ordering=ordering,
                                                                    origins=origins,
                                                                    axes=axes,
                                                                    **kwargs
                                                                    )

        if isinstance(res, tuple):
            zmcs, opts = res
        else:
            zmcs = res
            opts=res.converter_options
        # opts['origins'] = origins
        # opts['axes'] = axes
        opts['ordering'] = opts['ordering'][3:] - 3
        zmcs = zmcs[:, 2:]
        if 'derivs' in opts:
            derivs = opts['derivs']
            reshaped_derivs = [None] * len(derivs)
            for i, v in enumerate(derivs):
                # drop all terms relating to the embedding of the embedding
                sel_1 = (
                        (..., ) +
                        (slice(3, None, None), slice(None, None, None)) * (i + 1)
                        + (slice(2, None, None), slice(None, None, None))
                )
                reshaped_derivs[i] = v[sel_1]

                # print("> wat", v.shape, reshaped_derivs[i].shape)
            opts['derivs'] = reshaped_derivs
            # raise Exception(derivs.shape)
        return zmcs, opts

MolecularCartesianToMatrixConverter = MolecularCartesianToZMatrixConverter()
MolecularCartesianToMatrixConverter.register()

class MolecularCartesianToRegularCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """
    types = (MolecularCartesianCoordinateSystem, CartesianCoordinateSystem)
    def convert(self, coords, **kw):
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        return coords, kwargs
MolecularCartesianToRegularCartesianConverter = MolecularCartesianToRegularCartesianConverter()
MolecularCartesianToRegularCartesianConverter.register()

class MolecularZMatrixToCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """
    types = (MolecularZMatrixCoordinateSystem, MolecularCartesianCoordinateSystem)
    def convert(self, coords, **kw):
        total_points, opts = self.convert_many(coords[np.newaxis], **kw)
        return total_points[0], opts

    def convert_many(self, coords, molecule=None, origins=None, axes=None, ordering=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """

        n_sys = coords.shape[0]
        n_coords = coords.shape[1]
        n_atoms = len(molecule.atoms)
        if n_coords != n_atoms:
            raise ValueError('Embedding unclear when num_coords ({}) < num_atoms ({})'.format(
                n_coords,
                n_atoms
            ))

        x_ax = axes[..., 0, :]
        y_ax = axes[..., 1, :]
        extra_norms0 = nput.vec_norms(x_ax)
        extra_norms1 = nput.vec_norms(y_ax)
        extra_angles, _ = nput.vec_angles(x_ax, y_ax)
        extra_coords = np.zeros((n_sys, 2, 3))
        extra_coords[..., 0, 0] = extra_norms0
        extra_coords[..., 1, 0] = extra_norms1
        extra_coords[..., 1, 1] = extra_angles

        coords = np.concatenate([extra_coords, coords], axis=-2)
        if ordering is not None:
            ordering = np.array(ordering, dtype=int)
            ordering = ordering + 3
            ordering = np.concatenate([ [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], ordering])
        # raise Exception([ordering, coords])
        res = CoordinateSet(coords, ZMatrixCoordinates).convert(CartesianCoordinates3D,
                                                                        ordering=ordering,
                                                                        origins=origins,
                                                                        axes=axes,
                                                                        **kwargs)

        if isinstance(res, tuple):
            carts, opts = res
        else:
            carts = res
            opts = res.converter_options

        opts['origins'] = origins
        opts['axes'] = axes
        if ordering is not None:
            opts['ordering'] = ordering[3:] - 3
        carts = carts[..., 3:, :]
        if 'derivs' in opts:
            derivs = opts['derivs'][..., 2:, :, 3:, :]
            opts['derivs'] = derivs
        return carts, opts

MolecularZMatrixToCartesianConverter = MolecularZMatrixToCartesianConverter()
MolecularZMatrixToCartesianConverter.register()

class MolecularZMatrixToRegularZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """
    types = (MolecularZMatrixCoordinateSystem, ZMatrixCoordinateSystem)
    def convert(self, coords, **kw):
        return coords, kw

    def convert_many(self, coords, **kwargs):
        return coords, kwargs
MolecularZMatrixToRegularZMatrixConverter = MolecularZMatrixToRegularZMatrixConverter()
MolecularZMatrixToRegularZMatrixConverter.register()



