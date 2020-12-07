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

def _get_best_axes(first_pos, axes):
    """
    Determine the best pair of inertial axes so that we don't get large-scale breakdowns from the choice of embedding

    :param first_pos:
    :type first_pos:
    :param axes:
    :type axes:
    :return:
    :rtype:
    """
    fp_norm = np.linalg.norm(first_pos)
    if fp_norm > 1.0e-10:  # not chilling at the origin...
        first_pos = first_pos / fp_norm
        # check if it lies along an axis or is perpendicular to an axis
        a_proj = np.dot(first_pos, axes[0])
        b_proj = np.dot(first_pos, axes[1])
        c_proj = np.dot(first_pos, axes[2])
        if np.abs(b_proj) < .05: # lies in the A/C plane
            if np.abs(a_proj) > .95:
                axes = axes[1:]
                ax_names = ["B", "C"]
            else:
                axes = axes[:2]
                ax_names = ["A", "B"]
        elif np.abs(c_proj) < .05: # lies in the A/B plane
            if np.abs(a_proj) > .95:
                axes = axes[1:]
                ax_names = ["B", "C"]
            else:
                axes = axes[(0, 2),]
                ax_names = ["A", "C"]
        elif np.abs(a_proj) < .05:  # lies in the B/C plane
            if np.abs(b_proj) > .95:
                axes = axes[(0, 2),]
                ax_names = ["A", "C"]
            else:
                axes = axes[:2]
                ax_names = ["A", "B"]
        else: # not in any of the planes so no issues
            axes = axes[:2]
            ax_names = ["A", "B"]

    else:
        axes = axes[:2]
        ax_names = ["A", "B"]
    return axes, ax_names

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
        from .Molecule import Molecule
        self.molecule = molecule #type: Molecule
        if converter_options is None:
            converter_options = opts
            opts = {}
        nats = len(molecule.atoms)
        super().__init__(converter_options=converter_options, dimension=(nats, 3), coordinate_shape=(nats, 3), opts=opts)
        self.set_embedding()
    @property
    def origins(self):
        return self.converter_options['origins']
    @property
    def axes(self):
        return self.converter_options['axes']

    def set_embedding(self):
        molecule = self.molecule
        com = molecule.center_of_mass
        axes = molecule.inertial_axes
        converter_options = self.converter_options
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
        axes, ax_names = _get_best_axes(first_pos, axes)

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['axes_labels'] = ax_names
        converter_options['molecule'] = molecule

    def jacobian(self, *args, **kwargs):
        try:
            remb = self.converter_options['reembed']
        except KeyError:
            remb = None

        try:
            self.converter_options['reembed'] = True
            return super().jacobian(*args, **kwargs)
        finally:
            if remb is not None:
                self.converter_options['reembed'] = remb

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
        self.molecule = molecule #type: Molecule
        nats = len(self.molecule.atoms)
        if converter_options is None:
            converter_options = opts
            opts = {}
        super().__init__(converter_options=converter_options, dimension=(nats, 3), opts=opts)
        self.set_embedding()

    def set_embedding(self):
        """
        Sets up the embedding options...
        :return:
        :rtype:
        """
        molecule = self.molecule
        com = molecule.center_of_mass
        axes = molecule.inertial_axes
        converter_options = self.converter_options
        if 'ordering' in converter_options:
            ordering = np.array(converter_options['ordering'], dtype=int)
            ordering[0, 1] = -3; ordering[0, 2] = -2; ordering[0, 3] = -1
            ordering[1, 2] = -1; ordering[1, 3] = -2
            ordering[2, 3] = -2
            converter_options['ordering'] = ordering
            first = ordering[0, 0]
        else:
            first = 0

        first_pos = molecule.coords[first]
        axes, ax_names = _get_best_axes(first_pos, axes)

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['axes_labels'] = ax_names
        converter_options['molecule'] = molecule

class MolecularCartesianToZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """
    types = (MolecularCartesianCoordinateSystem, MolecularZMatrixCoordinateSystem)
    def convert(self, coords, molecule=None, origins=None, axes=None, ordering=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        :param coords:
        :type coords: CoordinateSet
        :param molecule:
        :type molecule:
        :param origins:
        :type origins:
        :param axes:
        :type axes:
        :param ordering:
        :type ordering:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
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

    def convert_many(self, coords, molecule=None, origins=None, axes=None, ordering=None, reembed=False, **kwargs):
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
            derivs = opts['derivs']
            reshaped_derivs = [None] * len(derivs)
            for i, v in enumerate(derivs):
                # drop all terms relating to the embedding of the embedding
                sel_1 = (
                        (..., ) +
                        (slice(2, None, None), slice(None, None, None)) * (i + 1)
                        + (slice(3, None, None), slice(None, None, None))
                )
                reshaped_derivs[i] = v[sel_1]
            opts['derivs'] = reshaped_derivs
            # raise Exception(derivs.shape)
        # reembed = False
        if reembed:
            from .Molecule import Molecule
            from .Transformations import MolecularTransformation

            if molecule is None:
                raise ValueError("can't reembed without a reference structure")
            # we calculate the Eckart embedding for the generated coords
            if carts.ndim != 3:
                raise ValueError("reembedding only currently implemented for stacks of structures")
            for i, crd in enumerate(carts): # lazy hack right now
                structures = Molecule(molecule.atoms, crd)
                ek = structures.eckart_frame(molecule) #type: MolecularTransformation
                # we chain on a rotation out of the inertial frame to keep stuff consistent
                # ek = MolecularTransformation(pax.transformation_function.transform, ek)
                carts[i] = ek.apply(crd)
                rot = ek.transformation_function.transform # ignore the shift when applying to derivs
                if 'derivs' in opts:
                    for v in opts['derivs']:
                        new = np.tensordot(v[i], rot, axes=[-1, 0])
                        v[i] = new
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



