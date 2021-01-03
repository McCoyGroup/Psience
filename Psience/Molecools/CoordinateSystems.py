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

    if axes.ndim > 2:
        axes = axes[..., (0, 1), :]
        ax_choice = (0, 1)
        ax_names = ["A", "B"]
    else:
        fp_norm = np.linalg.norm(first_pos)
        if fp_norm > 1.0e-10:  # not chilling at the origin...
            first_pos = first_pos / fp_norm
            # check if it lies along an axis or is perpendicular to an axis
            a_proj = np.dot(first_pos, axes[0])
            b_proj = np.dot(first_pos, axes[1])
            c_proj = np.dot(first_pos, axes[2])
            if np.abs(b_proj) < .05: # lies in the A/C plane
                if np.abs(a_proj) > .95:
                    ax_choice = (1, 2)
                    ax_names = ["B", "C"]
                else:
                    ax_choice = (0, 1)
                    ax_names = ["A", "B"]
            elif np.abs(c_proj) < .05: # lies in the A/B plane
                if np.abs(a_proj) > .95:
                    ax_choice = (1, 2)
                    ax_names = ["B", "C"]
                else:
                    ax_choice = (0, 2)
                    ax_names = ["A", "C"]
            elif np.abs(a_proj) < .05:  # lies in the B/C plane
                if np.abs(b_proj) > .95:
                    ax_choice = (0, 2)
                    ax_names = ["A", "C"]
                else:
                    ax_choice = (0, 1)
                    ax_names = ["A", "B"]
            else: # not in any of the planes so no issues
                ax_choice = (0, 1)
                ax_names = ["A", "B"]

        else:
            ax_choice = (0, 1)
            ax_names = ["A", "B"]
        axes = axes[ax_choice,]
    return axes, ax_names, ax_choice

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

    def pre_convert(self, system):
        self.set_embedding()

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
        axes, ax_names, ax_choice = _get_best_axes(first_pos, axes)

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['axes_labels'] = ax_names
        converter_options['axes_choice'] = ax_choice
        converter_options['molecule'] = molecule

    def jacobian(self, *args, **kwargs):
        try:
            remb = self.converter_options['reembed']
        except KeyError:
            remb = None

        try:
            self.converter_options['reembed'] = True if remb is None else remb
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

    def pre_convert(self, syystem):
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
        axes, ax_names, ax_choice = _get_best_axes(first_pos, axes)

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['axes_labels'] = ax_names
        converter_options['axes_choice'] = ax_choice
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

        zmcs, opts = self.convert_many(np.array([coords]), molecule=molecule, origins=origins, axes=axes, ordering=ordering, **kwargs)
        zmcs = zmcs[0]

        if 'derivs' in opts:
            derivs = opts['derivs']
            reshaped_derivs = [None] * len(derivs)
            for i, v in enumerate(derivs):
                reshaped_derivs[i] = v[0]
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

    def convert_many(self, coords, molecule=None, origins=None, axes=None, ordering=None,
                     reembed=False, axes_choice=None, return_derivs=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, attempting to preserve the embedding
        """
        from .Molecule import Molecule

        n_sys = coords.shape[0]
        n_coords = coords.shape[1]
        n_atoms = len(molecule.atoms)
        if n_coords != n_atoms + 2:
            # means we already added the embedding
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

        refuse_derivs = reembed and coords.squeeze().ndim != 2
        res = CoordinateSet(coords, ZMatrixCoordinates).convert(CartesianCoordinates3D,
                                                                        ordering=ordering,
                                                                        origins=origins,
                                                                        axes=axes,
                                                                        return_derivs=(return_derivs and not refuse_derivs),
                                                                        **kwargs)

        if isinstance(res, tuple):
            carts, opts = res
        else:
            carts = res
            opts = res.converter_options
        carts = carts[..., 3:, :]

        if reembed:
            if molecule is None:
                raise ValueError("can't reembed without a reference structure")
            reembed = not (carts.squeeze().ndim == 2 and np.allclose(molecule.coords, carts, atol=1.0e-5)) # agree to like a ten thousandth of an angstrom
            if reembed:
                if not return_derivs:
                    carts = molecule.embed_coords(carts)
                else:
                    print("garb")
                    inert_coords, coord_coms, coord_axes = Molecule(molecule.atoms, carts).principle_axis_data
                    print("parp")
                    if axes_choice is None:
                        axes_choice = (0, 1)
                    guh = self.convert_many(coords,
                                             origins=coord_coms,
                                             axes=coord_axes[:, axes_choice],
                                             molecule=molecule,
                                             reembed=False,
                                             ordering=ordering,
                                             return_derivs=return_derivs,
                                             axes_choice=axes_choice,
                                             **kwargs
                                             )
                    # print("...farts")
                    return guh

        opts['origins'] = origins
        opts['axes'] = axes
        if ordering is not None:
            opts['ordering'] = ordering[3:] - 3
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



