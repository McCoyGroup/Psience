"""
Defines useful extended internal coordinate frames
"""



import numpy as np
import McUtils.Numputils as nput
from McUtils.Coordinerds import (
    ZMatrixCoordinateSystem, CartesianCoordinateSystem, CoordinateSystemConverter,
    ZMatrixCoordinates, CartesianCoordinates3D, CoordinateSet, CoordinateSystemConverters
)
from .MoleculeInterface import AbstractMolecule

__all__ = [
    "MolecularZMatrixCoordinateSystem",
    "MolecularCartesianCoordinateSystem"
]

__reload_hook__ = [".MoleculeInterface"]

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
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        self.molecule = molecule
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

    def jacobian(self,
                 *args,
                 reembed=None,
                 strip_dummies=None,
                 converter_options=None,
                 **kwargs
                 ):
        if converter_options is None:
            converter_options = {}
        merged_convert_options = dict(self.converter_options, **converter_options)
        try:
            remb = merged_convert_options['reembed'] if reembed is None else reembed
        except KeyError:
            remb = None

        try:
            strip_dummies = merged_convert_options['strip_dummies'] if strip_dummies is None else strip_dummies
        except KeyError:
            strip_dummies = False

        if strip_dummies:
            dummies = self.molecule.dummy_positions
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(self.molecule.num_atoms),
                dummies
            )

        try:
            self.converter_options['reembed'] = True if remb is None else remb
            jacs = super().jacobian(*args, converter_options=converter_options, **kwargs)
            raw_jacs = []
            for j in jacs:
                ext_dim = j.ndim - 2
                shp = sum(
                    ((j.shape[i] // 3, 3) for i in range(ext_dim)),
                    ()
                ) +  j.shape[-2:]
                j = j.reshape(shp)
                if dummies is not None:
                    for i in range(ext_dim):
                        j = np.take(j, main_excludes, axis=2*i)

                #     j.shape[:i]
                #     + (j.shape[i] // 3, 3)
                #     + j.shape[i+1:]
                # )
                raw_jacs.append(j)
            jacs = raw_jacs
            return jacs
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
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        self.molecule = molecule #type: AbstractMolecule
        nats = len(self.molecule.atoms)
        if converter_options is None:
            converter_options = opts
            opts = {}
        super().__init__(converter_options=converter_options, dimension=(nats, 3), opts=opts)

    def pre_convert(self, system):
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

    def jacobian(self,
                 coords,
                 system,
                 strip_dummies=None,
                 converter_options=None,
                 analytic_deriv_order=None,
                 **kwargs
                 ):
        if converter_options is None:
            converter_options = {}
        merged_convert_options = dict(self.converter_options, **converter_options)
        try:
            strip_dummies = merged_convert_options['strip_dummies'] if strip_dummies is None else strip_dummies
        except KeyError:
            strip_dummies = False

        try:
            analytic_deriv_order = merged_convert_options['analytic_deriv_order'] if analytic_deriv_order is None else analytic_deriv_order
        except KeyError:
            analytic_deriv_order = 0

        if strip_dummies:
            dummies = self.molecule.dummy_positions
            if len(dummies) == 0:
                dummies = None
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(self.molecule.num_atoms),
                dummies
            )
        else:
            main_excludes = None

        jacs = super().jacobian(coords, system, analytic_deriv_order=analytic_deriv_order, converter_options=converter_options, **kwargs)
        raw_jacs = []
        for n,j in enumerate(jacs): # this expects a full filling of the jacobians which maybe I need to not expect...
            baseline = 2*analytic_deriv_order + len(coords.shape)
            ext_dim = j.ndim - baseline
            shp = sum(
                ((j.shape[i] // 3, 3) for i in range(ext_dim)),
                ()
            ) + j.shape[-baseline:]
            j = j.reshape(shp)
            if dummies is not None:
                for i in range(ext_dim):
                    j = np.take(j, main_excludes, axis=2*i)
                for i in range(analytic_deriv_order):
                    j = np.take(j, main_excludes, axis=-2*(i+2))

            if len(coords.shape) > 2:
                j = np.moveaxis(j, -3, 0)

            raw_jacs.append(j)
        jacs = raw_jacs
        return jacs

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

        zmcs, opts = self.convert_many(np.array([coords]),
                                       molecule=molecule, origins=origins, axes=axes, ordering=ordering, **kwargs)
        zmcs = zmcs[0]

        if 'derivs' in opts:
            derivs = opts['derivs']
            reshaped_derivs = [None] * len(derivs)
            for i, v in enumerate(derivs):
                reshaped_derivs[i] = v[0]
            opts['derivs'] = reshaped_derivs

        return zmcs, opts

    def convert_many(self, coords,
                     molecule=None,
                     origins=None, axes=None,
                     ordering=None,
                     strip_embedding=True,
                     strip_dummies=False,
                     **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding

        :param coords: coordinates in Cartesians to convert
        :type coords: np.ndarray
        :param molecule:
        :type molecule: AbstractMolecule
        :param origins: the origin for each individual structure
        :type origins: np.ndarray
        :param axes: the axes for each structure
        :type axes: np.ndarray
        :param ordering: the Z-matrix ordering spec
        :type ordering:
        :param strip_embedding: whether to strip the embedding coordinates
        :type strip_embedding:
        :param strip_dummies: whether to strip all dummy coordinates
        :type strip_dummies:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """

        n_sys = coords.shape[0]
        n_coords = coords.shape[1]
        n_atoms = len(molecule.atoms)

        # we add three dummy atoms at the origins and along the axes before doing the conversion
        if origins.ndim == 1:
            origins = np.broadcast_to(origins[np.newaxis, np.newaxis], (n_sys, 1, 3))
        elif origins.ndim == 2:
            origins = origins[:, np.newaxis, :]
        if axes.ndim == 2:
            axes = np.broadcast_to(axes[np.newaxis], (n_sys, 2, 3))
        if origins.shape[0] != n_sys:
            if n_sys % origins.shape[0] != 0:
                raise ValueError("inconsistent shapes; origins shape {} but coords shape {}".format(
                    origins.shape,
                    coords.shape
                ))
            num_coords = n_sys // origins.shape[0]
            origins = np.broadcast_to(origins[:, np.newaxis, :, :], (origins.shape[0], num_coords) + origins.shape[1:])
            origins = origins.reshape((n_sys,) + origins.shape[2:])
        if axes.shape[0] != n_sys:
            if n_sys % axes.shape[0] != 0:
                raise ValueError("inconsistent shapes; axes shape {} but coords shape {}".format(
                    axes.shape,
                    coords.shape
                ))
            num_coords = n_sys // axes.shape[0]
            axes = np.broadcast_to(axes[:, np.newaxis, :, :], (axes.shape[0], num_coords) + axes.shape[1:])
            axes = axes.reshape((n_sys,) + axes.shape[2:])
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
        # zmcs = zmcs[:, 2:]
        if strip_dummies:
            dummies = [0, 1, 2] + [x+3 for x in molecule.dummy_positions] # add on axes
        elif strip_embedding:
            dummies = [0, 1, 2]
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(len(molecule.atoms) + 3),
                dummies
            )
            sub_excludes = main_excludes - 1 # drop one fewer terms to drop I think...
            if 'derivs' in opts:
                derivs = opts['derivs']
                reshaped_derivs = [None] * len(derivs)
                deriv_excludes = np.arange(3, len(molecule.atoms) + 3)
                for i, v in enumerate(derivs):
                    # drop all terms relating to the embedding of the embedding
                    start_dim = v.ndim - 2*(i+2)
                    for j in range(start_dim, v.ndim-2, 2):
                        v = np.take(v, deriv_excludes, axis=j)
                    v = np.take(v, sub_excludes, axis=-2)
                    reshaped_derivs[i] = v

                opts['derivs'] = reshaped_derivs

            zmcs = zmcs[..., sub_excludes, :]
            # raise Exception(derivs.shape)
        return zmcs, opts

MolecularCartesianToZMatrixConverter = MolecularCartesianToZMatrixConverter()
MolecularCartesianToZMatrixConverter.register(CoordinateSystemConverters)

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
                     reembed=False, axes_choice=None, return_derivs=None,
                     strip_dummies=False,
                     strip_embedding=True,
                     planar_ref_tolerance=None,
                     **kwargs):
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

        if reembed:
            if molecule is None:
                raise ValueError("can't reembed without a reference structure")
            embed_carts = carts[..., 3:, :]
            reembed = not (
                    carts.squeeze().ndim == 2 and
                    np.allclose(molecule.coords, embed_carts, atol=1.0e-5)
            ) # agree to like a ten thousandth of an angstrom
            if reembed:
                if not return_derivs:
                    embed_carts = molecule.embed_coords(embed_carts, planar_ref_tolerance=planar_ref_tolerance)
                    carts = np.concatenate([
                        carts[..., :3, :],
                        embed_carts
                        ],
                    axis=-2
                    )
                else:
                    inert_coords, coord_coms, coord_axes = Molecule(molecule.atoms, embed_carts).principle_axis_data
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
                    return guh

        opts['origins'] = origins
        opts['axes'] = axes
        if ordering is not None:
            opts['ordering'] = ordering[3:] - 3
        if strip_dummies:
            # raise Exception("wwwwaaaaaaaaat")
            dummies = [0, 1, 2] + [x + 3 for x in molecule.dummy_positions]  # add on axes
        elif strip_embedding:
            dummies = [0, 1, 2]
        else:
            dummies = None
        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(len(molecule.atoms) + 3),
                dummies
            )
            sub_excludes = main_excludes - 1  # drop one fewer terms to drop I think...
            if 'derivs' in opts:
                derivs = opts['derivs']
                reshaped_derivs = [None] * len(derivs)
                deriv_excludes = np.arange(3, len(molecule.atoms) + 3)
                for i, v in enumerate(derivs):
                    # drop all terms relating to the embedding of the embedding
                    start_dim = v.ndim - i
                    for j in range(start_dim, v.ndim, 2):
                        v = np.take(v, deriv_excludes, axis=j)
                    v = np.take(v, sub_excludes, axis=-2)
                    reshaped_derivs[i] = v
                opts['derivs'] = reshaped_derivs

            carts = carts[..., main_excludes, :]

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


