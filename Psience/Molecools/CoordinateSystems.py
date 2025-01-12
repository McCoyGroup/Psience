"""
Defines useful extended internal coordinate frames
"""


import numpy as np
import McUtils.Numputils as nput
from McUtils.Parallelizers import Parallelizer
from McUtils.Coordinerds import (
    ZMatrixCoordinateSystem, CartesianCoordinateSystem, CoordinateSystemConverter,
    ZMatrixCoordinates, CartesianCoordinates3D, CoordinateSet, CompositeCoordinateSystem,
    GenericInternalCoordinateSystem, GenericInternalCoordinates,
    CartesianToGICSystemConverter, GICSystemToCartesianConverter
)

from ..Modes import RedundantCoordinateGenerator

from .Properties import StructuralProperties
# from .MoleculeInterface import AbstractMolecule

__all__ = [
    "MolecularEmbedding",
    "MolecularZMatrixCoordinateSystem",
    "MolecularCartesianCoordinateSystem"
]

__reload_hook__ = [".MoleculeInterface"]

class MolecularEmbedding:

    def __init__(self,
                 masses,
                 coords,
                 internals
                 ):

        self._coords = CoordinateSet(coords, MolecularCartesianCoordinateSystem(masses, coords))
        MolecularCartesianToRegularCartesianConverter(self.coords.system).register()
        RegularCartesianToMolecularCartesianConverter(self.coords.system).register()
        if isinstance(internals, CoordinateSet):
            self._int_spec = None
            self._ints = internals
        else:
            self._int_spec = self.canonicalize_internal_coordinate_spec(internals)
            self._ints = None
        self._jacobians = self._get_jacobian_storage()
        # self._int_spec = internals
        self._frame = None
        self._tr_modes = None

    @property
    def coords(self):
        return self._coords
    @coords.setter
    def coords(self, coords):
        if hasattr(coords, "system"):
            sys = coords.system
        else:
            sys = self.coords.system
            coords = CoordinateSet(coords, sys)
        if sys is not self.coords.system:
            MolecularCartesianToRegularCartesianConverter(sys).register()
            RegularCartesianToMolecularCartesianConverter(sys).register()
        self._jacobians = self._get_jacobian_storage()
        self._frame = None
        self._coords = coords
    @property
    def masses(self):
        return self._coords.system.masses

    @staticmethod
    def _wrap_conv(f):
        if f is None:
            return f

        def wrapped(*args, **kwargs):
            vals = f(*args, **kwargs)
            if not isinstance(vals, np.ndarray):
                vals, opts = vals
            else:
                opts = {}
            opts = dict(kwargs, **opts)
            return vals, opts

        return wrapped

    @classmethod
    def canonicalize_internal_coordinate_spec(cls, spec):
        if spec is not None:
            if hasattr(spec, 'items'):
                specs = spec.get('specs')
                if specs is None:
                    try:
                        zmatrix = spec['zmatrix']
                    except KeyError:
                        zmatrix = None
                    else:
                        zmatrix = np.asanyarray(zmatrix).astype(int)
                        if zmatrix.shape[1] != 4:
                            raise ValueError("can't understand Z-matrix {}".format(zmatrix))
                    spec['zmatrix'] = zmatrix
                    spec['conversion'] = cls._wrap_conv(spec.get('conversion'))
                    spec['inverse'] = cls._wrap_conv(spec.get('inverse'))
                    converter_options = spec.get('converter_options')
                    if converter_options is None:
                        converter_options = {}
                    if 'embedding_coords' not in converter_options:
                        if spec['zmatrix'] is not None:
                            converter_options['embedding_coords'] = MolecularZMatrixCoordinateSystem.embedding_coords
                    if 'jacobian_prep' not in converter_options:
                        if spec['zmatrix'] is not None:
                            converter_options['jacobian_prep'] = ZMatrixCoordinates.jacobian_prep_coordinates
                    spec['converter_options'] = converter_options
                else:
                    spec['zmatrix'] = None
                    spec['conversion'] = cls._wrap_conv(spec.get('conversion'))
                    spec['inverse'] = cls._wrap_conv(spec.get('inverse'))
            elif callable(spec):
                zmatrix = None
                specs = None
                conversion = cls._wrap_conv(spec)
                spec = {'specs': specs, 'zmatrix': zmatrix, 'conversion': conversion, 'inverse': None, 'converter_options': {}}
            else:
                if any(len(s) < 4 or hasattr(s, 'items') for s in spec):
                    conversion = None
                    specs = spec
                    zmatrix = None
                else:
                    conversion = None
                    specs = None
                    zmatrix = np.asanyarray(spec).astype(int)
                    if zmatrix.shape[1] != 4:
                        raise ValueError("can't understand Z-matrix {}".format(zmatrix))
                spec = {'specs': specs, 'zmatrix': zmatrix, 'conversion': conversion, 'inverse': None, 'converter_options': {}}
        return spec

    @property
    def internals(self):
        if self._int_spec is not None:
            return self._int_spec

    @internals.setter
    def internals(self, internals):
        self._int_spec = self.canonicalize_internal_coordinate_spec(internals)
        self._ints = None

    @property
    def zmatrix(self):
        if self._int_spec is not None:
            return self._int_spec['zmatrix']

    @zmatrix.setter
    def zmatrix(self, zmat):
        if zmat is not None:
            zmat = np.asanyarray(zmat).astype(int)
            if zmat.shape[1] != 4:
                raise ValueError("can't understand Z-matrix {}".format(zmat))
        if self._int_spec is None:
            self._int_spec = self.canonicalize_internal_coordinate_spec(zmat)
        else:
            self._int_spec['zmatrix'] = zmat
        self._ints = None

    @property
    def internal_coordinates(self):
        if self._ints is None and (
                self._int_spec is not None
                and any(self._int_spec[k] is not None for k in {'zmatrix', 'conversion', 'specs'})
        ):
            coords = self.coords
            if self._int_spec['specs'] is not None:
                ints = MolecularGenericInternalCoordinateSystem(self.masses, coords,
                                                                specs=self._int_spec['specs'],
                                                                redundant=self._int_spec.get('redundant', False),
                                                                untransformed_coordinates=self._int_spec.get('untransformed_coordinates'),
                                                                relocalize=self._int_spec.get('relocalize', False)
                                                                )
                MolecularCartesianToGICConverter(coords.system, ints).register()
                MolecularGICToCartesianConverter(ints, coords.system).register()
                coords = coords.convert(ints)
            elif self._int_spec['zmatrix'] is not None:
                zms = MolecularZMatrixCoordinateSystem(self.masses, coords, ordering=self._int_spec['zmatrix'])
                MolecularCartesianToZMatrixConverter(coords.system, zms).register()
                MolecularZMatrixToCartesianConverter(zms, coords.system).register()
                MolecularZMatrixToRegularZMatrixConverter(zms).register()
                RegularZMatrixToMolecularZMatrixConverter(zms).register()
                coords = coords.convert(zms)
            if self._int_spec['conversion'] is not None:
                conv = CompositeCoordinateSystem.register(
                    coords.system,
                    self._int_spec['conversion'],
                    inverse_conversion=self._int_spec['inverse'],
                    **self._int_spec['converter_options']
                )
                coords = coords.convert(conv)
            # print(zms)
            # print(zms, self.coords, self.coords.system.converter(zms))
            self._ints = coords
        return self._ints

    @internal_coordinates.setter
    def internal_coordinates(self, ics):
        if not isinstance(ics, CoordinateSet):
            raise ValueError("{} must be a {} to be valid internal coordinates".format(
                ics, CoordinateSet.__name__
            ))
        self._ints = ics

    def get_internals(self, strip_embedding=True):
        ics = self.internal_coordinates
        if ics is None:
            return None
        embedding_coords = self._get_embedding_coords()
        if embedding_coords is not None and strip_embedding:
            good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
            ics = ics.flatten()[good_coords]
        return ics

    @classmethod
    def _get_jacobian_storage(cls):
        return {
            'internals': [], #{"default":[], "fast":[], "generic":[]},
            'cartesian': []
        }
    def _get_int_jacobs(self,
                        jacs,
                        strip_dummies=False,
                        stencil=None, mesh_spacing=1.0e-3,
                        all_numerical=None, reembed=True,
                        planar_ref_tolerance=None,
                        parallelizer=None
                        ):
        """
        Gets the specified dX/dRs

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        intcds = self.internal_coordinates
        ccoords = self.coords
        carts = ccoords.system
        internals = intcds.system

        generics = 'GenericInternals' in internals.name
        zmatrix = 'ZMatrix' in internals.name

        if all_numerical is None:
            all_numerical = zmatrix

        converter_options = (
            dict(
                reembed=reembed,
                planar_ref_tolerance=planar_ref_tolerance,
                strip_dummies=strip_dummies
            )
                if zmatrix else
            {}
        )

        if isinstance(jacs, int):
            jacs = list(range(1, jacs + 1))

        exist_jacs = self._jacobians['internals']
        max_jac = max(jacs)
        need_jacs = [x + 1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
        if len(need_jacs) > 0:
            if generics:
                exist_jacs = [
                    x.squeeze() for x in
                    intcds.jacobian(carts, jacs,
                                    # odd behaves better
                                    mesh_spacing=mesh_spacing,
                                    stencil=stencil,
                                    all_numerical=all_numerical,
                                    converter_options=converter_options
                                    )
                ]
                self._jacobians['internals'] = exist_jacs
            else:
                stencil = (max(need_jacs) + 2 + (1 + max(need_jacs)) % 2) if stencil is None else stencil
                # odd behaves better
                with Parallelizer.lookup(parallelizer) as par:
                    new_jacs = [
                        x.squeeze() if isinstance(x, np.ndarray) else x
                        for x in intcds.jacobian(carts, need_jacs,
                                                 # odd behaves better
                                                 mesh_spacing=mesh_spacing,
                                                 stencil=stencil,
                                                 all_numerical=all_numerical,
                                                 converter_options=converter_options,
                                                 parallelizer=par
                                                 )
                    ]
                for j, v in zip(need_jacs, new_jacs):
                    for d in range(j - len(exist_jacs)):
                        exist_jacs.append(None)
                    exist_jacs[j - 1] = v

        return [exist_jacs[j - 1] for j in jacs]

    def _get_cart_jacobs(self, jacs,
                         strip_dummies=False,
                         stencil=None, mesh_spacing=1.0e-3,
                         all_numerical=None,
                         parallelizer=None
                         ):
        """
        Gets the specified dR/dXs

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        intcds = self.internal_coordinates
        ccoords = self.coords
        carts = ccoords.system
        internals = intcds.system
        generics = 'GenericInternals' in internals.name
        zmatrix = 'ZMatrix' in internals.name

        if all_numerical is None:
            all_numerical = zmatrix

        converter_options = (
            dict(strip_dummies=strip_dummies)
                if zmatrix else
            {}
        )

        if isinstance(jacs, int):
            jacs = list(range(1, jacs + 1))

        exist_jacs = self._jacobians['cartesian']
        max_jac = max(jacs)
        need_jacs = [x + 1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
        if len(need_jacs) > 0:
            if generics:
                exist_jacs = [
                    x.squeeze() for x in
                    ccoords.jacobian(internals, list(range(1, max_jac+1)),
                                     mesh_spacing=mesh_spacing,
                                     stencil=stencil,
                                     all_numerical=all_numerical,
                                     converter_options=converter_options
                                     )
                ]
            else:
                stencil = (max(need_jacs) + 2 + (1 + max(need_jacs)) % 2) if stencil is None else stencil
                # odd behaves better
                with Parallelizer.lookup(parallelizer) as par:
                    exist_jacs = [
                        x.squeeze() if isinstance(x, np.ndarray) else x
                        for x in ccoords.jacobian(internals,
                                                  order=list(range(1, max_jac+1)),
                                                  mesh_spacing=mesh_spacing,
                                                  stencil=stencil,
                                                  all_numerical=all_numerical,
                                                  converter_options=converter_options,
                                                  parallelizer=par
                                                  )
                    ]

        return [exist_jacs[j - 1] for j in jacs]

    @property
    def embedding_coords(self):
        return self._get_embedding_coords()
    def _get_embedding_coords(self):
        try:
            embedding = self.internal_coordinates.system.embedding_coords
        except AttributeError:
            try:
                embedding = self.internal_coordinates.system.converter_options['embedding_coords']
            except KeyError:
                embedding = None
        return embedding

    cartesian_by_internals_method = 'fast'
    def get_cartesians_by_internals(self, order=None, strip_embedding=False, reembed=True, method=None):
        if method is None:
            int_sys = self.internal_coordinates.system
            if "GenericInternals" in int_sys.name:
                method = 'classic'
            else:
                method = self.cartesian_by_internals_method

        if reembed and method == 'fast':
            L_base = self.get_translation_rotation_invariant_transformation(strip_embedding=strip_embedding, mass_weighted=False)
            jacs_1 = self.get_internals_by_cartesians(order, strip_embedding=strip_embedding)
            new_tf = nput.tensor_reexpand([L_base.T], jacs_1, order)
            inverse_tf = nput.inverse_transformation(new_tf, order, allow_pseudoinverse=True)
            return [
                np.tensordot(j, L_base, axes=[-1, -1])
                for j in inverse_tf
            ]
        elif not reembed and method == 'fast':
            wtf = self.get_internals_by_cartesians(order, strip_embedding=False) # faster to just do these derivs.
            base = nput.inverse_transformation(wtf, order-1, allow_pseudoinverse=True)
            # print(order, len(wtf), len(base))
        else:
            base = self._get_int_jacobs(order, reembed=reembed) if order is not None else self._jacobians['internals']
            if order is not None:
                if len(base) < order:
                    raise ValueError("insufficient {} (have {} but expected {})".format(
                        'CartesiansByInternals',
                        len(base),
                        order
                    ))
                base = base[:order]

            _ = []
            sh = self.coords.shape[:-2]
            nc = 3 * len(self.masses)
            nr = len(sh)
            n = -2 if base[0].shape[-2:] == (len(self.masses), 3) else -1
            for i, b in enumerate(base):
                rem = int(np.power(np.prod(b.shape[nr:n], dtype=int), 1 / (i+1)))
                b = b.reshape(sh + (rem,) * (i+1) + (nc,))
                _.append(b)
            base = _

        embedding_coords = self._get_embedding_coords() if strip_embedding else None
        if embedding_coords is not None and strip_embedding and base[0].shape[-2] == (3 * len(self.masses)):
            good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
            base = [t[np.ix_(*((good_coords,) * (t.ndim - 1)))] for t in base]
        return base

    def get_internals_by_cartesians(self, order=None, strip_embedding=False):
        base = self._get_cart_jacobs(order) if order is not None else self._jacobians['cartesian']
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'InternalsByCartesians',
                    len(base),
                    order
                ))
            base = base[:order]

        _ = []
        sh = self.coords.shape[:-2]
        nc = 3 * len(self.masses)
        for i, b in enumerate(base):
            b = b.reshape(sh + (nc,) * (i + 1) + (-1,))
            _.append(b)
        base = _

        if strip_embedding and base[0].shape[-1] == nc:
            embedding_coords = self._get_embedding_coords()
            if embedding_coords is not None:
                good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
                base = [t[..., good_coords] for t in base]
        return base

    def embed_coords(self, coords, sel=None, in_paf=False, planar_ref_tolerance=None):
        return StructuralProperties.get_eckart_embedded_coords(
            self.masses,
            self.coords,
            coords,
            sel=sel,
            in_paf=in_paf,
            planar_ref_tolerance=planar_ref_tolerance
        )

    @property
    def inertia_tensor(self):
        return StructuralProperties.get_prop_inertia_tensors(
            self.coords,
            self.masses
        )

    @property
    def inertial_frame(self):
        """
        Provides the inertial axis frame

        :return:
        :rtype:
        """

        if self._frame is None:
            # Need to put B in Hartree?
            #  I've got moments of inertia in amu * bohr^2 at the moment
            #  So we convert (amu * bohr^2) to (m_e * bohr^2) since hb^2/(m_e bohr^2) == E_h
            mom_i, eigs = StructuralProperties.get_prop_moments_of_inertia(
                self.coords,
                self.masses
            )
            B_e = 1 / (2 * mom_i)  # * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass"))
            self._inert_frame = B_e, eigs

        return self._inert_frame

    def inertial_frame_derivatives(self):

        return StructuralProperties.get_prop_inertial_frame_derivatives(
                self.coords,
                self.masses
            )

    @property
    def translation_rotation_modes(self):
        if self._tr_modes is None:
            self._tr_modes = StructuralProperties.get_prop_translation_rotation_eigenvectors(
                self.coords,
                self.masses
            )

        return self._tr_modes


    def get_translation_rotation_invariant_transformation(self,
                                                          order=0,
                                                          mass_weighted=True,
                                                          strip_embedding=True):
        L_tr = self.translation_rotation_modes[1]
        A = np.eye(L_tr.shape[0]) - (L_tr @ L_tr.T)
        evals, tf = np.linalg.eigh(A)
        zero_pos = np.abs(evals) < 1e-4 # the rest should be 1
        tf[:, zero_pos] = L_tr
        if strip_embedding:
            nzpos = np.abs(evals) > 1e-4
            tf = tf[:, nzpos]

        if not mass_weighted:
            tf = np.diag(np.repeat(1/np.sqrt(self.masses), 3)) @ tf

        return tf


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

class MolecularGenericInternalCoordinateSystem(GenericInternalCoordinateSystem):
    """
    Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding
    """
    name = "MolecularGenericInternals"
    def __init__(self, masses, coords, /, specs, converter_options=None,
                 redundant=False,
                 relocalize=False,
                 untransformed_coordinates=None,
                 angle_ordering='ijk',
                 **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        self.masses = masses
        self.coords = coords
        if converter_options is None:
            converter_options = dict(opts, reference_coordinates=coords, masses=masses, angle_ordering=angle_ordering, specs=specs)
            opts = {}
        if redundant:
            redundant_generator = RedundantCoordinateGenerator(
                specs,
                masses=self.masses,
                untransformed_coordinates=untransformed_coordinates,
                angle_ordering=angle_ordering,
                relocalize=relocalize
            )
            converter_options['redundant_generator'] = converter_options.get('redundant_generator', redundant_generator)
        nats = None#len(specs)
        super().__init__(converter_options=converter_options, dimension=(nats,), coordinate_shape=(nats,), opts=opts)

class MolecularZMatrixCoordinateSystem(ZMatrixCoordinateSystem):
    """
    Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding
    """
    name = "MolecularZMatrix"
    embedding_coords = [0, 1, 2, 4, 5, 8]
    def __init__(self, masses, coords, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        self.masses = masses
        self.coords = coords
        self.dummy_positions = [i for i, m in enumerate(masses) if m < 0]
        self.com = StructuralProperties.get_prop_center_of_mass(coords, masses)
        self.inertial_axes = StructuralProperties.get_prop_moments_of_inertia(coords, masses)[1]
        if converter_options is None:
            converter_options = opts
            opts = {}
        nats = len(masses)
        super().__init__(converter_options=converter_options, dimension=(nats, 3), coordinate_shape=(nats, 3), opts=opts)
        self.set_embedding()
    @property
    def origins(self):
        return self.converter_options['origins']
    @property
    def axes(self):
        return self.converter_options['axes']

    def pre_convert(self, system):
        # self.converter_options['molecule'] = self.molecule
        self.set_embedding()

    def set_embedding(self):
        com = self.com
        axes = self.inertial_axes
        converter_options = self.converter_options
        if 'ordering' in converter_options:
            ordering = np.array(converter_options['ordering'], dtype=int)
            ordering[0, 1] = -3; ordering[0, 2] = -1; ordering[0, 3] = -2
            ordering[1, 2] = -1; ordering[1, 3] = -2
            if len(ordering) > 2:
                ordering[2, 3] = -2
            converter_options['ordering'] = ordering
            first = ordering[0, 0]
        else:
            first = 0

        first_pos = self.coords[first]
        axes, ax_names, ax_choice = _get_best_axes(first_pos, axes)

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['axes_labels'] = ax_names
        converter_options['axes_choice'] = ax_choice
        converter_options['masses'] = self.masses
        converter_options['dummy_positions'] = self.dummy_positions
        converter_options['ref_coords'] = self.coords

    def jacobian(self,
                 coords,
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
            dummies = self.dummy_positions
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(self.masses),
                dummies
            )
        else:
            main_excludes = None

        if coords.ndim > 3:
            og_shape = coords.shape
            coords = coords.reshape((-1,) + coords.shape[-2:])
        else:
            og_shape = (coords.shape[0],) if coords.ndim == 3 else ()

        # print("...?", coords.shape)

        try:
            self.converter_options['reembed'] = True if remb is None else remb
            jacs = super().jacobian(coords, *args, converter_options=converter_options, **kwargs)
            if isinstance(jacs, np.ndarray):
                jacs = [jacs]
            raw_jacs = []
            for j in jacs:
                skip_dim = coords.ndim - 2
                if skip_dim > 0:
                    j = np.moveaxis(j, -3, 0) # skip_dim == 1 by construction so only need 1 move...
                ext_dim = j.ndim - 2 - skip_dim
                # print("????", j.shape, coords.shape, skip_dim)
                shp = og_shape + sum(
                    ((j.shape[i] // 3, 3) for i in range(skip_dim, ext_dim+skip_dim)),
                    ()
                ) + j.shape[-2:]
                j = j.reshape(shp)
                if main_excludes is not None:
                    for i in range(ext_dim):
                        j = np.take(j, main_excludes, axis=2*i)
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
    def __init__(self, masses, coords, dummy_positions=None, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        self.masses = masses
        self.coords = coords
        self.dummy_positions = [i for i,m in enumerate(masses) if m < 0] if dummy_positions is None else dummy_positions
        self.com = StructuralProperties.get_prop_center_of_mass(coords, masses)
        self.axes = StructuralProperties.get_prop_moments_of_inertia(coords, masses)[1]
        nats = len(masses)
        if converter_options is None:
            converter_options = opts
            opts = {}
        super().__init__(converter_options=converter_options, dimension=(nats, 3), coordinate_shape=(nats, 3), opts=opts)

    def pre_convert(self, system):
        self.converter_options['masses'] = self.masses
        if 'ZMatrix' in system.name:
            self.set_embedding()

    def set_embedding(self):
        """
        Sets up the embedding options...
        :return:
        :rtype:
        """
        com = self.com
        axes = self.axes
        converter_options = self.converter_options
        if 'ordering' in converter_options:
            ordering = np.array(converter_options['ordering'], dtype=int)
            ordering[0, 1] = -3; ordering[0, 2] = -2; ordering[0, 3] = -1
            ordering[1, 2] = -1; ordering[1, 3] = -2
            if len(ordering) > 2:
                ordering[2, 3] = -2
            converter_options['ordering'] = ordering
            first = ordering[0, 0]
        else:
            first = 0

        first_pos = self.coords[first]
        axes, ax_names, ax_choice = _get_best_axes(first_pos, axes)

        converter_options['origins'] = com
        converter_options['axes'] = axes
        converter_options['axes_labels'] = ax_names
        converter_options['axes_choice'] = ax_choice
        converter_options['masses'] = self.masses
        converter_options['ref_coords'] = self.coords
        converter_options['dummy_positions'] = self.dummy_positions

    def jacobian(self,
                 coords,
                 system,
                 order=None,
                 return_derivs=None,
                 strip_dummies=None,
                 converter_options=None,
                 analytic_deriv_order=None,
                 **kwargs
                 ):

        zmat_conv = 'ZMatrix' in system.name

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
            if zmat_conv:
                analytic_deriv_order = 0
            else:
                if return_derivs is None:
                    return_derivs = order
                analytic_deriv_order = return_derivs if nput.is_numeric(return_derivs) else max(return_derivs)

        if strip_dummies:
            dummies = self.dummy_positions
            if len(dummies) == 0:
                dummies = None
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(len(self.masses)),
                dummies
            )
        else:
            main_excludes = None

        jacs = super().jacobian(coords, system,
                                order=order,
                                analytic_deriv_order=analytic_deriv_order,
                                converter_options=converter_options,
                                **kwargs)
        if isinstance(jacs, np.ndarray):
            jacs = [jacs]

        if zmat_conv:
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
    def __init__(self, cart_system, zmat_system, **opts):
        self._types = (cart_system, zmat_system)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types
    def convert(self, coords, *,
                     masses,
                     dummy_positions, origins=None, axes=None, ordering=None, **kwargs):
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
                                       masses=masses,
                                       dummy_positions=dummy_positions,
                                       origins=origins, axes=axes, ordering=ordering, **kwargs)
        zmcs = zmcs[0]

        if 'derivs' in opts:
            derivs = opts['derivs']
            reshaped_derivs = [None] * len(derivs)
            for i, v in enumerate(derivs):
                reshaped_derivs[i] = v[0]
            opts['derivs'] = reshaped_derivs

        return zmcs, opts

    def convert_many(self,
                     coords,*,
                     masses,
                     dummy_positions,
                     origins=None, axes=None,
                     ordering=None,
                     strip_embedding=True,
                     strip_dummies=False,
                     return_derivs=None,
                     **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding

        :param coords: coordinates in Cartesians to convert
        :type coords: np.ndarray
        :param molecule:
        :type molecule: MolecularEmbedding
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
        return_derivs=False

        n_sys = coords.shape[0]
        n_coords = coords.shape[1]
        n_atoms = len(masses)

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
            if len(ordering) > 2:
                ordering[2, 3] = -1
            ordering = ordering + 3
            ordering = np.concatenate([ [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], ordering])
        res = CoordinateSet(coords, CartesianCoordinates3D).convert(ZMatrixCoordinates,
                                                                    ordering=ordering,
                                                                    origins=origins,
                                                                    axes=axes,
                                                                    return_derivs=return_derivs,
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
            dummies = [0, 1, 2] + [x+3 for x in dummy_positions] # add on axes
        elif strip_embedding:
            dummies = [0, 1, 2]
        else:
            dummies = None

        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(n_atoms + 3),
                dummies
            )
            sub_excludes = main_excludes - 1 # drop one fewer terms to drop I think...
            if 'derivs' in opts:
                derivs = opts['derivs']
                reshaped_derivs = [None] * len(derivs)
                deriv_excludes = np.arange(3, n_atoms + 3)
                for i, v in enumerate(derivs):
                    # drop all terms relating to the embedding of the embedding
                    start_dim = v.ndim - 2*(i+2)
                    for j in range(start_dim, v.ndim-2, 2):
                        v = np.take(v, deriv_excludes, axis=j)
                    v = np.take(v, sub_excludes, axis=-2)
                    reshaped_derivs[i] = v

                opts['derivs'] = reshaped_derivs

            zmcs = zmcs[..., sub_excludes, :]
        return zmcs, opts
# MolecularCartesianToZMatrixConverter = MolecularCartesianToZMatrixConverter()
# MolecularCartesianToZMatrixConverter.register(CoordinateSystemConverters)
class MolecularCartesianToRegularCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, cart_system, **opts):
        self._types = (CartesianCoordinates3D, cart_system)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types
    def convert(self, coords, **kw):
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        return coords, kwargs
# MolecularCartesianToRegularCartesianConverter = MolecularCartesianToRegularCartesianConverter()
# MolecularCartesianToRegularCartesianConverter.register()
class RegularCartesianToMolecularCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, cart_system, **opts):
        self._types = (cart_system, CartesianCoordinates3D)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types

    def convert(self, coords, **kw):
        return coords, kw

    def convert_many(self, coords, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        return coords, kwargs

class MolecularZMatrixToCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """
    def __init__(self, zmat_system, cart_system, **opts):
        self._types = (zmat_system, cart_system)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types

    def convert(self, coords, **kw):
        total_points, opts = self.convert_many(coords[np.newaxis], **kw)
        return total_points[0], opts

    def convert_many(self, coords, *,
                     masses, dummy_positions, ref_coords,
                     origins=None, axes=None, ordering=None,
                     reembed=False, axes_choice=None, return_derivs=None,
                     strip_dummies=False,
                     strip_embedding=True,
                     planar_ref_tolerance=None,
                     **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, attempting to preserve the embedding
        """
        from .Molecule import Molecule
        return_derivs = False

        n_sys = coords.shape[0]
        n_coords = coords.shape[1]
        n_atoms = len(masses)
        if n_coords != n_atoms + 2:
            # means we already added the embedding
            if n_coords != n_atoms:
                raise ValueError('Embedding unclear when num_coords ({}) < num_atoms ({})'.format(
                    n_coords,
                    n_atoms
                ))

            # TODO:
            # these need broadcasting when
            # doing jacobian calcs but this is currently pretty
            # hacky
            # print(axes.shape, n_sys)
            if axes.ndim == 2:
                axes = axes[np.newaxis]
            if axes.shape[0] > 1 and axes.shape[0] < n_sys:
                # print((n_sys // axes.shape[0], axes.shape[0]) + axes.shape)
                axes = np.reshape(
                    np.broadcast_to(
                        axes[np.newaxis],
                        (n_sys // axes.shape[0],) + axes.shape
                    ),
                    (n_sys, ) + axes.shape[1:]
                )

            if origins.ndim == 1:
                origins = origins[np.newaxis]
            # print(axes.shape, n_sys)
            if origins.shape[0] > 1 and origins.shape[0] < n_sys:
                # print((n_sys // axes.shape[0], axes.shape[0]) + axes.shape)
                origins = np.reshape(
                    np.broadcast_to(
                        origins[np.newaxis],
                        (n_sys // origins.shape[0],) + origins.shape
                    ),
                    (n_sys, ) + origins.shape[1:]
                )

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
            embed_carts = carts[..., 3:, :]
            reembed = not (
                    carts.squeeze().ndim == 2 and
                    np.allclose(ref_coords, embed_carts, atol=1.0e-5)
            ) # agree to like a ten thousandth of an angstrom
            if reembed:
                if not return_derivs:
                    embed_carts = StructuralProperties.get_eckart_embedded_coords(
                        masses,
                        ref_coords,
                        embed_carts,
                        planar_ref_tolerance=planar_ref_tolerance
                    )
                    carts = np.concatenate([
                        carts[..., :3, :],
                        embed_carts
                        ],
                    axis=-2
                    )
                else:
                    inert_coords, coord_coms, coord_axes = StructuralProperties.get_prop_principle_axis_data(
                        coords,
                        masses
                    ).principle_axis_data
                    if axes_choice is None:
                        axes_choice = (0, 1)
                    guh = self.convert_many(coords,
                                            origins=coord_coms,
                                            axes=coord_axes[:, axes_choice],
                                            masses=masses,
                                            ref_coords=ref_coords,
                                            dummy_positions=dummy_positions,
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
            dummies = [0, 1, 2] + [x + 3 for x in dummy_positions]  # add on axes
        elif strip_embedding:
            dummies = [0, 1, 2]
        else:
            dummies = None
        if dummies is not None:
            main_excludes = np.setdiff1d(
                np.arange(len(masses) + 3),
                dummies
            )
            sub_excludes = main_excludes - 1  # drop one fewer terms to drop I think...
            if 'derivs' in opts:
                derivs = opts['derivs']
                reshaped_derivs = [None] * len(derivs)
                deriv_excludes = np.arange(3, len(masses) + 3)
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
# MolecularZMatrixToCartesianConverter = MolecularZMatrixToCartesianConverter()
# MolecularZMatrixToCartesianConverter.register()

class MolecularZMatrixToRegularZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """
    def __init__(self, zmat_system, **opts):
        self._types = (zmat_system, ZMatrixCoordinateSystem)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types

    def convert(self, coords, **kw):
        return coords, kw

    def convert_many(self, coords, **kwargs):
        return coords, kwargs

class RegularZMatrixToMolecularZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, zmat_system, **opts):
        self._types = (ZMatrixCoordinateSystem, zmat_system)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types

    def convert(self, coords, **kw):
        return coords, kw

    def convert_many(self, coords, **kwargs):
        return coords, kwargs
# MolecularZMatrixToRegularZMatrixConverter = MolecularZMatrixToRegularZMatrixConverter()
# MolecularZMatrixToRegularZMatrixConverter.register()

class MolecularCartesianToGICConverter(CartesianToGICSystemConverter):
    """
    ...
    """
    def __init__(self, cart_system, zmat_system, **opts):
        self._types = (cart_system, zmat_system)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types

    def convert_many(self, coords, *,
                     order=0,
                     return_derivs=None,
                     redundant_generator:RedundantCoordinateGenerator=None,
                     **kw):
        """
        We'll implement this by having the ordering arg wrap around in coords?
        """
        internals, opts = super().convert_many(coords, order=order, return_derivs=return_derivs, **kw)
        if redundant_generator is not None:
            red_tf, red_exp = redundant_generator.get_redundant_transformation(
                opts['derivs'],
                untransformed_coordinates=redundant_generator.untransformed_coordinates,
                masses=redundant_generator.masses,
                relocalize=redundant_generator.relocalize
            )
            opts['reference_internals'] = internals
            internals = np.zeros(internals.shape[:-1] + red_tf.shape[-1:], dtype=internals.dtype)
            # internals = internals @ red_tf
            opts['derivs'] = red_exp
            opts['redundant_transformation'] = red_tf
            opts['redundant_inverse'] = np.moveaxis(red_tf, -1, -2)

        return internals, opts
# MolecularCartesianToZMatrixConverter = MolecularCartesianToZMatrixConverter()
# MolecularCartesianToZMatrixConverter.register(CoordinateSystemConverters)

class MolecularGICToCartesianConverter(GICSystemToCartesianConverter):
    """
    ...
    """
    def __init__(self, cart_system, zmat_system, **opts):
        self._types = (cart_system, zmat_system)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types
    def convert_many(self,
                     coords, *,
                     redundant_transformation=None,
                     redundant_inverse=None,
                     redundant_generator=None,
                     reference_internals=None,
                     **kw):

        if redundant_inverse is not None:
            coords = coords @ redundant_inverse + reference_internals
        carts, opts = super().convert_many(coords, **kw)
        if redundant_transformation is not None and 'derivs' in opts:
            opts['derivs'] = nput.tensor_reexpand([redundant_inverse], opts['derivs'])

        return carts, opts

class RegularGICToMolecularGICConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, **opts):
        self._types = (GenericInternalCoordinates, MolecularGenericInternalCoordinateSystem)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types

    def convert(self, coords, **kw):
        return coords, kw

    def convert_many(self, coords, **kwargs):
        return coords, kwargs
# MolecularZMatrixToRegularZMatrixConverter = MolecularZMatrixToRegularZMatrixConverter()
# MolecularZMatrixToRegularZMatrixConverter.register()

class MolecularGICConverterToRegularGIC(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, **opts):
        self._types = (MolecularGenericInternalCoordinateSystem, GenericInternalCoordinates)
        super().__init__(**opts)
    @property
    def types(self):
        return self._types

    def convert(self, coords, **kw):
        return coords, kw

    def convert_many(self, coords, **kwargs):
        return coords, kwargs
# MolecularZMatrixToRegularZMatrixConverter = MolecularZMatrixToRegularZMatrixConverter()
# MolecularZMatrixToRegularZMatrixConverter.register()