
import enum, numpy as np, functools as fp

from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer
from McUtils.Parallelizers import Parallelizer
from McUtils.Zachary import TensorDerivativeConverter, TensorExpansionTerms

class JacobianKeys(enum.Enum):
    CartesiansByInternals = "CartesiansByInternals"
    InternalsByCartesians = "InternalsByCartesians"
    InternalsByCartesianModes = "InternalsByModes"
    CartesianModesByInternals = "ModesByInternals"
    CartesiansByInternalModes = "CartesiansByModes"
    InternalModesByCartesians = "ModesByCartesians"
    CartesianModesByInternalModes = "CartesianModesByInternalModes"
    InternalModesByCartesianModes = "InternalModesByCartesianModes"
    InternalModesByInternals = "InternalModesByInternals"
    InternalsByInternalModes = "InternalsByInternalModes"
    CartesianModesByCartesians = "CartesianModesByCartesians"
    CartesiansByCartesianModes = "CartesiansByCartesianModes"


class DumbTensor:
    """
    A wrapper to make tensor algebra suck less
    """

    def __init__(self, tensor):
        self.t = tensor
    @property
    def shape(self):
        return self.t.shape
    @staticmethod
    def _dot(*t, axes=None):
        """
        Flexible tensordot
        """

        if len(t) == 1:
            return t[0]

        if any(isinstance(x, int) for x in t):
            return 0

        def tdot(a, b, **kw):
            if hasattr(a, "tensordot"):
                if 'axes' not in kw:
                    kw['axes'] = [-1, 0]
                td = a.tensordot(b, **kw)
            else:
                try:
                    td = np.tensordot(a, b, **kw)
                except ValueError:
                    if 'axes' not in kw:
                        axes = [-1, 0]
                    else:
                        axes = kw['axes']
                    raise ValueError("Shape-mismatch for sum: {} x {} along axes {}".format(a.shape, b.shape, axes))
            return td

        def td(a, b):
            if isinstance(a, int) or isinstance(b[0], int):
                res = 0
            else:
                res = tdot(a, b[0], axes=b[1])
            return res

        if axes is None:
            axes = [1] * (len(t) - 1)


        return fp.reduce(td, zip(t[1:], axes), t[0])

    def dot(self, b, *args, **kwargs):
        if isinstance(b, DumbTensor):
            b = b.t
        return type(self)(self._dot(self.t, b, *args, **kwargs))

    @staticmethod
    def _shift(a, *s):
        if isinstance(a, int):
            return a

        def shift_inds(n, i, j):
            if i < j:
                x = list(range(i)) + list(range(i + 1, j + 1)) + [i] + list(range(j + 1, n))
            else:
                x = list(range(j)) + [i] + list(range(j, i)) + list(range(i + 1, n))
            return x

        shiftIJ = lambda a, ij: np.transpose(a, shift_inds(a.ndim, *ij))
        return fp.reduce(shiftIJ, s, a)
    def shift(self, *args, **kwargs):
        return type(self)(self._shift(self.t, *args, **kwargs))
    def transpose(self, *perm):
        return type(self)(self.t.transpose(perm))

    @staticmethod
    def _contract_dim(R, targ_dim):
        # we figure out how much we're off by
        # and go from there, assuming that pairs of
        # dimensions to be contracted show up at the end
        for i in range(R.ndim - targ_dim):
            l_pos = R.ndim - (i + 2)
            gloobers = R.shape[:l_pos]
            if i > 0:
                r_pos = -i
                groobers = R.shape[r_pos:]
            else:
                groobers = ()
            R = R.reshape(gloobers + (-1,) + groobers)
        return R
    def contract_dim(self, targ_dim):
        return type(self)(self._contract_dim(self.t, targ_dim))

    def __add__(self, other):
        if isinstance(other, DumbTensor):
            other = other.t
        return type(self)(self.t+other)
    def __radd__(self, other):
        if isinstance(other, DumbTensor):
            other = other.t
        return type(self)(self.t+other)
    def __matmul__(self, other):
        return self.dot(other)
    def __getitem__(self, item):
        """
        :type item: slice
        """
        a = item.start
        b = item.stop
        return self.shift([a, b])

class MolecularFunctionExpander:
    """
    Handles all necessary details to be able to
    evaluate functions in molecular internal coordinates
    and then transform back to Cartesians, and to do the
    same with respect to normal dimensionless normal mode coordinates
    """

    def __init__(self,
                 cartesian_coordinates,
                 internal_coordinates,
                 transform=None,
                 masses=None, freqs=None,

                 parallelizer=None, logger=None,
                 checkpointer=None,

                 numerical_jacobians=True,
                 eckart_embed_derivatives=True,
                 eckart_embed_planar_ref_tolerance=None,
                 strip_dummies=False,
                 strip_embedding=True,
                 internal_fd_mesh_spacing=1.0e-3,
                 internal_fd_stencil=None,
                 cartesian_fd_mesh_spacing=1.0e-2,
                 cartesian_fd_stencil=None,
                 cartesian_analytic_deriv_order=0,
                 internal_by_cartesian_order=3,
                 cartesian_by_internal_order=4,
                 jacobian_warning_threshold=1e4
                 ):
        self._cached_transforms = {}

        self.cartesian_coordinates = cartesian_coordinates
        self.internal_coordinates = internal_coordinates
        self.modes = transform
        self.masses = np.ones(cartesian_coordinates.shape[-2]) if masses is None else masses
        self.freqs = np.ones(cartesian_coordinates.shape[-2]) if freqs is None else freqs
        if self.modes is not None:
            self.modes, self.inverse = self.undimensionalize(self.masses, self.freqs, transform)
        else:
            self.inverse = None

        self._internal_jacobians = []
        self._cartesian_jacobians = []
        self.parallelizer = parallelizer
        self.logger = logger

        self.internal_fd_mesh_spacing = internal_fd_mesh_spacing
        self.internal_fd_stencil = internal_fd_stencil
        self.cartesian_fd_mesh_spacing = cartesian_fd_mesh_spacing
        self.cartesian_fd_stencil = cartesian_fd_stencil
        self.cartesian_analytic_deriv_order = cartesian_analytic_deriv_order

        self.internal_by_cartesian_order = internal_by_cartesian_order
        self.cartesian_by_internal_order = cartesian_by_internal_order
        self.jacobian_warning_threshold = jacobian_warning_threshold

        self.strip_dummies = strip_dummies
        self.strip_embedding = strip_embedding

        self.reembed = eckart_embed_derivatives
        self.reembed_tol = eckart_embed_planar_ref_tolerance
        self.all_numerical = numerical_jacobians

        if logger is None:
            logger = NullLogger()
        self.logger = logger
        if parallelizer is None:
            parallelizer = Parallelizer.lookup(None)
        self.parallelizer = parallelizer
        if checkpointer is None:
            checkpointer = NullCheckpointer()
        self.checkpointer = checkpointer

    def _tripmass(self, masses):
        if self.strip_dummies:
            masses = masses[masses > 0]
        return np.broadcast_to(masses[np.newaxis, :], (3, len(masses))).T.flatten()
    def undimensionalize(self, masses, freqs, modes):
        """
        Removes units from normal modes

        :param masses:
        :type masses:
        :param modes:
        :type modes:
        :return:
        :rtype:
        """
        L = modes
        freq_conv = np.sqrt(np.broadcast_to(freqs[:, np.newaxis], L.shape))
        mass_conv = np.sqrt(np.broadcast_to(self._tripmass(masses)[np.newaxis, :], L.shape))
        conv = freq_conv * mass_conv
        L = L * conv
        Linv = (L / freq_conv**2)
        return L.T, Linv

    def get_int_jacobs(self, jacs):
        """
        Gets the specified Internal->Cartesian Jacobians

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        intcds = self.internal_coordinates
        ccoords = self.cartesian_coordinates
        carts = ccoords.system
        internals = intcds.system
        exist_jacs = self._internal_jacobians
        max_jac = max(jacs)
        need_jacs = [x+1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
        if len(need_jacs) > 0:
            stencil = (max(need_jacs) + 2 + (1+max(need_jacs))%2) if self.internal_fd_stencil is None else self.internal_fd_stencil
            # odd behaves better
            with Parallelizer.lookup(self.parallelizer) as par:
                new_jacs = [
                    x.squeeze() if isinstance(x, np.ndarray) else x
                    for x in intcds.jacobian(carts, need_jacs,
                                             # odd behaves better
                                             mesh_spacing=self.internal_fd_mesh_spacing,
                                             stencil=stencil,
                                             all_numerical=self.all_numerical,
                                             converter_options=dict(
                                                 reembed=self.reembed,
                                                 planar_ref_tolerance=self.reembed_tol,
                                                 strip_dummies=self.strip_dummies
                                             ),
                                             parallelizer=par
                                             )]
                # np.set_printoptions
                # with np.printoptions(linewidth=1e8, threshold=1e8, floatmode='fixed', precision=10):
                #     raise Exception(str(np.round(new_jacs[0].reshape(9, 9)[(3, 6, 7), :], 12)))
            for j,v in zip(need_jacs, new_jacs):
                for d in range(j-len(exist_jacs)):
                    exist_jacs.append(None)
                exist_jacs[j-1] = v

        return [exist_jacs[j-1] for j in jacs]
    def get_cart_jacobs(self, jacs):
        """
        Gets the specified Cartesian->Internal Jacobians

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        intcds = self.internal_coordinates
        ccoords = self.cartesian_coordinates
        carts = ccoords.system
        internals = intcds.system
        exist_jacs = self._cartesian_jacobians
        max_jac = max(jacs)
        need_jacs = [x+1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
        if len(need_jacs) > 0:
            stencil = (max(need_jacs) + 2 + (1+max(need_jacs))%2) if self.cartesian_fd_stencil is None else self.cartesian_fd_stencil
            # odd behaves better
            with Parallelizer.lookup(self.parallelizer) as par:
                new_jacs = [
                    x.squeeze() if isinstance(x, np.ndarray) else x
                    for x in ccoords.jacobian(internals, need_jacs,
                                                          mesh_spacing=self.cartesian_fd_mesh_spacing,
                                                          stencil=stencil,
                                                          # all_numerical=True,
                                                          analytic_deriv_order=self.cartesian_analytic_deriv_order,
                                                          converter_options=dict(strip_dummies=self.strip_dummies),
                                                          parallelizer=par
                                                          )
                ]
                if need_jacs[0] > self.cartesian_analytic_deriv_order:
                    new_jacs = new_jacs[self.cartesian_analytic_deriv_order:]

                for j, v in zip(need_jacs, new_jacs):
                    for d in range(j - len(exist_jacs)):
                        exist_jacs.append(None)
                    exist_jacs[j - 1] = v

        return [exist_jacs[j-1] for j in jacs]

    def _get_embedding_coords(self):
        try:
            embedding = self.internal_coordinates.system.embedding_coords
        except AttributeError:
            try:
                embedding = self.internal_coordinates.system.converter_options['embedding_coords']
            except KeyError:
                embedding = None
        return embedding
    def get_coordinate_transforms(self,
                                  internal_by_cartesian_order=None,
                                  cartesian_by_internal_order=None
                                  ):

        if internal_by_cartesian_order is None:
            internal_by_cartesian_order = self.internal_by_cartesian_order
        if cartesian_by_internal_order is None:
            cartesian_by_internal_order = self.cartesian_by_internal_order

        current_cache = self._cached_transforms

        if (
                current_cache is None
                or len(current_cache[JacobianKeys.CartesiansByInternals]) < cartesian_by_internal_order
                or len(current_cache[JacobianKeys.InternalsByCartesians]) < internal_by_cartesian_order
        ):

            if current_cache is None:
                current_cache = {}

            embedding_coords = self._get_embedding_coords() if self.strip_embedding else None
            # fill out
            if (
                    JacobianKeys.CartesiansByInternals not in current_cache
                    or len(current_cache[JacobianKeys.CartesiansByInternals]) < cartesian_by_internal_order
            ):
                cart_by_internal_jacobs = self.get_int_jacobs(list(range(1, cartesian_by_internal_order+1)))

                # The finite difference preserves too much shape by default
                _contract_dim = DumbTensor._contract_dim
                _ = []
                for i,x in enumerate(cart_by_internal_jacobs):
                    if isinstance(x, int) or x.ndim == 2+i:
                        _.append(x)
                    elif x.ndim > 2+i:
                        _.append(_contract_dim(x, 2+i))
                    else:
                        raise ValueError("bad shape for Cartesian by internal jacobian {} ({})".format(
                            i, x.shape
                        ))
                cart_by_internal_jacobs = _

                # we'll strip off the embedding coords just in case
                if embedding_coords is not None:
                    good_coords = np.setdiff1d(np.arange(3*len(self.masses)), embedding_coords)

                for i, x in enumerate(cart_by_internal_jacobs):
                    bad_spots = np.where(np.abs(x) > self.jacobian_warning_threshold)
                    bad_bad_spots = bad_spots  # so we don't lose it
                    if len(bad_spots) > 0:  # numpy fuckery
                        bad_spots = bad_spots[0]
                    if len(bad_spots) > 0:
                        m = np.max(np.abs(x[bad_bad_spots]))
                        self.logger.log_print('WARNING: maximum d^{i}X/dR^{i} term is {m}. '
                                              'This will likely mess up G-matrix terms and is probably coming from a planar structure. '
                                              'Setting to zero, but `jacobian_warning_threshold` can be increased if this is expected '
                                              'All terms >{t} (base shape:{s}): {b}',
                                              i=i+1,
                                              m=m,
                                              s=x.shape,
                                              b=np.array(bad_bad_spots).T.tolist(),
                                              t=self.jacobian_warning_threshold
                                              )
                        x[bad_bad_spots] = 0.
                        # raise Exception(";_;")


                # Need to then mass weight
                masses = self.masses
                mass_conv = np.sqrt(self._tripmass(masses))
                # mass weight the derivs w.r.t internals
                internal_weighting = mass_conv
                _ = []
                for i, x in enumerate(cart_by_internal_jacobs):
                    internal_weighting = np.expand_dims(internal_weighting, 0)
                    if isinstance(x, int):
                        _.append(x)
                    else:
                        x = x * internal_weighting
                        if embedding_coords is not None:
                            for j in range(i+1):
                                x = np.take(x, good_coords, axis=j)
                        _.append(x)
                cart_by_internal_jacobs = _

                current_cache[JacobianKeys.CartesiansByInternals] = cart_by_internal_jacobs

            else:
                cart_by_internal_jacobs = current_cache[JacobianKeys.CartesiansByInternals]

            if (
                    JacobianKeys.InternalsByCartesians not in current_cache
                    or len(current_cache[JacobianKeys.InternalsByCartesians]) < internal_by_cartesian_order
            ):
                int_by_cartesian_jacobs = self.get_cart_jacobs(list(range(1, internal_by_cartesian_order + 1)))

                _contract_dim = DumbTensor._contract_dim
                _ = []
                for i,x in enumerate(int_by_cartesian_jacobs):
                    if isinstance(x, int) or x.ndim == 2+i:
                        _.append(x)
                    elif x.ndim > 2+i:
                        _.append(_contract_dim(x, 2+i))
                    else:
                        raise ValueError("bad shape for internal by Cartesian jacobian {} ({})".format(
                            i, x.shape
                        ))
                int_by_cartesian_jacobs = _

                # we'll strip off the embedding coords just in case
                if embedding_coords is not None:
                    good_coords = np.setdiff1d(np.arange(3*len(self.masses)), embedding_coords)

                for i,x in enumerate(int_by_cartesian_jacobs):
                    bad_spots = np.where(np.abs(x) > self.jacobian_warning_threshold)
                    bad_bad_spots = bad_spots # so we don't lose it
                    if len(bad_spots) > 0: # numpy fuckery
                        bad_spots = bad_spots[0]
                    if len(bad_spots) > 0:
                        m = np.max(np.abs(x[bad_bad_spots]))
                        self.logger.log_print('WARNING: maximum d^{i}R/dX^{i} term is {m}. '
                                              'This will likely mess up G-matrix terms and is probably coming from a planar structure. '
                                              'Setting to zero, but `jacobian_warning_threshold` can be increased if this is expected. '
                                              'All terms >{t} (base shape:{s}): {b}',
                                              i=i+1,
                                              m=m,
                                              s=x.shape,
                                              b=np.array(bad_bad_spots).T,
                                              t=self.jacobian_warning_threshold
                                              )
                        x[bad_bad_spots] = 0.

                # Need to then mass weight
                masses = self.masses
                mass_conv = np.sqrt(self._tripmass(masses))
                # mass weight the derivs w.r.t cartesians
                cartesian_weighting = mass_conv
                mc = mass_conv
                _ = []
                for i, x in enumerate(int_by_cartesian_jacobs):
                    cartesian_weighting = np.expand_dims(cartesian_weighting, -1)#[..., np.newaxis]
                    if isinstance(x, int):
                        _.append(x)
                    else:
                        x = x / cartesian_weighting
                        if embedding_coords is not None:
                            x = np.take(x, good_coords, axis=-1)
                        _.append(x)
                    mc = np.expand_dims(mc, 0)
                    cartesian_weighting = cartesian_weighting * mc
                int_by_cartesian_jacobs = _

                current_cache[JacobianKeys.InternalsByCartesians] = int_by_cartesian_jacobs
            else:
                int_by_cartesian_jacobs = current_cache[JacobianKeys.InternalsByCartesians]

            QY = self.modes  # derivatives of Q with respect to the Cartesians
            YQ = self.inverse # derivatives of Cartesians with respect to Q

            if QY is not None:

                if (
                        JacobianKeys.InternalsByCartesianModes not in current_cache
                        or len(current_cache[JacobianKeys.InternalsByCartesianModes]) < internal_by_cartesian_order
                ):
                    RQ_derivs = TensorDerivativeConverter(
                        [YQ] + [0]*(len(int_by_cartesian_jacobs) - 1),
                        int_by_cartesian_jacobs
                    ).convert(order=len(int_by_cartesian_jacobs))#, check_arrays=True)
                    current_cache[JacobianKeys.InternalsByCartesianModes] = RQ_derivs
                else:
                    RQ_derivs = current_cache[JacobianKeys.InternalsByCartesianModes]

                if (
                        JacobianKeys.CartesianModesByInternals not in current_cache
                        or len(current_cache[JacobianKeys.CartesianModesByInternals]) < cartesian_by_internal_order
                ):
                    QR_derivs = TensorDerivativeConverter(
                        cart_by_internal_jacobs,
                        [QY] + [0]*(len(cart_by_internal_jacobs) - 1)
                    ).convert(order=len(cart_by_internal_jacobs))
                    current_cache[JacobianKeys.CartesianModesByInternals] = QR_derivs
                else:
                    QR_derivs = current_cache[JacobianKeys.CartesianModesByInternals]

                if (
                        JacobianKeys.CartesiansByInternalModes not in current_cache
                        or len(current_cache[JacobianKeys.CartesiansByInternalModes]) < len(cart_by_internal_jacobs)
                ):
                    x_derivs = cart_by_internal_jacobs#(YR, YRR, YRRR, YRRRR)
                    Q_derivs = RQ_derivs[:1] + [0]*(len(cart_by_internal_jacobs) - 1)
                    YQ_derivs = TensorDerivativeConverter(Q_derivs, x_derivs,
                                                          jacobians_name='Q',
                                                          values_name='X'
                                                          ).convert(order=len(cart_by_internal_jacobs))#, check_arrays=True)
                    # self._get_tensor_derivs(
                    #     YQ_derivs, (QY, 0, 0, 0),
                    #     mixed_XQ=False
                    # )

                    current_cache[JacobianKeys.CartesiansByInternalModes] = YQ_derivs

                if (
                        JacobianKeys.CartesianModesByInternalModes not in current_cache
                        or len(current_cache[JacobianKeys.CartesianModesByInternalModes]) < len(cart_by_internal_jacobs)
                ):
                    # modes = self.molecule.translation_rotation_modes[1]
                    # # raise Exception(
                    # #     modes @ modes.T,
                    # #     modes.T @ modes,
                    # #     self.molecule.translation_rotation_modes[1].shape,
                    # #     self.modes.inverse.shape
                    # # )
                    # YQ2 = np.concatenate([
                    #     modes,
                    #     self.modes.inverse  # derivatives of Cartesians with respect to Q
                    # ], axis=0)

                    YQ_derivs = current_cache[JacobianKeys.CartesiansByInternalModes]
                    qQ_derivs = TensorDerivativeConverter(YQ_derivs,
                                                          [QY] + [0] * (len(cart_by_internal_jacobs) - 1),
                                                          jacobians_name='YQ',
                                                          values_name='qY'
                                                          ).convert(order=len(cart_by_internal_jacobs))#, check_arrays=True)
                    # raise Exception(qQ_derivs[0][6:, 6:])
                    current_cache[JacobianKeys.CartesianModesByInternalModes] = qQ_derivs

                if (
                        JacobianKeys.InternalModesByCartesians not in current_cache
                        or len(current_cache[JacobianKeys.InternalModesByCartesians]) < len(int_by_cartesian_jacobs)
                ):
                    QR = QR_derivs[0]
                    QY_derivs = TensorDerivativeConverter(int_by_cartesian_jacobs,
                                                          [QR] + [0]*(len(int_by_cartesian_jacobs) - 1)
                                                          ).convert(order=len(int_by_cartesian_jacobs))#, check_arrays=True)
                    current_cache[JacobianKeys.InternalModesByCartesians] = QY_derivs

                if (
                        JacobianKeys.InternalModesByCartesianModes not in current_cache
                        or len(current_cache[JacobianKeys.InternalModesByCartesianModes]) < len(int_by_cartesian_jacobs)
                ):
                    RQ_derivs = current_cache[JacobianKeys.InternalsByCartesianModes]
                    QR = QR_derivs[0]
                    Qq_derivs = TensorDerivativeConverter(RQ_derivs,
                                                          [QR] + [0] * (len(RQ_derivs) - 1),
                                                          jacobians_name='Rq',
                                                          values_name='qR'
                                                          ).convert(order=len(RQ_derivs))
                    current_cache[JacobianKeys.InternalModesByCartesianModes] = Qq_derivs

                if (
                        JacobianKeys.CartesianModesByCartesians not in current_cache
                        or len(current_cache[JacobianKeys.CartesianModesByCartesians]) < len(cart_by_internal_jacobs)
                ):
                    current_cache[JacobianKeys.CartesianModesByCartesians] = [self.modes.matrix] + [0]*(len(cart_by_internal_jacobs)-1)
                if (
                        JacobianKeys.CartesiansByCartesianModes not in current_cache
                        or len(current_cache[JacobianKeys.CartesiansByCartesianModes]) < len(cart_by_internal_jacobs)
                ):
                    current_cache[JacobianKeys.CartesiansByCartesianModes] = [self.modes.inverse] + [0] * (len(cart_by_internal_jacobs) - 1)

                if (
                        JacobianKeys.InternalModesByInternals not in current_cache
                        or len(current_cache[JacobianKeys.InternalModesByInternals]) < len(int_by_cartesian_jacobs)
                ):
                    YR = current_cache[JacobianKeys.CartesiansByInternals][0]
                    QY = current_cache[JacobianKeys.InternalModesByCartesians][0]
                    current_cache[JacobianKeys.InternalModesByInternals] = [YR@QY] + [0]*(len(int_by_cartesian_jacobs)-1)

                if (
                        JacobianKeys.InternalsByInternalModes not in current_cache
                        or len(current_cache[JacobianKeys.InternalsByInternalModes]) < len(int_by_cartesian_jacobs)
                ):
                    RY = current_cache[JacobianKeys.InternalsByCartesians][0]
                    YQ = current_cache[JacobianKeys.CartesiansByInternalModes][0]
                    current_cache[JacobianKeys.InternalsByInternalModes] = [YQ@RY] + [0]*(len(int_by_cartesian_jacobs)-1)

        return current_cache#self._cached_transforms[self.molecule]


    @property
    def cartesian_L_matrix(self):
        return self.get_cartesians_by_cartesian_modes(1)[0]
    def get_cartesians_by_cartesian_modes(self, order=None):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=order,
            internal_by_cartesian_order=None if order is None else min(order, self.internal_by_cartesian_order)
        )[JacobianKeys.CartesiansByCartesianModes]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'CartesiansByInternalModes',
                    len(base),
                    order
                ))
            base = base[:order]
        return base
    @property
    def cartesian_L_inverse(self):
        return self.get_cartesian_modes_by_cartesians(1)[0]
    def get_cartesian_modes_by_cartesians(self, order=None):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=order,
            internal_by_cartesian_order=None if order is None else min(order, self.internal_by_cartesian_order)
        )[JacobianKeys.CartesianModesByCartesians]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'CartesiansByInternalModes',
                    len(base),
                    order
                ))
            base = base[:order]
        return base

    @property
    def internal_L_matrix(self):
        return self.get_internal_modes_by_internals(1)[0]
    def get_internal_modes_by_internals(self, order=None, strip_embedding=True):
        # print(dict(
        #     cartesian_by_internal_order=order,
        #     internal_by_cartesian_order=min(order, self.internal_by_cartesian_order)
        # ))
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=order,
            internal_by_cartesian_order=None if order is None else min(order, self.internal_by_cartesian_order)
        )[JacobianKeys.InternalModesByInternals]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'InternalModesByInternal',
                    len(base),
                    order
                ))
            base = base[:order]

        embedding_coords = self._get_embedding_coords() if (self.strip_embedding or strip_embedding) else None
        if embedding_coords is not None and (strip_embedding and not self.strip_embedding):
            good_coords = np.setdiff1d(np.arange(3 * self.num_atoms), embedding_coords)
            base = [t[good_coords,] if not isinstance(t, int) else t for t in base]
        return base
    @property
    def internal_L_inverse(self):
        return self.get_internals_by_internal_modes(1)[0]
    def get_internals_by_internal_modes(self, order=None, strip_embedding=True):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=order,
            internal_by_cartesian_order=None if order is None else min(order, self.internal_by_cartesian_order)
        )[JacobianKeys.InternalsByInternalModes]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'CartesiansByInternalModes',
                    len(base),
                    order
                ))
            base = base[:order]
        embedding_coords = self._get_embedding_coords() if (self.strip_embedding or strip_embedding) else None
        if embedding_coords is not None and (strip_embedding and not self.strip_embedding):
            good_coords = np.setdiff1d(np.arange(3 * self.num_atoms), embedding_coords)
            base = [t[..., good_coords] if not isinstance(t, int) else t for t in base]
        return base
    @property
    def cartesians_by_modes(self):
        return self.get_cartesians_by_modes()
    def get_cartesians_by_modes(self, order=None):
        # print(dict(
        #     cartesian_by_internal_order=order,
        #     internal_by_cartesian_order=min(order, self.internal_by_cartesian_order)
        # ))
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=order,
            internal_by_cartesian_order=None if order is None else min(order, self.internal_by_cartesian_order)
        )[JacobianKeys.CartesiansByInternalModes]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'CartesiansByInternalModes',
                    len(base),
                    order
                ))
            base = base[:order]
        return base
    @property
    def modes_by_cartesians(self):
        return self.get_coordinate_transforms()[JacobianKeys.InternalModesByCartesians]
    def get_modes_by_cartesians(self, order=None, strip_embedding=True):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=None if order is None else min(order, self.cartesian_by_internal_order),
            internal_by_cartesian_order=order
        )[JacobianKeys.InternalModesByCartesians]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'InternalModesByCartesians',
                    len(base),
                    order
                ))
            base = base[:order]
        return base
    @property
    def cartesians_by_internals(self):
        return self.get_coordinate_transforms()[JacobianKeys.CartesiansByInternals]
    def get_cartesians_by_internals(self, order=None, strip_embedding=False):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=order,
            internal_by_cartesian_order=None if order is None else min(order, self.internal_by_cartesian_order)
        )[JacobianKeys.CartesiansByInternals]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'CartesiansByInternals',
                    len(base),
                    order
                ))
            base = base[:order]

        embedding_coords = self._get_embedding_coords() if (self.strip_embedding or strip_embedding) else None
        if embedding_coords is not None and (strip_embedding and not self.strip_embedding):
            good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
            base = [t[np.ix_(*((good_coords,)*(t.ndim-1)))] for t in base]
        return base
    @property
    def internals_by_cartesians(self):
        return self.get_coordinate_transforms()[JacobianKeys.InternalsByCartesians]
    def get_internals_by_cartesians(self, order=None, strip_embedding=False):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=None if order is None else min(order, self.cartesian_by_internal_order),
            internal_by_cartesian_order=order
        )[JacobianKeys.InternalsByCartesians]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'InternalsByCartesians',
                    len(base),
                    order
                ))
            base = base[:order]
        embedding_coords = self._get_embedding_coords() if (self.strip_embedding or strip_embedding) else None
        if embedding_coords is not None and (strip_embedding and not self.strip_embedding):
            good_coords = np.setdiff1d(np.arange(3 * len(self.masses)), embedding_coords)
            base = [t[..., good_coords] for t in base]
        return base

    @property
    def cartesian_modes_by_internal_modes(self):
        return self.get_coordinate_transforms()[JacobianKeys.CartesianModesByInternalModes]
    def get_cartesian_modes_by_internal_modes(self, order=None):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=order,
            internal_by_cartesian_order=None if order is None else min(order, self.internal_by_cartesian_order)
        )[JacobianKeys.CartesianModesByInternalModes]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'CartesianModesByInternalModes',
                    len(base),
                    order
                ))
            base = base[:order]
        return base

    @property
    def internal_modes_by_cartesian_modes(self):
        return self.get_coordinate_transforms()[JacobianKeys.InternalModesByCartesianModes]
    def get_internal_modes_by_cartesian_modes(self, order=None):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=None if order is None else min(order, self.cartesian_by_internal_order),
            internal_by_cartesian_order=order
        )[JacobianKeys.InternalModesByCartesianModes]
        if order is not None:
            if len(base) < order:
                raise ValueError("insufficient {} (have {} but expected {})".format(
                    'InternalModesByCartesianModes',
                    len(base),
                    order
                ))
            base = base[:order]
        return base

    def evaluate(self, function, deriv_order=2):
        ...
