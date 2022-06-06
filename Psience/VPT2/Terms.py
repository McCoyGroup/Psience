"""
Stores all of the terms used inside the VPT2 representations
"""
import itertools

import numpy as np, functools as fp, itertools as ip, time, enum

from McUtils.Numputils import SparseArray, levi_cevita3, vec_tensordot, vec_outer
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer
from McUtils.Parallelizers import Parallelizer
from McUtils.Zachary import TensorDerivativeConverter, TensorExpansionTerms

from ..Molecools import Molecule, MolecularVibrations, MolecularNormalModes

from .Common import PerturbationTheoryException

__all__ = [
    "ExpansionTerms",
    "KineticTerms",
    "PotentialTerms",
    "DipoleTerms",
    "CoriolisTerm",
    "PotentialLikeTerm"
]

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

class MixedDerivativeHandlingModes(enum.Enum):
    Unhandled = "unhandled"
    Numerical = "numerical"
    Analytical = "analytical"
    Averaged = "averaged"

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

class ExpansionTerms:
    """
    Base class for kinetic, potential, and dipole derivative terms
    """

    # backpropagate_internals = False # just a flag that can be set to use Cartesian results _but_ do it with
    #                                # terms backpropagated from the internals
    # mixed_derivative_handling_mode = "unhandled"
    # undimensionalize_normal_modes = True
    # numerical_jacobians = True
    # eckart_embed_derivatives = True
    # strip_dummy_atoms = False
    # strip_embedding_coordinates = False

    # so they can be tracked/propagated up more easily
    __props__ = (
        "logger",
        "parallelizer",
        "checkpointer",
        "undimensionalize",
        "numerical_jacobians",
        "eckart_embed_derivatives",
        "eckart_embed_planar_ref_tolerance",
        "strip_dummies",
        "strip_embedding",
        "mixed_derivative_handling_mode",
        "backpropagate_internals",
        "direct_propagate_cartesians",
        "zero_mass_term",
        "internal_fd_mesh_spacing",
        "internal_fd_stencil",
        "cartesian_fd_mesh_spacing",
        "cartesian_fd_stencil",
        "cartesian_analytic_deriv_order",
        "internal_by_cartesian_order",
        "cartesian_by_internal_order",
        "jacobian_warning_threshold",
        "coordinate_transformations",
        "coordinate_derivatives",
    )
    _cached_jacobians = {}
    def __init__(self,
                 molecule,
                 modes=None,
                 mode_selection=None,
                 use_internal_modes=None,
                 logger=None,
                 parallelizer=None,
                 checkpointer=None,
                 undimensionalize=None,
                 numerical_jacobians=True,
                 eckart_embed_derivatives=True,
                 eckart_embed_planar_ref_tolerance=None,
                 strip_dummies=False,
                 strip_embedding=True,
                 mixed_derivative_handling_mode="unhandled",
                 backpropagate_internals=False,
                 direct_propagate_cartesians=False,
                 zero_mass_term=1e7,
                 internal_fd_mesh_spacing=1.0e-3,
                 internal_fd_stencil=9,
                 cartesian_fd_mesh_spacing=1.0e-3,
                 cartesian_fd_stencil=9,
                 cartesian_analytic_deriv_order=0,
                 internal_by_cartesian_order=3,
                 cartesian_by_internal_order=4,
                 jacobian_warning_threshold=1e4,
                 coordinate_transformations=None,
                 coordinate_derivatives=None
                 ):
        """
        :param molecule: the molecule we're doing the expansion for
        :type molecule: Molecule
        :param modes: normal modes in Cartesian coordinates
        :type modes: MolecularVibrations
        :param mode_selection: the selection of modes to use
        :type mode_selection: None | Iterable[int]
        :param undimensionalize: whether or not we need to do some units fuckery on the modes
        :type undimensionalize: bool
        """
        self._terms = None
        self.molecule = molecule

        self.strip_dummies = strip_dummies
        self.strip_embedding = strip_embedding
        self.backpropagate_internals = backpropagate_internals
        self.direct_propagate_cartesians = direct_propagate_cartesians

        self.zero_mass_term = zero_mass_term

        self.internal_fd_mesh_spacing = internal_fd_mesh_spacing
        self.internal_fd_stencil = internal_fd_stencil
        self.cartesian_fd_mesh_spacing = cartesian_fd_mesh_spacing
        self.cartesian_fd_stencil = cartesian_fd_stencil
        self.cartesian_analytic_deriv_order = cartesian_analytic_deriv_order

        self.internal_by_cartesian_order = internal_by_cartesian_order
        self.cartesian_by_internal_order = cartesian_by_internal_order
        self.jacobian_warning_threshold = jacobian_warning_threshold

        self.internal_coordinates = molecule.internal_coordinates
        self.coords = molecule.coords
        self.masses = molecule.masses * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        self.use_internal_modes = use_internal_modes
        if modes is None:
            modes = molecule.normal_modes.modes
        self._modes = modes.basis
        if undimensionalize is None:
            undimensionalize = not self._check_internal_modes(clean=False)
        if undimensionalize:
            self.raw_modes = modes
            modes = self.undimensionalize(self.masses, self._modes)
        else:
            self.raw_modes = None
            modes = self._modes
        if mode_selection is not None:
            modes = modes[mode_selection]
        self._modes = modes
        self.mode_sel = mode_selection
        self.freqs = self.modes.freqs
        self._inert_frame = None

        self.reembed=eckart_embed_derivatives
        self.reembed_tol=eckart_embed_planar_ref_tolerance
        self.all_numerical=numerical_jacobians

        if logger is None:
            logger = NullLogger()
        self.logger = logger
        if parallelizer is None:
            parallelizer = Parallelizer.lookup(None)
        self.parallelizer = parallelizer
        if checkpointer is None:
            checkpointer = NullCheckpointer()
        self.checkpointer = checkpointer

        if coordinate_derivatives is not None:
            self._cached_jacobians[molecule] = coordinate_derivatives

        if coordinate_transformations is not None:
            self._cached_transforms[molecule] = coordinate_transformations

        if not isinstance(mixed_derivative_handling_mode, MixedDerivativeHandlingModes):
            self.mixed_derivative_handling_mode = MixedDerivativeHandlingModes(mixed_derivative_handling_mode)

    @property
    def num_atoms(self):
        """
        Gets the number of atoms (excluding dummies if `strip_dummies` is `True`)

        :return:
        :rtype:
        """
        if self.strip_dummies:
            n = np.sum(self.masses > 0, dtype=int)
        else:
            n = len(self.masses)
        return n

    def _check_internal_modes(self, modes=None, clean=True):
        if self.use_internal_modes is not None:
            if clean and self.use_internal_modes:
                self._reshape_internal_modes()
            return self.use_internal_modes
        if modes is None:
            modes = self._modes
        mat = modes.matrix
        is_internal = mat.shape[0] == self.coords.shape[0] * self.coords.shape[1] - 6
        self.use_internal_modes = is_internal
        if clean and is_internal:
            self._reshape_internal_modes()
        return is_internal

    def _reshape_internal_modes(self):
        QR = self._modes.matrix  # derivatives of Q with respect to the internals
        # we need to add zeros for the orientation coordinates
        if not self.strip_embedding and QR.shape[0] != 3 * self.num_atoms:
            _QR = QR
            QR = np.zeros((3 * self.num_atoms, _QR.shape[1]))
            embedding_coords = [0, 1, 2, 4, 5, 8]
            good_coords = np.setdiff1d(np.arange(3 * self.num_atoms), embedding_coords)
            QR[good_coords, :] = _QR
            self._modes.matrix = QR

        RQ = self._modes.inverse  # derivatives of internals with respect to Q
        if not self.strip_embedding and RQ.shape[1] != 3 * self.num_atoms:
            _RQ = RQ
            # we need to add zeros for the orientation coordinates
            RQ = np.zeros((_RQ.shape[0], 3 * self.num_atoms))
            embedding_coords = [0, 1, 2, 4, 5, 8]
            good_coords = np.setdiff1d(np.arange(3 * self.num_atoms), embedding_coords)
            RQ[:, good_coords] = _RQ
            self._modes.inverse = RQ

    @property
    def modes(self):
        # if self._check_internal_modes():
        #     J, = self.get_cart_jacobs([1])
        #     return np.dot(J, self._modes)
        # else:
        #     # cartesian modes
        return self._modes

    def undimensionalize(self, masses, modes):
        """
        Removes units from normal modes

        :param masses:
        :type masses:
        :param modes:
        :type modes:
        :return:
        :rtype:
        """
        L = modes.matrix.T
        freqs = modes.freqs
        freq_conv = np.sqrt(np.broadcast_to(freqs[:, np.newaxis], L.shape))
        if self._check_internal_modes(clean=False):
            conv = freq_conv
            L = L * conv
            Linv = modes.inverse / conv
        else:
            mass_conv = np.sqrt(np.broadcast_to(self._tripmass(masses)[np.newaxis, :], L.shape))
            conv = freq_conv * mass_conv
            L = L * conv
            Linv = (L / freq_conv**2)
        modes = type(modes)(self.molecule, L.T, inverse=Linv, freqs=freqs)
        return modes

    def _tripmass(self, masses):
        if self.strip_dummies:
            masses = masses[masses > 0]
        else:
            masses = masses.copy()
            masses[masses < 0] = self.zero_mass_term
        return np.broadcast_to(masses[np.newaxis, :], (3, len(masses))).T.flatten()

    def get_terms(self, order=None):
        """
        Gets the terms up to the given order

        :param order:
        :type order:
        :return:
        :rtype:
        """
        raise NotImplementedError("base class")

    def get_term(self, t):
        """
        Provides the term at order `t`

        :param t:
        :type t:
        :return:
        :rtype:
        """
        if self._terms is None or len(self._terms) < t+1:
            self._terms = self.get_terms(order=t)
        return self._terms[t]

    @property
    def terms(self):
        if self._terms is None:
            self._terms = self.get_terms()
        return self._terms

    def __getitem__(self, item):
        return self.get_term(item)

    @staticmethod
    def _weight_derivatives(t, order = None):
        if isinstance(t, int):
            return t
        weighted = t
        if order is None:
            order = len(t.shape)
        if order > 1:
            s = t.shape
            weights = np.ones(s)
            all_inds = list(range(len(s)))
            for i in range(2, order + 1):
                for inds in ip.combinations(all_inds, i):
                    # define a diagonal slice through
                    sel = tuple(slice(None, None, None) if a not in inds else np.arange(s[a]) for a in all_inds)
                    weights[sel] = 1 / np.math.factorial(i)
            weighted = weighted * weights
            # print(weights, weighted.array)
        return weighted

    def get_int_jacobs(self, jacs):
        """
        Gets the specified Internal->Cartesian Jacobians

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        intcds = self.internal_coordinates
        ccoords = self.coords
        carts = ccoords.system
        internals = intcds.system
        if self.molecule not in self._cached_jacobians:
            self._cached_jacobians[self.molecule] = {
                'int': [],
                'cart': []
            }
        exist_jacs = self._cached_jacobians[self.molecule]['int']
        max_jac = max(jacs)
        need_jacs = [x+1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
        if len(need_jacs) > 0:
            with Parallelizer.lookup(self.parallelizer) as par:
                new_jacs = [
                    x.squeeze() if isinstance(x, np.ndarray) else x
                    for x in intcds.jacobian(carts, need_jacs,
                                             mesh_spacing=self.internal_fd_mesh_spacing,
                                             stencil=self.internal_fd_stencil,
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
        ccoords = self.coords
        carts = ccoords.system
        internals = intcds.system
        if self.molecule not in self._cached_jacobians:
            self._cached_jacobians[self.molecule] = {
                'int': [],
                'cart': []
            }
        exist_jacs = self._cached_jacobians[self.molecule]['cart']
        max_jac = max(jacs)
        need_jacs = [x+1 for x in range(0, max_jac) if x >= len(exist_jacs) or exist_jacs[x] is None]
        if len(need_jacs) > 0:
            with Parallelizer.lookup(self.parallelizer) as par:
                new_jacs = [
                    x.squeeze() if isinstance(x, np.ndarray) else x
                    for x in ccoords.jacobian(internals, need_jacs,
                                                          mesh_spacing=self.cartesian_fd_mesh_spacing,
                                                          stencil=self.cartesian_fd_stencil,
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

    @property
    def inertial_frame(self):
        """
        Provides the inertial axis frame

        :return:
        :rtype:
        """

        if self._inert_frame is None:
            # Need to put B in Hartree?
            #  I've got moments of inertia in amu * bohr^2 at the moment
            #  So we convert (amu * bohr^2) to (m_e * bohr^2) since hb^2/(m_e bohr^2) == E_h
            mom_i, eigs = self.molecule.inertial_eigensystem
            B_e = 1 / (2 * mom_i * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass"))
            # print(B_e * UnitsData.convert("Hartrees", "Wavenumbers") )
            self._inert_frame = B_e, eigs

        return self._inert_frame

    def inertial_frame_derivatives(self):

        if self.strip_dummies:
            real_pos = self.masses > 0
            mass = self.masses[real_pos]
            crds = self.molecule.coords[real_pos, :]
        else:
            mass = self.masses.copy()
            mass[mass < 0] = self.zero_mass_term
            crds = self.molecule.coords
        mass = np.sqrt(mass)
        carts = mass[:, np.newaxis] * crds  # mass-weighted Cartesian coordinates

        ### compute basic inertia tensor derivatives
        # first derivs are computed as a full (nAt, 3, I_rows (3), I_cols (3)) tensor
        # and then reshaped to (nAt * 3, I_rows, I_cols)
        eyeXeye = np.eye(9).reshape(3, 3, 3, 3).transpose((2, 0, 1, 3))
        I0Y_1 = np.tensordot(carts, eyeXeye, axes=[1, 0])

        nAt = carts.shape[0]
        nY = nAt * 3
        I0Y_21 = (
                np.reshape(np.eye(3), (9,))[np.newaxis, :, np.newaxis]
                * carts[:, np.newaxis, :]
        ) # a flavor of outer product
        I0Y_21 = I0Y_21.reshape((nAt, 3, 3, 3))
        I0Y_2 = (I0Y_21 + I0Y_21.transpose((0, 1, 3, 2)))
        I0Y = 2 * I0Y_1 - I0Y_2
        I0Y = I0Y.reshape((nY, 3, 3))

        # second derivatives are 100% independent of coorinates
        # only the diagonal blocks are non-zero, so we compute that block
        # and then tile appropriately
        keyXey = np.eye(9).reshape(3, 3, 3, 3)
        I0YY_nn = 2 * eyeXeye - (keyXey + keyXey.transpose((0, 1, 3, 2)))
        I0YY = np.zeros((nAt, 3, nAt, 3, 3, 3))
        for n in range(nAt):
            I0YY[n, :, n, :, :, :] = I0YY_nn
        I0YY = I0YY.reshape((nY, nY, 3, 3))

        return [I0Y, I0YY]

    def moment_of_inertia_derivs(self, order):

        B_e, _ = self.inertial_frame
        YQ = self.modes.inverse  # derivatives of Q with respect to the Cartesians
        u_0 = 2 * np.diag(B_e) # reconstruct inertia tensor

        if order > 0:
            IdY, _ = self.inertial_frame_derivatives()  # only ever two of these
            IdQ = np.tensordot(YQ, IdY, axes=[1, 0])
            A = np.tensordot(IdQ, u_0, axes=[2, 0])

        all_derivs = [u_0]
        for i in range(order):
            # take original term and multiply in a . u_0
            u = all_derivs[-1]
            u = np.moveaxis(np.tensordot(u, A, axes=[1, 1]), -1, 1)
            all_derivs.append(u)

        # then add in the binomial expansion temrs
        for i,u in enumerate(all_derivs):
            all_derivs[i] = (-1)**i * (i+1) / (2**i) * u

        return all_derivs

    _cached_transforms = {}
    def get_coordinate_transforms(self,
                                  internal_by_cartesian_order=None,
                                  cartesian_by_internal_order=None,
                                  current_cache=None
                                  ):

        if internal_by_cartesian_order is None:
            internal_by_cartesian_order = self.internal_by_cartesian_order
        if cartesian_by_internal_order is None:
            cartesian_by_internal_order = self.cartesian_by_internal_order

        if current_cache is None and self.molecule in self._cached_transforms:
            current_cache = self._cached_transforms[self.molecule]

        if (
                current_cache is None
                or len(current_cache[JacobianKeys.CartesiansByInternals]) < cartesian_by_internal_order
                or len(current_cache[JacobianKeys.InternalsByCartesians]) < internal_by_cartesian_order
        ):

            if current_cache is None:
                current_cache = {}

            if self.logger is not None:
                self.logger.log_print(
                    [
                        "Getting coordinate transforms for {m}",
                        "Embedding axes: {a}"
                    ],
                    m=self.molecule,
                    a=self.internal_coordinates.system.converter_options["axes_labels"]
                )

            # fill out
            if (
                    JacobianKeys.CartesiansByInternals not in current_cache
                    or len(current_cache[JacobianKeys.CartesiansByInternals]) < cartesian_by_internal_order
            ):
                # For speed reasons we've introduced class-level caching of these terms
                if self.logger is not None:
                    start = time.time()
                    self.logger.log_print(
                        "Getting d^nX/dR^n up to order {o}...",
                        o=cartesian_by_internal_order
                    )
                cart_by_internal_jacobs = self.get_int_jacobs(list(range(1, cartesian_by_internal_order+1)))
                if self.logger is not None:
                    end = time.time()
                    self.logger.log_print(
                        "took {t}s",
                        t=round(end-start, 3)
                    )

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
                if self.strip_embedding:
                    embedding_coords = [0, 1, 2, 4, 5, 8]
                    good_coords = np.setdiff1d(np.arange(3*self.num_atoms), embedding_coords)

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
                        if self.strip_embedding:
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
                if self.logger is not None:
                    start = time.time()
                    self.logger.log_print(
                        "Getting d^nR/dX^n up to order {o}...",
                        o=internal_by_cartesian_order
                    )
                int_by_cartesian_jacobs = self.get_cart_jacobs(list(range(1, internal_by_cartesian_order + 1)))
                # m = np.max([np.max(np.abs(x)) for x in int_by_cartesian_jacobs])
                if self.logger is not None:
                    end = time.time()
                    self.logger.log_print(
                        "took {t}s",
                        t=round(end-start, 3)
                    )

                # raise Exception([x.shape for x in int_by_cartesian_jacobs])

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
                if self.strip_embedding:
                    embedding_coords = [0, 1, 2, 4, 5, 8]
                    good_coords = np.setdiff1d(np.arange(3*self.num_atoms), embedding_coords)

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
                        if self.strip_embedding:
                            x = np.take(x, good_coords, axis=-1)
                        _.append(x)
                    mc = np.expand_dims(mc, 0)
                    cartesian_weighting = cartesian_weighting * mc
                int_by_cartesian_jacobs = _

                current_cache[JacobianKeys.InternalsByCartesians] = int_by_cartesian_jacobs
            else:
                int_by_cartesian_jacobs = current_cache[JacobianKeys.InternalsByCartesians]

            if self._check_internal_modes():
                self.use_internal_modes = True

                QR = self.modes.matrix  # derivatives of Q with respect to the internals
                # we need to add zeros for the orientation coordinates
                if not self.strip_embedding and QR.shape[0] != 3*self.num_atoms:
                    _QR = QR
                    QR = np.zeros((3*self.num_atoms, _QR.shape[1]))
                    embedding_coords = [0, 1, 2, 4, 5, 8]
                    good_coords = np.setdiff1d(np.arange(3*self.num_atoms), embedding_coords)
                    QR[good_coords, :] = _QR
                    self.modes.matrix = QR

                RQ = self.modes.inverse # derivatives of internals with respect to Q
                if not self.strip_embedding and RQ.shape[1] != 3 * self.num_atoms:
                    _RQ = RQ
                    # we need to add zeros for the orientation coordinates
                    RQ = np.zeros((_RQ.shape[0], 3*self.num_atoms))
                    embedding_coords = [0, 1, 2, 4, 5, 8]
                    good_coords = np.setdiff1d(np.arange(3*self.num_atoms), embedding_coords)
                    RQ[:, good_coords] = _RQ
                    self.modes.inverse = RQ

                if (
                        JacobianKeys.CartesiansByInternalModes not in current_cache
                        or len(current_cache[JacobianKeys.CartesiansByInternalModes]) < len(cart_by_internal_jacobs)
                ):
                    x_derivs = cart_by_internal_jacobs#(YR, YRR, YRRR, YRRRR)
                    Q_derivs = [RQ] + [0]*(len(cart_by_internal_jacobs) - 1)
                    YQ_derivs = TensorDerivativeConverter(Q_derivs, x_derivs,
                                                          jacobians_name='Q',
                                                          values_name='X'
                                                          ).convert(order=len(cart_by_internal_jacobs))#, check_arrays=True)

                    current_cache[JacobianKeys.CartesiansByInternalModes] = YQ_derivs

                if (
                        JacobianKeys.InternalModesByCartesians not in current_cache
                        or len(current_cache[JacobianKeys.InternalModesByCartesians]) < len(int_by_cartesian_jacobs)
                ):
                    QY_derivs = TensorDerivativeConverter(int_by_cartesian_jacobs,
                                                          [QR] + [0]*(len(int_by_cartesian_jacobs) - 1)
                                                          ).convert(order=len(int_by_cartesian_jacobs))#, check_arrays=True)
                    current_cache[JacobianKeys.InternalModesByCartesians] = QY_derivs
            else:
                # tr_modes = self.molecule.translation_rotation_modes[1].T
                QY = self.modes.matrix  # derivatives of Q with respect to the Cartesians
                YQ = self.modes.inverse # derivatives of Cartesians with respect to Q
                # YQ = np.concatenate([
                #     tr_modes,
                #     YQ  # derivatives of Cartesians with respect to Q
                # ], axis=0)
                # QY = np.concatenate([
                #     tr_modes.T,
                #     QY  # derivatives of Cartesians with respect to Q
                # ], axis=1)

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

            self._cached_transforms[self.molecule] = current_cache
            with self.checkpointer:
                try:
                    self.checkpointer['coordinate_transforms'] = {k.value:v for k,v in current_cache.items()}
                except (OSError, KeyError):
                    pass

        return current_cache#self._cached_transforms[self.molecule]

    @property
    def cartesian_L_matrix(self):
        return self.get_cartesians_by_cartesian_modes(1)[0]
    def get_cartesians_by_cartesian_modes(self, order=None):
        # print(dict(
        #     cartesian_by_internal_order=order,
        #     internal_by_cartesian_order=min(order, self.internal_by_cartesian_order)
        # ))
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
        if strip_embedding and not self.strip_embedding:
            embedding_coords = [0, 1, 2, 4, 5, 8]
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
        if strip_embedding and not self.strip_embedding:
            embedding_coords = [0, 1, 2, 4, 5, 8]
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

        if strip_embedding and not self.strip_embedding:
            embedding_coords = [0, 1, 2, 4, 5, 8]
            good_coords = np.setdiff1d(np.arange(3 * self.num_atoms), embedding_coords)
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
        if strip_embedding and not self.strip_embedding:
            embedding_coords = [0, 1, 2, 4, 5, 8]
            good_coords = np.setdiff1d(np.arange(3 * self.num_atoms), embedding_coords)
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

class PotentialTerms(ExpansionTerms):
    """
    A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates
    """
    __props__ = ExpansionTerms.__props__ + (
        "potential_derivatives",
        "check_input_force_constants",
        "hessian_tolerance",
        "grad_tolerance",
        "freq_tolerance"
    )
    def __init__(self,
                 molecule,
                 mixed_derivs=None,
                 modes=None,
                 potential_derivatives=None,
                 mode_selection=None,
                 logger=None,
                 parallelizer=None,
                 checkpointer=None,
                 check_input_force_constants=True,
                 allow_higher_potential_terms=False,
                 hessian_tolerance=1.0e-4,
                 grad_tolerance=1.0e-4,
                 freq_tolerance=2e-3,
                 **opts
                 ):
        """
        :param molecule: the molecule that will supply the potential derivatives
        :type molecule: Molecule
        :param mixed_derivs: whether or not the pulled derivatives are partially derivatives along the normal coords
        :type mixed_derivs: bool
        :param modes: the normal modes to use when doing calculations
        :type modes: None | MolecularVibrations
        :param mode_selection: the subset of normal modes to use
        :type mode_selection: None | Iterable[int]
        """

        super().__init__(molecule, modes, mode_selection=mode_selection,
                         logger=logger, parallelizer=parallelizer, checkpointer=checkpointer,
                         **opts
                         )

        self.check_input_force_constants=check_input_force_constants
        self.hessian_tolerance = hessian_tolerance
        self.grad_tolerance = grad_tolerance
        self.freq_tolerance = freq_tolerance

        self.mixed_derivs = mixed_derivs # we can figure this out from the shape in the future
        self._input_derivs = potential_derivatives
        self._v_derivs = None #
        self.allow_higher_potential_terms=allow_higher_potential_terms

    @property
    def v_derivs(self):
        if self._v_derivs is None:
            if self._input_derivs is None:
                self._input_derivs = self.molecule.potential_surface.derivatives
            self._v_derivs = self._canonicalize_derivs(self.freqs, self.masses, self._input_derivs)

        return self._v_derivs
    @v_derivs.setter
    def v_derivs(self, v):
        self._v_derivs = v


    def _check_mode_terms(self, derivs=None):
        modes_n = len(self.modes.freqs)
        if derivs is None:
            derivs = self.v_derivs
        for d in derivs:
            if d.shape != (modes_n,) * len(d.shape):
                return False
        return True
    def _canonicalize_derivs(self, freqs, masses, derivs):

        if self._check_mode_terms(derivs):
            return derivs

        if len(derivs) == 3:
            grad, fcs, fds = derivs
            try:
                fcs = fcs.array
            except AttributeError:
                fcs, thirds, fourths = derivs
                grad = None
            else:
                thirds = fds.third_deriv_array
                fourths = fds.fourth_deriv_array
        elif len(derivs) == 4:
            grad, fcs, thirds, fourths = derivs
        else:
            grad = derivs[0]
            fcs = derivs[1]
            thirds = derivs[2] if len(derivs) > 2 else None
            fourths = derivs[3] if len(derivs) > 3 else None

        n = self.num_atoms
        modes_n = len(self.modes.freqs)
        internals_n = 3 * n - 6
        coord_n = 3 * n

        if len(derivs) > 2 and self.mode_sel is not None and thirds.shape[0] == internals_n:
            thirds = thirds[(self.mode_sel,)]

        if len(derivs) > 3 and self.mode_sel is not None and fourths.shape[0] == internals_n:
            if not isinstance(self.mode_sel, slice):
                fourths = fourths[np.ix_(self.mode_sel, self.mode_sel)]
            else:
                fourths = fourths[self.mode_sel, self.mode_sel]

        if grad is not None:
            if grad.shape != (coord_n,) and grad.shape != (internals_n,):
                raise PerturbationTheoryException(
                    "{0}.{1}: length of gradient array {2[0]} is not {3[0]} or {4[0]}".format(
                        type(self).__name__,
                        "_canonicalize_force_constants",
                        grad.shape,
                        (coord_n,),
                        (internals_n,)
                    )
                )
        if (
                fcs.shape != (coord_n, coord_n)
                and fcs.shape != (internals_n, internals_n)
                and fcs.shape != (modes_n, modes_n)
        ):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of force constant array ({2[0]}x{2[1]}) is not in ({3})".format(
                    type(self).__name__,
                    "_canonicalize_force_constants",
                    fcs.shape,
                    ", ".join("({0[0]}x{0[1]})".format(x) for x in [
                        (coord_n, coord_n),
                        (internals_n, internals_n),
                        (modes_n, modes_n)
                    ])
                )
            )

        if (    len(derivs) > 2
                    and thirds.shape != (modes_n, coord_n, coord_n)
                    and thirds.shape != (modes_n, internals_n, internals_n)
                    and thirds.shape != (modes_n, modes_n, modes_n)
                    and thirds.shape != (coord_n, coord_n, coord_n)
                    and thirds.shape != (internals_n, internals_n, internals_n)
        ):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of third derivative array ({2[0]}x{2[1]}x{2[2]}) is not in ({3})".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    thirds.shape,
                    ", ".join("({0[0]}x{0[1]}x{0[2]})".format(x) for x in [
                        (modes_n, coord_n, coord_n),
                        (modes_n, internals_n, internals_n),
                        (modes_n, modes_n, modes_n)
                    ])
                )
            )
        # this might need to change in the future
        if (
                len(derivs) > 3
                    and fourths.shape != (modes_n, modes_n, coord_n, coord_n)
                    and fourths.shape != (modes_n, modes_n, internals_n, internals_n)
                    and fourths.shape != (modes_n, modes_n, modes_n, modes_n)
                    and fourths.shape != (coord_n, coord_n, coord_n, coord_n)
                    and fourths.shape != (internals_n, internals_n, internals_n, internals_n)
        ):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of fourth derivative array ({2[0]}x{2[1]}x{2[2]}x{2[3]}) is not ({3})".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    fourths.shape,
                    ", ".join("({0[0]}x{0[1]}x{0[2]}x{0[3]})".format(x) for x in [
                        (modes_n, modes_n, coord_n, coord_n),
                        (modes_n, modes_n, internals_n, internals_n),
                        (modes_n, modes_n, modes_n, modes_n),
                        (coord_n, coord_n, coord_n, coord_n),
                        (internals_n, internals_n, internals_n, internals_n)
                    ])
                )
            )

        for i in range(4, len(derivs)):
            if (
                    derivs[i].shape != (coord_n,) * (i+1)
                    and derivs[i].shape != (internals_n,) * (i+1)
            ):
                raise PerturbationTheoryException(
                    "{0}.{1}: dimension of {2}th derivative array {3} is not ({4})".format(
                        type(self).__name__,
                        "_canonicalize_derivs",
                        i+1,
                        derivs[i].shape,
                        ", ".join(str(x) for x in [
                            (coord_n,) * (i + 1),
                            (internals_n,) * (i+1)
                        ])
                    )
                )

        if self._check_mode_terms(derivs) or self.use_internal_modes:
            all_derivs = derivs
        else:
            # amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
            m_conv = np.sqrt(self._tripmass(masses))
            f_conv = np.sqrt(freqs)
            # f_conv = np.ones(f_conv.shape) # debugging
            if fcs.shape == (coord_n, coord_n):
                undimension_2 = np.outer(m_conv, m_conv)
            elif fcs.shape == (modes_n, modes_n):
                undimension_2 = f_conv[:, np.newaxis] * f_conv[np.newaxis, :]
            else:
                undimension_2 = 1
            fcs = fcs * (1 / undimension_2)

            if self.freq_tolerance is not None and self.check_input_force_constants:
                xQ2 = self.modes.inverse
                _, v2x = TensorDerivativeConverter((xQ2, 0), (grad, fcs)).convert(order=2)

                real_freqs = np.diag(v2x)
                nominal_freqs = self.modes.freqs
                # deviation on the order of a wavenumber can happen in low-freq stuff from numerical shiz
                if self.freq_tolerance is not None:
                    if np.max(np.abs(nominal_freqs - real_freqs)) > self.freq_tolerance:
                        raise PerturbationTheoryException(
                            "Input frequencies aren't obtained when transforming the force constant matrix;"
                            " this likely indicates issues with the input mode vectors"
                            " got \n{}\n but expected \n{}\n".format(
                                real_freqs * UnitsData.convert("Hartrees", "Wavenumbers"),
                                nominal_freqs * UnitsData.convert("Hartrees", "Wavenumbers")
                            )
                        )

            all_derivs = [grad, fcs]
            if len(derivs) > 2:
                if thirds.shape == (modes_n, coord_n, coord_n):
                    if self.mixed_derivs is None:
                        self.mixed_derivs = True
                    undimension_3 = (
                            f_conv[:, np.newaxis, np.newaxis]
                            * m_conv[np.newaxis, :, np.newaxis]
                            * m_conv[np.newaxis, np.newaxis, :]
                    )
                elif thirds.shape == (coord_n, coord_n, coord_n):
                    if self.mixed_derivs is None:
                        self.mixed_derivs = False
                    undimension_3 = (
                            m_conv[:, np.newaxis, np.newaxis]
                            * m_conv[np.newaxis, :, np.newaxis]
                            * m_conv[np.newaxis, np.newaxis, :]
                    )
                elif thirds.shape == (modes_n, modes_n, modes_n):
                    if self.mixed_derivs is None:
                        self.mixed_derivs = False
                    undimension_3 = (
                            f_conv[:, np.newaxis, np.newaxis]
                            * f_conv[np.newaxis, :, np.newaxis]
                            * f_conv[np.newaxis, np.newaxis, :]
                    )
                else:
                    if self.mixed_derivs is None:
                        self.mixed_derivs = False
                    undimension_3 = 1
                thirds = thirds * (1 / undimension_3)
                all_derivs.append(thirds)

            if len(derivs) > 3:
                if fourths.shape == (modes_n, modes_n, coord_n, coord_n):
                    undimension_4 = (
                            f_conv[:, np.newaxis, np.newaxis, np.newaxis]
                            * f_conv[np.newaxis, :, np.newaxis, np.newaxis]
                            * m_conv[np.newaxis, np.newaxis, :, np.newaxis]
                            * m_conv[np.newaxis, np.newaxis, np.newaxis, :]
                    )
                elif fourths.shape == (coord_n, coord_n, coord_n, coord_n):
                    undimension_4 = (
                            m_conv[:, np.newaxis, np.newaxis, np.newaxis]
                            * m_conv[np.newaxis, :, np.newaxis, np.newaxis]
                            * m_conv[np.newaxis, np.newaxis, :, np.newaxis]
                            * m_conv[np.newaxis, np.newaxis, np.newaxis, :]
                    )
                elif fourths.shape == (modes_n, modes_n, modes_n, modes_n):
                    undimension_4 = (
                            f_conv[:, np.newaxis, np.newaxis, np.newaxis]
                            * f_conv[np.newaxis, :, np.newaxis, np.newaxis]
                            * f_conv[np.newaxis, np.newaxis, :, np.newaxis]
                            * f_conv[np.newaxis, np.newaxis, np.newaxis, :]
                    )
                else:
                    undimension_4 = 1

                if isinstance(fourths, SparseArray):
                    fourths = fourths.asarray()
                fourths = fourths * (1 / undimension_4)

                all_derivs.append(fourths)

            for i in range(4, len(derivs)):
                term = derivs[i]
                if term.shape == (coord_n,) * (i + 1):
                    undimension = m_conv
                    mc = m_conv
                    for j in range(i):
                        mc = np.expand_dims(mc, 0)
                        undimension = np.expand_dims(undimension, -1) * mc
                elif term.shape != (internals_n,) * (i + 1):
                    undimension = f_conv
                    fc = f_conv
                    for j in range(i):
                        fc = np.expand_dims(fc, 0)
                        undimension = np.expand_dims(undimension, -1) * fc

                all_derivs.append(term / undimension)

        return all_derivs

    def get_terms(self, order=None, logger=None):

        if self._check_mode_terms():
            return self.v_derivs[1:]

        if logger is None:
            logger = self.logger

        logger.log_print('calculating potential derivatives')

        if order is None:
            order = len(self.v_derivs)
        else:
            order += 2

        if self.allow_higher_potential_terms and len(self.v_derivs) < order:
            self.v_derivs = tuple(self.v_derivs) + (0,) * order

        grad = self.v_derivs[0]
        if self.grad_tolerance is not None:
            if np.linalg.norm(grad) > self.grad_tolerance:
                # add some logger stuff...
                logger.log_print(
                    "WARNING: gradient norm is {n}",
                    n=np.linalg.norm(grad)
                )
        grad = np.zeros(grad.shape)
        V_derivs = [grad] + list(self.v_derivs[1:])
        # hess = self.v_derivs[1]

        # raise Exception([x.shape for x in self.v_derivs])

        # Use the Molecule's coordinates which know about their embedding by default
        intcds = self.internal_coordinates
        direct_prop = (
                self.direct_propagate_cartesians
                and not (isinstance(self.direct_propagate_cartesians, str) and self.direct_propagate_cartesians == 'dipoles')
        )
        mixed_derivs = self.mixed_derivs
        if intcds is None or direct_prop:
            # this is nice because it eliminates most of the terms in the expansion
            xQ = self.modes.inverse
            x_derivs = [xQ] + [0] * (order-1)
            # terms = self._get_tensor_derivs(x_derivs, V_derivs, mixed_terms=False, mixed_XQ=self.mixed_derivs)
            if mixed_derivs:
                terms = TensorDerivativeConverter(x_derivs, V_derivs, mixed_terms=[
                    [None, v] for v in V_derivs[2:]
                ]).convert(order=order)#, check_arrays=True)
            else:
                terms = TensorDerivativeConverter(x_derivs, V_derivs).convert(order=order)#, check_arrays=True)
        elif self._check_internal_modes() and not self._check_mode_terms():
            raise NotImplementedError("...")
            # It should be very rare that we are actually able to make it here
            terms = []
            RQ = self.modes.inverse
            for v in V_derivs:
                for j in range(v.ndim):
                    v = np.tensordot(RQ, v, axes=[1, -1])
                terms.append(v)
        else:
            if ( # handle mixed derivative resymmetrization
                    mixed_derivs
                    and self.mixed_derivative_handling_mode != MixedDerivativeHandlingModes.Unhandled
                    and order > 3
            ):
                ## TODO: figure out a way to make this work
                raise NotImplementedError("haven't included translation/rotation modes needed to make this work correctly")

                # raise NotImplementedError("different methods for handling mixed derivatives need patching")
                xQ = self.matrix.inverse
                x_derivs = [xQ] + [0] * (order - 1)
                cart_terms = TensorDerivativeConverter(x_derivs, V_derivs, mixed_terms=[
                    [None, v] for v in V_derivs[2:]
                ]).convert(order=order)
                # provides the derivatives expressed with respect to the Cartesian normal modes
                v3 = cart_terms[2]
                v4 = cart_terms[3]
                # transform, resymmetrize, and then go to internals
                for i in range(v4.shape[0]):
                    if (
                            self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Numerical
                            or self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Unhandled
                    ):
                        v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]
                    elif self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Analytical:
                        r = range(i, v4.shape[0])
                        v4[i, i, r, r] = v4[r, r, i, i]
                    elif self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Averaged:
                        r = range(i, v4.shape[0])
                        v4[i, i, r, r] = np.average(
                            [
                                v4[r, r, i, i],
                                v4[i, i, r, r]
                                ],
                            axis=0
                        )
                        # v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]
                    else:
                        raise ValueError("don't know what to do with `mixed_derivative_handling_mode` {} ".format(self.mixed_derivative_handling_mode))
                for i in range(v4.shape[0]):
                    for j in range(i+1, v4.shape[0]):
                        for k in range(j+1, v4.shape[0]):
                            for l in range(k+1, v4.shape[0]):
                                # if (i != j and i != k and i != l and j != k and j != l and k != l ): # all different
                                for p in itertools.permutations([i, j, k, l]):
                                    v4[p] = 0
                for i in range(v4.shape[0]):
                    v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]

                qQ_derivs = self.get_cartesian_modes_by_internal_modes(len(cart_terms)-1)
                terms = TensorDerivativeConverter(
                    qQ_derivs + [0],  # pad for the zeroed out gradient term
                    cart_terms
                ).convert(order=order)
                v4 = terms[3]
                for i in range(v4.shape[0]):
                    for j in range(i+1, v4.shape[0]):
                        for k in range(j+1, v4.shape[0]):
                            for l in range(k+1, v4.shape[0]):
                                # if (i != j and i != k and i != l and j != k and j != l and k != l ): # all different
                                for p in itertools.permutations([i, j, k, l]):
                                    v4[p] = 0
                for i in range(v4.shape[0]):
                    v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]
                terms[3] = v4
                mixed_derivs = False
                #   0 0 0 0 0 0     80.93408   -156.52080
                #   0 0 0 0 0 1    231.69991   -345.26155
                #
                #   0 0 0 0 0 0     89.28303   -156.52080
                #   0 0 0 0 0 1    255.14930   -345.26168
                #
                #   0 0 0 0 0 0    109.70479   -156.52080
                #   0 0 0 0 0 1    265.62694   -345.26168
            else:

                x_derivs = self.get_cartesians_by_modes(order=order-1)
                # raise Exception(x_derivs[1])
                x_derivs = list(x_derivs) + [0] # gradient term never matters

                if mixed_derivs:
                    if order > 4 and not self.allow_higher_potential_terms:
                        raise ValueError("don't currently have things tested for expansions beyond 4th V derivatives with mixed derivatives") #TODO: relax this once we have more flexible input determination
                    # terms = self._get_tensor_derivs(x_derivs, V_derivs, mixed_terms=True, mixed_XQ=self.mixed_derivs)

                    # since the normal modes are expressed over
                    # different sets of coordinates the fourth mixed deriv
                    # terms need to be partially corrected
                    qQ, qQQ = self.get_cartesian_modes_by_internal_modes(2)
                    f43 = np.tensordot(qQQ, V_derivs[2], axes=[2, 0])
                    fourths = V_derivs[3] + f43
                    V_derivs = V_derivs[:3] + [fourths] + V_derivs[4:]

                    terms = TensorDerivativeConverter(x_derivs, V_derivs,
                                                      mixed_terms=[
                                                          [None, v] for v in V_derivs[2:]
                                                      ]
                                                      ).convert(order=order)  # , check_arrays=True)
                else:
                    terms = TensorDerivativeConverter(x_derivs, V_derivs).convert(order=order)#, check_arrays=True)

            if self.hessian_tolerance is not None:
                    xQ2 = self.modes.inverse
                    _, v2x, = TensorDerivativeConverter((xQ2, 0), V_derivs).convert(order=2)
                    v2 = terms[1]
                    v2_diff = v2 - v2x
                    if np.max(np.abs(v2_diff)) > self.hessian_tolerance:
                        new_freqs = np.diag(v2)*UnitsData.convert("Hartrees", "Wavenumbers")
                        old_freqs = np.diag(v2x)*UnitsData.convert("Hartrees", "Wavenumbers")
                        zero_pos_new = np.where(np.abs(new_freqs) < 1.0e-10)
                        zero_pos_old = np.where(np.abs(old_freqs) < 1.0e-10)
                        if len(zero_pos_new) > 0 and (
                            len(zero_pos_old) == 0
                            or len(zero_pos_old[0]) != len(zero_pos_new[0])
                        ):
                            raise PerturbationTheoryException(
                                (
                                    "Encountered zero frequencies in internal normal mode Hessian that aren't in Cartesian normal mode Hessian."
                                    " Cartesian frequencies are \n{}\n but internals are \n{}\n"
                                    " This often indicates a planar dihedral angle where the derivatives are ill-defined.\n"
                                    " Try using dummy atoms to create a proper 3D structure.\n"
                                ).format(
                                    old_freqs,
                                    new_freqs
                                )
                            )
                        else:
                            raise PerturbationTheoryException(
                                (
                                    "Internal normal mode Hessian differs from Cartesian normal mode Hessian."
                                    " Cartesian frequencies are \n{}\n, internals are \n{}\n"
                                    " This often indicates issues with the derivatives.\n"
                                    " (YQ min/max: {} {} generally in the 10s for well-behaved systems)\n"
                                    " (YQQ min/max: {} {} generally in the 10s for well-behaved systems)"
                                 ).format(
                                    old_freqs,
                                    new_freqs,
                                    np.min(x_derivs[0]), np.max(x_derivs[0]),
                                    np.min(x_derivs[1]), np.max(x_derivs[1])
                                )
                            )

        if mixed_derivs and order > 3:
            v4 = terms[3]
            # we assume we only got second derivs in Q_i Q_i
            # at this point, then, we should be able to fill in the terms we know are missing
            if not isinstance(v4, np.ndarray):
                v4 = v4.asarray()
            for i in range(v4.shape[0]):
                for j in range(i+1, v4.shape[0]):
                    for k in range(j+1, v4.shape[0]):
                        for l in range(k+1, v4.shape[0]):
                            # if (i != j and i != k and i != l and j != k and j != l and k != l ): # all different
                            for p in itertools.permutations([i, j, k, l]):
                                v4[p] = 0
                            # v4[i, j, k, l] = 0
                            # v4[i, j, k, l] = 0
                            # v4[i, j, k, l] = 0
                            # v4[i, j, k, l] = 0
            for i in range(v4.shape[0]):
                v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]
            terms[3] = v4

        if intcds is not None and self.backpropagate_internals:
            # need to internal mode terms and
            # convert them back to Cartesian mode ones...
            Qq_derivs = self.get_internal_modes_by_cartesian_modes(len(terms) - 1)
            terms = TensorDerivativeConverter(
                Qq_derivs + [0], # pad for the zeroed out gradient term
                terms
            ).convert(order=order)
        elif intcds is not None and direct_prop:
            # need to internal mode terms and
            # convert them back to Cartesian mode ones...
            qQ_derivs = self.get_cartesian_modes_by_internal_modes(len(terms) - 1)
            terms = TensorDerivativeConverter(
                qQ_derivs + [0],  # pad for the zeroed out gradient term
                terms
            ).convert(order=order)

            if self.hessian_tolerance is not None:
                xQ2 = self.modes.inverse
                _, v2x, = TensorDerivativeConverter((xQ2, 0), V_derivs).convert(order=2)
                v2 = terms[1]
                v2_diff = v2 - v2x
                if np.max(np.abs(v2_diff)) > self.hessian_tolerance:
                    new_freqs = np.diag(v2) * UnitsData.convert("Hartrees", "Wavenumbers")
                    old_freqs = np.diag(v2x) * UnitsData.convert("Hartrees", "Wavenumbers")
                    zero_pos_new = np.where(np.abs(new_freqs) < 1.0e-10)
                    zero_pos_old = np.where(np.abs(old_freqs) < 1.0e-10)
                    if len(zero_pos_new) > 0 and (
                            len(zero_pos_old) == 0
                            or len(zero_pos_old[0]) != len(zero_pos_new[0])
                    ):
                        raise PerturbationTheoryException(
                            (
                                "Encountered zero frequencies in internal normal mode Hessian that aren't in Cartesian normal mode Hessian."
                                " Cartesian frequencies are \n{}\n but internals are \n{}\n"
                                " This often indicates a planar dihedral angle where the derivatives are ill-defined.\n"
                                " Try using dummy atoms to create a proper 3D structure.\n"
                            ).format(
                                old_freqs,
                                new_freqs
                            )
                        )
                    else:
                        raise PerturbationTheoryException(
                            (
                                "Internal normal mode Hessian differs from Cartesian normal mode Hessian."
                                " Cartesian frequencies are \n{}\n, internals are \n{}\n"
                                " This often indicates issues with the derivatives.\n"
                                " (YQ min/max: {} {} generally in the 10s for well-behaved systems)\n"
                                " (YQQ min/max: {} {} generally in the 10s for well-behaved systems)"
                            ).format(
                                old_freqs,
                                new_freqs,
                                np.min(x_derivs[0]), np.max(x_derivs[0]),
                                np.min(x_derivs[1]), np.max(x_derivs[1])
                            )
                        )

        # drop the gradient term as that is all zeros
        terms = terms[1:]

        try:
            self.checkpointer['potential_terms'] = terms
        except (OSError, KeyError):
            pass

        if self.hessian_tolerance is not None and np.linalg.norm(np.abs(terms[0] - np.diag(np.diag(terms[0])))) > self.hessian_tolerance:
            raise ValueError("F-matrix is not diagonal (got {})".format(terms[0]))

        new_freqs = np.diag(terms[0])
        old_freqs = self.modes.freqs
        # deviation on the order of a wavenumber can happen in low-freq stuff from numerical shiz
        if self.freq_tolerance is not None:
            if np.max(np.abs(new_freqs - old_freqs)) > self.freq_tolerance:
                raise PerturbationTheoryException(
                    "Force constants in normal modes don't return frequencies along diagonal;"
                    " this likely indicates issues with the mass-weighting"
                    " got \n{}\n but expected \n{}\n".format(
                        new_freqs*UnitsData.convert("Hartrees", "Wavenumbers"),
                        old_freqs*UnitsData.convert("Hartrees", "Wavenumbers")
                    )
                )

        return terms

class KineticTerms(ExpansionTerms):
    """Represents the KE coefficients"""

    __props__ = ExpansionTerms.__props__ + (
        'g_derivative_threshold',
        "gmatrix_tolerance",
        'use_cartesian_kinetic_energy'
    )
    def __init__(self,
                 molecule,
                 g_derivative_threshold=1e-3,
                 gmatrix_tolerance=1e-6,
                 use_cartesian_kinetic_energy=False,
                 **opts
                 ):
        super().__init__(molecule, **opts)
        self.g_derivative_threshold = g_derivative_threshold
        self.gmatrix_tolerance = gmatrix_tolerance
        self.use_cartesian_kinetic_energy = use_cartesian_kinetic_energy

    def get_terms(self, order=None, logger=None):

        if logger is None:
            logger = self.logger
        logger.log_print('calculating G-matrix derivatives')

        dot = DumbTensor._dot
        shift = DumbTensor._shift
        intcds = self.internal_coordinates
        if self.use_cartesian_kinetic_energy or intcds is None or self.backpropagate_internals:
            # this is nice because it eliminates a lot of terms in the expansion
            J = self.modes.matrix
            G = dot(J, J, axes=[[0, 0]])
            if order == 0:
                terms = [G]
            else:
                terms = [G] + [0]*(order)
        else:
            # should work this into the new layout
            uses_internal_modes = self._check_internal_modes()
            if uses_internal_modes:
                QY_derivs = self.get_internals_by_cartesians(order=order + 1) # really dRdY derivatives
                YQ_derivs = self.get_cartesians_by_internals(order=order + 1) # really dYdR derivatives
            else:
                QY_derivs = self.get_modes_by_cartesians(order=order+1)
                YQ_derivs = self.get_cartesians_by_modes(order=order+1)
            # RQ = dot(YQ, RY)

            term_getter = TensorDerivativeConverter(YQ_derivs, QY_derivs).terms
            term_getter.v_name = 'Y'
            J = term_getter.XV(1)
            G_terms = [J.dot(J, 1, 1)]
            for i in range(1, order+1):
                g_cur = G_terms[-1].dQ()#.simplify()
                G_terms.append(g_cur)
            terms = [x.array for x in G_terms]

            if uses_internal_modes:
                QR = self.modes.matrix
                RQ = self.modes.inverse
                for i,g in enumerate(terms):
                    for j in range(2):
                        g = np.tensordot(QR, g, axes=[0, -1])
                    for j in range(i):
                        g = np.tensordot(RQ, g, axes=[1, -1])
                    terms[i] = g


            for i,t in enumerate(terms):
                if i == 0:
                    if self.gmatrix_tolerance is not None and np.linalg.norm(np.abs(terms[0] - np.diag(np.diag(terms[0])))) > self.gmatrix_tolerance:
                        raise ValueError("G-matrix is not diagonal (got {})".format(terms[0]))
                    continue
                m = np.max(np.abs(t))
                if m > self.g_derivative_threshold:
                    # print(G_terms[i])
                    self.logger.log_print("WARNING: max of d^{i}G/dQ^{i} is {m}", i=i, m=m)

            # QY, QYY, QYYY = QY_derivs[:3]
            # YQ, YQQ, YQQQ = YQ_derivs[:3]
            # G = dot(QY, QY, axes=[[0, 0]])
            #
            # GQ = dot(YQ, QYY, QY, axes=[[1, 0], [1, 0]]) + dot(YQ, dot(QY, QYY, axes=[[0, 0]]), axes=[[-1, 1]])
            # GQQ = (
            #         dot(YQQ, QYY, QY, axes=[[-1, 0], [2, 0]])
            #         + dot(YQQ, dot(QY, QYY, axes=[[0, 0]]), axes=[[-1, 1]])
            #         + dot(YQ, dot(YQ, QYYY, QY, axes=[[-1, 0], [1, 0]]), axes=[[1, 1]])
            #         + dot(YQ, dot(YQ, dot(QY, QYYY, axes=[[0, 0]]), axes=[[1, 1]]), axes=[[1, 2]])
            #         + dot(YQ, dot(YQ, QYY, QYY, axes=[[-1, 0], [1, 1]]), axes=[[1, 2]])
            #         + dot(YQ, dot(YQ, QYY, QYY, axes=[[-1, 0], [1, 1]]), axes=[[1, 2]]).transpose((0, 1, 3, 2))
            #         # + dot(YQ, dot(YQ, dot(QYY, QYY, axes=[[0, 0]]), axes=[[-1, 0]]), axes=[[1, 2]])
            # )

        G_terms = terms
        try:
            self.checkpointer['gmatrix_terms'] = G_terms
        except (OSError, KeyError):
            pass

        return G_terms

class DipoleTerms(ExpansionTerms):
    __props__ = ExpansionTerms.__props__ + (
        "dipole_derivatives",
    )
    def __init__(self,
                 molecule,
                 dipole_derivatives=None,
                 mixed_derivs=None,
                 modes=None,
                 mode_selection=None,
                 logger=None,
                 parallelizer=None,
                 checkpointer=None,
                 **opts
                 ):
        """
        :param molecule: the molecule that will supply the dipole derivatives
        :type molecule: Molecule
        :param mixed_derivs: whether or not the pulled derivatives are partially derivatives along the normal coords
        :type mixed_derivs: bool
        :param modes: the normal modes to use when doing calculations
        :type modes: None | MolecularVibrations
        :param mode_selection: the subset of normal modes to use
        :type mode_selection: None | Iterable[int]
        """
        self.derivs = None
        super().__init__(molecule, modes=modes, mode_selection=mode_selection,
                         logger=logger, parallelizer=parallelizer, checkpointer=checkpointer,
                         **opts
                         )
        self.mixed_derivs = mixed_derivs
        if self.mixed_derivs is None:
            self.mixed_derivs = mixed_derivs
        if dipole_derivatives is None:
            dipole_derivatives = molecule.dipole_surface.derivatives
        self.derivs = self._canonicalize_derivs(self.freqs, self.masses, dipole_derivatives)

    def _canonicalize_derivs(self, freqs, masses, derivs):
        """
        Makes sure all of the dipole moments are clean and ready to rotate
        """

        if self._check_mode_terms(derivs):
            return derivs

        if len(derivs) == 4:
            mom, grad, seconds, thirds = derivs
            try:
                grad = grad.array
            except AttributeError:
                pass
        else:
            mom = derivs[0]
            grad = derivs[1]
            seconds = derivs[2] if len(derivs) > 2 else None
            thirds = derivs[3] if len(derivs) > 3 else None

        n = len(masses)
        modes_n = len(self.modes.freqs)
        internals_n = 3 * n - 6
        coord_n = 3 * n

        if len(derivs) > 2 and self.mode_sel is not None and seconds.shape[0] == internals_n:
            seconds = seconds[(self.mode_sel,)]

        if self.mode_sel is not None and thirds.shape[0] == internals_n:
            if not isinstance(self.mode_sel, slice):
                thirds = thirds[np.ix_(self.mode_sel, self.mode_sel)]
            else:
                thirds = thirds[self.mode_sel, self.mode_sel]

        if grad is not None:
            if (
                    grad.shape != (coord_n, 3)
                    and grad.shape != (internals_n, 3)
            ):
                raise PerturbationTheoryException(
                    "{0}.{1}: dimension of dipole derivative array ({2[0]}) is not {3[0]} or {4[0]}".format(
                        type(self).__name__,
                        "_canonicalize_derivs",
                        grad.shape,
                        (coord_n, 3),
                        (internals_n, 3)
                    )
                )

        if (
                len(derivs) > 2
                and seconds.shape != (modes_n, coord_n, 3)
                and seconds.shape != (modes_n, internals_n, 3)
                and seconds.shape != (modes_n, modes_n, 3)
                and seconds.shape != (coord_n, coord_n, 3)
                and seconds.shape != (internals_n, internals_n, 3)
        ):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of dipole second derivative array ({2[0]}x{2[1]}) not in {3}".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    seconds.shape,
                    ", ".join("({0[0]}x{0[1]})".format(x) for x in [
                        (modes_n, coord_n, 3),
                        (modes_n, internals_n, 3),
                        (modes_n, modes_n, 3),
                        (coord_n, coord_n, 3),
                        (internals_n, internals_n, 3)
                    ])
                )
            )

        if (

                len(derivs) > 3
                and thirds.shape != (modes_n, modes_n, coord_n, 3)
                and thirds.shape != (modes_n, modes_n, internals_n, 3)
                and thirds.shape != (modes_n, modes_n, modes_n, 3)
                and thirds.shape != (coord_n, coord_n, coord_n, 3)
                and thirds.shape != (internals_n, internals_n, internals_n, 3)
        ):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of dipole third derivative array ({2[0]}x{2[1]}x{2[2]}) not in ({3})".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    thirds.shape,
                    ", ".join("({0[0]}x{0[1]}x{0[2]})".format(x) for x in [
                        (modes_n, modes_n, coord_n, 3),
                        (modes_n, modes_n, internals_n, 3),
                        (modes_n, modes_n, modes_n, 3),
                        (coord_n, coord_n, coord_n, 3),
                        (internals_n, internals_n, internals_n, 3)
                    ])
                )
            )

        for i in range(4, len(derivs)):
            if (
                    derivs[i].shape != (coord_n,) * (i+1) + (3,)
                    and derivs[i].shape != (internals_n,) * (i+1) + (3,)
            ):
                raise PerturbationTheoryException(
                    "{0}.{1}: dimension of {2}th dipole derivative array {3} is not ({4})".format(
                        type(self).__name__,
                        "_canonicalize_derivs",
                        i+1,
                        derivs[i].shape,
                        ", ".join(str(x) for x in [
                            (coord_n,) * (i + 1) + (3,),
                            (internals_n,) * (i+1) + (3,)
                        ])
                    )
                )

        # We need to mass-weight the pure cartesian derivs
        # & undimensionalize the ones in terms of normal modes

        # amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        m_conv = np.sqrt(self._tripmass(masses))
        f_conv = np.sqrt(freqs)

        if grad.shape == (coord_n, 3):
            grad = grad / m_conv[:, np.newaxis]
        elif grad.shape == (modes_n, 3):
            grad = grad / f_conv[:, np.newaxis]

        all_derivs = [mom, grad]
        if len(derivs) > 2:
            if seconds.shape == (modes_n, coord_n, 3):
                if self.mixed_derivs is None:
                    self.mixed_derivs = True
                undimension_2 = (
                        f_conv[:, np.newaxis, np.newaxis]
                        * m_conv[np.newaxis, :, np.newaxis]
                )
            elif seconds.shape == (coord_n, coord_n, 3):
                if self.mixed_derivs is None:
                    self.mixed_derivs = False
                undimension_2 = (
                        m_conv[:, np.newaxis, np.newaxis]
                        * m_conv[np.newaxis, :, np.newaxis]
                        * m_conv[np.newaxis, np.newaxis, :]
                )
            else:
                undimension_2 = 1
            seconds = seconds / undimension_2
            all_derivs.append(seconds)

        if len(derivs) > 3:
            if thirds.shape == (modes_n, modes_n, coord_n, 3):
                if self.mixed_derivs is None:
                    self.mixed_derivs = True
                undimension_3 = (
                        f_conv[:, np.newaxis, np.newaxis, np.newaxis]
                        * f_conv[np.newaxis, :, np.newaxis, np.newaxis]
                        * m_conv[np.newaxis, np.newaxis, :, np.newaxis]
                )
            elif thirds.shape == (coord_n, coord_n, coord_n, 3):
                if self.mixed_derivs is None:
                    self.mixed_derivs = False
                undimension_3 = (
                        m_conv[:, np.newaxis, np.newaxis]
                        * m_conv[np.newaxis, :, np.newaxis]
                        * m_conv[np.newaxis, np.newaxis, :]
                )
            elif thirds.shape == (modes_n, modes_n, modes_n, 3):
                if self.mixed_derivs is None:
                    self.mixed_derivs = False
                undimension_3 = (
                        f_conv[:, np.newaxis, np.newaxis]
                        * f_conv[np.newaxis, :, np.newaxis]
                        * f_conv[np.newaxis, np.newaxis, :]
                )
            else:
                if self.mixed_derivs is None:
                    self.mixed_derivs = False
                undimension_3 = 1
            thirds = thirds / undimension_3
            all_derivs.append(thirds)

        for i in range(4, len(derivs)):
            term = derivs[i]
            if term.shape == (coord_n,) * (i + 1) + (3,):
                undimension = m_conv
                mc = m_conv
                for j in range(i):
                    mc = np.expand_dims(mc, 0)
                    undimension = np.expand_dims(undimension, -1) * mc
            elif term.shape != (internals_n,) * (i + 1) + (3,):
                undimension = f_conv
                fc = f_conv
                for j in range(i):
                    fc = np.expand_dims(fc, 0)
                    undimension = np.expand_dims(undimension, -1) * fc
            undimension = np.expand_dims(undimension, -1) # for the three components
            all_derivs.append(term / undimension)

        return all_derivs

    def _check_mode_terms(self, derivs=None):
        modes_n = len(self.modes.freqs)
        if derivs is None:
            derivs = self.derivs[1:]
        for d in derivs:
            if d.shape != (modes_n,) * len(d.shape):
                return False
        return True
    def get_terms(self, order=None):

        if self._check_mode_terms():
            return self.derivs[1:]

        if order is None:
            order = len(self.derivs) - 1
        else:
            order += 1

        v0 = self.derivs[0]
        # grad = self.derivs[1]
        # seconds = self.derivs[2]
        # thirds = self.derivs[3]

        # Use the Molecule's coordinates which know about their embedding by default
        intcds = self.internal_coordinates
        direct_prop = (
                self.direct_propagate_cartesians
                and not (isinstance(self.direct_propagate_cartesians, str) and self.direct_propagate_cartesians == 'potential')
        )
        mixed_derivs = self.mixed_derivs
        if intcds is None or direct_prop:# or not self.non_degenerate:
            # this is nice because it eliminates most of terms in the expansion
            xQ = self.modes.inverse
            x_derivs = [xQ] + [0] * (order-1)
            mu_derivs = self.derivs[1:]
        else:
            x_derivs = self.get_cartesians_by_modes(order=order)

            mu_derivs = self.derivs[1:]

            if len(mu_derivs) > 2:
                qQ, qQQ = self.get_cartesian_modes_by_internal_modes(2)
                f43 = np.tensordot(qQQ, mu_derivs[1], axes=[2, 0])
                mu_derivs = list(mu_derivs)
                mu_derivs[2] = mu_derivs[2] + f43

        mu = [None]*3
        for coord in range(3):

            u_derivs = [d[..., coord] for d in mu_derivs]
            mixed_terms = [
                [u_derivs[1]],  # dVdQXX
                [u_derivs[2]]  # dVdQQXX
            ]

            if mixed_derivs:
                if (  # handle mixed derivative resymmetrization
                        self.mixed_derivative_handling_mode != MixedDerivativeHandlingModes.Unhandled
                        and order > 1
                ):
                    # d^2X/dQ^2@dU/dX + dX/dQ@dU/dQdX
                    xQ = self.modes.inverse
                    v1, v2, v3 = TensorDerivativeConverter(
                        [xQ] + [0] * (order-1),
                        u_derivs,
                        mixed_terms=mixed_terms,
                        values_name="U"
                    ).convert(order=3)  # , check_arrays=True)
                    if self.mixed_derivs:  # and intcds is None:
                        # Gaussian gives slightly different constants
                        # depending on whether the analytic or numerical derivs
                        # were transformed
                        for i in range(v2.shape[0]):
                            if self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Numerical:
                                v2[i, :] = v2[:, i] = v2[i, :]
                            elif self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Analytical:
                                v2[i, :] = v2[:, i] = v2[:, i]
                            elif self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Averaged:
                                v2[i, :] = v2[:, i] = np.average(
                                    [
                                        v2[i, :], v2[:, i]
                                    ],
                                    axis=0
                                )
                            else:
                                raise ValueError(
                                    "don't know what to do with `mixed_derivative_handling_mode` {} ".format(
                                        self.mixed_derivative_handling_mode
                                    )
                                )
                    Qx = self.modes.matrix
                    v1, v2, v3 = TensorDerivativeConverter(
                        [Qx] + [0] * (order - 1),
                        [v1, v2, v3],
                        values_name="U"
                    ).convert(order=3)  # , check_arrays=True)
                    u_derivs = [v1, v2, v3]
                    mixed_terms = None

                # d^2X/dQ^2@dU/dX + dX/dQ@dU/dQdX
                terms = TensorDerivativeConverter(x_derivs, u_derivs,
                                                  mixed_terms=mixed_terms,
                                                  values_name="U"
                                                  ).convert(order=order)  # , check_arrays=True)
                # terms = list(terms)
                # if order > 1:
                #     v2 = terms[1]
                #     if self.mixed_derivs:  # and intcds is None:
                #         # Gaussian gives slightly different constants
                #         # depending on whether the analytic or numerical derivs
                #         # were transformed
                #         if self.mixed_derivative_handling_mode != MixedDerivativeHandlingModes.Unhandled:
                #             for i in range(v2.shape[0]):
                #                 if self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Numerical:
                #                     v2[i, :] = v2[:, i] = v2[i, :]
                #                 elif self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Analytical:
                #                     v2[i, :] = v2[:, i] = v2[:, i]
                #                 elif self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Averaged:
                #                     v2[i, :] = v2[:, i] = np.average(
                #                         [
                #                             v2[i, :], v2[:, i]
                #                         ],
                #                         axis=0
                #                     )
                #                 else:
                #                     raise ValueError(
                #                         "don't know what to do with `mixed_derivative_handling_mode` {} ".format(
                #                             self.mixed_derivative_handling_mode
                #                         )
                #                     )

                if order > 2:
                    v3 = terms[2]
                    if self.mixed_derivs:  # and intcds is None:
                        # Gaussian gives slightly different constants
                        # depending on whether the analytic or numerical derivs
                        # were transformed
                        # we assume we only got second derivs in Q_i Q_i
                        # at this point, then, we should be able to fill in the terms we know are missing
                        if not isinstance(v3, np.ndarray):
                            v3 = v3.asarray()
                            terms[2] = v3
                        for i in range(v3.shape[0]):
                            for j in range(i + 1, v3.shape[0]):
                                for k in range(j + 1, v3.shape[0]):
                                    # if (i != j and i != k and i != l and j != k and j != l and k != l ): # all different
                                    for p in itertools.permutations([i, j, k]):
                                        v3[p] = 0

                        # if self.mixed_derivative_handling_mode != MixedDerivativeHandlingModes.Unhandled:
                        #     for i in range(v3.shape[0]):
                        #         if self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Numerical:
                        #             r = np.arange(i, v3.shape[0])
                        #             v3[i, i, r] = v3[i, r, i] = v3[r, i, i] = v3[i, i, r]
                        #         elif self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Analytical:
                        #             r = np.arange(i, v3.shape[0])
                        #             # v3[i, :, :] = v3[:, i, :] = v3[:, :, i] = v3[i, :, :]
                        #             v3[i, i, r] = v3[i, r, i] = v3[r, i, i] = v3[r, i, i]
                        #             # v3[i, :, :] = v3[:, i, :] = v3[:, :, i] = v3[:, :, i]
                        #         elif self.mixed_derivative_handling_mode == MixedDerivativeHandlingModes.Averaged:
                        #             r = np.arange(i, v3.shape[0])
                        #             v3[i, i, r] = v3[i, r, i] = v3[r, i, i] = v3[r, i, i] = np.average(
                        #                 [
                        #                     v3[i, i, r],
                        #                     v3[r, i, i]
                        #                 ],
                        #                 axis=0
                        #             )
                        #         else:
                        #             raise ValueError(
                        #                 "don't know what to do with `mixed_derivative_handling_mode` {} ".format(
                        #                     self.mixed_derivative_handling_mode
                        #                 )
                        #             )

                if intcds is not None and self.backpropagate_internals:
                    # need to internal mode terms and
                    # convert them back to Cartesian mode ones...
                    Qq_derivs = self.get_internal_modes_by_cartesian_modes(len(terms))
                    terms = TensorDerivativeConverter(
                        Qq_derivs,
                        terms
                    ).convert(order=len(terms))
                elif intcds is not None and direct_prop:
                    qQ_derivs = self.get_cartesian_modes_by_internal_modes(len(terms))
                    terms = TensorDerivativeConverter(
                        qQ_derivs,
                        terms
                    ).convert(order=len(terms))

            else:
                terms = TensorDerivativeConverter(x_derivs, u_derivs).convert(order=order)  # , check_arrays=True)


            mu[coord] = (self.derivs[0][coord],) + tuple(terms)


        with self.checkpointer:
            try:
                self.checkpointer['dipole_terms'] = {'x':mu[0], 'y':mu[1], 'z':mu[2]}
            except (OSError, KeyError):
                pass

        return mu

class CoriolisTerm(ExpansionTerms):
    """
    Calculates the Coriolis coupling term
    """
    def get_zetas_and_momi(self):
        # mass-weighted mode matrix
        # (note that we want the transpose not the inverse for unit reasons)
        xQ = self.modes.matrix.T
        # remove the frequency dimensioning? -> this step makes me super uncomfortable but agrees with Gaussian
        freqs = self.freqs
        xQ = xQ / np.sqrt(freqs[:, np.newaxis])
        # reshape xQ so that it looks like atom x mode x Cartesian
        J = np.moveaxis(
            xQ.reshape((len(xQ), self.molecule.num_atoms, 3)),
            1, 0
        )

        # then rotate into the inertial frame
        B_e, eigs = self.inertial_frame
        # print(B_e * UnitsData.convert("Hartrees", "Wavenumbers"))
        # print(eigs)
        # print(J)
        J = np.tensordot(J, eigs, axes=[2, 0])

        # coriolis terms are given by zeta = sum(JeJ^T, n)
        ce = -levi_cevita3
        zeta = sum(
            np.tensordot(
                np.tensordot(ce, J[n], axes=[0, 1]),
                J[n],
                axes=[1, 1])
            for n in range(J.shape[0])
        )

        # raise Exception(np.round(zeta, 3))

        return zeta, B_e

    def get_zetas(self):

        z, m = self.get_zetas_and_momi()

        return z

    def get_terms(self, order=None, J=0):

        if J > 0:
            raise NotImplementedError("currently only have VibRot term for J=0")

        zeta_inert, _ = self.get_zetas_and_momi()
        inert_derivs = self.moment_of_inertia_derivs(order)

        # now we include the frequency dimensioning that comes from the q and p terms in Pi = Zeta*qipj
        freqs = self.freqs
        freq_term = np.sqrt(freqs[np.newaxis, :] / freqs[:, np.newaxis])
        zeta_inert = zeta_inert * freq_term[np.newaxis]

        terms = []
        coriolis = (
                           zeta_inert[np.newaxis, :, :, :, np.newaxis, np.newaxis] # ij
                           * zeta_inert[:, np.newaxis, np.newaxis, np.newaxis, :, :] # kl
        )
        for d in inert_derivs: # expansion of reciprocal inertia tensor
            d = np.expand_dims(d, [2, 3, -1, -2])
            term = d / 2 * coriolis
            terms.append(term)

            # add coordinates for `q`
            coriolis = np.expand_dims(coriolis, -3)


        try:
            self.checkpointer['coriolis_terms'] = terms
        except (OSError, KeyError):
            pass

        return terms

class PotentialLikeTerm(KineticTerms):
    """
    This accounts for the potential-like term.
    In Cartesian diplacement modes this is the Watson U.
    In proper internals, this is the V' term.
    """

    def get_terms(self, order=None, logger=None):

        ics = self.internal_coordinates
        if self.use_cartesian_kinetic_energy or self.backpropagate_internals or ics is None:

            wat_terms = self.moment_of_inertia_derivs(order)
            for i,d in enumerate(wat_terms):
                wat_terms[i] = -np.sum(d[(0, 1, 2), (0, 1, 2), ...], axis=0)

        else:
            ### transform inertia derivs into mode derivs
            YQ_derivs = self.get_cartesians_by_modes(order=2+order)

            I0_derivs = self.inertial_frame_derivatives() # only ever two of these
            if order > 0:
                I0_derivs = I0_derivs + [0]*order
            I0Q_derivs = TensorDerivativeConverter(YQ_derivs, I0_derivs).convert(order=2+order)#check_arrays=True)

            ### pull already computed G-matrix derivs
            # try:
            G_terms = super().get_terms(order=2+order, logger=NullLogger())
            # except:
            #     raise Exception(2+order)

            g_terms = TensorExpansionTerms(G_terms[1:], None, base_qx=G_terms[0], q_name='G')
            detG = g_terms.QX(0).det()

            # oooh this is a dangerous thing to have here
            amu2me = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
            I0 = amu2me * self.molecule.inertia_tensor
            I0_terms = [I0] + I0Q_derivs
            I_terms = TensorExpansionTerms(I0_terms[1:], None, base_qx=I0_terms[0], q_name='I')
            detI = I_terms.QX(0).det()

            # we skip the gamma term from Pickett altogether because it never directly
            # enters, instead only ever being treated as detIdQ - detGdQ
            gamdQ = (detI.dQ()/detI + -1*detG.dQ()/detG).simplify()#check_arrays=True)
            gamdQ.name = "dQ(gam)"
            gamdQQ = gamdQ.dQ().simplify()#check_arrays=True)

            v_term_1 = g_terms.QX(0).dot(gamdQQ, [1, 2], [1, 2])
            v_term_2 = g_terms.QX(1).dot(gamdQ, 3, 1).tr()
            v_term_3 = 1/4 * gamdQ.dot(gamdQ.dot(g_terms.QX(0), 1, 1), 1, 1)

            v0 = (
                    v_term_1
                    + v_term_2
                    + v_term_3
            )

            wat_exprs = [v0]
            for i in range(1, order+1):
                wat_exprs.append(wat_exprs[-1].dQ())#.simplify())

            wat_terms = []
            for i,w in enumerate(wat_exprs):
                try:
                    arr = w.array
                except TensorDerivativeConverter.TensorExpansionError:
                    raise ValueError("failed to construct U({})".format(i))
                else:
                    wat_terms.append(arr)

        try:
            self.checkpointer['psuedopotential_terms'] = wat_terms
        except (OSError, KeyError):
            pass

        return wat_terms