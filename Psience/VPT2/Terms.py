"""
Stores all of the terms used inside the VPT2 representations
"""

import numpy as np, functools as fp, itertools, time, enum, math

from McUtils.Numputils import SparseArray, levi_cevita3
import McUtils.Numputils as nput
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer
from McUtils.Parallelizers import Parallelizer
from McUtils.Zachary import TensorDerivativeConverter, TensorExpansionTerms
from McUtils.Combinatorics import UniquePermutations

from ..Molecools import Molecule, MolecularVibrations, MolecularNormalModes

from .Common import PerturbationTheoryException

__all__ = [
    "ExpansionTerms",
    "KineticTerms",
    "PotentialTerms",
    "CoriolisTerm",
    "PotentialLikeTerm",
    "DipoleTerms",
    "OperatorTerms"
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
    Old = 'old'

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
        "cartesian_by_internal_derivative_method",
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
                 mode_transformation=None,
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
                 mixed_derivative_handling_mode="old",
                 backpropagate_internals=False,
                 direct_propagate_cartesians=False,
                 zero_mass_term=1e7,
                 internal_fd_mesh_spacing=1.0e-2,
                 internal_fd_stencil=None,
                 cartesian_fd_mesh_spacing=1.0e-2,
                 cartesian_fd_stencil=None,
                 cartesian_analytic_deriv_order=0,
                 cartesian_by_internal_derivative_method='old',
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
        self.cartesian_by_internal_derivative_method = cartesian_by_internal_derivative_method
        self.cartesian_fd_mesh_spacing = cartesian_fd_mesh_spacing
        self.cartesian_fd_stencil = cartesian_fd_stencil
        self.cartesian_analytic_deriv_order = cartesian_analytic_deriv_order

        self.internal_by_cartesian_order = internal_by_cartesian_order
        self.cartesian_by_internal_order = cartesian_by_internal_order
        self.jacobian_warning_threshold = jacobian_warning_threshold

        self.internal_coordinates = molecule.internal_coordinates
        self.coords = molecule.coords
        self.masses = molecule._atomic_masses()# * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        self.use_internal_modes = use_internal_modes
        if modes is None:
            modes = molecule.normal_modes.modes
        if hasattr(modes, 'basis') and hasattr(modes.basis, 'to_new_modes'):
            modes = modes.basis
        if hasattr(modes, 'to_new_modes'):
            modes = modes.to_new_modes()
        if mode_transformation is not None:
            modes = modes.apply_transformation(mode_transformation)
        self._modes = modes#.basis
        if undimensionalize is None:
            undimensionalize = not self._check_internal_modes(clean=False)
        if undimensionalize:
            self.raw_modes = modes
            modes = self.modes.make_dimensionless()
            # modes = self.undimensionalize(self.masses, self._modes)
        else:
            self.raw_modes = None
            modes = self._modes
        self._presel_dim = len(self._modes.freqs)
        if mode_selection is not None:
            modes = modes[mode_selection]
        if mode_transformation is not None:
            self._pretf_dim = mode_transformation[0].shape[0]
        else:
            self._pretf_dim = len(modes.freqs)
        self._modes = modes
        self.mode_sel = mode_selection
        self.mode_tf = mode_transformation
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
        else:
            try:
                transf = self.checkpointer['coordinate_transforms']
            except (OSError, KeyError):
                pass
            else:
                self._cached_transforms[molecule] = {JacobianKeys(k):v for k,v in transf.items()}

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
        mat = modes.modes_by_coords
        is_internal = mat.shape[0] == self.coords.shape[0] * self.coords.shape[1] - 6
        self.use_internal_modes = is_internal
        if clean and is_internal:
            self._reshape_internal_modes()
        return is_internal

    def _reshape_internal_modes(self):
        # raise NotImplementedError("ordering has shifted")
        QR = self._modes.modes_by_coords  # derivatives of Q with respect to the internals
        # we need to add zeros for the orientation coordinates
        if not self.strip_embedding and QR.shape[0] != 3 * self.num_atoms:
            _QR = QR
            QR = np.zeros((3 * self.num_atoms, _QR.shape[1]))
            if hasattr(self.internal_coordinates.system, 'embedding_coords'):
                embedding_coords = self.internal_coordinates.system.embedding_coords
                good_coords = np.setdiff1d(np.arange(3 * self.num_atoms), embedding_coords)
                QR[good_coords, :] = _QR
                self._modes.modes_by_coords = QR

        RQ = self._modes.coords_by_modes  # derivatives of internals with respect to Q
        if not self.strip_embedding and RQ.shape[1] != 3 * self.num_atoms:
            _RQ = RQ
            # we need to add zeros for the orientation coordinates
            RQ = np.zeros((_RQ.shape[0], 3 * self.num_atoms))
            if hasattr(self.internal_coordinates.system, 'embedding_coords'):
                embedding_coords = self.internal_coordinates.system.embedding_coords
                good_coords = np.setdiff1d(np.arange(3 * self.num_atoms), embedding_coords)
                RQ[:, good_coords] = _RQ
                self._modes.coords_by_modes = RQ

    @property
    def modes(self):
        # if self._check_internal_modes():
        #     J, = self.get_cart_jacobs([1])
        #     return np.dot(J, self._modes)
        # else:
        #     # cartesian modes
        return self._modes

    # def undimensionalize(self, masses, modes):
    #     """
    #     Removes units from normal modes
    #
    #     :param masses:
    #     :type masses:
    #     :param modes:
    #     :type modes:
    #     :return:
    #     :rtype:
    #     """
    #     L = modes.modes_by_coords.T
    #     Linv = modes.coords_by_modes
    #
    #     freqs = modes.freqs
    #     freq_conv = np.sqrt(np.broadcast_to(freqs[:, np.newaxis], L.shape))
    #     if self._check_internal_modes(clean=False):
    #         conv = freq_conv
    #         inv_conv = 1 / freq_conv
    #     else:
    #         mass_conv = np.sqrt(np.broadcast_to(self._tripmass(masses)[np.newaxis, :], L.shape))
    #         conv = freq_conv * mass_conv
    #         inv_conv = 1 / freq_conv
    #     L = L * conv
    #     Linv = Linv / conv
    #     modes = type(modes)(self.molecule, L.T, inverse=Linv, freqs=freqs)
    #     return modes

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
                for inds in itertools.combinations(all_inds, i):
                    # define a diagonal slice through
                    sel = tuple(slice(None, None, None) if a not in inds else np.arange(s[a]) for a in all_inds)
                    weights[sel] = 1 / math.factorial(i)
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
            stencil = (max(need_jacs) + 2 + (1+max(need_jacs))%2) if self.internal_fd_stencil is None else self.internal_fd_stencil
            # odd behaves better
            with Parallelizer.lookup(self.parallelizer) as par:
                new_jacs = [
                    x.squeeze() if isinstance(x, np.ndarray) else x
                    for x in intcds.jacobian(carts, need_jacs,
                                             # odd behaves better
                                             mesh_spacing=self.internal_fd_mesh_spacing,
                                             stencil=stencil,
                                             # all_numerical=self.all_numerical,
                                             # all_numerical=True,
                                             # analytic_deriv_order=self.cartesian_analytic_deriv_order,
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
            # bad_moms = np.where(mom_i <= 1)
            # mom_i[bad_moms] = 1
            B_e = 1 / (2 * mom_i)  # * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass"))
            # B_e[bad_moms] = 0
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
        YQ = self.modes.coords_by_modes  # derivatives of Q with respect to the Cartesians
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

    def _get_embedding_coords(self):
        try:
            embedding = self.internal_coordinates.system.embedding_coords
        except AttributeError:
            try:
                embedding = self.internal_coordinates.system.converter_options['embedding_coords']
            except KeyError:
                embedding = None
        return embedding
    _cached_transforms = {}
    def get_coordinate_transforms(self,
                                  internal_by_cartesian_order=None,
                                  cartesian_by_internal_order=None,
                                  current_cache=None
                                  ):  #TODO: Cache this as a molecular property instead of some global store here...

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
                msg = [
                        "Getting coordinate transforms for {m}"
                    ]
                opt = dict(m=self.molecule)
                if 'axes_labels' in self.internal_coordinates.system.converter_options:
                    msg.append("Embedding axes: {a}")
                    opt['a'] = self.internal_coordinates.system.converter_options["axes_labels"]

                self.logger.log_print(msg, **opt)

            embedding_coords = self._get_embedding_coords() if self.strip_embedding else None


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
                # def_1 = self.molecule.embedding.cart_fd_defaults
                # try:
                    # self.molecule.embedding.cart_fd_defaults = def_1.copy()
                    # self.molecule.embedding.cart_fd_defaults['analytic_deriv_order'] = self.cartesian_analytic_deriv_order
                int_by_cartesian_jacobs = self.molecule.get_internals_by_cartesians(
                    internal_by_cartesian_order,
                    strip_embedding=self.strip_embedding,

                    mesh_spacing=self.cartesian_fd_mesh_spacing,
                    stencil=self.cartesian_fd_stencil,
                    # all_numerical=True,
                    analytic_deriv_order=self.cartesian_analytic_deriv_order,
                    strip_dummies=self.strip_dummies,
                    parallelizer=self.parallelizer
                )#self.get_cart_jacobs(list(range(1, internal_by_cartesian_order + 1)))
                # finally:
                #     self.molecule.embedding.cart_fd_defaults = def_1
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
                # if embedding_coords is not None:
                #     good_coords = np.setdiff1d(np.arange(3*self.num_atoms), embedding_coords)

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
                        # if embedding_coords is not None:
                        #     x = np.take(x, good_coords, axis=-1)
                        _.append(x)
                    mc = np.expand_dims(mc, 0)
                    cartesian_weighting = cartesian_weighting * mc
                int_by_cartesian_jacobs = _

                current_cache[JacobianKeys.InternalsByCartesians] = int_by_cartesian_jacobs
            else:
                int_by_cartesian_jacobs = current_cache[JacobianKeys.InternalsByCartesians]

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
                cart_by_internal_jacobs = self.molecule.get_cartesians_by_internals(
                    cartesian_by_internal_order,
                    strip_embedding=self.strip_embedding,
                    method=self.cartesian_by_internal_derivative_method,
                    mesh_spacing=self.internal_fd_mesh_spacing,
                    stencil=self.internal_fd_stencil,
                    analytic_deriv_order=self.cartesian_analytic_deriv_order,
                    reembed=self.reembed,
                    planar_ref_tolerance=self.reembed_tol,
                    strip_dummies=self.strip_dummies,
                    parallelizer=self.parallelizer
                )
                #self.get_int_jacobs(list(range(1, cartesian_by_internal_order+1)))

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
                # if embedding_coords is not None:
                #     good_coords = np.setdiff1d(np.arange(3*self.num_atoms), embedding_coords)

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
                        # if embedding_coords is not None:
                        #     for j in range(i+1):
                        #         x = np.take(x, good_coords, axis=j)
                        _.append(x)
                cart_by_internal_jacobs = _

                current_cache[JacobianKeys.CartesiansByInternals] = cart_by_internal_jacobs

            else:
                cart_by_internal_jacobs = current_cache[JacobianKeys.CartesiansByInternals]

            if self._check_internal_modes():
                self.use_internal_modes = True

                QR = self.modes.modes_by_coords  # derivatives of Q with respect to the internals
                # we need to add zeros for the orientation coordinates
                if not self.strip_embedding and QR.shape[0] != 3*self.num_atoms:
                    _QR = QR
                    QR = np.zeros((3*self.num_atoms, _QR.shape[1]))
                    if hasattr(self.internal_coordinates.system, 'embedding_coords'):
                        embedding_coords = self.internal_coordinates.system.embedding_coords
                        good_coords = np.setdiff1d(np.arange(3*self.num_atoms), embedding_coords)
                        QR[good_coords, :] = _QR
                        self.modes.modes_by_coords = QR

                RQ = self.modes.coords_by_modes # derivatives of internals with respect to Q
                if not self.strip_embedding and RQ.shape[1] != 3 * self.num_atoms:
                    _RQ = RQ
                    # we need to add zeros for the orientation coordinates
                    RQ = np.zeros((_RQ.shape[0], 3*self.num_atoms))
                    if hasattr(self.internal_coordinates.system, 'embedding_coords'):
                        embedding_coords = self.internal_coordinates.system.embedding_coords
                        good_coords = np.setdiff1d(np.arange(3*self.num_atoms), embedding_coords)
                        RQ[:, good_coords] = _RQ
                        self.modes.coords_by_modes = RQ

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
                QY = self.modes.modes_by_coords  # derivatives of Q with respect to the Cartesians
                YQ = self.modes.coords_by_modes # derivatives of Cartesians with respect to Q
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
                    # #     self.modes.coords_by_modes.shape
                    # # )
                    # YQ2 = np.concatenate([
                    #     modes,
                    #     self.modes.coords_by_modes  # derivatives of Cartesians with respect to Q
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
                    current_cache[JacobianKeys.CartesianModesByCartesians] = [self.modes.modes_by_coords] + [0]*(len(cart_by_internal_jacobs)-1)
                if (
                        JacobianKeys.CartesiansByCartesianModes not in current_cache
                        or len(current_cache[JacobianKeys.CartesiansByCartesianModes]) < len(cart_by_internal_jacobs)
                ):
                    current_cache[JacobianKeys.CartesiansByCartesianModes] = [self.modes.coords_by_modes] + [0] * (len(cart_by_internal_jacobs) - 1)

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
                    'CartesianModesByCartesians',
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
                    'InternalsByInternalModes',
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
        embedding_coords = self._get_embedding_coords() if (self.strip_embedding or strip_embedding) else None
        if embedding_coords is not None and (strip_embedding and not self.strip_embedding):
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
                 mode_transformation=None,
                 full_surface_mode_selection=None,
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

        self.full_mode_sel = full_surface_mode_selection

        super().__init__(molecule, modes,
                         mode_selection=mode_selection,
                         mode_transformation=mode_transformation,
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
            self._v_derivs = self._canonicalize_derivs(self.freqs, self.masses, self._input_derivs,
                                                       self.full_mode_sel,
                                                       self.mode_tf)

        return self._v_derivs
    @v_derivs.setter
    def v_derivs(self, v):
        self._v_derivs = v


    def _check_mode_terms(self, derivs=None):
        modes_n = len(self.freqs)
        if derivs is None:
            derivs = self.v_derivs
        for d in derivs:
            if d.shape != (modes_n,) * len(d.shape):
                return False
        return True
    def _canonicalize_derivs(self, freqs, masses, derivs, full_mode_sel, mode_transformation):

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
        full_modes_n = self._presel_dim
        modes_n = self._pretf_dim
        tf_n = len(freqs)
        internals_n = 3 * n - 6
        coord_n = 3 * n

        if len(derivs) > 2 and full_mode_sel is not None and thirds.shape[0] != modes_n:
            new_thirds = np.zeros((full_modes_n,) + fcs.shape)
            new_thirds[full_mode_sel,] = thirds
            thirds = new_thirds
        if len(derivs) > 2 and self.mode_sel is not None and thirds.shape[0] == self._presel_dim:
            thirds = thirds[(self.mode_sel,)]

        if len(derivs) > 2 and full_mode_sel is not None and fourths.shape[0] != modes_n:
            new_fourths = np.zeros((full_modes_n, full_modes_n) + fcs.shape)
            new_fourths[np.ix_(full_mode_sel, full_mode_sel)] = fourths
            fourths = new_fourths
        if len(derivs) > 3 and self.mode_sel is not None and fourths.shape[0] == self._presel_dim:
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


        if (
                mode_transformation is not None
                and (
                        len(derivs) > 2
                        and thirds.shape == (modes_n, coord_n, coord_n)
                )
        ):
            if len(derivs) > 3:
                thirds, fourths = nput.tensor_reexpand(
                    [mode_transformation[1]],
                    [thirds, fourths]
                )
            else:
                thirds, = nput.tensor_reexpand(
                    [mode_transformation[1]],
                    [thirds]
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
                xQ2 = self.modes.coords_by_modes
                _, v2x = TensorDerivativeConverter((xQ2, 0), (grad, fcs)).convert(order=2)

                real_freqs = np.diag(v2x)
                nominal_freqs = self.freqs
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
                if thirds.shape == (tf_n, coord_n, coord_n):
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
                elif thirds.shape == (tf_n, tf_n, tf_n):
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
                if fourths.shape == (tf_n, tf_n, coord_n, coord_n):
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
                elif fourths.shape == (tf_n, tf_n, tf_n, tf_n):
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

    @classmethod
    def _symmetrize_mixed_derivatives(cls, derivs, handling_mode, mode_axes, zero_rest=True,
                                      diagonal=True, restricted_diagonal=False, term_id=None, val_axes=0):

        if handling_mode == MixedDerivativeHandlingModes.Old:
            if term_id is None: raise ValueError(';_;')
            if term_id in {'v4_cart', 'v4_int'}:
                v4 = np.zeros_like(derivs)
                for i in range(v4.shape[0]):
                    v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = derivs[i, i, :, :]
                derivs = v4
            elif term_id in {'u3_cart', 'u3_int'}:
                v3 = derivs.copy()
                for i in range(v3.shape[0]):
                    for j in range(i + 1, v3.shape[0]):
                        for k in range(j + 1, v3.shape[0]):
                            for p in itertools.permutations([i, j, k]):
                                v3[p] = 0
                derivs = v3
        else:

            # v3 = terms[2]
            rem_axes = derivs.ndim - mode_axes - val_axes
            nmodes = derivs.shape[0]
            # s = slice(None)

            vals = derivs
            if zero_rest:
                derivs = np.zeros_like(derivs)

            ups = {}
            if restricted_diagonal:
                pos_spec = (
                    (i,) * mode_axes + (j,) * rem_axes
                    for i in range(nmodes)
                    for j in range(i, nmodes)
                )
            elif diagonal:
                pos_spec = (
                    (i,) * mode_axes + p
                    for i in range(nmodes)
                    for p in itertools.combinations_with_replacement(range(i, nmodes), rem_axes)
                )
            else:
                pos_spec = itertools.combinations_with_replacement(range(nmodes), derivs.ndim - val_axes)

            for pos in pos_spec:
                if handling_mode == MixedDerivativeHandlingModes.Unhandled:
                    _, counts = np.unique(pos, return_counts=True)
                    key = tuple(counts)
                    if key in ups:
                        inds = ups[key]
                    else:
                        inds, _ = UniquePermutations(pos).permutations(return_indices=True)
                        ups[key] = inds
                    for idx in inds:
                        new_pos = tuple(pos[p] for p in idx)
                        derivs[new_pos] = vals[new_pos]
                else:
                    if handling_mode == MixedDerivativeHandlingModes.Numerical:
                        val = vals[tuple(reversed(pos))]
                    elif handling_mode == MixedDerivativeHandlingModes.Analytical:
                        val = vals[pos]
                    elif handling_mode == MixedDerivativeHandlingModes.Averaged:
                        val = (vals[tuple(reversed(pos))] + vals[pos]) / 2
                    else:
                        raise ValueError("don't know what to do with `mixed_derivative_handling_mode` {} ".format(handling_mode))

                    _, counts = np.unique(pos, return_counts=True)
                    key = tuple(counts)
                    if key in ups:
                        inds = ups[key]
                    else:
                        inds, _ = UniquePermutations(pos).permutations(return_indices=True)
                        ups[key] = inds
                    for idx in inds:
                        new_pos = tuple(pos[p] for p in idx)
                        derivs[new_pos] = val

        return derivs

    def get_terms(self, order=None, logger=None):

        if self._check_mode_terms():
            return self.v_derivs[1:]

        if logger is None:
            logger = self.logger

        with logger.block(tag='calculating potential derivatives'):
            start = time.time()

            if order is None:
                order = len(self.v_derivs)
            else:
                order += 2

            if self.allow_higher_potential_terms and len(self.v_derivs) < order:
                self.v_derivs = tuple(self.v_derivs) + (0,) * order

            logger.log_print("prepping grad...")
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
                logger.log_print("Cartesian transformation...")
                # this is nice because it eliminates most of the terms in the expansion
                xQ = self.modes.coords_by_modes
                x_derivs = [xQ] + [0] * (order-1)
                # terms = self._get_tensor_derivs(x_derivs, V_derivs, mixed_terms=False, mixed_XQ=self.mixed_derivs)
                if mixed_derivs:
                    terms_base = nput.tensor_reexpand(
                        x_derivs,
                        V_derivs[:2]
                    )
                    terms_mixed = [
                        nput.tensor_reexpand(
                            x_derivs,
                            [0, v],
                            axes=[-1, n+1]
                        )[-1]
                        for n,v in enumerate(V_derivs[2:])
                    ]
                    terms = terms_base + terms_mixed
                    # terms = TensorDerivativeConverter(x_derivs, V_derivs, mixed_terms=[
                    #     [None, v] for v in V_derivs[2:]
                    # ]).convert(order=order)#, check_arrays=True)
                else:
                    terms = nput.tensor_reexpand(x_derivs, V_derivs)

                if mixed_derivs:
                    logger.log_print("handling mixed derivative symmetry ({mode})...", mode=self.mixed_derivative_handling_mode)
                    terms[2] = self._symmetrize_mixed_derivatives(terms[2],
                                                                  self.mixed_derivative_handling_mode,
                                                                  mode_axes=1,
                                                                  term_id='v3_cart'
                                                                  )
                    terms[3] = self._symmetrize_mixed_derivatives(terms[3],
                                                                  self.mixed_derivative_handling_mode,
                                                                  mode_axes=2,
                                                                  restricted_diagonal=True,
                                                                  term_id='v4_cart'
                                                                  )
            elif self._check_internal_modes() and not self._check_mode_terms():
                raise NotImplementedError("...")
                # It should be very rare that we are actually able to make it here
                terms = []
                RQ = self.modes.coords_by_modes
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
                    # raise NotImplementedError("haven't included translation/rotation modes needed to make this work correctly")

                    # raise NotImplementedError("different methods for handling mixed derivatives need patching")
                    QX = self.modes.coords_by_modes
                    x_derivs = [QX] + [0] * (order - 1)
                    cart_terms = TensorDerivativeConverter(x_derivs, V_derivs, mixed_terms=[
                        [None, v] for v in V_derivs[2:]
                    ]).convert(order=order)
                    # provides the derivatives expressed with respect to the Cartesian normal modes
                    v3 = cart_terms[2]
                    if mixed_derivs:
                        cart_terms[2] = self._symmetrize_mixed_derivatives(cart_terms[2],
                                                                           self.mixed_derivative_handling_mode,
                                                                           mode_axes=1,
                                                                           term_id='v3_cart')
                        cart_terms[3] = self._symmetrize_mixed_derivatives(cart_terms[3],
                                                                           self.mixed_derivative_handling_mode,
                                                                           mode_axes=2,
                                                                           restricted_diagonal=True,
                                                                           term_id='v4_cart')
                    qQ_derivs = self.get_cartesian_modes_by_internal_modes(len(cart_terms)-1)
                    terms = TensorDerivativeConverter(
                        qQ_derivs + [0],  # pad for the zeroed out gradient term
                        cart_terms
                    ).convert(order=order)
                    mixed_derivs = False
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

                        v_ders = [v.copy() for v in V_derivs]
                        v_ders[2] = self._symmetrize_mixed_derivatives(V_derivs[2],
                                                                       self.mixed_derivative_handling_mode,
                                                                       mode_axes=1,
                                                                       term_id='v3_int_pre')
                        v_ders[3] = self._symmetrize_mixed_derivatives(V_derivs[3],
                                                                       self.mixed_derivative_handling_mode,
                                                                       mode_axes=2,
                                                                       restricted_diagonal=True,
                                                                       term_id='v4_int_pre')

                        f43 = np.tensordot(qQQ, v_ders[2], axes=[2, 0])
                        fourths = v_ders[3] + f43
                        v_ders = v_ders[:3] + [fourths] + v_ders[4:]

                        terms = TensorDerivativeConverter(x_derivs, v_ders,
                                                          mixed_terms=[
                                                              [None, v] for v in v_ders[2:]
                                                          ]
                                                          ).convert(order=order)  # , check_arrays=True)


                        terms[2] = self._symmetrize_mixed_derivatives(terms[2],
                                                                       self.mixed_derivative_handling_mode,
                                                                       mode_axes=1,
                                                                       term_id='v3_int')
                        terms[3] = self._symmetrize_mixed_derivatives(terms[3],
                                                                       self.mixed_derivative_handling_mode,
                                                                       mode_axes=2,
                                                                       term_id='v4_int')

                    else:
                        terms = TensorDerivativeConverter(x_derivs, V_derivs).convert(order=order)#, check_arrays=True)

                if self.hessian_tolerance is not None:
                        xQ2 = self.modes.coords_by_modes
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
                    xQ2 = self.modes.coords_by_modes
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

            logger.log_print("checking Hessian...")
            if self.hessian_tolerance is not None and np.linalg.norm(np.abs(terms[0] - np.diag(np.diag(terms[0])))) > self.hessian_tolerance:
                raise ValueError(
                    "F-matrix<{}> is not diagonal (got {})".format(
                        terms[0].shape,
                        terms[0] * UnitsData.hartrees_to_wavenumbers
                    )
                )

            new_freqs = np.diag(terms[0])
            old_freqs = self.freqs
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

            end = time.time()
            logger.log_print('took {e:.3f}s...', e=end-start)

            return terms

    @classmethod
    def get_potential_optimized_coordinates(cls, V_expansion, order=2):
        V = [0] + V_expansion
        w = np.diag(V[1])
        forward_derivs = [np.eye(V[1].shape[0])]
        reverse_derivs = [np.eye(V[1].shape[0])]

        w = w[np.newaxis, :]
        for o in range(1, order + 1):
            V_rem = TensorDerivativeConverter.convert_fast(forward_derivs, V,
                                                           order=o+2, val_axis=0
                                                           )[-1]
            w = w[np.newaxis]
            new_Q = -V_rem / ((o+2)*w)
            forward_derivs = forward_derivs + [new_Q]
            new_R = -nput.tensordot_deriv(forward_derivs, reverse_derivs + [0], o)[-1]
            reverse_derivs = reverse_derivs + [new_R]

        return forward_derivs, reverse_derivs

    def optimize_coordinates(self, order=2):
        V = list(reversed([self[o] for o in range(order, -1, -1)]))
        RX = self.get_cartesians_by_internals(order=order, strip_embedding=True)
        XR = self.get_internals_by_cartesians(order=order, strip_embedding=True)
        QR, RQ = self.get_potential_optimized_coordinates(V, order=order)

        QX = TensorDerivativeConverter.convert_fast(QR, RX)
        XQ = TensorDerivativeConverter.convert_fast(XR, RQ)

        return (QR, RQ), (QX, XQ)

class KineticTerms(ExpansionTerms):
    """Represents the KE coefficients"""

    __props__ = ExpansionTerms.__props__ + (
        'g_derivative_threshold',
        "gmatrix_tolerance",
        'use_cartesian_kinetic_energy',
        "check_input_gmatrix"
        "freq_tolerance"
    )
    def __init__(self,
                 molecule,
                 g_derivative_threshold=1e-3,
                 gmatrix_tolerance=1e-6,
                 use_cartesian_kinetic_energy=False,
                 check_input_gmatrix=True,
                 freq_tolerance=2e-3,
                 **opts
                 ):
        super().__init__(molecule, **opts)
        self.g_derivative_threshold = g_derivative_threshold
        self.gmatrix_tolerance = gmatrix_tolerance
        self.use_cartesian_kinetic_energy = use_cartesian_kinetic_energy
        self.freq_tolerance = freq_tolerance
        self.check_input_gmatrix = check_input_gmatrix

    def get_terms(self, order=None, logger=None, return_expressions=False):

        if logger is None:
            logger = self.logger
        with logger.block(tag='calculating G-matrix derivatives'):
            start = time.time()

            dot = DumbTensor._dot
            shift = DumbTensor._shift
            intcds = self.internal_coordinates
            if self.use_cartesian_kinetic_energy or intcds is None or self.backpropagate_internals:
                # this is nice because it eliminates a lot of terms in the expansion
                J = self.modes.modes_by_coords
                G = dot(J, J, axes=[[0, 0]])
                if order == 0:
                    terms = [G]
                else:
                    terms = [G] + [0]*(order)
                G_terms = None
            else:
                # should work this into the new layout
                uses_internal_modes = self._check_internal_modes()
                if uses_internal_modes:
                    QY_derivs = self.get_internals_by_cartesians(order=order + 1) # really dRdY derivatives
                    YQ_derivs = self.get_cartesians_by_internals(order=order + 1) # really dYdR derivatives
                else:
                    QY_derivs = self.get_modes_by_cartesians(order=order+1)
                    YQ_derivs = self.get_cartesians_by_modes(order=order+1)
                    # QY_derivs = self.get_internals_by_cartesians(order=order + 1) # really dRdY derivatives
                    # YQ_derivs = self.get_cartesians_by_internals(order=order + 1) # really dYdR derivatives
                # RQ = dot(YQ, RY)

                term_getter = TensorDerivativeConverter(YQ_derivs, QY_derivs).terms
                # term_getter = TensorDerivativeConverter([np.eye(9,9), 0, 0], QY_derivs).terms
                term_getter.v_name = 'Y'
                J = term_getter.XV(1)
                G_terms = [J.dot(J, 1, 1)]
                for i in range(1, order+1):
                    g_cur = G_terms[-1].dQ()#.simplify(check_arrays=True)
                    G_terms.append(g_cur)
                terms = [x.array for x in G_terms]
                #
                # """
                # [[-0.0e+00  3.9e-07  8.8e-07]
                #  [ 3.9e-07 -0.0e+00  8.8e-07]
                #  [ 8.8e-07  8.8e-07  0.0e+00]]
                #  """
                #
                # """
                # [[-3.800e-07 -3.142e-05  0.000e+00]
                #  [-3.142e-05  1.843e-05 -0.000e+00]
                #  [ 0.000e+00 -0.000e+00 -1.198e-05]]
                # """
                #
                # QR_derivs = self.get_internals_by_internal_modes(order=order+1)
                # RQ = self.get_internal_modes_by_internals(order=1)[0]
                # # print(
                # #     [q.shape for q in QY_derivs],
                # #     [g.shape for g in terms])
                # terms = [terms[0]] + TensorDerivativeConverter.convert_fast(YQ_derivs, terms[1:], val_axis=0)
                # terms = [
                #     np.tensordot(np.tensordot(t, RQ, axes=[-2, 0]), RQ, axes=[-2, 0])
                #     for t in terms
                # ]
                # print(terms[0])

                if uses_internal_modes:
                    QR = self.modes.modes_by_coords
                    RQ = self.modes.coords_by_modes
                    for i,g in enumerate(terms):
                        for j in range(2):
                            g = np.tensordot(QR, g, axes=[0, -1])
                        for j in range(i):
                            g = np.tensordot(RQ, g, axes=[1, -1])
                        terms[i] = g

                for i,t in enumerate(terms):
                    if i == 0:
                        if self.gmatrix_tolerance is not None and np.linalg.norm(np.abs(terms[0] - np.diag(np.diag(terms[0])))) > self.gmatrix_tolerance:
                            raise ValueError("G-matrix is not diagonal (got {})".format(terms[0] * UnitsData.hartrees_to_wavenumbers))
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

            if return_expressions:
                G_terms = terms, G_terms
            else:
                G_terms = terms

            if self.freq_tolerance is not None and self.check_input_gmatrix:
                real_freqs = np.diag(G_terms[0])
                nominal_freqs = self.freqs
                # deviation on the order of a wavenumber can happen in low-freq stuff from numerical shiz
                if self.freq_tolerance is not None:
                    if np.max(np.abs(nominal_freqs - real_freqs)) > self.freq_tolerance:
                        raise PerturbationTheoryException(
                            "Input frequencies aren't obtained when transforming the G matrix;"
                            " this likely indicates issues with the input mode vectors"
                            " got \n{}\n but expected \n{}\n".format(
                                real_freqs * UnitsData.convert("Hartrees", "Wavenumbers"),
                                nominal_freqs * UnitsData.convert("Hartrees", "Wavenumbers")
                            )
                        )

            try:
                self.checkpointer['gmatrix_terms'] = G_terms
            except (OSError, KeyError):
                pass

            end = time.time()
            logger.log_print("took {e:.3f}s...", e=end-start)

            return G_terms

    @classmethod
    def _dRGQ_partition_contrib(cls, partition, R, G):
        r1, r2, s = partition
        if s - 1 >= len(G): return 0
        if r1 - 1 >= len(R) or r2 - 1 >= len(R): return 0

        base_term = G[s]
        if isinstance(base_term, (int, float, np.integer, np.floating)) and base_term == 0:
            return 0

        r_perm_counter = 1
        r_perm_idx = []
        # g_perm_idx = []

        a = base_term.ndim - 2
        for r in [r1, r2]:
            d = R[r]
            if isinstance(d, (int, float, np.integer, np.floating)) and d == 0:
                return 0
            if r == 0:
                base_term = np.moveaxis(base_term, a, -1)
            else:
                base_term = np.tensordot(base_term, d, axes=[a, 0])
                r_perm_idx.extend([r_perm_counter] * r)
                r_perm_counter -= 1
        if r2 > 0:
            # 1 1 -> (x, i, y, j) -> (-3 -> -2)
            base_term = np.moveaxis(base_term, -(r2 + 2), -2) # axis from the r1 contraction
        # if s > 0:
        g_perm_idx = ([1]*s) + ([0]*(r1+r2))

        if r1 > 0 or r2 > 0:
            # r1 and r2 need to be permuted to preserve symmetry, but all perms can overcount
            nterms = TensorDerivativeConverter.compute_partition_terms([p for p in partition[:2] if p > 0])
            r_perms, _ = UniquePermutations(r_perm_idx).permutations(return_indices=True)

            overcount = len(r_perms) / nterms
            base_term = base_term / overcount

            # if r1 != r2:
            base_term = base_term + np.moveaxis(base_term, -1, -2)


            # g indices need to be swapped into the r coords sometimes too
            g_perms, _ = UniquePermutations(g_perm_idx).permutations(return_indices=True)

            # take direct product of permutation indices
            perm_inds = []
            # print(overcount)
            padding = list(range(base_term.ndim - 2, base_term.ndim))
            for p in g_perms:
                for r in r_perms:
                    perm_inds.append(
                        list(p[:s]) + list(p[s:][r]) + padding
                    )

            base_term = sum(
                base_term.transpose(p)
                for p in perm_inds
            )

        return base_term

    @classmethod
    def _dRGQ_derivs(cls, R, G, o):
        # parts = IntegerPartitionPermutations(o, dim=3).get_partition_permutations(flatten=True)

        # print("="*50)
        total_cont = 0
        for g in range(o+1):
            rem = (o+1)-g
            # don't want to go smaller than half of rem, b.c. symmetry
            for r1 in range(rem//2, rem):
                r2 = o - (g + r1)
                # print("-"*10)
                # print(r1, r2, g)
                pc = cls._dRGQ_partition_contrib([r1, r2, g], R, G)
                # print(pc[0])
                total_cont += pc
        # print("_"*50)
        return total_cont

    @classmethod
    def reexpress_G(self,
                    G_expansion, forward_derivs, reverse_derivs=None, order=2
                    # G_subexpansion=None
                    ):
        """
        Apply a coordinate transformation to the G-matrix

        :param forward_derivs:
        :param reverse_derivs:
        :param order:
        :return:
        """

        if reverse_derivs is None:
            reverse_derivs = nput.inverse_transformation(forward_derivs, order)
        if forward_derivs is None:
            forward_derivs = nput.inverse_transformation(reverse_derivs, order)


        from ..Molecools.Hamiltonian import GMatrixExpansion

        return GMatrixExpansion.reexpress_G(
                    G_expansion, forward_derivs, reverse_derivs,
                    order
                    )


        R = reverse_derivs
        Q = forward_derivs
        G = G_expansion

        # if G_subexpansion is None:
        G_R = [self._dRGQ_derivs(R, G, o) for o in range(1, order+1)]
        # else:
        #     G_R = G_subexpansion + [self._dRGQ_derivs(R, G, o) for o in range(1, order + 1)]

        # print(G_R[1][0, 0])
        # raise Exception(...)

        return [G[0]] + TensorDerivativeConverter.convert_fast(Q, G_R, order=order, val_axis=0)
    def reexpress(self, forward_derivs, reverse_derivs=None, order=2):
        """
        Finds a coordinate transformation the give 0 contribution to the G-matrix

        :param forward_derivs:
        :param reverse_derivs:
        :param order:
        :return:
        """

        G = list(reversed([self[o] for o in range(order, -1, -1)]))
        return self.reexpress_G(G, forward_derivs, reverse_derivs)

    @classmethod
    def get_kinetic_optimized_coordinates(cls, G_expansion, order=2):
        # we do this order-by-order by noting that at each order we end up with a term
        # that includes the highest derivative of the new coordinates with respect to the old
        # multiplied by the G matrix, done two ways, so a transformation that eliminates the
        # coordinates is just (R^n)_a...bij = -1/2w_i [rem]_a...bij

        w = np.diag(G_expansion[0])
        ndim = len(w)
        R = np.eye(ndim)
        R2 = np.zeros((ndim,)*3)
        G1 = G_expansion[1]
        print(G1)
        for p in itertools.combinations_with_replacement(range(ndim), 3):
            for i,j,k in itertools.permutations(p):
                # print(i, j, k)
                R2[i,j,k] = -(w[i]*G1[i, j, k] + w[j]*G1[j, i, k] - w[k]*G1[k, j, i]) / (2 * w[i] * w[j])

        # raise Exception(R2)
        # raise Exception(R2 - np.transpose(R2, (1, 0, 2)))

        new_G = cls.reexpress_G(G_expansion, None, [R, R2], order=2)
        raise Exception(new_G[1])

        G = G_expansion
        w = np.diag(G[0])
        forward_derivs = [np.eye(G[0].shape[0])]
        reverse_derivs = [np.eye(G[0].shape[0])]

        w = 2 * w[:, np.newaxis]
        for o in range(1, order+1):
            G_rem = cls.reexpress_G(G, forward_derivs, reverse_derivs + [0], order=o)[-1]
            w = w[np.newaxis]
            new_R = -G_rem / w
            reverse_derivs = reverse_derivs + [new_R]
            new_Q = -nput.tensordot_deriv(forward_derivs + [0], reverse_derivs, o)[-1]
            forward_derivs = forward_derivs + [new_Q]

        return forward_derivs, reverse_derivs
    def optimize_coordinates(self, order=2):
        G = list(reversed([self[o] for o in range(order, -1, -1)]))
        RX = self.get_cartesians_by_internals(order=order, strip_embedding=True)
        XR = self.get_internals_by_cartesians(order=order, strip_embedding=True)
        QR, RQ = self.get_kinetic_optimized_coordinates(G, order=order)

        QX = TensorDerivativeConverter.convert_fast(QR, RX)
        XQ = TensorDerivativeConverter.convert_fast(XR, RQ)

        return (QR, RQ), (QX, XQ)


class CoriolisTerm(ExpansionTerms):
    """
    Calculates the Coriolis coupling term
    """
    def get_zetas_and_momi(self):
        # mass-weighted mode matrix
        # (note that we want the transpose not the inverse for unit reasons)
        xQ = self.modes.modes_by_coords.T
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
            # detG = g_terms.QX(0).det()

            # oooh this is a dangerous thing to have here
            # amu2me = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
            I0 = self.molecule.inertia_tensor
            I0_terms = [I0] + I0Q_derivs
            I_terms = TensorExpansionTerms(I0_terms[1:], None, base_qx=I0_terms[0], q_name='I')
            # detI = I_terms.QX(0).det()

            # we skip the gamma term from Pickett altogether because it never directly
            # enters, instead only ever being treated as ln(detI) - ln(detG)
            # and then we can make use of an identity from the Matrix Cookbook to say
            # ln(detX).dQ() = (X.inverse().dot(X.dQ)).tr()
            G0 = g_terms.QX(0)
            I0 = I_terms.QX(0)
            lndetdQ = lambda X:(X.inverse().dot(X.dQ(), 2, 2).shift(2, 1)).tr(2, 3)
            gamdQ = lndetdQ(I0)+-lndetdQ(G0)
            gamdQ.name = "dQ(gam)"

            # detG = G0.det()
            # detI = I0.det()
            # gamdQ_alt = (detI.dQ() / detI + -1 * detG.dQ() / detG).simplify()  # check_arrays=True)
            # raise Exception(gamdQ.array, gamdQ_alt.array)

            gamdQQ = gamdQ.dQ().simplify()

            # terms = [x.asarray() for x in gamdQQ.dQ().simplify().terms]
            # for t in gamdQQ.simplify().terms:
            #     print(t)
            # print("-=---==--==---")
            # for t in gamdQQ.dQ().simplify().terms:
            #     print(t)
            # raise Exception(gamdQQ.dQ().array)

            # print(terms)
            # print(gamdQQ.terms[0])
            # print(gamdQQ.dQ().simplify().terms[1])
            # sums = [terms[0]]
            # terms[0].transpose(1, 2, 0) + terms[2] +
            # b = terms[2]
            # for i in range(len(b)):
            #     for j in range(i+1, len(b)):
            #         for k in range(j + 1, len(b)):
            #             terms = []
            #             for p in itertools.permutations([i, j, k]):
            #                 terms.append(b[i, j, k]-b[p])
            #             if max(abs(x) for x in terms) > 1e-14:
            #                 raise ValueError("ugh")
            # for t in gamdQQ.dQ().terms:
            #     print(t)
            #     print(t.asarray(print_terms=True))
            # print(terms[1]+terms[3])
            # for s in terms:
            #     print(s)
            # for x in terms[1:]:
            #     sums.append(sums[-1]+x)
            # for s in sums:
            #     print(s)


            # raise Exception(...)
            #

            # a = gamdQQ.dQ().dQ().array
            # b = gamdQQ.dQ().array
            # for i in range(len(b)):
            #     for j in range(i+1, len(b)):
            #         for k in range(j + 1, len(b)):
            #             terms = []
            #             for p in itertools.permutations([i, j, k]):
            #                 terms.append(b[i, j, k]-b[p])
            #             if max(abs(x) for x in terms) > 1e-14:
            #                 raise ValueError("ugh")
            # ugh = gamdQQ.dQ().dQ().simplify()
            # raise Exception(ugh.array)
            # print(sum(u.array for u in ugh.terms))
            # for p in itertools.permutations([0, 1, 2, 3]):
            #     b = ugh.terms[0].array + ugh.terms[1].array.transpose(p)
            #     if abs(b[0, 0, 0, 1] - b[1, 0, 0, 0]) < 1e-6:
            #         print(b)
            # raise Exception(type(ugh.terms[0]))

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
                    if arr.shape == ():
                        arr = float(arr)
                    wat_terms.append(arr)
            # print(wat_terms)

        try:
            self.checkpointer['psuedopotential_terms'] = wat_terms
        except (OSError, KeyError):
            pass

        return wat_terms

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
                 mode_transformation=None,
                 full_surface_mode_selection=None,
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
        self.full_mode_sel = full_surface_mode_selection
        super().__init__(molecule, modes=modes, mode_selection=mode_selection, mode_transformation=mode_transformation,
                         logger=logger, parallelizer=parallelizer, checkpointer=checkpointer,
                         **opts
                         )
        self.mixed_derivs = mixed_derivs
        if self.mixed_derivs is None:
            self.mixed_derivs = mixed_derivs
        if dipole_derivatives is None:
            dipole_derivatives = molecule.dipole_derivatives
        self.derivs = self._canonicalize_derivs(self.freqs, self.masses, dipole_derivatives,
                                                self.full_mode_sel, self.mode_tf)

    def _canonicalize_derivs(self, freqs, masses, derivs, full_mode_sel, mode_transformation):
        """
        Makes sure all of the dipole moments are clean and ready to rotate
        """
        if (
                len(derivs) == 3
                and not (
                    all(
                        isinstance(d, (int, float, np.integer, np.floating))
                        or (isinstance(d, np.ndarray) and d.ndim > 0 and d.shape[-1] == 3)
                        for d in derivs
                    )
                )

        ):
            # need to effectively transpose things...
            tp_derivs = [[] for _ in range(len(derivs[0]))]
            for ders in derivs:
                for i,d in enumerate(ders):
                    tp_derivs[i].append(d)
            derivs = [
                np.moveaxis(np.array(x), 0, -1)
                for x in tp_derivs
            ]

        _ = list(derivs[:1])
        for n, w in enumerate(derivs[1:]):
            w = np.asanyarray(w)
            if w.shape == (3,):
                if np.allclose(w, [0, 0, 0]):
                    xn = 3*len(masses)
                    w = np.zeros((xn,) * (n + 1) + (3,), dtype=float)
                else:
                    raise ValueError("shape mismatch for dipole derivs")
            _.append(w)
        derivs = _

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
        full_mode_n = self._presel_dim
        modes_n = self._pretf_dim
        tf_n = len(freqs)
        internals_n = 3 * n - 6
        coord_n = 3 * n


        if len(derivs) > 2 and full_mode_sel is not None and seconds.shape[0] != modes_n:
            new_seconds = np.zeros((full_mode_n,) + grad.shape)
            new_seconds[full_mode_sel,] = seconds
            seconds = new_seconds
        if len(derivs) > 2 and self.mode_sel is not None and seconds.shape[0] == self._presel_dim:
            seconds = seconds[(self.mode_sel,)]

        if len(derivs) > 2 and full_mode_sel is not None and thirds.shape[0] != modes_n:
            new_thirds = np.zeros((full_mode_n,full_mode_n) + grad.shape)
            new_thirds[np.ix_(full_mode_sel, full_mode_sel)] = thirds
            thirds = new_thirds
        if self.mode_sel is not None and thirds.shape[0] == self._presel_dim:
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
                    "{0}.{1}: dimension of dipole derivative array ({2[0]}) is not {3} or {4}".format(
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

        if (
                mode_transformation is not None
                and (
                        len(derivs) > 2
                        and seconds.shape == (modes_n, coord_n, 3)
                )
        ):
            if len(derivs) > 3:
                seconds, thirds = nput.tensor_reexpand(
                    [mode_transformation[1]],
                    [seconds, thirds]
                )
            else:
                seconds, = nput.tensor_reexpand(
                    [mode_transformation[1]],
                    [seconds]
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
            if seconds.shape == (tf_n, coord_n, 3):
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
                )
            else:
                undimension_2 = 1
            seconds = seconds / undimension_2
            all_derivs.append(seconds)

        if len(derivs) > 3:
            if thirds.shape == (tf_n, tf_n, coord_n, 3):
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
                undimension_3 =(
                        m_conv[:, np.newaxis, np.newaxis, np.newaxis]
                        * m_conv[np.newaxis, :, np.newaxis, np.newaxis]
                        * m_conv[np.newaxis, np.newaxis, :, np.newaxis]
                )
            elif thirds.shape == (tf_n, tf_n, tf_n, 3):
                if self.mixed_derivs is None:
                    self.mixed_derivs = False
                undimension_3 = (
                        f_conv[:, np.newaxis, np.newaxis, np.newaxis]
                        * f_conv[np.newaxis, :, np.newaxis, np.newaxis]
                        * f_conv[np.newaxis, np.newaxis, :, np.newaxis]
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
                undimension = np.expand_dims(undimension, -1)
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
        modes_n = self._pretf_dim
        if derivs is None:
            derivs = self.derivs[1:]
        for d in derivs:
            if d.shape != (modes_n,) * len(d.shape):
                return False
        return True
    def get_terms(self, order=None, logger=None):

        if logger is None: logger = self.logger
        with logger.block(tag='calculating dipole derivatives'):
            start = time.time()

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
                xQ = self.modes.coords_by_modes
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
                if mixed_derivs:
                    mixed_terms = [
                        [u_derivs[1]],  # dVdQXX
                        [u_derivs[2]]  # dVdQQXX
                    ]
                    if intcds is not None:
                        if (  # handle mixed derivative resymmetrization
                                self.mixed_derivative_handling_mode != MixedDerivativeHandlingModes.Unhandled
                                and order > 1
                        ):
                            # d^2X/dQ^2@dU/dX + dX/dQ@dU/dQdX
                            xQ = self.modes.coords_by_modes
                            v1, v2, v3 = TensorDerivativeConverter(
                                [xQ] + [0] * (order-1),
                                u_derivs,
                                mixed_terms=mixed_terms,
                                values_name="U"
                            ).convert(order=3)  # , check_arrays=True)

                            # Gaussian gives slightly different constants
                            # depending on whether the analytic or numerical derivs
                            # were transformed
                            v2 = PotentialTerms._symmetrize_mixed_derivatives(v2,
                                                                              self.mixed_derivative_handling_mode, 1,
                                                                              term_id='u2_cart')
                            v3 = PotentialTerms._symmetrize_mixed_derivatives(v3,
                                                                              self.mixed_derivative_handling_mode, 2,
                                                                              term_id='u3_cart')

                            Qx = self.modes.modes_by_coords
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
                    terms = list(terms)
                    if (
                            self.mixed_derivative_handling_mode != MixedDerivativeHandlingModes.Unhandled
                            and order > 1
                    ):
                        terms[1] = PotentialTerms._symmetrize_mixed_derivatives(terms[1],
                                                                                self.mixed_derivative_handling_mode, 1,
                                                                                term_id='u2_int')
                        terms[2] = PotentialTerms._symmetrize_mixed_derivatives(terms[2],
                                                                                self.mixed_derivative_handling_mode, 2,
                                                                                term_id='u3_int')

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

            end = time.time()
            logger.log_print("took {e:.3f}s", e=end-start)

        return mu

class OperatorTerms(ExpansionTerms):
    """
    Literally as simple as it comes for an operator expansion.
    One dimensional, no mixed derivative stuff.
    """
    __props__ = ExpansionTerms.__props__ + (
        "operator_derivatives",
    )
    def __init__(self,
                 molecule,
                 operator_derivatives=None,
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
        if operator_derivatives is None:
            raise ValueError("can't transform derivatives without derivatives")
        super().__init__(molecule, modes=modes, mode_selection=mode_selection,
                         logger=logger, parallelizer=parallelizer, checkpointer=checkpointer,
                         **opts
                         )
        self.derivs = self._canonicalize_derivs(self.freqs, self.masses, operator_derivatives)

    def _check_mode_terms(self, derivs=None):
        modes_n = self._pretf_dim
        if derivs is None:
            derivs = self.derivs
        for d in derivs:
            if d.shape != (modes_n,) * len(d.shape):
                return False
        return True
    def _canonicalize_derivs(self, freqs, masses, derivs):

        if self._check_mode_terms(derivs):
            return derivs

        n = self.num_atoms
        modes_n = self._pretf_dim
        internals_n = 3 * n - 6
        coord_n = 3 * n

        for i in range(len(derivs)):
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
            all_derivs = []

            for i in range(len(derivs)):
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
            return self.derivs[1:]

        if logger is None:
            logger = self.logger

        logger.log_print('transformaing operator derivatives')

        if order is None:
            order = len(self.derivs)

        if len(self.derivs) < order:
            self.derivs = tuple(self.derivs) + (0,) * order

        # Use the Molecule's coordinates which know about their embedding by default
        intcds = self.internal_coordinates
        direct_prop = (
                self.direct_propagate_cartesians
                and not (isinstance(self.direct_propagate_cartesians, str) and self.direct_propagate_cartesians == 'dipoles')
        )
        if intcds is None or direct_prop:
            # this is nice because it eliminates most of the terms in the expansion
            xQ = self.modes.coords_by_modes
            x_derivs = [xQ] + [0] * (order-1)
            terms = TensorDerivativeConverter(x_derivs, self.derivs).convert(order=order)#, check_arrays=True)
        elif self._check_internal_modes() and not self._check_mode_terms():
            raise NotImplementedError("...")
            # It should be very rare that we are actually able to make it here
            terms = []
            RQ = self.modes.coords_by_modes
            for v in V_derivs:
                for j in range(v.ndim):
                    v = np.tensordot(RQ, v, axes=[1, -1])
                terms.append(v)
        else:
            x_derivs = self.get_cartesians_by_modes(order=order)
            # raise Exception(x_derivs[1])
            x_derivs = list(x_derivs) #+ [0] # gradient term never matters
            terms = TensorDerivativeConverter(x_derivs, self.derivs).convert(order=order)#, check_arrays=True)

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

        return terms