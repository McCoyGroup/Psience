"""
Stores all of the terms used inside the VPT2 representations
"""

import numpy as np, functools as fp, itertools as ip, time

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

class ExpansionTerms:
    """
    Base class for kinetic, potential, and dipole derivative terms
    """
    _cached_jacobians = {}
    def __init__(self,
                 molecule,
                 modes=None,
                 mode_selection=None,
                 undimensionalize=True,
                 logger=None,
                 parallelizer=None,
                 checkpointer=None,
                 numerical_jacobians=True,
                 eckart_embed=True
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
        dummies = self.molecule.dummy_positions
        dummy_comp = np.setdiff1d(np.arange(molecule.num_atoms), dummies)
        self.internal_coordinates = molecule.internal_coordinates
        self.coords = molecule.coords
        self.masses = molecule.masses[dummy_comp] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        if modes is None:
            modes = molecule.normal_modes.modes
        if undimensionalize:
            self.raw_modes = modes
            modes = self.undimensionalize(self.masses, modes.basis)
        else:
            self.raw_modes = None
        if mode_selection is not None:
            modes = modes[mode_selection]
        self.modes = modes
        self.mode_sel = mode_selection
        self.freqs = self.modes.freqs
        self._inert_frame = None

        self.reembed=eckart_embed
        self.all_numerical=numerical_jacobians

        if logger is None:
            logger = NullLogger()
        self.logger = logger
        if parallelizer is None:
            parallelizer = Parallelizer.get_default()
        self.parallelizer = parallelizer
        if checkpointer is None:
            checkpointer = NullCheckpointer()
        self.checkpointer = checkpointer

    def undimensionalize(self, masses, modes):
        L = modes.matrix.T
        freqs = modes.freqs
        freq_conv = np.sqrt(np.broadcast_to(freqs[:, np.newaxis], L.shape))
        mass_conv = np.sqrt(np.broadcast_to(self._tripmass(masses)[np.newaxis, :], L.shape))
        L = L * freq_conv * mass_conv
        Linv = (L / freq_conv**2)
        modes = type(modes)(self.molecule, L.T, inverse=Linv, freqs=freqs)
        return modes

    @staticmethod
    def _tripmass(masses):
        return np.broadcast_to(masses[np.newaxis, :], (3, len(masses))).T.flatten()

    def get_terms(self, order=None):
        raise NotImplemented

    def get_term(self, t):
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

    internal_fd_mesh_spacing = 1.0e-2
    internal_fd_stencil = 9
    def get_int_jacobs(self, jacs):
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
        if max_jac > len(exist_jacs):
            need_jacs = [x+1 for x in range(0, max_jac)]
            with Parallelizer.lookup(self.parallelizer) as par:
                new_jacs = [x.squeeze() for x in intcds.jacobian(carts, need_jacs,
                                                                 mesh_spacing=self.internal_fd_mesh_spacing,
                                                                 stencil=self.internal_fd_stencil,
                                                                 all_numerical=self.all_numerical,
                                                                 converter_options=dict(
                                                                     reembed=self.reembed,
                                                                     strip_dummies=True
                                                                 ),
                                                                 parallelizer=par
                                                                 )]
            self._cached_jacobians[self.molecule]['int'] = new_jacs
            exist_jacs = new_jacs
        return [exist_jacs[j-1] for j in jacs]

    cartesian_fd_mesh_spacing = 1.0e-5
    cartesian_fd_stencil = 5
    cartesian_analytic_deriv_order = 1
    def get_cart_jacobs(self, jacs):
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
        # print("C", jacs, max_jac, len(exist_jacs))
        if max_jac > len(exist_jacs):
            with Parallelizer.lookup(self.parallelizer) as par:
                need_jacs = [x + 1 for x in range(0, max_jac)]
                new_jacs = [
                    x.squeeze() for x in ccoords.jacobian(internals, need_jacs,
                                                          mesh_spacing=self.cartesian_fd_mesh_spacing,
                                                          stencil=self.cartesian_fd_stencil,
                                                          # all_numerical=True,
                                                          analytic_deriv_order=self.cartesian_analytic_deriv_order,
                                                          converter_options=dict(strip_dummies=True),
                                                          parallelizer=par
                                                          )
                ]

                self._cached_jacobians[self.molecule]['cart'] = new_jacs
                exist_jacs = new_jacs
        return [exist_jacs[j-1] for j in jacs]

    @property
    def inertial_frame(self):

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

        mass = np.sqrt(self.masses)
        carts = mass[:, np.newaxis] * self.molecule.coords  # mass-weighted Cartesian coordinates

        ### compute basic inertia tensor derivatives
        # first derivs are computed as a full (nAt, 3, I_rows (3), I_cols (3)) tensor
        # and then reshaped to (nAt * 3, I_rows, I_cols)
        eyeXeye = np.eye(9).reshape(3, 3, 3, 3).transpose((2, 0, 1, 3))
        I0Y_1 = np.tensordot(carts, eyeXeye, axes=[1, 0])

        nAt = carts.shape[0]
        nY = nAt * 3
        I0Y_21 = np.reshape(np.eye(3), (9,))[np.newaxis, :, np.newaxis] * carts[:, np.newaxis,
                                                                          :]  # a flavor of outer product
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

    @classmethod
    def _get_tensor_derivs(cls, x_derivs, V_derivs, order=4, mixed_XQ=False, mixed_terms=False):
        """
        Returns the derivative tensors of the potential with respect to the normal modes
        (note that this is fully general and the "cartesians" and "normal modes" can be any coordinate sets)
        :param x_derivs: The derivatives of the cartesians with respect to the normal modes
        :type x_derivs:
        :param V_derivs: The derivative of the potential with respect to the cartesians
        :type V_derivs:
        :param mixed_XQ: Whether the v_derivs[2] = V_Qxx and v_derivs[3] = V_QQxx or not
        :type mixed_XQ: bool
        """

        dot = DumbTensor._dot
        shift = DumbTensor._shift

        derivs = [None] * order

        # First Derivs
        xQ = x_derivs[0]
        Vx = V_derivs[0]
        V_Q = dot(xQ, Vx)

        derivs[0] = V_Q
        if order == 1:
            return tuple(derivs)

        # Second Derivs
        xQQ = x_derivs[1]
        Vxx = V_derivs[1]

        V_QQ_1 = dot(xQQ, Vx)
        V_QQ_2 = dot(xQ, dot(xQ, Vxx, axes=[[1, 0]]), axes=[[1, 1]])
        V_QQ_terms = (V_QQ_1, V_QQ_2)
        V_QQ = sum(x for x in V_QQ_terms if not isinstance(x, int))
        derivs[1] = V_QQ
        if order == 2:
            return tuple(derivs)

        # Third Derivs
        xQQQ = x_derivs[2]
        Vxxx = V_derivs[2]

        # If Q is just an expansion in X all of these terms will disappear except for V_QQQ_5

        # Gradient contribution
        V_QQQ_1 = dot(xQQQ, Vx)
        # Second deriv.
        # we generate the base arrangement
        Q32 = dot(xQQ, dot(xQ, Vxx, axes=[[1, 0]]), axes=[[2, 1]])
        # then we do the transpositions that put the xQ coordinate inside the xQQ
        if not isinstance(Q32, int):
            X = tuple(range(3, Q32.ndim))
            V_QQQ_2_terms = [
                Q32.transpose(0, 1, 2, *X),
                Q32.transpose(1, 2, 0, *X),
                Q32.transpose(2, 0, 1, *X)
                ]
            V_QQQ_2 = sum(V_QQQ_2_terms)
        else:
            V_QQQ_2 = 0

        # Third derivs.
        if not mixed_XQ:
            VQxx = dot(xQ, Vxxx, axes=[[1, 0]])
        else:
            VQxx = Vxxx

        V_QQQ_3 = dot(xQ, dot(xQ, VQxx, axes=[[1, 1]]), axes=[[1, 2]])

        V_QQQ_terms = (
            V_QQQ_1,
            V_QQQ_2,
            V_QQQ_3
        )
        V_QQQ = sum(x for x in V_QQQ_terms if not isinstance(x, int))

        derivs[2] = V_QQQ
        if order == 3:
            return tuple(derivs)

        # Fourth Derivs
        # For now we'll just generate everything rather than being particularly clever about it

        xQQQQ = x_derivs[3]
        Vxxxx = V_derivs[3]

        ## Gradient contribution
        V_QQQQ_1 = dot(xQQQQ, Vx)

        ## Hessian contribution
        #  All QQQ x Q permutations
        Q4231 = dot(xQQQ, dot(xQ, Vxx, axes=[[1, 0]]), axes=[[3, 1]])
        if not isinstance(Q4231, int):
            X = tuple(range(4, Q4231.ndim))
            V_QQQQ_21_terms =[
                Q4231,
                Q4231.transpose(3, 0, 1, 2, *X),
                Q4231.transpose(0, 3, 1, 2, *X),
                Q4231.transpose(0, 1, 3, 2, *X)
            ]
            V_QQQQ_21 = sum(V_QQQQ_21_terms)
        else:
            V_QQQQ_21_terms = 0
            V_QQQQ_21 = 0
        # QQ x QQ permutations
        Q4222 = dot(xQQ, dot(xQQ, Vxx, axes=[[2, 1]]), axes=[[2, 2]])
        if not isinstance(Q4222, int):
            X = tuple(range(4, Q4222.ndim))
            V_QQQQ_22_terms = [
                Q4222.transpose(2, 0, 1, 3, *X),
                Q4222.transpose(2, 3, 0, 1, *X),
                Q4222.transpose(2, 0, 3, 1, *X)
            ]
            V_QQQQ_22 = sum(V_QQQQ_22_terms)
        else:
            V_QQQQ_22 = 0

        V_QQQQ_2 = sum(x for x in [V_QQQQ_21, V_QQQQ_22] if not isinstance(x, int))

        Q4321 = dot(xQ, dot(xQQ, VQxx, axes=[[2, 2]]), axes=[[1, 3]])
        if not isinstance(Q4321, int):
            X = tuple(range(4, Q4321.ndim))
            V_QQQQ_3_terms = [
                Q4321.transpose(3, 1, 2, 0, *X),
                Q4321.transpose(1, 3, 2, 0, *X),
                Q4321.transpose(0, 1, 3, 2, *X),
                Q4321.transpose(0, 3, 1, 2, *X),
                Q4321.transpose(1, 0, 3, 2, *X)
                ]
            if mixed_terms:
                # bad name, but means that the Q in VQxx is different from the Q we're interested in,
                # but is equivalent to first order, so for numerical stability reasons rather than
                # use qQQ we just add on the appropriate transposition
                V_QQQQ_3_terms.append(Q4321.transpose(2, 1, 3, 0, *X))
            V_QQQQ_3 = sum(V_QQQQ_3_terms)
        else:
            V_QQQQ_3 = 0

        # fourth derivs
        if not mixed_XQ:
            VQQxx = dot(xQ, dot(xQ, Vxxxx), axes=[[1, 1]])
        else:
            VQQxx = Vxxxx

        if not isinstance(VQQxx, int):

            V_QQQQ_4 = dot(VQQxx, xQ, xQ, axes=[[3, 1], [2, 1]])

            N = V_QQQQ_4.ndim
            X = N - 4
            if X > 0:
                unroll = (0, 1) + tuple(range(2+X, N)) + tuple(range(2, 2+X))
                V_QQQQ_4 = V_QQQQ_4.transpose(unroll)
            X = tuple(range(4, V_QQQQ_4.ndim))
            # if mixed_XQ:
            #     # we need to zero out the elements we don't really have because Gaussian is mean
            #     import itertools
            #     nQ = V_QQQQ_4.shape[0]
            #     if nQ > 3:
            #         perms = np.array(list(itertools.permutations(range(nQ), 4))).T
            #         # print(V_QQQQ_4.shape[0], perms)
            #         V_QQQQ_4[perms] = 0.
        else:
            V_QQQQ_4 = 0

        V_QQQQ = (
                V_QQQQ_1 +
                V_QQQQ_2 +
                V_QQQQ_3 +
                V_QQQQ_4
        )

        return V_Q, V_QQ, V_QQQ, V_QQQQ

    _cached_transforms = {}
    internal_by_cartesian_order=3
    cartesian_by_internal_order=4
    jacobian_warning_threshold=1e12
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
                or len(current_cache["CartesiansByInternals"]) < cartesian_by_internal_order
                or len(current_cache["InternalsByCartesians"]) < internal_by_cartesian_order
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

            if (
                    "CartesiansByInternals" not in current_cache
                    or len(current_cache["CartesiansByInternals"]) < cartesian_by_internal_order
            ):
                # For speed reasons we've introduced class-level caching of these terms
                if self.logger is not None:
                    start = time.time()
                    self.logger.log_print(
                        "Getting d^nX/dR^n up to order {o}...",
                        o=cartesian_by_internal_order
                    )
                internal_jacobs = self.get_int_jacobs(list(range(1, cartesian_by_internal_order+1)))
                for i,x in enumerate(internal_jacobs):
                    bad_spots = np.where(np.abs(x) > self.jacobian_warning_threshold)
                    bad_bad_spots = bad_spots # so we don't lose it
                    if len(bad_spots) > 0: # numpy fuckery
                        bad_spots = bad_spots[0]
                    if len(bad_spots) > 0:
                        m = np.max(np.abs(x[bad_bad_spots]))
                        self.logger.log_print('WARNING: maximum d^{i}X/dR^{i} term is {m}. '
                                              'This will likely mess up G-matrix terms and is probably coming from a planar structure. '
                                              'Setting to zero, but `jacobian_warning_threshold` can be increased if this is expected',
                                              i=i,
                                              m=m
                                              )
                        x[bad_bad_spots] = 0.
                if self.logger is not None:
                    end = time.time()
                    self.logger.log_print(
                        "took {t}s",
                        t=round(end-start, 3)
                    )

                # The finite difference preserves too much shape by default
                _contract_dim = DumbTensor._contract_dim
                _ = []
                for i,x in enumerate(internal_jacobs):
                    if isinstance(x, int):
                        _.append(x)
                    elif x.ndim > 2+i:
                        _.append(_contract_dim(x, 2+i))
                internal_jacobs = _

                # we'll strip off the embedding coords just in case
                embedding_coords = [0, 1, 2, 4, 5, 8]
                good_coords = np.setdiff1d(np.arange(3*len(self.masses)), embedding_coords)
                # good_coords = np.arange(3*len(self.masses))

                # Need to then mass weight
                masses = self.masses
                mass_conv = np.sqrt(self._tripmass(masses))
                # mass weight the derivs w.r.t internals
                internal_weighting = mass_conv
                _ = []
                for i, x in enumerate(internal_jacobs):
                    internal_weighting = np.expand_dims(internal_weighting, 0)
                    if isinstance(x, int):
                        _.append(x)
                    else:
                        x = x * internal_weighting
                        for j in range(i+1):
                            x = np.take(x, good_coords, axis=j)
                        _.append(x)
                internal_jacobs = _

                current_cache["CartesiansByInternals"] = internal_jacobs
            else:
                internal_jacobs = current_cache["CartesiansByInternals"]

            if (
                    "InternalsByCartesians" not in current_cache
                    or len(current_cache["InternalsByCartesians"]) < internal_by_cartesian_order
            ):
                if self.logger is not None:
                    start = time.time()
                    self.logger.log_print(
                        "Getting d^nR/dX^n up to order {o}...",
                        o=internal_by_cartesian_order
                    )
                cartesian_jacobs = self.get_cart_jacobs(list(range(1, internal_by_cartesian_order + 1)))
                m = np.max([np.max(np.abs(x)) for x in cartesian_jacobs])
                for i,x in enumerate(cartesian_jacobs):
                    bad_spots = np.where(np.abs(x) > self.jacobian_warning_threshold)
                    bad_bad_spots = bad_spots # so we don't lose it
                    if len(bad_spots) > 0: # numpy fuckery
                        bad_spots = bad_spots[0]
                    if len(bad_spots) > 0:
                        m = np.max(np.abs(x[bad_bad_spots]))
                        self.logger.log_print('WARNING: maximum d^{i}R/dX^{i} term is {m}. '
                                              'This will likely mess up G-matrix terms and is probably coming from a planar structure. '
                                              'Setting to zero, but `jacobian_warning_threshold` can be increased if this is expected',
                                              i=i, m=m)
                        x[bad_bad_spots] = 0.
                if self.logger is not None:
                    end = time.time()
                    self.logger.log_print(
                        "took {t}s",
                        t=round(end-start, 3)
                    )

                _contract_dim = DumbTensor._contract_dim
                _ = []
                for i,x in enumerate(cartesian_jacobs):
                    if isinstance(x, int):
                        _.append(x)
                    elif x.ndim > 2+i:
                        _.append(_contract_dim(x, 2+i))
                cartesian_jacobs = _

                # we'll strip off the embedding coords just in case
                embedding_coords = [0, 1, 2, 4, 5, 8]
                good_coords = np.setdiff1d(np.arange(3*len(self.masses)), embedding_coords)
                # good_coords = np.arange(3*len(self.masses))

                # Need to then mass weight
                masses = self.masses
                mass_conv = np.sqrt(self._tripmass(masses))
                # mass weight the derivs w.r.t cartesians
                cartesian_weighting = mass_conv
                mc = mass_conv
                _ = []
                for i, x in enumerate(cartesian_jacobs):
                    cartesian_weighting = np.expand_dims(cartesian_weighting, -1)#[..., np.newaxis]
                    if isinstance(x, int):
                        _.append(x)
                    else:
                        x = x / cartesian_weighting
                        x = np.take(x, good_coords, axis=-1)
                        _.append(x)
                    mc = np.expand_dims(mc, 0)
                    cartesian_weighting = cartesian_weighting * mc
                cartesian_jacobs = _

                current_cache["InternalsByCartesians"] = cartesian_jacobs
            else:
                cartesian_jacobs = current_cache["InternalsByCartesians"]

            QY = self.modes.matrix  # derivatives of Q with respect to the Cartesians
            YQ = self.modes.inverse # derivatives of Cartesians with respect to Q

            if "InternalsByModes" not in current_cache:
                RQ, = TensorDerivativeConverter([YQ], cartesian_jacobs).convert(order=1, check_arrays=True)
                current_cache["InternalsByModes"] = [RQ]
            else:
                RQ, = current_cache["InternalsByModes"]

            if "ModesByInternals" not in current_cache:
                QR, = TensorDerivativeConverter(internal_jacobs, [QY]).convert(order=1, check_arrays=True)
                current_cache["ModesByInternals"] = [QR]
            else:
                QR, = current_cache["ModesByInternals"]

            if (
                    "CartesiansByModes" not in current_cache
                    or len(current_cache["CartesiansByModes"]) < len(internal_jacobs)
            ):
                x_derivs = internal_jacobs#(YR, YRR, YRRR, YRRRR)
                Q_derivs = [RQ] + [0]*(len(internal_jacobs) - 1)
                YQ_derivs = TensorDerivativeConverter(Q_derivs, x_derivs,
                                                      jacobians_name='Q',
                                                      values_name='X'
                                                      ).convert(order=len(internal_jacobs), check_arrays=True)

                qQ_derivs = TensorDerivativeConverter(YQ_derivs, [QY] + [0] * (len(internal_jacobs) - 1),
                                                      jacobians_name='Yq',
                                                      values_name='qY'
                                                      ).convert(order=len(internal_jacobs), check_arrays=True)
                # self._get_tensor_derivs(
                #     YQ_derivs, (QY, 0, 0, 0),
                #     mixed_XQ=False
                # )

                current_cache["CartesiansByModes"] = YQ_derivs
                current_cache["CartesianModesByInternalModes"] = qQ_derivs

                # "CartesiansByModes": [YQ, YQQ, YQQQ, YQQQQ],
                # "ModesByCartesians": [QY, QYY, QYYY],
                # "CartesianModesByInternalModes": [qQ, qQQ, qQQQ, qQQQQ]

            if (
                    "ModesByCartesians" not in current_cache
                    or len(current_cache["ModesByCartesians"]) < len(cartesian_jacobs)
            ):
                QY_derivs = TensorDerivativeConverter(cartesian_jacobs, [QR] + [0]*(len(cartesian_jacobs) - 1)).convert(order=len(cartesian_jacobs), check_arrays=True)
                current_cache["ModesByCartesians"] = QY_derivs

            # transf_data = {
            #     "CartesiansByInternals": internal_jacobs,
            #     "InternalsByCartesians": cartesian_jacobs,#[RY, RYY, RYYY],
            #     "InternalsByModes": [RQ],
            #     "CartesiansByModes": [YQ, YQQ, YQQQ, YQQQQ],
            #     "ModesByCartesians": [QY, QYY, QYYY],
            #     "CartesianModesByInternalModes": [qQ, qQQ, qQQQ, qQQQQ]
            # }

            self._cached_transforms[self.molecule] = current_cache
            self.checkpointer['coordinate_transforms'] = current_cache

        return current_cache#self._cached_transforms[self.molecule]

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
        )['CartesiansByModes']
        if order is not None:
            base = base[:order]
        return base

    @property
    def modes_by_cartesians(self):
        return self.get_coordinate_transforms()['ModesByCartesians']
    def get_modes_by_cartesians(self, order=None):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=None if order is None else min(order, self.cartesian_by_internal_order),
            internal_by_cartesian_order=order
        )['ModesByCartesians']
        if order is not None:
            base = base[:order]
        return base
    @property
    def cartesians_by_internals(self):
        return self.get_coordinate_transforms()['CartesiansByInternals']
    def get_cartesians_by_internals(self, order=None):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=order,
            internal_by_cartesian_order=None if order is None else min(order, self.internal_by_cartesian_order)
        )['CartesiansByInternals']
        if order is not None:
            base = base[:order]
        return base
    @property
    def internals_by_cartesians(self):
        return self.get_coordinate_transforms()['InternalsByCartesians']
    def get_internals_by_cartesians(self, order=None):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=None if order is None else min(order, self.cartesian_by_internal_order),
            internal_by_cartesian_order=order
        )['InternalsByCartesians']
        if order is not None:
            base = base[:order]
        return base

    @property
    def cartesian_modes_by_internal_modes(self):
        return self.get_coordinate_transforms()['CartesianModesByInternalModes']
    def get_cartesian_modes_by_internal_modes(self, order=None):
        base = self.get_coordinate_transforms(
            cartesian_by_internal_order=None if order is None else min(order, self.cartesian_by_internal_order),
            internal_by_cartesian_order=order
        )['CartesianModesByInternalModes']
        if order is not None:
            base = base[:order]
        return base

class PotentialTerms(ExpansionTerms):
    """
    A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates
    """
    def __init__(self,
                 molecule,
                 mixed_derivs=None,
                 modes=None,
                 potential_derivatives=None,
                 mode_selection=None,
                 logger=None,
                 parallelizer=None,
                 checkpointer=None
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
                         logger=logger, parallelizer=parallelizer, checkpointer=checkpointer)
        if potential_derivatives is None:
            potential_derivatives = molecule.potential_surface.derivatives
        self.mixed_derivs = mixed_derivs # we can figure this out from the shape in the future
        self.v_derivs = self._canonicalize_derivs(self.freqs, self.masses, potential_derivatives)

    def _canonicalize_derivs(self, freqs, masses, derivs):

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

        n = len(masses)
        modes_n = len(self.modes.freqs)
        internals_n = 3 * n - 6
        coord_n = 3 * n

        if len(derivs) > 2 and self.mode_sel is not None and thirds.shape[0] == internals_n:
            # TODO: need to handle more cases of input formats...
            thirds = thirds[(self.mode_sel,)]

        if len(derivs) > 3 and self.mode_sel is not None and fourths.shape[0] == internals_n:
            # TODO: need to handle more cases of input formats...
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

    hessian_tolerance=1.0e-4
    grad_tolerance=1.0e-4
    freq_tolerance=2e-3
    def old_get_terms(self, order=None, logger=None):

        if logger is None:
            logger = self.logger

        logger.log_print('calculating potential derivatives')

        if order is None:
            order = len(self.v_derivs)
        else:
            order += 2

        grad = self.v_derivs[0]
        # hess = self.v_derivs[1]

        # raise Exception([x.shape for x in self.v_derivs])

        # Use the Molecule's coordinates which know about their embedding by default
        intcds = self.internal_coordinates
        if intcds is None:
            # this is nice because it eliminates most of the terms in the expansion
            xQ = self.modes.inverse

            x_derivs = [xQ] + [0] * (order-1)
            V_derivs = self.v_derivs

            # terms = self._get_tensor_derivs(x_derivs, V_derivs, mixed_terms=False, mixed_XQ=self.mixed_derivs)

            if self.mixed_derivs:
                terms = TensorDerivativeConverter(x_derivs, V_derivs, mixed_terms=[
                    [None, self.v_derivs[2]], # dVdQXX
                    [None, self.v_derivs[3]]  # dVdQQXX
                ]).convert(order=order, check_arrays=True)
            else:
                terms = TensorDerivativeConverter(x_derivs, V_derivs).convert(order=order, check_arrays=True)

        else:
            x_derivs = self.get_cartesians_by_modes(order=order-1)
            x_derivs = list(x_derivs) + [0] # gradient term never matters

            if self.grad_tolerance is not None:
                if np.linalg.norm(grad) > self.grad_tolerance:
                    # add some logger stuff...
                    logger.log_print(
                        "WARNING: gradient norm is {n}",
                        n = np.linalg.norm(grad)
                    )
                    # grad = np.zeros(grad.shape)

            V_derivs = self.v_derivs

            if self.mixed_derivs:
                if order > 4:
                    raise ValueError("don't currently have things tested for expansions beyond 4th V derivatives with mixed derivatives") #TODO: relax this once we have more flexible input determination
                # terms = self._get_tensor_derivs(x_derivs, V_derivs, mixed_terms=True, mixed_XQ=self.mixed_derivs)

                # since the normal modes are expressed over
                # different sets of coordinates the fourth mixed deriv
                # terms need to be partially corrected
                QY, = self.get_modes_by_cartesians(1)
                qQQ = np.tensordot(x_derivs[1], QY, axes=[2, 0])
                f43 = np.tensordot(qQQ, V_derivs[2], axes=[2, 0])
                fourths = V_derivs[3] + f43
                V_derivs = V_derivs[:3] + [fourths]

                terms = TensorDerivativeConverter(x_derivs, V_derivs,
                                                  mixed_terms=[
                                                      [None, V_derivs[2]],  # dVdQXX
                                                      [None, V_derivs[3]]  # dVdQQXX
                                                  ]
                                                  ).convert(order=order, check_arrays=True)
            else:
                terms = TensorDerivativeConverter(x_derivs, V_derivs).convert(order=order, check_arrays=True)

            xQ2 = self.modes.inverse
            _, v2x,  =  TensorDerivativeConverter((xQ2, 0), V_derivs).convert(order=2, check_arrays=True)#self._get_tensor_derivs((xQ2, 0, 0, 0), V_derivs, order=2, mixed_XQ=False)

            if self.hessian_tolerance is not None:
                v2 = terms[1]
                v2_diff = v2 - v2x
                if True or np.max(np.abs(v2_diff)) > self.hessian_tolerance:
                    raise PerturbationTheoryException(
                        (
                            "Internal normal mode Hessian differs from Cartesian normal mode Hessian."
                            " Cartesian frequencies are {}, internals are {}.\n"
                            " This often indicates issues with the derivatives.\n"
                            " (YQ min/max: {} {} generally in the 10s for well-behaved systems)\n"
                            " (YQQ min/max: {} {} generally in the 10s for well-behaved systems)"
                         ).format(
                            np.diag(v2x)*UnitsData.convert("Hartrees", "Wavenumbers"),
                            np.diag(v2)*UnitsData.convert("Hartrees", "Wavenumbers"),
                            np.min(x_derivs[0]), np.max(x_derivs[0]),
                            np.min(x_derivs[1]), np.max(x_derivs[1])
                        )
                    )

        terms = terms[1:]

        if order > 3:
            v4 = terms[2]
            if self.mixed_derivs:# and intcds is None:
                # we assume we only got second derivs in Q_i Q_i
                # at this point, then, we should be able to fill in the terms we know are missing
                if not isinstance(v4, np.ndarray):
                    v4 = v4.asarray()
                for i in range(v4.shape[0]):
                    v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]

        self.checkpointer['potential_terms'] = terms

        new_freqs = np.diag(terms[0])
        old_freqs = self.modes.freqs
        # deviation on the order of a wavenumber can happen in low-freq stuff from numerical shiz
        if self.freq_tolerance is not None:
            if np.max(np.abs(new_freqs - old_freqs)) > self.freq_tolerance:
                raise PerturbationTheoryException(
                    "Force constants in normal modes don't return frequencies along diagonal;"
                    " this likely indicates issues with the mass-weighting"
                    " got {} but expected {}".format(new_freqs, old_freqs)
                )

        return terms

    get_terms = old_get_terms

class KineticTerms(ExpansionTerms):
    """Represents the KE coefficients"""

    g_derivative_threshold = 1e-3
    def get_terms(self, order=None, logger=None):

        if logger is None:
            logger = self.logger
        logger.log_print('calculating G-matrix derivatives')

        dot = DumbTensor._dot
        shift = DumbTensor._shift
        intcds = self.internal_coordinates
        if intcds is None:
            # this is nice because it eliminates a lot of terms in the expansion
            J = self.modes.matrix
            G = dot(J, J, axes=[[0, 0]])
            if order == 0:
                terms = [G]
            else:
                terms = [G] + [0]*(order)

        else:
            # should work this into the new layout
            QY_derivs = self.get_modes_by_cartesians(order=order+1)
            YQ_derivs = self.get_cartesians_by_modes(order=order+1)

            # RQ = dot(YQ, RY)

            term_getter = TensorDerivativeConverter(YQ_derivs, QY_derivs).terms
            term_getter.v_name = 'Y'
            J = term_getter.XV(1)
            G_terms = [J.dot(J, 1, 1)]
            for i in range(1, order+1):
                g_cur = G_terms[-1].dQ().simplify()
                G_terms.append(g_cur)
            terms = [x.array for x in G_terms]

            for i,t in enumerate(terms):
                if i == 0:
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

            # raise Exception([G, GQ, GQQ])

        G_terms = terms
        self.checkpointer['gmatrix_terms'] = G_terms

        return G_terms

class DipoleTerms(ExpansionTerms):
    def __init__(self,
                 molecule,
                 derivatives=None,
                 mixed_derivs=None,
                 modes=None,
                 mode_selection=None,
                 logger=None,
                 parallelizer=None,
                 checkpointer=None
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
                         logger=logger, parallelizer=parallelizer, checkpointer=checkpointer)
        self.mixed_derivs = mixed_derivs
        if self.mixed_derivs is None:
            self.mixed_derivs = mixed_derivs
        if derivatives is None:
            derivatives = molecule.dipole_surface.derivatives
        self.derivs = self._canonicalize_derivs(self.freqs, self.masses, derivatives)

    def _canonicalize_derivs(self, freqs, masses, derivs):
        """
        Makes sure all of the dipole moments are clean and ready to rotate
        """

        mom, grad, seconds, thirds = derivs
        try:
            grad = grad.array
        except AttributeError:
            pass

        n = len(masses)
        modes_n = len(self.modes.freqs)
        internals_n = 3 * n - 6
        coord_n = 3 * n

        if self.mode_sel is not None and thirds.shape[0] == internals_n:
            # TODO: need to handle more cases of input formats...
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
                seconds.shape != (modes_n, coord_n, 3)
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
                thirds.shape != (modes_n, modes_n, coord_n, 3)
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

        # We need to mass-weight the pure cartesian derivs
        # & undimensionalize the ones in terms of normal modes

        # amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        m_conv = np.sqrt(self._tripmass(masses))
        f_conv = np.sqrt(freqs)

        if grad.shape == (coord_n, 3):
            grad = grad / m_conv[:, np.newaxis]
        elif grad.shape == (modes_n, 3):
            grad = grad / f_conv[:, np.newaxis]

        if seconds.shape == (modes_n, coord_n, 3):
            # raise Exception('wat')
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

        return mom, grad, seconds, thirds

    def get_terms(self, order=None):
        if order is None:
            order = 2
        if order > 2:
            raise ValueError("only have this up order 2...")
        v0 = self.derivs[0]
        grad = self.derivs[1]
        seconds = self.derivs[2]
        thirds = self.derivs[3]

        # Use the Molecule's coordinates which know about their embedding by default
        intcds = self.internal_coordinates
        if intcds is None:# or not self.non_degenerate:
            # this is nice because it eliminates most of terms in the expansion
            xQ = self.modes.inverse
            xQQ = 0
            xQQQ = 0
            xQQQQ = 0
        else:

            QY, QYY, QYYY = self.modes_by_cartesians
            xQ, xQQ, xQQQ, xQQQQ = self.cartesians_by_modes

        x_derivs = (xQ, xQQ, xQQQ, xQQQQ)
        mu = [None]*3
        for coord in range(3):
            u_derivs = (grad[..., coord], seconds[..., coord], thirds[..., coord])
            # raise Exception(intcds, self.mixed_derivs)
            if intcds is not None and self.mixed_derivs:
                xQ, xQQ, xQQQ, xQQQQ = [DumbTensor(x) for x in x_derivs]
                u1, u2, u3 = [DumbTensor(u) for u in u_derivs]

                v1 = (xQ@u1).t
                v2 = xQ.dot(u2, axes=[[1, 1]]).t + xQQ.dot(u1, axes=[[2, 0]]).t

                # rather than do the naive total transformation, for numerical stability reasons
                # we'll want to directly apply xQQ since we have YQ QY _should_ be the identity
                # but is not since we're working in a subspace
                v3_1 = xQQQ.dot(u1, axes=[[3, 0]]).t

                v3_2 = xQQ.dot(u2, axes=[[2, 1]]).t
                v3_2 = 1/6 * sum(v3_2.transpose(p) for p in ip.permutations([0, 1, 2]))

                # raise Exception(v3_2 - v3_2.transpose([1, 2, 0]))

                v3_3 = xQ.dot(u3, axes=[[1, 2]]).t
                # for p in ip.permutations([0, 1, 2]):
                #     # we don't have enough data to actually determine anything where i!=j!=k
                #     v3_3[p] = 0.
                v3 = v3_1 + v3_2 + v3_3

            else:
                u1, u2, u3 = u_derivs
                v1 = np.tensordot(xQ, u1, axes=[1, 0])
                v2 = np.tensordot(xQ, u2, axes=[1, 1])
                v3 = np.tensordot(xQ, u3, axes=[1, 2])

            # print(">>>>>", v2)

            # raise Exception(v1, v2, v3)

            mu[coord] = (v0[coord], v1, v2, v3)#(v1, v2, 0)#(v1, 0, 0)


        self.checkpointer['dipole_terms'] = mu

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
        J = xQ.reshape((len(xQ), self.molecule.num_atoms, 3)).transpose(1, 0, 2)

        # then rotate into the inertial frame
        B_e, eigs = self.inertial_frame
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

    def get_terms(self, order=None):

        if order > 0:
            raise ValueError("Only have coriolis up to order 0 for now")

        zeta_inert, B_e = self.get_zetas_and_momi()

        # now we include the frequency dimensioning that comes from the q and p terms in Pi = Zeta*qipj
        freqs = self.freqs
        freq_term = np.sqrt(freqs[np.newaxis, :] / freqs[:, np.newaxis])
        zeta_inert = zeta_inert * freq_term[np.newaxis]

        coriolis = (
                           zeta_inert[:, :, :, np.newaxis, np.newaxis] # ij
                           * zeta_inert[:, np.newaxis, np.newaxis, :, :] # kl
        )

        corr = B_e[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis] * coriolis

        self.checkpointer['coriolis_terms'] = (corr[0], corr[1], corr[2])

        return [[corr[0], corr[1], corr[2]]]

class PotentialLikeTerm(KineticTerms):
    """
    This accounts for the potential-like term.
    In Cartesian diplacement modes this is the Watson U.
    In proper internals, this is the V' term.
    """

    def get_terms(self, order=None, logger=None):

        ics = self.internal_coordinates
        if ics is None:

            if order > 0:
                raise ValueError("currently only have Watson term up to 0 order")
            B_e, eigs = self.inertial_frame
            wat = -2*sum(B_e)
            wat_terms = [wat]

        else:
            ### transform inertia derivs into mode derivs
            YQ_derivs = self.get_cartesians_by_modes(order=2+order)

            I0_derivs = self.inertial_frame_derivatives() # only ever two of these
            if order > 0:
                I0_derivs = I0_derivs + [0]*order
            I0Q_derivs = TensorDerivativeConverter(YQ_derivs, I0_derivs).convert(check_arrays=True)

            ### pull already computed G-matrix derivs
            G_terms = super().get_terms(order=2+order, logger=NullLogger())

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
            gamdQ = (detI.dQ()/detI + -1*detG.dQ()/detG).simplify(check_arrays=True)
            gamdQQ = gamdQ.dQ().simplify(check_arrays=True)

            v0 = (
                    g_terms.QX(0).dot(gamdQQ, [1, 2], [1, 2])
                    + g_terms.QX(1).dot(gamdQ, 3, 1).tr()
                    + 1/4 * gamdQ.dot(gamdQ.dot(g_terms.QX(0), 1, 1), 1, 1)
            )

            wat_terms = [v0]
            for i in range(1, order+1):
                wat_terms.append(wat_terms[-1].dQ())
            wat_terms = [w.array for w in wat_terms]

        self.checkpointer['psuedopotential_terms'] = wat_terms

        return wat_terms
