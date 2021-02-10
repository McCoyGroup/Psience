"""
Stores all of the terms used inside the VPT2 representations
"""

import numpy as np, functools as fp, itertools as ip, scipy, time

from McUtils.Numputils import SparseArray, levi_cevita3, vec_tensordot, vec_outer
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer
from McUtils.Parallelizers import Parallelizer

from ..Molecools import MolecularVibrations, MolecularNormalModes

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
                 checkpointer=None
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
        self.internal_coordinates = molecule.internal_coordinates
        self.coords = molecule.coords
        self.masses = molecule.masses * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        if modes is None:
            modes = molecule.normal_modes
        if undimensionalize:
            self.raw_modes = modes
            modes = self.undimensionalize(self.masses, modes.basis)
        else:
            self.raw_modes = None
        if mode_selection is not None:
            modes = modes[(mode_selection,)]
        self.modes = modes
        self.mode_sel = mode_selection
        self.freqs = self.modes.freqs
        self._inert_frame = None

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

    def get_terms(self):
        raise NotImplemented

    @property
    def terms(self):
        if self._terms is None:
            self._terms = self.get_terms()
        return self._terms

    def __getitem__(self, item):
        return self.terms[item]

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
            new_jacs = [x.squeeze() for x in intcds.jacobian(carts, need_jacs, mesh_spacing=1.0e-2,
                                                             all_numerical=True,
                                                             converter_options=dict(reembed=False)
                                                             )]
            self._cached_jacobians[self.molecule]['int'] = new_jacs
            exist_jacs = new_jacs
        return [exist_jacs[j-1] for j in jacs]

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
            need_jacs = [x + 1 for x in range(0, max_jac)]
            new_jacs = [
                x.squeeze() for x in ccoords.jacobian(internals, need_jacs, mesh_spacing=1.0e-5, analytic_deriv_order=1)
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
    def get_coordinate_transforms(self):

        if self.molecule in self._cached_transforms:
            return self._cached_transforms[self.molecule]

        # ncoords = self.molecule.coords.shape[0]
        # gc = np.delete(
        #     np.arange(ncoords*3),
        #     [0, 1, 2, 4, 5, 8]
        # )

        # raise Exception(gc)

        if self.logger is not None:
            self.logger.log_print(
                [
                    "Getting coordinate transforms for {m}",
                    "Embedding axes: {a}"
                    ],
                m=self.molecule,
                a=self.internal_coordinates.system.converter_options["axes_labels"]
            )

        # For speed reasons we've introduced class-level caching of these terms
        if self.logger is not None:
            start = time.time()
            self.logger.log_print(
                "Getting d^nX/dR^n up to order {o}...",
                o=3
            )
        XR, XRR, XRRR = self.get_int_jacobs([1, 2, 3])
        XRRRR = 0
        if self.logger is not None:
            end = time.time()
            self.logger.log_print(
                "took {t}s",
                t=round(end-start, 3)
            )

        # The finite difference preserves too much shape by default
        _contract_dim = DumbTensor._contract_dim
        if XR.ndim > 2:
            XR = _contract_dim(XR, 2)
            # XR = XR[gc]
        if not isinstance(XRR, int):
            if XRR.ndim > 3:
                XRR = _contract_dim(XRR, 3)
            # XRR = XRR[np.ix_(gc, gc)]
        if not isinstance(XRRR, int):
            if XRRR.ndim > 4:
                XRRR = _contract_dim(XRRR, 4)
            # XRRR = XRRR[np.ix_(gc, gc, gc)]
        if not isinstance(XRRRR, int):
            if XRRRR.ndim > 5:
                XRRRR = _contract_dim(XRRRR, 5)
            # XRRRR = XRRRR[np.ix_(gc, gc, gc, gc)]

        if self.logger is not None:
            start = time.time()
            self.logger.log_print(
                "Getting d^nR/dX^n up to order {o}...",
                o=3
            )
        RX, RXX, RXXX = self.get_cart_jacobs([1, 2, 3])
        if self.logger is not None:
            end = time.time()
            self.logger.log_print(
                "took {t}s",
                t=round(end-start, 3)
            )
        # print(self.molecule, self.molecule.internal_coordinates, np.round(RX, 3), np.round(RXX, 3))
        # raise Exception(RXX.shape)

        if RX.ndim > 2:
            RX = _contract_dim(RX, 2)
        # RX = RX[:, gc]
        if RXX.ndim > 3:
            RXX = _contract_dim(RXX, 3)
        # RXX = RXX[:, :, gc]
        if RXXX.ndim > 4:
            RXXX = _contract_dim(RXXX, 4)
        # RXXX = RXXX[:, :, :, gc]

        # Need to then mass weight
        masses = self.masses
        mass_conv = np.sqrt(self._tripmass(masses))
        YR = XR * mass_conv[np.newaxis, :]
        if isinstance(XRR, int):
            YRR = 0
        else:
            YRR = XRR * mass_conv[np.newaxis, np.newaxis, :]
        if isinstance(XRRR, int):
            YRRR = 0
        else:
            YRRR = XRRR * mass_conv[np.newaxis, np.newaxis, np.newaxis, :]
        if isinstance(XRRRR, int):
            YRRRR = 0
        else:
            YRRRR = XRRRR * mass_conv[np.newaxis, np.newaxis, np.newaxis, np.newaxis, :]

        # We need to compute all these terms then mass weight them
        RY = RX / mass_conv[:, np.newaxis]
        RYY = RXX / (mass_conv[:, np.newaxis, np.newaxis] * mass_conv[np.newaxis, :, np.newaxis])
        RYYY = RXXX / (
                mass_conv[:, np.newaxis, np.newaxis, np.newaxis]
                * mass_conv[np.newaxis, :, np.newaxis, np.newaxis]
                * mass_conv[np.newaxis, np.newaxis, :, np.newaxis]
        )

        QY = self.modes.matrix  # derivatives of Q with respect to the Cartesians
        YQ = self.modes.inverse # derivatives of Cartesians with respect to Q

        RQ, = self._get_tensor_derivs((YQ,), (RY,), order=1, mixed_XQ=False)
        QR = np.tensordot(YR, QY, axes=[-1, 0])

        x_derivs = (YR, YRR, YRRR, YRRRR)
        Q_derivs = (RQ, 0, 0, 0)
        YQ_derivs = self._get_tensor_derivs(Q_derivs, x_derivs, mixed_XQ=False)
        YQ, YQQ, YQQQ, YQQQQ = YQ_derivs

        QYY = np.tensordot(RYY, QR, axes=[-1, 0])
        QYYY = np.tensordot(RYYY, QR, axes=[-1, 0])

        qQ, qQQ, qQQQ, qQQQQ = self._get_tensor_derivs(
            YQ_derivs, (QY, 0, 0, 0),
            mixed_XQ=False
        )

        transf_data = {
            "InternalsByModes": [RQ],
            "CartesiansByModes": [YQ, YQQ, YQQQ, YQQQQ],
            "ModesByCartesians": [QY, QYY, QYYY],
            "CartesiansByInternals": [YR, YRR, YRRR, YRRRR],
            "InternalsByCartesians": [RY, RYY, RYYY],
            "CartesianModesByInternalModes": [qQ, qQQ, qQQQ, qQQQQ]
        }

        self._cached_transforms[self.molecule] = transf_data
        self.checkpointer['coordinate_transforms'] = transf_data

        return self._cached_transforms[self.molecule]

    @property
    def cartesians_by_modes(self):
        return self.get_coordinate_transforms()['CartesiansByModes']

    @property
    def modes_by_cartesians(self):
        return self.get_coordinate_transforms()['ModesByCartesians']

    @property
    def cartesians_by_internals(self):
        return self.get_coordinate_transforms()['CartesiansByInternals']

    @property
    def internals_by_cartesians(self):
        return self.get_coordinate_transforms()['InternalsByCartesians']

    @property
    def cartesian_modes_by_internal_modes(self):
        return self.get_coordinate_transforms()['CartesianModesByInternalModes']

class PotentialTerms(ExpansionTerms):
    """
    A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates
    """
    def __init__(self,
                 molecule,
                 mixed_derivs=True,
                 modes=None,
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
        :type modes: None | MolecularNormalModes
        :param mode_selection: the subset of normal modes to use
        :type mode_selection: None | Iterable[int]
        """
        super().__init__(molecule, modes, mode_selection=mode_selection,
                         logger=logger, parallelizer=parallelizer, checkpointer=checkpointer)
        self.v_derivs = self._canonicalize_derivs(self.freqs, self.masses, molecule.potential_derivatives)
        self.mixed_derivs = mixed_derivs # we can figure this out from the shape in the future

    def _canonicalize_derivs(self, freqs, masses, derivs):

        if len(derivs) == 3:
            grad, fcs, fds = derivs
            fcs = fcs.array
            thirds = fds.third_deriv_array
            fourths = fds.fourth_deriv_array
        else:
            grad, fcs, thirds, fourths = derivs

        n = len(masses)
        # modes_matrix = self.modes.inverse
        # modes_n = len(modes_matrix)
        modes_n = 3*n - 6
        # if modes_n == 3*n:
        #     modes_n = modes_n - 6
        #     modes_matrix = modes_matrix[6:]
        #     freqs = freqs[6:]
        coord_n = 3*n
        if grad.shape != (coord_n,):
            raise PerturbationTheoryException(
                "{0}.{1}: length of gradient array ({2[0]}) is not {3[0]}".format(
                    type(self).__name__,
                    "_canonicalize_force_constants",
                    grad.shape,
                    (coord_n,)
                )
            )
        if fcs.shape != (coord_n, coord_n):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of force constant array ({2[0]}x{2[1]}) is not {3[0]}x{3[1]}".format(
                    type(self).__name__,
                    "_canonicalize_force_constants",
                    fcs.shape,
                    (coord_n, coord_n)
                )
            )
        if thirds.shape != (modes_n, coord_n, coord_n):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of third derivative array ({2[0]}x{2[1]}x{2[2]}) is not ({3[0]}x{3[1]}x{3[2]})".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    thirds.shape,
                    (modes_n, coord_n, coord_n)
                )
            )
        # this might need to change in the future
        if fourths.shape != (modes_n, modes_n, coord_n, coord_n):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of fourth derivative array ({2[0]}x{2[1]}x{2[2]}x{2[3]}) is not ({3[0]}x{3[1]}x{3[2]}x{3[3]})".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    fourths.shape,
                    (modes_n, modes_n, coord_n, coord_n)
                )
            )

        amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        m_conv = np.sqrt(self._tripmass(masses))
        f_conv = np.sqrt(freqs)
        # f_conv = np.ones(f_conv.shape) # debugging

        undimension_2 = np.outer(m_conv, m_conv)
        fcs = fcs / undimension_2

        if self.mode_sel is not None:
            thirds = thirds[(self.mode_sel,)]
        undimension_3 = np.outer(m_conv, m_conv)[np.newaxis, :, :] * f_conv[:, np.newaxis, np.newaxis]
        thirds = thirds * (1 / undimension_3 / np.sqrt(amu_conv))

        wat = np.outer(m_conv, m_conv)[np.newaxis, :, :] * (f_conv ** 2)[:, np.newaxis, np.newaxis]
        undimension_4 = SparseArray.from_diag(1 / wat / amu_conv)
        if self.mode_sel is not None:
            if not isinstance(self.mode_sel, slice):
                fourths = fourths[np.ix_(self.mode_sel, self.mode_sel)]
            else:
                fourths = fourths[self.mode_sel, self.mode_sel]
        fourths = fourths * undimension_4

        return grad, fcs, thirds, fourths

    def new_get_terms(self):
        # I'd have liked this to work...but it looks like numerics is fucking it up?
        # or lack of invertibility in QY/YQ?
        grad = self.v_derivs[0]
        hess = self.v_derivs[1]
        thirds = self.v_derivs[2]
        fourths = self.v_derivs[3]

        # transform into proper Cartesian mode derivatives first
        x_derivs = (self.modes.inverse, 0, 0, 0)
        V_derivs = (grad, hess, thirds, fourths)

        v1, v2, v3, v4 = self._get_tensor_derivs(x_derivs, V_derivs, mixed_XQ=self.mixed_derivs)
        if self.mixed_derivs:
            # we assume we only got second derivs in Q_i Q_i
            # at this point, then, we should be able to fill in the terms we know are missing
            for i in range(v4.shape[0]):
                v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]

        # Now if we've got an internal spec, transform into internal modes
        intcds = self.internal_coordinates
        if intcds is not None:
            YQ, YQQ, _, _ = self.cartesians_by_modes
            qQ_terms = self.cartesian_modes_by_internal_modes
            v1, v2, v3, v4 = self._get_tensor_derivs(qQ_terms, (v1, v2, v3, v4), mixed_XQ=False)

        return v2, v3, v4

    def old_get_terms(self):
        grad = self.v_derivs[0]
        hess = self.v_derivs[1]
        thirds = self.v_derivs[2]
        fourths = self.v_derivs[3]

        # Use the Molecule's coordinates which know about their embedding by default
        intcds = self.internal_coordinates
        if intcds is None:
            # this is nice because it eliminates most of terms in the expansion
            xQ = self.modes.inverse
            xQQ = 0
            xQQQ = 0
            xQQQQ = 0

            x_derivs = (xQ, xQQ, xQQQ, xQQQQ)
            V_derivs = (grad, hess, thirds, fourths)

            v1, v2, v3, v4 = self._get_tensor_derivs(x_derivs, V_derivs, mixed_XQ=self.mixed_derivs)
        else:

            xQ, xQQ, xQQQ, xQQQQ = self.cartesians_by_modes

            if np.linalg.norm(grad) > 1.0e-4:
                # add some logger stuff...
                if self.logger is not None:
                    self.logger.log_print("WARNING: gradient norm is {}".format(
                        np.linalg.norm(grad)
                    ))
                grad = np.zeros(grad.shape)

            x_derivs = (xQ, xQQ, xQQQ, xQQQQ)
            V_derivs = (grad, hess, thirds, fourths)

            v1, v2, v3, v4 = self._get_tensor_derivs(x_derivs, V_derivs, mixed_terms=True, mixed_XQ=self.mixed_derivs)

            xQ2 = self.modes.inverse
            _, v2x,  = self._get_tensor_derivs((xQ2, 0, 0, 0), V_derivs, order=2, mixed_XQ=self.mixed_derivs)

            v2_diff = v2 - v2x

            if np.max(np.abs(v2_diff)) > 1.0e-4:
                raise PerturbationTheoryException(
                    "Internal normal mode Hessian differs from Cartesian normal mode Hessian;"
                    " this likely indicates issues with the second derivatives"
                    " (YQQ min/max: {} {} generally in the 10s for well-behaved systems)".format(np.min(xQQ), np.max(xQQ))
                )

        if self.mixed_derivs:# and intcds is None:
            # we assume we only got second derivs in Q_i Q_i
            # at this point, then, we should be able to fill in the terms we know are missing
            for i in range(v4.shape[0]):
                v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]

        self.checkpointer['potential_terms'] = (v2, v3, v4)

        return v2, v3, v4


    get_terms = old_get_terms

class KineticTerms(ExpansionTerms):
    """Represents the KE coefficients"""

    def get_terms(self):

        dot = DumbTensor._dot
        shift = DumbTensor._shift
        intcds = self.internal_coordinates
        if intcds is None:
            # this is nice because it eliminates a lot of terms in the expansion
            J = self.modes.matrix
            G = dot(J, J, axes=[[0, 0]])
            GQ = 0
            GQQ = 0
        else:
            QY, QYY, QYYY = self.modes_by_cartesians
            YQ, YQQ, YQQQ, YQQQQ = self.cartesians_by_modes
            # YR, YRR, YRRR, YRRRRR = self.cartesians_by_internals
            RY, RYY, RYYY = self.internals_by_cartesians

            # RQ = dot(YQ, RY)

            G = dot(QY, QY, axes=[[0, 0]])

            # J = DumbTensor(QY)
            # Jd = DumbTensor(YQ)
            # K = DumbTensor(dot(RYY, YR, QY))
            # U = K.dot(J, axes=[[0, 0]])
            #
            # GQ = Jd@(U + U[2:1])
            # GQ = GQ.t
            #
            # L = DumbTensor(dot(RYYY, YR, QY))
            # H = DumbTensor(dot(RQ, dot(RQ, YRR, axes=[[1, 0]]), axes=[[1, 1]]))
            # K22 = K.dot(K, axes=[[1, 1]])
            # V = L[3:2]@J + K22[2:0]
            #
            # GQQ = (H@(U + U[2:1])).t + (Jd@(Jd@(V+V[3:2]))[0:1]).t
            # raise Exception(QYY.shape)

            GQ = dot(YQ, QYY, QY, axes=[[1, 0], [1, 0]]) + dot(YQ, dot(QY, QYY, axes=[[0, 0]]), axes=[[-1, 1]])
            GQQ = (
                    dot(YQQ, QYY, QY, axes=[[-1, 0], [2, 0]])
                    + dot(YQQ, dot(QY, QYY, axes=[[0, 0]]), axes=[[-1, 1]])
                    + dot(YQ, dot(YQ, QYYY, QY, axes=[[-1, 0], [1, 0]]), axes=[[1, 1]])
                    + dot(YQ, dot(YQ, dot(QY, QYYY, axes=[[0, 0]]), axes=[[1, 1]]), axes=[[1, 2]])
                    + dot(YQ, dot(YQ, QYY, QYY, axes=[[-1, 0], [1, 1]]), axes=[[1, 2]])
                    + dot(YQ, dot(YQ, QYY, QYY, axes=[[-1, 0], [1, 1]]), axes=[[1, 2]]).transpose((0, 1, 3, 2))
                    # + dot(YQ, dot(YQ, dot(QYY, QYY, axes=[[0, 0]]), axes=[[-1, 0]]), axes=[[1, 2]])
            )
            # for l, t in (
            #         ("G", G),
            #         ("dGdQ", GQ),
            #         ("dGdQQ", GQQ)
            # ):
            #     t = np.abs(t)
            #     M = np.max(t)
            #     m = np.min(t)
            #     print(l + ": ", np.round(M, 4), np.where(t == M), np.round(m, 4), np.average(t))

        G_terms = (G, GQ, GQQ)

        self.checkpointer['gmatrix_terms'] = G_terms

        return G_terms

class DipoleTerms(ExpansionTerms):
    def __init__(self,
                 molecule,
                 mixed_derivs=True,
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
        :type modes: None | MolecularNormalModes
        :param mode_selection: the subset of normal modes to use
        :type mode_selection: None | Iterable[int]
        """
        super().__init__(molecule, modes=modes, mode_selection=mode_selection,
                         logger=logger, parallelizer=parallelizer, checkpointer=checkpointer)
        self.derivs = self._canonicalize_derivs(self.freqs, self.masses, molecule.dipole_derivatives)
        self.mixed_derivs = mixed_derivs # we can figure this out from the shape in the future

    def _canonicalize_derivs(self, freqs, masses, derivs):
        """
        Makes sure all of the dipole moments are clean and ready to rotate
        """

        # TODO: this will need major clean up now that I've improved
        #       dipole handling in Molecule
        mom, grad, higher = derivs
        grad = grad.array
        seconds = higher.second_deriv_array
        thirds = higher.third_deriv_array

        n = len(masses)
        modes_matrix = self.modes.inverse
        modes_n = len(modes_matrix)
        if modes_n == 3*n:
            modes_n = modes_n - 6
            modes_matrix = modes_matrix[6:]
            freqs = freqs[6:]
        coord_n = modes_n + 6
        if grad.shape != (coord_n, 3):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of dipole derivative array ({2[0]}) is not {3[0]}".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    grad.shape,
                    (coord_n, 3)
                )
            )
        if seconds.shape != (modes_n, coord_n, 3):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of dipole second derivative array ({2[0]}x{2[1]}) is not {3[0]}x{3[1]}".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    seconds.shape,
                    (modes_n, coord_n, 3)
                )
            )
        if thirds.shape != (modes_n, modes_n, coord_n, 3):
            raise PerturbationTheoryException(
                "{0}.{1}: dimension of dipole third derivative array ({2[0]}x{2[1]}x{2[2]}) is not ({3[0]}x{3[1]}x{3[2]})".format(
                    type(self).__name__,
                    "_canonicalize_derivs",
                    thirds.shape,
                    (modes_n, modes_n, coord_n, 3)
                )
            )

        # We need to mass-weight the pure cartesian derivs
        # & undimensionalize the ones in terms of normal modes

        amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        m_conv = np.sqrt(self._tripmass(masses))
        f_conv = np.sqrt(freqs*amu_conv)

        grad = grad / m_conv[:, np.newaxis]
        seconds = seconds / (
            f_conv[:, np.newaxis, np.newaxis]
            * m_conv[np.newaxis, :, np.newaxis]
        )
        thirds = thirds / (
                f_conv[:, np.newaxis, np.newaxis, np.newaxis]
                * f_conv[np.newaxis, :, np.newaxis, np.newaxis]
                * m_conv[np.newaxis, np.newaxis, :, np.newaxis]
        )

        # undimension_2 = np.outer(m_conv, m_conv)
        # fcs = fcs / undimension_2
        #
        # undimension_3 = np.outer(m_conv, m_conv)[np.newaxis, :, :] * f_conv[:, np.newaxis, np.newaxis]
        # thirds = thirds * (1 / undimension_3 / np.sqrt(amu_conv))
        #
        # wat = np.outer(m_conv, m_conv)[np.newaxis, :, :] * (f_conv ** 2)[:, np.newaxis, np.newaxis]
        # undimension_4 = SparseArray.from_diag(1 / wat / amu_conv)
        # fourths = fourths
        # fourths = fourths * undimension_4

        return mom, grad, seconds, thirds

    def get_terms(self):
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

                # raise Exception([
                #     np.max(np.abs(v3_2.t)),
                #     np.max(np.abs(v3_3)),
                #     np.max(np.abs(v2)),
                #     np.max(np.abs(v1))
                # ])

                # ds = dict(plot_style={'vmin':-1.0e-5, 'vmax':1.0e-5})
                # import McUtils.Plots as plt
                # plt.TensorPlot(v3_21.t).show()
                # raise Exception("...")
                # plt.TensorPlot(u_derivs[2]).show()
                # # plt.ArrayPlot(u2.t)
                # # plt.ArrayPlot(v2.t)
                # # plt.TensorPlot(u_derivs[2])
                # v3_3 = v3_2.t
                # plt.TensorPlot(v3_3 - v3_3.transpose(1, 0, 2), **ds)
                # plt.TensorPlot(v3_3 - v3_3.transpose(0, 2, 1), **ds)
                # plt.TensorPlot(v3_3 - v3_3.transpose(2, 0, 1), **ds)
                # plt.TensorPlot(v3_3 - v3_3.transpose(2, 1, 0), **ds)
                # plt.TensorPlot(v3_3 - v3_3.transpose(1, 2, 0), **ds).show()
                #
                # raise Exception("...?")

                # for i in range(v3.shape[0]):
                #     v3[i, :, i] = v3[:, i, i] = v3[i, i, :]

            else:
                u1, u2, u3 = u_derivs
                v1 = np.tensordot(xQ, u1, axes=[1, 0])
                v2 = np.tensordot(xQ, u2, axes=[1, 1])
                v3 = np.tensordot(xQ, u3, axes=[1, 2])

            # print(">>>>>", v2)

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
        J = np.tensordot(J, eigs, axes=1)

        # coriolis terms are given by zeta = sum(JeJ^T, n)
        ce = -levi_cevita3
        zeta = sum(
            np.tensordot(
                np.tensordot(J[n], ce, axes=[1, 0]),
                J[n],
                axes=[2, 1]).transpose(1, 0, 2)
            for n in range(J.shape[0])
        )

        return zeta, B_e

    def get_zetas(self):

        z, m = self.get_zetas_and_momi()

        return z

    def get_terms(self):

        zeta_inert, B_e = self.get_zetas_and_momi()

        # new we include the frequency dimensioning that comes from the q and p terms in Pi = Zeta*qipj
        freqs = self.freqs
        freq_term = np.sqrt(freqs[np.newaxis, :] / freqs[:, np.newaxis])
        zeta_inert = zeta_inert * freq_term[np.newaxis]

        coriolis = (
                           zeta_inert[:, :, :, np.newaxis, np.newaxis] # ij
                           * zeta_inert[:, np.newaxis, np.newaxis, :, :] # kl
        )

        corr = B_e[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis] * coriolis

        self.checkpointer['coriolis_terms'] = (corr[0], corr[1], corr[2])

        return corr[0], corr[1], corr[2]

class PotentialLikeTerm(KineticTerms):
    """
    This accounts for the potential-like term.
    In Cartesian diplacement modes this is the Watson U.
    In proper internals, this is the V' term.
    """

    def get_terms(self):


        ics = self.internal_coordinates
        if ics is None:

            B_e, eigs = self.inertial_frame
            wat = -2*sum(B_e)

        else:
            # much more complicated, but we have
            # wat = sum(dGdQ_ii . dgdq) + G . dgdQQ + 1/4 (G.dgdq)^T . dgdq
            # where g = det(I_0) / get(G)

            mass = np.sqrt(self.masses)
            carts = mass[:, np.newaxis] * self.molecule.coords # mass-weighted Cartesian coordinates


            amu2me = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
            I0 = amu2me * self.molecule.inertia_tensor

            ### compute basic inertia tensor derivatives
            # first derivs are computed as a full (nAt, 3, I_rows (3), I_cols (3)) tensor
            # and then reshaped to (nAt * 3, I_rows, I_cols)
            eyeXeye = np.eye(9).reshape(3, 3, 3 ,3).transpose((2, 0, 1, 3))
            I0Y_1 = np.tensordot(carts, eyeXeye, axes=[1, 0])

            nAt = carts.shape[0]
            nY = nAt * 3
            I0Y_21 = np.reshape(np.eye(3), (9,))[np.newaxis, :, np.newaxis] * carts[:, np.newaxis, :] # a flavor of outer product
            I0Y_21 = I0Y_21.reshape((nAt, 3, 3, 3))
            I0Y_2 = (I0Y_21 + I0Y_21.transpose((0, 1, 3, 2)))
            I0Y = 2 * I0Y_1 - I0Y_2
            I0Y = I0Y.reshape((nY, 3, 3))

            # second derivatives are 100% independent of coorinates
            # only the diagonal blocks are non-zero, so we compute that block
            # and then tile appropriately
            keyXey = np.eye(9).reshape(3, 3, 3 ,3)
            I0YY_nn = 2 * eyeXeye - (keyXey + keyXey.transpose((0, 1, 3, 2)))
            I0YY = np.zeros((nAt, 3, nAt, 3, 3, 3))
            for n in range(nAt):
                I0YY[n, :, n, :, :, :] = I0YY_nn
            I0YY = I0YY.reshape((nY, nY, 3, 3))

            # import McUtils.Plots as plt
            # plt.TensorPlot(
            #     I0YY[:5, :5]
            # ).show()

            ### transform inertia derivs into mode derivs
            YQ, YQQ, YQQQ, YQQQQ = self.cartesians_by_modes

            I0Q = np.tensordot(YQ, I0Y, axes=[-1, 0])
            I0QQ = np.tensordot(YQQ, I0Y, axes=[-1, 0]) + np.tensordot(
                YQ,
                np.tensordot(
                    YQ,
                    I0YY,
                    axes=[-1, 0]
                ),
                axes=[-1, 1]
            )

            ### pull already computed G-matrix derivs
            G, GQ, GQQ = super().get_terms()


            # now build the actual dg/dQ terms

            detI = np.linalg.det(I0)
            detG = np.linalg.det(G)
            gam = detI / detG

            invI = np.linalg.inv(I0)
            invG = np.linalg.inv(G)

            adjI = invI*detI
            adjG = invG*detG

            invIdQ = - np.tensordot(np.tensordot(invI, I0Q, axes=[-1, 1]), invI, axes=[-1, 0]).transpose(1, 0, 2)
            invGdQ = - np.tensordot(np.tensordot(invG, GQ, axes=[-1, 1]), invG, axes=[-1, 0]).transpose(1, 0, 2)

            # not quite enough terms to want to be clever here...
            nQ = GQ.shape[0]
            ## First derivatives of the determinant
            detIdQ = np.array([
                np.trace(np.dot(adjI, I0Q[i]))
                for i in range(nQ)
            ])
            detGdQ = np.array([
                np.trace(np.dot(adjG, GQ[i]))
                for i in range(nQ)
            ])

            adjIdQ = detI * invIdQ + detIdQ[:, np.newaxis, np.newaxis] * invI[np.newaxis, :, :]
            adjGdQ = detG * invGdQ + detGdQ[:, np.newaxis, np.newaxis] * invG[np.newaxis, :, :]

            ## Second derivatives of the determinant
            detIdQQ = np.array([
                [
                    np.tensordot(I0Q[i], adjIdQ[j], axes=2)
                    + np.tensordot(adjI, I0QQ[i, j], axes=2)
                    for i in range(nQ)
                ]
                for j in range(nQ)
            ])
            detGdQQ = np.array([
                [
                    np.tensordot(GQ[i], adjGdQ[j], axes=2)
                    + np.tensordot(adjG, GQQ[i, j], axes=2)
                    for i in range(nQ)
                ]
                for j in range(nQ)
            ])

            ## Derivatives of Gamma
            gamdQ_I = 1/detI * detIdQ
            gamdQ_G = 1/detG * detGdQ
            gamdQ = gamdQ_I - gamdQ_G

            gamdQQ_I = -1 / detI**2 * np.outer(detIdQ, detIdQ) + 1 / detI * detIdQQ
            gamdQQ_G = -1 / detG**2 * np.outer(detGdQ, detGdQ) + 1 / detG * detGdQQ
            gamdQQ = gamdQQ_I - gamdQQ_G

            # Build out the proper Watson term
            wat_diag = sum(
                    np.dot(GQ[i, i], gamdQ)
                    for i in range(nQ)
                )
            wat_QQ = np.tensordot(G, gamdQQ, axes=2)
            wat_Q = np.tensordot(
                    np.tensordot(G, gamdQ, axes=[1, 0]),
                    gamdQ, axes=[0, 0]
                )
            wat = (wat_diag + wat_QQ + 1/4 * wat_Q)

            # import McUtils.Plots as plt
            #
            # plt.ArrayPlot(gamdQQ).show()

            # raise Exception([
            #     wat,
            #     [
            #         np.max(np.abs(x)) for x in [
            #         # wat_diag,
            #         # wat_Q,
            #         # wat_QQ,
            #         # # gamdQ_G,
            #         # gamdQQ_G,
            #         gam,
            #         detI,
            #         detIdQ,
            #         detIdQQ,
            #         I0,
            #         I0Q,
            #         I0QQ,
            #         gamdQ_I,
            #         gamdQQ_I,
            #         # detGdQQ
            #     ]
            #     ]
            #     # np.max(np.abs(I0)),
            #     # np.max(np.abs(I0Q)),
            #     # np.max(np.abs(gamdQ_I)),
            #     # np.max(np.abs(gamdQQ_G)),
            #     # np.max(np.abs(gamdQQ_I)),
            #
            # ])

            # wat = 0. # turn off for the moment

        # print(wat)

        self.checkpointer['psuedopotential_terms'] = [wat]

        return [wat]
