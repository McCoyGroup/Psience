"""
Stores all of the terms used inside the VPT2 representations
"""

import numpy as np, functools as fp, itertools as ip
from McUtils.Numputils import SparseArray
from McUtils.Data import UnitsData

from ..Molecools import Molecule, MolecularNormalModes
from .Common import PerturbationTheoryException

import McUtils.Plots as plt
import McUtils.Coordinerds as crds

__all__ = [
    "ExpansionTerms",
    "KineticTerms",
    "PotentialTerms"
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
    def __init__(self):
        self._terms = None

    @staticmethod
    def _tripmass(masses):
        return np.broadcast_to(masses, (len(masses), 3)).T.flatten()

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

    @classmethod
    def _get_tensor_derivs(cls, x_derivs, V_derivs, order=4, mixed_XQ=False):
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

        xQ_Vxx = dot(xQ, Vxx)
        V_QQ = dot(xQQ, Vx) + dot(xQ_Vxx, xQ, axes=[[1, 1]])
        derivs[1] = V_QQ
        if order == 2:
            return tuple(derivs)

        # Third Derivs
        xQQQ = x_derivs[2]
        Vxxx = V_derivs[2]

        # If Q is just an expansion in X all of these terms will disappear except for V_QQQ_5
        xQQ_Vxx = dot(xQQ, Vxx)

        V_QQQ_1 = dot(xQQQ, Vx)
        V_QQQ_2 = dot(xQ_Vxx, xQQ, axes=[[-1, -1]])
        V_QQQ_3 = dot(xQQ_Vxx, xQ, axes=[[-1, -1]])
        V_QQQ_4 = shift(V_QQQ_2, (1, 0))

        if not mixed_XQ:
            VQxx = dot(xQ, Vxxx)
        else:
            VQxx = Vxxx

        V_QQQ_5 = dot(VQxx, xQ, xQ, axes=[[2, 1], [1, 1]])

        V_QQQ_terms = (
            V_QQQ_1,
            V_QQQ_2,
            V_QQQ_3,
            V_QQQ_4,
            V_QQQ_5
        )
        # if mixed_XQ:
        #     ehhh = [
        #         UnitsData.convert("Hartrees", "Wavenumbers") * x for x in V_QQQ_terms
        #     ]
        #     np.savetxt("/Users/Mark/Desktop/eh.dat",
        #                np.array([e.flatten() for e in ehhh])
        #                )
        #     raise Exception(...)
        # print(np.reshape(VQxx, (3, 81)))
        V_QQQ = sum(x for x in V_QQQ_terms if not isinstance(x, int))
        derivs[2] = V_QQQ
        if order == 3:
            return tuple(derivs)

        # Fourth Derivs
        # For now we'll just generate everything rather than being particularly clever about it

        xQQQQ = x_derivs[3]
        Vxxxx = V_derivs[3]
        V_QQQQ_1 = dot(xQQQQ, Vx) + dot(xQ, Vxx, xQQQ, axes=[[-1, 0], [-1, -1]])

        xQQQ_Vxx_xQ = dot(xQQQ, Vxx, xQ, axes=[[-1, 0], [-1, -1]])
        xQ_22_Vxxx = dot(xQ, Vxxx, axes=[[-1, 1]])
        xQQ_Vxx_xQQ = dot(xQQ_Vxx, xQQ, axes=[[-1, -1]])

        V_QQQQ_2 = (
                xQQ_Vxx_xQQ +
                dot(Vxxx, xQ, xQQ, axes=[[1, -1], [1, -1]]) +
                shift(xQQQ_Vxx_xQ, (0, 1), (2, 3))
        )

        V_QQQQ_3 = (
                xQQQ_Vxx_xQ +
                dot(Vxxx, xQQ, xQ, axes=[[-1, -1], [1, -1]]) +
                shift(xQQ_Vxx_xQQ, (0, 3))
        )

        V_QQQQ_4 = (
                shift(xQQ_Vxx_xQQ, (1, 2)) +
                shift(dot(xQ, VQxx, xQQ, axes=[[1, 1], [2, 2]]), (2, 3)) +
                shift(xQQQ_Vxx_xQ, (0, 1), (3, 1))
        )

        if not mixed_XQ:
            VQQxx = dot(xQ, dot(xQ, Vxxxx), axes=[[1, 1]])
        else:
            VQQxx = Vxxxx

        V_QQQQ_5 = (
                dot(xQQ, VQxx, xQ, axes=[[2, 1], [3, 1]]) +
                shift(dot(xQQ, VQxx, xQ, axes=[[2, 2], [3, 1]]), (3, 1)) +
                shift(dot(xQ, VQxx, xQQ, axes=[[1, 1], [2, 2]]), (2, 0)) +
                dot(VQQxx, xQ, xQ, axes=[[3, 1], [2, 1]])
        )

        V_QQQQ = (
                V_QQQQ_1 +
                V_QQQQ_2 +
                V_QQQQ_3 +
                V_QQQQ_4 +
                V_QQQQ_5
        )

        return V_Q, V_QQ, V_QQQ, V_QQQQ


class PotentialTerms(ExpansionTerms):
    def __init__(self, molecule, mixed_derivs=True, non_degenerate=False):
        self.molecule = molecule
        self.internal_coordinates = molecule.internal_coordinates
        self.coords = molecule.coords
        self.masses = molecule.masses*UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        self.modes = self.undimensionalize(self.masses, molecule.normal_modes.basis)
        self.freqs = self.modes.freqs
        self.v_derivs = self._canonicalize_derivs(self.freqs, self.masses, molecule.potential_derivatives)
        self.non_degenerate=non_degenerate
        self.mixed_derivs = mixed_derivs # we can figure this out from the shape in the future
        super().__init__()

    def undimensionalize(self, masses, modes):
        L = modes.matrix.T
        freqs = modes.freqs
        freq_conv = np.sqrt(np.broadcast_to(freqs[:, np.newaxis], L.shape))
        mass_conv = np.sqrt(np.broadcast_to(self._tripmass(masses)[np.newaxis, :], L.shape))
        L = L * freq_conv * mass_conv
        Linv = (L / freq_conv**2)
        modes = type(modes)(self.molecule, L.T, inverse=Linv, freqs=freqs)
        return modes

    def _canonicalize_derivs(self, freqs, masses, derivs):

        if len(derivs) == 3:
            grad, fcs, fds = derivs
            fcs = fcs.array
            thirds = fds.third_deriv_array
            fourths = fds.fourth_deriv_array
        else:
            grad, fcs, thirds, fourths = derivs

        n = len(masses)
        modes_matrix = self.modes.inverse
        modes_n = len(modes_matrix)
        if modes_n == 3*n:
            modes_n = modes_n - 6
            modes_matrix = modes_matrix[6:]
            freqs = freqs[6:]
        coord_n = modes_n + 6
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

        undimension_2 = np.outer(m_conv, m_conv)
        fcs = fcs / undimension_2

        undimension_3 = np.outer(m_conv, m_conv)[np.newaxis, :, :] * f_conv[:, np.newaxis, np.newaxis]
        thirds = thirds * (1 / undimension_3 / np.sqrt(amu_conv))

        wat = np.outer(m_conv, m_conv)[np.newaxis, :, :] * (f_conv ** 2)[:, np.newaxis, np.newaxis]
        undimension_4 = SparseArray.from_diag(1 / wat / amu_conv)
        fourths = fourths
        fourths = fourths * undimension_4

        return grad, fcs, thirds, fourths

    def get_terms(self):
        grad = self.v_derivs[0]
        hess = self.v_derivs[1]
        thirds = self.v_derivs[2]
        fourths = self.v_derivs[3]

        # Use the Molecule's coordinates which know about their embedding by default
        intcds = self.internal_coordinates
        if intcds is None or not self.non_degenerate:
            # this is nice because it eliminates most of terms in the expansion
            xQ = self.modes.inverse
            xQQ = 0
            xQQQ = 0
            xQQQQ = 0
        else:
            YQ = self.modes.inverse
            QY = self.modes.matrix
            # We need to compute all these terms then mass weight them
            ccoords = self.coords
            carts = ccoords.system
            internals = intcds.system
            # XR, XRR, XRRR, XRRRR = [x.squeeze() for x in intcds.jacobian(carts, [1, 2, 3, 4])]
            # the 3rd and fourth derivative tensors will probably need to be optimized out

            # XR, XRR = [x.squeeze() for x in intcds.jacobian(carts, [1, 2])]
            # XRRR = XRRRR = 0

            XR, XRR, XRRR = [x.squeeze() for x in intcds.jacobian(carts, [1, 2, 3])]
            XRRRR = 0

            # The finite difference preserves too much shape by default
            _contract_dim = DumbTensor._contract_dim
            if XR.ndim > 2:
                XR = _contract_dim(XR, 2)
            if XRR.ndim > 3:
                XRR = _contract_dim(XRR, 3)
            if XRRR.ndim > 4:
                XRRR = _contract_dim(XRRR, 4)
            # if XRRRR.ndim > 5:
            #     XRRRR = _contract_dim(XRRRR, 5)

            RX = ccoords.jacobian(internals, 1)
            if RX.ndim > 2:
                RX = _contract_dim(RX, 2)

            # Need to then mass weight
            masses = self.masses
            mass_conv = np.sqrt(np.broadcast_to(masses[np.newaxis], (len(masses), 3)).flatten())
            YR = XR * mass_conv[np.newaxis]
            YRR = XRR * mass_conv[np.newaxis, np.newaxis]
            YRRR = XRRR * mass_conv[np.newaxis, np.newaxis, np.newaxis]
            YRRRR = XRRRR * mass_conv[np.newaxis, np.newaxis, np.newaxis, np.newaxis]
            RY = RX / mass_conv[:, np.newaxis]

            # transform into proper dXdQ space
            dot = DumbTensor._dot
            shift = DumbTensor._shift
            RQ = dot(YQ, RY)
            x_derivs = (YR, YRR, YRRR, YRRRR)
            Q_derivs = (RQ, 0, 0, 0)
            # xQ, xQQ, xQQQ, xQQQQ = self._get_tensor_derivs(Q_derivs, x_derivs, mixed_XQ=False)
            # xQ, xQQ = self._get_tensor_derivs(Q_derivs, x_derivs, mixed_XQ=False, order=2)
            xQ, xQQ, xQQQ = self._get_tensor_derivs(Q_derivs, x_derivs, mixed_XQ=False, order=3)
            # raise Exception(xQ-YQ)
            # xQ = xQ.transpose(0, 1)
            xQQ = xQQ.transpose(0, 2, 1)
            xQQQ = xQQQ.transpose(0, 2, 3, 1)
            # xQQQQ = xQQQQ.transpose(0, 1, 3, 4, 2)

            RQ = DumbTensor(RQ)
            YRR = DumbTensor(YRR)
            YRRR = DumbTensor(YRRR)
            YQQ = (RQ@(RQ@YRR)[0:1])[0:1]
            YQQQ = YRRR
            for i in range(3):
                YQQQ = RQ.dot(YQQQ, axes=[[1, i]])
            # plt.TensorPlot(YQQQ.t)
            # plt.TensorPlot(xQQQ)
            # plt.TensorPlot(xQQQ-YQQQ.t, plot_style=dict(vmin=-1e-8, vmax=1e-8)).show()

            # raise Exception(xQQ)
            # xQQ = 0
            # xQQQ = 0
            xQQQQ = 0

        x_derivs = (xQ, xQQ, xQQQ, xQQQQ)
        V_derivs = (grad, hess, thirds, fourths)

        v1, v2, v3, v4 = self._get_tensor_derivs(x_derivs, V_derivs, mixed_XQ=self.mixed_derivs)

        if self.mixed_derivs:
            for i in range(v4.shape[0]):
                v4[i, :, i, :] = v4[i, :, :, i] = v4[:, i, :, i] = v4[:, i, i, :] = v4[:, :, i, i] = v4[i, i, :, :]

        # test = UnitsData.convert("Hartrees", "Wavenumbers") * np.array([
        #     v4[2, 2, 2, 2],
        #     v4[1, 1, 2, 2],
        #     v4[1, 1, 1, 1],
        #     v4[0, 0, 2, 2],
        #     v4[0, 0, 1, 1],
        #     v4[0, 0, 0, 0]
        # ]).T

        # raise Exception([v2*UnitsData.convert("Hartrees", "Wavenumbers"),
        #                 v3*UnitsData.convert("Hartrees", "Wavenumbers")])

        v3 = UnitsData.convert("Wavenumbers", "Hartrees")*np.loadtxt("/Users/Mark/Desktop/test.dat")
        v3 = v3.reshape(3, 3, 3)

        return v2, v3, v4


class KineticTerms(ExpansionTerms):
    def __init__(self, molecule):
        """Represents the KE coefficients

        :param molecule: the molecule these modes are valid for
        :type molecule: Molecule
        :param internals: Optional internal coordinate set to rexpress in
        :type internals: CoordinateSystem | None
        """
        self.molecule = molecule
        self.masses = molecule.masses*UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        self.modes = self.undimensionalize(self.masses, molecule.normal_modes.basis)
        self.internal_coordinates = molecule.internal_coordinates
        self.coords = molecule.coords
        super().__init__()

    def undimensionalize(self, masses, modes):
        L = modes.matrix.T
        freqs = modes.freqs
        freq_conv = np.sqrt(np.broadcast_to(freqs[:, np.newaxis], L.shape))
        mass_conv = np.sqrt(np.broadcast_to(self._tripmass(masses)[np.newaxis, :], L.shape))
        L = L * freq_conv * mass_conv
        Linv = (L / freq_conv**2)
        modes = type(modes)(self.molecule, L.T, inverse=Linv, freqs=freqs)
        return modes

    def get_terms(self):

        dot = DumbTensor._dot
        shift = DumbTensor._shift
        # got_the_terms = self.masses is not None and \
        #                 len(self.masses) == 3 and \
        #                 not isinstance(self.masses[0], (int, np.integer, float, np.floating))
        #
        # if got_the_terms:
        #     G, GQ, GQQ = self.masses
        intcds = self.internal_coordinates
        if intcds is None:
            # this is nice because it eliminates a lot of terms in the expansion
            J = self.modes.matrix
            G = dot(J, J, axes=[[1, 1]])
            GQ = 0
            GQQ = 0
        else:
            ccoords = self.coords
            carts = ccoords.system
            internals = intcds.system

            # First we take derivatives of internals with respect to Cartesians
            RX, RXX, RXXX= ccoords.jacobian(internals, [1, 2, 3])
            # FD tracks too much shape

            _contract_dim = DumbTensor._contract_dim
            if RX.ndim > 2:
                RX = _contract_dim(RX, 2)
            if RXX.ndim > 3:
                RXX = _contract_dim(RXX, 3)
            if RXXX.ndim > 4:
                RXXX = _contract_dim(RXXX, 4)

            # Now we take derivatives of Cartesians with respect to internals
            XR, XRR = [x.squeeze() for x in intcds.jacobian(carts, [1, 2])]
            if XR.ndim > 2:
                XR = _contract_dim(XR, 2)
            if XRR.ndim > 3:
                XRR = _contract_dim(XRR, 3)

            # take only the well-defined coordinates

            sp = [x for x in np.arange(XR.shape[0]) if x not in (0, 1, 2, 4, 5, 8)]
            # XR = XR[sp, :]
            # plt.ArrayPlot(XR)
            # plt.ArrayPlot(XR@RX, plot_style=dict(vmin=-2, vmax=2))
            # plt.ArrayPlot(RX@XR, plot_style=dict(vmin=-2, vmax=2)).show()
            # XRR = XRR[sp, sp, :]
            # RX = RX[:, sp]
            # RXX = RXX[:, :, sp]
            # RXXX = RXXX[:, :, :, sp]
            # xr = xr[tuple(s-3 for s in sp), :]
            # rx = rx[:, tuple(s-3 for s in sp)]

            # next we need to mass-weight
            masses = self.masses
            mass_conv = np.sqrt(np.broadcast_to(masses[:, np.newaxis], (3, len(masses))).flatten())
            RY = RX / mass_conv[:, np.newaxis]
            RYY = RXX / (mass_conv[:, np.newaxis, np.newaxis] * mass_conv[np.newaxis, :, np.newaxis])
            RYYY = RXXX / (
                    mass_conv[:, np.newaxis, np.newaxis,   np.newaxis]
                    * mass_conv[np.newaxis, :, np.newaxis, np.newaxis]
                    * mass_conv[np.newaxis, np.newaxis, :, np.newaxis]
            )
            YR = XR * mass_conv[np.newaxis]
            YRR = XRR * mass_conv[np.newaxis, np.newaxis]

            QY = self.modes.matrix  # derivatives of Q with respect to the mass-weighted Cartesians
            YQ = self.modes.inverse
            G = dot(QY, QY, axes=[[0, 0]])

            J = DumbTensor(QY)
            Jd = DumbTensor(YQ)
            H = YQQ = 0
            K = DumbTensor(dot(RYY, YR, QY))
            U = K.dot(J, axes=[[0, 0]])
            L = DumbTensor(dot(RYYY, YR, QY))
            K22 = K.dot(K, axes=[[1, 1]])
            V = L[3:2]@J + K22[2:0]

            GQ = Jd@(U+U[2:1])
            GQ = GQ.t

            dRdYY = DumbTensor(RYY)
            dRdY = DumbTensor(RY)
            dYdR = DumbTensor(YR)
            U = dYdR@(dRdYY[2:1]@dRdY)
            U2 = U[2:1]
            intdGdR = U.t + U2.t

            dRdYYY = DumbTensor(RYYY)
            dYdRR = DumbTensor(YRR)
            dUdY = (
                    dRdY@dYdRR@dRdYY[2:1]@dRdY +
                    (dYdR@(dRdYYY[3:2]@dRdY)[0:1])[0:1] +
                    (dYdR@(dRdYY[2:1]@dRdYY[0:1])[2:1])[0:1]
            )
            intdGdRR = dYdR@(dUdY+dUdY[3:2])
            intdGdRR = intdGdRR.t

            QR = dot(YR, QY)
            RQ = dot(YQ, RY)
            intdGdQ = dot(RQ, intdGdR)
            GQ_2 = dot(intdGdQ, QR, QR, axes=[[1, 0], [1, 0]])
            GQ_2w = GQ_2 * UnitsData.convert("Hartrees", "Wavenumbers")

            intdGdQQ = dot(RQ, dot(RQ, intdGdRR), axes=[[1, 1]])
            GQQ_2 = dot(intdGdQQ, QR, QR, axes=[[2, 0], [2, 0]])

            CoordinateSet = crds.CoordinateSet
            CartesianCoordinates3D = crds.CartesianCoordinates3D
            ZMatrixCoordinateSystem = crds.ZMatrixCoordinateSystem
            icsys = ZMatrixCoordinateSystem(ordering=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
            ccs = CoordinateSet(ccoords, CartesianCoordinates3D)
            rx = _contract_dim(ccs.jacobian(icsys), 2)
            ics = ccs.convert(icsys)

            ic = np.asarray(intcds).flatten()[sp]
            r1 = ic[0]
            r2 = ic[1]
            q = ic[2]
            mO = masses[0]; mH = masses[1]; mD = masses[2]
            G_analytic = np.array([
                [(mH + mO) / (mH * mO), np.cos(q) / mO, -(np.sin(q) / (mO * r2))],
                [np.cos(q) / mO, (mD + mO) / (mD * mO), -(np.sin(q) / (mO * r1))],
                [-(np.sin(q) / (mO * r2)), -(np.sin(q) / (mO * r1)),
                 (mH + mO) / (mH * mO * r1 ** 2) + (mD + mO) / (mD * mO * r2 ** 2) - (2 * np.cos(q)) / (mO * r1 * r2)]
            ])
            GR_analytic = np.array([
                [
                    [0, 0, 0],
                    [0, 0, np.sin(q) / (mO * r1 ** 2)],
                    [0, np.sin(q) / (mO * r1 ** 2), (-2 * (mH + mO)) / (mH * mO * r1 ** 3) + (2 * np.cos(q)) / (mO * r1 ** 2 * r2)]
                ],
                [
                    [0, 0, np.sin(q) / (mO * r2 ** 2)],
                    [0, 0, 0],
                    [np.sin(q) / (mO * r2 ** 2), 0, (-2 * (mD + mO)) / (mD * mO * r2 ** 3) + (2 * np.cos(q)) / (mO * r1 * r2 ** 2)]
                ],
                [
                    [0, -(np.sin(q) / mO), -(np.cos(q) / (mO * r2))],
                    [-(np.sin(q) / mO), 0, -(np.cos(q) / (mO * r1))],
                    [-(np.cos(q) / (mO * r2)), -(np.cos(q) / (mO * r1)), (2 * np.sin(q)) / (mO * r1 * r2)]
                ]
            ])
            GRR_analytic = np.array([
                [
                    [
                        [0, 0, 0],
                        [0, 0, (-2*np.sin(q))/(mO*r1**3)],
                        [0, (-2*np.sin(q))/(mO*r1**3), (6*(mH + mO))/(mH*mO*r1**4) - (4*np.cos(q))/(mO*r1**3*r2)]
                    ],
                    [
                        [0, 0, 0],
                        [0, 0, 0],
                        [0, 0, (-2*np.cos(q))/(mO*r1**2*r2**2)]
                    ],
                    [
                        [0, 0, 0],
                        [0, 0, np.cos(q)/(mO*r1**2)],
                        [0, np.cos(q)/(mO*r1**2), (-2*np.sin(q))/(mO*r1**2*r2)]
                    ]
                ],
                [
                    [
                        [0, 0, 0],
                        [0, 0, 0],
                        [0, 0, (-2*np.cos(q))/(mO*r1**2*r2**2)]
                    ],
                    [
                        [0, 0, (-2*np.sin(q))/(mO*r2**3)],
                        [0, 0, 0],
                        [(-2*np.sin(q))/(mO*r2**3), 0, (6*(mD + mO))/(mD*mO*r2**4) - (4*np.cos(q))/(mO*r1*r2**3)]
                    ],
                    [
                        [0, 0, np.cos(q)/(mO*r2**2)],
                        [0, 0, 0],
                        [np.cos(q)/(mO*r2**2), 0, (-2*np.sin(q))/(mO*r1*r2**2)]
                    ]
                ],
                [
                    [
                        [0, 0, 0],
                        [0, 0, np.cos(q)/(mO*r1**2)],
                        [0, np.cos(q)/(mO*r1**2), (-2*np.sin(q))/(mO*r1**2*r2)]
                    ],
                    [
                        [0, 0, np.cos(q)/(mO*r2**2)],
                        [0, 0, 0],
                        [np.cos(q)/(mO*r2**2), 0, (-2*np.sin(q))/(mO*r1*r2**2)]
                    ],
                    [
                        [0, -(np.cos(q)/mO), np.sin(q)/(mO*r2)],
                        [-(np.cos(q)/mO), 0, np.sin(q)/(mO*r1)],
                        [np.sin(q)/(mO*r2), np.sin(q)/(mO*r1), (2*np.cos(q))/(mO*r1*r2)]
                    ]
                ]
            ])
            intdGdRR = intdGdRR[sp, :, :, :][:, sp, :, :][:, :, sp, :][:, :, :, sp]
            plt.TensorPlot(GRR_analytic, plot_style=dict(vmin=-1e-4, vmax=1e-4))
            plt.TensorPlot(intdGdRR, plot_style=dict(vmin=-1e-4, vmax=1e-4))
            plt.TensorPlot(GRR_analytic-intdGdRR, plot_style=dict(vmin=-1e-4, vmax=1e-4)).show()

            GQQ = (Jd@(Jd@(V+V[3:2]))[0:1]).t

            GQQ_w = GQQ * UnitsData.convert("Hartrees", "Wavenumbers")

            # plt.TensorPlot(GQQ)
            # plt.TensorPlot(GQQ_2)
            # plt.TensorPlot(GQQ-GQQ_2).show()

            # GQQ = 0#self._weight_derivatives(GQQ_2, 2)
            # GQQ = GQQ_2
            # GQQ = 0

            GQQ_2w = GQQ_2 * UnitsData.convert("Hartrees", "Wavenumbers")

        G_terms = (G, GQ, GQQ)
        return G_terms