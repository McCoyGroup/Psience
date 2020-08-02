"""
Stores all of the terms used inside the VPT2 representations
"""

import numpy as np, functools as fp, itertools as ip
from McUtils.Numputils import SparseArray
from McUtils.Data import UnitsData

from ..Molecools import Molecule, NormalModeCoordinates
from .Common import PerturbationTheoryException

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
                    raise Exception([a.shape, b.shape, kw])
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
        Linv = L / freq_conv / mass_conv
        modes_V = type(modes)(self.molecule, Linv, freqs=freqs)
        return modes_V

    def _canonicalize_derivs(self, freqs, masses, derivs):

        if len(derivs) == 3:
            grad, fcs, fds = derivs
            fcs = fcs.array
            thirds = fds.third_deriv_array
            fourths = fds.fourth_deriv_array
        else:
            grad, fcs, thirds, fourths = derivs

        modes_n = len(self.modes.matrix)
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
        fcs = fcs * undimension_2

        undimension_3 = np.outer(m_conv, m_conv)[np.newaxis, :, :] / f_conv[:, np.newaxis, np.newaxis]
        thirds = thirds * (undimension_3 / np.sqrt(amu_conv))

        wat = np.outer(m_conv, m_conv)[np.newaxis, :, :] / (f_conv ** 2)[:, np.newaxis, np.newaxis]
        undimension_4 = SparseArray.from_diag(wat / amu_conv)
        fourths = fourths
        fourths = fourths * undimension_4

        return grad, fcs, thirds, fourths

    def get_terms(self):
        # Use the Molecule's coordinates which know about their embedding by default
        intcds = self.internal_coordinates
        if intcds is None or not self.non_degenerate:
            # this is nice because it eliminates most of terms in the expansion
            xQ = self.modes.matrix
            xQQ = 0
            xQQQ = 0
            xQQQQ = 0
        else:
            # We need to compute all these terms then mass weight them
            ccoords = self.coords
            carts = ccoords.system
            internals = intcds.system
            xQ = intcds.jacobian(carts, 1)
            xQQ = intcds.jacobian(carts, 2)
            xQQQ = intcds.jacobian(carts, 3)
            xQQQQ = intcds.jacobian(carts, 4, stencil=6)

            # The finite difference preserves too much shape by default
            if xQ.ndim == 3:  # means we have NX3-type coords
                ncoords = np.product(xQ.shape[1:]).astype(int)
                xQ = xQ.reshape((xQ.shape[0], ncoords))
            if xQQ.ndim == 4:  # means we have NX3-type coords
                ncoords = np.product(xQQ.shape[2:]).astype(int)
                xQQ = xQQ.reshape(xQQ.shape[:2] + (ncoords,))
            if xQQQ.ndim == 5:  # means we have NX3-type coords
                ncoords = np.product(xQQQ.shape[3:]).astype(int)
                xQQQ = xQQQ.reshape(xQQQ.shape[:3] + (ncoords,))
            if xQQQQ.ndim == 6:  # means we have NX3-type coords
                ncoords = np.product(xQQQQ.shape[4:]).astype(int)
                xQQQQ = xQQQQ.reshape(xQQQQ.shape[:4] + (ncoords,))

            # Need to then mass weight
            masses = self.masses
            mass_conv = np.sqrt(np.broadcast_to(masses[np.newaxis], (len(masses), 3)).flatten())
            # xQ /= mass_conv[np.newaxis]
            # xQQ /= mass_conv[np.newaxis, np.newaxis]
            # xQQQ /= mass_conv[np.newaxis, np.newaxis, np.newaxis]
            # xQQQQ /= mass_conv[np.newaxis, np.newaxis, np.newaxis, np.newaxis]

            # Put everything into proper normal modes
            dot = self._dot
            RX = ccoords.jacobian(internals, 1)
            if RX.ndim == 3:
                ncrd = np.prod(RX.shape[1:]).astype(int)
                RX = RX.reshape((RX.shape[0], ncrd))
            YQ = self.modes.matrix
            L = dot(YQ, RX)
            xQ = dot(L, xQ, axes=[[1, 0]])
            for i in range(2):
                xQQ = dot(L, xQQ, axes=[[1, i]])
            for i in range(3):
                xQQQ = dot(L, xQQQ, axes=[[1, i]])
            for i in range(4):
                xQQQQ = dot(L, xQQQQ, axes=[[1, i]])

        x_derivs = (xQ, xQQ, xQQQ, xQQQQ)

        grad = self.v_derivs[0]
        hess = self.v_derivs[1]
        thirds = self.v_derivs[2]
        fourths = self.v_derivs[3]
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
        modes = type(modes)(self.molecule, L.T, freqs=freqs)
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
            # raise Exception([ccoords.convert(internals) - intcds])

            # First we take derivatives of internals with respect to Cartesians
            RX = ccoords.jacobian(internals, 1)
            RXX = ccoords.jacobian(internals, 2)
            RXXX = ccoords.jacobian(internals, 3)
            # FD tracks too much shape
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
            if RX.ndim > 2:
                RX = _contract_dim(RX, 2)
            if RXX.ndim > 3:
                RXX = _contract_dim(RXX, 3)
            if RXXX.ndim > 4:
                RXXX = _contract_dim(RXXX, 4)

            # Now we take derivatives of Cartesians with respect to internals
            XR = intcds.jacobian(carts, 1).squeeze()
            if XR.ndim > 2:
                XR = _contract_dim(XR, 2)

            XRR = intcds.jacobian(carts, 2)
            if XRR.ndim > 3:
                XRR = _contract_dim(XRR, 3)

            # non_garbage_coords = spec = (0, 3, 4)
            # RX   =   RX[:, spec]
            # RXX  =  RXX[:, :, spec]
            # RXXX = RXXX[:, :, :, spec]
            # XR   =  XR[spec, :]
            # XRR  =  XRR[spec, spec, :]

            # next we need to mass-weight
            masses = self.masses
            mass_conv = np.sqrt(np.broadcast_to(masses[:, np.newaxis], (3, len(masses))).flatten())
            RY = RX / mass_conv[:, np.newaxis]
            RYY = RXX / (mass_conv[:, np.newaxis]*mass_conv[:, np.newaxis])
            RYYY = RXXX / (
                    mass_conv[:, np.newaxis, np.newaxis, np.newaxis]
                    * mass_conv[np.newaxis, :, np.newaxis, np.newaxis]
                    * mass_conv[np.newaxis, np.newaxis, :, np.newaxis]
            )
            YR = XR * mass_conv[np.newaxis]
            YRR = XRR * mass_conv[np.newaxis, np.newaxis]

            QY = self.modes.matrix
            YQ = self.modes.matrix.T / self.modes.freqs[:, np.newaxis]
            G = dot(QY, QY, axes=[[0, 0]])

            QYY = dot(RYY, YR, QY)
            J = DumbTensor(QY)
            Jd = DumbTensor(YQ)
            H = YQQ = 0
            K = DumbTensor(QYY)
            U = K.dot(J, axes=[[1, 0]])
            L = DumbTensor(dot(RYYY, YR, QY))
            K22 = K.dot(K, axes=[[1, 1]])
            V = L[3:2]@J + K22[2:0]

            GQ = Jd@(U+U[2:1])
            GQ = GQ.t

            GQQ = Jd@(Jd@(V+V[3:2]))[0:1]
            GQQ = GQQ.t

            GQQ = self._weight_derivatives(GQQ, 2)

            # import McUtils.Plots as plt
            #
            # QR = dot(YR, QY)
            # RQ = dot(YQ, RY)

            # plot_style = dict(vmin=-5e-4, vmax=.5e-4, cmap='cividis')
            # toot = plt.ArrayPlot(dot(YR, RY))
            # plt.ArrayPlot(dot(RY, YR))
            # plt.ArrayPlot(dot(QY, YQ))
            # plt.ArrayPlot(dot(QR, RQ))
            # plt.ArrayPlot(dot(RQ, QR))
            # toot = plt.ArrayPlot(YQ, image_size=(900, 300), aspect_ratio=1/3)
            # plt.ArrayPlot(np.round(QR.T, 3), image_size=(900, 300), aspect_ratio=1 / 3)
            # plt.ArrayPlot(dot(YQ, QY))
            # toot.show()
            # tp = plt.TensorPlot(GQ)
            # wat = dot(RYY, QR)
            # woop = np.average(wat, axis=0)
            # plt.TensorPlot(wat.reshape(3, 3, 9, 9))
            # plt.ArrayPlot(woop)
            # diff_woop = np.round((wat - woop[np.newaxis, :, :])/np.max(wat), 2)
            # plt.TensorPlot(diff_woop.reshape(3, 3, 9, 9)).show()
            # plt.TensorPlot(QYY.reshape(3, 3, 9, 3)).show()
            # tp.colorbar = True
            # plt.ArrayPlot(GQ[0], colorbar={'tick_padding':50}, padding=((1, 80), (1, 1))).show()
            # meh = np.round(G, 4)
            # plt.ArrayPlot(meh).show()
            # toot.show()
            # raise Exception([np.average(G), np.average(GQ), np.average(GQQ)])

            # GQ = GQ * self.modes.freqs[:, np.newaxis, np.newaxis]

        G_terms = (G, GQ, GQQ)
        return G_terms