import numpy as np, scipy.sparse as sp, functools as fp, itertools as ip
from ..Wavefun import Wavefunctions, Wavefunction
from ..Molecools import NormalModeCoordinates
from McUtils.Numputils import SparseArray
from McUtils.Coordinerds import CoordinateSystem, CoordinateSet, CartesianCoordinates3D
from McUtils.Data import UnitsData

__all__ = [
    'PerturbationTheoryWavefunctions',
    'PerturbationTheoryException',
    'PerturbationTheoryHamiltonian'
]

class PerturbationTheoryException(Exception):
    pass

class PerturbationTheoryWavefunction(Wavefunction):
    @property
    def coeffs(self):
        return self.data

class PerturbationTheoryWavefunctions(Wavefunctions):
    def __init__(self, basis, energies, coeffs, hamiltonian):
        """Energies for the wavefunctions generated

        :param energies:
        :type energies:
        :param coeffs:
        :type coeffs:
        :param hamiltonian:
        :type hamiltonian:
        """
        self.hamiltonian = hamiltonian
        self.basis = basis
        super().__init__(
            energies=energies,
            wavefunctions=coeffs,
            wavefunction_class=PerturbationTheoryWavefunction
        )

class PerturbationTheoryHamiltonian:
    """Represents the main handler used in the perturbation theory calculation
    I probably want it to rely on a Molecule object...

    :param coords: The current coordinates of the system (used for computing the Jacobians)
    :type coords: np.ndarray | CoordinateSet
    :param masses: The masses of the coordinate system (used for computing the G-matrix and friends)
    :param pot_derivs: The Cartesian derivatives of the potential
    :type pot_derivs:
    :param modes: The Cartesian coordinate normal modes and the corresponding inverse
    :type modes: NormalModeCoordinates
    :param internals: The internal coordinate normal modes
    :type internals: CoordinateSystem
    :param numbers_of_quanta: The numbers of quanta of excitation to use for every mode
    :type numbers_of_quanta: np.ndarray | Iterable[int]
    """
    def __init__(self,
                 *ignore,
                 coords = None,
                 masses = None,
                 pot_derivs = None,
                 modes = None,
                 internals = None,
                 n_quanta = 3,
                 undimensionalize = True
                 ):
        if len(ignore) > 0:
            raise PerturbationTheoryException("{} takes no positional arguments".format(
                type(self).__name__
            ))
        mode_n = modes.matrix.shape[0]
        self.mode_n = mode_n
        self.n_quanta = np.full((mode_n,), n_quanta) if isinstance(n_quanta, (int, np.int)) else tuple(n_quanta)
        # we assume that our modes are already in mass-weighted coordinates, so now we need to make them dimensionless
        # to make life easier, we include this undimensionalization differently for the kinetic and potential energy
        # the kinetic energy needs to be weighted by a sqrt(omega) term while the PE needs a 1/sqrt(omega)
        self.modes = modes
        if isinstance(coords, CoordinateSet) and internals is None:
            internals = coords.system
        if undimensionalize:
            modes_T, modes_V = self.undimensionalize(masses, modes)
        else:
            modes_T, modes_V = modes, modes
        self.internals = internals
        self.V_terms = self.PotentialTerms(pot_derivs, coords, masses, modes_V, internals)
        self.G_terms = self.KineticTerms(coords, masses, modes_T, internals)
    def undimensionalize(self, masses, modes):
        L = modes.matrix
        freqs = modes.freqs
        freq_conv = np.sqrt(np.broadcast_to(freqs[:, np.newaxis], L.shape))
        mass_conv = np.sqrt(np.broadcast_to(self._tripmass(masses)[np.newaxis, :], L.shape))
        L, Linv = L * freq_conv * mass_conv, L / freq_conv / mass_conv
        modes_T = NormalModeCoordinates(L, freqs=freqs)
        modes_V = NormalModeCoordinates(Linv, freqs=freqs)
        return modes_T, modes_V
    @staticmethod
    def _tripmass(masses):
        return np.broadcast_to(masses, (len(masses), 3)).T.flatten()
    @classmethod
    def from_fchk(cls, file, internals = None, n_quanta = 3):
        from McUtils.GaussianInterface import GaussianFChkReader

        with GaussianFChkReader(file) as gr:
            parse = gr.parse(["Coordinates", "Gradient", "AtomicMasses",
                              "ForceConstants", "ForceDerivatives", "VibrationalModes", "VibrationalData"])

        amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        coords = UnitsData.convert("Angstroms", "AtomicUnitOfLength") * parse["Coordinates"]

        masses = amu_conv * parse["AtomicMasses"]
        modes = parse["VibrationalModes"]
        freqs = parse["VibrationalData"]["Frequencies"] * UnitsData.convert("Wavenumbers", "Hartrees")
        fcs = parse["ForceConstants"].array

        raw = np.sqrt(np.diag(np.dot(np.dot(modes, fcs), modes.T)))
        reweight = freqs / raw
        modes = modes * reweight[:, np.newaxis]
        # new = np.sqrt(np.diag(np.dot(np.dot(modes, fcs), modes.T)))
        # raise Exception(UnitsData.convert("Hartrees", "Wavenumbers")*new)

        fds = parse["ForceDerivatives"]

        m_conv = np.sqrt(cls._tripmass(masses))
        f_conv = np.sqrt(freqs)

        undimension_2 = np.outer(m_conv, m_conv)
        fcs = fcs * undimension_2

        undimension_3 = np.outer(m_conv, m_conv)[np.newaxis, :, :] / f_conv[:, np.newaxis, np.newaxis]
        thirds = fds.third_deriv_array * (undimension_3 / np.sqrt(amu_conv))

        wat = np.outer(m_conv, m_conv)[np.newaxis, :, :] / (f_conv**2)[:, np.newaxis, np.newaxis]
        undimension_4 = SparseArray.from_diag(wat/ amu_conv)
        fourths = fds.fourth_deriv_array
        fourths = fourths * undimension_4
        #symm4=np.tensordot(fourths.tensordot(modes, axes=[2, 1]), modes, axes=[2, 1])

        # test4 = np.array([
        #     [symm4[1, 2, 0], symm4[1, 0, 2]],
        #     [symm4[0, 2, 1], symm4[2, 0, 1]]
        # ])

        return cls(
            coords=coords,
            masses=masses,
            pot_derivs=[parse["Gradient"], fcs, thirds, fourths],
            modes=NormalModeCoordinates(modes, freqs=freqs),
            internals=internals,
            n_quanta=n_quanta # need a way to pull the internals from the FChk or provide an alternate way to get them...
        )

    class ExpansionTerms:
        def __init__(self):
            self._terms = None
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
        def _dot(*t, axes=None):
            """
            Flexible tensordot
            """

            if len(t) == 1:
                return t[0]

            if any(isinstance(x, int) for x in t):
                return 0

            tdot = lambda a, b, **kw: (a.tensordot(b, **kw) if hasattr(a, "tensordot") else np.tensordot(a, b, **kw))
            td = lambda a, b: (0 if isinstance(a, int) or isinstance(b[0], int) else tdot(a, b[0], axes=b[1]))
            if axes is None:
                axes = [1] * (len(t) - 1)

            return fp.reduce(td, zip(t[1:], axes), t[0])

        @staticmethod
        def _weight_derivatives(t):
            if isinstance(t, int):
                return t
            weighted = t
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

            dot = cls._dot
            shift = cls._shift

            derivs = [None]*order


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
                    dot(xQ_22_Vxxx, xQ, xQQ, axes=[[0, -1], [1, -1]]) +
                    shift(xQQQ_Vxx_xQ, (0, 1), (2, 3))
            )

            V_QQQQ_3 = (
                    xQQQ_Vxx_xQ +
                    dot(xQ_22_Vxxx, xQQ, xQ, axes=[[0, 2], [1, 1]]) +
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
                    shift(dot(xQQ, VQxx, xQ, axes=[[2, 2], [2, 1]]), (3, 1)) +
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
        def __init__(self, v_derivs, coords, masses, modes, internals, mixed_derivs = True):
            self.coords = coords
            self.modes = modes
            self.masses = masses
            self.internals = internals
            self.v_derivs = self._canonicalize_derivs(v_derivs)
            self.mixed_derivs = mixed_derivs
            super().__init__()

        def _canonicalize_derivs(self, derivs):

            grad, fcs, thirds, fourths = derivs

            modes_n = len(self.modes.matrix)
            coord_n = modes_n+6
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

            return grad, fcs, thirds, fourths

        def get_terms(self):
            if self.internals is None:
                # this is nice because it eliminates most of terms in the expansion
                xQ = self.modes.matrix
                xQQ = 0
                xQQQ = 0
                xQQQQ = 0
            else:
                # I think this could do with some sprucing up...
                ccoords = CoordinateSet(self.coords, CartesianCoordinates3D)
                xQ = self.internals.jacobian(ccoords, CartesianCoordinates3D, self.modes, 1)
                xQQ = self.internals.jacobian(ccoords, CartesianCoordinates3D, self.modes, 2)
                xQQQ = self.internals.jacobian(ccoords, CartesianCoordinates3D, self.modes, 3)
                xQQQQ = self.internals.jacobian(ccoords, CartesianCoordinates3D, self.modes, 4)

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
        def __init__(self, coords, masses, modes, internals = None):
            """Represents the KE coefficients

            :param masses: masses of the atoms in the modes
            :type masses: np.ndarray | None
            :param modes: Cartesian displacement normal modes
            :type modes: NormalModeCoordinates
            :param internals: Optional internal coordinate set to rexpress in
            :type internals: CoordinateSystem | None
            """
            self.coords = coords
            self.masses = masses
            self.modes = modes
            self.internals = internals
            super().__init__()

        def get_terms(self):

            dot = self._dot
            shift = self._shift
            got_the_terms = self.masses is not None and \
                            len(self.masses) == 3 and \
                            not isinstance(self.masses[0], (int, np.integer, float, np.floating))

            if got_the_terms:
                G, GQ, GQQ = self.masses
            elif self.internals is None:
                # this is nice because it eliminates a lot of terms in the expansion
                J = self.modes.matrix
                G = dot(J, J, axes=[[1, 1]])
                GQ = 0
                GQQ = 0
            else:
                # this is like _probably_ correct
                L = self.modes.matrix

                # we're dropping this because we want to go dimensionless
                J = L#dot(L, m)

                carts = CartesianCoordinates3D
                ccoords = CoordinateSet(self.coords, CartesianCoordinates3D)
                XR = ccoords.jacobian(self.internals, 1)
                XRR = ccoords.jacobian(self.internals, 2)
                intcds = ccoords.convert(self.internals)
                RXX = intcds.jacobian(carts, 2)
                RXXX = intcds.jacobian(carts, 3)
                RQ = np.linalg.pinv(J)

                YQ = dot(RQ, XR)
                YQQ = dot(RQ, shift(dot(RQ, XRR, axes=[[1, 1]]), (1, 0)))
                JY = dot(L, RXX)
                JYY = dot(L, RXXX)

                JQ = dot(YQ, JY)
                JQQ = dot(YQQ, JY) + dot(YQ, shift(dot(YQ, JYY, axes=[[1, 1]]), (1, 0)))

                G = dot(J, shift(J, (0, 1)))
                GQ_1 = dot(JQ, J, axes=[[2, 1]])
                GQ = GQ_1 + shift(GQ_1, (2, 3))
                GQQ_1 = dot(JQQ, J, axes=[[3, 1]])
                GQQ_2 = dot(JQ, JQ, axes=[[2, 2]])
                GQQ = GQQ_1 + shift(GQQ_2, (2, 0)) + shift(GQQ_2, (2, 1)) + shift(GQQ_1, (2, 3))

            G_terms = (G, GQ, GQQ)
            return G_terms

    class SubHamiltonian:
        def __init__(self, compute, n_quanta):
            self.compute = compute
            self.dims = n_quanta
        @property
        def diag(self):
            ndims = int(np.prod(self.dims))
            return self[np.arange(ndims), np.arange(ndims)]
        def get_element(self, n, m):
            """Pulls elements of the Hamiltonian, first figures out if it's supposed to be pulling blocks...

            :param n:
            :type n:
            :param m:
            :type m:
            :return:
            :rtype:
            """

            dims = self.dims
            ndims = int(np.prod(dims))
            idx = (n, m)

            pull_elements = True
            if isinstance(n, np.ndarray) and isinstance(m, np.ndarray):
                if len(n.shape) > 1 and len(m.shape) > 1:
                    pull_elements = False
            if pull_elements:
                pull_elements = all(isinstance(x, (int, np.integer)) for x in idx)
                if not pull_elements:
                    pull_elements = all(not isinstance(x, (int, np.integer, slice)) for x in idx)
                    if pull_elements:
                        e1 = len(idx[0])
                        pull_elements = all(len(x) == e1 for x in idx)
            if not isinstance(n, int):
                if isinstance(n, np.ndarray):
                    n = n.flatten()
                if not isinstance(n, slice):
                    n = np.array(n)
                n = np.arange(ndims)[n]
            else:
                n = [n]
            if not isinstance(m, int):
                if isinstance(m, np.ndarray):
                    m = m.flatten()
                if not isinstance(m, slice):
                    m = np.array(m)
                m = np.arange(ndims)[m]
            else:
                m = [m]
            if pull_elements:
                n = np.unravel_index(n, dims)
                m = np.unravel_index(m, dims)
            else:
                blocks = np.array(list(ip.product(n, m)))
                n = np.unravel_index(blocks[:, 0], dims)
                m = np.unravel_index(blocks[:, 1], dims)
            def pad_lens(a, b):
                if isinstance(a, (int, np.integer)) and not isinstance(b, (int, np.integer)):
                    a = np.full((len(b),), a)
                if isinstance(b, (int, np.integer)) and not isinstance(a, (int, np.integer)):
                    b = np.full((len(a),), b)
                return a,b
            i = tuple(pad_lens(a, b) for a,b in zip(n,m))
            els = self.compute(i)
            if pull_elements:
                return els
            else:
                shp = (len(np.unique(blocks[:, 0])), len(np.unique(blocks[:, 1])))
                return els.reshape(shp).squeeze()

        def __getitem__(self, item):
            if not isinstance(item, tuple):
                item = (item,)
            if len(item) == 1:
                item = item + (slice(None, None, None),)
            return self.get_element(*item)

    @property
    def H0(self):
        def compute_H1(inds,
                       G=self.G_terms[0],
                       V=self.V_terms[0],
                       pp=self.ProductOperator.pp(self.n_quanta),
                       QQ=self.ProductOperator.QQ(self.n_quanta),
                       H=self._compute_h0
                       ):
            return H(inds, G, V, pp, QQ)

        return self.SubHamiltonian(compute_H1, self.n_quanta)

    def _compute_h0(self, inds, G, F, pp, QQ):
        """

        :param inds: which elements of H0 to compute
        :param G: The G-matrix
        :type G: np.ndarray
        :param F: The Hessian
        :type F: np.ndarray
        :param pp: Matrix representation of pp
        :type pp:
        :param QQ: Matrix representation of QQ
        :type QQ:
        :return:
        :rtype:
        """

        # print(type(gmatrix_derivs))
        if not isinstance(G, int):
            # takes a 5-dimensional SparseTensor and turns it into a contracted 2D one
            subKE = pp[inds]
            if isinstance(subKE, np.ndarray):
                ke = np.tensordot(subKE.squeeze(), -G, axes=[[0, 1], [0, 1]])
            else:
                ke = subKE.tensordot(-G, axes=[[0, 1], [0, 1]]).squeeze()
        else:
            ke = 0

        if not isinstance(F, int):
            subPE = QQ[inds]
            if isinstance(subPE, np.ndarray):
                pe = np.tensordot(subPE.squeeze(), F, axes=[[0, 1], [0, 1]])
            else:
                pe = subPE.tensordot(F, axes=[[0, 1], [0, 1]]).squeeze()

        else:
            pe = 0

        return 1/2*ke + 1/2*pe

    @property
    def H1(self):
        def compute_H1(inds,
                       G=self.G_terms[1],
                       V=self.V_terms[1],
                       pQp=self.ProductOperator.pQp(self.n_quanta),
                       QQQ=self.ProductOperator.QQQ(self.n_quanta),
                       H=self._compute_h1
                       ):
            return H(inds, G, V, pQp, QQQ)
        return self.SubHamiltonian(compute_H1, self.n_quanta)

    def _compute_h1(self, inds, gmatrix_derivs, V_derivs, pQp, QQQ):
        """

        :param inds: which elements of H1 to compute
        :param gmatrix_derivs: The derivatives of the G-matrix with respect to Q
        :type gmatrix_derivs: np.ndarray
        :param V_derivs: The derivatives of V with respect to QQQ
        :type V_derivs: np.ndarray
        :param pQp: Matrix representation of pQp
        :type pQp:
        :param QQQ: Matrix representation of QQQ
        :type QQQ:
        :return:
        :rtype:
        """

        if not isinstance(gmatrix_derivs, int):
            subpQp = pQp[inds]
            if isinstance(subpQp, np.ndarray):
                subpQp = subpQp.squeeze()
                ke = np.tensordot(subpQp, -gmatrix_derivs, axes=[[0, 1, 2], [0, 1, 2]])
            else:
                ke = subpQp.tensordot(-gmatrix_derivs, axes=[[0, 1, 2], [0, 1, 2]]).squeeze()
        else:
            ke = 0

        if not isinstance(V_derivs, int):
            subQQQ = QQQ[inds]
            if isinstance(subQQQ, np.ndarray):
                subQQQ = subQQQ.squeeze()
                pe = np.tensordot(subQQQ, V_derivs, axes=[[0, 1, 2], [0, 1, 2]])
            else:
                pe = subQQQ.tensordot(V_derivs, axes=[[0, 1, 2], [0, 1, 2]]).squeeze()

        else:
            pe = 0

        # test = UnitsData.convert("Hartrees", "Wavenumbers")/6*subQQQ*V_derivs

        return 1/2*ke + 1/6 * pe

    @property
    def H2(self):
        def compute_H2(inds,
                       G=self.G_terms[2],
                       V=self.V_terms[2],
                       KE=self.ProductOperator.pQQp(self.n_quanta),
                       PE=self.ProductOperator.QQQQ(self.n_quanta),
                       H=self._compute_h2
                       ):
            return H(inds, G, V, KE, PE)

        return self.SubHamiltonian(compute_H2, self.n_quanta)

    def _compute_h2(self, inds, gmatrix_derivs, V_derivs, KE, PE):
        """

        :param inds: which elements of H1 to compute
        :param gmatrix_derivs: The derivatives of the G-matrix with respect to QQ
        :type gmatrix_derivs: np.ndarray
        :param V_derivs: The derivatives of V with respect to QQQQ
        :type V_derivs: np.ndarray
        :param KE: Matrix representation of pQQp
        :type KE:
        :param PE: Matrix representation of QQQQ
        :type PE:
        :return:
        :rtype:
        """

        # print(type(gmatrix_derivs))
        if not isinstance(gmatrix_derivs, int):
            keTens = KE[inds]
            if isinstance(keTens, np.ndarray):
                ke = np.tensordot(keTens.squeeze(), -gmatrix_derivs, axes=[[0, 1, 2, 3], [0, 1, 2, 3]])
            else:
                ke = keTens.tensordot(-gmatrix_derivs, axes=[[0, 1, 2, 3], [0, 1, 2, 3]]).squeeze()
        else:
            ke = 0

        # print(inds)

        if not isinstance(V_derivs, int):
            peTens = PE[inds]
            if isinstance(peTens, np.ndarray):
                pe = np.tensordot(peTens.squeeze(), V_derivs, axes=[[0, 1, 2, 3], [0, 1, 2, 3]])
            else:
                pe = peTens.tensordot(V_derivs, axes=[[0, 1, 2, 3], [0, 1, 2, 3]]).squeeze()
        else:
            pe = 0

        return 1/4*ke + 1/24*pe

    class ProductOperator:
        """
        Provides a (usually) _lazy_ representation of an operator, which allows things like
        QQQ and pQp to be calculated block-by-block
        """
        def __init__(self, funcs, quanta):
            """

            :param funcs:
            :type funcs:
            :param quanta:
            :type quanta:
            """
            self.funcs = funcs
            self.quanta = tuple(quanta)
            self.mode_n = len(quanta)
            self._tensor = None

        @property
        def ndim(self):
            return len(self.funcs) + len(self.quanta)
        @property
        def shape(self):
            return (self.mode_n,)*len(self.funcs) + self.quanta
        @property
        def tensor(self):
            if self._tensor is None:
                self._tensor = self.product_operator_tensor()
            return self._tensor

        def get_inner_indices(self):
            """Gets the n-dimensional array of ijkl (e.g.) indices that functions will map over

            :return:
            :rtype:
            """
            funcs = self.funcs
            dims = len(self.funcs)
            shp = (self.mode_n,) * dims
            inds = np.indices(shp, dtype=int)
            tp = np.roll(np.arange(len(funcs) + 1), -1)
            base_tensor = np.transpose(inds, tp)
            return base_tensor

        def __getitem__(self, item):
            return self.get_elements(item)
        def get_individual_elements(self, idx):
            if len(idx) != len(self.quanta):
                raise ValueError("number of indices requested must be the same as the number of modes")
            inds = self.get_inner_indices()
            idx = tuple(tuple(np.array([i]) if isinstance(i, (int, np.integer)) else i for i in j) for j in idx)
            funcs = self.funcs
            quants = self.quanta
            def pull(inds, f=funcs, x=idx, qn = quants):
                uinds = np.unique(inds)
                mats = self._operator_submatrix(f, qn, inds, return_kron=False)
                els = [m[x[i]] for m,i in zip(mats, uinds)]
                if isinstance(els[0], np.matrix):
                    els = [np.asarray(e).squeeze() for e in els]
                res = np.prod(els, axis=0)
                # print(inds, [len(x[i][0]) for i in inds])
                return res
            res = np.apply_along_axis(pull, -1, inds)
            return res
        def get_elements(self, idx):
            if len(idx) != len(self.quanta):
                raise ValueError("number of indices requested must be the same as the number of quanta")
            inds = self.get_inner_indices()
            idx = tuple(tuple(np.array([i]) if isinstance(i, (int, np.integer)) else i for i in j) for j in idx)
            tens = self.tensor
            quants = self.quanta
            def pull(inds, t=tens, x=idx, qn = quants):
                sly = t[tuple(inds)]
                uinds = np.unique(inds)
                sub = tuple(tuple(j) for i in uinds for j in x[i])
                try:
                    res = sly[sub]
                except:
                    print(sub)
                    raise

                missing = [i for i in range(len(x)) if i not in inds]
                equivs = [x[i][0] == x[i][1] for i in missing]

                orthog = np.prod(equivs, axis=0).astype(int)
                return res * orthog
            res = np.apply_along_axis(pull, -1, inds)
            return SparseArray(res.squeeze())

        def product_operator_tensor(self):
            """Generates the tensor created from the product of funcs over the dimensions dims, except for the fact that it
            makes a _ragged_ tensor in the final dimensions

            :param funcs:
            :type funcs:
            :param dims:
            :type dims:
            :return:
            :rtype:
            """

            dims = self.quanta
            funcs = self.funcs
            base_tensor = self.get_inner_indices()
            news_boy = lambda inds, f=funcs, d=dims: self._operator_submatrix(f, d, inds)
            news_boys = np.apply_along_axis(news_boy, -1, base_tensor)

            return news_boys

        def _operator_submatrix(self, funcs, dims, inds, padding = 3, return_kron = True):
            """Returns the operator submatrix for a product operator like piQjpk or whatever

            :param funcs: the functions that take a dimension size and return a matrix for the suboperator
            :type funcs:
            :param dims: dimensions of each coordinate (e.g. (5, 8, 2, 9))
            :type dims: tuple | np.ndarray
            :param inds: the list of indices
            :type inds: tuple | np.ndarray
            :param padding: the representation can be bad if too few terms are used so we add a padding
            :type padding: int
            :return:
            :rtype:
            """

            uinds = np.unique(inds)
            mm = {k:i for i,k in enumerate(uinds)}
            ndim = len(uinds)
            pieces = [None] * ndim
            for f, i in zip(funcs, inds):
                n = mm[i]
                if pieces[n] is None:
                    pieces[n] = f(dims[i]+padding)
                else:
                    pieces[n] = pieces[n].dot(f(dims[i]+padding))

            # for j in np.setdiff1d(totinds, inds):
            #     pieces[j] = sp.identity(dims[j])

            if return_kron:
                mat = sp.csr_matrix(fp.reduce(sp.kron, pieces))
                sub_shape = tuple(dims[i]+padding for i in np.unique(inds) for j in range(2))
                trans = tuple(j for i in zip(range(ndim), range(ndim, 2*ndim)) for j in i)
                mat = SparseArray(mat, shape=sub_shape).transpose(trans)
            else:
                mat = pieces
            return mat

        @staticmethod
        def pmatrix_ho(n):
            """

            :param n:
            :type n:
            :return:
            :rtype: sp.csr_matrix
            """
            # the imaginary terms pull out and just become a negative sign
            ar = 1/np.sqrt(2)*np.sqrt(np.arange(1, n))
            bands = [
                [ ar,  1],
                [-ar, -1]
            ]
            return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

        @staticmethod
        def qmatrix_ho(n):
            """

            :param n:
            :type n:
            :return:
            :rtype: sp.csr_matrix
            """

            ar = 1/np.sqrt(2)*np.sqrt(np.arange(1, n))
            bands = [
                [ar, 1],
                [ar, -1]
            ]
            return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

        @classmethod
        def QQ(cls, n_quanta, qmatrix=None):
            if qmatrix is None:
                qmatrix = cls.qmatrix_ho
            return cls((qmatrix, qmatrix), n_quanta)

        @classmethod
        def pp(cls, n_quanta, pmatrix=None):
            if pmatrix is None:
                pmatrix = cls.pmatrix_ho
            return cls((pmatrix, pmatrix), n_quanta)

        @classmethod
        def QQQ(cls, n_quanta, qmatrix=None):
            if qmatrix is None:
                qmatrix = cls.qmatrix_ho
            return cls((qmatrix, qmatrix, qmatrix), n_quanta)

        @classmethod
        def pQp(cls, n_quanta, pmatrix=None, qmatrix=None):
            if pmatrix is None:
                pmatrix = cls.pmatrix_ho
            if qmatrix is None:
                qmatrix = cls.qmatrix_ho
            return cls((pmatrix, qmatrix, pmatrix), n_quanta)

        @classmethod
        def QQQQ(cls, n_quanta, qmatrix=None):
            if qmatrix is None:
                qmatrix = cls.qmatrix_ho
            return cls((qmatrix, qmatrix, qmatrix, qmatrix), n_quanta)

        @classmethod
        def pQQp(cls, n_quanta, pmatrix=None, qmatrix=None):
            if pmatrix is None:
                pmatrix = cls.pmatrix_ho
            if qmatrix is None:
                qmatrix = cls.qmatrix_ho
            return cls((pmatrix, qmatrix, qmatrix, pmatrix), n_quanta)

    def get_state_indices(self, states):
        if isinstance(states, (int, np.integer)):
            states = np.arange(min([np.prod(self.n_quanta), states]))
        if not isinstance(states, slice):
            if not isinstance(states[0], (int, np.integer)):
                states = np.ravel_multi_index(np.array(states).T, self.n_quanta)
            if isinstance(states, tuple):  # numpy is weird
                states = np.array(states)
        return states

    def get_state_quantum_numbers(self, states):
        if isinstance(states, slice):
            states = np.arange(np.prod(self.n_quanta))[states]
        elif isinstance(states, int):
            states = np.arange(min([np.prod(self.n_quanta), states]))
        qns = tuple(np.array(np.unravel_index(states, self.n_quanta)).T)
        return qns

    def _get_corrections(self, states=15, coupled_states=None, coeff_threshold=None, energy_threshold=None):
        if states is None:
            states = np.prod(self.n_quanta)
        states = self.get_state_indices(states)
        if coupled_states is None:
            coupled_states = slice(None, None, None)
        if isinstance(coupled_states, slice):
            coupled_states = np.arange(np.prod(self.n_quanta))
        coupled_states = self.get_state_indices(coupled_states)

        H0 = self.H0
        H1 = self.H1
        H2 = self.H2

        energies = H0.diag
        state_E = energies[states]
        energies = energies[coupled_states]
        if isinstance(coupled_states, slice):
            H1_blocks = H1[states, coupled_states]
        else:
            ixes = np.ix_(states, coupled_states)
            H1_blocks = H1[ixes]

        e_blocks = state_E[:, np.newaxis] - np.broadcast_to(energies[np.newaxis], (len(states), len(energies)))

        for n,s in enumerate(states):
            w = np.where(coupled_states==s)[0]
            if len(w) > 0:
                e_blocks[n, w[0]] = 1 # gotta prevent blowups

        if energy_threshold is not None:
            if isinstance(energy_threshold, (int, float, np.integer, np.floating)):
                energy_threshold = (energy_threshold, 1)
            dropped = np.abs(e_blocks) < energy_threshold[0]
            e_blocks[dropped] = np.sign(e_blocks[dropped]) * energy_threshold[1]

        coeffs = H1_blocks / e_blocks

        if coeff_threshold is not None:
            if isinstance(coeff_threshold, (int, float, np.integer, np.floating)):
                coeff_threshold = (coeff_threshold, 0)
            dropped = np.abs(coeffs) > coeff_threshold[0]
            coeffs[dropped] = np.sign(coeffs[dropped]) * coeff_threshold[1]

        for n,s in enumerate(states):
            w = np.where(coupled_states==s)[0]
            if len(w) > 0:
                coeffs[n, w[0]] = 0 # we don't want to double count this

        e_co = np.expand_dims(coeffs, axis=1)
        e_H1 = np.expand_dims(H1_blocks, axis=2)
        e1 = np.matmul(e_co, e_H1).squeeze()
        e2s = H2[states, states]

        # print(self._fmt_corr2_matrix(states, coupled_states, H1_blocks, coeffs))

        return coeffs, (state_E , e1,  e2s), e_blocks, H1_blocks, states, coupled_states

    def _fmt_corr_matrix(self, states, coupled_states, H1_blocks, coeffs):
        """A useful debug function"""
        h1_corr = UnitsData.convert("Hartrees", "Wavenumbers") * H1_blocks * coeffs
        return "\n".join(
            [
                "Corr:",
                "    " + (" ".join(" {2}{1}{0}".format(*x) for x in self.get_state_quantum_numbers(coupled_states))),
                *("{2}{1}{0} ".format(*s) + (" ".join("{:<+4.0f}".format(x) for x in h)) for s, h in
                  zip(self.get_state_quantum_numbers(states), h1_corr))
                ]
        )

    def _fmt_corr2_matrix(self, states, coupled_states, H1_blocks, coeffs):
        """A useful debug function"""
        h1_corr = UnitsData.convert("Hartrees", "Wavenumbers") * H1_blocks
        return "\n".join(
            [
                "Corr:",
                "    " + (" ".join(" {2}{1}{0}".format(*x) for x in self.get_state_quantum_numbers(coupled_states))),
                *("{2}{1}{0} ".format(*s) + (" ".join("{:<+4.0f}".format(x) for x in h)) for s, h in
                  zip(self.get_state_quantum_numbers(states), h1_corr))
                ]
        )

    def get_corrections(self, states=15, coupled_states=None, coeff_threshold=None, energy_threshold=None):
        """

        :param states:
        :type states:
        :param coeff_threshold: a hack for ditching near degeneracies
        :type coeff_threshold: float | Iterable[float]
        :return:
        :rtype:
        """

        return self._get_corrections(states=states, coupled_states=coupled_states,
                                     coeff_threshold=coeff_threshold, energy_threshold=energy_threshold)[:2]

    def get_wavefunctions(self, states=15, coupled_states=None, coeff_threshold=None, energy_threshold=None):
            """Computes perturbation expansion of the wavefunctions and energies

            :param states: the states to target
            :type states: int | iterable[int] | None
            :return: coeffs
            :rtype: np.ndarray
            """

            coeffs, corrs, e_blocks, H1_blocks, states, coupled_states = self._get_corrections(
                states, coeff_threshold=coeff_threshold, coupled_states=coupled_states
            )
            energies = sum(corrs)

            # basis = np.array(np.unravel_index(np.arange(num_ens), self.n_quanta)).T

            return PerturbationTheoryWavefunctions(
                None, # need to get the basis back int, but using the coupled_states args
                # (UnitsData.convert("AtomicUnitOfEnergy", "Wavenumbers")*energies, basis),
                energies,
                coeffs,
                self
            )

    def martin_test(self, states=15, coupled_states=None):
        """Applies the Martin Test to all of the specified states and returns the resulting correlation matrix

        :param states:
        :type states:
        :return:
        :rtype:
        """
        states = self.get_state_indices(states)
        if coupled_states is None:
            coupled_states = states#slice(None, None, None)
        else:
            coupled_states = self.get_state_indices(states)
        if isinstance(coupled_states, slice):
            H1_blocks = self.H1[states, coupled_states]
        else:
            H1_blocks = self.H1[np.ix_(states, coupled_states)]

        energies = self.H0.diag
        state_energies = energies[states]
        coupled_energies = energies[coupled_states]
        diffs = state_energies[:, np.newaxis] - coupled_energies[np.newaxis, :]
        for n,s in enumerate(states):
            w = np.where(coupled_states==s)[0]
            if len(w) > 0:
                diffs[n, w[0]] = 1

        return (H1_blocks**4)/(diffs**3)



