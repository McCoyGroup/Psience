
__all__ = [
    "MolecularHamiltonian"
]

import numpy as np
import McUtils.Numputils as nput
from McUtils.Combinatorics import UniquePermutations
from McUtils.Zachary import TensorDerivativeConverter

from .CoordinateSystems import MolecularEmbedding
from .Properties import PotentialSurfaceManager, DipoleSurfaceManager, NormalModesManager

class MolecularHamiltonian:
    """
    Provides a generic molecular Hamiltonian with support for expansions (for VPT),
    construction of molecular models, and other wave function generators
    """

    def __init__(self,
                 embedding:MolecularEmbedding,
                 potential_manager:PotentialSurfaceManager=None,
                 modes_manager:NormalModesManager=None,
                 dipole_manager:DipoleSurfaceManager=None
                 ):
        self.embedding = embedding
        self.modes = modes_manager
        self.potential = ScalarOperatorManager(potential_manager)
        self.dipole = ScalarOperatorManager(dipole_manager)

    def get_VPT_expansions(self,
                           order=2,
                           coordinate_transformation=None
                           ):
        carts = self.embedding.internals is None
        embedding = ModeEmbedding(self.embedding, self.modes, dimensionless=True)
        if isinstance(order, (int, np.integer)):
            n = order
            order = {
                'potential':n,
                'kinetic':n,
                'pseudopotential':n-2,
            }
            if carts:
                order['coriolis'] = n-2

        if coordinate_transformation is not None:
            if isinstance(coordinate_transformation[0][0], (int, float, np.integer, np.floating)):
                coordinate_transformation = (
                    coordinate_transformation,
                    nput.matinv_deriv(coordinate_transformation, len(coordinate_transformation))
                )

        expansions = {}
        for k,o in order.items():
            if k == 'kinetic':
                generator = GMatrixExpansion(embedding)
            elif k == 'potential':
                generator = ScalarExpansion(embedding, self.potential.derivatives)
            elif k == 'coriolis':
                generator = CoriolisRotationExpansion(embedding)
            elif k == 'pseudopotential':
                generator = PseudopotentialExpansion(GMatrixExpansion(embedding))
            else:
                raise ValueError("don't know what to do with key {}".format(k))
            expansions[k] = generator.get_terms(o, transformation=coordinate_transformation)

        return expansions

    @property
    def G_terms(self):
        return GMatrixExpansion(ModeEmbedding(self.embedding, self.modes, dimensionless=True))

class ScalarOperatorManager:
    def __init__(self, manager_or_derivs):
        self.manager = manager_or_derivs
        self._derivs = None
    @property
    def derivatives(self):
        if hasattr(self.manager, 'derivatives'):
            return self.manager.derivatives
        else:
            return self.manager

class ModeEmbedding:
    """
    Provides a specialization on a `MoleculaEmbedding` to express all properties
    in terms of the attendant normal modes
    """
    def __init__(self,
                 embedding:MolecularEmbedding,
                 modes:NormalModesManager,
                 mass_weight=False,
                 dimensionless=False
                 ):
        self.embedding = embedding
        modes = modes.modes
        self.mass_weighted = mass_weight or dimensionless
        if modes is not None:
            modes = modes.basis.to_new_modes()
            if dimensionless:
                modes = modes.make_dimensionless()
            elif mass_weight:
                modes = modes.make_mass_weighted()
        self.modes = modes

    def mw_conversion(self, strip_dummies=None):
        masses = self.embedding.masses
        if strip_dummies:
            masses = masses[masses > 0]
        mvec = np.broadcast_to(
                np.asanyarray(masses)[:, np.newaxis],
                (len(masses), 3)
            ).flatten()
        return np.diag(np.sign(mvec) * np.sqrt(np.abs(mvec)))
    def mw_inverse(self, strip_dummies=None):
        masses = self.embedding.masses
        if strip_dummies:
            masses = masses[masses > 0]
        mvec = np.broadcast_to(
                np.asanyarray(masses)[:, np.newaxis],
                (len(masses), 3)
            ).flatten()
        return np.diag(np.sign(mvec) / np.sqrt(np.abs(mvec)))

    def get_mw_cartesians_by_internals(self, order=None, strip_embedding=True):
        RX = self.embedding.get_cartesians_by_internals(
                order=order,
                strip_embedding=strip_embedding
            )
        if self.mass_weighted:
            XY = self.mw_conversion()
            RX = [
                np.tensordot(tf_X, XY, axes=[-1, 0])
                for tf_X in RX
            ]
        return RX
    def get_internals_by_mw_cartesians(self, order=None, strip_embedding=True):
        XR = self.embedding.get_internals_by_cartesians(
                order=order,
                strip_embedding=strip_embedding
            )
        if self.mass_weighted:
            YX = self.mw_inverse()
            YR = []
            for X_tf in XR:
                for d in range(X_tf.ndim-1):
                    X_tf = np.tensordot(YX, X_tf, axes=[1, d])
                YR.append(X_tf)
        else:
            YR = XR
        return YR

    def get_internals_by_cartesians(self, order=None, strip_embedding=True):
        """
        expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians

        :param order:
        :param strip_embedding:
        :return:
        """
        if self.embedding.internals is None:
            if self.modes is None:
                raise NotImplementedError("not sure what's most consistent for just...plain Cartesians")
            return [self.modes.inverse]
        else:
            YR = self.get_internals_by_mw_cartesians(
                order=order,
                strip_embedding=strip_embedding
            )
            if self.modes is not None:
                YQ = self.modes.matrix
                RY = self.get_mw_cartesians_by_internals(
                    order=1,
                    strip_embedding=strip_embedding
                )[0]
                RQ = RY @ YQ

                _ = []
                for Y_tf in YR:
                    Y_tf = np.tensordot(Y_tf, RQ, axes=[-1, 0])
                    _.append(Y_tf)
                YR = _
            return YR

    def get_cartesians_by_internals(self, order=None, strip_embedding=True):
        """
        expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians

        :param order:
        :param strip_embedding:
        :return:
        """
        if self.embedding.internals is None:
            if self.modes is None:
                raise NotImplementedError("not sure what's most consistent for just...plain Cartesians")
            return [self.modes.inverse]
        else:
            if self.modes is not None: strip_embedding = False
            RY = self.get_mw_cartesians_by_internals(
                order=order,
                strip_embedding=strip_embedding
            )
            if self.modes is not None:
                QY = self.modes.inverse
                YR = self.get_internals_by_mw_cartesians(
                    order=1,
                    strip_embedding=strip_embedding
                )[0]
                QR = QY @ YR

                _ = []
                for R_tf in RY:
                    for d in range(R_tf.ndim - 1):
                        R_tf = np.tensordot(QR, R_tf, axes=[1, d])
                    _.append(R_tf)
                RY = _
            return RY

    def get_inertia_tensor_expansion(self, order=None, strip_embedding=True):
        YI0 = self.embedding.inertial_frame_derivatives()
        QY = self.get_cartesians_by_internals(order=order, strip_embedding=strip_embedding)
        return [self.embedding.inertia_tensor] + nput.tensor_reexpand(QY, YI0, order=order)

    def get_inertial_frame(self):
        return self.embedding.inertial_frame

class ScalarExpansion:
    """
    A helper class that can transform scalar derivatives from Cartesians to internals
    """
    def __init__(self,
                 embedding:ModeEmbedding,
                 derivs
                 ):
        self.embedding = embedding
        self.derivs = derivs

    def get_terms(self, order=None, *, derivs=None):

        if derivs is None:
            derivs = self.derivs
        if order is None:
            order = len(self.derivs)
        if len(derivs) < order:
            derivs = tuple(derivs) + (0,) * (order - len(derivs))

        zero_derivs = 0
        for d in derivs:
            if isinstance(d, (int, np.integer, float, np.floating)) and d == 0:
                d+=1
            else:
                break

        QX = self.embedding.get_cartesians_by_internals(order=order-zero_derivs)
        terms = TensorDerivativeConverter.convert_fast(QX, derivs).convert(order=order)  # , check_arrays=True)

        return terms

    # @classmethod
    # def get_optimized_coordinates(cls, V_expansion, order=2):
    #     V = [0] + V_expansion
    #     w = np.diag(V[1])
    #     forward_derivs = [np.eye(V[1].shape[0])]
    #     reverse_derivs = [np.eye(V[1].shape[0])]
    #
    #     w = w[np.newaxis, :]
    #     for o in range(1, order + 1):
    #         V_rem = TensorDerivativeConverter.convert_fast(forward_derivs, V,
    #                                                        order=o+2, val_axis=0
    #                                                        )[-1]
    #         w = w[np.newaxis]
    #         new_Q = -V_rem / ((o+2)*w)
    #         forward_derivs = forward_derivs + [new_Q]
    #         new_R = -nput.tensordot_deriv(forward_derivs, reverse_derivs + [0], o)[-1]
    #         reverse_derivs = reverse_derivs + [new_R]
    #
    #     return forward_derivs, reverse_derivs
    #
    # def optimize_coordinates(self, order=2):
    #     V = list(reversed([self[o] for o in range(order, -1, -1)]))
    #     RX = self.get_cartesians_by_internals(order=order, strip_embedding=True)
    #     XR = self.get_internals_by_cartesians(order=order, strip_embedding=True)
    #     QR, RQ = self.get_potential_optimized_coordinates(V, order=order)
    #
    #     QX = TensorDerivativeConverter.convert_fast(QR, RX)
    #     XQ = TensorDerivativeConverter.convert_fast(XR, RQ)
    #
    #     return (QR, RQ), (QX, XQ)

class GMatrixExpansion:
    # I thought about having this be a generic "tensor product" expansion, but this is cleaner
    def __init__(self, embedding:ModeEmbedding):
        self.embedding = embedding

    def get_terms(self, order, transformation=None):
        if transformation is not None:
            forward_derivs, reverse_derivs = transformation
            return self.reexpress(order, forward_derivs, reverse_derivs)
        XQ = self.embedding.get_internals_by_cartesians(order=order+1)
        XG = nput.tensordot_deriv(XQ, XQ, order, axes=[0, 0])#, identical=False)

        if order > 0:
            QX = self.embedding.get_cartesians_by_internals(order=order+1)
            # XG = self.reexpress_G(XG, QX, XQ, order=order)
            XG = [XG[0]] + nput.tensor_reexpand(QX, XG[1:], order)
        return XG

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
    def reexpress_G(self, G_expansion, forward_derivs, reverse_derivs, order=2):
        """
        Apply a coordinate transformation to the G-matrix

        :param forward_derivs:
        :param reverse_derivs:
        :param order:
        :return:
        """
        R = reverse_derivs
        Q = forward_derivs
        G = G_expansion

        G_R = [self._dRGQ_derivs(R, G, o) for o in range(1, order+1)]

        return [G[0]] + nput.tensor_reexpand(Q, G_R, order=order, axes=[-1, 0])

    def reexpress(self, order, forward_derivs, reverse_derivs):
        return self.reexpress_G(
            self.get_terms(order),
            forward_derivs,
            reverse_derivs,
            order=order
        )

class PseudopotentialExpansion:
    def __init__(self, g_expansion:GMatrixExpansion):
        self.g_expansion = g_expansion

    @classmethod
    def lambda_expansion(cls,
                        G_expansion, G_inv,
                        I0_expansion, I0_inv,
                        order):

        g_exp = [np.trace(t, axis1=-1, axis2=-2) for t in nput.tensordot_deriv(G_expansion[1:], G_inv, order)]
        i_exp = [np.trace(t, axis1=-1, axis2=-2) for t in nput.tensordot_deriv(I0_expansion[1:], I0_inv, order)]
        return [ti - tg for ti,tg in zip(i_exp, g_exp)]

    @classmethod
    def get_U(cls, G_expansion, I0_expansion):

        QG = G_expansion
        QI0 = I0_expansion
        n = len(QG) - 1
        order = n - 2

        QU = nput.matinv_deriv(QI0, n)
        QiG = nput.matinv_deriv(G_expansion, n)

        L = cls.lambda_expansion(
            QG, QiG,
            QI0, QU,
            order + 1
        )

        LxL = nput.tensorprod_deriv(L, L, order)

        inner_exp = [
            l1 + 1/4*l2
            for l1, l2 in zip(L[1:], LxL)
        ]
        G_contract = nput.tensordot_deriv(G_expansion, inner_exp, order, axes=[[0, 1], [0, 1]])
        QG_contract = [
            np.trace(g, axis1=-1, axis2=-2) for g in
            nput.tensordot_deriv(G_expansion[1:], L, order, axes=[2, 0])
        ]

        return [g1 + g2 for g1, g2 in zip(G_contract, QG_contract)]
        # return nput.scalarprod_deriv(Qig, t, order=order)

    def get_terms(self, order, transformation=None):
        if transformation is not None:
            forward_derivs, reverse_derivs = transformation
            return self.reexpress(order, forward_derivs, reverse_derivs)

        # build setup terms, need to go 2 orders past the nominal order
        n = order + 2
        QG = self.g_expansion.get_terms(n)
        QI0 = self.g_expansion.embedding.get_inertia_tensor_expansion(order=n)

        return self.get_U(QG, QI0)

    def reexpress(self, order, forward_derivs, reverse_derivs):
        """
        Apply a coordinate transformation to the G-matrix

        :param forward_derivs:
        :param reverse_derivs:
        :param order:
        :return:
        """

        # build setup terms, need to go 2 orders past the nominal order
        n = order + 2
        QG = self.g_expansion.reexpress(n, forward_derivs, reverse_derivs)
        QI0 = nput.tensor_reexpand(
            forward_derivs,
            self.g_expansion.embedding.get_inertia_tensor_expansion(order=n),
            n
        )

        return self.get_U(QG, QI0)

class ZetaExpansion:
    def __init__(self, embedding:ModeEmbedding):
        self.embedding = embedding
    def get_terms(self, order, transformation=None):
        QY = self.embedding.get_cartesians_by_internals(order=order)
        if transformation is not None:
            forward_derivs, reverse_derivs = transformation
            QY = nput.tensor_reexpand(forward_derivs, QY, order)

        B_e, eigs = self.embedding.get_inertial_frame()

        _ = []
        for J in QY:
            J = J.reshape(J.shape[:-2] + (-1, 3))
            J = np.tensordot(eigs, J, axes=[1, 2])
            _.append(np.moveaxis(J, -2, -1))
        QY = _

        QYQY = [
            np.sum(qyqy, axis=-1)
            for qyqy in nput.tensorprod_deriv(QY, QY, axes=[[0, 1], [0, 1]], order=order)
        ]

        zeta_expansion = []
        for qyqy in QYQY:
            for a in range(3):
                b = (a + 1) % 3
                c = (a + 2) % 3
                if b % 2 == 1:
                    c, b = b, c
                zeta_expansion.append(qyqy[..., b, :, c] - qyqy[..., c, :, b])

        return zeta_expansion

class ReciprocalInertiaExpansion:
    """
    Pulled from Watson
    """
    def __init__(self, embedding:ModeEmbedding):
        self.embedding = embedding

    def get_terms(self, order, transformation=None):
        I0, a = self.embedding.get_inertia_tensor_expansion(order=1)
        if transformation is not None:
            forward_derivs, reverse_derivs = transformation
            a = nput.tensor_reexpand(forward_derivs, [a], 1)[0]
        u = np.linalg.inv(I0)

        expansion = [u]
        for _ in range(order):
            prev = expansion[-1]
            new = np.tensordot(
                np.tensordot(a, prev, axes=[-1, -1]),
                u,
                axes=[-1, 0]
            )
            expansion.append(new)

        return expansion

class CoriolisRotationExpansion:
    """
    Provides an expansion of the Coriolis rotation operator
    """
    def __init__(self, embedding:ModeEmbedding):
        self.embedding = embedding
        self.zetas = ZetaExpansion(self.embedding)
        self.us = ReciprocalInertiaExpansion(self.embedding)

    def get_terms(self, order, transformation=None):
        Z = self.zetas.get_terms(order, transformation=transformation)
        ZZ = nput.tensorprod_deriv(Z, Z, order, axes=[1, 1]) # x/y/z shared
        u = self.us.get_terms(order, transformation=transformation)

        return nput.tensordot_deriv(u, ZZ, order, axes=[1, 0])
