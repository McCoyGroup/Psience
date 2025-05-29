
__all__ = [
    "MolecularHamiltonian"
]

import numpy as np, itertools
import McUtils.Numputils as nput
import McUtils.Iterators as itut
from McUtils.Combinatorics import UniquePermutations
from McUtils.Zachary import TensorDerivativeConverter

from .CoordinateSystems import MolecularEmbedding, ModeEmbedding
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
            }
            if carts:
                order['coriolis'] = n-2
                order['watson'] = n-2
            else:
                order['pseudopotential'] = n-2

        coordinate_transformation = self.prep_transformation(coordinate_transformation)

        expansions = {}
        for k,o in order.items():
            post_processor = None
            if k == 'kinetic':
                generator = self.gmatrix_expansion(embedding=embedding)
            elif k == 'potential':
                # ders = self.potential.derivatives
                # if np.sum(np.abs(ders[0])) > 1e-8: ...
                generator = self.potential_expansion(embedding=embedding, exclude_gradient=True)
                post_processor = lambda exp: exp[1:]
            elif k == 'coriolis':
                generator = self.coriolis_expansion(embedding=embedding)
            elif k == 'pseudopotential':
                generator = self.pseudopotential_expansion(embedding)
            elif k == 'watson':
                generator = self.watson_expansion(embedding)
            else:
                raise ValueError("don't know what to do with key {}".format(k))
            terms = generator.get_terms(o, transformation=coordinate_transformation)
            if post_processor is not None:
                terms = post_processor(terms)
            expansions[k] = terms

        return expansions

    @classmethod
    def prep_transformation(cls, coordinate_transformation):
        if coordinate_transformation is not None:
            if nput.is_numeric(coordinate_transformation[0][0][0]):
                coordinate_transformation = (
                    coordinate_transformation,
                    nput.inverse_transformation(coordinate_transformation, len(coordinate_transformation))
                )

        return coordinate_transformation

    def get_embedding(self, masses=None, modes=True, mass_weight=True, dimensionless=False):
        if modes is True:
            modes = self.modes
        elif modes is False:
            modes = None
        return ModeEmbedding(self.embedding, modes, masses=masses, mass_weight=mass_weight, dimensionless=dimensionless)

    def potential_expansion(self,
                            order=None,
                            *,
                            dimensionless=False,
                            embedding=None,
                            exclude_gradient=False,
                            **embedding_opts
                            ):
        if embedding is None:
            embedding = self.get_embedding(dimensionless=dimensionless, **embedding_opts)

        ders = self.potential.derivatives
        if exclude_gradient:
            ders = [0] + ders[1:]

        exp = ScalarExpansion(embedding, ders)
        if order is not None:
            exp = exp.get_terms(order)
        return exp
    def dipole_expansion(self,
                            order=None,
                            *,
                            expansion=None,
                            dimensionless=False,
                            embedding=None,
                            exclude_gradient=False,
                            **embedding_opts
                            ):
        if embedding is None:
            embedding = self.get_embedding(dimensionless=dimensionless, **embedding_opts)

        if expansion is None:
            expansion = self.dipole.derivatives

        exp = DipoleExpansion(embedding, expansion)
        if order is not None:
            exp = exp.get_terms(order, shared=1)
        return exp
    def get_potential_optmizing_transformation(self,
                                               order,
                                               *,
                                               dimensionless=False,
                                               embedding=None,
                                               exclude_gradient=True,
                                               **embedding_opts
                                               ):
        return nput.optimizing_transformation(
            self.potential_expansion(order,
                                     dimensionless=dimensionless,
                                     embedding=embedding,
                                     exclude_gradient=exclude_gradient,
                                     **embedding_opts
                                     ),
            order
        )

    def _get_ke_expansion(self, cls,
                          order=None,
                          *,
                          coords=None,
                          masses=None,
                          dimensionless=False,
                          embedding=None,
                          **embedding_opts
                          ):
        if embedding is None:
            embedding = self.get_embedding(masses=masses, dimensionless=dimensionless, **embedding_opts)

        exp = cls(embedding)
        if order is not None:
            exp = exp.get_terms(order, coords=coords)
        return exp
    def gmatrix_expansion(self,
                          order=None,
                          *,
                          coords=None,
                          dimensionless=False,
                          embedding=None,
                          **embedding_opts
                          ) -> "GMatrixExpansion|np.ndarray":
        return self._get_ke_expansion(GMatrixExpansion,
                                      coords=coords,
                                      order=order,
                                      dimensionless=dimensionless,
                                      embedding=embedding,
                                      **embedding_opts
                                      )
    def get_kinetic_optmizing_transformation(self,
                                             order,
                                             *,
                                             dimensionless=False,
                                             embedding=None,
                                             **embedding_opts
                                             ):
        g = self.gmatrix_expansion(dimensionless=dimensionless, embedding=embedding, **embedding_opts) #type:GMatrixExpansion
        G_r = g.get_terms(order)

        w = np.diag(G_r[0])
        # Q = [v]
        # R = [v.T]

        # Q0 = Q[0]
        # R0 = R[0]
        # G0 = G_r[0]
        Q = [np.eye(len(w))]
        R = [np.eye(len(w))]

        # w = w
        G_partial = None
        for o in range(1, order + 1):
            G_rem, G_partial = g.reexpress_G(G_r, Q, R,
                                             order=[o],
                                             # GR_expansion=G_partial,
                                             return_GR=True
                                             )
            # G_partial = G_partial[:-1]
            # w = w[np.newaxis]
            ggg = (
                          G_rem[-1] + G_rem[-1].transpose(1, 0, 2)
                   ) / 4
            # print(ggg - ggg.transpose(1, 0, 2))
            # raise ValueError("wtf")
            new_R = -ggg / (w[:, np.newaxis, np.newaxis] + w[np.newaxis, :, np.newaxis])
            # for i in range(o+1):
            #     new_R = np.tensordot(w, new_R)
            # print(w.shape)
            # print(new_R - np.transpose(new_R, (1, 0, 2)))
            # raise Exception(...)
            R = R + [new_R]

            print(new_R - new_R.transpose(1, 0, 2))

            Q = nput.inverse_transformation(R, [o], reverse_expansion=Q)
            print("="*50)
            G_rem2, G_partial = g.reexpress_G(G_r[:1], Q, R,
                                             order=[o],
                                             # GR_expansion=G_partial,
                                             return_GR=True
                                             )
            # print(G_rem[-1])
            print(G_rem2[-1] + G_rem2[-1].transpose(1, 0, 2))

            raise Exception("???")

        return Q, R

    def pseudopotential_expansion(self,
                                  order=None,
                                  *,
                                  dimensionless=False,
                                  embedding=None,
                                  **embedding_opts
                                  ):
        return self._get_ke_expansion(PseudopotentialExpansion,
                                      order=order,
                                      dimensionless=dimensionless,
                                      embedding=embedding,
                                      **embedding_opts
                                      )
    def coriolis_expansion(self,
                           order=None,
                           *,
                           dimensionless=False,
                           embedding=None,
                           **embedding_opts
                           ):
        return self._get_ke_expansion(CoriolisRotationExpansion,
                                      order=order,
                                      dimensionless=dimensionless,
                                      embedding=embedding,
                                      **embedding_opts
                                      )
    def watson_expansion(self,
                         order=None,
                         *,
                         dimensionless=False,
                         embedding=None,
                         **embedding_opts
                         ):
        return self._get_ke_expansion(ReciprocalInertiaExpansion,
                                      order=order,
                                      dimensionless=dimensionless,
                                      embedding=embedding,
                                      **embedding_opts
                                      )

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

class PotentialManager(ScalarOperatorManager):
    def __init__(self, manager_or_derivs, modes:NormalModesManager):
        super().__init__(manager_or_derivs)
        modes = modes.modes
        if modes is not None:
            modes = modes.basis.to_new_modes().make_dimensionless()
        self.modes = modes

    @property
    def derivatives(self):
        if hasattr(self.manager, 'derivatives'):
            derivs = self.manager.derivatives
        else:
            derivs = self.manager
        return self.canonicalize_derivs(derivs)


class ScalarExpansion:
    """
    A helper class that can transform scalar derivatives from Cartesians (or Cartesian modes) to internals
    """
    def __init__(self,
                 embedding:ModeEmbedding,
                 derivs
                 ):
        self.embedding = embedding
        self.derivs = derivs

    def canonicalize_derivs(self, derivs, shared=0):
        if self.embedding.modes is None: return derivs
        # in_modes = False
        # for d in derivs:
        #     if not nput.is_zero(d):
        #         in_modes = d.shape[0] != d.shape[-1]
        #         if in_modes: break
        YX = self.embedding.mw_inverse()
        _ = []
        for d in derivs:
            if not nput.is_zero(d):
                for i,x in enumerate(d.shape[shared:]):
                    if x == YX.shape[-1]:
                        d = np.moveaxis(
                            np.tensordot(YX, d, axes=[-1, shared+i]),
                            0,
                            shared + i
                        )
            _.append(d)
        # if in_modes:
        #     QX = self.embedding.modes.coords_by_modes @ YX
        #     Qq = np.diag(1/np.sqrt(self.embedding.modes.freqs))
        #     _ = []
        #     for d in derivs:
        #         if not nput.is_zero(d):
        #             for j in range(d.ndim):
        #                 if d.shape[j] == QX.shape[1]:
        #                     d = np.tensordot(QX, d, axes=[1, j])
        #                 else:
        #                     d = np.tensordot(Qq, d, axes=[1, j])
        #                 d = np.moveaxis(d, 0, j)
        #             # # symmetrize d
        #             if symmetrization_mode is not None and symmetrization_mode != 'unhandled':
        #                 nput.symmetrize_array(d, symmetrization_mode=symmetrization_mode)
        #                 # n = 0
        #                 # t = 0
        #                 # for p in itertools.permutations(range(d.ndim-1)):
        #                 #     n += 1
        #                 #     t += d.transpose(p + (d.ndim-1,))
        #                 # d = t / n
        #         _.append(d)
        #     derivs = _
        # else:
        #     derivs = nput.tensor_reexpand([YX], derivs, len(derivs))
        return derivs

    @classmethod
    def check_mixed_expansion(self, derivs, shared=0):
        real_ders = [d for d in derivs if not nput.is_zero(d)]
        if any(d.shape[shared] != d.shape[-1] for d in real_ders):
            shp = real_ders[-1].shape
            mixed_axes = itut.counts(shp)[shp[-1]]
            return True, mixed_axes
        else:
            return False, None

    def get_scalar_transformation(self, order, transformation=None):
        base_tf = self.embedding.get_internals_by_cartesians(order)
        rev_tf = self.embedding.get_cartesians_by_internals(order)
        if transformation is not None:
            forward_derivs, reverse_derivs = MolecularHamiltonian.prep_transformation(transformation)
            base_tf = nput.tensor_reexpand(forward_derivs, base_tf)
            rev_tf = nput.tensor_reexpand(reverse_derivs, rev_tf)
        return base_tf, rev_tf

    def get_terms(self,
                  order=None,
                  *,
                  derivs=None,
                  transformation=None,
                  mixed_transformation=None,
                  mixed_derivative_handling_mode='high',
                  modes=None,
                  shared=0
                  ):
        if transformation is not None:
            forward_derivs, reverse_derivs = MolecularHamiltonian.prep_transformation(transformation)
            V = self.get_terms(order, derivs=derivs, shared=shared)
            return nput.tensor_reexpand(forward_derivs, V, order, axes=[-1, shared])

        if derivs is None:
            derivs = self.derivs
        if order is None:
            order = len(derivs)
        if len(derivs) < order:
            derivs = tuple(derivs) + (0,) * (order - len(derivs))

        derivs = self.canonicalize_derivs(derivs, shared=shared)
        is_mixed, num_base = self.check_mixed_expansion(derivs)

        zero_derivs = 0
        for d in derivs:
            if nput.is_numeric(d) and d == 0:
                d+=1
            else:
                break

        QX = self.embedding.get_cartesians_by_internals(order=order - zero_derivs)
        real_ders = [d for d in derivs if not nput.is_zero(d)]
        if self.embedding.modes is not None and (
                real_ders[0].shape[0] == self.embedding.modes.modes_by_coords.shape[1]
        ):
            XQ = self.embedding.modes.modes_by_coords
            QX = nput.tensor_reexpand(QX, [XQ], order - zero_derivs, axes=[-1, 0])
        else:
            YX = self.embedding.mw_inverse()
            QX = [np.tensordot(q, YX, axes=[-1, 0]) for q in QX]


        if is_mixed:
            canonical_derivs = nput.tensor_reexpand(
                QX[:num_base], derivs[:num_base], num_base,
                axes=[-1, -1]
            ) + [
                nput.tensor_reexpand(QX, [0] * (num_base-1) + [d], num_base, axes=[-1, -1])[-1]
                for d in derivs[num_base:]
            ]
            if mixed_transformation is None:
                if modes is None:
                    modes = self.embedding.modes
                if modes.is_cartesian:
                    if self.embedding.embedding.internals is not None:
                        mixed_transformation = nput.tensor_reexpand(
                            QX,
                            [modes.modes_by_coords]
                        )
            if mixed_transformation is None:
                return canonical_derivs
            else:
                tf_base = nput.tensor_reexpand(mixed_transformation, canonical_derivs, order, axes=[-1, -1])
                return [
                    nput.symmetrize_array(t, symmetrization_mode=mixed_derivative_handling_mode)
                    for t in tf_base
                ]
        else:
            # return TensorDerivativeConverter(QX, derivs).convert()
            return nput.tensor_reexpand(QX, derivs, order, axes=[-1, shared])

class DipoleExpansion(ScalarExpansion):
    def __init__(self,
                 embedding: ModeEmbedding,
                 derivs
                 ):
        ders = [
            np.moveaxis(d, -1, 0)
                if not nput.is_zero(d) and d.shape[-1] == 3 else
            d
            for d in derivs
        ]
        super().__init__(embedding, ders)
    def get_terms(self,
                  order=None,
                  *,
                  derivs=None,
                  transformation=None,
                  mixed_transformation=None,
                  mixed_derivative_handling_mode='high',
                  modes=None,
                  shared=0
                  ):
        if derivs is not None:
            derivs = [
                np.moveaxis(d, -1, 0)
                if not nput.is_zero(d) and d.shape[-1] == 3 else
                d
                for d in derivs
            ]
        else:
            derivs = self.derivs
        ref = derivs[0]
        derivs = derivs[1:]
        terms = super().get_terms(
            order=order,
            derivs=derivs,
            transformation=transformation,
            mixed_transformation=mixed_transformation,
            mixed_derivative_handling_mode=mixed_derivative_handling_mode,
            modes=modes,
            shared=shared
        )
        expansion = [ref] + terms
        new_exp = [
            np.moveaxis(d, 0, -1)
                if not nput.is_zero(d) and d.shape[0] == 3 else
            d
            for d in expansion
        ]
        return new_exp


class GMatrixExpansion:
    # I thought about having this be a generic "tensor product" expansion, but this is cleaner
    def __init__(self, embedding:ModeEmbedding):
        self.embedding = embedding

    def get_terms(self, order, coords=None, transformation=None):
        if transformation is not None:
            forward_derivs, reverse_derivs = MolecularHamiltonian.prep_transformation(transformation)
            return self.reexpress(order, forward_derivs, reverse_derivs)
        XQ = self.embedding.get_internals_by_cartesians(order=order+1, coords=coords)
        XG = nput.tensordot_deriv(XQ, XQ, order, axes=[-2, -2], shared=XQ[0].ndim-2)#, identical=False)

        if order > 0:
            QX = self.embedding.get_cartesians_by_internals(order=order+1, coords=coords)
            # XG = self.reexpress_G(XG, QX, XQ, order=order)
            XG = [XG[0]] + nput.tensor_reexpand(QX, XG[1:], order)
        return XG

    @classmethod
    def _dRGQ_partition_contrib(cls, partition, R, G):
        # print(partition)
        r1, r2, s = partition
        if s >= len(G): return 0
        if r1 >= len(R) or r2 >= len(R): return 0

        base_term = G[s]
        if isinstance(base_term, (int, float, np.integer, np.floating)) and base_term == 0:
            return 0

        r_perm_counter = 1
        r_perm_idx = []
        # g_perm_idx = []

        a = base_term.ndim - 2
        for r in [r1, r2]:
            d = R[r]
            if nput.is_zero(d): return 0

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
                    G_expansion, forward_derivs, reverse_derivs,
                    order,
                    *,
                    GR_expansion=None,
                    return_GR=False
                    ):
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

        if nput.is_numeric(order): order = list(range(order+1))
        m = max(order)

        G_R = GR_expansion
        if G_R is None:
            G_R = [self._dRGQ_derivs(R, G, o) for o in range(1, m+1)]

        if len(G_R) < m:
            G_R = G_R + [self._dRGQ_derivs(R, G, o) for o in range(len(G_R)+1, m+1)]

        rexpansion = nput.tensor_reexpand(Q,
                                          G_R,
                                          order=[o for o in order if o > 0],
                                          axes=[-1, 0]
                                          )
        G_Q = []
        n = 0
        G0 = None
        for o in order:
            if o == 0:
                if G0 is None:
                    G0 = reverse_derivs[0] @ G[0] @ reverse_derivs[0].T
                G_Q.append(G0)
            else:
                G_Q.append(rexpansion[n])
                n += 1

        if return_GR:
            return G_Q, G_R
        else:
            return G_Q

    def reexpress(self, order, forward_derivs, reverse_derivs):
        return self.reexpress_G(
            self.get_terms(order),
            forward_derivs,
            reverse_derivs,
            order=order
        )

class PseudopotentialExpansion:
    def __init__(self, g_expansion:"ModeEmbedding|GMatrixExpansion"):
        if not isinstance(g_expansion, GMatrixExpansion):
            g_expansion = GMatrixExpansion(g_expansion)
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
            forward_derivs, reverse_derivs = MolecularHamiltonian.prep_transformation(transformation)
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

        QY = self.embedding.get_internals_by_cartesians(order=order)
        if transformation is not None:
            forward_derivs, reverse_derivs = MolecularHamiltonian.prep_transformation(transformation)
            QY = nput.tensor_reexpand(forward_derivs, QY, order)
        QY = [
            np.moveaxis(yq, -1, 0) if not nput.is_zero(yq) else 0
            for yq in QY
        ]

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
            forward_derivs, reverse_derivs = MolecularHamiltonian.prep_transformation(transformation)
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
