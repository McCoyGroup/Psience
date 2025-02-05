
import numpy as np, scipy as sp, itertools as it, collections, math, time, enum
import scipy.linalg

import McUtils.Numputils as nput
from McUtils.Scaffolding import Logger
from McUtils.Formatters import TableFormatter
from McUtils.Zachary import DensePolynomial
from McUtils.Combinatorics import IntegerPartitioner, IntegerPartitioner2D, UniquePermutations

from ..BasisReps import StateMaker
from ..Modes import NormalModes

__all__ = [
    'FranckCondonModel'
]

class State(enum.Enum):
    GroundState = 'gs'
    ExcitedState = 'es'

    @classmethod
    def get_aliases(cls):
        return {
            'ground-state': 'gs',
            'excited-state': 'es'
        }

    @classmethod
    def resolve(cls, method_name):
        if isinstance(method_name, State):
            return method_name
        method_name = method_name.lower()
        method_name = cls.get_aliases().get(method_name, method_name)
        return cls(method_name)

class RotationMethods(enum.Enum):
    LeastSquares = 'least-squares'
    Duschinsky = 'duschinsky'

    @classmethod
    def get_aliases(cls):
        return {
            'default': 'duschinsky',
            'lstsq': 'least-squares'
        }
    @classmethod
    def resolve(cls, method_name):
        if isinstance(method_name, RotationMethods):
            return method_name
        method_name = method_name.lower()
        method_name = cls.get_aliases().get(method_name, method_name)
        return cls(method_name)


class FranckCondonModel:

    # @classmethod
    # def get_simple_poly_evaluator(cls, exponents):
    #     """
    #     Provides a simple function for evaluating based on the exponents without rotation
    #     (i.e. product of 1D integrals)
    #
    #     :param exponents:
    #     :return:
    #     """
    #

    default_rotation_method = 'default'
    default_rotation_order = State.GroundState
    default_rotation_center = State.ExcitedState
    default_embedding_ref = State.ExcitedState
    default_include_rotation = True
    def __init__(self, gs_nms: NormalModes, es_nms, atoms=None, *,
                 logger=None,
                 embed=True, embedding_ref=None,
                 masses=None,
                 mass_weight=True,
                 dimensionless=False,
                 mode_selection=None,
                 mode_reordering=None,
                 rotation_method=None,
                 rotation_order=None,
                 rotation_center=None,
                 include_rotation=None,
                 rotation_blocks=None
                 ):

        self.gs_nms, self.es_nms, _ = self.prep_modes(
            gs_nms, es_nms,
            embed=embed, embedding_ref=embedding_ref,
            masses=masses,
            mass_weight=mass_weight,
            dimensionless=dimensionless,
            mode_selection=mode_selection,
            mode_reordering=mode_reordering,
        )

        self.rotation_method = rotation_method
        self.rotation_order = rotation_order
        self.rotation_center = rotation_center
        self.include_rotation = include_rotation
        self.rotation_blocks = rotation_blocks
        self.atoms = atoms

        self.logger = Logger.lookup(logger)

    @classmethod
    def prep_modes(cls, gs_nms, es_nms,
                   embed=True, embedding_ref=None,
                   masses=None,
                   mass_weight=False,
                   dimensionless=False,
                   mode_selection=None,
                   mode_reordering=None,
                   **rotation_opts
                   ):
        if embed:
            if embedding_ref is None:
                embedding_ref = cls.default_embedding_ref
            gs_nms, es_nms = cls.embed_modes(gs_nms, es_nms, masses=masses, ref=embedding_ref)

        if mass_weight:
            gs_nms = cls.mass_weight_nms(gs_nms, masses=masses)
            es_nms = cls.mass_weight_nms(es_nms, masses=masses)

        if dimensionless:
            gs_nms = cls.make_dimensionless(gs_nms, masses=masses)
            es_nms = cls.make_dimensionless(es_nms, masses=masses)

        if mode_selection is not None:
            gs_nms = gs_nms[mode_selection]
            es_nms = es_nms[mode_selection]

        if mode_reordering is not None:
            # gs_nms = gs_nms[mode_selection]
            es_nms = es_nms[mode_reordering]

        return gs_nms, es_nms, rotation_opts

    @classmethod
    def from_files(cls, gs_file, es_file,
                   logger=None, mode_selection=None,
                   mode_reordering=None,
                   internals=None, internals_ref='gs',
                   **rotation_embedding_opts):
        from ..Molecools import Molecule
        internals_ref = State.resolve(internals_ref)
        if internals_ref == State.GroundState:
            gs_mol = Molecule.from_file(gs_file, internals=internals)
            ics = gs_mol.internal_coordinates  # set up the reference
            es_mol = Molecule.from_file(es_file, internals=gs_mol.internals)
        else:
            es_mol = Molecule.from_file(es_file, internals=internals)
            ics = es_mol.internal_coordinates  # set up the reference
            gs_mol = Molecule.from_file(gs_file, internals=es_mol.internals)
        return cls.from_mols(
            gs_mol,
            es_mol,
            logger=logger,
            mode_selection=mode_selection,
            mode_reordering=mode_reordering,
            **rotation_embedding_opts
        )

    @classmethod
    def convert_internal_modes(cls, mol, nms):
        raise NotImplementedError("ignored")
        RX = mol.get_cartesians_by_internals(1, strip_embedding=True, reembed=True)[0]
        XR = mol.get_internals_by_cartesians(1, strip_embedding=True)[0]

        mat = RX @ nms.modes_by_coords # RX @ XQ
        inv = nms.coords_by_modes @ XR # QX @ XR

        # print(inv @ mat) #, nms.matrix @ nms.inverse)
        # raise Exception(...)
        freqs = nms.freqs

        return NormalModes(
            mol.internal_coordinates.system,
            mat,
            inverse=inv,
            freqs=freqs,
            origin=mol.internal_coordinates
        )

    @classmethod
    def from_mols(cls, gs, es,
                  logger=None,
                  remove_transrot=True,
                  use_internals=True,
                  embed=True,
                  mass_weight=True,
                  **rotation_embedding_opts):
        use_internals = use_internals and (gs.internals is not None)
        if use_internals:
            gg = NormalModes.from_molecule(gs, use_internals=True)
            ee = NormalModes.from_molecule(es, use_internals=True)
        elif remove_transrot:
            gg = NormalModes.from_molecule(gs, project_transrot=True, mass_weighted=True, use_internals=False)
            ee = NormalModes.from_molecule(es, project_transrot=True, mass_weighted=True, use_internals=False)
        else:
            gg = gs.normal_modes.modes.basis.to_new_modes()
            ee = es.normal_modes.modes.basis.to_new_modes()
        # if use_internals:
        #     # es = gs.get_embedded_molecule(ref=gs)
        #     # gg, ee = cls.embed_modes(gg, ee, masses=gs.atomic_masses, ref='gs')
        #     gg = cls.convert_internal_modes(gs, gg)
        #     ee = cls.convert_internal_modes(es, ee)
        #
        #     raise Exception(gg.inverse @ gg.matrix)

        return cls(
            gg, ee,
            atoms=gs.atoms,
            embed=not use_internals and embed,
            logger=logger,
            **rotation_embedding_opts
        )

    def get_overlaps(self, excitations, *, duschinsky_cutoff=None, ground_states=None,
                     return_states=True,
                     **rotation_opts
                     ):
        gs_nms, es_nms, rotation_opts = self.prep_modes(
            self.gs_nms, self.es_nms,
            **rotation_opts
        )

        excitations = self.prep_state_space(excitations, es_nms)
        ground_states = self.prep_state_space(ground_states, gs_nms)

        ret = self.get_fcfs(
            gs_nms, es_nms,
            excitations,
            **self.prep_opts(
                ground_states=ground_states,
                duschinsky_cutoff=duschinsky_cutoff,
                **rotation_opts
            )
        )

        if return_states:
            ret = (ret, (ground_states, excitations))

        return ret


    OverlapData = collections.namedtuple('OverlapData', ['alphas', 'scaling', 'center', 'gs', 'es'])
    Embedding = collections.namedtuple('Embedding', ['modes', 'center'])
    def get_overlap_data(self, **rotation_embedding_opts) -> OverlapData:
        opts = self.prep_opts(**rotation_embedding_opts)
        gs_nms, es_nms, rotation_opts = self.prep_modes(self.gs_nms, self.es_nms, **opts)
        gs_data, es_data = self.prep_overlap_args(self.gs_nms, self.es_nms)

        (alphas, modes_c, center, scaling), (L_gs, c_gs), (L_es, c_es) = self.get_overlap_gaussian_data(
            *gs_data,
            *es_data,
            **rotation_opts
        )

        return self.OverlapData(
            alphas,
            scaling,
            self.Embedding(modes_c, center),
            self.Embedding(L_gs, c_gs),
            self.Embedding(L_es, c_es)
        )

    _rot_opts = {'rotation_center', 'rotation_order', 'include_rotation', 'rotation_method'}
    @classmethod
    def prep_overlap_args(self, gs_nms, es_nms):
        # gs_nms = gs_nms.make_mass_weighted()
        # es_nms = es_nms.make_mass_weighted()

        freqs_gs = np.asanyarray(gs_nms.freqs)
        if gs_nms.frequency_scaled:
            freqs_gs = np.ones_like(freqs_gs)
        freqs_es = np.asanyarray(es_nms.freqs)
        if es_nms.frequency_scaled:
            freqs_es = np.ones_like(freqs_es)
        modes_gs = np.asanyarray(gs_nms.modes_by_coords)
        modes_es = np.asanyarray(es_nms.modes_by_coords)
        inv_gs = np.asanyarray(gs_nms.coords_by_modes)
        inv_es = np.asanyarray(es_nms.coords_by_modes)
        center_gs = np.asanyarray(gs_nms.origin).flatten()
        center_es = np.asanyarray(es_nms.origin).flatten()

        return (
            (freqs_gs, modes_gs, inv_gs, center_gs),
            (freqs_es, modes_es, inv_es, center_es),
        )

    @classmethod
    def get_poly_evaluation_plan(self, exponents, alphas=None, zpe_prod=None):
        """
        Provides a function that can take a set of indices and rotation matrices
        from the gs and es bases to the shared basis of the central Gaussian and compute
        the corresponding term contributions by considering every even
        permutation of indices that could lead to a non-zero contribution

        :param tg:
        :param te:
        :return:
        """
        key = exponents
        exponents = np.asanyarray(exponents)
        n = np.sum(exponents, dtype=int)

        max_gamma = self.df_weights(n) # we only need to compute this once

        evaluators = []
        if n == 0:
            def evaluate_zero(alphas, zpe_prod, exped_poly_coeffs, batch_size=int(1e7)):
                base_val = self.zero_point_alpha_contrib(alphas) if zpe_prod is None else zpe_prod
                if not (isinstance(exped_poly_coeffs[0], np.ndarray) and exped_poly_coeffs[0].ndim == 2):
                    base_val = np.full(len(exped_poly_coeffs[0]), base_val)
                return base_val
            evaluators.append(evaluate_zero)
        else:

            for parts in IntegerPartitioner.partitions(n): # arrays of partitions with the same total order
                good_parts = np.all(parts % 2 == 0, axis=1) # only even order partitions make sense
                parts = parts[good_parts]
                if len(parts) == 0: continue

                # find the different numbers of ways to split the terms up
                term_classes = [
                    IntegerPartitioner2D.get_partitions(exponents, terms)#[:1]
                    for terms in parts
                ]

                term_class_inds = [
                    [
                        (slice(None),) + np.where(np.ones(t.shape)) + (t.flatten(),)
                        for t in term_class
                    ]
                    for term_class in term_classes
                ]

                weights = [
                    np.array([
                        np.prod([UniquePermutations.count_permutations(counts) for counts in term])
                        for term in term_class
                        ])
                    for term_class in term_classes
                ]

                evaluators.append(
                    self.term_evaluator(parts, term_classes, term_class_inds, weights, max_gamma, alphas=alphas, zpe_prod=zpe_prod)
                )

        def evaluate_contrib(alphas, zpe_prod, poly_coeffs, batch_size=int(1e7)):
            return sum(
                [
                    ev(alphas, zpe_prod, poly_coeffs, batch_size=batch_size)
                    for ev in evaluators
                 ]
            )

        return evaluate_contrib

    @classmethod
    def zero_point_alpha_contrib(cls, alphas):
        return (np.sqrt(2 * np.pi) ** len(alphas)) / np.prod(np.sqrt(alphas))
    @classmethod
    def term_evaluator(self, exponents_list, splits_list, splits_inds_list, weights_list, gammas, alphas=None, zpe_prod=None):
        ndim = len(exponents_list[0])
        prefacs = np.array([np.prod(gammas[exponents]) for exponents in exponents_list])
        if alphas is not None:
            prefacs = prefacs * (self.zero_point_alpha_contrib(alphas) if zpe_prod is None else zpe_prod)
        prealphas = alphas

        def evaluate_contrib(alphas, zpe_prod, exped_poly_coeffs, batch_size=int(1e7)):
            contrib = [0] * len(exponents_list)

            scaling = prefacs * (
                1
                    if prealphas is not None else
                self.zero_point_alpha_contrib(alphas)
                    if zpe_prod is None else
                zpe_prod
            )

            ncoords = len(alphas)
            nperms = math.comb(ncoords, ndim)
            # try to split evenly into k batches
            num_batches = 1 + (nperms // batch_size)
            real_batch_size = nperms // num_batches
            batch_remainder = nperms % real_batch_size
            # # see how well we can split the remainder over the batches
            # subremainder = batch_remainder % num_batches
            # real_batch_size += batch_remainder // num_batches
            # batch_remainder = subremainder

            # print(":::::", exponents_list, ":::::")
            # print("    scaling:", scaling)

            batch_generator = it.combinations(range(ncoords), ndim)
            for batch_number in range(num_batches):
                num_terms = real_batch_size
                if batch_remainder > 0: # exhaust remainder bit by bit
                    num_terms += 1
                    batch_remainder -= 1
                ind_array = np.empty((num_terms, ndim), dtype=int)
                for k,ind in enumerate(it.islice(batch_generator, num_terms)):
                    ind_array[k] = ind

                new_contribs = self.evaluate_poly_contrib_chunk(
                    ind_array,
                    exponents_list,
                    splits_list,
                    splits_inds_list,
                    weights_list,
                    alphas,
                    exped_poly_coeffs
                )
                for i,c in enumerate(new_contribs):
                    contrib[i] += c#sum(cc*ss for cc,ss in zip(c, scaling))
                    # if isinstance(c, list):
                    #     contrib[i] += sum(c)
                    #     # if not isinstance(contrib[i], list):
                    #     #     if contrib[i] != 0: raise ValueError("mismatch in contribs...")
                    #     #     contrib[i] = c
                    #     # else:
                    #     #     for j,cc in enumerate(c):
                    #     #         contrib[i][j] += cc
                    # else:
                    #     contrib[i] += c

            cont = sum(c * s for c,s in zip(contrib, scaling))
            # print("     contrib:", cont)
            return cont

        return evaluate_contrib

    @staticmethod
    def _expand(pc, splits):
        return pc[:, np.newaxis, :, :] ** splits[np.newaxis, :, :, :]  # c x l x k x d
    @staticmethod
    def _contract(poly_exps, weights):
        return np.dot(
                np.prod(np.prod(poly_exps, axis=-1), axis=-1),
                weights
            )
    @classmethod
    def _index_contract(cls, pc, weights, split_inds):
        return sum(
            w * np.prod(pc[s], axis=-1)
            for w, s in zip(weights, split_inds)
        )
    @classmethod
    def evaluate_poly_chunks(cls, poly_coeffs, exps, splits, split_inds, weights, alphas, include_baseline=False):
        if isinstance(poly_coeffs, np.ndarray) and poly_coeffs.ndim == 2:
            return cls.evaluate_poly_chunks([poly_coeffs], exps, splits, split_inds, weights, alphas)[0]

        if include_baseline:
            raise NotImplementedError("need to add")

        alpha_contrib = np.prod(1 / (alphas ** (exps[np.newaxis] // 2)), axis=-1)  # c array of values
        # compute contribution from polynomial factors, can add them up immediately
        multi_contrib = np.zeros(len(poly_coeffs))
        for ii,pc in enumerate(poly_coeffs):
            # poly_exps = cls._expand(pc[:, :, :, 1], splits)
            # poly_contrib = cls._contract(poly_exps, weights)
            poly_contrib = cls._index_contract(pc, weights, split_inds)
            multi_contrib[ii] = np.dot(alpha_contrib, poly_contrib)

        return multi_contrib

    @classmethod
    def evaluate_poly_contrib_chunk(self, inds, exponents_list, splits_list, splits_inds_list,
                                    weights_list, alphas, coeffs):
        # k is the number of initial coords, d is the number of final coords
        # n is the total number of coords
        # exps is a d vector of exponents
        # splits is an l x k x d array of the way the exps divide over the initial coords
        # coeffs is an k x n x e array of exponentiated polynomial coeffs for each coord
        # inds is a c x d array where c is the number of inds in this chunk
        # alphas is a n vector

        multicoeff = not isinstance(coeffs, np.ndarray) and isinstance(coeffs[0], np.ndarray)

        e = coeffs.shape[-1] if not multicoeff else coeffs[0].shape[-1]
        c = inds.shape[0]
        k = coeffs.shape[1] if not multicoeff else coeffs[0].shape[1]
        d = inds.shape[1]
        if multicoeff:
            poly_coeffs = []
            for cc in coeffs:
                pcs = np.empty((c, k, d, e), dtype=cc.dtype)
                for i in range(k): pcs[:, i, :, :] = cc[inds, i, :]
                cc = pcs
                poly_coeffs.append(cc)
        else:
            poly_coeffs = np.empty((c, k, d, e), dtype=coeffs.dtype)
            for i in range(k): poly_coeffs[:, i, :, :] = coeffs[inds, i, :]

        # compute extra alpha contrib to integrals
        subalphas = alphas[inds] # c x d array of non-zero alphas

        contribs = []
        for exps, splits, splits_inds, weights in zip(exponents_list, splits_list, splits_inds_list, weights_list):
            poly_chunks = self.evaluate_poly_chunks(
                poly_coeffs, exps, splits, splits_inds, weights, subalphas
            )
            contribs.append(poly_chunks)

        return contribs


    duschinsky_cutoff = 1e-20
    evaluator_plans = {}
    @classmethod
    def evaluate_shifted_poly_overlap(self,
                                      poly:'HermiteProductPolynomial',
                                      Q,
                                      alphas,
                                      zpe_prod,
                                      duschinsky_cutoff=None
                                      ):
        if duschinsky_cutoff is None: duschinsky_cutoff = self.duschinsky_cutoff

        contrib = 0

        work_queues = {}
        max_exponent = max([0] + [
            np.max(e) if len(e) > 0 else 0
            for e, _, _ in poly.term_iter(filter=lambda x,nzi,nzc:sum(x[i] for i in nzi)%2==0)
        ])
        if Q is not None:
            Q_sels = Q[:, :, np.newaxis] ** np.arange(max_exponent+1)[np.newaxis, np.newaxis, :]
        else:
            Q_sels = None
        for exponents, scaling, inds in poly.term_iter(filter=lambda x,nzi,nzc:sum(x[i] for i in nzi)%2==0): # even/odd
            # scaling = scaling
            exponents = np.asanyarray(exponents)
            ordering = np.argsort(-exponents)
            exponents = exponents[ordering]
            inds = [inds[o] for o in ordering]
            poly_coeffs = Q_sels[:, inds, :] if Q is not None else None

            # put coeffs on the queue so we can batch the calls
            key = tuple(exponents)
            queue = work_queues.get(key, None)
            if queue is None:
                queue = []
                work_queues[key] = queue
            queue.append([poly_coeffs, scaling])

        for key,queue in work_queues.items():
            exponents = np.array(key)
            # print("="*10, exponents, "="*10)
            plan = self.evaluator_plans.get(key, None)
            if plan is None:
                if Q is None:
                    plan = self.get_simple_poly_evaluator(exponents)
                else:
                    plan = self.get_poly_evaluation_plan(exponents)
                self.evaluator_plans[key] = plan
            coeff_stack = [poly for poly,scaling in queue if abs(scaling) > duschinsky_cutoff]
            coeff_weights = [scaling for poly,scaling in queue if abs(scaling) > duschinsky_cutoff]
            if len(coeff_weights) > 0:
                subcontribs = plan(alphas, zpe_prod, coeff_stack)
                # print("--->", len(coeff_stack), subcontribs, coeff_weights)
                wtf = sum(w*c for w,c in zip(coeff_weights, subcontribs))
                contrib += wtf
                # print("   >", wtf, contrib)

        return contrib

    @classmethod
    def df_weights(cls, n):
        # gives double factorial weights up to order n
        weights = np.zeros(n+1, dtype=float)
        inds = np.arange(0, n+1, 2)
        double_fac = np.concatenate([[1], np.cumprod(np.arange(1, n, 2))])
        weights[inds] = double_fac
        return weights

    @classmethod
    def get_overlap_gaussian_data(cls,
                                  freqs_gs, modes_gs, inv_gs, center_gs,
                                  freqs_es, modes_es, inv_es, center_es,
                                  rotation_method=None,
                                  rotation_order=None,
                                  rotation_center=None,
                                  include_rotation=None,
                                  rotation_blocks=None
                                  ):

        if include_rotation is None:
            include_rotation = cls.default_include_rotation
        if rotation_order is None:
            rotation_order = cls.default_rotation_order
        if rotation_method is None:
            rotation_method = cls.default_rotation_method
        if rotation_center is None:
            rotation_center = cls.default_rotation_center

        gs_order = State(rotation_order) == State.GroundState
        L_e = np.eye(inv_gs.shape[0])
        L_g = np.eye(inv_gs.shape[0])

        if include_rotation:
            if rotation_blocks is None:
                rotation_blocks = [np.arange(L_e.shape[0])]
            rotation_method = RotationMethods.resolve(rotation_method)
            for block in rotation_blocks:
                if rotation_method == RotationMethods.LeastSquares:
                    # use least squares to express ground state modes as LCs of excited state modes
                    if gs_order:
                        m1 = modes_es[:, block]
                        m2 = modes_gs[:, block]
                    else:
                        m1 = modes_gs[:, block]
                        m2 = modes_es[:, block]
                    ls_tf = np.linalg.inv(m1.T @ m1) @ m1.T @ m2
                    U, s, V = np.linalg.svd(ls_tf)
                    L = U @ V  # Unitary version of a Duschinsky matrix
                    L = L.T
                    if gs_order:
                        l_g = L
                        l_e = np.eye(L.shape[0])
                    else:
                        l_e = L
                        l_g = np.eye(L.shape[0])
                elif rotation_method == RotationMethods.Duschinsky:
                    if gs_order:
                        # del_G E
                        l_g = inv_gs[block, :] @ modes_es[:, block]
                        # l_g = (inv_gs[block, :] @ modes_es[:, block]).T
                        l_e = np.eye(l_g.shape[0])
                    else:
                        # del_E G
                        l_e = inv_es[block, :] @ modes_gs[:, block]
                        # l_e = (inv_gs[block, :] @ modes_es[:, block]).T
                        l_g = np.eye(l_e.shape[0])
                elif callable(rotation_method):
                    l_g, l_e = rotation_method(rotation_order, block)
                else:
                    raise ValueError(f"don't understand rotation method {rotation_method}")
                L_e[np.ix_(block, block)] = l_e
                L_g[np.ix_(block, block)] = l_g

        # find displacement vector (center of gs ground state in this new basis)
        # modes_gs = modes_es @ L

        es_center = State(rotation_center) == State.ExcitedState
        if gs_order:
            tf = modes_es
        else:
            tf = modes_gs
        c_g = center_gs
        c_e = center_es
        if es_center:
            dx_g = (c_g - c_e) @ tf
            dx_e = np.zeros(dx_g.shape)
        else:
            dx_e = (c_e - c_g) @ tf
            dx_g = np.zeros(dx_e.shape)

        # inverse covariance matrices in this coordinate system

        # del_E G or del_G E
        Z_gs = (L_g @ np.diag(freqs_gs) @ L_g.T)
        Z_es = (L_e @ np.diag(freqs_es) @ L_e.T)
        Z_c = Z_gs + Z_es

        # X_g = nput.fractional_power(L_g, -1)
        # print(L_g)
        # print(X_g)
        # X_e = nput.fractional_power(L_e, -1)
        X_g = np.linalg.inv(L_g)
        X_e = np.linalg.inv(L_e)
        S_gs = (X_g.T @ np.diag(1/freqs_gs) @ X_g)
        S_es = (X_e.T @ np.diag(1/freqs_es) @ X_e)
        S_c = S_gs + S_es

        norm_gs = np.power(np.linalg.eigvalsh(Z_gs), 1/4)
        norm_es = np.power(np.linalg.eigvalsh(Z_es), 1/4)
        # prefactor = (
        #         np.prod(norm_gs * norm_es) / np.power(np.pi, len(freqs_es) / 2)
        #         * np.exp(-1/2 * np.dot(dx_g+dx_e, np.dot(np.linalg.inv(S_c), dx_g+dx_e))) # one of these is always zero
        # )
        # prefactor = 1
        decay = (
                np.exp(-1/2 * np.dot(dx_g+dx_e, np.dot(np.linalg.inv(S_c), dx_g+dx_e))) # one of these is always zero
        )
        # print(np.linalg.norm(dx_g+dx_e))
        # print(np.dot(dx_g + dx_e, np.dot(np.linalg.inv(S_c), dx_g + dx_e)))

        alphas_c, modes_c = np.linalg.eigh(Z_c)
        center = np.linalg.inv(Z_c) @ (Z_es @ dx_e + Z_gs @ dx_g)

        zpe_prod = decay * np.prod(norm_gs*norm_es/np.sqrt(alphas_c/2))
        # zpe_prod = prefactor * self.zero_point_alpha_contrib(alphas_c)
        # zpe_prod = decay * self.zero_point_alpha_contrib(norm_gs*norm_es*alphas_c / np.power(np.pi, len(freqs_es) / 2))

        return (alphas_c, modes_c, center, zpe_prod), (L_g, dx_g), (L_e, dx_e)

    integral_block_size = 1000
    @classmethod
    def eval_fcf_overlaps(self,
                          excitations_gs, freqs_gs, modes_gs, inv_gs, center_gs,
                          excitations_es, freqs_es, modes_es, inv_es, center_es,
                          duschinsky_cutoff=None,
                          logger=None,
                          **rotation_opts
                          ):
        """
        Evaluates the Gaussian overlaps between two H.O. wave functions defined by
        a set of polynomial coefficients, broadening factors, and centers, assuming
        the modes and centers are in an Eckart fream
        """

        logger = Logger.lookup(logger)

        # The steps to this process
        # 0. Express ground state nms. w.r.t excited state frame
        # 1. find the overlap Gaussian parameters
        # 2. express the shifts w.r.t the corresponding Gaussian frames
        # 3. shift the polynomials so they are centered at the overlap Gaussian (but don't rotate)
        # 4. express the normal modes (transformation coeffs) w.r.t. the central Gaussian coords
        # 5. for every pair of polynomial terms, if even, evaluate polynomial contrib, reusing when possible

        freqs_gs = np.asanyarray(freqs_gs)
        freqs_es = np.asanyarray(freqs_es)
        modes_gs = np.asanyarray(modes_gs)
        modes_es = np.asanyarray(modes_es)
        inv_gs = np.asanyarray(inv_gs)
        inv_es = np.asanyarray(inv_es)
        center_gs = np.asanyarray(center_gs)
        center_es = np.asanyarray(center_es)

        (alphas, modes_c, center, zpe_prod), (L_gs, c_gs), (L_es, c_es) = self.get_overlap_gaussian_data(
            freqs_gs, modes_gs, inv_gs, center_gs,
            freqs_es, modes_es, inv_es, center_es,
            **rotation_opts
        )

        with logger.block(tag='Duschinksy Transformation'):
            logger.log_print(
                L_gs,
                message_prepper=logger.prep_array
            )
            logger.log_print(
                "disp: {av}",
                av=c_es-c_gs
            )

        # We now have the space to define the parameters that go into the overlap calculation
        shift_gs = center - c_gs #- L_gs @ modes_c.T @ center
        shift_es = center - c_es #- modes_c.T @ center

        Q_gs = (modes_c @ L_gs)
        Q_es = (modes_c @ L_es)

        start = time.time()
        with logger.block(tag="Computing base polynomials"):
            polys1 = [
                HermiteProductPolynomial.from_quanta(eg, freqs_gs).shift(shift_gs)
                for eg in excitations_gs
            ]
            polys2 = [
                HermiteProductPolynomial.from_quanta(ee, freqs_es).shift(shift_es)
                for ee in excitations_es
            ]
            end = time.time()
            logger.log_print("   took: {e:.3f}s", e=(end-start))

        Q = np.concatenate([Q_gs, Q_es], axis=1)


        with logger.block(tag="Evaluating Overlaps"):
            overlaps = []
            n = 0
            block_counter = 1
            block_size = self.integral_block_size
            full_blocks = len(excitations_gs) * len(excitations_es)
            num_blocks = int(np.ceil(full_blocks / block_size))
            block_start = None
            logger.log_print("num integrals: {n}", n=full_blocks)
            for eg,spoly_gs in zip(excitations_gs,polys1):
                subovs = []
                subcounter = 0
                with logger.block(
                        tag="Overlaps from: {gs}",
                        gs=eg,
                        preformatter=lambda **kw:dict(kw, gs=StateMaker.parse_state(kw['gs']))
                    ):
                    for spoly_es in polys2:
                        if n % block_size == 0:
                            if block_start is not None:
                                block_end = time.time()
                                logger.log_print(
                                    "{fmt_table}",
                                    es=excitations_es,
                                    ov=subovs[(subcounter - 1)*block_size:(subcounter)*block_size],
                                    preformatter=lambda **kw: dict(kw,
                                                                   fmt_table=self.format_overlap_tables(
                                                                       kw['es'],
                                                                       kw['ov'],
                                                                       include_headers=False
                                                                   )
                                                                   )
                                )
                                logger.log_print("   took: {e:.3f}s", e=block_end-block_start)
                            block_start = time.time()
                            logger.log_print("evaluating block {} of {}".format(
                                block_counter, num_blocks
                            ))
                            block_counter += 1
                            subcounter += 1
                        n += 1
                        poly = spoly_gs.concat(spoly_es)
                        ov = self.evaluate_shifted_poly_overlap(
                            poly, Q,
                            alphas,
                            zpe_prod,
                            duschinsky_cutoff=duschinsky_cutoff
                        )
                        # print("<::", scaling, ov/scaling, ov)
                        subovs.append(ov)
                logger.log_print(
                    "{fmt_table}",
                    es=excitations_es,
                    ov=subovs,
                    preformatter=lambda **kw:dict(kw,
                                                  fmt_table=self.format_overlap_tables(kw['es'], kw['ov'])
                                                  )
                )

                overlaps.extend(subovs)

        return overlaps

    @classmethod
    def embed_modes(cls, gs_nms: 'NormalModes', es_nms, ref=None, masses=None):
        from ..Molecools import Molecule

        if masses is None:
            masses = gs_nms.masses
        if masses is None:
            masses = es_nms.masses
        if masses is None:
            raise ValueError("need masses to reembed normal modes")
        if ref is None:
            ref = cls.default_embedding_ref

        gs_embedding = State(ref) == State.GroundState

        if not gs_embedding:
            es_nms, gs_nms = gs_nms, es_nms

        mw_gs = gs_nms.mass_weighted
        mw_es = es_nms.mass_weighted
        gs_nms = gs_nms.remove_mass_weighting()
        es_nms = es_nms.remove_mass_weighting()

        gs = Molecule(['H'] * len(masses), gs_nms.origin, masses=masses)

        embedding = gs.get_embedding_data(es_nms.origin)
        ref_coords = embedding.reference_data.coords # embedded in PAF
        emb_ax = 0
        cart_ax = 1
        mat = gs_nms.matrix.reshape(ref_coords.shape + (-1,))
        ref_mat = np.moveaxis(
            np.tensordot(embedding.reference_data.axes, mat, axes=[emb_ax, cart_ax]),
            0, cart_ax
        ).reshape(gs_nms.matrix.shape)
        # emb_ax = (emb_ax + 1) % 2
        # cart_ax = 2
        # mat = gs_nms.inverse.reshape((-1,) + ref_coords.shape)
        # ref_inv = np.moveaxis(
        #     np.tensordot(embedding.reference_data.axes, mat,
        #                  axes=[emb_ax, cart_ax]),
        #     0, cart_ax
        # ).reshape(gs_nms.inverse.shape)
        gs_modes = NormalModes(
            gs_nms.basis,
            ref_mat,
            origin=ref_coords,
            inverse=ref_mat.T / np.repeat(gs_nms.masses, 3)[np.newaxis, :],
            freqs=gs_nms.freqs,
            masses=gs_nms.masses,
            mass_weighted=False
        )

        rot = embedding.rotations[0]
        emb_coords = (embedding.coord_data.coords[0]) @ rot.T # embedded in PAF,
        # raise Exception(emb_coords, ref_coords, rot)
        emb_ax = 0
        cart_ax = 1
        mat = es_nms.matrix.reshape(emb_coords.shape + (-1,))
        rot = embedding.coord_data.axes[0].T
        emb_mat = np.moveaxis(
            np.tensordot(rot, mat, axes=[emb_ax, cart_ax]),
            0, cart_ax
        ).reshape(es_nms.matrix.shape)
        # emb_ax = (emb_ax + 1) % 2
        # cart_ax = 2
        # mat = es_nms.inverse.reshape((-1,) + ref_coords.shape)
        # emb_inv = np.moveaxis(
        #     np.tensordot(rot, mat, axes=[emb_ax, cart_ax]),
        #     0, cart_ax
        # ).reshape(es_nms.inverse.shape)
        es_modes = NormalModes(
            gs_modes.basis,
            emb_mat,
            inverse=emb_mat.T / np.repeat(es_nms.masses, 3)[np.newaxis, :],
            origin=emb_coords,
            freqs=es_nms.freqs,
            masses=es_nms.masses,
            mass_weighted=False
        )

        if mw_gs:
            gs_modes = gs_modes.make_mass_weighted()
        if mw_es:
            es_modes = es_modes.make_mass_weighted()

        if not gs_embedding:
            gs_modes, es_modes = es_modes, gs_modes

        return gs_modes, es_modes

    @classmethod
    def mass_weight_nms(cls, nms, masses=None):
        # if masses is None: masses = nms.masses
        new = nms.make_mass_weighted(masses=masses)
        return new

    @classmethod
    def make_dimensionless(cls, nms, freqs=None, masses=None):
        # if masses is None: masses = nms.masses
        new = nms.make_dimensionless(freqs=freqs, masses=masses)
        return new

    @classmethod
    def prep_states_from_threshold_and_quanta(cls, nms, *,
                                              threshold=None, min_freq=None, max_state=None,
                                              min_quanta=None, max_quanta=None
                                              ):
        from ..BasisReps import BasisStateSpace

        if threshold is None and max_state is None and max_quanta is None:
            raise ValueError("need at least one filter criterion out of `threshold`, `max_state`, or `max_quanta`")

        freqs = nms.freqs
        if threshold is None:
            threshold = freqs[-1] * (np.sum(max_state) if max_state is not None else max_quanta)

        return BasisStateSpace.states_under_freq_threshold(
            freqs,
            threshold,
            min_freq=min_freq,
            max_state=max_state,
            min_quanta=min_quanta,
            max_quanta=max_quanta
        )

    @classmethod
    def prep_states_from_excitations(cls, nms, *, states, **opts):
        from ..BasisReps import StateMaker
        state = StateMaker(len(nms.freqs), **opts)
        return np.array([
            state(*exc)
            for exc in states
        ])

    state_space_prep_registry = {}
    @classmethod
    def state_space_prep_dispatchers(cls):
        return dict(
            {
                'threshold': cls.prep_states_from_threshold_and_quanta,
                'states': cls.prep_states_from_excitations
            },
            **cls.state_space_prep_registry
        )
    @classmethod
    def dispatch_state_space_prep(cls, spec, nms):
        method = spec.get('method', None)
        if method is None:
            # infer method from spec
            if 'threshold' in spec or 'quanta' in spec:
                method = 'threshold'
            elif len(spec.keys() & cls.state_space_prep_dispatchers()) > 0:
                method = next(iter(spec.keys() & cls.state_space_prep_dispatchers()))
            else:
                raise ValueError("can't infer state space method from {}".format(spec))
        spec = spec.copy()
        if 'method' in spec: del spec['method']

        prepper = cls.state_space_prep_dispatchers().get(method)

        return prepper(nms, **spec)

    @classmethod
    def _check_listable(cls, excitations):
        try:
            return all(
                isinstance(exc, np.ndarray) or all(
                    nput.is_numeric(e) or all(nput.is_numeric(ee) for ee in e)
                    for e in exc
                )
                for exc in excitations
            )
        except TypeError:
            return False
    @classmethod
    def prep_state_space(cls, excitations, nms, check=True):
        """
        Dispatcher to get appropriate state spaces
        :param excitations:
        :param check:
        :return:
        """
        if excitations is None: return excitations

        refilter = not check
        if check:
            try:
                exc0 = excitations[0][0]
            except (IndexError, KeyError, TypeError):
                refilter = True
            else:
                refilter = (
                    not isinstance(exc0, int)
                    or any(len(e) < len(nms.freqs) for e in excitations)
                )

        if refilter:
            if isinstance(excitations, dict):
                # dispatch on methods
                excitations = cls.dispatch_state_space_prep(excitations, nms)
            elif isinstance(excitations, (int, np.integer)) or isinstance(excitations[0], (int, np.integer)):
                from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis as HO
                excitations = BasisStateSpace.from_quanta(
                    HO(len(nms.freqs)), excitations
                ).excitations
            elif cls._check_listable(excitations):
                excitations = cls.dispatch_state_space_prep(
                    {
                        "states": excitations
                    }, nms
                )

            else:
                raise ValueError("can't build a state space from {}".format(excitations))

        return excitations


    @classmethod
    def get_fcfs(cls,
                 gs_nms: 'NormalModes', es_nms: 'NormalModes',
                 excitations, ground_states=None,
                 duschinsky_cutoff=None,
                 logger=None,
                 **rotation_embedding_opts
                 ):

        gs_nms, es_nms, rotation_opts = cls.prep_modes(
            gs_nms, es_nms,
            **rotation_embedding_opts
        )
        gs_args, es_args = cls.prep_overlap_args(gs_nms, es_nms)

        ndim = len(gs_nms.freqs)
        if ground_states is None:
            ground_states = [[0] * ndim]

        return cls.eval_fcf_overlaps(
            ground_states, *gs_args,
            excitations, *es_args,
            duschinsky_cutoff=duschinsky_cutoff,
            logger=logger,
            **rotation_opts
        )
    @classmethod
    def format_overlap_tables(cls, es, overlaps, include_headers=True):

        fmt_table = TableFormatter(
            [StateMaker.parse_state, ".3e"],
            headers=["State", "Overlap"] if include_headers else None
        ).format(
            [[s, o] for s, o in zip(es, overlaps)]
        )

        return fmt_table
    @classmethod
    def get_fcf_spectrum(self,
                         gs_nms: 'NormalModes', es_nms: 'NormalModes',
                         excitations, ground_states=None,
                         logger=None,
                         duschinsky_cutoff=None,
                         return_states=False,
                         **rotation_embedding_opts
                         ):
        from ..Spectra import DiscreteSpectrum
        from McUtils.Data import UnitsData

        gs_nms, es_nms, rotation_opts = self.prep_modes(
            gs_nms, es_nms,
            **rotation_embedding_opts
        )


        excitations = self.prep_state_space(excitations, es_nms)
        ground_states = self.prep_state_space(ground_states, gs_nms)

        overlaps = self.get_fcfs(
            gs_nms, es_nms,
            excitations, ground_states=ground_states,
            embed=False, masses=None, mass_weight=False,
            duschinsky_cutoff=duschinsky_cutoff,
            logger=logger,
            **rotation_opts
        )

        freqs = np.dot(np.asanyarray(excitations), es_nms.freqs) * UnitsData.convert("Hartrees", "Wavenumbers")
        ints = np.power(overlaps, 2)

        if ground_states is not None:
            ret = []
            bs = len(freqs)
            for i in range(len(ground_states)):
                ret.append(
                    DiscreteSpectrum(freqs, ints[i*bs:(i+1)*bs])
                )
        else:
            ret = DiscreteSpectrum(freqs, np.power(overlaps, 2))

        if return_states:
            ret = (ret, (ground_states, excitations))

        return ret

    def prep_opts(self, **opts):
        return dict(
            dict(
                dict(
                    logger=self.logger,
                    rotation_order=self.rotation_order if self.rotation_order is not None else self.default_rotation_order,
                    rotation_center=self.rotation_center if self.rotation_center is not None else self.rotation_center,
                    rotation_method=self.rotation_method if self.rotation_method is not None else self.rotation_method,
                    rotation_blocks=self.rotation_blocks,
                    include_rotation=self.include_rotation if self.include_rotation is not None else self.default_include_rotation
                ),
                **opts,
            ),
            embed=False,
            mass_weight=False,
            dimensionless=False,
            mode_reordering=None
        )

    def get_spectrum(self,
                     excitations, *, ground_states=None,
                     return_states=False,
                     duschinsky_cutoff=None,
                     **rotation_embedding_opts
                     ):

        return self.get_fcf_spectrum(
            self.gs_nms, self.es_nms,
            excitations,
            **self.prep_opts(
                return_states=return_states,
                ground_states=ground_states,
                duschinsky_cutoff=duschinsky_cutoff,
                **rotation_embedding_opts
            )
        )

    def get_ezFCF_input(self, excitations, atoms=None, ground_states=None, **rotation_embedding_opts):
        from .ezFCFInterface import ezFCFInterface

        return ezFCFInterface(
            self.atoms if atoms is None else atoms,
            self.gs_nms,
            self.es_nms,
            excitations,
            **{
                k:(v.value if hasattr(v, 'value') else v)
                for k,v in
                self.prep_opts(
                    ground_states=ground_states,
                    # duschinsky_cutoff=duschinsky_cutoff,
                    **rotation_embedding_opts
                ).items()
            }
        )


class HermiteProductPolynomial:

    def __init__(self, poly_dict: dict, ndim):
        self.polys = poly_dict
        self.ndim = ndim

    def shift(self, s):
        return HermiteProductPolynomial(
            {
                k: p.shift(s[k])
                for k, p in self.polys.items()
            },
            self.ndim
        )

    def concat(self, other:'HermiteProductPolynomial'):
        new_poly = self.polys.copy()
        for k,p in other.polys.items():
            new_poly[k + self.ndim] = p
        return HermiteProductPolynomial(
            new_poly,
            self.ndim + other.ndim
        )

    def term_iter(self, filter=None):
        term_coords = list(self.polys.keys())
        term_orders = [
            self.polys[k].coeffs.shape[0]
            for k in term_coords
        ]
        term_inds = [
            k for k, n in enumerate(term_orders)
            if n > 1
        ]
        term_coords = [term_coords[i] for i in term_inds]
        # term_orders = [term_orders[i] for i in term_inds]

        terms_iter = it.product(*[range(o) for o in term_orders])
        for term in terms_iter:
            if filter is None or filter(term, term_inds, term_coords):
                # term_list = sum(([i]*k for i,k in zip(term_coords, term)), [])
                sub_inds = [i for i in term_inds if term[i] > 0]
                term_prod = np.prod([
                    self.polys[i].coeffs[k]
                    for i,k in zip(term_coords, term)
                ])
                term = [term[i] for i in sub_inds]
                yield term, term_prod, [term_coords[i] for i in sub_inds]

    @classmethod
    def from_quanta(cls, quanta, alphas):
        return cls(
            {
                k: cls.get_1D_hermite_poly(n, a)
                for k, (n, a) in enumerate(zip(quanta, alphas))
                if n > 0
            },
            len(alphas)
        )

    _hermite_cache = {}
    @classmethod
    def _hermite_coeffs(cls, n):
        if n not in cls._hermite_cache:
            cls._hermite_cache[n] = np.flip(np.asarray(sp.special.hermite(n, monic=False)))
        return cls._hermite_cache[n]
    @classmethod
    def get_1D_hermite_poly(cls, n, a):
        return DensePolynomial(
            cls._hermite_coeffs(n) * np.sqrt(
                a ** np.arange(n + 1) / (2 ** (n) * math.factorial(n))
            )
        )