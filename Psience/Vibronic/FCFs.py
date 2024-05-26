
import numpy as np, scipy as sp, itertools as it
from McUtils.Zachary import DensePolynomial
from McUtils.Combinatorics import IntegerPartitioner, IntegerPartitioner2D, UniquePermutations

from ..Modes import NormalModes

__all__ = [
    'FranckCondonModel'
]

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



    @classmethod
    def get_poly_evaluation_plan(self, exponents, alphas=None):
        """
        Provides a function that can take a set of indices and rotation matrices
        from the gs and es bases to the shared basis of the central Gaussian and compute
        the corresponding term contributions by considering every even
        permutation of indices that could lead to a non-zero contribution

        :param tg:
        :param te:
        :return:
        """
        exponents = np.asanyarray(exponents)
        n = np.sum(exponents, dtype=int)

        max_gamma = self.df_weights(n) # we only need to compute this once

        evaluators = []
        if n == 0:
            def evaluate_zero(alphas, poly_coeffs, batch_size=int(1e7)):
                base_val = np.sqrt(2 * np.pi) ** (len(alphas)) / np.prod(np.sqrt(alphas))
                if not (isinstance(poly_coeffs, np.ndarray) and poly_coeffs.ndim == 2):
                    base_val = np.full(len(poly_coeffs), base_val)
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

                weights = [
                    np.array([
                        np.prod([UniquePermutations.count_permutations(counts) for counts in term])
                        for term in term_class
                        ])
                    for term_class in term_classes
                ]

                # print(">---- exps -------:", exponents)
                # for e,s,w in zip(parts, term_classes, weights):
                #     print(e)
                #     print(s)
                #     print(w)
                #     print(":---------------<")

                evaluators.append(
                    self.term_evaluator(parts, term_classes, weights, max_gamma, alphas=alphas)
                )

        def evaluate_contrib(alphas, poly_coeffs, batch_size=int(1e7)):
            return sum(
                [
                    ev(alphas, poly_coeffs, batch_size=batch_size)
                    for ev in evaluators
                 ]
            )

        return evaluate_contrib

    @classmethod
    def zero_point_alpha_contrib(cls, alphas):
        return (np.sqrt(2 * np.pi)**len(alphas))/np.prod(np.sqrt(alphas))
    @classmethod
    def term_evaluator(self, exponents_list, splits_list, weights_list, gammas, alphas=None):
        ndim = len(exponents_list[0])
        prefacs = np.array([np.prod(gammas[exponents]) for exponents in exponents_list])
        if alphas is not None:
            prefacs = prefacs * self.zero_point_alpha_contrib(alphas)
        prealphas = alphas

        def evaluate_contrib(alphas, poly_coeffs, batch_size=int(1e7)):
            contrib = [0] * len(exponents_list)

            scaling = prefacs * (1 if prealphas is not None else self.zero_point_alpha_contrib(alphas))

            ncoords = len(alphas)
            nperms = np.math.comb(ncoords, ndim)
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
                    weights_list,
                    alphas,
                    poly_coeffs
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

    @classmethod
    def evaluate_poly_chunks(cls, poly_coeffs, exps, splits, weights, alphas, include_baseline=False):
        if isinstance(poly_coeffs, np.ndarray) and poly_coeffs.ndim == 2:
            return cls.evaluate_poly_chunks([poly_coeffs], exps, splits, weights, alphas)[0]

        if include_baseline:
            raise NotImplementedError("need to add")

        alpha_contrib = np.prod(1 / (alphas ** (exps[np.newaxis] // 2)), axis=-1)  # c array of values
        # compute contribution from polynomial factors, can add them up immediately
        multi_contrib = np.zeros(len(poly_coeffs))
        for ii,pc in enumerate(poly_coeffs):
            poly_exps = pc[:, np.newaxis, :, :] ** splits[np.newaxis, :, :, :]  # c x l x k x d
            poly_contrib = np.dot(
                np.prod(np.prod(poly_exps, axis=-1), axis=-1),
                weights
            )
            # print("??", alpha_contrib)
            # print("??", poly_contrib)
            multi_contrib[ii] = np.dot(alpha_contrib, poly_contrib)

        return multi_contrib

    @classmethod
    def evaluate_poly_contrib_chunk(self, inds, exponents_list, splits_list, weights_list, alphas, coeffs):
        # k is the number of initial coords, d is the number of final coords
        # n is the total number of coords
        # exps is a d vector of exponents
        # splits is an l x k x d array of the way the exps divide over the initial coords
        # coeffs is a k x n array of polynomial coeffs for each coord
        # inds is a c x d array where c is the number of inds in this chunk
        # alphas is a n vector

        multicoeff = not isinstance(coeffs, np.ndarray) and isinstance(coeffs[0], np.ndarray)

        c = inds.shape[0]
        k = coeffs.shape[1] if not multicoeff else coeffs[0].shape[1]
        d = inds.shape[1]
        if multicoeff:
            poly_coeffs = []
            for cc in coeffs:
                pcs = np.empty((c, k, d), dtype=cc.dtype)
                for i in range(k): pcs[:, i, :] = cc[:, i][inds]
                cc = pcs
                poly_coeffs.append(cc)
        else:
            poly_coeffs = np.empty((c, k, d), dtype=coeffs.dtype)
            for i in range(k): poly_coeffs[:, i, :] = coeffs[:, i][inds]

        # compute extra alpha contrib to integrals
        subalphas = alphas[inds] # c x d array of non-zero alphas

        contribs = []
        for exps, splits, weights in zip(exponents_list, splits_list, weights_list):
            poly_chunks = self.evaluate_poly_chunks(
                poly_coeffs, exps, splits, weights, subalphas
            )
            contribs.append(poly_chunks)

        return contribs



    evaluator_plans = {}
    @classmethod
    def evaluate_shifted_poly_overlap(self,
                                      poly:'HermiteProductPolynomial',
                                      Q,
                                      alphas
                                      ):

        contrib = 0

        work_queues = {}
        for exponents, scaling, inds in poly.term_iter(filter=lambda x,nzi,nzc:sum(x[i] for i in nzi)%2==0):
            exponents = np.array(exponents)
            ordering = np.argsort(-exponents)
            exponents = exponents[ordering]
            poly_coeffs = Q[:, inds][:, ordering] if Q is not None else None

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
            coeff_stack = [poly for poly,scaling in queue if abs(scaling) > 1e-15]
            coeff_weights = [scaling for poly,scaling in queue if abs(scaling) > 1e-15]
            if len(coeff_weights) > 0:
                subcontribs = plan(alphas, coeff_stack)
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
    def get_overlap_gaussian_data(self,
                                  freqs_gs, modes_gs, center_gs,
                                  freqs_es, modes_es, center_es,
                                  rotation_method='duschinsky'
                                  ):

        if rotation_method == 'least-squares':
            # use least squares to express ground state modes as LCs of excited state modes
            ls_tf = np.linalg.inv(modes_es.T @ modes_es) @ modes_es.T @ modes_gs
            U, s, V = np.linalg.svd(ls_tf)
            L = U @ V  # Effectively a Duschinsky matrix
        else:
            L = modes_es.T @ modes_gs
            # L = L.T

        # raise Exception(L)


        # find displacement vector (center of gs ground state in this new basis)
        # modes_gs = modes_es @ L
        dx = modes_es.T @ (center_gs - center_es)
        c_gs = L.T @ dx

        # print(L)
        # raise Exception(c_gs)

        # print(center_es)
        # print(center_gs)
        # print(c_gs)
        #
        # raise Exception(c_gs)

        # inverse covariance matrices in this coordinate system
        Z_gs = (L.T @ np.diag(freqs_gs) @ L)
        Z_es = np.diag(freqs_es)
        Z_c = Z_gs + Z_es

        S_c = (L.T @ np.diag(1/ (freqs_gs)) @ L) + np.diag(1/(freqs_es))

        norm_gs = np.prod(np.power(freqs_gs, 1/4))
        norm_es = np.prod(np.power(freqs_es, 1/4))
        prefactor = (
                norm_gs * norm_es / np.power(np.pi, len(freqs_es) / 2)
                * np.exp(-1/2 * np.dot(dx, np.dot(np.linalg.inv(S_c), dx)))
        )
        # prefactor = 1

        alphas_c, modes_c = np.linalg.eigh(Z_c)
        center = np.linalg.inv(Z_c) @ (Z_gs @ c_gs)

        return (alphas_c, modes_c, center, prefactor), (L, c_gs), (np.eye(L.shape[0]), np.zeros(c_gs.shape))

    @classmethod
    def eval_fcf_overlaps(self,
                           excitations_gs, freqs_gs, modes_gs, center_gs,
                           excitations_es, freqs_es, modes_es, center_es
                           ):
        """
        Evaluates the Gaussian overlaps between two H.O. wave functions defined by
        a set of polynomial coefficients, broadening factors, and centers, assuming
        the modes and centers are in an Eckart fream
        """

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
        center_gs = np.asanyarray(center_gs)
        center_es = np.asanyarray(center_es)

        (alphas, modes_c, center, scaling), (L_gs, c_gs), (L_es, c_es) = self.get_overlap_gaussian_data(
            freqs_gs, modes_gs, center_gs,
            freqs_es, modes_es, center_es
        )

        # We now have the space to define the parameters that go into the overlap calculation
        shift_gs = center - c_gs #- L_gs @ modes_c.T @ center
        shift_es = center - c_es #- modes_c.T @ center

        # raise Exception(shift_es, shift_gs, c_gs)

        Q_gs = (modes_c.T @ L_gs)
        Q_es = (modes_c.T @ L_es)

        polys1 = [
            HermiteProductPolynomial.from_quanta(eg, freqs_gs).shift(shift_gs)
            for eg in excitations_gs
        ]
        polys2 = [
            HermiteProductPolynomial.from_quanta(ee, freqs_es).shift(shift_es)
            for ee in excitations_es
        ]

        Q = np.concatenate([Q_gs, Q_es], axis=1)

        overlaps = []
        for spoly_gs in polys1:
            for spoly_es in polys2:
                poly = spoly_gs.concat(spoly_es)
                ov = scaling * self.evaluate_shifted_poly_overlap(
                        poly, Q,
                        alphas
                    )
                # print("<::", scaling, ov/scaling, ov)
                overlaps.append(ov)

        return overlaps

    @classmethod
    def embed_modes(cls, gs_nms: 'NormalModes', es_nms, masses=None):
        from ..Molecools import Molecule

        if masses is None:
            masses = gs_nms.masses
        if masses is None:
            masses = es_nms.masses
        if masses is None:
            raise ValueError("need masses to reembed normal modes")

        gs = Molecule(['H'] * len(masses), gs_nms.origin, masses=masses)

        embedding = gs.get_embedding_data(es_nms.origin)
        # raise Exception(embedding)
        ref_coords = embedding.reference_data.coords
        ref_mat = np.tensordot(
            gs_nms.matrix.reshape((-1,) + ref_coords.shape), embedding.reference_data.axes,
            axes=[-1, 0]
        ).reshape(gs_nms.matrix.shape)
        gs_modes = NormalModes(
            gs_nms.basis,
            ref_mat,
            origin=ref_coords,
            freqs=gs_nms.freqs,
            masses=gs_nms.masses
        )

        emb_coords = embedding.coord_data.coords[0] @ embedding.rotations[0].T
        emb_mat = np.tensordot(
            np.tensordot(
                es_nms.matrix.reshape((-1,) + emb_coords.shape), embedding.coord_data.axes[0],
                axes=[-1, 0]
            ),
            embedding.rotations[0],
            axes=[-1, 1]
        ).reshape(es_nms.matrix.shape)
        es_modes = NormalModes(
            gs_modes.basis,
            emb_mat,
            origin=emb_coords,
            freqs=es_nms.freqs,
            masses=es_nms.masses
        )

        return gs_modes, es_modes

    @classmethod
    def mass_weight_nms(cls, nms, masses=None):
        if masses is None: masses = nms.masses

        mat = nms.matrix
        mat = np.reshape(
            np.sqrt(masses)[:, np.newaxis, np.newaxis] * mat.reshape((-1, 3, mat.shape[-1])),
            mat.shape
        )

        inv = nms.inverse
        inv = np.reshape(
            np.sqrt(masses)[np.newaxis, np.newaxis, :] * inv.reshape((inv.shape[0], 3, -1)),
            inv.shape
        )

        origin = np.sqrt(masses)[:, np.newaxis] * nms.origin

        return NormalModes(
            nms.basis,
            mat,
            inverse=inv,
            origin=origin,
            freqs=nms.freqs
        )


    @classmethod
    def get_fcfs(self,
                 gs_nms: 'NormalModes', es_nms: 'NormalModes',
                 excitations, ground_states=None,
                 embed=True, masses=None, mass_weight=True
                 ):

        if embed:
            gs_nms, es_nms = self.embed_modes(gs_nms, es_nms, masses=masses)

        if mass_weight:
            gs_nms = self.mass_weight_nms(gs_nms, masses=masses)
            es_nms = self.mass_weight_nms(es_nms, masses=masses)

        gs_freqs = gs_nms.freqs
        gs_center = gs_nms.origin.flatten()
        gs_basis = gs_nms.matrix
        ndim = len(gs_freqs)
        if ground_states is None:
            ground_states = [[0] * ndim]

        es_freqs = es_nms.freqs
        es_center = es_nms.origin.flatten()
        es_basis = es_nms.matrix
        return self.eval_fcf_overlaps(
            ground_states, gs_freqs, gs_basis, gs_center,
            excitations, es_freqs, es_basis, es_center
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

    @classmethod
    def get_1D_hermite_poly(self, n, a):
        return DensePolynomial(
            np.flip(np.asarray(sp.special.hermite(n, monic=False)))
            * np.sqrt(
                a ** np.arange(n + 1) / (2 ** (n) * np.math.factorial(n))
            )
        )