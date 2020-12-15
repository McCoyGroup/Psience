"""
Provides the operator representations needed when building a Hamiltonian representation.
I chose to only implement direct product operators. Not sure if we'll need a 1D base class...
"""

import numpy as np, scipy.sparse as sp, os, tempfile as tf
from collections import OrderedDict
from McUtils.Numputils import SparseArray
from McUtils.Scaffolding import Logger, NullLogger

from .StateSpaces import BraKetSpace

__all__ = [
    "Operator",
    "ContractedOperator"
]

class Operator:
    """
    Provides a (usually) _lazy_ representation of an operator, which allows things like
    QQQ and pQp to be calculated block-by-block.
    Crucially, the underlying basis for the operator is assumed to be orthonormal.
    """
    def __init__(self, funcs, quanta, prod_dim=None, symmetries=None,
                 selection_rules=None,
                 zero_threshold=1.0e-14):
        """
        :param funcs: The functions use to calculate representation
        :type funcs: callable | Iterable[callable]
        :param quanta: The number of quanta to do the deepest-level calculations up to (also tells us dimension)
        :type quanta: int | Iterable[int]
        :param prod_dim: The number of functions in `funcs`, if `funcs` is a direct term generator
        :type prod_dim: int | None
        :param symmetries: Labels for the funcs where if two funcs share a label they are symmetry equivalent
        :type symmetries: Iterable[int] | None
        """
        if isinstance(quanta, int):
            quanta = [quanta]
            # funcs = [funcs]
        try:
            funcs = tuple(funcs)
            single_func = False
        except TypeError:  # can't iterate
            single_func = True
        if prod_dim is None:
            if single_func:
                funcs = (funcs,)
            prod_dim = len(funcs)
        self.fdim = prod_dim
        self.funcs = funcs
        self.sel_rules = selection_rules
        self.symmetry_inds = symmetries
        self.quanta = tuple(quanta)
        self.mode_n = len(quanta)
        self.zero_threshold = zero_threshold
        # self._tensor = None
        self._parallelizer = None
            # in the future people will be able to supply this so that they can fully control how
            # the code gets parallelized

    @property
    def ndim(self):
        return self.fdim + self.mode_n
    @property
    def shape(self):
        return (self.mode_n, ) * self.fdim + self.quanta

    def get_inner_indices(self):
        """
        Gets the n-dimensional array of ijkl (e.g.) indices that functions will map over
        Basically returns the indices of the inner-most tensor

        :return:
        :rtype:
        """
        dims = self.fdim
        if dims == 0:
            return None
        shp = (self.mode_n,) * dims
        inds = np.indices(shp, dtype=int)
        tp = np.roll(np.arange(dims + 1), -1)
        base_tensor = np.transpose(inds, tp)
        return base_tensor

    def __getitem__(self, item):
        # check to see if we were actually given a bunch of bra and ket states
        if (
                isinstance(item, tuple)
                and len(item) == 2
                and all(len(x) == self.mode_n for x in item[0])
        ):
            # expected to look like num_modes X [bra, ket] X quanta
            # so we reshape it
            bras, kets = item
            modes = [None] * self.mode_n
            for i in range(self.mode_n):
                i_bras = [x[i] for x in bras]
                i_kets = [x[i] for x in kets]
                modes[i] = (i_bras, i_kets)
            item = modes
        return self.get_elements(item)

    def __getstate__(self):
        state = self.__dict__.copy()
        state['_parallelizer'] = None
        state['_tensor'] = None
        return state

    def _get_eye_tensor(self, states):
        """

        :param states:
        :type states: BraKetSpace
        :param quants:
        :type quants:
        :return:
        :rtype:
        """

        ntotal = len(states)
        states, non_orthog = states.apply_non_orthogonality([])
        return sp.csr_matrix(
            (
                np.ones(len(non_orthog)),
                (
                    np.zeros(len(non_orthog)),
                    non_orthog
                )
            ),
            shape=(1, ntotal)
        )

    def _mat_prod_operator_terms(self, inds, funcs, states, sel_rules):
        """
        Evaluates product operator terms based on 1D representation matrices coming from funcs

        :param inds:
        :type inds:
        :param funcs: functions that generate 1D representation matrices
        :type funcs:
        :param states:
        :type states: BraKetSpace
        :param sel_rules: selection rules as 1D arrays of rules for each dimension
        :type sel_rules:
        :return:
        :rtype:
        """
        from collections import OrderedDict

        # We figure out which indices are actually unique; this gives us a way to map
        # indices onto operators
        # We assume we always have as many indices as dimensions in our operator
        # and so the more unique indices the _lower_ the dimension of the operator
        # since it acts on the same index more times
        # The indices will be mapped onto ints for storage in `pieces`

        if sel_rules is not None:
            sel_rules = states.get_sel_rules_from1d(inds, sel_rules)
            states, all_sels = states.apply_sel_rules(sel_rules)

        if len(states) == 0:
            return None, None

        bras, kets = states.state_pairs
        max_dim = max(np.max(bras), np.max(kets))
        padding = 3  # for approximation reasons we need to pad our representations...

        uinds = OrderedDict((k, None) for k in inds)
        uinds = np.array(list(uinds.keys()))

        mm = {k: i for i, k in enumerate(uinds)}
        subdim = len(uinds)
        # set up a place to hold onto the pieces we'll calculate
        pieces = [None] * subdim
        # now we construct the reps from 1D ones
        for f, i in zip(funcs, inds):
            n = mm[i]  # makes sure that we fill in in the same order as uinds
            if pieces[n] is None:
                pieces[n] = f(max_dim + padding)  # type: sp.spmatrix
            else:
                # QQ -> Q.Q & QQQ -> Q.Q.Q
                og = pieces[n]  # type: sp.spmatrix
                pieces[n] = og.dot(f(max_dim + padding))

        # now we take the requisite products of the chunks for the indices that are
        # potentially non-orthogonal
        chunk = None
        for i, j, o in zip(bras, kets, pieces):
            op = o  # type: sp.spmatrix
            blob = np.asarray(op[i, j]).squeeze()
            if chunk is None:
                if isinstance(blob, (int, np.integer, float, np.floating)):
                    chunk = np.array([blob])
                else:
                    chunk = np.asarray(blob)
            else:
                chunk = chunk * blob

        # weird stuff can happen with sp.spmatrix
        if (
                isinstance(chunk, (int, np.integer, float, np.floating))
                or (isinstance(chunk, np.ndarray) and chunk.ndim == 0)
        ):
            chunk = np.array([chunk])

        return chunk, all_sels

    def _direct_prod_operator_terms(self, inds, func, states):
        """
        Evaluates product operator terms based on 1D representation matrices,
        but using a direct function to generate them

        :param inds:
        :type inds:
        :param func: a function that will take a set of indices and term-generator
        :type func:
        :return:
        :rtype:
        """

        # next we get the term generator defined by inds
        # this will likely end up calling uinds again, but ah well, that operation is cheap
        gen = func(inds)
        try:
            g, sel_rules = gen
            r = sel_rules[0][0]
            if not isinstance(r, (int, np.integer)):
                raise TypeError("Bad selection rule")
        except (TypeError, IndexError):
            sel_rules = None
        else:
            gen = g
            states, sel_rules = states.apply_sel_rules(sel_rules)
            if len(states) == 0:
                return None, None

        chunk = gen(states)

        return chunk, sel_rules

    def _calculate_single_pop_elements(self, inds, funcs, states, sel_rules):
        """
        Calculates terms for a single product operator.
        Assumes orthogonal product bases.

        :param inds: the index in the total operator tensor (i.e. which of the funcs to use)
        :type inds:
        :param funcs: the functions to use when generating representations, must return matrices
        :type funcs:
        :param states: the states to compute terms between stored like ((s_1l, s_1r), (s_2l, s_2r), ...)
                        where `s_il` is the set of quanta for mode i for the bras and
                        `s_ir` is the set of quanta for mode i for the kets
        :type states: BraKetSpace
        :param quants: the total quanta to use when building representations
        :type quants:
        :param sel_rules: the selection rules to use when building representations
        :type sel_rules:
        :return:
        :rtype:
        """

        # logger = Logger.lookup('debug')
        # if isinstance(logger, NullLogger):
        #     logger = Logger()
        #     logger.register('debug')
        #
        # with logger.block(tag="{} - {}".format(inds, isinstance(funcs, (list, tuple)))):


        # determine how many states aren't potentially coupled by the operator
        # & then determine which of those are non-orthogonal
        nstates = len(states)
        states, non_orthog = states.apply_non_orthogonality(inds)

        # logger.log_print('nort: {n}', n=len(non_orthog))

        # if none of the states are non-orthogonal...just don't calculate anything
        if len(non_orthog) == 0:
            return sp.csr_matrix((1, nstates), dtype='float')
        else:
            if isinstance(funcs, (list, tuple)):
                # if sel_rules is not None:
                #     if all_sels is not None:
                #         non_orthog = non_orthog[all_sels,]
                chunk, all_sels = self._mat_prod_operator_terms(inds, funcs, states, sel_rules)
            else:
                chunk, all_sels = self._direct_prod_operator_terms(inds, funcs, states)

            # if chunk is None:
            #     return sp.csr_matrix((1, nstates), dtype='float')

            if all_sels is not None:
                if not (isinstance(all_sels, np.ndarray) and all_sels.dtype == np.dtype(bool)):
                    raise ValueError("bad selection rule application result {}".format(all_sels))
                # logger.log_print('sels: {ns} {sels}', ns=len(all_sels), sels=np.sum(np.where(all_sels)))
                non_orthog = non_orthog[all_sels,]

            if chunk is None:
                return sp.csr_matrix((1, nstates), dtype='float')

            # finally we make sure that everything we're working with is
            # non-zero because it'll buy us time on dot products later
            non_zero = np.where(np.abs(chunk) >= self.zero_threshold)[0]
            # logger.log_print('nz: {nz}', nz=len(non_zero))

            if len(non_zero) == 0:
                return sp.csr_matrix((1, nstates), dtype='float')
            chunk = chunk[non_zero,]
            non_orthog = non_orthog[non_zero,]

            wat = sp.csr_matrix(
                (
                    chunk,
                    (
                        np.zeros(len(non_zero)),
                        non_orthog
                    )
                ), shape=(1, nstates))

            return wat

    def _get_pop_sequential(self, inds, idx, save_to_disk=False):
        """
        Sequential method for getting elements of our product operator tensors

        :param inds:
        :type inds:
        :param idx:
        :type idx:
        :param save_to_disk: whether to save to disk or not; used by the parallelizer
        :type save_to_disk: bool
        :return:
        :rtype:
        """
        res = np.apply_along_axis(self._calculate_single_pop_elements,
                                  -1, inds, self.funcs, idx, self.sel_rules
                                  )
        if save_to_disk:
            new = sp.vstack([x for y in res.flatten() for x in y])
            with tf.NamedTemporaryFile() as tmp:
                res = tmp.name + ".npz"
            # print(res)
            # os.remove(tmp.name)
            sp.save_npz(res, new, compressed=False)

        return res

    def _get_pop_parallel(self, inds, idx, parallelizer=None):
        if parallelizer != "multiprocessing":
            raise NotImplementedError("More parallelization methods are coming--just not yet")

        # in the future this will be like Parallizer.num_cores
        # and we can just have Parallelizer.apply_along_axis(func, ...)

        import multiprocessing as mp
        cores = mp.cpu_count()

        if self._parallelizer is None:
            self._parallelizer = mp.Pool()

        ind_chunks = np.array_split(inds, min(cores, inds.shape[0]))

        chunks = [(sub_arr, idx, True) for sub_arr in ind_chunks]
        res = self._parallelizer.starmap(self._get_pop_sequential, chunks)
        if isinstance(res[0], str):
            try:
                sparrays = np.array([sp.load_npz(f) for f in res])
            finally: # gotta clean up after myself
                try:
                    for f in res:
                        os.remove(f)
                except Exception as e:
                    # print(e)
                    pass
            res = sparrays

        # raise Exception(res)

        return res

    def filter_symmetric_indices(self, inds):
        """
        Determines which inds are symmetry unique.
        For something like `qqq` all permutations are equivalent, but for `pqp` we have `pi qj pj` distinct from `pj qj pi`.
        This means for `qqq` we have `(1, 0, 0) == (0, 1, 0)` but for `pqp` we only have stuff like `(2, 0, 1) == (1, 0, 2)` .

        :param inds: indices to filter symmetric bits out of
        :type inds: np.ndarray
        :return: symmetric indices & inverse map
        :rtype:
        """


        symm_labels = self.symmetry_inds
        flat = np.reshape(inds, (-1, inds.shape[-1]))
        if symm_labels is None:
            return flat, None

        ind_grps = OrderedDict()
        for i, l in enumerate(symm_labels):
            if l in ind_grps:
                ind_grps[l].append(i)
            else:
                ind_grps[l] = [i]
        ind_grps = list(ind_grps.values())

        if len(ind_grps) == len(symm_labels):
            # short-circuit case if no symmetry
            return flat, None
        elif len(ind_grps) == 1:
            # totally symmetric so we can just sort and delete the dupes
            flat = np.sort(flat, axis=1)
            return np.unique(flat, axis=0, return_inverse=True)
        else:
            # if indices are shared between groups we can't do anything about them
            # so we figure out where there are overlaps in indices
            # we do this iteratively by looping over groups, figuring out pairwise if they
            # share elements, and using an `or` op to update our list of indep terms
            indep_grps = np.full((len(flat),), False, dtype=bool)
            # this is 4D, but over relatively small numbers of inds (maybe 6x6x6x6 at max)
            # and so hopefully still fast, esp. relative to computing actual terms
            for i in range(len(ind_grps)):
                for j in range(i+1, len(ind_grps)):
                    for t1 in ind_grps[i]:
                        for t2 in ind_grps[j]:
                            indep_grps = np.logical_or(
                                indep_grps,
                                flat[:, t1] == flat[:, t2]
                            )
            indep_grps = np.logical_not(indep_grps)
            symmable = flat[indep_grps]

            # pull out the groups of indices and sort them
            symmetriz_grps = [np.sort(symmable[:, g], axis=1) for g in ind_grps]
            # stack them together to build a full array
            # then reshuffle to get the OG ordering
            reordering = np.argsort(np.concatenate(ind_grps))
            symmetrized = np.concatenate(symmetriz_grps, axis=1)[:, reordering]
            flat[indep_grps] = symmetrized

            return np.unique(flat, axis=0, return_inverse=True)

    def get_elements(self, idx, parallelizer=None):#'multiprocessing'):
        """
        Calculates a subset of elements

        :param idx: bra and ket states as tuples of elements
        :type idx: Iterable[(Iterable[int], Iterable[int])]
        :return:
        :rtype:
        """

        # we expect this to be an iterable object that looks like
        # num_modes X [bra, ket] X quanta
        if not isinstance(idx, BraKetSpace):
            idx = BraKetSpace.from_indices(idx, quanta=self.quanta)

        inds = self.get_inner_indices()

        if inds is None: # just a number
            new = self._get_eye_tensor(idx)
        else:
            mapped_inds, inverse = self.filter_symmetric_indices(inds)

            if parallelizer is not None:
                res = self._get_pop_parallel(mapped_inds, idx,
                                             parallelizer=parallelizer
                                             )
            else:
                res = self._get_pop_sequential(mapped_inds, idx)

            if inverse is not None:
                res = res.flatten()[inverse]

            wat = [x for y in res for x in y]
            new = sp.vstack(wat)

        shp = inds.shape if inds is not None else ()
        res = SparseArray(new)
        res = res.reshape(shp[:-1] + res.shape[-1:])

        return res

    def __repr__(self):
        return "{}(<{}>, {})".format(
            type(self).__name__,
            ", ".join(str(s) for s in self.shape),
            self.funcs
        )

class ContractedOperator(Operator):
    """
    Provides support for terms that look like `pGp` or `p(dG/dQ)Qp` by
    expanding them out as the pure operator component that depends on the basis states (i.e. `pp` or `pQp`)
    and doing the appropriate tensor contractions with the expansion coefficients (i.e. `G` or `dG/dQ`)
    """

    def __init__(self, coeffs, funcs, quanta, prod_dim=None, axes=None, symmetries=None,
                 selection_rules=None,
                 zero_threshold=1.0e-14
                 ):
        """
        :param coeffs: The tensor of coefficients contract with the operator representation (`0` means no term)
        :type coeffs: np.ndarray | int
        :param funcs: The functions use to calculate representation
        :type funcs: callable | Iterable[callable]
        :param quanta: The number of quanta to do the deepest-level calculations up to
        :type quanta: int | Iterable[int]
        :param axes: The axes to use when doing the contractions
        :type axes: Iterable[int] | None
        :param symmetries: The symmetries to pass through to `Operator`
        :type symmetries: Iterable[int] | None
        """
        self.coeffs = coeffs
        self.axes = axes
        super().__init__(funcs, quanta, symmetries=symmetries, prod_dim=prod_dim,
                         selection_rules=selection_rules,
                         zero_threshold=zero_threshold
                         )

    def get_elements(self, idx, parallelizer=None):
        """
        Computes the operator values over the specified indices

        :param idx: which elements of H0 to compute
        :type idx: Iterable[int]
        :return:
        :rtype:
        """

        c = self.coeffs
        if not isinstance(c, (int, np.integer, float, np.floating)):
            # takes an (e.g.) 5-dimensional SparseTensor and turns it into a contracted 2D one
            axes = self.axes
            if axes is None:
                axes = (tuple(range(c.ndim)), )*2
            subTensor = super().get_elements(idx, parallelizer=parallelizer)
            if isinstance(subTensor, np.ndarray):
                contracted = np.tensordot(subTensor.squeeze(), c, axes=axes)
            else:
                contracted = subTensor.tensordot(c, axes=axes).squeeze()
        elif c == 0:
            contracted = 0  # a short-circuit
        else:
            subTensor = super().get_elements(idx)
            if c == 1:
                return subTensor
            contracted = c * subTensor

        return contracted

    def __repr__(self):
        return "{}(opdim=<{}>, cdim=<{}>, {})".format(
            type(self).__name__,
            ", ".join(str(s) for s in self.shape),
            None if not hasattr(self.coeffs, 'shape') else ", ".join(str(s) for s in self.shape),
            self.funcs
        )

