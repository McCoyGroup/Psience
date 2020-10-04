"""
Provides the operator representations needed when building a Hamiltonian representation.
I chose to only implement direct product operators. Not sure if we'll need a 1D base class...
"""

import numpy as np, scipy.sparse as sp, functools as fp, os, tempfile as tf
from McUtils.Numputils import SparseArray

#TODO: abstract the VPT2 stuff out so we can use it for a general operator too

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
    def __init__(self, funcs, quanta, symmetries=None):
        """
        :param funcs: The functions use to calculate representation
        :type funcs: callable | Iterable[callable]
        :param quanta: The number of quanta to do the deepest-level calculations up to
        :type quanta: int | Iterable[int]
        :param symmetry_inds: Labels for the funcs where if two funcs share a label they are symmetry equivalent
        :type symmetry_inds: Iterable[int] | None
        """
        if isinstance(quanta, int):
            quanta = [quanta]
            # funcs = [funcs]
        self.funcs = tuple(funcs)
        self.symmetry_inds = symmetries
        self.quanta = tuple(quanta)
        self.mode_n = len(quanta)
        # self._tensor = None
        self._parallelizer = None
            # in the future people will be able to supply this so that they can fully control how
            # the code gets parallelized

    @property
    def ndim(self):
        return len(self.funcs) + len(self.quanta)
    @property
    def shape(self):
        return (self.mode_n, ) *len(self.funcs) + self.quanta
    # @property
    # def tensor(self):
    #     if self._tensor is None:
    #         self._tensor = self.product_operator_tensor()
    #     return self._tensor

    def get_inner_indices(self):
        """
        Gets the n-dimensional array of ijkl (e.g.) indices that functions will map over
        Basically returns the indices of the inner-most tensor

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

    def _calculate_single_pop_elements(self, inds, funcs, states, quants):
        """
        Calculates terms for a single product operator.
        Assumes orthogonal product bases.

        :param inds: the index in the total operator tensor (i.e. which of the funcs to use)
        :type inds:
        :param funcs: the functions to use when generating representations
        :type funcs:
        :param states: the states to compute terms between stored like ((s_1l, s_1r), (s_2l, s_2r), ...)
                        where `s_il` is the set of quanta for mode i for the bras and
                        `s_ir` is the set of quanta for mode i for the kets
        :type states: Iterable[Iterable[Iterable[int]]]
        :param quants: the total quanta to use when building representations
        :type quants:
        :return:
        :rtype:
        """
        # determine how many states are properly non-orthogonal
        dims = quants
        ndim = len(quants)
        nstates = len(states[0][0])
        uinds = np.unique(inds)
        missing = [i for i in range(ndim) if i not in uinds]
        equivs = [states[i][0] == states[i][1] for i in missing]
        orthog = np.prod(equivs, axis=0).astype(int) # taking advantage of the fact that bools act like ints
        single_state = isinstance(orthog, (int, np.integer)) # means we got a single state to calculate over
        if single_state:
            orthog = np.array([orthog])

        non_orthog = np.where(orthog != 0)[0]

        # if none of the states are non-orthogonal...just don't calculate anything
        if len(non_orthog) == 0:
            return sp.csr_matrix((1, nstates), dtype='float')
        else:
            non_orthog_states = [(states[i][0][non_orthog], states[i][1][non_orthog]) for i in uinds]
            # otherwise we calculate only the necessary terms (after building the 1D reps)
            padding = 3 # for approximation reasons we need to pad our representations, actually...
            # we figure out which indices are actually unique
            uinds = np.unique(inds)
            mm = {k: i for i, k in enumerate(uinds)}
            # then set up a place to hold onto the pieces we'll calculate
            subdim = len(uinds)
            pieces = [None] * subdim
            # TODO: if it turns out that we only need a small number of non-orthogonal terms
            #       we should directly compute terms for cheapness
            for f, i in zip(funcs, inds):
                n = mm[i] # makes sure that we fill in in the same order as uinds
                if pieces[n] is None:
                    pieces[n] = f(dims[i] + padding) #type: sp.spmatrix
                else:
                    # QQ -> Q.Q & QQQ -> Q.Q.Q
                    og = pieces[n] #type: sp.spmatrix
                    pieces[n] = og.dot(f(dims[i] + padding))

            # now we take the requisite products of the chunks
            chunk = None
            for s, o in zip(non_orthog_states, pieces):
                op = o #type: sp.spmatrix
                blob = np.asarray(op[s]).squeeze()
                if chunk is None:
                    if isinstance(blob, (int, np.integer, float, np.floating)):
                        chunk = np.array([blob])
                    else:
                        chunk = np.asarray(blob)
                else:
                    chunk *= blob

            # all sorts of weird shit can happen with sp.spmatrix
            if (
                    isinstance(chunk, (int, np.integer, float, np.floating))
                    or isinstance(chunk, np.ndarray) and chunk.shape==()
            ):
                chunk = np.array([chunk])

            non_zero = np.where(chunk != 0.)[0] # by filtering again we save on computation in the dot products

            if len(non_zero) == 0:
                return sp.csr_matrix((1, nstates), dtype='float')

            non_orthog = non_orthog[non_zero,]
            chunk = chunk[non_zero,]

            # try:
            wat = sp.csr_matrix(
                (
                    chunk,
                    (
                        np.zeros(len(non_orthog)),
                        non_orthog
                    )
                ), shape=(1, nstates))

            # if wat.shape[1] == 2 or wat.nnz == 2:
            #     raise Exception(wat, nstates, non_orthog, inds)

            return wat
            # except:
            #     raise Exception([
            #         (
            #             non_zero,
            #             chunk,
            #             (
            #                 np.zeros(len(non_orthog)),
            #                 non_orthog
            #             )
            #         )
            #     ])

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
        res = np.apply_along_axis(self._calculate_single_pop_elements, -1, inds, self.funcs, idx, self.quanta)
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
    def get_elements(self, idx, parallelizer=None):#parallelizer=None):#'multiprocessing'):
        """
        Calculates a subset of elements

        :param idx: bra and ket states as tuples of elements
        :type idx: Iterable[(Iterable[int], Iterable[int])]
        :return:
        :rtype:
        """
        inds = self.get_inner_indices()

        # we expect this to be an iterable object that looks like
        # num_modes X [bra, ket] X quanta
        idx = tuple(
            tuple(np.array([i]) if isinstance(i, (int, np.integer)) else np.asarray(i) for i in j)
            for j in idx
        )

        # raise Exception(idx)
        if len(idx) != self.mode_n:
            raise ValueError("BraKet spec {} isn't valid for Operator with dimension {}".format(
                idx,
                (self.ndim, self.mode_n)
            ))
        shp = inds.shape

        if parallelizer is not None:
            # parallelizer = self._parallelizer
            inds = inds.reshape((-1, inds.shape[-1]))
            res = self._get_pop_parallel(inds, idx, parallelizer=parallelizer)
        else:
            res = self._get_pop_sequential(inds, idx)

        wat = [x for y in res.flatten() for x in y]
        # try:
        new = sp.vstack(wat)
        # except:
        #     raise Exception(wat, res[0, 1, -1].toarray(), parallelizer)

        res = SparseArray(new)
        res = res.reshape(shp[:-1] + res.shape[-1:])
        # print("=======>", res)

        # print(self.symmetry_inds)
        # if tuple(self.symmetry_inds) == (2, 2, 2):
        #     raise Exception(
        #         list(
        #             zip(
        #                 np.array([
        #                     [s[0][(0, 1),] for s in idx],
        #                     [s[1][(0, 1),] for s in idx]
        #                 ]).T,
        #                 res[:, :, :, (0, 1)].toarray().transpose(3, 0, 1, 2)
        #             )
        #         )
        #     )
        return res


class DebugOperator:
    """
    Old implementation, here for debug purposes
    """
    def __init__(self, funcs, quanta, symmetries=None):
        """
        :param funcs: The functions use to calculate representation
        :type funcs: callable | Iterable[callable]
        :param quanta: The number of quanta to do the deepest-level calculations up to
        :type quanta: int | Iterable[int]
        """
        if isinstance(quanta, int):
            quanta = [quanta]
            # funcs = [funcs]
        self.funcs = funcs
        self.quanta = tuple(quanta)
        self.mode_n = len(quanta)
        self._tensor = None

    @property
    def ndim(self):
        return len(self.funcs) + len(self.quanta)
    @property
    def shape(self):
        return (self.mode_n, ) *len(self.funcs) + self.quanta
    @property
    def tensor(self):
        if self._tensor is None:
            self._tensor = self.product_operator_tensor()
        return self._tensor

    def get_inner_indices(self):
        """
        Gets the n-dimensional array of ijkl (e.g.) indices that functions will map over
        Basically returns the indices of the inner-most tensor

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
    def get_individual_elements(self, idx):
        """
        TBH I can't remember what this function is supposed to do ?_?
        :param idx:
        :type idx:
        :return:
        :rtype:
        """
        if len(idx) != len(self.quanta):
            raise ValueError("number of indices requested must be the same as the number of modes")
        inds = self.get_inner_indices()
        idx = tuple(tuple(np.array([i]) if isinstance(i, (int, np.integer)) else i for i in j) for j in idx)
        funcs = self.funcs
        quants = self.quanta
        def pull(inds, f=funcs, x=idx, qn = quants):
            uinds = np.unique(inds)
            mats = self._operator_submatrix(f, qn, inds, return_kron=False)
            els = [m[x[i]] for m ,i in zip(mats, uinds)]
            if isinstance(els[0], np.matrix):
                els = [np.asarray(e).squeeze() for e in els]
            res = np.prod(els, axis=0)

            return res
        res = np.apply_along_axis(pull, -1, inds)
        return res
    def get_elements(self, idx, parallelizer=None):
        if len(idx) != len(self.quanta):
            raise ValueError("number of indices requested must be the same as the number of quanta")
        inds = self.get_inner_indices()
        idx = tuple(tuple(np.array([i]) if isinstance(i, (int, np.integer)) else i for i in j) for j in idx)
        tens = self.tensor
        quants = self.quanta

        pull = lambda inds, t=tens,x=idx,qn=quants,f=self._take_subtensor: f(inds, t, x, qn)
        res = np.apply_along_axis(pull, -1, inds)
        return SparseArray(res.squeeze())
    @staticmethod
    def _take_subtensor(inds, t, x, qn):
        """
        Takes the subtensor of `t` defined by `inds` given a total set of indices `x`
        Then applies orthonormality conditions, i.e. _this assumes an orthonormal basis_
        """
        # finds the appropriate indices of t to sample
        sly = t[tuple(inds)]
        uinds = np.unique(inds)
        sub = tuple(tuple(j) for i in uinds for j in x[i])
        res = sly[sub]

        # compute orthonormality indices
        missing = [i for i in range(len(x)) if i not in inds]
        equivs = [x[i][0] == x[i][1] for i in missing]
        orthog = np.prod(equivs, axis=0).astype(int)

        return res * orthog

    def product_operator_tensor(self):
        """
        Generates the tensor created from the product of funcs over the dimensions dims,
        Note that this isn't a totally legit tensor since it's ragged

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

    def _operator_submatrix(cls, funcs, dims, inds, padding = 3, return_kron = True):
        """
        Returns the operator submatrix for a product operator like piQjpk or whatever

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
        mm = {k:i for i ,k in enumerate(uinds)}
        ndim = len(uinds)
        pieces = [None] * ndim
        for f, i in zip(funcs, inds):
            n = mm[i]
            if pieces[n] is None:
                pieces[n] = f(dims[i] +padding)
            else:
                pieces[n] = pieces[n].dot(f(dims[i] +padding))

        if return_kron:
            mat = sp.csr_matrix(fp.reduce(sp.kron, pieces))
            sub_shape = tuple(dims[i ] +padding for i in np.unique(inds) for j in range(2))
            trans = tuple(j for i in zip(range(ndim), range(ndim, 2* ndim)) for j in i)
            mat = SparseArray(mat, shape=sub_shape).transpose(trans)
        else:
            mat = pieces
        return mat

# Operator = DebugOperator

class ContractedOperator(Operator):
    """
    Provides support for terms that look like `pGp` or `p(dG/dQ)Qp` by
    expanding them out as the pure operator component that depends on the basis states (i.e. `pp` or `pQp`)
    and doing the appropriate tensor contractions with the expansion coefficients (i.e. `G` or `dG/dQ`)
    """

    def __init__(self, coeffs, funcs, quanta, axes=None, symmetries=None):
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
        super().__init__(funcs, quanta, symmetries=symmetries)

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

