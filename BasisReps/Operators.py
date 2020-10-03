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
        nstates = len(states[0])
        uinds = np.unique(inds)
        missing = [i for i in range(ndim) if i not in uinds]
        equivs = [states[i][0] == states[i][1] for i in missing]
        orthog = np.prod(equivs, axis=0).astype(int) # taking advantage of the fact that bools act like ints
        non_orthog = np.where(orthog != 0)[0]

        # if none of the states are non-orthogonal...just don't calculate anything
        if len(non_orthog) == 0:
            return sp.csr_matrix([1, nstates], dtype='float')
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

            chunk = None
            for s, o in zip(non_orthog_states, pieces):
                op = o #type: sp.spmatrix
                blob = np.asarray(op[s]).squeeze()
                if chunk is None:
                    chunk = np.asarray(blob)
                else:
                    chunk *= blob

            return sp.csr_matrix(
                (
                    np.asarray(chunk),
                    (
                        np.zeros(len(non_orthog)),
                        non_orthog
                    )
                ), shape=(1, len(orthog)))

    def _get_pop_sequential(self, inds, idx, save_to_disk=False):
        """
        Sequential method for getting elements of our product operator tensors

        :param inds:
        :type inds:
        :param idx:
        :type idx:
        :param save_io: whether to save to disk or not
        :type save_io:
        :return:
        :rtype:
        """
        tensor = np.apply_along_axis(self._calculate_single_pop_elements, -1, inds, self.funcs, idx, self.quanta)

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


        return np.apply_along_axis(self._calculate_single_pop_elements, -1, inds, self.funcs, idx, self.quanta)
    def get_elements(self, idx, parallelizer=None):
        """
        Calculates a subset of elements

        :param idx: bra and ket states as tuples of elements
        :type idx: Iterable[Iterable[int]]
        :return:
        :rtype:
        """
        inds = self.get_inner_indices()
        idx = tuple(tuple(np.array([i]) if isinstance(i, (int, np.integer)) else i for i in j) for j in idx)
        shp = inds.shape

        if parallelizer is None:
            # parallelizer = self._parallelizer
            parallelizer = 'multiprocessing'

        chunks = [(sub_arr, idx, True)
                  for sub_arr in ind_chunks]
        # if len(inds[0]) == 4:
        #     raise Exception(chunks)
        res = self._parallelizer.starmap(self._pull_sequential, chunks)
        res =

        new = sp.vstack([x for y in res.flatten() for x in y])

        res = SparseArray(new)
        res = res.reshape(shp[:-1] + res.shape[-1:])
        print("=======>", res)

        return res

    def _pull_parallel(self,
                          inds, tens, idx, quants,
                          parallel_method = "multiprocessing"
                          ):
        """
        Allows us to take the puller function and evaluate it in parallel.
        The hope is that once we've got a more general `Parallelizer` architecture
        up and running, we'll be able to distribute over more and more cores.
        """
        # adapted from https://stackoverflow.com/a/45555516/5720002

        if parallel_method != "multiprocessing":
            raise NotImplementedError("More parallelization methods are coming--just not yet")

        # in the future this will be like Parallizer.num_cores
        # and we can just have Parallelizer.apply_along_axis(func, ...)

        import multiprocessing as mp
        cores = mp.cpu_count()

        if self._parallelizer is None:
            self._parallelizer = mp.Pool()

        ind_chunks = np.array_split(inds, min(cores, inds.shape[0]))

        # save tensors to disk and let pull_sequential reload them on the other side...
        tens_files = np.full(tens.shape, None, dtype=object)
        # if len(inds[0]) == 4:
        #     raise Exception(len(idx), len(idx[0]), [len(x) for x in idx[0]], idx[0][0].shape)
        try:
            flit = tens.flat
            i = flit.coords
            for y in flit:
                if not isinstance(y, SparseArray):
                        y = SparseArray(y)
                with tf.NamedTemporaryFile() as tmp:
                    res = tmp.name + ".npz"
                tens_files[i] = y.savez(res, compressed=False)
                i = flit.coords

            chunks = [(sub_arr, tens_files, idx, quants, True)
                      for sub_arr in ind_chunks]
            # if len(inds[0]) == 4:
            #     raise Exception(chunks)
            res = self._parallelizer.starmap(self._pull_sequential, chunks)
        finally:
            # clean up all these tensors after ourselves
            for x in tens_files.flat:
                try:
                    os.remove(x)
                except:# Exception as e: # sometimes shit happens but we stil need to clean up...
                    # print(e)
                    pass

        if not isinstance(res[0], str):
            res = np.concatenate(res, axis=0)

        return res # convenient list-join
    def _pull_sequential(self, inds, tens, idx, quants, save_io=False):
        # save_io is so we don't have to pass huge arrays through pickle
        if tens.dtype == np.dtype(object):
            # to avoid passing shit through the pipe we saved stuff on the other end
            tens_files = tens
            tens = np.full(tens_files.shape, None, dtype=object)
            flit = tens_files.flat
            i = flit.coords
            for y in flit:
                tens[i] = SparseArray.loadz(y)
                i = flit.coords
        # raise Exception(tens, tens.dtype, np.dtype(object), tens.dtype == np.dtype(object))
        res = np.apply_along_axis(self._take_subtensor, -1, inds, tens, idx, quants)
        if save_io:
            if isinstance(res[0], sp.spmatrix):
                new = sp.vstack([x for y in res for x in y])
                with tf.NamedTemporaryFile() as tmp:
                    res = tmp.name + ".npz"
                # os.remove(tmp.name)
                sp.save_npz(res, new, compressed=False)

        return res
    def get_elements_old(self, idx, method='parallel'):
        if len(idx) != len(self.quanta):
            raise ValueError("number of indices requested must be the same as the number of quanta")
        inds = self.get_inner_indices()
        idx = tuple(
            tuple(np.array([i]) if isinstance(i, (int, np.integer)) else i for i in j)
            for j in idx
        )
        tens = self.tensor
        quants = self.quanta

        flat_inds = inds.reshape((-1, inds.shape[-1]))

        if method == 'sequential':
            res = self._pull_sequential(flat_inds, tens, idx, quants)
        else:
            res = self._pull_parallel(flat_inds, tens, idx, quants)

        test = res[0]
        shp = inds.shape
        if isinstance(test, sp.spmatrix):
            new = sp.vstack([x for y in res.flatten() for x in y])
            res = SparseArray(new)
        elif isinstance(test, str):
            # print(res)
            # this means we wrote out some temporary npz files to use IO as a data-transfer mechanism
            try:
                sp_arrs = [ sp.load_npz(f) for f in res ]
            finally:
                for f in res:
                    os.remove(f)
            new = sp.vstack(sp_arrs)
            res = SparseArray(new)
        else:
            # raise Exception(test)
            res = SparseArray(np.concatenate(res, axis=0).squeeze())

        return res

    @staticmethod
    def _take_subtensor(inds, t, x, qn):
        """
        Takes the subtensor of `t` defined by `inds` given a total set of indices `x`
        Then applies orthonormality conditions, i.e. _this assumes an orthonormal basis_
        """

        # compute orthonormality indices
        missing = [i for i in range(len(x)) if i not in inds]
        equivs = [x[i][0] == x[i][1] for i in missing]
        orthog = np.prod(equivs, axis=0).astype(float)

        if isinstance(orthog, np.ndarray):
            sly = t[tuple(inds)]
            # finds the appropriate indices of t to sample
            uinds = np.unique(inds)
            non_zero_orthog = np.where(orthog != 0)[0]
            nz_sub = tuple(tuple(j[non_zero_orthog]) for i in uinds for j in x[i])
            non_zero_vals = sly[nz_sub][0] # we get a spurious `1` in the shape from NumPy or SparseArray...
            res = sp.csr_matrix(
                (
                    non_zero_vals,
                    (
                        np.zeros(len(non_zero_vals)),
                        non_zero_orthog
                    )
                ),
                shape=(1, len(orthog))
            )
        elif orthog != 0.:
            sly = t[tuple(inds)]
            # finds the appropriate indices of t to sample
            uinds = np.unique(inds)
            # non_zero_orthog = np.where(orthog != 0)[0]
            sub = tuple(tuple(j) for i in uinds for j in x[i])
            res = sly[sub] * orthog
        else:
            # raise Exception(orthog)
            res = np.zeros(len(x[0][0]))

        return res

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

    _op_mat_cache = {} # we try to cache equivalent things
    def _operator_submatrix(self, funcs, dims, inds, padding=3, return_kron=True):
        """
        Returns the operator submatrix for a product operator like piQjpk or whatever...
        Used in the initial implementation that directly constructed the entire product operator tensor

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

        # we apply caching so that we only compute the symmetry distinct ones
        if return_kron and self.symmetry_inds is not None:
            if funcs not in self._op_mat_cache:
                symm_inds = self.symmetry_inds
                # symm_groups = [np.where(symm_inds==i)[0] for i in np.unique(symm_inds)]
                self._op_mat_cache[funcs] = {}
                self._op_mat_cache[funcs]["symm_inds"] = np.array(symm_inds)
                # self._op_mat_cache[funcs]["symm_groups"] = symm_groups
            f_cache = self._op_mat_cache[funcs]
            if dims not in f_cache:
                f_cache[dims] = {}
            symm_cache = f_cache[dims]
            symm_inds = f_cache["symm_inds"]
            # determine which of our distinct indices go with which symmetry-distinct terms
            ind_groups = {
                u: np.array([i for i, k in enumerate(inds) if k == u]) for u in uinds
            }
            symm_groups = {
                u: symm_inds[i] for u,i in ind_groups.items()
            }
            # now we sort _this_, since these bits can be transposed at will
            symm_key = tuple(sorted(
                tuple((u, tuple(v)) for u,v in symm_groups.items()),
                key=lambda v:v[0]
            ))
            if symm_key in symm_cache:
                return symm_cache[symm_key]

        # we figure out how many unique indices we have so that we can figure out our object dimension
        uinds = np.unique(inds)
        mm = {k:i for i, k in enumerate(uinds)}
        ndim = len(uinds)
        pieces = [None] * ndim
        for f, i in zip(funcs, inds):
            n = mm[i]
            if pieces[n] is None:
                pieces[n] = f(dims[i] + padding)
            else:
                pieces[n] = pieces[n].dot(f(dims[i] + padding))

        if return_kron:
            # if we want to return the Kronecker product, we build it using sparse methods
            # this is also the only branch of this for which we do caching...
            mat = sp.csr_matrix(fp.reduce(sp.kron, pieces))
            sub_shape = tuple(dims[i] + padding for i in np.unique(inds) for j in range(2))
            trans = tuple(j for i in zip(range(ndim), range(ndim, 2* ndim)) for j in i)
            mat = SparseArray(mat, shape=sub_shape).transpose(trans)
            if self.symmetry_inds is not None:
                symm_cache[symm_key] = mat
        else:
            mat = pieces

        return mat

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

    def get_elements(self, idx):
        """
        Computes the operator values over the specified indices

        :param idx: which elements of H0 to compute
        :type idx: Iterable[int]
        :return:
        :rtype:
        """

        c = self.coeffs
        if not isinstance(c, int):
            # takes an (e.g.) 5-dimensional SparseTensor and turns it into a contracted 2D one
            axes = self.axes
            if axes is None:
                axes = (tuple(range(c.ndim)), )*2
            subTensor = super().get_elements(idx)
            if isinstance(subTensor, np.ndarray):
                contracted = np.tensordot(subTensor.squeeze(), c, axes=axes)
            else:
                contracted = subTensor.tensordot(c, axes=axes).squeeze()
        else:
            contracted = 0 # a short-circuit

        return contracted

