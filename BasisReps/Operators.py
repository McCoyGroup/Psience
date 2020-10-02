"""
Provides the operator representations needed when building a Hamiltonian representation.
I chose to only implement direct product operators. Not sure if we'll need a 1D base class...
"""

import numpy as np, scipy.sparse as sp, functools as fp
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
    def get_elements(self, idx):
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
            res = np.zeros(len(orthog))
            non_zero_vals = sly[nz_sub]
            res[non_zero_orthog] = non_zero_vals
        elif orthog != 0.:
            sly = t[tuple(inds)]
            # finds the appropriate indices of t to sample
            uinds = np.unique(inds)
            # non_zero_orthog = np.where(orthog != 0)[0]
            sub = tuple(tuple(j) for i in uinds for j in x[i])
            res = sly[sub] * orthog
        else:
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

