import numpy as np, scipy.sparse as sp, functools as fp, itertools as ip
from .Common import PerturbationTheoryException
from McUtils.Numputils import SparseArray

__all__ = [
    "SubHamiltonian",
    "ProductOperator"
]

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
            return a ,b
        i = tuple(pad_lens(a, b) for a ,b in zip(n ,m))
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
        return (self.mode_n, ) *len(self.funcs) + self.quanta
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
            els = [m[x[i]] for m ,i in zip(mats, uinds)]
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
        mm = {k :i for i ,k in enumerate(uinds)}
        ndim = len(uinds)
        pieces = [None] * ndim
        for f, i in zip(funcs, inds):
            n = mm[i]
            if pieces[n] is None:
                pieces[n] = f(dims[i ] +padding)
            else:
                pieces[n] = pieces[n].dot(f(dims[i ] +padding))

        # for j in np.setdiff1d(totinds, inds):
        #     pieces[j] = sp.identity(dims[j])

        if return_kron:
            mat = sp.csr_matrix(fp.reduce(sp.kron, pieces))
            sub_shape = tuple(dims[i ] +padding for i in np.unique(inds) for j in range(2))
            trans = tuple(j for i in zip(range(ndim), range(ndim, 2* ndim)) for j in i)
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
        ar = 1 / np.sqrt(2) * np.sqrt(np.arange(1, n))
        bands = [
            [ar, 1],
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

        ar = 1 / np.sqrt(2) * np.sqrt(np.arange(1, n))
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