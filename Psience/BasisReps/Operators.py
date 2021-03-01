"""
Provides the operator representations needed when building a Hamiltonian representation.
I chose to only implement direct product operators. Not sure if we'll need a 1D base class...
"""

import numpy as np, scipy.sparse as sp, os, tempfile as tf, time, gc
from collections import OrderedDict
from McUtils.Numputils import SparseArray
from McUtils.Scaffolding import Logger, NullLogger
from McUtils.Parallelizers import Parallelizer, SerialNonParallelizer

from .StateSpaces import BraKetSpace

__all__ = [
    "Operator",
    "ContractedOperator"
]

__reload_hook__ = [ '.StateSpaces' ]

class Operator:
    """
    Provides a (usually) _lazy_ representation of an operator, which allows things like
    QQQ and pQp to be calculated block-by-block.
    Crucially, the underlying basis for the operator is assumed to be orthonormal.
    """
    def __init__(self, funcs, quanta,
                 prod_dim=None,
                 symmetries=None,
                 selection_rules=None,
                 parallelizer=None,
                 logger=None,
                 zero_threshold=1.0e-14,
                 chunk_size=None
                 ):
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
        if logger is None:
            logger = NullLogger()
        self.logger = logger
        self._parallelizer = parallelizer
        self.chunk_size = chunk_size

    @property
    def ndim(self):
        return self.fdim + self.mode_n
    @property
    def shape(self):
        return (self.mode_n, ) * self.fdim + self.quanta

    def load_parallelizer(self):
        """
        Loads the held parallelizer if needed

        :return:
        :rtype:
        """
        from McUtils.Parallelizers import Parallelizer

        if not isinstance(self._parallelizer, Parallelizer):
            return Parallelizer.lookup(self._parallelizer)
        else:
            return self._parallelizer

    @property
    def parallelizer(self):
        """
        Loads a parallelizer that can be used to speed up various bits of the calculation

        :return:
        :rtype:
        """
        if self._parallelizer is not None:
            return self.load_parallelizer()
        else:
            return self._parallelizer

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
        if isinstance(item, BraKetSpace): # common case, so we use it to skip further fuckery
            pass
        elif (
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
        state['logger'] = None
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

        uinds = OrderedDict((k, None) for k in inds)
        uinds = np.array(list(uinds.keys()))
        states = states.take_subdimensions(uinds)

        if sel_rules is not None:
            sel_rules = states.get_sel_rules_from1d(inds, sel_rules)
            states, all_sels = states.apply_sel_rules(sel_rules)

        if len(states) == 0:
            return None, None

        bras, kets = states.state_pairs
        max_dim = max(np.max(bras), np.max(kets))
        padding = 3  # for approximation reasons we need to pad our representations...

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
        Evaluates product operator terms using a function to directly generate them

        :param inds:
        :type inds:
        :param func: a function that will take a set of indices and term-generator
        :type func:
        :return:
        :rtype:
        """

        # next we get the term generator defined by inds
        # this will likely end up calling uinds again, but ah well, that operation is cheap
        uinds = OrderedDict((k, None) for k in inds)
        uinds = np.array(list(uinds.keys()))
        # the issue here is that _this_ operation is not necessarily cheap...
        # the subdimension states still can be quite expensive to get indices from
        states = states.take_subdimensions(uinds)

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

        # bras, kets = states.state_pairs
        # states = [bk for i, bk in enumerate(zip(bras, kets)) if i in inds]
        chunk = gen(zip(*states.state_pairs))

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

        # determine how many states aren't potentially coupled by the operator
        # & then determine which of those are non-orthogonal
        nstates = len(states)
        states, non_orthog = states.apply_non_orthogonality(inds)

        # if none of the states are non-orthogonal...just don't calculate anything
        if len(non_orthog) == 0:
            return sp.csr_matrix((1, nstates), dtype='float')
        else:
            if isinstance(funcs, (list, tuple)):
                chunk, all_sels = self._mat_prod_operator_terms(inds, funcs, states, sel_rules)
            else:
                chunk, all_sels = self._direct_prod_operator_terms(inds, funcs, states)

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
            flattened = [x for y in res.flatten() for x in y]
            res = []
            try:
                for array in flattened:
                    with tf.NamedTemporaryFile() as tmp:
                        dump = tmp.name + ".npz"
                    # print(res)
                    # os.remove(tmp.name)
                    sp.save_npz(dump, array, compressed=False)
                    res.append(dump)
            except:
                for file in res:
                    try:
                        os.remove(file)
                    except:
                        pass
                raise

        return res

    @Parallelizer.main_restricted
    def _load_arrays(self, res, parallelizer=None):
        if isinstance(res[0], str):
            res = sum(res, [])

            start = time.time()
            try:
                sparrays = [sp.load_npz(f) for f in res]
            finally: # gotta clean up after myself
                try:
                    for f in res:
                        os.remove(f)
                except Exception as e:
                    # print(e)
                    pass
            self.logger.log_print("reloading arrays from disk: took {e:.3f}s", e=(time.time() - start))
            res = np.array(sparrays, dtype=object)

        return res

    def _get_pop_parallel(self, inds, idx, parallelizer):
        """
        :param inds:
        :type inds:
        :param idx:
        :type idx:
        :param parallelizer:
        :type parallelizer: Parallelizer
        :return:
        :rtype:
        """

        # parallelizer.print("calling...")
        inds = parallelizer.scatter(inds)
        # parallelizer.print(inds.shape)
        # parallelizer.print(inds.shape)
        dump = self._get_pop_sequential(inds, idx, False)
        res = parallelizer.gather(dump)
        # parallelizer.print(res)
        if parallelizer.on_main:
            res = self._load_arrays(res)

        return res

    @Parallelizer.main_restricted
    def _main_get_elements(self, inds, idx, parallelizer=None):
        """
        Implementation of get_elements to be run on the main process...

        :param idx:
        :type idx:
        :param parallelizer:
        :type parallelizer:
        :return:
        :rtype:
        """

        # we expect this to be an iterable object that looks like
        # num_modes X [bra, ket] X quanta
        if not isinstance(idx, BraKetSpace):
            self.logger.log_print("constructing BraKet space")
            idx = BraKetSpace.from_indices(idx, quanta=self.quanta)

        if inds is None: # just a number
            self.logger.log_print("returning identity tensor")
            new = self._get_eye_tensor(idx)
        else:
            mapped_inds, inverse = self.filter_symmetric_indices(inds)
            if parallelizer is not None and not isinstance(parallelizer, SerialNonParallelizer):
                self.logger.log_print("evaluating {nel} elements over {nind} unique indices using {cores} processes".format(
                    nel = len(idx),
                    nind = len(mapped_inds),
                    cores = parallelizer.nproc
                ))
                res = self._get_pop_parallel(mapped_inds, idx, parallelizer)
            else:
                self.logger.log_print("evaluating {nel} elements over {nind} unique indices sequentially".format(
                    nel = len(idx),
                    nind = len(mapped_inds)
                ))
                res = self._get_pop_sequential(mapped_inds, idx)

            if inverse is not None:
                res = res.flatten()[inverse]

            wat = [x for y in res for x in y]
            new = sp.vstack(wat)

        shp = inds.shape if inds is not None else ()
        res = SparseArray.from_data(new)
        res = res.reshape(shp[:-1] + res.shape[-1:])

        return res

    @Parallelizer.worker_restricted
    def _worker_get_elements(self, idx, parallelizer=None):
        """

        :param idx:
        :type idx:
        :param parallelizer:
        :type parallelizer:
        :return:
        :rtype:
        """
        # worker threads only need to do a small portion of the work
        # and actually inherit most of their info from the parent process
        self._get_pop_parallel(None, idx, parallelizer)

    def _get_elements(self, inds, idx, parallelizer=None):
        """
        Runs the regular element getting algorithm in parallel

        :param idx:
        :type idx:
        :param parallelizer:
        :type parallelizer:
        :return:
        :rtype:
        """
        self._worker_get_elements(idx, parallelizer=parallelizer)
        return self._main_get_elements(inds, idx, parallelizer=parallelizer)

    def get_elements(self, idx, parallelizer=None):
        """
        Calculates a subset of elements

        :param idx: bra and ket states as tuples of elements
        :type idx: Iterable[(Iterable[int], Iterable[int])]
        :return:
        :rtype:
        """

        if self.chunk_size is not None:
            if isinstance(idx, BraKetSpace):
                if len(idx) > self.chunk_size:
                    idx_splits = idx.split(self.chunk_size)
                else:
                    idx_splits = [idx]
            else:
                idx_splits = [idx]
        else:
            idx_splits = [idx]

        chunks = []
        for idx in idx_splits:

            inds = self.get_inner_indices()
            if inds is None:
                return self._get_elements(inds, idx)

            if parallelizer is None:
                parallelizer = self.parallelizer

            if parallelizer is not None and not isinstance(parallelizer, SerialNonParallelizer):
                parallelizer.printer = self.logger.log_print
                elem_chunk = parallelizer.run(self._get_elements, inds, idx)
            else:
                elem_chunk = self._main_get_elements(inds, idx, parallelizer=None)

            chunks.append(elem_chunk)

        if all(isinstance(x, np.ndarray) for x in chunks):
            elems = np.concatenate(chunks, axis=0)
        else:
            from functools import reduce
            elems = reduce(lambda a, b: a.concatenate(b), chunks[1:], chunks[0])

        return elems

    @staticmethod
    def _get_dim_string(dims):
        if len(dims) < 7:
            return ", ".join(str(s) for s in dims)
        else:
            return "{}, ({} more), {}".format(
                ", ".join(str(s) for s in dims[:2]),
                len(dims) - 4,
                ", ".join(str(s) for s in dims[-2:])
            )
    def __repr__(self):
        return "{}(<{}>, {})".format(
            type(self).__name__,
            self._get_dim_string(self.shape),
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
                 zero_threshold=1.0e-14,
                 chunk_size=None,
                 parallelizer=None, logger=None
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
        :param prod_dim:
        :type prod_dim:
        :param selection_rules:
        :type selection_rules:
        :param parallelizer:
        :type parallelizer:
        :param logger:
        :type logger:
        :param zero_threshold:
        :type zero_threshold:
        :param chunk_size: number of elements that can be evaluated at once (for memory reasons)
        :type chunk_size: int | None
        """
        self.coeffs = coeffs
        self.axes = axes
        super().__init__(funcs, quanta, symmetries=symmetries, prod_dim=prod_dim,
                         selection_rules=selection_rules,
                         zero_threshold=zero_threshold,
                         parallelizer=parallelizer,
                         logger=logger
                         )
        self.chunk_size = chunk_size

    def _get_element_block(self, idx, parallelizer=None):
        c = self.coeffs
        if not isinstance(c, (int, np.integer, float, np.floating)):
            # takes an (e.g.) 5-dimensional SparseTensor and turns it into a contracted 2D one
            axes = self.axes
            if axes is None:
                axes = (tuple(range(c.ndim)),) * 2
            subTensor = super().get_elements(idx, parallelizer=parallelizer)

            if isinstance(subTensor, np.ndarray):
                contracted = np.tensordot(subTensor.squeeze(), c, axes=axes)
            else:
                contracted = subTensor.tensordot(c, axes=axes).squeeze()

            # if self.fdim > 3:
            #     raise RuntimeError("wwwwooooof")

        elif c == 0:
            contracted = 0  # a short-circuit
        else:
            subTensor = super().get_elements(idx)
            if c == 1:
                return subTensor
            contracted = c * subTensor

        gc.collect()

        return contracted

    def get_elements(self, idx, parallelizer=None):
        """
        Computes the operator values over the specified indices

        :param idx: which elements of H0 to compute
        :type idx: Iterable[int]
        :return:
        :rtype:
        """
        c = self.coeffs
        if isinstance(c, (int, np.integer, float, np.floating)) and c == 0:
            return 0

        if self.chunk_size is not None:
            if isinstance(idx, BraKetSpace):
                if len(idx) > self.chunk_size:
                    idx_splits = idx.split(self.chunk_size)
                else:
                    idx_splits = [idx]
            else:
                idx_splits = [idx]
        else:
            idx_splits = [idx]

        chunks = [self._get_element_block(idx) for idx in idx_splits]
        if all(isinstance(x, np.ndarray) for x in chunks):
            contracted = np.concatenate(chunks, axis=0)
        else:
            from functools import reduce
            contracted = reduce(lambda a,b: a.concatenate(b), chunks[1:], chunks[0])
        return contracted

    def __repr__(self):
        return "{}(<{}>x<{}>, {})".format(
            type(self).__name__,
            None if not hasattr(self.coeffs, 'shape') else self._get_dim_string(self.coeffs.shape),
            self._get_dim_string(self.shape),
            self.funcs
        )

