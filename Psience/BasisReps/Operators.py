"""
Provides the operator representations needed when building a Hamiltonian representation.
I chose to only implement direct product operators. Not sure if we'll need a 1D base class...
"""

import tracemalloc

import numpy as np, scipy.sparse as sp, os, tempfile as tf, time, gc
from collections import OrderedDict
from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Scaffolding import Logger, NullLogger
from McUtils.Parallelizers import Parallelizer, SerialNonParallelizer

from .StateSpaces import BasisStateSpace, SelectionRuleStateSpace, PermutationallyReducedStateSpace, BraKetSpace

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
                 selection_rule_steps=None, # for product operators
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
        self.selection_rules = selection_rules
        self.selection_rule_steps = selection_rule_steps
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

    def clear_cache(self):
        """
        :return:
        :rtype:
        """

        funcs = self.funcs
        if not isinstance(funcs, tuple):
            funcs = (funcs,)
        for f in funcs:
            if hasattr(f, 'clear_cache'):
                f.clear_cache()

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

    def get_inner_indices(self, reduced_inds=False):
        """
        Gets the n-dimensional array of ijkl (e.g.) indices that functions will map over
        Basically returns the indices of the inner-most tensor

        :return:
        :rtype:
        """
        dims = self.fdim
        if dims == 0:
            return None
        n_terms = max(len(x) for x in self.selection_rules) if reduced_inds else self.mode_n
        shp = (n_terms,) * dims
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

        # next pull the unique indices
        uinds_basic, inverse = np.unique(flat, axis=0, return_inverse=True)

        # Finally apply a lexical sorting so that we
        # can do as little work as possible later since work can be shared between (0,0,0) and (0,0,1) but not
        # (1, 1, 0)

        indsort = np.lexsort(uinds_basic.T)
        # raise Exception(indsort.shape, uinds_basic.shape)
        uinds = uinds_basic[indsort]
        # if len(uinds) > 500:
        #     raise Exception(inds, uinds)

        return uinds, (inverse, np.argsort(indsort))

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
        else:
            all_sels = None

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
            sel_rules = None

        # bras, kets = states.state_pairs
        # states = [bk for i, bk in enumerate(zip(bras, kets)) if i in inds]
        chunk = gen(zip(*states.state_pairs))

        return chunk, sel_rules

    def _calculate_single_pop_elements(self, inds, funcs, states, sel_rules,
                                       check_orthogonality=True,
                                       use_sel_rule_filtering=False,
                                       use_sel_sum_filtering=False
                                       ):
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
        if check_orthogonality:
            # TODO: selection rules are actually _cheaper_ to apply than this in general, esp. if we focus
            #       only on the number of quanta that can change within the set of indices
            #       so we should support applying them first & then only doing this for the rest
            if use_sel_rule_filtering and sel_rules is not None:
                _, idx = np.unique(inds, return_index=True)
                uinds = inds[np.sort(idx)]
                states, non_orthog = states.apply_sel_rules_along(sel_rules, uinds)
                if len(non_orthog) > 0:
                    states, non_orthog_2 = states.apply_non_orthogonality(inds)#, max_inds=self.fdim)
                    non_orthog = non_orthog[non_orthog_2]
            elif use_sel_sum_filtering and sel_rules is not None:
                states, non_orthog = states.apply_sel_sums(sel_rules, inds)
                if len(non_orthog) > 0:
                    states, non_orthog_2 = states.apply_non_orthogonality(inds)  # , max_inds=self.fdim)
                    non_orthog = non_orthog[non_orthog_2]
            else:
                states, non_orthog = states.apply_non_orthogonality(inds)#, max_inds=self.fdim)
        else:
            non_orthog = np.arange(nstates)


        def return_empty():
            return sp.csr_matrix((1, nstates), dtype='float')

        # if none of the states are non-orthogonal...just don't calculate anything
        if len(non_orthog) == 0:
            return return_empty()
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
                return return_empty()

            # if check_orthogonality:
            # finally we make sure that everything we're working with is
            # non-zero because it'll buy us time on dot products later
            non_zero = np.where(np.abs(chunk) >= self.zero_threshold)[0]
            if len(non_zero) == 0:
                return return_empty()
            chunk = chunk[non_zero,]
            non_orthog = non_orthog[non_zero,]

            wat = sp.csr_matrix(
                (
                    chunk,
                    (
                        np.zeros(len(non_zero), dtype='int8'),
                        non_orthog
                    )
                ), shape=(1, nstates))
            # else:
            #     wat = chunk

            return wat

    def _get_pop_sequential(self, inds, idx, check_orthogonality=True, save_to_disk=False):
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
                                  -1, inds, self.funcs, idx,
                                  self.selection_rules,
                                  check_orthogonality
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

    def _get_pop_parallel(self, inds, idx, parallelizer, check_orthogonality=True):
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

        inds = parallelizer.scatter(inds)
        # print(f"Num inds={len(idx)}")
        parallelizer.print("evaluating over inds of shape {s}", s=inds.shape, log_level=parallelizer.logger.LogLevel.Debug)
        start = time.time()
        dump = self._get_pop_sequential(inds, idx, save_to_disk=False, check_orthogonality=check_orthogonality)
        end = time.time()
        parallelizer.print("took {t:.3f}s", t=end-start, log_level=parallelizer.logger.LogLevel.Debug)
        res = parallelizer.gather(dump)
        parallelizer.print("gather took {t:.3f}s", t=time.time()-end, log_level=parallelizer.logger.LogLevel.Debug)
        if isinstance(res, np.ndarray):
            parallelizer.print("got res of shape {s}", s=res.shape, log_level=parallelizer.logger.LogLevel.Debug)
        if parallelizer.on_main:
            res = self._load_arrays(res)

        return res

    @Parallelizer.main_restricted
    # @profile
    def _main_get_elements(self, inds, idx,
                           parallelizer=None,
                           check_orthogonality=True,
                           ):
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
        try:
            idx = idx.obj # shared proxy deref
        except AttributeError:
            pass

        if not isinstance(idx, BraKetSpace):
            self.logger.log_print("constructing BraKet space")
            idx = BraKetSpace.from_indices(idx, quanta=self.quanta)

        if inds is None: # just a number
            self.logger.log_print("evaluating identity tensor over {} elements".format(len(idx)))
            new = self._get_eye_tensor(idx)
        else:
            mapped_inds, inv_dat = self.filter_symmetric_indices(inds)
            if parallelizer is not None and not isinstance(parallelizer, SerialNonParallelizer):
                self.logger.log_print("evaluating {nel} elements over {nind} unique indices using {cores} processes".format(
                    nel=len(idx),
                    nind=len(mapped_inds),
                    cores=parallelizer.nproc
                ))
                res = self._get_pop_parallel(mapped_inds, idx, parallelizer, check_orthogonality=check_orthogonality)
            else:
                self.logger.log_print("evaluating {nel} elements over {nind} unique indices sequentially".format(
                    nel = len(idx),
                    nind = len(mapped_inds)
                ))
                res = self._get_pop_sequential(mapped_inds, idx, check_orthogonality=check_orthogonality)

            if inv_dat is not None:
                inverse, preinv = inv_dat
                res = res.flatten()[preinv][inverse]

            wat = [x for y in res for x in y]
            new = sp.vstack(wat)

        shp = inds.shape if inds is not None else ()
        res = SparseArray.from_data(new, cache_block_data=False)
        res = res.reshape(shp[:-1] + res.shape[-1:])

        return res

    @Parallelizer.worker_restricted
    def _worker_get_elements(self, idx, parallelizer=None, check_orthogonality=True):
        """

        :param idx:
        :type idx:
        :param parallelizer:
        :type parallelizer:
        :return:
        :rtype:
        """
        try:
            idx = idx.obj # shared proxy deref
        except AttributeError:
            pass
        # worker threads only need to do a small portion of the work
        # and actually inherit most of their info from the parent process
        self._get_pop_parallel(None, idx, parallelizer, check_orthogonality=True)

    def _get_elements(self, inds, idx, full_inds=None, parallelizer=None, check_orthogonality=True):
        """
        Runs the regular element getting algorithm in parallel

        :param idx:
        :type idx:
        :param parallelizer:
        :type parallelizer:
        :return:
        :rtype:
        """
        if full_inds is not None:
            inds = full_inds
        self._worker_get_elements(idx, parallelizer=parallelizer, check_orthogonality=check_orthogonality)
        return self._main_get_elements(inds, idx, parallelizer=parallelizer, check_orthogonality=check_orthogonality)

    def _split_idx(self, idx):
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
        return idx_splits

    # from memory_profiler import profile
    # @profile
    def _eval_chunks(self, inds, idx_splits,
                     parallelizer=None,
                     check_orthogonality=True,
                     memory_constrained=False
                     ):

        # import tracemalloc
        #
        # tracemalloc.start()

        chunks = []

        if parallelizer is None:
            parallelizer = self.parallelizer
        if parallelizer is not None and not isinstance(parallelizer, SerialNonParallelizer):
            with parallelizer:

                for n, idx in enumerate(idx_splits):
                    idx: BraKetSpace
                    # parallelizer.printer = self.logger.log_print
                    max_nproc = len(inds.flatten())
                    # idx.load_space_diffs()
                    idx = parallelizer.share(idx)
                    elem_chunk = parallelizer.run(self._get_elements, None, idx,
                                                  check_orthogonality=check_orthogonality,
                                                  main_kwargs={'full_inds': inds},
                                                  comm=None if max_nproc >= parallelizer.nprocs else list(range(max_nproc))
                                                  )
                    if memory_constrained:
                        idx = idx.obj
                        if len(idx_splits) > 1:
                            idx.free()
                        else:
                            idx.clear_cache()
                        idx_splits[n] = None
                        del idx
                        gc.collect()

                    chunks.append(elem_chunk)
        else:

            for n, idx in enumerate(idx_splits):
                idx: BraKetSpace
                elem_chunk = self._main_get_elements(inds, idx,
                                                     parallelizer=None,
                                                     check_orthogonality=check_orthogonality
                                                     )
                if memory_constrained:
                    # self.logger.log_print("mem usage post\n{m}", m="\n".join(str(x) for x in tracemalloc.take_snapshot().statistics('lineno')))
                    if len(idx_splits) > 1:
                        idx.free()
                    else:
                        idx.clear_cache()
                    idx_splits[n] = None
                    del idx
                    gc.collect()

                    # self.logger.log_print("mem usage freed\n{m}", m="\n".join(str(x) for x in tracemalloc.take_snapshot().statistics('lineno')))

                chunks.append(elem_chunk)

        return chunks

    # @profile
    def _construct_elem_array(self, chunks):

        if all(isinstance(x, np.ndarray) for x in chunks):
            elems = np.concatenate(chunks, axis=0)
        else:
            elems = chunks[0].concatenate(*chunks[1:])

        return elems

    # @profile
    def get_elements(self, idx,
                     parallelizer=None,
                     check_orthogonality=True,
                     memory_constrained=False
                     # need to work in further cache clearing to account for mem. constraints
                     ):
        """
        Calculates a subset of elements

        :param idx: bra and ket states as tuples of elements
        :type idx: BraKetSpace
        :return:
        :rtype:
        """

        idx_splits = self._split_idx(idx)

        inds = self.get_inner_indices()
        # print(">>>>>", inds.size)
        if inds is None:
            return self._get_elements(inds, idx, check_orthogonality=check_orthogonality)

        chunks = self._eval_chunks(inds, idx_splits,
                                   parallelizer=parallelizer,
                                   check_orthogonality=check_orthogonality,
                                   memory_constrained=memory_constrained
                                   )

        arr = self._construct_elem_array(chunks)
        return arr

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

    def get_transformed_space(self, base_space, rules=None, parallelizer=None, logger=None, **opts):
        """
        Returns the space one would get from applying
        the selection rules from this operator

        :param base_space:
        :type base_space: BasisStateSpace
        :return:
        :rtype: SelectionRuleStateSpace
        """
        if rules is None:
            rules = self.selection_rules
        if parallelizer is None:
            parallelizer = self.parallelizer
        return base_space.apply_selection_rules(rules, parallelizer=parallelizer, logger=logger, **opts)

    def _calculate_single_transf(self, inds, funcs, base_space, sel_rules):
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

        # once we're doing this over all possible indices I think I've lost any benefit

        # TODO: edit direct sum generator to be a true generator where I can directly resume
        #       the calculation over a different set of target dimensions if I need to or something
        #       although I'm low-key not sure that's feasible
        #       _Alternately_ allow for splitting by which indices were touched rather than what the
        #       input state was

        # trying to not calculate anything unnecessary
        uinds = np.unique(inds)
        sr = [r for r in sel_rules if len(r) == len(uinds) or (len(uinds) == 1 and len(r) == 0)]
        states = base_space.apply_selection_rules(sr, target_dimensions=uinds)
        brakets = states.get_representation_brakets()
        vals = self._calculate_single_pop_elements(inds, funcs, brakets, sel_rules, check_orthogonality=False)

        return vals, brakets#, states
    def _get_transf_sequential(self, inds, base_space):#, save_to_disk=False):
        """
        Sequential method for getting elements of our product operator tensors

        :param inds:
        :type inds:
        :param base_space:
        :type base_space:
        :param save_to_disk: whether to save to disk or not; used by the parallelizer
        :type save_to_disk: bool
        :return:
        :rtype:
        """

        res = np.apply_along_axis(self._calculate_single_transf,
                                  -1, inds, self.funcs, base_space, self.selection_rules
                                  )

        # spaces = [r[2] for r in res]
        brakets = [r[1] for r in res]
        res = [r[0] for r in res]

        # total_space = res[0][2]
        # total_brakets = res[0][1]
        # for s in res[1:]:
        #     total_brakets = total_brakets.concatenate(s[1]) # to preserve ordering since union will destroy that
        #     total_space = total_space.union(s[2])
        # res = np.concatenate([r[0] for r in res])

        # if save_to_disk:
        #     flattened = [x for y in res.flatten() for x in y]
        #     res = []
        #     try:
        #         for array in flattened:
        #             with tf.NamedTemporaryFile() as tmp:
        #                 dump = tmp.name + ".npz"
        #             # print(res)
        #             # os.remove(tmp.name)
        #             sp.save_npz(dump, array, compressed=False)
        #             res.append(dump)
        #     except:
        #         for file in res:
        #             try:
        #                 os.remove(file)
        #             except:
        #                 pass
        #         raise
        #
        # raise Exception(res, total_space, total_brakets, total_brakets.bras.excitations[5], inds, total_brakets.kets.excitations[5])

        return res, brakets
    @Parallelizer.main_restricted
    def _apply_transformations_sequential(self, inds, base_space, perm_class_map, parallelizer=None):
        """
        Implementation of get_elements to be run on the main process...

        :param idx:
        :type idx:
        :param parallelizer:
        :type parallelizer:
        :return:
        :rtype:
        """

        raise NotImplementedError("no")

        if inds is None:  # just a number
            self.logger.log_print("returning identity tensor")
            new = self._get_eye_tensor(base_space)
        else:
            mapped_inds, inverse = self.filter_symmetric_indices(inds) # I can't actually reduce here...
            inds = np.reshape(inds, (-1, inds.shape[-1]))
            if parallelizer is not None and not isinstance(parallelizer, SerialNonParallelizer):
                raise NotImplementedError("still working up parallelism")
                self.logger.log_print(
                    "evaluating {nel} elements over {nind} unique indices using {cores} processes".format(
                        nel=len(idx),
                        nind=len(mapped_inds),
                        cores=parallelizer.nproc
                    ))
                res = self._get_pop_parallel(mapped_inds, idx, parallelizer)
            else:
                self.logger.log_print("transforming space of size {nel} over {nind} unique indices sequentially".format(
                    nel=len(base_space),
                    nind=len(mapped_inds)
                ))

                res, brakets = self._get_transf_sequential(mapped_inds, base_space)

            if inverse is not None: # fix the reductions by symmetry
                res = [res[i] for i in inverse]
                brakets = [brakets[i] for i in inverse]

            # TODO: this will need to be _a lot_ faster if this is going to become the default
            #       method for getting representations...
            # now we need to apply the permutations which
            # we'll do in a slightly dumb/inefficient way
            # by determining which brakets correspond to which
            # input classes and using this to do the necessary
            # permutations of all relevant components
            full_res = []
            full_bras = []
            full_kets = []
            full_inds = []
            for r, b, i in zip(res, brakets, inds):
                res_groups, sorting = nput.group_by(r, b.bras.excitations)
                ket_groups, _ = nput.group_by(b.kets.excitations, b.bras.excitations, sorting=sorting)
                res_keys, res_vals = res_groups
                # raise Exception(res_groups[0], ket_groups[0])
                for key, vals, kets in zip(res_keys, res_vals, ket_groups[1]):
                    for test_key, perms, inv_perms in perm_class_map:
                        if np.all(key == test_key):
                            # permute kets
                            perm_kets = np.array([k[perms] for k in kets]).transpose((1, 0, 2))
                            perm_bras = np.broadcast_to(key[perms][:, np.newaxis, :], (len(perms), len(vals), key.shape[-1]))
                            perm_inds = np.broadcast_to(inv_perms[:, np.newaxis, i], (len(perms), len(vals), len(i)))
                            perm_vals = np.broadcast_to(vals[np.newaxis], (len(perms), len(vals)))
                            full_res.append(perm_vals.reshape(-1))
                            full_bras.append(perm_bras.reshape(-1, perm_bras.shape[-1]))
                            full_kets.append(perm_kets.reshape(-1, perm_kets.shape[-1]))
                            full_inds.append(perm_inds.reshape(-1, perm_inds.shape[-1]))
                            break
                    else:
                        raise ValueError("couldn't find key {} in input keys?".format(key))

            full_bras = np.concatenate(full_bras, axis=0)
            full_kets = np.concatenate(full_kets, axis=0)
            full_inds = np.concatenate(full_inds)
            full_res = np.concatenate(full_res)

            new_brakets = BraKetSpace(
                BasisStateSpace(base_space.basis, full_bras, mode=BasisStateSpace.StateSpaceSpec.Excitations),
                BasisStateSpace(base_space.basis, full_kets, mode=BasisStateSpace.StateSpaceSpec.Excitations),
            )
            _, _, bk_inds = nput.unique(np.concatenate([full_bras, full_kets], axis=1), return_inverse=True)
            real_full_inds = np.concatenate([full_inds, bk_inds[:, np.newaxis]], axis=1).T

            res = SparseArray.from_data(
                (
                    full_res,
                    real_full_inds
                ),
                shape=(self.mode_n, self.mode_n, np.max(bk_inds)+1)
            )

            # raise Exception("oookay")
            #
            # res = np.concatenate(res)
            # full_brakets = brakets[0]
            # for s in brakets[1:]:
            #     full_brakets = full_brakets.concatenate(s)

        return res, new_brakets

    def apply_reduced(self, base_space, parallelizer=None, logger=None):
        """
        Takes a base space as input and applies the held selection rules in semi-efficient
        fashion only on the indices that can change and then uses this to compute all matrix
        elements, returning then the final generated space

        :param base_space:
        :type base_space: BasisStateSpace | PermutationallyReducedStateSpace
        :return:
        :rtype: tuple[SparseArray, BraKetSpace]
        """

        if not isinstance(base_space, PermutationallyReducedStateSpace):
            base_space = base_space.permutationally_reduce()

        rep_space = base_space.representative_space()
        if self.chunk_size is not None and len(rep_space) > self.chunk_size:
            raise NotImplementedError("not handling chunking yet...")
            idx_splits = rep_space.split(self.chunk_size)
        else:
            idx_splits = [rep_space]

        inds = self.get_inner_indices(reduced_inds=False)
        perm_class_map = [(e, p, np.argsort(p, axis=1)) for e,p in zip(rep_space.excitations, base_space.equivalence_classes)]
        chunks = []
        brakets = []
        for idx in idx_splits:
            if inds is None:
                return self._apply_transformations_sequential(inds, idx)

            if parallelizer is None:
                parallelizer = self.parallelizer

            if parallelizer is not None and not isinstance(parallelizer, SerialNonParallelizer):
                raise NotImplementedError("Don't have parallelism yet")
                parallelizer.printer = self.logger.log_print
                elem_chunk = parallelizer.run(self._get_elements, None, idx,
                                              main_kwargs={'full_inds':inds},
                                              comm=None if len(inds) >= parallelizer.nprocs else list(range(len(inds)))
                                              )
            else:
                elem_chunk, elem_brakets = self._apply_transformations_sequential(inds, idx, perm_class_map, parallelizer=None)

            brakets.append(elem_brakets)
            chunks.append(elem_chunk)

        brakets = brakets[0]
        chunks = chunks[0]

        # raise Exception(chunks, brakets)

        # if all(isinstance(x, np.ndarray) for x in chunks):
        #     elems = np.concatenate(chunks, axis=0)
        # else:
        #     from functools import reduce
        #     elems = reduce(lambda a, b: a.concatenate(b), chunks[1:], chunks[0])

        return chunks, brakets

class ContractedOperator(Operator):
    """
    Provides support for terms that look like `pGp` or `p(dG/dQ)Qp` by
    expanding them out as the pure operator component that depends on the basis states (i.e. `pp` or `pQp`)
    and doing the appropriate tensor contractions with the expansion coefficients (i.e. `G` or `dG/dQ`)
    """

    def __init__(self, coeffs, funcs, quanta,
                 prod_dim=None, axes=None, symmetries=None,
                 selection_rules=None,
                 selection_rule_steps=None,
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
                         selection_rule_steps=selection_rule_steps,
                         zero_threshold=zero_threshold,
                         parallelizer=parallelizer,
                         logger=logger
                         )
        self.chunk_size = chunk_size

    def _get_element_block(self, idx, parallelizer=None, check_orthogonality=True, memory_constrained=False):
        c = self.coeffs
        if not isinstance(c, (int, np.integer, float, np.floating)):
            # takes an (e.g.) 5-dimensional SparseTensor and turns it into a contracted 2D one
            axes = self.axes
            if axes is None:
                axes = (tuple(range(c.ndim)),) * 2
            subTensor = super().get_elements(idx, parallelizer=parallelizer, check_orthogonality=check_orthogonality, memory_constrained=memory_constrained)

            # we collect here to minimize the effect of memory spikes if possible
            # self.clear_cache()
            # SparseArray.clear_cache()
            gc.collect()
            if isinstance(subTensor, np.ndarray):
                if len(axes[1]) > 0:
                    contracted = np.tensordot(subTensor.squeeze(), c, axes=axes)
                else:
                    # TODO: make this broadcasting more robust
                    contracted = c[np.newaxis, :] * subTensor.squeeze()[:, np.newaxis]
            else:
                with subTensor.cache_options(enabled=False):
                    if len(axes[1]) > 0:
                        contracted = subTensor.tensordot(c, axes=axes).squeeze()
                    else:
                        # TODO: make this broadcasting more robust
                        subTensor = subTensor.expand_dims(-1)
                        contracted = subTensor * c[np.newaxis, :]

            # if self.fdim > 3:
            #     raise RuntimeError("wwwwooooof")

        elif c == 0:
            contracted = 0  # a short-circuit
        else:
            subTensor = super().get_elements(idx, check_orthogonality=check_orthogonality)
            if c == 1:
                return subTensor
            contracted = c * subTensor

        gc.collect()

        return contracted

    def get_elements(self, idx, parallelizer=None, check_orthogonality=True, memory_constrained=False):
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

        idx_splits = self._split_idx(idx)

        chunks = [self._get_element_block(idx, check_orthogonality=check_orthogonality, memory_constrained=memory_constrained) for idx in idx_splits]
        if all(isinstance(x, np.ndarray) for x in chunks):
            subchunks = [np.array([y]) if y.shape == () else y for y in chunks]
            contracted = np.concatenate(subchunks, axis=0)
        else:
            contracted = chunks[0].concatenate(*chunks[1:])
            # print([x.shape for x in chunks])
            # print(contracted.shape)
            #
            # raise Exception(contracted.asarray() -
            #                 np.concatenate([x.asarray() for x in chunks], axis=0)
            #                 )
            # assert np.all(contracted == np.concatenate([x.asarray() for x in chunks], axis=0))
        # print(contracted.shape)
        return contracted

    def apply_reduced(self, base_space, parallelizer=None, logger=None):

        c = self.coeffs
        if isinstance(c, (int, np.integer, float, np.floating)) and c == 0:
            return 0

        c = self.coeffs
        if not isinstance(c, (int, np.integer, float, np.floating)):
            # takes an (e.g.) 5-dimensional SparseTensor and turns it into a contracted 2D one
            axes = self.axes
            if axes is None:
                axes = (tuple(range(c.ndim)),) * 2

            subTensor, brakets = super().apply_reduced(base_space, parallelizer=parallelizer, logger=logger)

            # we collect here to minimize the effect of memory spikes if possible
            # self.clear_cache()
            # SparseArray.clear_cache()
            gc.collect()
            if isinstance(subTensor, np.ndarray):
                contracted = np.tensordot(subTensor.squeeze(), c, axes=axes)
            else:
                contracted = subTensor.tensordot(c, axes=axes).squeeze()

            # if self.fdim > 3:
            #     raise RuntimeError("wwwwooooof")

        elif c == 0:
            contracted = 0  # a short-circuit
            brakets = None
        else:
            subTensor, brakets = super().apply_reduced(base_space, parallelizer=parallelizer, logger=logger)
            if c == 1:
                return subTensor
            contracted = c * subTensor

        gc.collect()

        return contracted, brakets

        # return chunks, brakets


    def __repr__(self):
        return "{}(<{}>x<{}>, {})".format(
            type(self).__name__,
            None if not hasattr(self.coeffs, 'shape') else self._get_dim_string(self.coeffs.shape),
            self._get_dim_string(self.shape),
            self.funcs
        )

