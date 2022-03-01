"""
Provides an abstract Hamiltonian object that can be used when building representations
"""

__all__ = [
    "Representation",
    "ExpansionRepresentation"
]

import tracemalloc

import numpy as np, itertools as ip, scipy.sparse as sp, time, gc, tempfile as tf

from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Scaffolding import Logger, NullLogger

from .Operators import Operator
from .StateSpaces import BraKetSpace, SelectionRuleStateSpace, StateSpaceMatrix

# for interactive work
__reload_hook__ = [ '.StateSpaces', ".Operators" ]

class Representation:
    """
    A `Representation` provides a simple interface to build matrix representations of operators expressed
    in high-dimensional spaces.

    """

    def __init__(self, compute, basis, name=None, logger=None,
                 selection_rules=None,
                 selection_rule_steps=None,
                 memory_constrained=False
                 ):
        """
        :param compute: the function that turns indices into values
        :type compute: callable | Operator
        :param basis: the basis quanta used in the representations
        :type basis: RepresentationBasis
        :param logger: logger for printing out debug info
        :type logger: None | Logger
        """
        if isinstance(compute, Operator):
            operator = compute
            compute = self._compute_op_els
        else:
            operator = None
        self.operator = operator
        self._compute = compute
        self.basis = basis
        self.dims = basis.dimensions
        self._diminds = None
        if logger is None:
            logger = NullLogger()
        self.logger = logger
        self.array = StateSpaceMatrix(basis)
        self._selection_rules = selection_rules
        self._selection_rule_steps = selection_rule_steps
        self.name=name
        self.memory_constrained = memory_constrained

    @property
    def parallelizer(self):
        if self.operator is not None:
            return self.operator.parallelizer
        else:
            return None

    def compute(self, inds, **kwargs):
        # if isinstance(inds, BraKetSpace):
        #     # allows us to use cached stuff
        #     return self.array.compute_values(self._compute, inds)
        # else:
        return self._compute(inds, **kwargs)
    def compute_cached(self, inds):
        if isinstance(inds, BraKetSpace):
            # allows us to use cached stuff
            self.array._compute_uncached_values(self._compute, inds)
        else:
            raise ValueError("Can only compute cached values when given explicit BraKets")
    @property
    def chunk_size(self):
        if self.operator is None:
            return None
        else:
            return self.operator.chunk_size

    def clear_cache(self):
        # print(">>>>>>>>>>>>> wat", self.compute, self.operator)
        if self.operator is not None:
            # print("?????")
            self.operator.clear_cache()
        elif hasattr(self.compute, 'clear_cache'):
            self.compute.clear_cache()
        gc.collect()
    def _compute_op_els(self, inds, check_orthogonality=True, memory_constrained=None):
        res = self.operator.get_elements(inds,
                                          check_orthogonality=check_orthogonality,
                                          memory_constrained=self.memory_constrained if memory_constrained is None else memory_constrained
                                          )
        return res

    @property
    def diag(self):
        return self[self.dim_inds, self.dim_inds]
    @property
    def ndims(self):
        return int(np.prod(self.dims))
    @property
    def dim_inds(self):
        # we probably want to avoid using this...
        if self._diminds is None:
            self._diminds = np.arange(self.ndims)
        return self._diminds
    def _get_inds(self, n):
        if isinstance(n, slice):
            start = n.start
            stop = n.stop
            step = n.step
            if start is None:
                start = 0
            elif start < 0:
                start = self.ndims + start
            if stop is None:
                stop = self.ndims
            elif stop < 0:
                stop = self.ndims + stop
            if step is None:
                step = 1
            return np.arange(start, stop, step)
        elif isinstance(n, np.ndarray):
            negs = np.where(n < 0)
            if len(negs[0]) > 0:
                n[negs] += self.ndims
            return n
    # @profile
    def get_brakets(self, states, check_orthogonality=True, memory_constrained=False):
        """
        Computes term elements based on getting a BraKetSpace.
        Can directly pass element specs through, since the shape management shit
        is managed by the BraKetSpace

        :param states:
        :type states: BraKetSpace | Tuple[np.ndarray, np.ndarray]
        :return:
        :rtype:
        """
        if not isinstance(states, BraKetSpace):
            states = BraKetSpace.from_indices(states, basis=self.basis)

        self.logger.log_print(
            "evaluating in BraKet space {states}",
            states=states
        )

        return self.compute(states, check_orthogonality=check_orthogonality, memory_constrained=memory_constrained)

    def get_element(self, n, m):
        """
        Computes term elements.
        Determines first whether it needs to pull single elements or blocks of them.

        :param n:
        :type n:
        :param m:
        :type m:
        :return:
        :rtype:
        """

        dims = self.dims
        # ndims = int(np.prod(dims))
        idx = (n, m)

        # There are two possible modes for this pulling individual elements or pulling blocks
        # the block pulling can be quite a bit faster, so we try to detect first off if we want to do that
        pull_elements = True
        # first we check if we've got something like from `np.ix_(n, m)`
        if isinstance(n, np.ndarray) and isinstance(m, np.ndarray):
            if len(n.shape) > 1 and len(m.shape) > 1:
                pull_elements = False
        if pull_elements:
            # next we check to see if idx is really just a single element
            pull_elements = all(isinstance(x, (int, np.integer)) for x in idx)
            if not pull_elements:
                # if not a single element, we make sure there are no slices
                pull_elements = all(not isinstance(x, (int, np.integer, slice)) for x in idx)
                # print(">?>", pull_elements, idx)
                if pull_elements:
                    e1 = len(idx[0])
                    pull_elements = all(len(x) == e1 for x in idx)
                    # print(">??", pull_elements)

        # This manages figuring out which indices we pull...

        # We figure out the row spec
        if not isinstance(n, int):
            if isinstance(n, np.ndarray):
                n = n.flatten()
            if not isinstance(n, slice):
                n = np.array(n, dtype=int)
            n = self._get_inds(n)
        else:
            n = [n]
        # Then the column spec
        if not isinstance(m, int):
            if isinstance(m, np.ndarray):
                m = m.flatten()
            if not isinstance(m, slice):
                m = np.array(m)
            m = self._get_inds(m)
        else:
            m = [m]

        if pull_elements:
            # If we're just pulling elements we need only unravel those indices
            n = np.unravel_index(n, dims)
            m = np.unravel_index(m, dims)
        else:
            # If we're pulling blocks we need to compute the product of the row
            #  and column indices to get the total index spec
            blocks = np.array(list(ip.product(n, m)))
            n = np.unravel_index(blocks[:, 0], dims)
            m = np.unravel_index(blocks[:, 1], dims)

        # we define a temporary helper to pad repeat the indices if necessary
        def pad_lens(a, b):
            if isinstance(a, (int, np.integer)) and not isinstance(b, (int, np.integer)):
                a = np.full((len(b),), a)
            if isinstance(b, (int, np.integer)) and not isinstance(a, (int, np.integer)):
                b = np.full((len(a),), b)
            return a, b

        i = tuple(pad_lens(a, b) for a, b in zip(n, m))
        els = self.compute(i)
        if isinstance(els, int) and els == 0:
            # short-circuited :|
            return els

        elif not pull_elements:
            shp = (len(np.unique(blocks[:, 0])), len(np.unique(blocks[:, 1])))
            extra_shp = els.shape[:-1] # some terms will return higher-dimensional results?
            # for sparse arrays this happens in-place :|
            els = els.reshape(extra_shp + shp).squeeze()

        return els

    def __getitem__(self, item):
        if isinstance(item, BraKetSpace): # short circuit here
            return self.get_brakets(item)

        if not isinstance(item, tuple):
            item = (item,)
        if len(item) == 1:
            item = item + (slice(None, None, None),)
        if len(item) > 2:
            raise Exception("index spec '{}' must be of dimension 2".format(item))
        else:
            inds = np.array(item)
            if inds.dtype == np.dtype(int):
                return self.get_brakets(inds)
            else:
                return self.get_element(*item)

    def __rmul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return ExpansionRepresentation([other], [self], self.basis,
                                           name=self.name,
                                           logger=self.logger,
                                           memory_constrained=self.memory_constrained
            )
        else:
            raise TypeError("operator * not defined for objects of type {0} and {1} (only numbers are supported with {0})".format(
                type(self).__name__,
                type(other).__name__
            ))
    def __mul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return ExpansionRepresentation([other], [self],
                                           self.basis,
                                           name=self.name,
                                           logger=self.logger,
                                           memory_constrained=self.memory_constrained
                                           )
        else:
            raise TypeError("operator * not defined for objects of type {0} and {1} (only numbers are supported with {0})".format(
                type(self).__name__,
                type(other).__name__
            ))
    def __add__(self, other):
        if isinstance(other, Representation):
            if other.dims != self.dims:
                if isinstance(other, ExpansionRepresentation):
                    return other + self
                raise ValueError("Can't combine TermComputer objects with dim {} and {}".format(
                    self.dims,
                    other.dims
                ))
            return ExpansionRepresentation([1, 1],
                                           [self, other],
                                           self.basis,
                                           name=self.name + other.name if self.name is not None and other.name is not None else None,
                                           logger=self.logger,
                                           memory_constrained=self.memory_constrained
            )
        else:
            raise TypeError("operator + not defined for objects of type {0} and {1} (only subclasses of TermComputer are supported with {0})".format(
                type(self).__name__,
                type(other).__name__
            ))

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
        if self.name is None:
            return "{}(<{}>, {})".format(
                type(self).__name__,
                self._get_dim_string(self.dims) if hasattr(self, 'dims') else '???',
                self.operator if self.operator is not None else self.compute
            )
        else:
            return "{}<{}>".format(
                type(self).__name__, self.name
            )

    @property
    def selection_rules(self):
        """

        :return:
        :rtype:
        """
        if self._selection_rules is None:
            if self.operator is not None:
                return self.operator.selection_rules

            else:
                return None
        else:
            return self._selection_rules
    @selection_rules.setter
    def selection_rules(self, rules):
        if isinstance(rules[0], (int, np.integer)):
            raise ValueError('selection rules expected to be a list of lists')
        self._selection_rules = rules

    @property
    def selection_rule_steps(self):
        """
        :return:
        :rtype:
        """
        if self._selection_rule_steps is None:
            if self.operator is not None:
                return self.operator.selection_rule_steps

            else:
                return None
        else:
            return self._selection_rule_steps

    def get_transformed_space(self, space, parallelizer=None, logger=None, **opts):
        """
        Returns the state space obtained by using the
        held operator to transform `space`

        :param space:
        :type space:
        :return: connected state spaces
        :rtype: SelectionRuleStateSpace
        """

        if parallelizer is None:
            parallelizer = self.parallelizer

        if self.operator is not None:
            return self.operator.get_transformed_space(space, rules=self.selection_rules, parallelizer=parallelizer, logger=logger,
                                                       **opts
                                                       )
        elif self.selection_rules is not None:
            return space.apply_selection_rules(self.selection_rules, parallelizer=parallelizer, logger=logger,
                                               **opts
                                               )
        else:
            raise ValueError("can't get a transformed space without a held operator or selection rules")


    def apply(self, other):
        raise NotImplementedError("This code path is temporarily abandoned")

        if self.operator is None:
            raise ValueError("")

        if not isinstance(other, (Representation, StateSpaceMatrix)):
            raise TypeError("{} doesn't support application to objects that don't provide a state space".format(
                type(self).__name__
            ))

        if isinstance(other, Representation):
            other = other.array

        other_states = other.brakets.bras
        new_states = self.operator.get_transformed_space(other_states)
        # try:
        brakets = new_states.get_representation_brakets()
        # raise Exception(brakets.bras.excitations[:5], brakets.kets.excitations[:5])
        # except:
        #     raise Exception(other_states, other)

        self.compute_cached(brakets)

        return self.array.dot(other)

    # from memory_profiler import profile
    def _build_rep_matrix(self, row_inds, col_inds, N, sub, method='unique'):

        if method == 'low_mem':
            # logger.log_print("building upper/triangle indices...")
            # upper triangle of indices
            utri_pos = np.where(row_inds < col_inds)[0]
            diag_pos = np.where(row_inds == col_inds)
            if len(diag_pos) > 0:
                diag_pos = diag_pos[0]

            utri_rows = row_inds[utri_pos]
            utri_cols = col_inds[utri_pos]
            if len(utri_pos) + len(diag_pos) < len(row_inds):
                ltri_pos = np.setdiff1d(np.arange(len(row_inds)),
                                        np.concatenate([utri_pos, diag_pos])
                                        )

                # ltri_rows = row_inds[ltri_pos]

                # this was originally just np.unique but the memory
                # for constructing the sort arrays got large?
                # remove dupes by focusing on col inds first
                ltri_cols = col_inds[ltri_pos]
                col_grps = nput.group_by(ltri_pos, ltri_cols)[0]
                unduped_pos = []
                utri_rem_rows = np.full(len(utri_rows), True) # we build a mask to cut down on duplicate tests
                for k,idx in zip(*col_grps):
                    possible_matches = np.where(utri_rows[utri_rem_rows] == k)
                    if len(possible_matches) > 0:
                        possible_matches = possible_matches[0]
                        # now we match corresponding upper-triangle cols against lower-tri rows
                        utri_test_cols = utri_cols[utri_rem_rows][possible_matches]
                        # we mask out the matched positions so we don't test
                        # them on later iterations
                        utri_rem_rows[utri_rem_rows][possible_matches] = False

                        # we assume no duplicates at this point
                        # so we just test the utri_test_cols against
                        # the different column vals to see if we have
                        # any dupes

                        good_spots = np.where(
                            np.all(
                                utri_test_cols[np.newaxis, :] != row_inds[idx][:, np.newaxis],
                                axis=1
                                )
                        )
                        # print(good_spots.shape, len(idx))
                        # print(idx)
                        if len(good_spots) > 0:
                            good_spots = good_spots[0]
                            unduped_pos.append(idx[good_spots])
                    else: # no upper-triangle rows matched this lower-triangle column
                        unduped_pos.append(idx)

                ltri_pos = np.concatenate(unduped_pos)

                tri_pos = np.concatenate([utri_pos, ltri_pos], axis=0)

            else:
                tri_pos = utri_pos

            tri_rows = row_inds[tri_pos]
            tri_cols = col_inds[tri_pos]
            tri_vals = sub[tri_pos]

            up_tri = np.array([tri_rows, tri_cols]).T
            lo_tri = np.array([tri_cols, tri_rows]).T

            if len(diag_pos) > 0: # diags can get duped...somehow
                _, wtf_pos = np.unique(row_inds[diag_pos], return_index=True)
                diag_pos = diag_pos[np.sort(wtf_pos)]

            diag_inds = row_inds[diag_pos]
            diag_vals = sub[diag_pos]
            diag = np.array([diag_inds, diag_inds]).T

            full_inds = np.concatenate([up_tri, lo_tri, diag], axis=0)
            full_dat = np.concatenate([tri_vals, tri_vals, diag_vals], axis=0)

        elif method == 'unique':

            # upper triangle of indices
            up_tri = np.array([row_inds, col_inds]).T
            # lower triangle is made by transposition
            low_tri = np.array([col_inds, row_inds]).T
            # but now we need to remove the duplicates, because many sparse matrix implementations
            # will sum up any repeated elements...
            full_inds = np.concatenate([up_tri, low_tri])
            full_dat = np.concatenate([sub, sub], axis=0)

            _, idx = np.unique(full_inds, axis=0, return_index=True)
            sidx = np.sort(idx)
            full_inds = full_inds[sidx]
            full_dat = full_dat[sidx]

        else:
            raise ValueError('bad method {}'.format(method))

        # N = len(total_space)
        if full_dat.ndim > 1:
            # need to tile the inds enough times to cover the shape of full_dat
            # _, inds = np.unique(full_inds, axis=0, return_index=True)
            # sidx = np.sort(inds)
            # full_inds = full_inds[sidx]
            # full_dat = full_dat[sidx]

            ext_shape = full_dat.shape[1:]
            shape = (N, N) + ext_shape
            full_dat = np.moveaxis(full_dat, 0, -1).flatten()

            rows, cols = full_inds.T
            full_inds = np.empty(
                (2 + len(ext_shape), len(full_dat)),
                dtype=int
            )
            block_size = len(rows)
            # create a tensor of indices that will be used to appropriately tile
            # the system
            ind_tensor = np.moveaxis(np.indices(ext_shape), 0, -1).reshape(-1, len(ext_shape))
            for n, idx in enumerate(ind_tensor):
                s, e = n * block_size, (n + 1) * block_size
                full_inds[0, s:e] = rows
                full_inds[1, s:e] = cols
                for j, x in enumerate(idx):
                    full_inds[2 + j, s:e] = x
        else:
            full_inds = full_inds.T
            shape = (N, N)

        return SparseArray.from_data((full_dat, full_inds), shape=shape, cache_block_data=False)

    # from memory_profiler import profile
    # @profile
    def get_representation_matrix(self,
                                  coupled_space,
                                  total_space,
                                  filter_space=None,
                                  diagonal=False,
                                  logger=None,
                                  zero_element_warning=True,
                                  clear_sparse_caches=True,
                                  clear_operator_caches=True,
                                  assume_symmetric=True,
                                  remove_duplicates=True,
                                  memory_constrained=None
                                  ):
        """
        Actively constructs a perturbation theory Hamiltonian representation

        :param h:
        :type h:
        :param cs:
        :type cs:
        :return:
        :rtype:
        """

        if not assume_symmetric:
            raise NotImplementedError("don't have non-symmetric representations implemented")

        if memory_constrained is None:
            memory_constrained = self.memory_constrained

        if logger is None:
            logger = self.logger

        if not isinstance(coupled_space, BraKetSpace):
            if diagonal:
                m_pairs = BraKetSpace(coupled_space, coupled_space)
            else:
                m_pairs = coupled_space.get_representation_brakets(other=filter_space)
        else:
            m_pairs = coupled_space

        if remove_duplicates:
            m_pairs = m_pairs.remove_duplicates(assume_symmetric=assume_symmetric)

        if len(m_pairs) > 0:
            # logger.log_print(["coupled space dimension {d}"], d=len(m_pairs))
            sub = self.get_brakets(m_pairs, check_orthogonality=not diagonal, memory_constrained=memory_constrained)
            if isinstance(sub, SparseArray):
                sub = sub.asarray()
            if clear_sparse_caches:
                SparseArray.clear_cache()
            if clear_operator_caches:
                self.clear_cache()
        else:
            logger.log_print('no states to couple!', log_level=logger.LogLevel.Debug)
            sub = 0

        logger.log_print("constructing sparse representation...", log_level=logger.LogLevel.Debug)

        N = len(total_space)
        if isinstance(sub, (int, np.integer, np.floating, float)):
            if sub == 0:
                sub = SparseArray.empty((N, N), dtype=float)
            else:
                raise ValueError("Using a constant shift of {} will force representation to be dense...".format(sub))
                sub = np.full((N, N), sub)
        else:
            if zero_element_warning:
                zeros = np.where(sub == 0.)
                if len(zeros) > 0 and len(zeros[0]) > 0:
                    zeros = zeros[0]
                    bad_pairs = (m_pairs.bras.take_subspace(zeros), m_pairs.kets.take_subspace(zeros))
                    raise ValueError(
                        ('got {0} zero elements from {1}; if you expect zeros, set `zero_element_warning=False`. '
                         'First zero element: <|{2}|{1}|{3}>').format(len(zeros), self,
                                                                      bad_pairs[0].excitations[0],
                                                                      bad_pairs[1].excitations[0]
                                                                      )
                    )
            if diagonal and coupled_space is total_space: # fast shortcut
                sub = SparseArray.from_diag(sub)
            else:
                # we assume duplicates were removed earlier

                utri, full_inds = self._pull_rep_mat_inds(logger, total_space, m_pairs)

                logger.log_print("getting vals array...", log_level=logger.LogLevel.Debug)
                full_dat = np.concatenate([sub, sub[utri]], axis=0)
                del sub, utri

                sorting = np.argsort(full_inds[:, 0])
                np.take(full_dat, sorting, out=full_dat, axis=0, mode='clip')
                np.take(full_inds, sorting, out=full_inds, axis=0, mode='clip')
                del sorting

                full_dat, full_inds, shape = self._handle_multidim_rep_mat(logger, full_dat, full_inds, N)

                logger.log_print("building sparse tensor...", log_level=logger.LogLevel.Debug)
                # another attempt
                # if False and memory_constrained:
                #     with tf.NamedTemporaryFile() as full_dat_dump:
                #         with tf.NamedTemporaryFile() as full_inds_dump:
                #             full_dat_map = np.memmap(full_dat_dump.name, dtype=full_dat.dtype, shape=full_dat.shape)
                #             full_dat_map[:] = full_dat[:]
                #             del full_dat
                #
                #             full_inds_map = np.memmap(full_inds_dump.name, dtype=full_inds.dtype, shape=full_inds.shape)
                #             full_inds_map[:] = full_inds[:]
                #             del full_inds
                #
                #             sub = SparseArray.from_data((full_dat_map, full_inds_map),
                #                                         shape=shape,
                #                                         cache_block_data=False,
                #                                         logger=logger
                #                                         # layout='coo'
                #                                         )
                #             # sub2 = SparseArray.from_data((full_dat, full_inds), shape=shape, cache_block_data=False, logger=logger)
                #             #
                #             # raise Exception(np.sum(sub.asarray() - sub2.asarray()))
                #
                # else:
                cs = self.chunk_size
                if cs is not None:
                    blocks = int(np.ceil(len(full_dat) / cs))
                    sub = None
                    for i in range(blocks):
                        s = i*cs; e = min((i+1)*cs, len(full_dat))
                        init = SparseArray.initializer_list(full_dat[s:e], full_inds[:, s:e])
                        subub = SparseArray.from_data(init,
                                                shape=shape,
                                                cache_block_data=False,
                                                logger=logger
                                                , init_kwargs=dict(assume_sorted=True)
                                                )
                        if sub is None:
                            sub = subub
                        else:
                            sub += subub

                    del full_dat
                    del full_inds

                else:
                    init = SparseArray.initializer_list(full_dat, full_inds)
                    del full_dat
                    del full_inds
                    sub = SparseArray.from_data(init,
                                                shape=shape,
                                                cache_block_data=False,
                                                logger=logger
                                                , init_kwargs=dict(assume_sorted=True)
                                                )
                    del init


        return sub

    # @profile
    def _pull_rep_mat_inds(self, logger, total_space, m_pairs):
        logger.log_print("finding row/column indices...", log_level=logger.LogLevel.Debug)
        # figure out the appropriate inds for this data in the sparse representation

        # I should probably add some level of chunking on this?
        # or force `bras` to keep its data on disk?
        row_inds = total_space.find(m_pairs.bras, minimal_dtype=True)
        col_inds = total_space.find(m_pairs.kets, minimal_dtype=True)

        logger.log_print("constructing full row/col index array...", log_level=logger.LogLevel.Debug)
        # ninds = len(row_inds)
        # first_bits = sidx[sidx < ninds]
        # last_bits = sidx[sidx >= ninds] - ninds

        utri = row_inds != col_inds
        full_inds = np.array([
            np.concatenate([row_inds, col_inds[utri]]),
            np.concatenate([col_inds, row_inds[utri]])
        ]).T
        # _, idx = np.unique(full_inds, axis=0, return_index=True)
        # sidx = np.sort(idx)

        # nddup = len(np.unique(diag_inds))
        # logger.log_print("pre/pos: {nfull}/{nun} diags/undiags:{ndig}/{nddup} ??? {ndiff}",
        #                  nfull=len(full_inds), nun=len(sidx),
        #                  ndig=len(diag_inds), nddup=nddup,
        #                  ndiff = len(full_inds) - nddup
        #                  )
        # full_inds = full_inds[sidx]

        return utri, full_inds

    # @profile
    def _handle_multidim_rep_mat(self, logger, full_dat, full_inds, N):

        if full_dat.ndim > 1:
            # need to tile the inds enough times to cover the shape of full_dat

            logger.log_print("making data multi-dim...", log_level=logger.LogLevel.Debug)
            ext_shape = full_dat.shape[1:]
            shape = (N, N) + ext_shape
            full_dat = full_dat.flatten()

            logger.log_print("populating index array...", log_level=logger.LogLevel.Debug)
            rows, cols = full_inds.T
            full_inds = np.empty(
                (len(shape), len(full_dat)),
                dtype=nput.infer_inds_dtype(len(full_dat))
            )

            # create a tensor of indices that will be used to appropriately tile
            # the system
            stride = len(shape)
            # print(len(shape), len(full_dat)/len(rows))
            # row_repeats = np.broadcast_to(rows[:, np.newaxis], (len(rows), len(ext_shape))).flatten()
            # del rows
            # col_repeats = np.broadcast_to(cols[:, np.newaxis], (len(cols), len(ext_shape))).flatten()
            # del cols
            ind_tensor = np.moveaxis(np.indices(ext_shape), 0, -1).reshape(-1, len(ext_shape))
            for n, idx in enumerate(ind_tensor):
                # we want this to remain sorted, so
                # we want to tile this in such a way that
                # all of the things associated with row 0 come first
                # then with row 1
                # then row 2, etc
                # To do so, we fill with an offset of num extra columns
                full_inds[0, n::stride] = rows
                full_inds[1, n::stride] = cols
                for j, x in enumerate(idx):
                    full_inds[2 + j, n::stride] = x
            del rows
            del cols

        else:
            full_inds = full_inds.T
            shape = (N, N)

        return full_dat, full_inds, shape

    def get_diagonal_representation(self,
                                  coupled_space,
                                  total_space,
                                  logger=None,
                                  zero_element_warning=True,
                                  clear_sparse_caches=True
                                  ):
        """
        Actively constructs a perturbation theory Hamiltonian representation

        :param h:
        :type h:
        :param cs:
        :type cs:
        :return:
        :rtype:
        """

        if logger is None:
            logger = self.logger
        m_pairs = coupled_space.get_representation_brakets()

        if len(m_pairs) > 0:
            logger.log_print(["coupled space dimension {d}"], d=len(m_pairs), log_level=logger.LogLevel.Debug)
            sub = self[m_pairs]
            if clear_sparse_caches:
                SparseArray.clear_cache()
        else:
            logger.log_print('no states to couple!', log_level=logger.LogLevel.Debug)
            sub = 0

        logger.log_print("constructing sparse representation...", log_level=logger.LogLevel.Debug)

        N = len(total_space)
        if isinstance(sub, (int, np.integer, np.floating, float)):
            if sub == 0:
                sub = SparseArray.empty((N, N), dtype=float)
            else:
                raise ValueError("Using a constant shift of {} will force representation to be dense...".format(sub))
                sub = np.full((N, N), sub)
        else:
            # figure out the appropriate inds for this data in the sparse representation
            row_inds = total_space.find(m_pairs.bras)
            col_inds = total_space.find(m_pairs.kets)

            if zero_element_warning:
                zeros = np.where(sub == 0.)
                if len(zeros) > 0 and len(zeros[0]) > 0:
                    zeros = zeros[0]
                    bad_pairs = (m_pairs.bras.take_subspace(zeros), m_pairs.kets.take_subspace(zeros))
                    raise ValueError(
                        ('got zero elements from {0}; if you expect zeros, set `zero_element_warning=False`. '
                         'First zero element: <|{1}|{0}|{2}>').format(self,
                                                                      bad_pairs[0].excitations[0],
                                                                      bad_pairs[1].excitations[0]
                                                                      )
                    )

            # upper triangle of indices
            up_tri = np.array([row_inds, col_inds]).T
            # lower triangle is made by transposition
            low_tri = np.array([col_inds, row_inds]).T
            # but now we need to remove the duplicates, because many sparse matrix implementations
            # will sum up any repeated elements
            full_inds = np.concatenate([up_tri, low_tri])
            full_dat = np.concatenate([sub, sub])

            _, idx = np.unique(full_inds, axis=0, return_index=True)
            sidx = np.sort(idx)
            full_inds = full_inds[sidx]
            full_dat = full_dat[sidx]
            sub = SparseArray.from_data((full_dat, full_inds.T), shape=(N, N))

        return sub

class ExpansionRepresentation(Representation):
    """
    Provides support for terms that look like `1/2 pGp + 1/2 dV/dQdQ QQ` by computing each term on its own
    """
    def __init__(self, coeffs, computers, basis, name=None, logger=None, memory_constrained=False):
        """
        :param coeffs: The expansion coefficients
        :type coeffs: Iterable[float]
        :param compute: the functions that turns indices into values
        :type compute: Iterable[callable | Operator]
        :param n_quanta: the total quanta used in the representations (necessary for shape reasons)
        :type n_quanta: tuple[int]
        """
        self.coeffs = np.array(coeffs)
        self.computers = [Representation(c, basis) if not isinstance(c, Representation) else c for c in computers]
        super().__init__(None, basis, name=name, logger=logger, memory_constrained=memory_constrained)

    def clear_cache(self):
        for c in self.computers:
            c.clear_cache()

    def __rmul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return type(self)(self.coeffs * other, self.computers, self.basis,
                              logger=self.logger
                              )
        else:
            raise TypeError(
                "operator * not defined for objects of type {0} and {1} (only numbers are supported with {0})".format(
                    type(self).__name__,
                    type(other).__name__
                ))
    def __mul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return type(self)(self.coeffs*other, self.computers, self.basis,
                              name=self.name,
                              logger=self.logger
                              )
        else:
            raise TypeError(
                "operator * not defined for objects of type {0} and {1} (only numbers are supported with {0})".format(
                    type(self).__name__,
                    type(other).__name__
                ))
    def __add__(self, other):
        if isinstance(other, Representation):
            if other.dims != self.dims:
                raise ValueError("Can't combine TermComputer objects with dim {} and {}".format(
                    self.dims,
                    other.dims
                ))
            if isinstance(other, ExpansionRepresentation):
                return type(self)(
                    np.concatenate([self.coeffs, other.coeffs]),
                    self.computers + other.computers,
                    self.basis,
                    name=self.name + other.name if (self.name is not None and other.name is not None) else None,
                    logger=self.logger
                )
            elif isinstance(other, Representation):
                return type(self)(
                    np.concatenate([self.coeffs, [1]]),
                    self.computers + [other],
                    self.basis,
                    name=self.name + other.name if (self.name is not None and other.name is not None) else None,
                    logger=self.logger
                )
            else:
                raise TypeError(
                    "operator + not defined for objects of type {0} and {1} (only subclasses of TermComputer are supported with {0})".format(
                        type(self).__name__,
                        type(other).__name__
                    ))

        else:
            raise TypeError(
                "operator + not defined for objects of type {0} and {1} (only subclasses of TermComputer are supported with {0})".format(
                    type(self).__name__,
                    type(other).__name__
                ))

    def _dispatch_over_expansion(self, attr, *args, **kwargs):
        els = None
        for c, t in zip(self.coeffs, self.computers):
            if not (isinstance(c, (int, float, np.integer, np.floating)) and c == 0):
                with self.logger.block(tag="in {}".format(t)):
                    start = time.time()
                    bits = getattr(t, attr)(*args, **kwargs)
                    scaled = bits * c
                    if self.ndims == 4:
                        raise ValueError("woop")
                    if isinstance(scaled, (int, float, np.integer, np.floating)) and scaled != 0:
                        # raise Exception(bits, c, scaled, n,m, t)
                        raise ValueError(" ".join([
                            "Adding a constant ({}) to a sparse operator ({}) would cast to dense.",
                            "Error likely occurred in getting elements for {}.",
                            "Explicitly subclass {} if you truly need the constant shift.",
                        ]).format(
                            scaled,
                            els,
                            t,
                            type(self).__name__
                        ))
                    else:
                        if els is None:
                            els = scaled
                        else:
                            if isinstance(scaled, (SparseArray,)):
                                scaled = scaled.asarray()
                            elif isinstance(scaled, (sp.spmatrix,)):
                                scaled = scaled.toarray()
                                # import McUtils.Plots as plt
                                #
                                # plt.ArrayPlot(scaled).show()
                            # print(scaled.shape, els.shape)
                            els += scaled
                    self.logger.log_print("took {e:.3f}s", e=(time.time() - start))
                    # this can be memory intensive so we collect each step
                    gc.collect()

        return els

    def get_brakets(self, states, check_orthogonality=True, memory_constrained=False):
        if not isinstance(states, BraKetSpace):
            states = BraKetSpace.from_indices(states, basis=self.basis)
        return self._dispatch_over_expansion('get_brakets', states, check_orthogonality=check_orthogonality, memory_constrained=memory_constrained)
    def get_element(self, n, m):
        return self._dispatch_over_expansion('get_element', n, m)
    @property
    def chunk_size(self):
        for c in self.computers:
            if c.chunk_size is not None:
                return c.chunk_size
        else:
            return None

    @property
    def selection_rules(self):
        """

        :return:
        :rtype:
        """
        if self._selection_rules is None:
            if all(hasattr(x, 'selection_rules') for x in self.computers):
                _ = []
                for x in self.computers:
                    for r in x.selection_rules:
                        if r not in _:
                            _.append(r)
                self._selection_rules = _
        return self._selection_rules

    @property
    def selection_rule_steps(self):
        """

        :return:
        :rtype:
        """
        if self._selection_rule_steps is None:
            if all(hasattr(x, 'selection_rule_steps') for x in self.computers):
                _ = []
                for x in self.computers:
                    r = x.selection_rule_steps
                    if r not in _:
                        _.append(r)
                self._selection_rule_steps = _
        return self._selection_rule_steps

    def get_transformed_space(self, space, rules=None, parallelizer=None, logger=None, **opts):
        """
        Returns the state space obtained by using the
        held operators to transform `space`

        :param space:
        :type space: BasisStateSpace
        :return:
        :rtype: SelectionRuleStateSpace
        """
        import functools

        if parallelizer is None:
            parallelizer = self.parallelizer

        # we take a union of all transformation rules and just apply that
        # if possible
        if rules is not None:
            ooooh_shiz = space.apply_selection_rules(rules, parallelizer=parallelizer, logger=logger, **opts)
        elif self.selection_rules is not None:
            ooooh_shiz = space.apply_selection_rules(self.selection_rules, parallelizer=parallelizer, logger=logger, **opts)
        else:
            spaces = [r.get_transformed_space(space, parallelizer=parallelizer, logger=logger, **opts) for r in self.computers]
            ooooh_shiz = functools.reduce(lambda s1,s2: s1.union(s2), spaces[1:], spaces[0])

        return ooooh_shiz

    def __repr__(self):
        if self.name is None:
            return "{}(<{}>, ({}), ({}))".format(
                type(self).__name__,
                self._get_dim_string(self.dims) if hasattr(self, 'dims') else '???',
                self.coeffs,
                self.computers
            )
        else:
            return "{}<{}>".format(
                type(self).__name__, self.name
            )