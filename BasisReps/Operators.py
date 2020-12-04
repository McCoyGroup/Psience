"""
Provides the operator representations needed when building a Hamiltonian representation.
I chose to only implement direct product operators. Not sure if we'll need a 1D base class...
"""

import numpy as np, scipy.sparse as sp, os, tempfile as tf
from collections import OrderedDict
from McUtils.Numputils import SparseArray

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

    def _get_eye_tensor(self, states, quants):
        orthog = self._get_orthogonality([], states, quants)
        non_orthog = np.where(orthog != 0.)[0]
        return sp.csr_matrix(
            (
                orthog[non_orthog],
                (
                    np.zeros(len(non_orthog)),
                    non_orthog
                )
            ),
            shape=(1, len(orthog))
        )

    def _get_orthogonality(self, inds, states, quants):
        #TODO: this can mostly be delegated to BraKetSpace

        nstates = len(states[0][0])
        ndim = len(quants)
        uinds = np.unique(inds)
        if len(uinds) == ndim:  # means we've got no states that are unused (can only rarely happen, but important edge case)
            orthog = np.ones(nstates, dtype=int)
        else:
            # determine which of the unused states (`missing`) have different numbers of
            # quanta in the bra (`states[i][0]`) than the ket (`states[i][1]`)
            missing = [i for i in range(ndim) if i not in uinds]
            equivs = [states[i][0] == states[i][1] for i in missing]
            orthog = np.prod(equivs, axis=0).astype(int)  # taking advantage of the fact that bools act like ints
        return orthog

    def _apply_sel_rules(self, states, rules):
        #TODO: this can be delegated to BraKetSpace
        import functools as fp

        # we have a set of rules for every dimension in states...
        # so we define a function that will tell us where any of the possible selection rules are satisfied
        # and then we apply that in tandem to the quanta changes between bra and ket and the rules
        apply_rules = lambda diff, rule: fp.reduce(
            lambda d, i: np.logical_or(d, diff==i),
            rule[1:],
            diff==rule[0]
        )
        sels = [
            apply_rules(x[1] - x[0], r) for x, r in zip(states, rules)
        ]
        # then we figure out which states by satisfied all the rules
        all_sels = fp.reduce(np.logical_and, sels[1:], sels[0])
        # and then filter the OG lists by those
        states = [(s[0][all_sels], s[1][all_sels]) for s in states]

        return states, all_sels

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

        # We figure out which indices are actually unique; this gives us a way to map
        # indices onto operators
        # We assume we always have as many indices as dimensions in our operator
        # and so the more unique indices the _lower_ the dimension of the operator
        # since it acts on the same index more times
        # The indices will be mapped onto ints for storage in `pieces`
        uinds = np.unique(inds)
        mm = {k: i for i, k in enumerate(uinds)}

        non_orthog_states, non_orthog_spec = states.apply_non_orthogonality(uinds, assume_unique=True)

        subdim = len(uinds)
        # and then apply selection rules if we have them
        if sel_rules is not None:
            # we need selection rules for every dimension in the product
            sel_bits = [None] * subdim
            for s, i in zip(sel_rules, inds):
                n = mm[i]  # makes sure that we fill in in the same order as uinds
                if sel_bits[n] is None:
                    sel_bits[n] = s
                else:
                    sel_bits[n] = np.unique(np.add.outer(s, sel_bits[n])).flatten() # make the next set of rules
            # print(sel_bits)
            # apply full selection rules
            non_orthog_states, all_sels = non_orthog_states.apply_sel_rules(sel_bits)
            if len(non_orthog_states) == 0:
                return None # short-circuit because there's nothing to calculate

        else:
            all_sels = None

        # determine what size we need for the 1D reps
        # max_dim = np.max([
        #     np.max(non_orthog_states.bras.excitations),
        #     np.max(non_orthog_states.kets.excitations)
        #     ])
        max_dim = np.max(np.array(non_orthog_states))
        padding = 3  # for approximation reasons we need to pad our representations...

        # then set up a place to hold onto the pieces we'll calculate
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
        for s, o in zip(non_orthog_states, pieces):
            op = o  # type: sp.spmatrix
            blob = np.asarray(op[s]).squeeze()
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

    def _direct_prod_operator_terms(self, inds, func, states, non_orthog):
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

        # We figure out which indices are actually unique; this gives us a way to map
        # indices onto operators for our generator
        # non_orthog_states = states.apply_non_orthog(inds)
        _, idx = np.unique(inds, return_index=True)
        uinds = inds[np.sort(idx)]
        # print(inds, uinds)
        # mm = {k: i for i, k in enumerate(uinds)}

        # then we determine which states are potentially non-orthogonal based on precomputed info
        non_orthog_states = [
            (states[i][0][non_orthog], states[i][1][non_orthog])
            for i in uinds
        ]

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

        if sel_rules is not None:
            # apply full selection rules
            non_orthog_states, all_sels = self._apply_sel_rules(non_orthog_states, sel_rules)
            if len(non_orthog_states[0][0]) == 0:
                return None  # short-circuit because there's nothing to calculate
        else:
            all_sels = None

        chunk = gen(non_orthog_states)

        return chunk, all_sels

    def _calculate_single_pop_elements(self, inds, funcs, states, quants, sel_rules):
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
        :type states: Iterable[Iterable[Iterable[int]]]
        :param quants: the total quanta to use when building representations
        :type quants:
        :param sel_rules: the selection rules to use when building representations
        :type sel_rules:
        :return:
        :rtype:
        """

        # determine how many states aren't potentially coupled by the operator
        # & then determine which of those are non-orthogonal

        # nstates = len(states)
        #
        # ndim = len(quants)
        # uinds = np.unique(inds)

        # raise NotImplementedError("Currently adding support for much better orthogonality handling")
        nstates = len(states[0][0])
        orthog = self._get_orthogonality(inds, states, quants)
        single_state = isinstance(orthog, (int, np.integer)) # means we got a single state to calculate over
        if single_state:
            orthog = np.array([orthog])

        non_orthog = np.where(orthog != 0)[0]

        # if none of the states are non-orthogonal...just don't calculate anything
        if len(non_orthog) == 0:
            return sp.csr_matrix((1, nstates), dtype='float')
        else:

            if isinstance(funcs, (list, tuple)):
                chunk = self._mat_prod_operator_terms(inds, funcs, states, non_orthog, sel_rules)
            else:
                chunk = self._direct_prod_operator_terms(inds, funcs, states, non_orthog)

            if chunk is None: # after application of sel rules we get nothing
                return sp.csr_matrix((1, nstates), dtype='float')
            else:
                chunk, all_sels = chunk # secondary application of selection rules
                if all_sels is not None:
                    non_orthog = non_orthog[all_sels,]

            non_zero = np.where(np.abs(chunk) >= self.zero_threshold)[0] # by filtering again we save on computation in the dot products
            # print(len(non_zero), len(non_orthog))

            if len(non_zero) == 0:
                return sp.csr_matrix((1, nstates), dtype='float')

            non_orthog = non_orthog[non_zero,]
            # print(chunk)
            chunk = chunk[non_zero,]

            wat = sp.csr_matrix(
                (
                    chunk,
                    (
                        np.zeros(len(non_orthog)),
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
                                  -1, inds, self.funcs, idx, self.quanta, self.sel_rules
                                  )
        if save_to_disk:
            new = sp.vstack([x for y in res.flatten() for x in y])
            with tf.NamedTemporaryFile() as tmp:
                res = tmp.name + ".npz"
            # print(res)
            # os.remove(tmp.name)
            sp.save_npz(res, new, compressed=False)

        return res

    def _get_pop_parallel(self, inds, idx, parallelizer="multiprocessing"):
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

    def get_elements(self, idx, parallelizer=None):#parallelizer='multiprocessing'):
        """
        Calculates a subset of elements

        :param idx: bra and ket states as tuples of elements
        :type idx: Iterable[(Iterable[int], Iterable[int])]
        :return:
        :rtype:
        """

        # we expect this to be an iterable object that looks like
        # num_modes X [bra, ket] X quanta
        idx = tuple(
            tuple(np.array([i]) if isinstance(i, (int, np.integer)) else np.array(i) for i in j)
            for j in idx
        )

        if len(idx) != self.mode_n:
            raise ValueError("BraKet spec {} isn't valid for Operator with dimension {}".format(
                idx,
                (self.ndim, self.mode_n)
            ))


        # if len(idx) != self.mode_n:
        #     raise ValueError("BraKet spec {} isn't valid for Operator with dimension {}".format(
        #         idx,
        #         (self.ndim, self.mode_n)
        #     ))

        inds = self.get_inner_indices()

        if inds is None: # just a number
            new = self._get_eye_tensor(idx, self.quanta)
            # raise Exception(new[0], type(new), new.shape)
        else:

            # self.symmetry_inds = None
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

            # if inverse is not None:
            #     # raise Exception(res)
            #
            #     new = new.reshape((1, int(np.prod(new.shape))))
            #     print(new, inverse)
            #     new = new[inverse]




        shp = inds.shape if inds is not None else ()
        # print(type(new), new.shape)
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

