
import numpy as np
from ..BasisReps import PermutationStateIndexer
from .Wavefunctions import Wavefunctions, Wavefunction

__all__ = [
    "DirectProductWavefunction",
    "DirectProductWavefunctions",
]
class DirectProductWavefunction(Wavefunction):
    def __init__(self, wfns:'Iterable[Wavefunction]', parent=None, index=None, **opts):
        """
        **LLM Docstring**

        Build a single direct-product wavefunction (a product state across multiple independent degrees of freedom) from a collection of component wavefunctions, with its energy given by the sum of the components' energies.

        :param wfns: the component wavefunctions making up this product state
        :type wfns: Iterable[Wavefunction]
        :param parent: the parent `Wavefunctions` collection this wavefunction belongs to
        :type parent: Wavefunctions | None
        :param index: this wavefunction's index within its parent collection
        :type index: Any | None
        :param opts: extra options forwarded to the base `Wavefunction.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        energies = [w.energy for w in wfns]
        self.energies = energies
        super().__init__(
            sum(energies),
            wfns,
            parent=parent,
            index=index,
            **opts
        )

class DirectProductWavefunctions(Wavefunctions):
    wavefunction_class = DirectProductWavefunction
    def __init__(self,
                 wfns,
                 indices=None, wavefunction_class=None, **opts):
        """
        **LLM Docstring**

        Build a collection of direct-product wavefunctions spanning the full tensor-product state space of a set of per-degree-of-freedom wavefunction collections, indexed via a `PermutationStateIndexer`.

        :param wfns: the per-degree-of-freedom wavefunction collections making up the product space
        :type wfns: Iterable[Wavefunctions]
        :param indices: explicit state indices this collection represents, if a subset
        :type indices: Any | None
        :param wavefunction_class: the class used to build individual product wavefunctions; defaults to `cls.wavefunction_class`
        :type wavefunction_class: type | None
        :param opts: extra options forwarded to the base `Wavefunctions.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        self.base_wfns = wfns
        self.indexer = PermutationStateIndexer(len(wfns))
        super().__init__(
            None,
            wfns,
            indices=indices,
            wavefunction_class=wavefunction_class,
            **opts
        )

    def get_slice(self, n):
        """
        **LLM Docstring**

        Intended to retrieve a slice of product wavefunctions by an integer count or other index specification

        :param n: the slice specification
        :type n: int | object
        """
        raise NotImplementedError("simple slice ill-defined fro a direct product wavefunction")

    def _get_slice_inds(self, n):
        """
        **LLM Docstring**

        Resolve a flat index/slice specification into the corresponding per-degree-of-freedom quantum number indices, via the `PermutationStateIndexer`.

        :param n: an integer (treated as `arange(n)`), `None` (treated as all states, `arange(len(self))`), a `slice`, or an explicit array of flat state indices
        :type n: int | slice | np.ndarray | None
        :return: the per-degree-of-freedom index array, one row per selected state
        :rtype: np.ndarray
        """
        if isinstance(n, int):
            idx = np.arange(n)
        elif n is None:
            idx = np.arange(len(self))
        elif isinstance(n, slice):
            stop = n.stop
            if stop is None:
                stop = len(self)
            idx = np.arange(stop)[n]
        else:
            idx = n
        return self.indexer.from_indices(idx)

    def __len__(self):
        """
        **LLM Docstring**

        The total number of product states, the product of the sizes of each component wavefunction collection.

        :return: the number of product states
        :rtype: int
        """
        return np.prod([len(w) for w in self.base_wfns]).astype(int)

    def get_energies(self, idx):
        """
        **LLM Docstring**

        Compute the total energies of a set of product states by summing the corresponding component energies.

        :param idx: the state index/indices to compute energies for, in any form accepted by `_get_slice_inds`
        :type idx: Any
        :return: the summed energies
        :rtype: np.ndarray
        """
        base_idx = self._get_slice_inds(idx).T
        return sum(w.energies[i] for w,i in zip(self.base_wfns, base_idx))
    @property
    def energies(self):
        """
        **LLM Docstring**

        The energies of every product state in this collection.

        :return: the full set of product-state energies
        :rtype: np.ndarray
        """
        return self.get_energies(None)

    def get_wavefunctions(self, which):
        """
        **LLM Docstring**

        Build the requested product wavefunction(s) by selecting the corresponding component wavefunctions from each base collection and either combining them into a single `DirectProductWavefunction` (if a single state was requested) or wrapping the resulting sub-collections in a new `DirectProductWavefunctions`.

        :param which: the state index/indices to retrieve, in any form accepted by `_get_slice_inds`
        :type which: Any
        :return: the requested product wavefunction, or a `DirectProductWavefunctions` collection if multiple states were requested
        :rtype: DirectProductWavefunction | DirectProductWavefunctions
        """
        base_idx = self._get_slice_inds(which).T
        wfn_bits = [w.wavefunctions[i] for w, i in zip(self.base_wfns, base_idx)]
        if isinstance(wfn_bits[0], Wavefunction):
            return self.wavefunction_class(wfn_bits, parent=self, index=which)
        else:
            return type(self)(wfn_bits, indices=base_idx.T)