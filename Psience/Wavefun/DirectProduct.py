
import numpy as np
from ..BasisReps import PermutationStateIndexer
from .Wavefunctions import Wavefunctions, Wavefunction

__all__ = [
    "DirectProductWavefunction",
    "DirectProductWavefunctions",
]
class DirectProductWavefunction(Wavefunction):
    def __init__(self, wfns:'Iterable[Wavefunction]', parent=None, index=None, **opts):
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
        if isinstance(n, int):
            ...
        else:
            np.arange()
        self.indexer.from_indices(len)

    def _get_slice_inds(self, n):
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
        return np.prod([len(w) for w in self.base_wfns]).astype(int)

    def get_energies(self, idx):
        base_idx = self._get_slice_inds(idx).T
        return sum(w.energies[i] for w,i in zip(self.base_wfns, base_idx))
    @property
    def energies(self):
        return self.get_energies(None)

    def get_wavefunctions(self, which):
        base_idx = self._get_slice_inds(which).T
        wfn_bits = [w.wavefunctions[i] for w, i in zip(self.base_wfns, base_idx)]
        if isinstance(wfn_bits[0], Wavefunction):
            return self.wavefunction_class(wfn_bits, parent=self, index=which)
        else:
            return type(self)(wfn_bits, indices=base_idx.T)