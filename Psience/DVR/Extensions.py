
import numpy as np
from ..VSCF import GridSCF
from McUtils.Scaffolding import ParameterManager

from .BaseDVR import BaseDVR
from .Wavefunctions import DVRWavefunctions
from .DirectProduct import DirectProductDVR
from .FiniteBasisDVR import WavefunctionBasisDVR

__all__ = [
    "SelfConsistentDVR",
    "PotentialOptimizedDVR"
]

class SCFWavefunctionGenerator:
    def __init__(self, dvr_1D:BaseDVR):
        self.dvr = dvr_1D
        self.prev = None
    def __call__(self, pot, **kwargs):
        if self.prev is None:
            res = self.prev = self.dvr.run(potential_values=pot)
        else:
            res = self.dvr.run(
                potential_values=pot,
                grid=self.prev.grid,
                kinetic_energy=self.prev.kinetic_energy
            )
        return res.wavefunctions

class SelfConsistentDVR(GridSCF):
    def __init__(self, base_dvr:"DirectProductDVR", **opts):
        props = ParameterManager(**opts)
        self.base_dvr = base_dvr
        generators = [SCFWavefunctionGenerator(d) for d in self.base_dvr.dvrs]
        pot_data = self.base_dvr.run(result='potential_energy')
        grid = pot_data.grid
        pe = pot_data.potential_energy
        if hasattr(pe, 'diagonal'):
            pe = pe.diagonal()
        else:
            pe = np.diag(pe)
        pe = pe.reshape(grid.shape[:-1])
        super().__init__(grid, pe, generators, **props.filter(GridSCF))
    # def initialize(self):
    #     d = super().initialize()
    #     for w in d.wavefunctions:
    #         w[:1].plot()
    #     import McUtils.Plots as plt
    #     plt.DensityPlot(*self.grid.transpose(2, 0, 1), self.vals, plot_style=dict(vmin=0, vmax=1)).show()
    #     return d
    def __repr__(self):
        return "{}({})".format(
            type(self).__name__,
            self.base_dvr
        )

class PotentialOptimizedDVR(DirectProductDVR):
    def __init__(self,
                 wfns_1D:'Iterable[DVRWavefunctions]',
                 **base_opts
                 ):
        base_opts = {k:base_opts[k] for k in base_opts.keys() - {'mass', 'g', 'g_deriv', 'include_kinetic_coupling'}}
        super().__init__(
            [WavefunctionBasisDVR(w) for w in wfns_1D],
            mass=1,
            g=None, #
            g_deriv=None,
            include_kinetic_coupling=False,
            **base_opts
        )

    @classmethod
    def from_scf(cls, scf_dvr:SelfConsistentDVR, wfns=None, **opts):
        if wfns is None:
            wfns = scf_dvr.run().wavefunctions
        return cls(
            wfns,
            **dict(
                scf_dvr.base_dvr.opts,
                **opts
            )
        )