
# from .DVR import DVRConstructor
from .Wavefunctions import DVRWavefunctions
from .DirectProduct import DirectProductDVR
# from .SCF import SelfConsistentDVR

import numpy as np

class PotentialOptimizedDVR(DirectProductDVR):
    def __init__(self,
                 dvrs_1D,
                 wfns_1D:'Iterable[DVRWavefunctions]',
                 **base_opts
                 ):
        super().__init__(dvrs_1D, **base_opts)
        self.presolutions = wfns_1D
        # self.opt_divs = [opt_divs] * len(dvrs_1D) if isinstance(opt_divs, (int, np.integer)) else opt_divs

    def get_grid(self, domain=None, divs=None, **kwargs):
        subgrids = []
        orderings = []
        for wfs in self.presolutions:
            x = wfs.expectation(lambda w:wfs.grid, wfs)
            print(x)
            x = np.diag(x)
            ordering = np.argsort(x)
            orderings.append(ordering)
            subgrids.append(x[ordering])

        mesh = np.array(np.meshgrid(*subgrids, indexing='ij'))
        MEHSH = np.moveaxis(mesh, 0, len(subgrids))
        raise Exception("...")

        return MEHSH, orderings
    def get_kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, include_kinetic_coupling=True, **kwargs):
        raise NotImplementedError("...")

    # @classmethod
    # def from_scf(cls, scf_dvr):
    #     base_nd = DVRConstructor.construct(**opts)
    #     res = base_nd.run(result='potential_energy')


