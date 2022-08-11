import numpy as np

from .BaseDVR import BaseDVR
import typing

class InitialBasis(typing.Protocol):
    @property
    def dimensions(self)->'Iterable[int]':
        raise NotImplementedError("abstract")
    def x(self, n:int)->np.ndarray:
        raise NotImplementedError("abstract")
    def p2(self, n:int)->np.ndarray:
        raise NotImplementedError("abstract")
    def p(self, n:int)->np.ndarray:
        raise NotImplementedError("abstract")

class FiniteBasisDVR(BaseDVR):
    def __init__(self, basis:InitialBasis, domain=None, divs=None, **opts):
        self.basis = basis
        super().__init__(domain=None, divs=self.basis.dimensions[0], **opts)

    def get_grid(self, domain=None, divs=None, **kwargs):
        x = self.basis.x(self.divs)
        if hasattr(x, 'toarray'):
            x = x.toarray()
        return np.linalg.eigh(x)

    def real_momentum(self, grid=None, mass=None, hb=1, **kwargs):
        gps, tf = grid
        p = self.basis.p(self.divs)
        if hasattr(p, 'toarray'):
            p = p.toarray()
        return np.dot(np.dot(tf.T, p), tf)

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        gps, tf = grid
        # try:
        #     base_ke = self.basis.kinetic_energy(self.divs)
        # except (AttributeError, NotImplementedError):
        try:
            p2 = self.basis.p2(self.divs)
        except (AttributeError, NotImplementedError):
            p = self.basis.p(self.divs)
            if hasattr(p, 'toarray'):
                p = p.toarray()
            p2 = np.dot(p, p)
        base_ke = -hb / (2 * mass) * p2
        return np.dot(np.dot(tf.T, base_ke), tf)

    def potential_energy(self,
                         grid=None,
                         potential_function=None,
                         potential_values=None,
                         potential_grid=None,
                         logger=None,
                         **pars
                         ):
        return super().potential_energy(
            grid[0],
            potential_function=potential_function,
            potential_values=potential_values,
            potential_grid=potential_grid,
            logger=logger,
            **pars
        )

class HarmonicDVR(FiniteBasisDVR):
    def __init__(self, divs=None, **opts):
        from ..BasisReps import HarmonicOscillatorBasis
        super().__init__(
            HarmonicOscillatorBasis(divs),
            **opts
        )

class WavefunctionBasisDVR(FiniteBasisDVR):
    def __init__(self, wfns=None, **opts):
        from ..BasisReps import WavefunctionBasis
        super().__init__(
            WavefunctionBasis(wfns),
            **opts
        )

