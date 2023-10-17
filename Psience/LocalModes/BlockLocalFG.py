import numpy as np, scipy.linalg as slag, collections

from McUtils.Scaffolding import Logger
from McUtils.Data import UnitsData
from ..Molecools import Molecule

__all__ = [
    'BlockLocalFGOrthogonalizer',
    'BlockLocalFGOrthogonalizerIterative'
]

class BlockLocalFGOrthogonalizer:

    def __init__(self, f, g, dimensionless=False, sel=None):
        f = np.asanyarray(f)
        g = np.asanyarray(g)
        if sel is not None:
            f = f[sel, :][:, sel]
            g = g[sel, :][:, sel]
        if dimensionless:
            a = np.diag(np.power(np.diag(f) / np.diag(g), 1 / 4))
            ai = np.diag(np.power(np.diag(g) / np.diag(f), 1 / 4))
            f = ai @ f @ ai
            g = a @ g @ a
        self.f = f
        self.g = g
        self.ndim = f.shape[0]

        ## Alternate formulation that makes U by default, not U^-1
        # freq2, modes = slag.eigh(f, g, type=2)
        # freqs = np.sqrt(freq2)
        # modes = modes * np.sqrt(freqs)[np.newaxis, :]
        # rU, a, rM = np.linalg.svd(modes)
        # R = (rU @ rM)
        # print(modes @ R.T)

        freq2, modes = slag.eigh(f, g, type=3)
        freqs = np.sqrt(freq2)
        modes = modes / np.sqrt(freqs)[np.newaxis, :]
        rU, a, rM = np.linalg.svd(modes)
        R = (rU @ rM)

        # with np.printoptions(linewidth=1e8, suppress=True):
        #     print(np.round(r, 8))

        self.modes = modes
        self.rotation = R
        self.scaling = modes @ R.T # wow...

    @classmethod
    def from_molecule(cls, mol: Molecule, dimensionless=True, sel=None):
        return cls(
            mol.get_internal_potential_derivatives(2)[1],
            mol.g_matrix,
            dimensionless=dimensionless,
            sel=sel
        )

    def run(self):
        ui = self.scaling
        u = np.linalg.inv(ui)
        g = u @ self.g @ u
        f = ui @ self.f @ ui
        return f, g, u, ui


class BlockLocalFGOrthogonalizerIterative:

    default_tolerance=1e-10
    default_max_iterations=250
    default_scaling=1
    default_damping=1
    def __init__(self, f, g, tolerance=None, max_iterations=None, scaling=None, damping=1, logger=True):
        f = np.asanyarray(f)
        g = np.asanyarray(g)
        self.f = f
        self.g = g
        self.ndim = f.shape[0]
        a = np.diag(np.power(np.diag(f) / np.diag(g), 1/4))
        ai = np.diag(np.power(np.diag(g) / np.diag(f), 1/4))
        self._f = ai @ f @ ai
        self._g = a @ g @ a
        self.u = a
        self.damping = damping
        self.scaling = scaling
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.logger = Logger.lookup(logger)

    @staticmethod
    def _step(f, g, scaling=1, rescale=True):
        wf = np.diag(f)
        wg = np.diag(g)
        u = (f - g) / (
                wf[:, np.newaxis] + wf[np.newaxis, :]
                    + wg[np.newaxis, :] + wg[:, np.newaxis]
        ) / scaling
        np.fill_diagonal(u, 1)
        ui = np.linalg.inv(u)
        g = u @ g @ u
        f = ui @ f @ ui
        if rescale:
            a = np.diag(np.power(np.diag(f) / np.diag(g), 1/4))
            ai = np.diag(np.power(np.diag(g) / np.diag(f), 1/4))
            f = ai @ f @ ai
            g = a @ g @ a
            u = u @ a

        return f, g, u

    @staticmethod
    def off_diag_norm(f, g):
        d = f - g
        inds = np.triu_indices_from(d, 1)
        return np.linalg.norm(d[inds])
    @classmethod
    def iterate(cls, f, g, rescale=True, scaling=None, damping=None, tolerance=None, max_iterations=None):
        if tolerance is None:
            tolerance = cls.default_tolerance
        if scaling is None:
            scaling = cls.default_scaling
        if damping is None:
            damping = cls.default_damping
        if max_iterations is None:
            max_iterations = cls.default_max_iterations
        cur_tol = cls.off_diag_norm(f, g)
        u = np.eye(f.shape[0])
        conv = collections.deque()
        conv.append(cur_tol)
        for _ in range(max_iterations):
            if cur_tol < tolerance:
                break
            f, g, u_ = cls._step(f, g, scaling=scaling, rescale=rescale)
            scaling = scaling * damping
            u = u @ u_
            cur_tol = cls.off_diag_norm(f, g)
            conv.append(cur_tol)
        num_its = len(conv) - 1
        return f, g, u, cur_tol, num_its, conv

    def run(self):
        evs = slag.eigvalsh(self._f, self._g, type=2)
        if np.min(evs) <= 0:
            raise ValueError("need all real frequencies")
        self.logger.log_print("Frequencies: {f}".format(f=np.sqrt(evs)*UnitsData.convert("Hartrees", "Wavenumbers")))
        self.logger.log_print("Condition Number: {f}".format(f=np.max(evs)/np.min(evs)))
        f, g, u, cur_tol, num_its, conv = self.iterate(
            self._f,
            self._g,
            damping=self.damping,
            scaling=self.scaling,
            tolerance=self.tolerance,
            max_iterations=self.max_iterations
        )
        tol = self.default_tolerance if self.tolerance is None else self.tolerance
        if cur_tol > tol:
            self.logger.log_print("Local block FG failed to converge after {i} iterations (off-diag norm {n} greater than tolerance {t})",
                                  i=num_its,
                                  n=cur_tol,
                                  t=tol
                                  )
        u = self.u @ u
        return f, g, u, num_its, conv


    @classmethod
    def from_molecule(cls, mol:Molecule, **opts):
        return cls(
            mol.get_internal_potential_derivatives(2)[1],
            mol.g_matrix,
            **opts
        )

