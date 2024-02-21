import numpy as np, scipy.linalg as slag, collections

from McUtils.Scaffolding import Logger
from McUtils.Data import UnitsData
from .NormalModes import NormalModes

__all__ = [
    'ObliqueModeGenerator'
]

class ObliqueModeGenerator:

    def __init__(self, f, g, dimensionless=False, sel=None, frequency_scaled=True):
        f = np.asanyarray(f)
        g = np.asanyarray(g)
        if g.ndim == 1: # vector of masses
            g = np.broadcast_to(1/g[:, np.newaxis], (len(g), len(f) // len(g))).flatten()
            g = np.diag(g)
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

        freqs, modes, inv = NormalModes.get_normal_modes(
            f, g,
            remove_transrot=True,
            dimensionless=frequency_scaled
        )

        rU, a, rM = np.linalg.svd(modes)
        good_a = np.where(a > 1e-8)[0]
        R = (rU[:, good_a] @ rM[good_a, :])

        self.modes = modes
        self.rotation = R
        self.inverse_scaling = modes @ R.T
        self.scaling = R @ inv

    @classmethod
    def from_molecule(cls, mol, dimensionless=True, sel=None, use_internals=None, frequency_scaled=True):
        from ..Molecools import Molecule
        mol = mol # type:Molecule
        if use_internals is None:
            use_internals = mol.internal_coordinates is not None
        if use_internals:
            return cls(
                mol.get_internal_potential_derivatives(2)[1],
                mol.g_matrix,
                dimensionless=dimensionless,
                sel=sel,
                frequency_scaled=frequency_scaled
            )
        else:
            return cls(
                mol.potential_derivatives[1],
                mol.atomic_masses,
                dimensionless=dimensionless,
                sel=sel,
                frequency_scaled=frequency_scaled
            )

    def run(self, scaling_type='normal'):
        u = self.scaling
        ui = self.inverse_scaling
        if scaling_type == 'inverse':
            u, ui = ui, u
        g = u @ self.g @ u.T
        f = ui.T @ self.f @ ui
        return f, g, u, ui

