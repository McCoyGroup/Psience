"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""

import numpy as np

from McUtils.Zachary import Mesh

from Psience.Wavefun import Wavefunction, Wavefunctions

__all__ = ["DGBWavefunctions", "DGBWavefunction"]

class DGBWavefunction(Wavefunction):
    def __init__(self, energy, data, centers=None, alphas=None, inds=None, **opts):
        super().__init__(energy, data, **opts)
        if centers is None:
            centers = self.parent.centers
        self.centers = centers
        if alphas is None:
            alphas = self.parent.alphas
        self.alphas = alphas
        if inds is None and self.parent is not None:
            inds = self.parent.hamiltonian.inds
        self.inds = inds

    def plot(self,
             figure=None,
             domain=None,
             **opts
             ):

        if domain is None:
            domain = Mesh(self.centers).bounding_box

        return super().plot(figure=figure, domain=domain, **opts)

    def evaluate(self, points):
        points = np.asanyarray(points)
        reshape = None
        if points.ndim == 1:
            points = points[np.newaxis]
        elif points.ndim > 2:
            reshape = points.shape[:-1]
            points.reshape(-1, points.shape[-1])

        # actual basis evaluation
        c_disps = points[np.newaxis, :, :] - self.centers[:, np.newaxis, :]
        alphas = self.alphas[:, np.newaxis]
        if self.inds is not None:
            c_disps = c_disps[..., self.parent.inds]
            alphas = alphas[..., self.parent.inds]
        normas = (2 * alphas / np.pi) ** (1/4)
        exp_evals = np.exp(-alphas * c_disps**2)

        vals = np.dot(
            self.data,
            np.prod(normas * exp_evals, axis=-1)
        )
        if reshape is not None:
            vals = vals.reshape(reshape)
        return vals
    def expectation(self, op, other=None):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """
        #
        raise NotImplementedError("operator evaluation for DGBs not implemented yet")
    def project(self, dofs):
        """
        Computes the projection of the current wavefunction onto a set of degrees
        of freedom, returning a projected wave function object

        :return:
        :rtype: Wavefunction
        """
        if isinstance(dofs, (int, np.integer)):
            dofs = [dofs]

        # assert self.alphas.ndim == 1

        proj_dofs = dofs if self.inds is None else np.intersect1d(dofs, self.inds)
        scaling = np.power(2 * np.pi, len(dofs)/4) / np.power(np.prod(self.alphas[:, proj_dofs], axis=-1), 1/4)
        remaining = np.setdiff1d(np.arange(self.centers.shape[-1]), dofs)
        centers = self.centers[:, remaining]
        alphas = self.alphas[:, remaining]

        return type(self)(
            self.energy,
            self.data * scaling,
            centers=centers,
            alphas=alphas
        )

class DGBWavefunctions(Wavefunctions):
    wavefunction_class = DGBWavefunction
    def __init__(self, energies=None, wavefunctions=None, hamiltonian=None, **opts):
        super().__init__(energies=energies, wavefunctions=wavefunctions, hamiltonian=hamiltonian, **opts) # add all opts
        self.hamiltonian = hamiltonian
        self.centers = self.hamiltonian.centers
        self.alphas = self.hamiltonian.alphas
    def __repr__(self):
        return "{}(num={}, DVR={})".format(
            type(self).__name__,
            len(self),
            self.hamiltonian
        )