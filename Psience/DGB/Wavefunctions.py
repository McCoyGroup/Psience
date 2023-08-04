"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""

import numpy as np

from McUtils.Zachary import Mesh
import McUtils.Plots as plt

from Psience.Wavefun import Wavefunction, Wavefunctions

__all__ = ["DGBWavefunctions", "DGBWavefunction"]

class DGBWavefunction(Wavefunction):
    def __init__(self, energy, data, centers=None, alphas=None, inds=None, transformations=None,
                 mass_weighted=None, masses=None,
                 **opts):
        super().__init__(energy, data, **opts)
        if centers is None:
            centers = self.parent.centers
        self.centers = centers
        if alphas is None:
            alphas = self.parent.alphas
        self.alphas = alphas
        if transformations is None and self.parent is not None:
            transformations = self.parent.transformations
        self.transformations = transformations
        if inds is None and self.parent is not None:
            inds = self.parent.hamiltonian.inds
        self.inds = inds
        if mass_weighted is None and self.parent is not None:
            mass_weighted = self.parent.mass_weighted
        self.mass_weighted = mass_weighted
        if self.mass_weighted:
            if masses is None:
                masses = self.parent.masses
            masses = np.asanyarray(masses)
        self.masses = masses

    def get_dimension(self):
        return self.centers.shape[-1]

    def plot(self,
             figure=None,
             domain=None,
             plot_centers=False,
             **opts
             ):

        if domain is None:
            domain = Mesh(self.centers).bounding_box

        fig = super().plot(figure=figure, domain=domain, **opts)
        if plot_centers:
            if self.ndim == 1:
                plt.ScatterPlot(
                    self.centers[:, 0],
                    np.zeros(len(self.centers)),
                    figure=fig
                )
            elif self.ndim == 2:
                plt.ScatterPlot(
                    self.centers[:, 0],
                    self.centers[:, 1],
                    figure=fig
                )
            else:
                raise ValueError("can't plot centers for more than 2D data...?")
        return fig

    def evaluate(self, points):
        points = np.asanyarray(points)
        reshape = None
        if points.ndim == 1:
            points = points[np.newaxis]
        elif points.ndim > 2:
            reshape = points.shape[:-1]
            points.reshape(-1, points.shape[-1])

        centers = self.centers
        if self.mass_weighted:
            points = points * np.sqrt(self.masses[np.newaxis])
            centers = centers * np.sqrt(self.masses[np.newaxis])


        # actual basis evaluation
        c_disps = points[np.newaxis, :, :] - centers[:, np.newaxis, :]
        if self.transformations is not None:
            inv, tfs = self.transformations # not sure if I need to transpose this or not...
            tfs = np.broadcast_to(tfs[:, np.newaxis, :, :], (len(self.centers), len(points)) + tfs.shape[1:])
            c_disps = tfs @ c_disps[:, :, :, np.newaxis]
            c_disps = c_disps.reshape(c_disps.shape[:-1])
        alphas = self.alphas[:, np.newaxis]
        if self.inds is not None:
            c_disps = c_disps[..., self.inds]
            alphas = alphas[..., self.inds]
        normas = (2 * alphas / np.pi) ** (1/4)
        exp_evals = np.exp(-alphas * c_disps**2)

        vals = np.dot(
            self.data,
            np.prod(normas * exp_evals, axis=-1)
        )
        if reshape is not None:
            vals = vals.reshape(reshape)
        return vals

    def marginalize_out(self, dofs):
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

        remaining = np.setdiff1d(np.arange(self.centers.shape[-1]), dofs)
        centers = self.centers[:, remaining]

        masses = self.masses[remaining] if isinstance(self.masses, np.ndarray) else self.masses

        if self.transformations is None:
            scaling = np.power(2 * np.pi, len(dofs)/4) / np.power(np.prod(self.alphas[:, proj_dofs], axis=-1), 1/4)
            alphas = self.alphas[:, remaining]

            return type(self)(
                self.energy,
                self.data * scaling,
                centers=centers,
                alphas=alphas,
                transformations=self.transformations,
                mass_weighted=self.mass_weighted,
                masses=masses
            )
        else:
            alphas = self.alphas
            tfs, inv = self.transformations
            if self.inds is not None:
                alphas = alphas[:, self.inds]
                tfs = tfs[:, :, self.inds]
                inv = inv[:, self.inds, :]
            n = alphas.shape[-1]
            npts = len(alphas)
            diag_covs = np.zeros((npts, n, n))
            diag_inds = (slice(None, None, None),) + np.diag_indices(n)
            diag_covs[diag_inds] = 1 / (2 * alphas)


            scaling = np.power(2 * np.pi, len(dofs)/4) / np.power(np.prod(self.alphas[:, proj_dofs], axis=-1), 1/4)

            tfs = tfs[:, remaining, :]
            inv = inv[:, :, remaining]
            covs = tfs @ diag_covs @ inv
            # covs = covs[:, remaining[np.newaxis, :], remaining[:, np.newaxis]]

            # raise Exception(covs[0])
            # covs[np.abs(covs) < 1e-12] = 0  # numerical garbage can be an issue...

            if np.allclose(covs, covs.transpose(0, 2, 1)):
                two_a_inv, tfs = np.linalg.eigh(covs)  # eigenvalues of inverse tensor...
                inv = tfs.transpose(0, 2, 1)
            else:
                two_a_inv, tfs = np.linalg.eig(covs)  # eigenvalues of inverse tensor...

                comp_part = np.imag(two_a_inv)
                if np.sum(comp_part) > 0:
                    raise ValueError(
                        "complex alphas obtained... ({})".format(np.sum(comp_part))
                    )

                two_a_inv = np.real(two_a_inv)
                tfs = np.real(tfs)
                inv_sort = np.argsort(np.abs(two_a_inv), axis=-1)
                idx = np.arange(two_a_inv.shape[0])[:, np.newaxis]
                two_a_inv = two_a_inv[idx, inv_sort]

                idx = np.arange(two_a_inv.shape[0])[:, np.newaxis, np.newaxis]
                cdx = np.arange(two_a_inv.shape[1])[np.newaxis, :, np.newaxis]
                inv_sort = inv_sort[:, np.newaxis, :]

                tfs = tfs[idx, cdx, inv_sort]
                inv = np.linalg.inv(tfs)


            return type(self)(
                self.energy,
                self.data * scaling,
                centers=centers,
                alphas=1/(2*two_a_inv),
                transformations=(tfs, inv),
                mass_weighted=self.mass_weighted,
                masses=masses
            )

class DGBWavefunctions(Wavefunctions):
    wavefunction_class = DGBWavefunction
    def __init__(self, energies=None, wavefunctions=None, hamiltonian=None, **opts):
        super().__init__(energies=energies, wavefunctions=wavefunctions, hamiltonian=hamiltonian, **opts) # add all opts
        self.hamiltonian = hamiltonian
        self.centers = self.hamiltonian.centers
        self.alphas = self.hamiltonian.alphas
        self.transformations = self.hamiltonian.transformations
        self.mass_weighted = self.hamiltonian.mass_weighted
        self.masses = self.hamiltonian.masses

    def expectation(self, op,
                    expansion_degree=None,
                    quadrature_degree=None,
                    expansion_type=None,
                    other=None
                    ):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction | np.ndarray
        :param op:
        :type op:
        :return:
        :rtype:
        """

        if other is None:
            other = self

        if self.hamiltonian is None:
            from .DGB import DGB
            self.hamiltonian = DGB(
                self.centers,
                potential_function=None,
                alphas=self.alphas,
                transformations=self.transformations,
                min_singular_value=None,
                num_svd_vectors=None,
                clustering_radius=None,
                expansion_type=expansion_type,
                expansion_degree=expansion_degree,
                quadrature_degree=quadrature_degree
            )

        if other.hamiltonian is not self.hamiltonian:
            raise ValueError("mismatch in DGBs between {} and {}".format(
                self, other
            ))

        op_mat = self.hamiltonian.evaluate_multiplicative_operator(
            op,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            quadrature_degree=quadrature_degree
        )
        w1 = self.wavefunctions.T
        w2 = other.wavefunctions
        if op_mat.ndim > 2: # lots of shape fuckery to get broacasting right
            op_mat = np.moveaxis(np.moveaxis(op_mat, 0, -1), 0, -1)
            for _ in range(op_mat.ndim - 2):
                w1 = np.expand_dims(w1, 0)
                w2 = np.expand_dims(w2, 0)
            res = np.moveaxis(np.moveaxis(w1 @ op_mat @ w2, -1, 0), -1, 0)
        else:
            res = w1 @ op_mat @ w2

        return res

    def __repr__(self):
        return "{}(num={}, DVR={})".format(
            type(self).__name__,
            len(self),
            self.hamiltonian
        )