"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""

import numpy as np

from McUtils.Zachary import Mesh
import McUtils.Plots as plt

from Psience.Wavefun import Wavefunction, Wavefunctions
from .Components import DGBGaussians, DGBCartesians

__all__ = ["DGBWavefunctions", "DGBWavefunction"]

class DGBWavefunction(Wavefunction):
    def __init__(self, energy, data, gaussians:DGBGaussians = None, **opts):
        super().__init__(energy, data, **opts)
        if gaussians is None:
            gaussians = self.parent.gaussians
        self.gaussians = gaussians

    def get_dimension(self):
        return self.gaussians.coords.shape[-1]

    def plot(self,
             figure=None,
             domain=None,
             plot_centers=False,
             **opts
             ):

        centers = self.gaussians.coords.centers
        if domain is None:
            domain = Mesh(centers).bounding_box

        fig = super().plot(figure=figure, domain=domain, **opts)
        if plot_centers:
            if self.ndim == 1:
                plt.ScatterPlot(
                    centers[:, 0],
                    np.zeros(len(centers)),
                    figure=fig
                )
            elif self.ndim == 2:
                plt.ScatterPlot(
                    centers[:, 0],
                    centers[:, 1],
                    figure=fig
                )
            else:
                raise ValueError("can't plot centers for more than 2D data...?")
        return fig

    def plot_cartesians(self,
                        xyz_sel=None,
                        atom_sel=None,
                        figure=None,
                        plot_centers=False,
                        atom_styles=None,
                        **plot_styles
                        ):
        if not isinstance(self.gaussians.coords, DGBCartesians):
            raise ValueError("can't plot projections onto Cartesian coordinates for coords of type {}".format(
                type(self.gaussians.coords).__name__
            ))

        _, natoms, nxyz = self.gaussians.coords.cart_shape

        full_spec = np.arange(natoms*nxyz).reshape((natoms, nxyz))
        if (xyz_sel is None and nxyz > 2) or (xyz_sel is not None and len(xyz_sel) > 2):
            raise ValueError("can't plot 3D Gaussians")
        all_ats = np.arange(natoms)
        if atom_sel is None:
            atom_sel = all_ats

        if atom_styles is None:
            atom_styles = [{} for _ in atom_sel]

        for i,atom in enumerate(atom_sel):
            rem = np.setdiff1d(all_ats, [atom])
            if xyz_sel is None:
                subspec = full_spec[rem, :].flatten()
            else:
                subspec = full_spec[np.ix_(rem, xyz_sel)].flatten()

            proj = self.marginalize_out(subspec)  # what we're projecting _out_

            ps = dict(plot_styles, **atom_styles[i])

            figure = proj.plot(
                figure=figure,
                plot_centers=plot_centers,
                plotter=plt.TriContourLinesPlot if proj.ndim == 2 else None,
                # levels=np.linspace(-max_val, max_val, 16),
                # domain=[[-2, 2], [-2, .2]],
                # cmap='RdBu'
                **ps
            )

        return figure

    def evaluate(self, points):
        # raise NotImplementedError(...)
        points = np.asanyarray(points)
        reshape = None
        if points.ndim == 1:
            points = points[np.newaxis]
        elif points.ndim > 2:
            reshape = points.shape[:-1]
            points.reshape(-1, points.shape[-1])

        centers = self.gaussians.coords.centers
        # if self.mass_weighted:
        #     points = points * np.sqrt(self.masses[np.newaxis])
        #     centers = centers * np.sqrt(self.masses[np.newaxis])


        # actual basis evaluation
        c_disps = points[np.newaxis, :, :] - centers[:, np.newaxis, :]
        tfs = self.gaussians.transformations
        if tfs is not None:
            tfs = np.broadcast_to(tfs[:, np.newaxis, :, :], (len(tfs), len(points)) + tfs.shape[1:])
            c_disps = tfs @ c_disps[:, :, :, np.newaxis]
            c_disps = c_disps.reshape(c_disps.shape[:-1])
        alphas = self.gaussians.alphas[:, np.newaxis]
        # if self.inds is not None:
        #     c_disps = c_disps[..., self.inds]
        #     alphas = alphas[..., self.inds]
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

        scaling, subgaussians = self.gaussians.marginalize_out(dofs)

        return type(self)(
            self.energy,
            self.data * scaling,
            subgaussians
            )

class DGBWavefunctions(Wavefunctions):
    wavefunction_class = DGBWavefunction
    def __init__(self, energies=None, wavefunctions=None, hamiltonian=None,  **opts):
        super().__init__(energies=energies, wavefunctions=wavefunctions, hamiltonian=hamiltonian, **opts) # add all opts
        self._hamiltonian = hamiltonian
        self.gaussians = self.hamiltonian.gaussians

    @property
    def hamiltonian(self):
        # if self._hamiltonian is None:
        #     from .DGB import DGB
        #     self._hamiltonian = DGB(
        #         self.centers,
        #         potential_function=None,
        #         alphas=self.alphas,
        #         transformations=self.transformations,
        #         min_singular_value=None,
        #         num_svd_vectors=None,
        #         clustering_radius=None
        #     )
        return self._hamiltonian
    def operator_representation(self,
                                op,
                                expansion_degree=None,
                                quadrature_degree=None,
                                expansion_type=None,
                                ):
        return self.hamiltonian.evaluate_multiplicative_operator(
            op,
            expansion_degree=expansion_degree,
            expansion_type=expansion_type,
            quadrature_degree=quadrature_degree
        )

    def expectation(self, op,
                    expansion_degree=None,
                    quadrature_degree=None,
                    expansion_type=None,
                    embed=True,
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

        if other.hamiltonian is not self.hamiltonian:
            raise ValueError("mismatch in DGBs between {} and {}".format(
                self, other
            ))

        if embed:
            op = self.gaussians.coords.embed_function(op)

        op_mat = self.operator_representation(
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

    def localize(self,
                 criterion,
                 which=None
                 ):
        """
        Find a transformation that maximally localizes the wavefunctions in the Boys' sense
        by minimizing <r^2> - <r>^2 over unitary transformations

        :param criterion:
        :param which:
        :return:
        """
        raise NotImplementedError('woof...')

        ndim = self.centers.shape[0]
        def r2_func(centers, deriv_order=2):
            c2 = centers**2

        r2 = self.operator_representation(

        )
        r = self.hamiltonian.S


    def __repr__(self):
        return "{}(num={}, DVR={})".format(
            type(self).__name__,
            len(self),
            self.hamiltonian
        )