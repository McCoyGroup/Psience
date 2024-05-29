"""
Provides a DVRWavefunction class that inherits from the base Psience wavefunction
"""

import numpy as np

from McUtils.Zachary import Mesh
import McUtils.Plots as plt

from Psience.Wavefun import Wavefunction, Wavefunctions
from .Gaussians import DGBGaussians
from .Coordinates import DGBCartesians, DGBWatsonModes

__all__ = ["DGBWavefunctions", "DGBWavefunction"]

__reload_hook__ = ['.Components']

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
        if isinstance(plot_centers, dict):
            plot_centers_opts = plot_centers.copy()
            plot_centers = True
        else:
            plot_centers_opts = {}
        if 'figure' not in plot_centers_opts:
            plot_centers_opts['figure'] = fig
        if plot_centers:
            if self.ndim == 1:
                plt.ScatterPlot(
                    centers[:, 0],
                    np.zeros(len(centers)),
                    **plot_centers_opts
                )
            elif self.ndim == 2:
                plt.ScatterPlot(
                    centers[:, 0],
                    centers[:, 1],
                    **plot_centers_opts
                )
            else:
                raise ValueError("can't plot centers for more than 2D data...?")
        return fig

    def to_cartesian_wavefunction(self):
        """
        Projects the wavefunction back to Cartesians
        :return:
        """
        if isinstance(self.gaussians.coords, DGBCartesians):
            return self

        new_gauss = self.gaussians.as_cartesians()
        return type(self)(
            self.energy,
            self.data,
            new_gauss,
            **self.opts
        )

    def plot_cartesians(self,
                        xyz_sel=None,
                        *,
                        atom_sel=None,
                        figure=None,
                        plot_centers=False,
                        atom_styles=None,
                        **plot_styles
                        ):

        if isinstance(self.gaussians.coords, DGBWatsonModes):
            return self.to_cartesian_wavefunction().plot_cartesians(
                xyz_sel=xyz_sel,
                atom_sel=atom_sel,
                figure=figure,
                plot_centers=plot_centers,
                atom_styles=atom_styles,
                **plot_styles
            )
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

        if xyz_sel is not None:
            xyz_remv = np.setdiff1d(np.arange(nxyz), xyz_sel)
        else:
            xyz_remv = None
        for i,atom in enumerate(atom_sel):
            rem = np.setdiff1d(all_ats, [atom])
            subspec = full_spec[rem, :].flatten()
            if xyz_sel is not None:
                subspec = np.sort(
                    np.concatenate([subspec, full_spec[atom][xyz_remv]])
                )
            proj = self.marginalize_out(subspec)  # what we're projecting _out_

            ps = dict(plot_styles, **atom_styles[i])
            ps['plotter'] = ps.get(
                'plotter',
                plt.TriContourLinesPlot if proj.ndim == 2 else None
            )

            figure = proj.plot(
                figure=figure,
                plot_centers=plot_centers,
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
            tfs, inv = tfs
            tfs = np.broadcast_to(tfs[:, np.newaxis, :, :], (len(tfs), len(points)) + tfs.shape[1:])
            c_disps = tfs.transpose((0, 1, 3, 2)) @ c_disps[:, :, :, np.newaxis]
            c_disps = c_disps.reshape(c_disps.shape[:-1])

        alphas = self.gaussians.alphas[:, np.newaxis, :]
        # if self.inds is not None:
        #     c_disps = c_disps[..., self.inds]
        #     alphas = alphas[..., self.inds]
        normas = np.prod((2 * alphas / np.pi) ** (1/4), axis=-1)
        exp_evals = np.prod(np.exp(-alphas * c_disps**2), axis=-1)

        momenta = self.gaussians.momenta
        if momenta is not None:
            points = points[np.newaxis, :, :, np.newaxis]
            if tfs is not None:
                points = tfs.transpose((0, 1, 3, 2)) @ points
            momenta = momenta[:, np.newaxis, np.newaxis, :]
            correlation = momenta @ points
            correlation = correlation.reshape(correlation.shape[:2])
            exp_evals *= np.cos(correlation)
            sum_term = np.sum(
                (momenta.reshape(alphas.shape)**2) / (2 * alphas),
                axis=-1
            )
            cos_corr = momenta @ centers[:, np.newaxis, :, np.newaxis]
            cos_corr = np.reshape(cos_corr, cos_corr.shape[:2])
            normas /= np.sqrt((1 + np.exp(-sum_term)*np.cos(cos_corr))/2)

        vals = np.dot(
            self.data,
            normas * exp_evals
        )
        if reshape is not None:
            vals = vals.reshape(reshape)
        return vals

    def marginalize_out(self, dofs, rescale=True) -> 'DGBWavefunction':
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
    def as_cartesian_wavefunction(self):
        """
        Projects the wavefunction back to Cartesians
        :return:
        """
        if isinstance(self.gaussians.coords, DGBCartesians):
            return self
        new_ham = self.hamiltonian.as_cartesian_dgb()
        new = type(self)(
            self.energies,
            self.wavefunctions,
            hamiltonian=new_ham
        )

        return new

        # raise Exception(
        #     self[0].evaluate(
        #         self.gaussians.coords.centers[:5]
        #     ),
        #     new[0].evaluate(
        #         new.gaussians.coords.centers[:5]
        #     )
        # )

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
                                embed=True,
                                expansion_degree=None,
                                quadrature_degree=None,
                                expansion_type=None,
                                ):
        return self.hamiltonian.evaluate_multiplicative_operator(
            op,
            embed=embed,
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

        op_mat = self.operator_representation(
            op,
            embed=embed,
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