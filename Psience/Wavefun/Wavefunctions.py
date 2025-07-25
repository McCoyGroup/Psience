"""
Provides very general support for an abstract wavefunction object
Allows different methods to provide their own concrete implementation details
"""
import abc
from abc import *
import numpy as np

import McUtils.Devutils as dev
import McUtils.Numputils as nput
from McUtils.Plots import Graphics, Plot, TriContourPlot
from ..Spectra import DiscreteSpectrum

__all__ = [
    "Wavefunction",
    "Wavefunctions",
    "MatrixWavefunction",
    "MatrixWavefunctions",
    "WavefunctionException"
]

__reload_hook__ = ["..Spectra"]

class WavefunctionException(Exception):
    pass

class Wavefunction:
    """Represents a single wavefunction object"""
    def __init__(self, energy, data, parent=None, index=None, **opts):
        self.energy = energy
        self.data   = data
        self.parent = parent
        self.index = index
        self.opts   = opts

    @abc.abstractmethod
    def get_dimension(self):
        raise NotImplementedError("abstract base method")
    @property
    def ndim(self):
        return self.get_dimension()

    @classmethod
    def prep_plot_grid(cls,
                       domain,
                       plot_points=100,
                       domain_padding=None
                       ):
        if isinstance(domain[0], (int, np.integer, float, np.floating)):
            domain = [domain]
        if isinstance(plot_points, (int, np.integer)):
            plot_points = [plot_points] * len(domain)

        if domain_padding is not None:
            if isinstance(domain_padding, (int, float, np.integer, np.floating)):
                domain_padding = [-domain_padding, domain_padding]
            domain_padding = np.asanyarray(domain_padding)
            if domain_padding.ndim == 1:
                domain_padding = domain_padding[np.newaxis, :]
            domain = np.asanyarray(domain) + domain_padding

        grids = []
        for dom, pts in zip(domain, plot_points):
            grids.append(np.linspace(*dom, pts))
        grid = np.moveaxis(np.array(np.meshgrid(*grids, indexing='xy')), 0, -1).reshape(-1, len(domain))  # vector of points

        return grid

    def plot(self,
             figure=None, domain=None, *, domain_padding=None, grid=None, values=None, plot_points=100,
             which=None,
             index=0, scaling=1, shift='auto', plotter=None, plot_density=False, return_values=False,
             zero_tol=1e-8, contour_levels=None,
             **opts
             ):
        """
        Plots a single wave function on the grid

        :param figure:
        :type figure:
        :param grid:
        :type grid:
        :param index:
        :type index:
        :param scaling:
        :type scaling:
        :param shift:
        :type shift:
        :param opts:
        :type opts:
        :return:
        :rtype:
        """
        if which is None:
            which = index

        if grid is None and domain is None:
            raise ValueError("can't plot a wave function without a specified domain")

        if grid is None:
            grid = self.prep_plot_grid(domain=domain, domain_padding=domain_padding, plot_points=plot_points)

        grid = np.asanyarray(grid)
        if grid.ndim == 1:
            grid = grid[:, np.newaxis]
        elif grid.ndim > 2:
            grid = grid.reshape(-1, grid.shape[-1])
        dim = grid.shape[-1]

        if dim > 2 and plotter is None: # if people want to try, let 'em
            raise ValueError("can't plot data with dimension higher than 2, take a projection first")

        # allows us to scale wave functions independently
        if not isinstance(scaling, (int, float, np.integer, np.floating)):
            scaling = scaling[which]
        if shift is None: shift = 0
        if isinstance(shift, str) and shift == 'auto':
            shift = self.energy
        if not isinstance(shift, (int, float, np.integer, np.floating, str)):
            shift = shift[which]

        if values is None:
            if plot_density:
                values = self.probability_density(grid)
            else:
                values = self.evaluate(grid)
        values[np.abs(values) < zero_tol] = 0.

        values = values * scaling
        if contour_levels is not None and 'levels' not in opts:
            if np.std(values) > 1e-8:
                max_val = np.max(np.abs(values))
                levels = np.linspace(-max_val + shift, max_val + shift, contour_levels)
                opts['levels'] = levels
        values = values + shift

        if plotter is None:
            if dim == 1:
                plotter = Plot
            else:
                plotter = TriContourPlot

        if return_values:
            return {
                'plotter':plotter,
                'grid':grid,
                'values':values,
                'figure':figure,
                'opts':opts
            }
        else:
            return plotter(*np.moveaxis(grid, -1, 0), values, figure=figure, **opts)

    def projection_plot(self,
                        coords,
                        figure=None,
                        **plot_options
                        ):
        """
        A convenience function to plot multiple projections
        on the same set of axes

        Deprecated in favor of `plot_cartesians` for its primary use case (`DGBWavefunctions`)

        :param coords:
        :type coords:
        :param figure:
        :type figure:
        :param plot_options:
        :type plot_options:
        :return:
        :rtype:
        """
        if isinstance(coords, (int, np.integer)):
            coords = [[coords]]
        elif isinstance(coords[0], (int, np.integer)):
            coords = [[c] for c in coords]

        for proj_inds in coords:
            fig = self.project(proj_inds).plot(
                figure=figure,
                **plot_options
            )
            if figure is None:
                figure = fig

        return figure

    @abstractmethod
    def expectation(self, op, other=None):
        """Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunction
        :param op:
        :type op:
        :return:
        :rtype:
        """
        pass
    def overlap(self, other):
        return self.expectation(lambda w:w, other=other)
    @abstractmethod
    def evaluate(self, points):
        """
        Evaluates the current wavefunction

        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base method")
    @property
    def probability_density(self):
        """
        Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        return lambda pts:self.evaluate(pts)**2 # we're assuming real I guess...
    @abstractmethod
    def marginalize_out(self, dofs):
        """
        Integrates out the contributions from the degrees of freedom `dofs`

        :return:
        :rtype: Wavefunction
        """
        raise NotImplementedError("abstract base method")
    def project(self, dofs):
        """
        Computes the projection of the current wavefunction onto a set of degrees
        of freedom, returning a projected wave function object

        :return:
        :rtype: Wavefunction
        """

        dof_complement = np.setdiff1d(np.arange(self.ndim), dofs)
        return self.marginalize_out(dof_complement)

class Wavefunctions:
    """
    An object representing a set of wavefunctions.
    Provides concrete, but potentially inefficient methods for doing all the wavefunction ops.

    """
    wavefunction_class = Wavefunction
    def __init__(self,
                 energies=None, wavefunctions=None,
                 indices=None, wavefunction_class=None,
                 dipole_function=None,
                 **opts):
        self.wavefunctions = wavefunctions
        self.energies = energies
        self.wavefunction_class = self.wavefunction_class if wavefunction_class is None else wavefunction_class
        self.indices = indices
        self.dipole_function = dipole_function
        self.opts = opts

    def get_modification_dict(self,
                              *,
                              energies=None,
                              wavefunctions=None,
                              wavefunction_class=None,
                              indices=None,
                              dipole_function=None,
                              opts=None,
                              **rem_opts
                              ):
        if opts is None:
            opts =dev.merge_dicts(self.opts, rem_opts)
        return dict(
            energies=energies if energies is not None else self.energies,
            wavefunctions=wavefunctions if wavefunctions is not None else self.wavefunctions,
            wavefunction_class=wavefunction_class if wavefunction_class is not None else self.wavefunction_class,
            indices=indices if indices is not None else self.indices,
            dipole_function=dipole_function if dipole_function is not None else self.dipole_function,
            **opts
        )
    def modify(self,
               *,
               energies=None,
               wavefunctions=None,
               wavefunction_class=None,
               indices=None,
               dipole_function=None,
               opts=None,
               **rem_opts
               ):
        mod_dict = self.get_modification_dict(
            energies=energies,
            wavefunctions=wavefunctions,
            wavefunction_class=wavefunction_class,
            indices=indices,
            dipole_function=dipole_function,
            opts=opts,
            **rem_opts
        )
        return type(self)(
            mod_dict.pop('energies'),
            mod_dict.pop('wavefunctions'),
            **mod_dict
        )

    def get_wavefunctions(self, which):
        inds = self.indices
        if inds is None:
            inds = np.arange(len(self.wavefunctions))
        if not isinstance(which, (int, np.integer)):
            return self.modify(
                energies=self.energies[which,],
                wavefunctions=self.wavefunctions[:, which],
                indices=inds[which,]
            )
        else:
            return self.wavefunction_class(
                self.energies[which],
                self.wavefunctions[:, which],
                parent=self,
                index=inds[which],
                **self.opts
            )
    def __getitem__(self, item):
        """Returns a single Wavefunction object"""
        # iter comes for free with this
        return self.get_wavefunctions(item)
    def __len__(self):
        return len(self.energies)
    def __iter__(self):
        for i in range(len(self)):
            yield self.__getitem__(i)

    def frequencies(self, start_at=0):
        return np.concatenate([self.energies[:start_at], self.energies[1+start_at:]]) - self.energies[start_at]

    def get_spectrum(self,
                     dipole_function=None,
                     *,
                     start_at=0,
                     **options
                     ):
        if dipole_function is None: # it's just so convenient to have this on the object...
            dipole_function = self.dipole_function
        if dipole_function is None:
            raise ValueError("a dipole function is required to get a spectrum (none stored in wavefunctions)")
        freqs = self.frequencies(start_at=start_at)
        transition_moments = self.expectation(dipole_function,
                                              other=self[(start_at,)],
                                              **options,
                                              ).reshape(-1, 3)
        transition_moments = np.concatenate([
            transition_moments[:start_at],
            transition_moments[start_at+1:]
        ])
        return DiscreteSpectrum.from_transition_moments(
            freqs,
            transition_moments
        )

    def plot(self, figure=None, graphics_class=None, **opts):
        """Plots all of the wavefunctions on one set of axes

        :param opts:
        :type opts:
        :return:
        :rtype:
        """

        k = "plot_defaults"
        opts = dict(self.opts[k] if k in self.opts else (), **opts)

        if figure == None:
            dim = self.opts['dimension'] if 'dimension' in self.opts else 1
            if graphics_class is None:
                if dim ==1:
                    graphics_class = Graphics
                elif dim == 2:
                    graphics_class = Graphics#Graphics3D
                else:
                    raise WavefunctionException(
                        "{}.{}: don't know how to plot wavefunctions of dimension {}".format(
                            type(self).__name__, 'plot', dim
                        )
                    )
            figure = graphics_class(strict=False, **opts)

        for i, wfn in enumerate(self):
            ind = wfn.index
            if ind is None:
                ind = i
            wfn.plot(figure, which=i, index=ind, **opts)

        return figure

    def expectation(self, op, other=None):
        """
        Computes the expectation value of operator op over the wavefunction other and self

        :param other:
        :type other: Wavefunctions
        :param op:
        :type op:
        :return:
        :rtype:
        """
        if other is None:
            other = self

        res = []
        for wfn in self:
            subres = []
            for ofn in other:
                subres.append(wfn.expectation(op, ofn))
            res.append(subres)
        return np.array(res)
    def overlap(self, other):
        return self.expectation(lambda w:w, other=other)

    def coordinate(self):
        """
        Provides the coordinate operator in the wavefunction basis

        :return:
        :rtype:
        """
        raise NotImplementedError("no coordinate rep implemented for {}".format(self))
    def momentum(self):
        """
        Provides the real part of the representation of the momentum operator in the wavefunction basis

        :return:
        :rtype:
        """
        raise NotImplementedError("no momentum implemented for {}".format(self))
    def laplacian(self):
        """
        Provides the representation of the laplacian in the wavefunction basis

        :return:
        :rtype:
        """
        raise NotImplementedError("no momentum implemented for {}".format(self))
    def kinetic_energy(self):
        """
        Provides the representation of the KE in the wavefunction basis

        :return:
        :rtype:
        """
        raise NotImplementedError("no KE implemented for {}".format(self))

class MatrixWavefunction(Wavefunction):
    """
    Simple wave function that takes a set of expansion coefficients alongside its basis.
    Technically this should be called a _linear expansion wave function_, but
    that was too long for my taste.
    """
    def __init__(self, energy, coefficient_vector, basis=None, dipole_matrix=None, dipole_function=None, **etc):
        """
        :param energy: energy of the wavefunction
        :type energy: float
        :param coefficients: expansion coefficients
        :type coefficients: Iterable[float]
        :param basis_wfns: basis functions for the expansion
        :type basis_wfns: Wavefunctions
        """
        if dipole_function is None:
            dipole_function = self._get_dipoles
        self.dipole_matrix = dipole_matrix
        super().__init__(
            energy,
            {
                'coeffs':coefficient_vector,
                'basis':basis
            },
            dipole_function=dipole_function
            **etc
        )

    @property
    def coeffs(self):
        return self.data['coeffs']
    @property
    def basis(self):
        return self.data['basis']

    def evaluate(self, *args, **kwargs):
        """
        Evaluates the wavecfunction as any other linear expansion.

        :param args: coordinates + any other args the basis takes
        :type args:
        :param kwargs: any keyword arguments the basis takes
        :type kwargs:
        :return: values of the wavefunction
        :rtype:
        """
        if self.basis is None:
            raise ValueError("can't evaluate without basis")
        return np.dot(self.data['coeffs'], self.basis(*args, **kwargs))

    def _get_dipoles(self, _):
        return self.dipole_matrix

    def expect(self, operator):
        """
        Provides the expectation value of the operator `op`.
        Uses the basis to compute the reps and then expands with the expansion coeffs.

        :param operator:
        :type operator:
        :return:
        :rtype:
        """
        return np.dot(np.dot(self.data['coeffs'], operator), self.data['coeffs'])

    def expectation(self, operator, other):
        """
        Computes the expectation value of operator `op` over the wavefunction `other` and `self`.
        **Note**: _the basis of `other`, `self`, and `op` are assumed to be the same_.

        :param op: an operator represented in the basis of the expansion
        :type op: Operator
        :param other: the other wavefunction to expand over
        :type other: ExpansionWavefunction
        :return:
        :rtype:
        """

        return np.dot(np.dot(self.data['coeffs'], operator), other.data['coeffs'])

    def probability_density(self):
        """Computes the probability density of the current wavefunction

        :return:
        :rtype:
        """
        raise NotImplementedError("expansion wave function probability densities not yet implemented")

    def project(self, dofs):
        """
        Computes the projection of the current wavefunction onto a set of degrees
        of freedom

        :return:
        :rtype:
        """
        raise NotImplementedError("expansion wave function projections not yet implemented")

class MatrixWavefunctions(Wavefunctions):
    wavefunctions: np.ndarray
    def __init__(self, energies, coefficients,
                 basis=None,
                 hamiltonian=None,
                 dipole_matrix=None,
                 dipole_function=None, wavefunction_class=None, **ops):
        # self._coeffs = coefficients
        # self._energies = energies
        self.basis = basis
        self.hamiltonian = hamiltonian
        self.dipole_matrix = dipole_matrix
        if dipole_function is None:
            dipole_function = self._get_dipoles
        if wavefunction_class is None:
            wavefunction_class = MatrixWavefunction
        super().__init__(energies, coefficients,
                         dipole_function=dipole_function,
                         wavefunction_class=wavefunction_class,
                         **ops
                         )

    def _get_dipoles(self, _):
        return self.dipole_matrix

    def expectation(self, op, other=None):
        if other is None:
            other = self
        if not nput.is_numeric_array_like(op):
            op = op(self)
        op = np.asanyarray(op)
        return np.tensordot(
            self.wavefunctions,
            np.tensordot(
                other.wavefunctions,
                op,
                axes=[0, 0]
            ),
            axes=[0, 1]
        )

    # def expect(self, op, other):
    #     """
    #     Provides expectation values of the wavefunctions o
    #     :param op:
    #     :type op:
    #     :param other:
    #     :type other:
    #     :return:
    #     :rtype:
    #     """
    #     return NotImplemented


    # def get_wavefunctions(self, which):
    #     energy = self.energies[which]
    #     wfn = self.wavefunctions[which]
    #     return self.wavefunction_class(energy, wfn, self._basis)