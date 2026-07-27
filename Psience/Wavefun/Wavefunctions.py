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
        """
        **LLM Docstring**

        Store the core data defining a single wavefunction: its energy, the underlying representation data (basis-specific), and optional linkage back to a parent `Wavefunctions` collection and index.

        :param energy: the wavefunction's energy
        :type energy: float
        :param data: the representation-specific data defining the wavefunction (e.g. expansion coefficients, grid values)
        :type data: Any
        :param parent: the `Wavefunctions` collection this wavefunction belongs to, if any
        :type parent: Wavefunctions | None
        :param index: this wavefunction's index within its parent collection
        :type index: Any | None
        :param opts: extra representation-specific options
        :type opts: dict
        :return: None
        :rtype: None
        """
        self.energy = energy
        self.data   = data
        self.parent = parent
        self.index = index
        self.opts   = opts

    @abc.abstractmethod
    def get_dimension(self):
        """
        **LLM Docstring**

        Abstract hook for the number of degrees of freedom (dimensionality) this wavefunction is defined over. Concrete subclasses must implement this.

        :return: never returns on the base class
        :rtype: int
        :raises NotImplementedError: always, on the base class
        """
        raise NotImplementedError("abstract base method")
    @property
    def ndim(self):
        """
        **LLM Docstring**

        The number of degrees of freedom this wavefunction is defined over, via `get_dimension`.

        :return: the dimensionality
        :rtype: int
        """
        return self.get_dimension()

    @classmethod
    def prep_plot_grid(cls,
                       domain,
                       plot_points=100,
                       domain_padding=None
                       ):
        """
        **LLM Docstring**

        Build a flattened grid of evaluation points spanning one or more coordinate domains, for use in plotting a wavefunction, optionally padding each domain first.

        :param domain: the coordinate domain(s) to grid over, either a single `(min, max)` pair or a list of such pairs (one per dimension)
        :type domain: tuple | list[tuple]
        :param plot_points: the number of grid points per dimension, either a single integer (applied to every dimension) or a per-dimension list
        :type plot_points: int | list[int]
        :param domain_padding: extra padding to add to each domain before gridding; a single number is applied symmetrically to every dimension, or an explicit `[low, high]` (or per-dimension array of such) can be given
        :type domain_padding: float | Iterable | None
        :return: the flattened grid of evaluation points, shape `(npoints, ndim)`
        :rtype: np.ndarray
        """
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
        """
        **LLM Docstring**

        Compute the overlap `<self|other>` between this wavefunction and another, via `expectation` with the identity operator.

        :param other: the wavefunction to compute the overlap with
        :type other: Wavefunction
        :return: the overlap value
        :rtype: float | complex
        """
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
        """
        **LLM Docstring**

        Build a collection of wavefunctions sharing a common representation, storing their energies, representation data, optional state indices, the concrete `Wavefunction` subclass used to build individual wavefunctions on demand, and an optional dipole function used for spectra.

        :param energies: the energies of each wavefunction in the collection
        :type energies: np.ndarray | None
        :param wavefunctions: the representation-specific data for the collection (e.g. a matrix of expansion coefficients)
        :type wavefunctions: Any | None
        :param indices: explicit state indices this collection represents, if a subset of a larger space
        :type indices: Any | None
        :param wavefunction_class: the `Wavefunction` subclass used to build individual wavefunctions; defaults to `cls.wavefunction_class`
        :type wavefunction_class: type | None
        :param dipole_function: a function providing the dipole operator representation, used by `get_spectrum`
        :type dipole_function: callable | None
        :param opts: extra representation-specific options
        :type opts: dict
        :return: None
        :rtype: None
        """
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
        """
        **LLM Docstring**

        Build the full keyword-argument dict that would be used to reconstruct this collection with the given fields overridden, defaulting each unspecified field to this object's own current value.

        :param energies: replacement energies; defaults to `self.energies`
        :type energies: np.ndarray | None
        :param wavefunctions: replacement representation data; defaults to `self.wavefunctions`
        :type wavefunctions: Any | None
        :param wavefunction_class: replacement wavefunction class; defaults to `self.wavefunction_class`
        :type wavefunction_class: type | None
        :param indices: replacement state indices; defaults to `self.indices`
        :type indices: Any | None
        :param dipole_function: replacement dipole function; defaults to `self.dipole_function`
        :type dipole_function: callable | None
        :param opts: an explicit options dict to use instead of merging `self.opts` with `rem_opts`
        :type opts: dict | None
        :param rem_opts: extra options merged with `self.opts` when `opts` isn't given directly
        :type rem_opts: dict
        :return: the full modification dict
        :rtype: dict
        """
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
        """
        **LLM Docstring**

        Build a new `Wavefunctions` collection (of the same concrete subclass) with the given fields overridden, via `get_modification_dict`.

        :param energies: replacement energies
        :type energies: np.ndarray | None
        :param wavefunctions: replacement representation data
        :type wavefunctions: Any | None
        :param wavefunction_class: replacement wavefunction class
        :type wavefunction_class: type | None
        :param indices: replacement state indices
        :type indices: Any | None
        :param dipole_function: replacement dipole function
        :type dipole_function: callable | None
        :param opts: an explicit options dict to use instead of merging `self.opts` with `rem_opts`
        :type opts: dict | None
        :param rem_opts: extra options merged with `self.opts`
        :type rem_opts: dict
        :return: the new, modified collection
        :rtype: Wavefunctions
        """
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
        """
        **LLM Docstring**

        Retrieve one or more wavefunctions from the collection: a scalar index returns a single `Wavefunction` object, while any other index returns a new `Wavefunctions` collection restricted to the selected states.

        :param which: the state index/indices to retrieve
        :type which: int | object
        :return: the selected wavefunction, or a restricted collection
        :rtype: Wavefunction | Wavefunctions
        """
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
        """
        **LLM Docstring**

        The number of wavefunctions in the collection.

        :return: the number of states
        :rtype: int
        """
        return len(self.energies)
    def __iter__(self):
        """
        **LLM Docstring**

        Iterate over the individual `Wavefunction` objects in the collection, in index order.

        :return: an iterator over the collection's wavefunctions
        :rtype: Iterator[Wavefunction]
        """
        for i in range(len(self)):
            yield self.__getitem__(i)

    def frequencies(self, start_at=0):
        """
        **LLM Docstring**

        Compute the transition frequencies from a reference state to every other state in the collection, as energy differences relative to the reference.

        :param start_at: the index of the reference state
        :type start_at: int
        :return: the transition frequencies (energy differences), excluding the reference state itself
        :rtype: np.ndarray
        """
        return np.concatenate([self.energies[:start_at], self.energies[1+start_at:]]) - self.energies[start_at]

    def get_spectrum(self,
                     dipole_function=None,
                     *,
                     start_at=0,
                     **options
                     ):
        """
        **LLM Docstring**

        Build a discrete IR spectrum for transitions out of a reference state, computing transition dipole moments via `expectation` against the configured (or given) dipole function.

        :param dipole_function: the dipole operator function to use; defaults to `self.dipole_function`
        :type dipole_function: callable | None
        :param start_at: the index of the reference (typically ground) state to compute transitions from
        :type start_at: int
        :param options: extra options forwarded to `expectation`
        :type options: dict
        :return: the constructed discrete spectrum
        :rtype: DiscreteSpectrum
        :raises ValueError: if no dipole function is available (neither given nor stored on the collection)
        """
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
        """
        **LLM Docstring**

        Compute the overlap matrix `<self_i|other_j>` between every wavefunction in this collection and every wavefunction in `other`, via `expectation` with the identity operator.

        :param other: the collection to compute overlaps with
        :type other: Wavefunctions
        :return: the overlap matrix
        :rtype: np.ndarray
        """
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
            dipole_function=dipole_function,
            **etc
        )

    @property
    def coeffs(self):
        """
        **LLM Docstring**

        The expansion coefficients defining this wavefunction in its basis.

        :return: the coefficient vector
        :rtype: np.ndarray
        """
        return self.data['coeffs']
    @property
    def basis(self):
        """
        **LLM Docstring**

        The basis wavefunctions this wavefunction's expansion is defined over.

        :return: the basis
        :rtype: Wavefunctions | None
        """
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
        """
        **LLM Docstring**

        Return the stored dipole matrix, used as the default dipole function for this wavefunction (ignoring the argument, since the dipole matrix is precomputed).

        :param _: unused (accepted for interface compatibility as a dipole-function callback)
        :type _: Any
        :return: the stored dipole matrix
        :rtype: np.ndarray
        """
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
        """
        **LLM Docstring**

        Build a collection of `MatrixWavefunction`s sharing a common basis, Hamiltonian, and dipole matrix.

        :param energies: the energies of each wavefunction in the collection
        :type energies: np.ndarray
        :param coefficients: the matrix of expansion coefficients, one column per wavefunction
        :type coefficients: np.ndarray
        :param basis: the basis the coefficients are expressed in
        :type basis: Wavefunctions | None
        :param hamiltonian: the Hamiltonian matrix in the same basis
        :type hamiltonian: np.ndarray | None
        :param dipole_matrix: the dipole operator matrix in the same basis, used as the default dipole function
        :type dipole_matrix: np.ndarray | None
        :param dipole_function: an explicit dipole function to use instead of `_get_dipoles`
        :type dipole_function: callable | None
        :param wavefunction_class: the `Wavefunction` subclass to use; defaults to `MatrixWavefunction`
        :type wavefunction_class: type | None
        :param ops: extra options forwarded to the base `Wavefunctions.__init__`
        :type ops: dict
        :return: None
        :rtype: None
        """
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
        """
        **LLM Docstring**

        Return the stored dipole matrix, used as the default dipole function for this collection (ignoring the argument, since the dipole matrix is precomputed).

        :param _: unused (accepted for interface compatibility as a dipole-function callback)
        :type _: Any
        :return: the stored dipole matrix
        :rtype: np.ndarray
        """
        return self.dipole_matrix

    def expectation(self, op, other=None):
        """
        **LLM Docstring**

        Compute the expectation-value matrix of an operator between every wavefunction in this collection and every wavefunction in `other`, by contracting the coefficient matrices against the operator matrix (evaluating `op` against `self` first, if it isn't already given as a raw array).

        :param op: the operator, either a raw matrix (in the shared basis) or a callable that builds one when given `self`
        :type op: np.ndarray | callable
        :param other: the other collection to compute expectation values against; defaults to `self`
        :type other: MatrixWavefunctions | None
        :return: the expectation-value matrix
        :rtype: np.ndarray
        """
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