"""
Redoes what was originally PyDVR but in the _right_ way using proper subclassing and abstract properties
"""

import abc, numpy as np, scipy.sparse as sp, scipy.interpolate as interp

from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger

__all__ = ["BaseDVR", "DVRResults", "DVRException"]

class BaseDVR(metaclass=abc.ABCMeta):
    """
    Provides the abstract interface for creating a
    convenient runnable DVR that can be cleanly subclassed to provide
    extensions
    """

    def __init__(self,
                 domain=None,
                 divs=None,
                 potential_function=None,
                 logger=None,
                 **base_opts
                 ):
        """
        :param base_opts: base opts to use when running
        :type base_opts:
        """
        self.domain = domain
        base_opts['domain'] = domain
        self.divs = divs
        base_opts['divs'] = divs
        self.potential_function = potential_function
        base_opts['potential_function'] = potential_function
        self.opts = base_opts

        if isinstance(logger, Logger):
            self.logger = logger
        elif logger is True:
            self.logger = Logger()
        elif logger is False or logger is None:
            self.logger = NullLogger()
        else:
            self.logger = Logger(logger)

    def __repr__(self):
        if self.potential_function is not None:
            return "{}({}, pts={}, pot={})".format(
                type(self).__name__,
                self.domain,
                self.divs,
                self.potential_function
            )
        else:
            return "{}({}, pts={}, pot={})".format(
                type(self).__name__,
                self.domain,
                self.divs,
                self.potential_function
            )

    def _logger(self, logger):
        return self.logger if logger is None else logger

    @abc.abstractmethod
    def get_grid(self, domain=None, divs=None, **kwargs):
        raise NotImplementedError("abstract interface")
    def grid(self, domain=None, divs=None, **kwargs):
        if domain is None:
            domain = self.domain
        if divs is None:
            divs = self.divs

        if domain is None:
            raise ValueError("need a value for `domain`")
        if divs is None:
            raise ValueError("need a value for `divs`")

        return self.get_grid(domain=domain, divs=divs, **kwargs)

    @abc.abstractmethod
    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        raise NotImplementedError("abstract interface")
    def handle_kinetic_coupling(self, grid, ke_1D, g, g_deriv, hb=1, logger=None, **kwargs):
        logger = self._logger(logger)
        if g is not None:
            with logger.block(tag="handling kinetic coupling"):
                if g_deriv is None:
                    raise ValueError(
                        "if a function for `g` is supplied, also need a function, `g_deriv` for the second derivative of `g`"
                    )
                # add the average value of `g` across the grid points
                if isinstance(g, (int, float, np.integer, np.floating)):
                    logger.log_print("constant G-matrix element")
                    g_vals = np.full(grid.shape, g)
                    g_deriv_vals = np.zeros(grid.shape)
                else:
                    logger.log_print("variable G-matrix element")
                    try:
                        iter(g)
                    except TypeError:
                        g_vals = g(grid)
                    else:
                        g_vals = np.asanyarray(g)

                    try:
                        iter(g_deriv)
                    except TypeError:
                        g_deriv_vals = g_deriv(grid)
                    else:
                        g_deriv_vals = np.asanyarray(g_deriv)

                g_vals = 1 / 2 * (g_vals[:, np.newaxis] + g_vals[np.newaxis, :])
                g_deriv_vals = (hb ** 2) / 2 * np.diag(g_deriv_vals)
                ke_1D = ke_1D * g_vals + g_deriv_vals
        return ke_1D
    def kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs):

        if grid is None:
            grid = self.grid()

        if g is not None:
            mass = 1

        if mass is None:
            raise ValueError("need a value for the mass")

        ke_1D = self.get_kinetic_energy(grid=grid, mass=mass, g=g, g_deriv=g_deriv, hb=hb, **kwargs)
        ke_1D = self.handle_kinetic_coupling(grid, ke_1D, g, g_deriv, hb=hb, **kwargs)

        return ke_1D

    def real_momentum(self, grid=None, mass=None, hb=1, **kwargs):
        raise NotImplementedError("real momentum needs to be implemented")

    def potential_energy(self,
                         grid=None,
                         potential_function=None,
                         potential_values=None,
                         potential_grid=None,
                         logger=None,
                         **pars
                         ):
        """
        Calculates the potential energy at the grid points based
        on dispatching on the input form of the potential

        :param grid: the grid of points built earlier in the DVR
        :type grid:
        :param potential_function: a function to evaluate the potential energy at the points
        :type potential_function:
        :param potential_values: the values of the potential at the DVR points
        :type potential_values:
        :param potential_grid: a grid of points and values to be interpolated
        :type potential_grid:
        :param pars: ignored keyword arguments
        :type pars:
        :return:
        :rtype:
        """

        if grid is None:
            grid = self.grid()

        if isinstance(grid, tuple) and len(grid) == 2: # for FBR DVRs
            grid, tfs = grid

        if potential_function is None and potential_grid is None and potential_values is None:
            potential_function = self.potential_function

        logger = self._logger(logger)
        if potential_function is not None:
            logger.log_print('evaluating potential function over grid')
            # explicit potential function passed; map over coords
            pf=potential_function

            dim = len(grid.shape)
            if dim > 1:
                npts = np.prod(grid.shape[:-1], dtype=int)
                grid = np.reshape(grid, (npts, grid.shape[-1]))
                pot = sp.diags([pf(grid)], [0])
            else:
                pot = np.diag(pf(grid))
        elif potential_values is not None:
            logger.log_print('constructing potential matrix from diagonal values')
            # array of potential values at coords passed
            dim = len(grid.shape)
            if dim > 1:
                pot = sp.diags([potential_values], [0])
            else:
                pot = np.diag(potential_values)
        elif potential_grid is not None:
            logger.log_print('interpolating potential over grid')
            # TODO: extend to include ND, scipy.griddata

            dim = len(grid.shape)
            if dim > 1:
                dim -= 1
                npts = npts = np.prod(grid.shape[:-1], dtype=int)
                grid = np.reshape(grid, (npts, grid.shape[-1]))

            if dim == 1:
                # use a cubic spline interpolation
                interpolator = lambda g1, g2: interp.interp1d(g1[:, 0], g1[:, 1], kind='cubic')(g2)
            else:
                # use griddata to do a general purpose interpolation
                def interpolator(g, g2):
                    # g is an np.ndarray of potential points and values
                    # g2 is the set of grid points to interpolate them over

                    shape_dim = len(g.shape)
                    if shape_dim == 2:
                        points = g[:, :-1]
                        vals = g[:, -1]
                        return interp.griddata(points, vals, g2)
                    else:
                        # assuming regular structured grid
                        mesh = np.moveaxis(g, 0, shape_dim)
                        points = tuple(np.unique(x) for x in mesh[:-1])
                        vals = mesh[-1]
                        return interp.interpn(points, vals, g2)
            wtf = np.nan_to_num(interpolator(potential_grid, grid))
            pot = sp.diags([wtf], [0])
        else:
            raise DVRException("couldn't construct potential matrix")

        return pot

    def hamiltonian(self, kinetic_energy=None, potential_energy=None, potential_threshold=None, **pars):
        """
        Calculates the total Hamiltonian from the kinetic and potential matrices

        :param kinetic_energy:
        :type kinetic_energy:
        :param potential_energy:
        :type potential_energy: np.ndarray | sp.spmatrix
        :param potential_threshold:
        :type potential_threshold:
        :param pars:
        :type pars:
        :return:
        :rtype:
        """

        if potential_threshold is not None:
            diag = potential_energy.diagonal()
            chops = np.where(diag > 0)
            if len(chops) == 0:
                return kinetic_energy + potential_energy
            chops = chops[0]

            ham = kinetic_energy + potential_energy
            ham[chops, :] = 0
            ham[:, chops] = 0

            return ham

        else:
            return kinetic_energy + potential_energy

    def wavefunctions(self, hamiltonian=None, num_wfns=25, nodeless_ground_state=False, diag_mode=None, logger=None, **pars):
        """
        Calculates the wavefunctions for the given Hamiltonian.
        Doesn't support any kind of pruning based on potential values although that might be a good feature
        to support explicitly in the future

        :param hamiltonian:
        :type hamiltonian:
        :param num_wfns:
        :type num_wfns:
        :param nodeless_ground_state:
        :type nodeless_ground_state:
        :param diag_mode:
        :type diag_mode:
        :param pars:
        :type pars:
        :return:
        :rtype:
        """


        logger = self._logger(logger)
        with logger.block(tag="diagonalizing Hamiltonian"):
            logger.log_print([
                "dimension={dimension}",
                "density={density:.3f}%",
                "mode={mode}"
            ],
                dimension=hamiltonian.shape,
                density=100 * hamiltonian.nnz/np.prod(hamiltonian.shape) if isinstance(hamiltonian, sp.spmatrix) else 100,
                mode=diag_mode
            )
            if isinstance(hamiltonian, sp.spmatrix) and diag_mode == 'dense':
                hamiltonian = hamiltonian.toarray()
            if isinstance(hamiltonian, sp.spmatrix):
                import scipy.sparse.linalg as la
                engs, wfns = la.eigsh(hamiltonian, num_wfns, which='SM')
            else:
                engs, wfns = np.linalg.eigh(hamiltonian)
                if num_wfns is not None:
                    engs = engs[:num_wfns]
                    wfns = wfns[:, :num_wfns]

        if nodeless_ground_state:
            logger.log_print('correcting phases to construct nodeless ground state')
            s = np.sign(wfns[:, 0])
            wfns *= s[:, np.newaxis]
        return engs, wfns

    def run(self, result='wavefunctions', logger=None, grid=None, potential_energy=None, kinetic_energy=None, hamiltonian=None, **opts):
        """
        :return:
        :rtype: DVRResults
        """
        from .Wavefunctions import DVRWavefunctions

        opts = dict(self.opts, **opts)

        try:
            self._opts = self.opts
            self.opts = opts

            logger = self._logger(logger)
            with logger.block(tag="Running DVR"):
                logger.log_print("{dvr}", dvr=self)
                logger.log_print(opts, message_prepper=logger.prep_dict)
                opts['logger'] = logger

                res = DVRResults(parent=self, **self.opts)

                with logger.block(tag="constructing grid"):
                    if grid is None:
                        grid = self.grid(**self.opts)
                    res.grid = grid
                    if result == 'grid':
                        return res


                with logger.block(tag="constructing potential matrix"):
                    if potential_energy is None:
                        potential_energy = self.potential_energy(grid=res.grid, **opts)
                    res.potential_energy = potential_energy
                    if result == 'potential_energy':
                        return res


                with logger.block(tag="constructing kinetic matrix"):
                    if kinetic_energy is None:
                        kinetic_energy = self.kinetic_energy(grid=res.grid, **opts)
                    res.kinetic_energy = kinetic_energy
                    if result == 'kinetic_energy':
                        return res


                with logger.block(tag="building Hamiltonian"):
                    if hamiltonian is None:
                        hamiltonian = self.hamiltonian(
                            kinetic_energy=res.kinetic_energy,
                            potential_energy=res.potential_energy,
                            **opts
                        )
                    res.hamiltonian = hamiltonian
                    if result == 'hamiltonian':
                        return res


                with logger.block(tag="evaluating wavefunctions"):
                    energies, wfn_data = self.wavefunctions(
                        hamiltonian=res.hamiltonian,
                        **opts
                    )
                    wfns = DVRWavefunctions(energies=energies, wavefunctions=wfn_data, grid=res.grid, results=res, **opts)
                    res.wavefunctions = wfns

                return res
        finally:
            self.opts = self._opts

class DVRResults:
    """
    A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension
    """
    def __init__(self,
                 grid=None,
                 kinetic_energy=None,
                 potential_energy=None,
                 hamiltonian=None,
                 wavefunctions=None,
                 parent=None,
                 **opts
                 ):

        # self.parent = None
        self.grid = grid
        self.kinetic_energy = kinetic_energy
        self.potential_energy = potential_energy
        self.parent = parent
        self.wavefunctions = wavefunctions
        self.hamiltonian = hamiltonian
        self.opts = opts

    @property
    def dimension(self):
        dim = len(self.grid.shape)
        if dim > 1:
            dim -= 1
        return dim

    def plot_potential(self, plot_class=None, figure=None, plot_units=None, energy_threshold=None, zero_shift=False, **opts):
        """
        Simple plotting function for the potential.
        Should be updated to deal with higher dimensional cases

        :param plot_class: the graphics class to use for the plot
        :type plot_class: McUtils.Plots.Graphics
        :param opts: plot styling options
        :type opts:
        :return:
        :rtype: McUtils.Plots.Graphics
        """
        from McUtils.Plots import Plot, ContourPlot

        # get the grid for plotting
        MEHSH = self.grid
        dim = self.dimension
        if dim == 1:
            mesh = [MEHSH]
        else:
            mesh = np.moveaxis(MEHSH, dim, 0)
        if plot_class is None:
            if dim == 1:
                plot_class = Plot
            elif dim == 2:
                plot_class = ContourPlot
            else:
                raise DVRException("{}.{}: don't know how to plot {} dimensional potential".format(
                    type(self).__name__,
                    'plot',
                    dim
                ))

        pot = self.potential_energy.diagonal()
        if isinstance(plot_units, str) and plot_units == 'wavenumbers':
            pot = pot * UnitsData.convert("Hartrees", "Wavenumbers")

        if zero_shift:
            pot = pot - np.min(pot)

        if energy_threshold:
            pot = pot.copy()
            pot[pot > energy_threshold] = energy_threshold

        return plot_class(*mesh, pot.reshape(mesh[0].shape), figure=figure, **opts)

class DVRException(Exception):
    """
    Base exception class for working with DVRs
    """


