"""
Redoes what was originally PyDVR but in the _right_ way using proper subclassing and abstract properties
"""

import abc, numpy as np, scipy.sparse as sp, scipy.interpolate as interp

class AbstractDVR(metaclass=abc.ABCMeta):
    """
    Provides the abstract interface for creating a
    convenient runnable DVR that can be cleanly subclassed to provide
    extensions
    """

    @abc.abstractmethod
    def grid(self, **kwargs):
        raise NotImplementedError("abstact interface")

    @abc.abstractmethod
    def kinetic_energy(self, **kwargs):
        raise NotImplementedError("abstact interface")

    def potential_energy(self, grid=None,
                         potential_function=None,
                         potential_values=None,
                         **pars
                         ):
        """ A default ND potential implementation for reuse"""

        if potential_function is not None:
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
            # array of potential values at coords passed
            dim = len(grid.shape)
            if dim > 1:
                pot = sp.diags([potential_values], [0])
            else:
                pot = np.diag(potential_values)
        elif potential_values is not None:
            # TODO: extend to include ND, scipy.griddata

            grid = pars['grid']
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
            wtf = np.nan_to_num(interpolator(pars['potential_grid'], grid))
            pot = sp.diags([wtf], [0])
        else:
            raise DVRException("couldn't construct potential matrix")

        return pot

    def hamiltonain(self, kinetic_energy=None, potential_energy=None, potential_threshold=None, **pars):
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

    def wavefunctions(self, hamiltonian=None, num_wfns=None, nodeless_ground_state=False, diag_mode=None, **pars):
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

        if isinstance(hamiltonian, sp.spmatrix) and diag_mode == 'dense':
            hamiltonian = hamiltonian.toarray()
        if isinstance(hamiltonian, sp.spmatrix):
            import scipy.sparse.linalg as la
            engs, wfns = la.eigsh(hamiltonian, num_wfns, which='SM')
        else:
            engs, wfns = np.linalg.eigh(hamiltonian)
            if num_wfns is not None:
                engs = engs[:num_wfns]
                wfns = wfns[:, num_wfns]

        if nodeless_ground_state:
            s = np.sign(wfns[:, 0])
            wfns *= s[np.newaxis, :]
        return engs, wfns

    def run(self, result='wavefunctions', **opts):
        """
        :return:
        :rtype: DVR.Results
        """
        from .Wavefunctions import DVRWavefunctions

        res = DVRResults(parent=self, **opts)

        grid = self.grid(**opts)
        res.grid = grid
        if result == 'grid':
            return res

        pe = self.potential_energy(grid=res.grid, **opts)
        res.potential_energy = pe
        if result == 'potential_energy':
            return res

        ke = self.kinetic_energy(grid=res.grid, **opts)
        res.kinetic_energy = ke
        if result == 'kinetic_energy':
            return res

        h = self.hamiltonian(
            kinetic_energy=res.kinetic_energy,
            potential_energy=res.potential_energy,
            **opts
        )
        res.hamiltonian = h
        if result == 'hamiltonian':
            return res

        energies, wfn_data = self.wavefunctions(
            hamiltonian=res.hamiltonian,
            **opts
        )
        wfns = DVRWavefunctions(energies=energies, wavefunctions=wfn_data, **opts)
        res.wavefunctions = wfns

        return res

class DVRException(Exception):
    """
    Base exception class for working with DVRs
    """

class DVRResults:
    """
    A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension
    """
    def __init__(self,
                 grid = None,
                 kinetic_energy=None,
                 potential_energy=None,
                 hamiltonian=None,
                 wavefunctions=None,
                 parent=None,
                 **opts
                 ):

        self.parent=None,
        self.grid=grid
        self.kinetic_energy=kinetic_energy
        self.potential_energy=potential_energy
        self.parent=parent
        self.wavefunctions=wavefunctions
        self.hamiltonian=hamiltonian
        self.opts = opts

    @property
    def dimension(self):
        dim = len(self.grid.shape)
        if dim > 1:
            dim -= 1
        return dim

    def plot_potential(self, plot_class=None, **opts):
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

        return plot_class(*mesh, pot.reshape(mesh[0].shape), **opts)


