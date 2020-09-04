import os, numpy as np

class DVR:
    '''
    This is a manager class for working with DVRs.
    It uses files from which it loads the spec data and methods, which are placed in the `Classes` subfolder.
    '''

    loaded_DVRs = {}  # for storing DVRs loaded from file
    dvr_dir = os.path.join(os.path.dirname(__file__), "Classes")

    def __init__(self, dvr_file="ColbertMiller1D", **kwargs):

        self.params = kwargs  # these are the global parameters passed to all methods
        if dvr_file is None:
            self._dvr = None
        else:
            self._dvr = self.load_dvr(dvr_file)

    @classmethod
    def dvr_file(self, dvr):
        """
        Locates the file used to initialize a DVR

        :param dvr: name of the DVR type
        :type dvr: str
        :return: config file to load
        :rtype: str
        """

        if os.path.exists(dvr):
            dvr_file = dvr
        else:
            dvr_file = os.path.join(self.dvr_dir, dvr+".py")

        if not os.path.exists(dvr_file):
            raise DVRException("couldn't load DVR "+dvr)

        return dvr_file

    @classmethod
    def _load_env(cls, file, env):
        """
        Fills in necessary parameters in an env from a module
        :param env:
        :type env:
        :return:
        :rtype:
        """

        cls.loaded_DVRs[file] = env
        for k in ('grid', 'kinetic_energy', 'potential_energy', 'hamiltonian', 'wavefunctions'):
            if 'grid' not in env:
                raise DVRException("{}.{}: DVR class '{}' didn't export property '{}' (in file {})".format(
                    cls.__name__,
                    'load_dvr',
                    env['class'],
                    k,
                    env['file']
                ))

        return env

    @classmethod
    def load_dvr(self, dvr):
        """
        Loads a DVR object from a class file or name

        :param dvr: file or DVR name
        :type dvr: str
        :return: new DVR instance
        :rtype: DVR
        """

        dvr_file = self.dvr_file(dvr)
        if dvr_file not in self.loaded_DVRs:
            # defaults and parameters passed as global
            dvr_class = os.path.splitext(os.path.basename(dvr_file))[0]
            _load_env = {
                'DVR': self,
                'class': dvr_class,
                'file': dvr_file,
                'potential_energy': self._potential_energy,
                'hamiltonian': self._hamiltonian,
                'wavefunctions': self._wavefunctions
                }

            import sys

            dvr_dir = os.path.dirname(dvr_file)
            sys.path.insert(0, os.path.dirname(dvr_dir))
            sys.path.insert(0, dvr_dir)
            try:
                exec("from .Classes.{} import *".format(dvr_class), globals(), _load_env)
            except ImportError:
                exec("from .Classes.{} import *".format(dvr_class), globals(), _load_env)
            finally:  # want to preserve error handling but still close the file
                sys.path.pop(0)

            load = self._load_env(dvr_file, _load_env)

        else:

            load = self.loaded_DVRs[dvr_file]

        return load

    def _dvr_prop(self, prop_name):
        """
        Gets a DVR property by delegating to the loaded class info.

        :param prop_name: the name of the property we want to load
        :type prop_name: str
        :return:
        :rtype:
        """

        if prop_name in self.params:
            prop = self.params[prop_name]
        elif prop_name in self._dvr:
            prop = self._dvr[prop_name]
            if callable(prop):
                prop = prop(**self.params)
        else:
            raise DVRException("no property "+prop_name)

        return prop

    def domain(self):
        """
        :return: the domain for the DVR
        :rtype:
        """
        return self._dvr_prop('domain')

    def divs(self):
        """
        :return: the number of DVR points
        :rtype: tuple(int)
        """
        return self._dvr_prop('divs')

    def grid(self):
        """
        :return: the grid for the DVR
        :rtype: np.ndarray
        """
        return self._dvr_prop('grid')

    def kinetic_energy(self):
        """
        :return: the kinetic energy matrix
        :rtype: np.ndarray
        """
        return self._dvr_prop('kinetic_energy')

    def potential_energy(self):
        """
        :return: the potential energy matrix
        :rtype: np.ndarray
        """
        return self._dvr_prop('potential_energy')

    def hamiltonian(self):
        """
        :return: the total Hamiltonian matrix
        :rtype: np.ndarray
        """
        return self._dvr_prop('hamiltonian')

    def wavefunctions(self):
        """
        :return: the DVR wavefunctions
        :rtype: (np.ndarray, np.ndarray)
        """
        return self._dvr_prop('wavefunctions')

    def _run(self):
        """
        :return:
        :rtype: DVR.Results
        """
        from .Wavefunctions import DVRWavefunctions

        def get_res():
            return self.Results(parent=self, **self.params)
        grid = self.grid()
        self.params['grid'] = grid
        if self.params['result'] == 'grid':
            return get_res()

        pe = self.potential_energy()
        self.params['potential_energy'] = pe
        if self.params['result'] == 'potential_energy':
            return get_res()

        ke = self.kinetic_energy()
        self.params['kinetic_energy'] = ke
        if self.params['result'] == 'kinetic_energy':
            return get_res()

        h = self.hamiltonian()
        self.params['hamiltonian'] = h
        if self.params['result'] == 'hamiltonian':
            return get_res()

        energies, wfn_data = self.wavefunctions()
        self.params['wavefunctions'] = DVRWavefunctions(energies=energies, wavefunctions=wfn_data, **self.params)
        if self.params['result'] == 'wavefunctions':
            return get_res()

        return get_res()

    def run(self, **runpars):
        """
        Runs the DVR with the passed parameters, delegating most calls to the loaded class
        :param runpars:
        :type runpars:
        :return:
        :rtype: DVR.Results
        """

        par = self.params.copy()
        try:
            if 'result' not in self.params:
                self.params['result'] = 'wavefunctions'
            self.params.update(runpars)
            res = self._run()
        finally:
            self.params = par

        return res

    @staticmethod
    def _potential_energy(**pars):
        """
        A default potential energy implementation for reuse
        :param pars: parameters; important keys are potential_values, potential_grid, potential_function
        :type pars:
        :return: potential energy matrix
        :rtype: np.ndarray
        """

        import numpy as np

        if 'potential_function' in pars:
            # explicit potential function passed; map over coords
            pf=pars['potential_function']
            if isinstance(pf, str):
                # these mostly just exist as a few simple test cases
                if pf == 'harmonic_oscillator':
                    k = pars['k'] if 'k' in pars else 1
                    re = pars['re'] if 're' in pars else 0
                    pf = lambda x, k=k, re=re: 1/2*k*(x-re)**2
                elif pf == 'morse_oscillator':
                    de = pars['De'] if 'De' in pars else 10
                    a = pars['alpha'] if 'alpha' in pars else 1
                    re = pars['re'] if 're' in pars else 0
                    pf = lambda x, a=a, de=de, re=re: de*(1-np.exp((-a*(x-re)))**2)
                else:
                    raise DVRException("unknown potential "+pf)
            pot = np.diag([pf(x) for x in pars['grid']])
        elif 'potential_values' in pars:
            # array of potential values at coords passed
            pot = np.diag(pars['potential_values'])
        elif 'potential_grid' in pars:
            # TODO: extend to include ND, scipy.griddata
            import scipy.interpolate as interp, scipy.sparse as sp

            dim = len(pars['grid'].shape)
            if dim > 1:
                dim -= 1

            if dim == 1:
                interpolator = lambda g1, g2: interp.interp1d(g1[:, 0], g1[:, 1], kind='cubic')(g2)
            else:
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
                        mesh = g.transpose(np.roll(np.arange(len(g.shape)), 1))
                        points = tuple(np.unique(x) for x in mesh[:-1])
                        vals = mesh[-1]
                        return interp.interpn(points, vals, g2)
            interp_vals = interpolator(pars['potential_grid'], pars['grid']).flatten()
            pot = sp.diags([interp_vals], [0])
        else:
            raise DVRException("couldn't construct potential matrix")

        return pot

    @staticmethod
    def _hamiltonian(**pars):
        """
        A default Hamiltonian matrix implementation for reuse

        :param pars: parameters; important keys are kinetic_energy and potential_energy
        :type pars:
        :return: Hamiltonian matrix
        :rtype: np.ndarray
        """
        ke = pars['kinetic_energy'] #type:np.ndarray
        pe = pars['potential_energy']  #type:np.ndarray
        return ke + pe

    @staticmethod
    def _wavefunctions(**pars):
        """
        A default wavefunction implementation for reuse
        :param pars: parameters; important keys are hamiltonian, and nodeless_ground_state
        :type pars:
        :return: potential energy matrix
        :rtype: np.ndarray
        """

        engs, wfns = np.linalg.eigh(pars['hamiltonian'])
        if 'nodeless_ground_state' in pars and pars['nodeless_ground_state']:
            s = np.sign(wfns[:, 0])
            wfns *= s[np.newaxis, :]
        return engs, wfns.T

    class Results:
        """
        A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension
        """
        def __init__(self,
                     grid = None,
                     kinetic_energy = None,
                     potential_energy = None,
                     hamiltonian = None,
                     wavefunctions = None,
                     parent = None,
                     **opts
                     ):

            self.parent = None,
            self.grid = grid
            self.kinetic_energy = kinetic_energy
            self.potential_energy = potential_energy
            self.parent = parent
            self.wavefunctions = wavefunctions
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
            from McUtils.Plots import Plot, Plot3D
            import numpy as np

            # get the grid for plotting
            MEHSH = self.grid
            dim = self.dimension
            if dim == 1:
                mesh = [MEHSH]
            else:
                unrolly_polly_OLLY = np.roll(np.arange(len(MEHSH.shape)), 1)
                mesh = MEHSH.transpose(unrolly_polly_OLLY)
            if plot_class is None:
                if dim == 1:
                    plot_class = Plot
                elif dim == 2:
                    plot_class = Plot3D
                else:
                    raise DVRException("{}.{}: don't know how to plot {} dimensional potential".format(
                        type(self).__name__,
                        'plot',
                        dim
                    ))

            pot = self.potential_energy.diagonal()

            return plot_class(*mesh, pot.reshape(mesh[0].shape), **opts)

            #
            #
            #
            # # Plot3D(mesh[0], mesh[1], np.diag(res["potential_energy"].todense()).reshape(mesh[0].shape) ).show()
            # for i in range(10):
            #     ContourPlot(mesh[0], mesh[1], res["wavefunctions"][1][:, i].reshape(mesh[0].shape)).show()
            # # self.assertIsInstance(res[0], np.ndarray)

class DVRException(Exception):
    '''An Exception in a DVR '''
    pass
