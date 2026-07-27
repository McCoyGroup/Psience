import functools

import numpy as np
from ..VSCF import GridSCF
from McUtils.Scaffolding import ParameterManager

from .BaseDVR import BaseDVR
from .Wavefunctions import DVRWavefunctions
from .DirectProduct import DirectProductDVR
from .FiniteBasisDVR import WavefunctionBasisDVR

__all__ = [
    "SelfConsistentDVR",
    "PotentialOptimizedDVR"
]

class SCFWavefunctionGenerator:
    def __init__(self, dvr_1D:BaseDVR):
        """
        **LLM Docstring**

        Wrap a 1D DVR so it can be re-solved repeatedly with updated potential values (reusing the cached grid/kinetic-energy operator after the first solve), for use as a per-dimension wavefunction generator in a self-consistent-field iteration.

        :param dvr_1D: the 1D DVR to wrap
        :type dvr_1D: BaseDVR
        :return: None
        :rtype: None
        """
        self.dvr = dvr_1D
        self.prev = None
    def __call__(self, pot, **kwargs):
        """
        **LLM Docstring**

        Solve the wrapped 1D DVR for the given potential values, reusing the cached grid and kinetic-energy operator from the previous call (if any) to avoid recomputing them each SCF iteration.

        :param pot: the potential values (on this dimension's DVR grid) to solve with
        :type pot: np.ndarray
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the resulting 1D wavefunctions
        :rtype: DVRWavefunctions
        """
        if self.prev is None:
            res = self.prev = self.dvr.run(potential_values=pot)
        else:
            res = self.dvr.run(
                potential_values=pot,
                grid=self.prev.grid,
                kinetic_energy=self.prev.kinetic_energy
            )
        return res.wavefunctions

class SelfConsistentDVR(GridSCF):
    def __init__(self, base_dvr:"DirectProductDVR", **opts):
        """
        **LLM Docstring**

        Set up a self-consistent-field (SCF) treatment of a multi-dimensional `DirectProductDVR`, precomputing the grid and potential-energy values and building the per-dimension SCF wavefunction generators before delegating to `GridSCF.__init__`.

        :param base_dvr: the multi-dimensional direct-product DVR to run the SCF procedure on
        :type base_dvr: DirectProductDVR
        :param opts: extra options, filtered and forwarded to the base `GridSCF.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        props = ParameterManager(**opts)
        self.base_dvr = base_dvr
        grid, pe = self.create_grid_vals()
        generators = self.create_solvers(grid, pe)
        super().__init__(grid, pe, generators, **props.filter(GridSCF))
    # def initialize(self):
    #     d = super().initialize()
    #     for w in d.wavefunctions:
    #         w[:1].plot()
    #     import McUtils.Plots as plt
    #     plt.DensityPlot(*self.grid.transpose(2, 0, 1), self.vals, plot_style=dict(vmin=0, vmax=1)).show()
    #     return d
    def create_grid_vals(self):
        """
        **LLM Docstring**

        Compute the full grid and potential-energy values for the underlying multi-dimensional DVR, used to initialize the SCF procedure.

        :return: `(grid, pe)` -- the DVR grid and the potential-energy values reshaped to match the grid's spatial shape
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        pot_data = self.base_dvr.run(result='potential_energy')
        grid = pot_data.grid
        pe = pot_data.potential_energy
        if hasattr(pe, 'diagonal'):
            pe = pe.diagonal()
        else:
            pe = np.diag(pe)
        pe = pe.reshape(grid.shape[:-1])
        return grid, pe
    def create_solvers(self, grid, pe):
        """
        **LLM Docstring**

        Build the per-dimension SCF wavefunction generators, first rebinding each 1D DVR's `g`/`g_deriv` kinetic-coupling functions (if callable) so they can be evaluated at the fixed "other-dimension" grid point (from the initial guess) while varying only their own coordinate.

        :param grid: the full multi-dimensional grid
        :type grid: np.ndarray
        :param pe: the potential-energy values on the grid
        :type pe: np.ndarray
        :return: the list of per-dimension `SCFWavefunctionGenerator` objects
        :rtype: list[SCFWavefunctionGenerator]
        """
        # rebind the
        init = self.get_initial_point(pe)
        gp = grid[init]
        for i, d in enumerate(self.base_dvr.dvrs):
            g = d.opts.get('g', None)
            if g is not None:
                try:
                    iter(g)
                except TypeError:
                    g = self._wrap_g(g, gp, i)
                else:
                    g = [
                        [self._wrap_g(el, gp, i) for el in _]
                        for _ in g
                    ]
                d.opts['g'] = g

                gd = d.opts.get('g_deriv', None)
                if gd is not None:
                    try:
                        iter(gd)
                    except TypeError:
                        gd = self._wrap_g(gd, gp, i)
                    else:
                        gd = [self._wrap_g(el, gp, i) for el in gd]
                    d.opts['g_deriv'] = gd

        generators = [
            SCFWavefunctionGenerator(d)
            for d in self.base_dvr.dvrs
        ]

        return generators
    @staticmethod
    def _wrap_g(g, gp, i):
        """
        **LLM Docstring**

        Wrap a per-dimension kinetic-coupling function `g` so that, when evaluated during the SCF procedure for dimension `i`, the other dimensions' coordinates are held fixed at the reference point `gp` while only dimension `i`'s coordinate is allowed to vary over the passed-in 1D grid.

        :param g: the kinetic-coupling function to wrap; passed through unchanged if not callable
        :type g: callable | object
        :param gp: the fixed reference point (coordinates for every dimension) to hold the other dimensions at
        :type gp: np.ndarray
        :param i: the dimension index that should vary
        :type i: int
        :return: the wrapped function, or `g` unchanged if it wasn't callable
        :rtype: callable | object
        """
        if not callable(g):
            return g
        @functools.wraps(g)
        def eval_g(grid, g=g):
            """
            **LLM Docstring**

            Evaluate the enclosing kinetic-coupling function `g` at a full multi-dimensional grid built by broadcasting the fixed reference point `gp` across every dimension except `i`, which is set to the given 1D grid.

            :param grid: the 1D grid of values for dimension `i`
            :type grid: np.ndarray
            :param g: the underlying kinetic-coupling function, captured from the enclosing scope
            :type g: callable
            :return: the evaluated kinetic-coupling values
            :rtype: np.ndarray
            """
            full_ars = [np.full_like(grid, x) for x in gp]
            full_ars[i] = grid
            grid = np.concatenate(
                [np.expand_dims(a, -1) for a in full_ars],
                axis=1
            )
            return g(grid)
        return eval_g

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation showing the class name and the wrapped base DVR.

        :return: string of the form `ClassName(base_dvr)`
        :rtype: str
        """
        return "{}({})".format(
            type(self).__name__,
            self.base_dvr
        )

class PotentialOptimizedDVR(DirectProductDVR):
    def __init__(self,
                 wfns_1D:'Iterable[DVRWavefunctions]',
                 **base_opts
                 ):
        """
        **LLM Docstring**

        Build a multi-dimensional DVR whose per-dimension bases are potential-optimized (built from a given set of 1D wavefunctions rather than a fixed grid), via `WavefunctionBasisDVR` and `DirectProductDVR`.

        :param wfns_1D: the per-dimension 1D wavefunctions to use as the optimized basis
        :type wfns_1D: Iterable[DVRWavefunctions]
        :param base_opts: extra options forwarded to the base `DirectProductDVR.__init__`
        :type base_opts: dict
        :return: None
        :rtype: None
        """
        # base_opts = {k:base_opts[k] for k in base_opts.keys() - {'mass', 'g', 'g_deriv', 'include_kinetic_coupling'}}
        super().__init__(
            [WavefunctionBasisDVR(w) for w in wfns_1D],
            **base_opts
        )

    @classmethod
    def from_minimum(cls, base_dvr:"DirectProductDVR|SelfConsistentDVR", **opts):
        """
        **LLM Docstring**

        Build a `PotentialOptimizedDVR` using the SCF wavefunctions computed from an initial (unconverged, single-iteration) guess at the potential minimum, wrapping `base_dvr` in a `SelfConsistentDVR` first if it isn't already one.

        :param base_dvr: the base multi-dimensional DVR (or an already-built `SelfConsistentDVR`) to optimize the basis from
        :type base_dvr: DirectProductDVR | SelfConsistentDVR
        :param opts: extra options, overriding the base DVR's own stored options, forwarded to the constructor
        :type opts: dict
        :return: the potential-optimized DVR
        :rtype: PotentialOptimizedDVR
        """
        if not isinstance(base_dvr, SelfConsistentDVR):
            base_dvr = SelfConsistentDVR(base_dvr)
        wfns = base_dvr.initialize().wavefunctions
        return cls(
            wfns,
            **dict(
                base_dvr.base_dvr.opts,
                **opts
            )
        )

    @classmethod
    def from_scf(cls, scf_dvr:"DirectProductDVR|SelfConsistentDVR", wfns=None, **opts):
        """
        **LLM Docstring**

        Build a `PotentialOptimizedDVR` using the wavefunctions from a fully-converged SCF run (or an explicitly supplied set of wavefunctions), wrapping `scf_dvr` in a `SelfConsistentDVR` first if it isn't already one.

        :param scf_dvr: the base multi-dimensional DVR (or an already-built `SelfConsistentDVR`) to run/use for the optimized basis
        :type scf_dvr: DirectProductDVR | SelfConsistentDVR
        :param wfns: explicit wavefunctions to use instead of running the SCF procedure
        :type wfns: Iterable[DVRWavefunctions] | None
        :param opts: extra options, overriding the base DVR's own stored options, forwarded to the constructor
        :type opts: dict
        :return: the potential-optimized DVR
        :rtype: PotentialOptimizedDVR
        """
        if not isinstance(scf_dvr, SelfConsistentDVR):
            scf_dvr = SelfConsistentDVR(scf_dvr)
        if wfns is None:
            wfns = scf_dvr.run().wavefunctions
        return cls(
            wfns,
            **dict(
                scf_dvr.base_dvr.opts,
                **opts
            )
        )