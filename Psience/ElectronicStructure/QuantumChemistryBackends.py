"""
Quantum-chemistry backends: GFN-xTB (through `xtb-python`) and PySCF
"""

import re
import warnings

import numpy as np

import McUtils.Devutils as dev

from .Backends import ElectronicStructureBackend, import_or_explain

__all__ = [
    "XTBBackend",
    "PySCFBackend",
]


@ElectronicStructureBackend.register('xtb')
class XTBBackend(ElectronicStructureBackend):
    """
    GFN-xTB energies, gradients, Mulliken charges, and dipoles through the `xtb` python bindings
    (conda-forge: `xtb-python`). Every requested property comes from one single point.

    :param atoms: the atom symbols
    :param method: ``'GFN2-xTB'`` (default), ``'GFN1-xTB'`` or ``'GFN0-xTB'``
    :param charge: the total charge (default 0)
    :param multiplicity: the spin multiplicity (default 1); passed to xtb as ``multiplicity - 1`` unpaired electrons
    :param quiet: silence xtb's output
    :param accuracy: xtb's numerical accuracy parameter (smaller is tighter), if not the default
    :param max_iterations: the SCF iteration limit, if not the default
    :param electronic_temperature: the electronic temperature in Kelvin, if not the default
    :param solvent: an ALPB solvent name (e.g. ``'water'``)
    :param retry_iterations: if the SCF fails to converge, try once more with this iteration limit
        (the electronic temperature and accuracy are never changed behind your back); `None` disables the retry
    """
    distance_units = "BohrRadius"
    default_multiplicity = 1
    property_units = {
        'energy': 'Hartrees',
        'charges': 'ElementaryCharge',
        'dipole': ('ElementaryCharge', 'BohrRadius'),
    }
    analytic_derivative_orders = {
        'energy': 1,
        'charges': 0,
        'dipole': 0,
    }
    property_shapes = {
        'charges': lambda n: (n,),
        'dipole': (3,),
    }

    def __init__(self,
                 atoms,
                 method="GFN2-xTB",
                 charge=None,
                 multiplicity=None,
                 quiet=True,
                 accuracy=None,
                 max_iterations=None,
                 electronic_temperature=None,
                 solvent=None,
                 retry_iterations=1000):
        super().__init__(atoms, charge=charge, multiplicity=multiplicity, quiet=quiet)
        self.method = method
        self.accuracy = accuracy
        self.max_iterations = max_iterations
        self.electronic_temperature = electronic_temperature
        self.solvent = solvent
        self.retry_iterations = retry_iterations
        self.check_electron_count(self.charge, self.multiplicity)

    @classmethod
    def handle_specialization(cls, tag):
        return {'method': tag}

    @classmethod
    def discover_dependencies(cls, **opts):
        """
        The `xtb` python bindings come from conda-forge (``xtb-python``); the PyPI ``xtb`` package needs a
        separately installed libxtb.
        """
        from McUtils.ExternalPrograms import JobDependencies
        deps = JobDependencies()
        deps.add('conda', 'xtb-python', import_name='xtb', reason='xtb python bindings')
        deps.add('conda', 'tblite', '<0.6', reason='xtb-python links tblite; tblite 0.6.0 is ABI-incompatible with the '
                 'current dftd4 build on conda-forge (undefined symbol `smooth_cutoff` when importing xtb)')
        deps.add('conda', 'typing_extensions', import_name='typing_extensions',
                 reason='imported by xtb.interface, but missing from the xtb-python conda package')
        return deps

    def get_calculator(self, coords, charge=None, multiplicity=None):
        """
        Build an xtb `Calculator` at the given Cartesians (Bohr) with this backend's settings.
        """
        charge, multiplicity = self.resolve_state(charge, multiplicity)
        interface = import_or_explain('xtb.interface', 'xtb-python')
        utils = import_or_explain('xtb.utils', 'xtb-python')
        param = utils.get_method(self.method)
        if param is None:
            raise ValueError(f"unknown xTB method {self.method!r} (try 'GFN2-xTB', 'GFN1-xTB', 'GFN0-xTB')")
        calc = interface.Calculator(
            param, self.numbers, coords,
            charge=charge,
            uhf=multiplicity - 1
        )
        if self.quiet:
            calc.set_verbosity('muted')
        if self.accuracy is not None:
            calc.set_accuracy(self.accuracy)
        if self.max_iterations is not None:
            calc.set_max_iterations(self.max_iterations)
        if self.electronic_temperature is not None:
            calc.set_electronic_temperature(self.electronic_temperature)
        if self.solvent is not None:
            solvent = utils.get_solvent(self.solvent)
            if solvent is None:
                raise ValueError(f"unknown xTB solvent {self.solvent!r}")
            calc.set_solvent(solvent)
        return calc

    def get_single_point(self, coords, charge=None, multiplicity=None):
        """
        Run one xTB single point, retrying once with a larger iteration limit if the SCF doesn't converge.
        """
        calc = self.get_calculator(coords, charge=charge, multiplicity=multiplicity)
        try:
            return calc.singlepoint()
        except Exception as e:  # xtb raises `XTBException`, which we can't import without xtb
            if type(e).__name__ != 'XTBException' or self.retry_iterations is None:
                raise
            if self.max_iterations is not None and self.max_iterations >= self.retry_iterations:
                raise RuntimeError(f"xTB ({self.method}) failed: {e}") from e
            calc = self.get_calculator(coords, charge=charge, multiplicity=multiplicity)
            calc.set_max_iterations(self.retry_iterations)
            try:
                return calc.singlepoint()
            except Exception as e2:
                raise RuntimeError(
                    f"xTB ({self.method}) failed even with {self.retry_iterations} iterations: {e2}"
                ) from e2

    def evaluate_properties(self, coords, orders, charge, multiplicity, **opts):
        if (charge, multiplicity) != (self.charge, self.multiplicity):
            self.check_electron_count(charge, multiplicity)
        nstruct, n = coords.shape[:2]
        out = {p: [np.empty((nstruct,) + self.get_property_shape(p))] for p in orders}
        if orders.get('energy', 0) > 0:
            out['energy'].append(np.empty((nstruct, 3 * n)))
        for i, coord in enumerate(coords):
            res = self.get_single_point(coord, charge=charge, multiplicity=multiplicity)
            for p, o in orders.items():
                if p == 'energy':
                    out[p][0][i] = res.get_energy()
                    if o > 0:
                        out[p][1][i] = res.get_gradient().reshape(-1)
                elif p == 'charges':
                    out[p][0][i] = res.get_charges()
                elif p == 'dipole':
                    out[p][0][i] = res.get_dipole()
                else:
                    raise ValueError(f"unknown property {p!r}")
        return out


@ElectronicStructureBackend.register('pyscf')
class PySCFBackend(ElectronicStructureBackend):
    """
    PySCF energies, gradients and (for HF/DFT) Hessians, plus SCF-level charges, dipoles, and
    static polarizabilities. Every requested property comes from one converged calculation.

    :param atoms: the atom symbols
    :param level_of_theory: ``'hf'``/``'rhf'``/``'uhf'``/``'rohf'``, ``'mp2'``, or a DFT functional name
        (``'b3lyp'``, ``'wb97x-d3bj'``, ...; a trailing ``-d3bj``/``-d3zero``/``-d4`` selects a dispersion correction,
        which needs ``pyscf-dispersion``); or a callable ``caller(molecule, **opts)`` returning a finished calculation
    :param basis_set: the basis set (``basis`` is accepted as an alias)
    :param charge: the total charge
    :param multiplicity: the spin multiplicity; PySCF's own ``spin`` (2S) is set to ``multiplicity - 1``.
        Open shells use unrestricted references (``'rohf'`` for restricted open-shell HF).
    :param quiet: silence PySCF's output
    :param molecule_options: extra options for `gto.Mole` (coordinates are always passed in Bohr)
    :param level_of_theory_options: extra options for the calculation (``functional_options``, ``restricted``,
        ``reference`` for MP2, and any SCF attribute such as ``conv_tol`` or ``max_cycle``)
    :param disp: a dispersion correction (``'d3bj'``, ``'d3zero'``, ``'d4'``) for DFT
    :param scf_options: SCF attributes to set on every calculation (``{'conv_tol': 1e-10, 'max_cycle': 200}``)
    :param grid_level: the DFT integration grid level
    :param density_fit: `True`, or an auxiliary basis name, to use density fitting
    :param threads: the number of OpenMP threads to use for each calculation
    :param max_memory: PySCF's memory limit in MB
    :param soscf_fallback: if the first-order SCF doesn't converge, continue with a second-order (Newton) SCF
    :param allow_unconverged: warn instead of raising when the SCF still hasn't converged
    :param reuse_guess: start each SCF from the previous geometry's density (faster along scans and
        optimizations; fresh guesses are used if that fails)
    :param analytic_hessian: treat analytic Hessians as available (HF/DFT only)
    :param charge_model: ``'mulliken'`` (default) or ``'meta_lowdin'`` population charges
    """
    distance_units = "BohrRadius"
    property_units = {
        'energy': 'Hartrees',
        'charges': 'ElementaryCharge',
        'dipole': ('ElementaryCharge', 'BohrRadius'),
        'polarizability': None,  # atomic units, Bohr^3
    }
    analytic_derivative_orders = {
        'energy': 1,
        'charges': 0,
        'dipole': 0,
        'polarizability': 0,
    }
    property_shapes = {
        'charges': lambda n: (n,),
        'dipole': (3,),
        'polarizability': (3, 3),
    }

    def __init__(self,
                 atoms, *,
                 level_of_theory,
                 basis_set=None,
                 basis=None,
                 caller=None,
                 charge=None,
                 multiplicity=None,
                 quiet=True,
                 molecule_options=None,
                 level_of_theory_options=None,
                 disp=None,
                 scf_options=None,
                 grid_level=None,
                 density_fit=None,
                 threads=None,
                 max_memory=None,
                 soscf_fallback=True,
                 allow_unconverged=False,
                 reuse_guess=False,
                 analytic_hessian=False,
                 charge_model='mulliken'):
        if basis_set is None:
            basis_set = basis
        elif basis is not None and basis != basis_set:
            raise ValueError(f"got two different basis sets: {basis_set!r} and {basis!r}")
        if basis_set is None:
            raise ValueError("PySCFBackend needs a `basis_set`")

        molecule_options = dict(molecule_options or {})
        if molecule_options.get('unit', 'Bohr').lower() not in ('bohr', 'b', 'au', 'a.u.'):
            raise ValueError("coordinates are supplied in Bohr; `molecule_options['unit']` can't be changed")
        molecule_options['unit'] = 'Bohr'
        # explicit `molecule_options` charge/spin take precedence, as before
        if 'charge' in molecule_options:
            charge = molecule_options.pop('charge')
        if 'spin' in molecule_options:
            multiplicity = molecule_options.pop('spin') + 1  # PySCF's spin is 2S
        if quiet:
            molecule_options['verbose'] = molecule_options.get('verbose', 0)
        if max_memory is not None:
            molecule_options['max_memory'] = max_memory
        self.molecule_options = molecule_options

        super().__init__(atoms, charge=charge, multiplicity=multiplicity, quiet=quiet)
        self.basis_set = basis_set
        self.level_of_theory = level_of_theory
        self.level_of_theory_options = dict(level_of_theory_options or {})
        if disp is not None:
            self.level_of_theory_options['disp'] = disp
        self.scf_options = dict(scf_options or {})
        self.grid_level = grid_level
        self.density_fit = density_fit
        self.threads = threads
        self.soscf_fallback = soscf_fallback
        self.allow_unconverged = allow_unconverged
        self.reuse_guess = reuse_guess
        self.charge_model = charge_model
        self._last_dm = None
        if analytic_hessian:
            self.analytic_derivative_orders = dict(self.analytic_derivative_orders, energy=2)
        self._pyscf_caller_dispatch = dev.uninitialized
        self.caller = self.resolve_caller(caller, level_of_theory)

    @classmethod
    def handle_specialization(cls, tag):
        if '/' not in tag:
            raise ValueError(f"expected a 'method/basis' tag for pyscf (e.g. 'b3lyp/cc-pvdz'), got {tag!r}")
        lot, basis = tag.split("/", 1)
        return {'level_of_theory': lot, 'basis_set': basis}

    @classmethod
    def discover_dependencies(cls, *, disp=None, level_of_theory=None, **opts):
        """
        PySCF, plus ``pyscf-dispersion`` if a dispersion correction is requested (via ``disp`` or a ``-d3bj``-style
        suffix on the level of theory).
        """
        from McUtils.ExternalPrograms import JobDependencies
        deps = JobDependencies()
        deps.add('pip', 'pyscf', import_name='pyscf', reason='PySCF itself')
        deps.add('conda', 'pyscf', import_name='pyscf', reason='PySCF itself')
        lot = level_of_theory if isinstance(level_of_theory, str) else ''
        if disp is not None or cls._split_dispersion(lot)[1] is not None:
            deps.add('pip', 'pyscf-dispersion', import_name='pyscf.dispersion', reason='dispersion corrections')
        return deps

    # region calculation setup

    _default_caller = 'dft'
    @property
    def caller_dispatch(self) -> dev.OptionsMethodDispatch:
        self._pyscf_caller_dispatch = dev.handle_uninitialized(
            self._pyscf_caller_dispatch,
            dev.OptionsMethodDispatch,
            args=(self.get_callers,),
            kwargs=dict(
                default_method=self._default_caller,
            )
        )
        return self._pyscf_caller_dispatch

    def get_callers(self):
        return {
            'mp2': self._call_mp2,
            'dft': self._call_dft,
            'hf': self._call_hf,
            'rhf': self._call_rhf,
            'uhf': self._call_uhf,
            'rohf': self._call_rohf,
        }

    def resolve_caller(self, caller, lot):
        if caller is None:
            if callable(lot):
                caller = lot
            else:
                caller, opts = self.caller_dispatch.resolve(lot)
                self.level_of_theory_options.update(opts)
        return caller

    def get_molecule(self, coords, charge=None, multiplicity=None):
        gto = import_or_explain('pyscf.gto', 'pyscf')
        charge, multiplicity = self.resolve_state(charge, multiplicity)
        opts = dict(self.molecule_options)
        if charge is not None:
            opts['charge'] = charge
        if multiplicity is not None:
            opts['spin'] = multiplicity - 1  # PySCF wants 2S, not 2S+1

        molecule = gto.Mole(
            basis=self.basis_set,
            **opts
        )
        molecule.atom = [
            (a, crd)
            for a,crd in zip(self.atoms, coords)
        ]
        molecule.build()

        return molecule

    def _setup_scf(self, mf, attrs=None):
        if self.density_fit:
            mf = mf.density_fit(auxbasis=None if self.density_fit is True else self.density_fit)
        for k, v in dict(self.scf_options, **(attrs or {})).items():
            setattr(mf, k, v)
        return mf

    def _converge(self, mf):
        """
        Run the SCF, falling back to a fresh guess and then a second-order solver if it doesn't converge, and
        refuse to hand back an unconverged result (unless `allow_unconverged`).
        """
        dm0 = self._last_dm if self.reuse_guess else None
        if dm0 is not None and np.shape(dm0)[-1] != mf.mol.nao:
            dm0 = None
        mf.kernel(dm0=dm0)
        if not mf.converged and dm0 is not None:
            mf.kernel()
        if not mf.converged and self.soscf_fallback:
            mf = mf.newton()
            mf.kernel()
        if not mf.converged:
            msg = f"PySCF SCF ({self.level_of_theory}/{self.basis_set}) did not converge"
            if self.allow_unconverged:
                warnings.warn(msg)
            else:
                raise RuntimeError(msg)
        if self.reuse_guess:
            self._last_dm = mf.make_rdm1()
        return mf

    def _is_restricted(self, molecule, restricted):
        return (molecule.spin == 0) if restricted is None else restricted

    _dispersion_suffix = re.compile(r"^(.+?)[-_]?(d3bj|d3zero|d3|d4)$", re.I)
    @classmethod
    def _split_dispersion(cls, lot):
        m = cls._dispersion_suffix.match(lot or '')
        if m is None:
            return lot, None
        disp = m.group(2).lower()
        return m.group(1), 'd3zero' if disp == 'd3' else disp

    def _call_dft(self, molecule, restricted=None, functional_options=None,  disp=None, **opts):
        dft = import_or_explain('pyscf.dft', 'pyscf')

        functional_options = dict(functional_options or {})
        xc, suffix_disp = self._split_dispersion(self.level_of_theory)
        functional_options['xc'] = xc
        d = functional_options.pop('disp', None)
        if disp is None:
            disp = d if d is not None else suffix_disp
        if self._is_restricted(molecule, restricted):
            mf = dft.RKS(molecule, **functional_options)
        else:
            mf = dft.UKS(molecule, **functional_options)
        if self.grid_level is not None:
            mf.grids.level = self.grid_level
        if disp is not None:
            mf.disp = disp
        return self._converge(self._setup_scf(mf, opts))

    def _call_hf(self, molecule, restricted=None, **opts):
        scf = import_or_explain('pyscf.scf', 'pyscf')
        if self._is_restricted(molecule, restricted):
            mf = scf.RHF(molecule)
        else:
            mf = scf.UHF(molecule)
        return self._converge(self._setup_scf(mf, opts))
    def _call_rhf(self, molecule, **opts):
        return self._call_hf(molecule, restricted=True, **opts)
    def _call_uhf(self, molecule, **opts):
        return self._call_hf(molecule, restricted=False, **opts)
    def _call_rohf(self, molecule, **opts):
        scf = import_or_explain('pyscf.scf', 'pyscf')
        return self._converge(self._setup_scf(scf.ROHF(molecule), opts))

    def _call_mp2(self, molecule, reference='hf', **etc):
        mp = import_or_explain('pyscf.mp', 'pyscf')
        hf_caller, opts = self.caller_dispatch.resolve(reference)
        mf = hf_caller(molecule, **opts)
        return mp.MP2(mf, **etc).run()

    def run_calculation(self, coords, charge=None, multiplicity=None, **opts):
        molecule = self.get_molecule(coords, charge=charge, multiplicity=multiplicity)
        if self.threads is not None:
            lib = import_or_explain('pyscf.lib', 'pyscf')
            with lib.with_omp_threads(self.threads):
                return self.caller(molecule, **(self.level_of_theory_options | opts))
        return self.caller(molecule, **(self.level_of_theory_options | opts))

    # endregion

    # region results

    def get_energy_from_calc(self, calc) -> float:
        if hasattr(calc, 'e_tot'):
            return float(calc.e_tot)
        elif hasattr(calc, 'tot_energy'):
            return float(calc.tot_energy())
        else:
            return float(calc.kernel())

    def get_gradient_from_calc(self, calc):
        grad = calc.nuc_grad_method()
        grad.verbose = 0
        return np.asanyarray(grad.kernel(), dtype=float)

    def get_hessian_from_calc(self, calc):
        if not hasattr(calc, 'Hessian'):
            raise NotImplementedError(f"PySCF has no analytic Hessian for {type(calc).__name__}")
        hess = calc.Hessian()
        hess.verbose = 0
        h = np.asanyarray(hess.kernel(), dtype=float)  # (natoms, natoms, 3, 3)
        n = h.shape[0]
        return h.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)

    @staticmethod
    def _scf_only(calc, prop):
        # the SCF density of a correlated calculation would silently be the wrong level of theory
        if not hasattr(calc, 'make_rdm1') or not hasattr(calc, 'mo_coeff') or not hasattr(calc, 'dip_moment'):
            raise NotImplementedError(f"{prop} are only available for SCF (HF/DFT) calculations, not {type(calc).__name__}")
        return calc

    def get_charges_from_calc(self, calc):
        mf = self._scf_only(calc, "charges")
        if self.charge_model == 'mulliken':
            pop, chg = mf.mulliken_pop(verbose=0)
        elif self.charge_model in ('meta_lowdin', 'meta-lowdin', 'mulliken_meta'):
            pop, chg = mf.mulliken_meta(verbose=0)
        else:
            raise ValueError(f"unknown charge model {self.charge_model!r} (use 'mulliken' or 'meta_lowdin')")
        return np.asanyarray(chg, dtype=float)

    def get_dipole_from_calc(self, calc):
        mf = self._scf_only(calc, "dipoles")
        return np.asanyarray(mf.dip_moment(unit='AU', verbose=0), dtype=float)

    def get_polarizability_from_calc(self, calc):
        mf = self._scf_only(calc, "polarizabilities")
        scf = import_or_explain('pyscf.scf', 'pyscf')
        dft = import_or_explain('pyscf.dft', 'pyscf')
        unrestricted = isinstance(mf, scf.uhf.UHF)
        if isinstance(mf, dft.rks.KohnShamDFT):
            kind = 'uks' if unrestricted else 'rks'
        else:
            kind = 'uhf' if unrestricted else 'rhf'
        polar = import_or_explain(f'pyscf.prop.polarizability.{kind}', 'pyscf-properties')
        p = polar.Polarizability(mf)
        p.verbose = 0
        return np.asanyarray(p.polarizability(), dtype=float)

    # endregion

    def evaluate_properties(self, coords, orders, charge, multiplicity, **opts):
        nstruct, n = coords.shape[:2]
        out = {p: [np.empty((nstruct,) + self.get_property_shape(p))] for p in orders}
        e_order = orders.get('energy', 0)
        if e_order > 0:
            out['energy'].append(np.empty((nstruct, 3 * n)))
        if e_order > 1:
            out['energy'].append(np.empty((nstruct, 3 * n, 3 * n)))
        getters = {
            'charges': self.get_charges_from_calc,
            'dipole': self.get_dipole_from_calc,
            'polarizability': self.get_polarizability_from_calc,
        }
        for i, coord in enumerate(coords):
            res = self.run_calculation(coord, charge=charge, multiplicity=multiplicity, **opts)
            for p, o in orders.items():
                if p == 'energy':
                    out[p][0][i] = self.get_energy_from_calc(res)
                    if o > 0:
                        out[p][1][i] = self.get_gradient_from_calc(res).reshape(-1)
                    if o > 1:
                        out[p][2][i] = self.get_hessian_from_calc(res)
                else:
                    out[p][0][i] = getters[p](res)
        return out
