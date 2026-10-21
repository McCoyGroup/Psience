"""
Models that predict atomic charge (and spin) densities for a given total charge and multiplicity
"""

import abc

import numpy as np

from McUtils.Data import AtomData, UnitsData

from .Densities import AtomCenteredDensity
from .Backends import ElectronicStructureBackend
from .MLIPBackends import MACEPolarBackend

__all__ = [
    "ElectronicState",
    "ElectronicStateModel",
    "MACEPolarStateModel",
    "BackendChargeStateModel",
    "ChargeEvaluatorStateModel",
]


class ElectronicState:
    """
    The densities a model predicts for one charge/spin state at one geometry.

    `charge_density` is the net charge density (nuclei included, so its populations are atomic
    partial charges); `spin_density` is majority-minus-minority spin (integrating to
    `multiplicity - 1`), or `None` if the model has no spin resolution.
    """

    def __init__(self, atoms, coords, charge, multiplicity, charge_density, spin_density=None, energy=None, model=None):
        """
        :param coords: nuclear coordinates in Bohr
        :type coords: np.ndarray
        :type charge_density: AtomCenteredDensity
        :type spin_density: AtomCenteredDensity | None
        :param energy: total energy in Hartree, if the model provides it
        :type energy: float | np.ndarray | None
        """
        self.atoms = list(atoms)
        self.coords = np.asanyarray(coords)
        self.charge = charge
        self.multiplicity = multiplicity
        self.charge_density = charge_density
        self.spin_density = spin_density
        self.energy = energy
        self.model = model

    @property
    def partial_charges(self):
        return self.charge_density.populations

    @property
    def electron_density(self):
        """
        The electron-number density implied by the charge density: populations `Z_i - q_i`,
        higher moments `-p_ilm` (the nuclei carry no higher moments)
        """
        Z = np.array([AtomData[a, "Number"] for a in self.atoms], dtype=float)
        m = -self.charge_density.multipoles.copy()
        m[..., 0] += Z
        return self.charge_density.modify(multipoles=m, label="electron density")

    def __repr__(self):
        return "{}(charge={}, multiplicity={}, natoms={}{})".format(
            type(self).__name__, self.charge, self.multiplicity, len(self.atoms),
            "" if self.spin_density is None else ", spin-resolved"
        )


class ElectronicStateModel(metaclass=abc.ABCMeta):
    """
    Predicts `ElectronicState`s. Concrete models take coordinates in their `distance_units`;
    `evaluate_state_bohr` takes Bohr.
    """
    distance_units = 'Angstroms'
    supports_spin = False

    def __init__(self, atoms):
        self.atoms = list(atoms)

    @abc.abstractmethod
    def evaluate_state(self, coords, charge=0, multiplicity=1):
        """
        :param coords: nuclear coordinates, `(natoms, 3)` in `distance_units`
        :type coords: np.ndarray
        :rtype: ElectronicState
        """
        ...

    def evaluate_state_bohr(self, coords, charge=0, multiplicity=1):
        return self.evaluate_state(
            np.asanyarray(coords) * UnitsData.convert("BohrRadius", self.distance_units),
            charge=charge, multiplicity=multiplicity
        )

    def evaluate_states(self, coords, states):
        """
        :param states: `(charge, multiplicity)` pairs
        :type states: Iterable[tuple[int, int]]
        :rtype: list[ElectronicState]
        """
        return [self.evaluate_state(coords, charge=c, multiplicity=m) for c, m in states]

    @classmethod
    def get_models(cls):
        return {
            'mace-polar': MACEPolarStateModel,
            'mace_polar': MACEPolarStateModel,
            'charges': BackendChargeStateModel,
            'charge-evaluator': ChargeEvaluatorStateModel,
            'charge_evaluator': ChargeEvaluatorStateModel,
        }

    @classmethod
    def from_molecule(cls, mol, model='mace-polar', **opts):
        """
        Resolve a model for a `Molecule`.

        :param model: a model instance, `'mace-polar[:size]'`, `'charges:backend[:tag]'` (any
            `ElectronicStructureBackend` that provides charges, e.g. `'charges:aimnet2'`,
            `'charges:xtb'`, `'charges:pyscf:b3lyp/def2-svp'`), or `'charge-evaluator:name'` (a
            charge evaluator from `mol.get_charge_evaluator`)
        :rtype: ElectronicStateModel
        """
        if isinstance(model, ElectronicStateModel):
            return model
        tag = None
        if ':' in model:
            model, tag = model.split(':', 1)
        model_type = cls.get_models()[model.lower()]
        return model_type.from_molecule(mol, tag=tag, **opts)


def _molecule_state_opts(mol, opts):
    opts = dict(opts)
    if mol.charge is not None:
        opts.setdefault('charge', mol.charge)
    if mol.spin is not None:
        opts.setdefault('multiplicity', mol.spin)
    return opts


class MACEPolarStateModel(ElectronicStateModel):
    """
    MACE-POLAR atomic multipoles (`density_coefficients`) and spin-channel multipoles
    (`spin_charge_density`, spin density = up channel minus down channel, which integrates to
    `multiplicity - 1`), Gaussian-smeared with the model's width.

    Benchmark note (MACE-POLAR-1-S, 2026-10-20): the spin channels carry per-atom "spin" of order
    +-1 even for closed-shell singlets (it integrates to zero) and miss the textbook spin patterns
    of allyl, benzyl and phenoxyl radicals, so they should not be read as a physical spin density
    without validation. Total-charge differences (Fukui functions) behave sensibly.
    """
    supports_spin = True

    def __init__(self, backend=None, atoms=None, **backend_opts):
        """
        :param backend: a `MACEPolarBackend` (built from `atoms` and `backend_opts` if not given)
        :type backend: MACEPolarBackend | None
        """
        if backend is None:
            if atoms is None:
                raise ValueError("need `atoms` or a `backend`")
            backend = MACEPolarBackend(atoms, **backend_opts)
        super().__init__(backend.atoms if atoms is None else atoms)
        self.backend = backend
        self.distance_units = backend.distance_units

    @classmethod
    def from_molecule(cls, mol, tag=None, **opts):
        if tag is not None:
            opts = dict(opts, **MACEPolarBackend.handle_specialization(tag))
        return cls(MACEPolarBackend(mol.atoms, **_molecule_state_opts(mol, opts)))

    @property
    def smearing(self):
        s = self.backend.smearing
        if s is None:
            raise ValueError("couldn't determine the MACE-POLAR smearing width")
        return s

    def evaluate_state(self, coords, charge=0, multiplicity=1):
        coords = np.asanyarray(coords, dtype=float)
        res = self.backend.evaluate_polar_properties(
            coords,
            properties=['energy', 'density_coefficients', 'spin_charge_density'],
            charge=charge, multiplicity=multiplicity
        )
        conv = UnitsData.convert(self.distance_units, "BohrRadius")
        dens = np.asanyarray(res['density_coefficients'])
        L = int(round(np.sqrt(dens.shape[-1]))) - 1
        scale = conv ** np.concatenate([np.full(2 * l + 1, l) for l in range(L + 1)])
        centers = coords * conv
        sigma = self.smearing * conv
        scd = np.asanyarray(res['spin_charge_density'])  # (..., natoms, 2, (L+1)^2)
        charge_density = AtomCenteredDensity(dens * scale, centers, smearing=sigma, atoms=self.atoms, label="charge density")
        spin_density = AtomCenteredDensity((scd[..., 0, :] - scd[..., 1, :]) * scale, centers, smearing=sigma,
                                           atoms=self.atoms, label="spin density")
        energy = np.asanyarray(res['energy']) * UnitsData.convert(self.backend.get_property_units('energy'), "Hartrees")
        return ElectronicState(self.atoms, centers, charge, multiplicity, charge_density,
                               spin_density=spin_density, energy=energy, model=self)


class BackendChargeStateModel(ElectronicStateModel):
    """
    Partial charges (and energies) from any `ElectronicStructureBackend` that provides
    `'charges'` (AIMNet2, xTB, PySCF, ASE calculators that report charges), evaluated per state
    with the backend's per-call charge/multiplicity overrides. No spin resolution; sites are
    point charges unless a `smearing` (Bohr) is supplied for visualization.
    """

    def __init__(self, backend, smearing=None, with_energy=True):
        """
        :type backend: ElectronicStructureBackend
        :param with_energy: also evaluate the state energy
        """
        super().__init__(backend.atoms)
        self.backend = backend
        self.smearing = smearing
        self.with_energy = with_energy and 'energy' in backend.supported_properties
        self.distance_units = backend.distance_units

    @classmethod
    def from_molecule(cls, mol, tag=None, backend=None, smearing=None, **opts):
        if backend is None:
            backend = 'aimnet2' if tag is None else tag
        backend = ElectronicStructureBackend.from_spec(backend, mol.atoms, **_molecule_state_opts(mol, opts))
        return cls(backend, smearing=smearing)

    def evaluate_state(self, coords, charge=0, multiplicity=1):
        coords = np.asanyarray(coords, dtype=float)
        props = ('charges', 'energy') if self.with_energy else ('charges',)
        res = self.backend.evaluate(coords, order=0, properties=props, charge=charge, multiplicity=multiplicity)
        q = res['charges'][0] * UnitsData.convert(self.backend.get_property_units('charges'), "ElementaryCharge")
        energy = None
        if self.with_energy:
            energy = res['energy'][0] * UnitsData.convert(self.backend.get_property_units('energy'), "Hartrees")
        centers = coords * UnitsData.convert(self.distance_units, "BohrRadius")
        charge_density = AtomCenteredDensity.from_populations(q, centers, smearing=self.smearing,
                                                              atoms=self.atoms, label="charge density")
        return ElectronicState(self.atoms, centers, charge, multiplicity, charge_density, energy=energy, model=self)


class ChargeEvaluatorStateModel(ElectronicStateModel):
    """
    Partial charges from a factory of charge evaluators, one per `(charge, multiplicity)`: any
    object with `evaluate(coords, order=0) -> [charges, ...]`, `distance_units`, and
    `target_property_units` (e.g. a `Psience.Molecools` `ChargeEvaluator`). No spin resolution;
    sites are point charges unless a `smearing` (Bohr) is supplied for visualization.
    """

    def __init__(self, evaluator_factory, atoms, smearing=None, distance_units=None):
        """
        :param evaluator_factory: `(charge, multiplicity) -> ChargeEvaluator`
        :type evaluator_factory: callable
        """
        super().__init__(atoms)
        self.evaluator_factory = evaluator_factory
        self.smearing = smearing
        self._evaluators = {}
        if distance_units is not None:
            self.distance_units = distance_units

    @classmethod
    def from_molecule(cls, mol, tag=None, evaluator=None, smearing=None, **opts):
        if evaluator is None:
            evaluator = 'aimnet2' if tag is None else tag
        factory = lambda charge, multiplicity: mol.get_charge_evaluator(
            evaluator, charge=charge, multiplicity=multiplicity, **opts
        )
        model = cls(factory, mol.atoms, smearing=smearing)
        model.distance_units = factory(mol.charge or 0, mol.spin or 1).distance_units
        return model

    def get_evaluator(self, charge, multiplicity):
        key = (charge, multiplicity)
        if key not in self._evaluators:
            self._evaluators[key] = self.evaluator_factory(charge, multiplicity)
        return self._evaluators[key]

    def evaluate_state(self, coords, charge=0, multiplicity=1):
        coords = np.asanyarray(coords, dtype=float)
        ev = self.get_evaluator(charge, multiplicity)
        q = ev.evaluate(coords, order=0)[0]
        q = q * UnitsData.convert(ev.target_property_units, "ElementaryCharge")
        centers = coords * UnitsData.convert(self.distance_units, "BohrRadius")
        charge_density = AtomCenteredDensity.from_populations(q, centers, smearing=self.smearing,
                                                              atoms=self.atoms, label="charge density")
        return ElectronicState(self.atoms, centers, charge, multiplicity, charge_density, model=self)
