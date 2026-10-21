"""
Spin densities and Fukui functions from `ElectronicStateModel`s
"""

import numpy as np

from McUtils.Data import UnitsData

from .StateModels import ElectronicStateModel

__all__ = [
    "SpinDensity",
    "FukuiFunctions",
]


def _molecule_coords(mol, model):
    return mol.coords * UnitsData.convert("BohrRadius", model.distance_units)


def _default_ion_multiplicity(multiplicity):
    # one electron more or fewer: the low-spin neighbor (singlet -> doublet, doublet -> singlet, ...)
    return multiplicity + 1 if multiplicity == 1 else multiplicity - 1


class SpinDensity:
    """
    The spin density (majority minus minority) of one state, optionally relative to a reference
    state at the same geometry.
    """

    def __init__(self, state, reference=None):
        """
        :type state: ElectronicState
        :param reference: a state whose spin density is subtracted (e.g. the closed-shell neutral),
            a heuristic for models whose spin channels carry a closed-shell partition; the
            difference still integrates to `multiplicity - 1` if the reference is a singlet
        :type reference: ElectronicState | None
        """
        if state.spin_density is None:
            raise ValueError(f"{type(state.model).__name__} does not resolve spin")
        self.state = state
        self.reference = reference

    @classmethod
    def compute(cls, model, coords, charge=0, multiplicity=2, reference=None):
        """
        :param model: the state model
        :type model: ElectronicStateModel
        :param coords: coordinates in the model's distance units
        :param reference: `(charge, multiplicity)` of a reference state to subtract, or `None`
        :type reference: tuple[int, int] | None
        :rtype: SpinDensity
        """
        state = model.evaluate_state(coords, charge=charge, multiplicity=multiplicity)
        ref = None if reference is None else model.evaluate_state(coords, charge=reference[0], multiplicity=reference[1])
        return cls(state, reference=ref)

    @classmethod
    def from_molecule(cls, mol, model='mace-polar', charge=None, multiplicity=None, reference=None, **opts):
        """
        :param mol: the molecule (its charge/spin are the defaults)
        :type mol: Molecule
        """
        model = ElectronicStateModel.from_molecule(mol, model=model, **opts)
        if charge is None:
            charge = 0 if mol.charge is None else mol.charge
        if multiplicity is None:
            multiplicity = 1 if mol.spin is None else mol.spin
        return cls.compute(model, _molecule_coords(mol, model), charge=charge, multiplicity=multiplicity,
                           reference=reference)

    @property
    def density(self):
        """the spin density as an `AtomCenteredDensity`"""
        s = self.state.spin_density
        if self.reference is not None:
            s = s - self.reference.spin_density
        return s.modify(label="spin density" if self.reference is None else "spin density difference")

    @property
    def populations(self):
        return self.density.populations

    def participation_ratio(self, heavy_atoms_only=False):
        return self.density.participation_ratio(heavy_atoms_only=heavy_atoms_only)

    def summary(self, top=8):
        d = self.density
        p = d.populations
        order = np.argsort(-np.abs(p))[:top]
        return {
            "charge": self.state.charge, "multiplicity": self.state.multiplicity,
            "reference": None if self.reference is None else (self.reference.charge, self.reference.multiplicity),
            "total": float(d.total),
            "participation_ratio": float(d.participation_ratio()),
            "participation_ratio_heavy_atoms": float(d.participation_ratio(heavy_atoms_only=True)),
            "largest": [(int(i), self.state.atoms[i], float(p[i])) for i in order],
        }


class FukuiFunctions:
    """
    Finite-difference (vertical) Fukui functions in the electron-density convention,

        f+ = n(N+1) - n(N)        (electron addition, LUMO-like)
        f- = n(N) - n(N-1)        (electron removal, HOMO-like)
        f0 = (f+ + f-) / 2,  dual descriptor = f+ - f-

    from the predicted charge densities of the reference state and its anion/cation at the same
    geometry. Each integrates to 1 (the dual descriptor to 0) and is an `AtomCenteredDensity`
    (condensed Fukui indices are its populations; the continuous form uses the model's smearing).
    """

    def __init__(self, reference, cation=None, anion=None):
        """
        :type reference: ElectronicState
        :type cation: ElectronicState | None
        :type anion: ElectronicState | None
        """
        self.reference = reference
        self.cation = cation
        self.anion = anion

    @classmethod
    def compute(cls, model, coords, charge=0, multiplicity=1,
                cation_multiplicity=None, anion_multiplicity=None, which=('+', '-')):
        """
        :param model: the state model
        :type model: ElectronicStateModel
        :param coords: coordinates in the model's distance units
        :param which: which ions to evaluate (`'+'` needs the anion, `'-'` the cation)
        :rtype: FukuiFunctions
        """
        ref = model.evaluate_state(coords, charge=charge, multiplicity=multiplicity)
        cation = anion = None
        if '-' in which:
            m = _default_ion_multiplicity(multiplicity) if cation_multiplicity is None else cation_multiplicity
            cation = model.evaluate_state(coords, charge=charge + 1, multiplicity=m)
        if '+' in which:
            m = _default_ion_multiplicity(multiplicity) if anion_multiplicity is None else anion_multiplicity
            anion = model.evaluate_state(coords, charge=charge - 1, multiplicity=m)
        return cls(ref, cation=cation, anion=anion)

    @classmethod
    def from_molecule(cls, mol, model='mace-polar', charge=None, multiplicity=None, **opts):
        """
        :param mol: the molecule (its charge/spin are the defaults)
        :type mol: Molecule
        """
        compute_opts = {k: opts.pop(k) for k in ('cation_multiplicity', 'anion_multiplicity', 'which') if k in opts}
        model = ElectronicStateModel.from_molecule(mol, model=model, **opts)
        if charge is None:
            charge = 0 if mol.charge is None else mol.charge
        if multiplicity is None:
            multiplicity = 1 if mol.spin is None else mol.spin
        return cls.compute(model, _molecule_coords(mol, model), charge=charge, multiplicity=multiplicity,
                           **compute_opts)

    # n = Z - (charge density), so electron-density differences are reversed charge-density differences
    @property
    def f_plus(self):
        if self.anion is None:
            raise ValueError("f+ needs the anion")
        return (self.reference.charge_density - self.anion.charge_density).modify(label="f+")

    @property
    def f_minus(self):
        if self.cation is None:
            raise ValueError("f- needs the cation")
        return (self.cation.charge_density - self.reference.charge_density).modify(label="f-")

    @property
    def f_zero(self):
        return ((self.f_plus + self.f_minus) / 2).modify(label="f0")

    @property
    def dual_descriptor(self):
        return (self.f_plus - self.f_minus).modify(label="dual descriptor")

    def summary(self, top=8):
        out = {"charge": self.reference.charge, "multiplicity": self.reference.multiplicity}
        for name in ('f_plus', 'f_minus'):
            try:
                f = getattr(self, name)
            except ValueError:
                continue
            p = f.populations
            order = np.argsort(-np.abs(p))[:top]
            out[name] = {
                "total": float(f.total),
                "participation_ratio": float(f.participation_ratio()),
                "participation_ratio_heavy_atoms": float(f.participation_ratio(heavy_atoms_only=True)),
                "largest": [(int(i), self.reference.atoms[i], float(p[i])) for i in order],
            }
        return out
