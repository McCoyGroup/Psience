"""
Atom-centered densities as Gaussian-smeared multipole expansions
"""

import numpy as np

from McUtils.Data import AtomData
from McUtils.Zachary import MultipoleExpansion

__all__ = [
    "AtomCenteredDensity",
]


class AtomCenteredDensity:
    """
    A density written as a sum of atom-centered Gaussian-smeared multipoles,

        rho(r) = sum_i sum_lm p_ilm phi_lm(r - R_i; sigma_i),

    in the GTO basis of `graph_longrange`/MACE-POLAR, in which `p_ilm` are the exact multipole
    moments of each site (so the `l = 0` coefficients are condensed atomic populations).
    Moments use the real Racah-normalized layout `m = -l..l` per `l` (`[q, mu_y, mu_z, mu_x, ...]`).

    Everything is stored in atomic units: centers and widths in Bohr, moments in `e Bohr^l`
    (or the corresponding units of whatever is being represented, e.g. unpaired electrons for a
    spin density). Atoms are the last batch axis.
    """

    def __init__(self, multipoles, centers, smearing=None, atoms=None, label=None):
        """
        :param multipoles: real moments per site, `(..., natoms, (L+1)^2)` (or `(..., natoms)` for populations)
        :type multipoles: np.ndarray
        :param centers: site positions in Bohr, `(..., natoms, 3)`
        :type centers: np.ndarray
        :param smearing: Gaussian width(s) in Bohr (`None` for point sites: populations only)
        :type smearing: float | np.ndarray | None
        :param atoms: atom symbols
        :type atoms: list[str] | None
        :param label: a name for the density
        :type label: str | None
        """
        multipoles = np.asanyarray(multipoles, dtype=float)
        centers = np.asanyarray(centers, dtype=float)
        if multipoles.ndim == centers.ndim - 1:
            multipoles = multipoles[..., np.newaxis]
        n = multipoles.shape[-1]
        L = int(round(np.sqrt(n))) - 1
        if (L + 1) ** 2 != n:
            raise ValueError(f"trailing dimension {n} is not a real (L+1)^2 multipole size")
        self.multipoles = multipoles
        self.centers = centers
        self.smearing = None if smearing is None else np.asanyarray(smearing, dtype=float)
        self.atoms = None if atoms is None else list(atoms)
        self.label = label
        self._expansion = None

    # region construction / arithmetic

    @classmethod
    def from_populations(cls, populations, centers, smearing=None, atoms=None, label=None):
        """
        :param populations: condensed (monopole) values per site, `(..., natoms)`
        :type populations: np.ndarray
        :rtype: AtomCenteredDensity
        """
        return cls(np.asanyarray(populations, dtype=float)[..., np.newaxis], centers,
                   smearing=smearing, atoms=atoms, label=label)

    def modify(self, multipoles=None, label=None, smearing=None):
        return type(self)(
            self.multipoles if multipoles is None else multipoles,
            self.centers,
            smearing=self.smearing if smearing is None else smearing,
            atoms=self.atoms,
            label=self.label if label is None else label
        )

    @property
    def order(self):
        return int(round(np.sqrt(self.multipoles.shape[-1]))) - 1

    def _aligned(self, other):
        if not isinstance(other, AtomCenteredDensity):
            return NotImplemented
        if self.centers.shape != other.centers.shape or not np.allclose(self.centers, other.centers):
            raise ValueError("densities can only be combined on the same atomic centers")
        if (self.smearing is None) != (other.smearing is None) or (
                self.smearing is not None and not np.allclose(self.smearing, other.smearing)):
            raise ValueError("densities can only be combined with the same smearing")
        L = max(self.order, other.order)
        def pad(m, L):
            k = (L + 1) ** 2
            if m.shape[-1] == k:
                return m
            out = np.zeros(m.shape[:-1] + (k,))
            out[..., :m.shape[-1]] = m
            return out
        return pad(self.multipoles, L), pad(other.multipoles, L)

    def __add__(self, other):
        res = self._aligned(other)
        if res is NotImplemented:
            return res
        a, b = res
        return self.modify(multipoles=a + b, label=None)

    def __sub__(self, other):
        res = self._aligned(other)
        if res is NotImplemented:
            return res
        a, b = res
        return self.modify(multipoles=a - b, label=None)

    def __neg__(self):
        return self.modify(multipoles=-self.multipoles)

    def __mul__(self, scalar):
        return self.modify(multipoles=self.multipoles * np.asanyarray(scalar, dtype=float)[..., np.newaxis, np.newaxis])

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self * (1 / np.asanyarray(scalar, dtype=float))

    # endregion

    # region condensed quantities

    @property
    def populations(self):
        """the condensed (monopole) values per site, `(..., natoms)`"""
        return self.multipoles[..., 0]

    @property
    def total(self):
        """the integral of the density"""
        return np.sum(self.populations, axis=-1)

    @property
    def dipoles(self):
        """site dipoles (Cartesian, `e Bohr`), `(..., natoms, 3)`"""
        if self.order < 1:
            return np.zeros(self.multipoles.shape[:-1] + (3,))
        return self.multipoles[..., [3, 1, 2]]

    def _mask(self, heavy_atoms_only):
        if not heavy_atoms_only:
            return slice(None)
        if self.atoms is None:
            raise ValueError("`heavy_atoms_only` needs `atoms`")
        return np.array([AtomData[a, "Number"] > 1 for a in self.atoms])

    def participation_ratio(self, heavy_atoms_only=False):
        """
        The participation ratio `(sum_i |p_i|)^2 / sum_i p_i^2` of the condensed values: the
        effective number of sites the density is spread over (1 for a single site, N for
        an even spread over N sites). Negative populations count by magnitude.

        :param heavy_atoms_only: restrict to non-hydrogen sites
        :type heavy_atoms_only: bool
        :rtype: np.ndarray
        """
        p = self.populations[..., self._mask(heavy_atoms_only)]
        return np.sum(np.abs(p), axis=-1) ** 2 / np.sum(p ** 2, axis=-1)

    def localization(self, heavy_atoms_only=False):
        """the inverse participation ratio, between `1/N` (delocalized) and 1 (one site)"""
        return 1 / self.participation_ratio(heavy_atoms_only=heavy_atoms_only)

    def condensed(self, atoms=None):
        """the condensed values as a list of `(index, atom, value)`"""
        atoms = self.atoms if atoms is None else atoms
        p = self.populations
        if p.ndim != 1:
            raise ValueError("`condensed` is defined for a single structure")
        return [(i, None if atoms is None else atoms[i], float(v)) for i, v in enumerate(p)]

    # endregion

    # region continuous representation

    @property
    def expansion(self):
        """the underlying distributed `MultipoleExpansion` (atomic units)"""
        if self._expansion is None:
            self._expansion = MultipoleExpansion.from_real(self.multipoles, center=self.centers, smearing=self.smearing)
        return self._expansion

    def density(self, points):
        """
        The density at `points` (Bohr), per Bohr^3. Requires smeared sites.

        :param points: `(N, 3)` or `(..., N, 3)`
        :type points: np.ndarray
        :rtype: np.ndarray
        """
        if self.smearing is None:
            raise ValueError("point-site densities only have condensed values; give a `smearing` (`modify(smearing=...)`)")
        return self.expansion.density(points, sum_axis=-1)

    def potential(self, points, deriv_order=0):
        """
        The Coulomb potential of the density at `points` (Bohr), in atomic units, with optional
        Cartesian derivatives (as `MultipoleExpansion.evaluate`).
        """
        return self.expansion.evaluate(points, deriv_order=deriv_order, sum_axis=-1)

    # endregion

    def __repr__(self):
        return "{}({}natoms={}, order={}, total={})".format(
            type(self).__name__,
            "" if self.label is None else f"{self.label!r}, ",
            self.centers.shape[-2], self.order, np.round(self.total, 6)
        )
