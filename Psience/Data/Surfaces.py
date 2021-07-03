"""
Provides concrete tools for dealing with two of the most useful types of surfaces we have
"""

import numpy as np
from collections import namedtuple
from McUtils.GaussianInterface import GaussianLogReader, GaussianFChkReader
from McUtils.Zachary import Surface, MultiSurface, InterpolatedSurface, TaylorSeriesSurface

__all__=[
    "DipoleSurface",
    "PotentialSurface"
]

class DipoleSurface(MultiSurface):
    """
    Provides a unified interface to working with dipole surfaces.
    Currently basically no fancier than a regular surface (although with convenient loading functions), but dipole-specific
    stuff could come
    """
    def __init__(self, mu_x, mu_y, mu_z):
        """

        :param mu_x: X-component of dipole moment
        :type mu_x: Surface
        :param mu_y: Y-component of dipole moment
        :type mu_y: Surface
        :param mu_z: Z-component of dipole moment
        :type mu_z: Surface
        """
        if isinstance(mu_x.base, TaylorSeriesSurface):
            self.mode = "taylor"
        else:
            self.mode = "interp"
        super().__init__(
            mu_x,
            mu_y,
            mu_z
        )

    @staticmethod
    def get_log_values(log_file, keys=("StandardCartesianCoordinates", "DipoleMoments")):

        with GaussianLogReader(log_file) as parser:
            parse_data = parser.parse(keys)

        carts = parse_data[keys[0]][1]
        dipoles = np.array(parse_data[keys[1]])

        return namedtuple('dipole_log_values', ['cartesians', 'dipoles'])(carts, dipoles)
    @classmethod
    def from_log_file(cls, log_file, coord_transf, keys=("StandardCartesianCoordinates", "DipoleMoments"), tol = .001, **opts):
        """
        Loads dipoles from a Gaussian log file and builds a dipole surface by interpolating.
        Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
        to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
        Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions

        :param log_file: a Gaussian log file to pull from
        :type log_file: str
        :return:
        :rtype:
        """

        carts, dipoles = cls.get_log_values(log_file, keys=keys)

        scan_coords = coord_transf(carts)
        if len(dipoles) != len(scan_coords):
            raise ValueError(
                "mismatch between number of dipoles ({}) and number of coordinates ({})".format(
                    len(dipoles),
                    len(scan_coords)
                )
            )

        if scan_coords.ndim == 1:
            scan_sort = np.argsort(scan_coords)
        else:
            scan_sort = np.lexsort(tuple(reversed(tuple(scan_coords.T))))
        scan_coords = scan_coords[scan_sort]
        dipoles = dipoles[scan_sort]

        # this is particularly relevant for optimization scans...but we pull all the "unique" indices
        # then we pull the indices right _before_ each unique one since that's the final one in the block of "uniques"
        # finally we do a "roll" to make sure the order from the sort is preserved
        tol_coords = np.floor(scan_coords/tol)
        if tol_coords.ndim == 1:
            diffs = np.diff(tol_coords)
        else:
            diffs = np.sum(abs(np.diff(tol_coords, axis=0)), axis=1)
        inds = np.where(diffs != 0)[0]
        inds = np.concatenate((inds, [len(inds)]))
        scan_coords = scan_coords[inds]
        dipoles = dipoles[inds]

        dipoles = list(np.transpose(dipoles))

        return DipoleSurface(*(
            Surface(
                ((scan_coords, d), opts),
                base = InterpolatedSurface,
                dipole_component = "x" if i == 0 else "y" if i == 1 else "z"
            ) for i,d in enumerate(dipoles)
        ))

    @staticmethod
    def get_fchk_values(fchk_file):

        with GaussianFChkReader(fchk_file) as parser:
            parse_data = parser.parse(["Coordinates", "Dipole Moment", "Dipole Derivatives"])

        center = parse_data["Coordinates"]
        const_dipole = parse_data["Dipole Moment"]
        derivs = parse_data["Dipole Derivatives"]
        derivs = np.reshape(derivs, (int(len(derivs) / 3), 3))

        return namedtuple('dipole_fchk_values', ['center', 'const', 'derivs'])(center, const_dipole, derivs)
    @classmethod
    def from_fchk_file(cls, fchk_file, **opts):
        """
        Loads dipoles from a Gaussian formatted checkpoint file and builds a dipole surface via a linear approximation

        :param fchk_file: a Gaussian fchk file to pull from
        :type log_file: str
        :return:
        :rtype:
        """

        center, const_dipole, derivs = cls.get_fchk_values(fchk_file)
        derivs = list(np.transpose(derivs))

        opts['center'] = center.flatten()
        surfs = [None]*3
        for i, d in enumerate(zip(derivs, list(const_dipole))):
            d, r = d
            opts = opts.copy()
            opts["ref"] = r
            surfs[i] = Surface(
                ((d,), opts),
                base = TaylorSeriesSurface,
                dipole_component="x" if i == 0 else "y" if i == 1 else "z"
            )

        return cls(*surfs)

    def __call__(self, gridpoints, **opts):
        """
        Explicitly overrides the Surface-level evaluation because we know the Taylor surface needs us to flatten our gridpoints

        :param gridpoints:
        :type gridpoints:
        :param opts:
        :type opts:
        :return:
        :rtype:
        """

        gps = np.asarray(gridpoints)
        if self.mode == "taylor":
            if gps.ndim == 2:
                gps = gps.flatten()
            elif gps.ndim > 2:
                gps = np.reshape(gps, gps.shape[:-2] + (np.product(gps.shape[-2:]),))

        return super().__call__(gps, **opts)

class PotentialSurface(Surface):
    """
    A potential surface structure to go along with the DipoleSurface.
    Provides convenient access to dipole data + a unified interface to things like energy minimization
    """

    @staticmethod
    def get_log_values(log_file, keys=("StandardCartesianCoordinates", "ScanEnergies")):

        with GaussianLogReader(log_file) as parser:
            parse_data = parser.parse(keys)

        # Need to be smarter about this. At some point we might be able to infer what type of log file we have...
        coord_key = keys[0]
        coords = parse_data[coord_key][1]
        eng_key = keys[1]
        if eng_key == "ScanEnergies":
            energies = np.array(parse_data[eng_key].energies[:, -1])
        else:
            raise Exception("Haven't dealt with scan types beyond rigid ones...")

        return namedtuple('potential_log_values', ['coords', 'energies'])(coords, energies)
    @classmethod
    def from_log_file(cls, log_file, coord_transf, keys=("StandardCartesianCoordinates", "ScanEnergies"), tol = .001, **opts):
        """
        Loads dipoles from a Gaussian log file and builds a potential surface by interpolating.
        Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
        to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
        Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions.

        :param log_file: a Gaussian log file to pull from
        :type log_file: str
        :return:
        :rtype:
        """

        dat = cls.get_log_values(log_file, keys=keys)
        carts = dat.coords
        pots = dat.energies

        # raise Exception(carts, pots)
        scan_coords = coord_transf(carts)
        if len(pots) != len(scan_coords):
            raise ValueError(
                "mismatch between number of potential values ({}) and number of coordinates ({})".format(
                    len(pots),
                    len(scan_coords)
                )
            )

        if scan_coords.ndim == 1:
            scan_sort = np.argsort(scan_coords)
        else:
            scan_sort = np.lexsort(tuple(reversed(tuple(scan_coords.T))))
        scan_coords = scan_coords[scan_sort]
        pots = pots[scan_sort]

        # this is particularly relevant for optimization scans...but we pull all the "unique" indices
        # then we pull the indices right _before_ each unique one since that's the final one in the block of "uniques"
        # finally we do a "roll" to make sure the order from the sort is preserved
        tol_coords = np.floor(scan_coords/tol)
        if tol_coords.ndim == 1:
            diffs = np.diff(tol_coords)
        else:
            diffs = np.sum(abs(np.diff(tol_coords, axis=0)), axis=1)
        inds = np.where(diffs != 0)[0]
        inds = np.concatenate((inds, [len(inds)]))
        scan_coords = scan_coords[inds].squeeze()
        pots = pots[inds]

        return cls(
                ((scan_coords, pots), opts),
                base=InterpolatedSurface
        )

    @staticmethod
    def get_fchk_values(fchk_file):
        # TODO: I know I probably didn't do this right but I'm just getting a thing out for now
        with GaussianFChkReader(fchk_file) as parser:
            parse_data = parser.parse(["Coordinates", "Total Energy", "Gradient", "ForceConstants", "ForceDerivatives"])

        center = parse_data["Coordinates"]
        eng = parse_data["Total Energy"]
        derivs = [parse_data['Gradient'], parse_data["ForceConstants"], parse_data["ForceDerivatives"]]

        return namedtuple('potential_fchk_values', ['center', 'energy', 'derivs'])(
            center, eng, derivs
        )
    @classmethod
    def from_fchk_file(cls, fchk_file, **opts):
        """
        Loads potential from a Gaussian formatted checkpoint file and builds a potential surface via a quartic approximation

        :param fchk_file: a Gaussian fchk file to pull from
        :type log_file: str
        :return:
        :rtype:
        """

        center, energy, derivs = cls.get_fchk_values(fchk_file)

        return cls((derivs, dict(ref=energy, center=center.flatten())), base=TaylorSeriesSurface, **opts)