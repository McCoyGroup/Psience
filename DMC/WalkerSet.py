from McUtils.Data import UnitsData, AtomData
from McUtils.Coordinerds.CoordinateSystems import CoordinateSet
import numpy as np

__all__ = ["WalkerSet"]

##################################################################################################################
##
##                                                  WalkerSet
##
##
class WalkerSet:
    def __init__(self,
                 atoms,
                 coords,
                 *miscargs,
                 initial_walkers=5000,
                 masses = None,
                 sigmas = None,
                 weights = None,
                 _initialize = True
                 ):

        if _initialize:
            if len(miscargs) > 0:
                raise ValueError("{}: only the `atoms` and `coords` arguments should be passed as non-keyword args".format(
                    type(self).__name__
                ))

            atoms = [ AtomData[atom, "Symbol"] for atom in atoms ] # take whatever atom spec and symbolize it
            if masses is None:
                masses = np.array([ UnitsData.convert(AtomData[atom, "Mass"], "AtomicUnitOfMass") for atom in atoms ])

            if sigmas is None:
                sigmas = np.sqrt(masses)

            coords = CoordinateSet(coords) # by default CartesianCoordinates3D

            if isinstance(initial_walkers, int):
                #out of laziness we'll just duplicate our original a bunch
                initial_walkers =  coords*initial_walkers
            else:
                initial_walkers = CoordinateSet(initial_walkers)

        # I could add some fancy property tricks to prevent this stuff from being mangled, but it's not worth it...
        self.atoms = atoms
        self.masses = masses
        self.sigmas = sigmas
        self.coords = initial_walkers
        self.weights = weights

    def initialize(self, time_step, D):
        """Sets up necessary parameters for use in calculating displacements and stuff

        :param deltaT:
        :type deltaT:
        :param D:
        :type D:
        :return:
        :rtype:
        """
        self.time_step = time_step
        self.sigmas = np.sqrt((2 * D * time_step) / self.masses)

    def get_displacements(self, n=1, sigmas = None):
        """Computes n random Gaussian displacements from the sigmas and the timestep

        :param n:
        :type n:
        :param sigmas:
        :type sigmas:
        :return:
        :rtype:
        """
        if sigmas is None:
            sigmas = self.sigmas
        shape = (n, ) + self.coords.shape[:-2] + self.coords.shape[-1:]
        disps = np.array([
            np.random.normal(0.0, sig, size = shape) for sig in sigmas
        ])
        # transpose seems to be somewhat broken (?)
        # disp_roll = np.arange(len(disps.shape))
        # disp_roll = np.concatenate((np.roll(disp_roll[:-1], 1), disp_roll[-1:]))
        # print(disp_roll)
        # disps = disps.transpose(disp_roll)

        for i in range(len(shape) - 1):
            disps = disps.swapaxes(i, i + 1)
        return disps

    def get_displaced_coords(self, n=1):
        """Computes n new displaced coordinates based on the original ones from the sigmas

        :param n:
        :type n:
        :return:
        :rtype:
        """
        accum_disp = np.cumsum(self.get_displacements(n), axis=1)
        return np.broadcast_to(self.coords, (n,) + self.coords.shape) + accum_disp

    def displace(self, n=1):
        """Generates multiple displacements and binds the last one

        :param n:
        :type n:
        :return:
        :rtype:
        """
        coords = self.get_displaced_coords(n)
        self.coords = coords[-1]
        return coords