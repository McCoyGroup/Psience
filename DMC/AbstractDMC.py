from abc import *
import numpy as np
from .WalkerSet import WalkerSet
from .PotentialEvaluator import PotentialEvaluator

__all__ = ["AbstractDMC"]

## needs updates from the RynLib implementation

##################################################################################################################
##
##                                                  AbstractDMC
##
##
class AbstractDMC(metaclass=ABCMeta):
    """A completely abstract implementation of DMC"""

    def __init__(self,
                 name, description,
                 walker_set = None,
                 D = None, time_step = None,
                 steps_per_propagation = None,
                 num_time_steps = None,
                 equilibration = None,
                 potential = None,
                 descendent_weighting_delay = None,
                 log_file = None
                 ):
        from collections import deque

        self.name = name
        self.description = description

        if not isinstance(walker_set, WalkerSet):
            walker_set=WalkerSet(*walker_set)
        self.walkers = walker_set
        self.D = D
        self.time_step = time_step
        self.walkers.initialize(self.time_step, self.D)

        if not isinstance(potential, PotentialEvaluator):
            potential = PotentialEvaluator(potential)
        self.potential = potential

        self.log = log_file # we'll dump to binary here
        self.reference_potentials = deque()

        self.step_num = 0
        self.num_time_steps = num_time_steps
        self.steps_per_propagation = steps_per_propagation

        self.descendent_weighting_delay = descendent_weighting_delay
        self.descendent_weights = None

        self._equilibrated = False
        if isinstance(equilibration, int):
            equilibration = lambda s, e=equilibration: s.step_num > e #classic lambda parameter binding
        self.equilibration_check = equilibration

    @property
    def zpe(self):
        return self.get_zpe()
    def get_zpe(self, n = 30):
        import itertools
        return np.average(np.array(itertools.islice(self.reference_potentials, -n, 1)))

    @property
    def equilibrated(self):
        if not self._equilibrated:
            self._equilibrated = self.equilibration_check(self)
        return self._equilibrated

    def get_energy(self, coords=None, atoms=None):
        """Handles potential calls

        :return: list of potentials
        :rtype:
        """
        if coords is None:
            coords = self.walkers.coords
        coords = np.asarray(coords)
        if atoms is None:
            atoms = self.walkers.atoms

        if len(coords.shape) > 3:
            return np.array([ self.potential(atoms, crds) for crds in coords ])
        else:
            return self.potential(atoms, coords)

    @abstractmethod
    def branch(self):
        """Handles branching in DMC. Returns the new energy array after the branching occurs.

        :return:
        :rtype:
        """
        pass

    @abstractmethod
    def update_weights(self, energies, weights):
        """Handles weights in DMC"""
        pass

    @abstractmethod
    def weight_descendants(self):
        """Applies descendant weighting"""
        pass

    def propagate(self, nsteps = None):
        """Propagates the system forward n steps

        :param nsteps: number of steps to propagate for; None means automatic
        :type nsteps:
        :return:
        :rtype:
        """
        if nsteps is None:
            nsteps = self.steps_per_propagation

        coord_sets = self.walkers.displace(nsteps)
        energies = self.potential(self.walkers.atoms, coord_sets)
        self.step_num += nsteps

        weights = self.update_weights(energies, self.walkers.weights)
        self.walkers.weights = weights

        self.branch()
        self.weight_descendants()

    def snapshot(self, file=None):
        """Saves a snapshot of the simulation to file

        :param file:
        :type file:
        :return:
        :rtype:
        """
        import pickle

        if file is None:
            file = self.log
        if isinstance(file, str):
            with open(file, "w+") as dump:
                pickle.dump(self, dump)
        else:
            pickle.dump(self, file)

