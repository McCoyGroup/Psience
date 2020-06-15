from .AbstractDMC import AbstractDMC
import numpy as np

__all__ = ["ContinuousWeightingDMC"]

## needs updates with everything we've learned from the RynLib implementation
class ContinuousWeightingDMC(AbstractDMC):

    def __init__(self, *args, alpha = None, **kwargs):
        super().__init__(*args, **kwargs)

        self.parents = np.arange(len(self.walkers))
        self._parents = self.walkers.coords.copy()
        self._parent_weights = self.walkers.weights.copy()
        self.alpha = alpha

    def branch(self):
        weights = self.walkers.weights
        walkers = self.walkers.coords
        parents = self.parents
        threshold = 1.0 / len(self.walkers)
        eliminated_walkers = np.argwhere(weights < threshold)

        for dying in eliminated_walkers: # gotta do it iteratively to get the max_weight_walker right..
            cloning = np.argmax(weights)
            parents[dying] = parents[cloning]
            walkers[dying] = walkers[cloning]
            weights[dying] = weights[cloning] / 2.0
            weights[cloning] /= 2.0

    def _compute_vref(self, energies, weights):
        """Takes a single set of energies and weights and computes the average potential

        :param energies: single set of energies
        :type energies:
        :param weights: single set of weights
        :type weights:
        :return:
        :rtype: float
        """
        Vbar = np.average(energies, weights=weights, axis = 0)
        num_walkers = len(weights)
        correction=np.sum(weights-np.ones(num_walkers), axis = 0)/num_walkers
        vref = Vbar - (self.alpha * correction)
        return vref

    def update_weights(self, energies, weights):
        """Iteratively updates the weights over a set of vectors of energies

        :param energies:
        :type energies: np.ndarray
        :param weights:
        :type weights: np.ndarray
        :return:
        :rtype: np.ndarray
        """
        for e in energies: # this is basically a reduce call, but there's no real reason not to keep it like this
            Vref = self._compute_vref(e, weights)
            self.reference_potentials.append(Vref) # a constant time operation
            new_wts = np.exp(-1.0 * (e - Vref) * self.time_step)
            weights *= new_wts
        return weights

    def weight_descendants(self):
        do_it = self.equilibrated and (self.step_num - self._last_desc_weighting_step >= self.descendent_weighting_delay)
        if do_it:
            self._last_desc_weighting_step = self.step_num
            num_walkers = len(self.walkers.coords)
            weights = np.array( [ np.sum(self.walkers.weights[ self.parents == i ]) for i in range(num_walkers) ] )
            self.descendent_weights = (self._parents, weights, self._parent_weights)
            self._parents = self.walkers.coords.copy()
            self._parent_weights = self.walkers.weights.copy()
            self.parents = np.arange(num_walkers)