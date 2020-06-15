from Psience import Wavefunction
from McUtils.Plots import Plot
import numpy as np

__all__ = ["DMCWavefunction"]

class DMCWavefunction(Wavefunction):
    @property
    def coords(self):
        return self.data[0]
    @property
    def weights(self):
        return self.data[1]
    @property
    def descendant_weights(self):
        return self.data[2]

    def plot(self, figure = None, probability_density = False, bins = 20, **opts):
        #we'll assume 1D for now...
        coords = self.coords.flatten()
        if probability_density:
            weights = self.weights * self.descendant_weights
        else:
            weights = self.weights
        hist, bins = np.histogram(coords, weights=weights, bins = bins, density = True)
        bins -= (bins[1] - bins[0]) / 2
        return Plot(bins[:-1], hist, figure = figure)






