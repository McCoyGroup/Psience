from Peeves.TestUtils import *
from unittest import TestCase
from Psience.Spectra import *
import sys, os, numpy as np

class SpectrumTests(TestCase):

    @debugTest
    def test_LorentzSpectrum(self):

        spec = DiscreteSpectrum(
            [ 3685.825401, 2706.131282, 1383.456813, 7202.853894, 5323.915794, 2749.223062, 6377.965209, 5044.828562, 4072.476551],
            [ 35.6790758, 13.99412984, 54.64570094, 1.177683293, 0.544362793, 2.17610184, 0.219226886, 1.617797499, 2.331812556]
        )

        ploot = spec.plot()

        brood = spec.broaden('lorentzian', breadth=[500/x for x in spec.intensities])

        brood.plot(figure=ploot, axes_labels=["Frequency cm$^{-1}$", "Intensity km mol$^{-1}$"], image_size=500)

        ploot.show()
