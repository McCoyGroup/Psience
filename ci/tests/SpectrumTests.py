from Peeves.TestUtils import *
from unittest import TestCase
from Psience.Spectra import *
import McUtils.Plots as plt
import sys, os, numpy as np

class SpectrumTests(TestCase):

    @validationTest
    def test_LorentzSpectrum(self):

        spec = DiscreteSpectrum(
            [ 3685.825401, 2706.131282, 1383.456813, 7202.853894, 5323.915794, 2749.223062, 6377.965209, 5044.828562, 4072.476551],
            [ 35.6790758, 13.99412984, 54.64570094, 1.177683293, 0.544362793, 2.17610184, 0.219226886, 1.617797499, 2.331812556]
        )

        ploot = spec.plot()

        brood = spec.broaden('lorentzian', breadth=[500/x for x in spec.intensities])

        brood.plot(figure=ploot, axes_labels=["Frequency cm$^{-1}$", "Intensity km mol$^{-1}$"], image_size=500)

        ploot.show()

    @debugTest
    def test_ExtractSpec(self):
        extractor = SpectrumExtractor.from_file(
            TestManager.test_data('spec_split.png')
        )
        color, specs = extractor.extract_spectra(
            'red',
            x_range=[2800, 3000],
            tolerances=[25, 25, 25],
            spectrum_direction='down',
            extract_lines=True,
            smoothing=False,
            max_pixel_distance=.02
        )

        # print(specs)
        fig = None
        for s in specs:
            fig = plt.Plot(
                *s,
                color=plt.ColorPalette.rgb_code(
                    plt.ColorPalette.color_convert(color, 'lab', 'rgb')
                ),
                figure=fig
                # marker_size=1
            )
        fig.show()
