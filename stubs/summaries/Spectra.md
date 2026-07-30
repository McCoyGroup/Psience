### `BaseSpectrum.py` — Provides a general base spectrum class that can be extended to new fancy spectral forms
  - **class `BaseSpectrum`**
    > Base class to support spectral operation
    - `__init__(frequencies, intensities, **meta)`
    - `take_subspectrum(pos)` — Takes a subset of frequencies/intensities specified by `pos`
    - `frequency_filter(freq_min, freq_max)` — Filters by frequencies >= `freq_min` and <= `freq_max`
    - `intensity_filter(int_min, int_max)` — Filters by intensities >= `int_min` and <= `int_max`
    - `save(file)` — Saves the spectrum in JSON format
    - `load(file)` — Saves a spectrum from a JSON file
    - `plot(**opts)` — A stub so that subclasses can implement their own `plot` methods
  - **class `DiscreteSpectrum`** (BaseSpectrum)
    > Concrete implementation of `BaseSpectrum` that exists
    > solely to allow for plotting and broadening.
    - `from_transition_moments(frequencies, transition_moments, **meta)` — Assumes frequencies and transition moments in a.u.
    - `from_raman_moments(frequencies, transition_polarizabilities, pump_frequency=0, **meta)` — Build a discrete Raman spectrum from transition frequencies and polarizability transition moments (…
    - `normalize(which=None)` — Build a copy of the spectrum with intensities normalized by either the overall maximum intensity or…
    - `plot(figure=None, plot_style=None, **opts)` — :param figure: figure to plot the spectrum on
    - `broaden(breadth=10, *, broadening_type='gaussian', noising_function=None)` — Applies a broadening to the spectrum
  - **class `ContinuousSpectrum`** (BaseSpectrum)
    > Concrete implementation of `BaseSpectrum` that exists
    > solely to allow for plotting & maybe some day for interchange with like experimental formats
    - `plot(figure=None, filled=False, plot_style=None, **opts)` — :param figure: figure to plot the spectrum on
  - **class `BroadenedSpectrum`** (BaseSpectrum)
    > A stick spectrum with associated broadening function
    - `__init__(frequencies, intensities, broadening_type='gaussian', breadth=10, noising_function=None, **meta)`
    - `plot(step_size=0.5, freq_min=None, freq_max=None, figure=None, plot_style=None, filled=False, adjust_width=True, renormalize=True, **opts)` — Applies the broadening then plots it using `McUtils.Plots.Plot`
  - **class `MultiSpectrum`**
    > A wrapper for multiple spectra, really just for the plotting convenience
    - `__init__(spectra, **meta)`
    - `frequency_filter(freq_min, freq_max)` — Filters by frequencies >= `freq_min` and <= `freq_max`
    - `intensity_filter(int_min, int_max)` — Filters by intensities >= `int_min` and <= `int_max`
    - `plot(figure=None, **opts)` — A just plots all the spectra on the same figure

### `HarmonicSpectrum.py`
  - **class `HarmonicSpectrum`** (DiscreteSpectrum)
    - `from_normal_modes(nms, dipole_derivatives, **opts)` — Build a harmonic (double-harmonic-approximation) IR spectrum from a set of normal modes and their C…
    - `raman_from_modes(nms, polarizability_derivatives, **opts)` — Build a harmonic Raman spectrum from a set of normal modes and their Cartesian polarizability deriv…
    - `from_mol(mol, modes=None, dipole_derivatives=None, **opts)` — Build a harmonic IR spectrum for a molecule, using its normal modes and Cartesian dipole derivative…
    - `raman_from_mol(mol, modes=None, polarizability_derivatives=None, **opts)` — Build a harmonic Raman spectrum for a molecule, using its normal modes and Cartesian polarizability…

### `Multidimensional.py`
  - **class `TwoDimensionalSpectrum`**
    > Base class to support spectral operation
    - `__init__(freq1, freq2, intensities, **meta)`
    - `take_subspectrum(sample_x, sample_y)` — Takes a subset of frequencies/intensities specified by `pos`
    - `frequency_filter(freq_span_x, freq_span_y)` — Restrict the spectrum to a rectangular window of frequency values along both axes, via `take_subspe…
    - `intensity_filter(int_min, int_max)` — Restrict the spectrum to the rectangular index range spanning the points whose intensity falls with…
    - `clip(int_min, int_max, clip_abs=True)` — Clip the spectrum's intensities to a `[int_min, int_max]` range, either clipping the absolute magni…
    - `plot(plot_filled=True, contour_line_style=None, figure=None, symmetric_range=True, remove_baseline=True, vmin=None, vmax=None, levels=None, **opts)` — Render the 2D spectrum as a filled and/or line contour plot, optionally removing the median baselin…

### `SpectrumExtractor.py`
  - **class `SpectrumExtractor`**
    - `__init__(image_data, color_space='rgb', operation_color_space='lab')`
    - `from_pil(pil_image, color_space=None, **etc)` — Build a `SpectrumExtractor` from a PIL image, inferring the source color space from the image's mod…
    - `from_file(file, color_space=None, **etc)` — Build a `SpectrumExtractor` by loading an image from a file, via `PILInterface.from_file` and `from…
    - `from_url(file, color_space=None)` — Build a `SpectrumExtractor` by downloading an image from a URL, via `PILInterface.from_url` and `fr…
    - `find_pixels(color, distance_tol=None, tolerances=None, color_space='rgb', search_color_space=None, image_range=None)` — Find the pixel coordinates in the image whose color is within a per-channel tolerance of a target c…
    - `find_spectrum_lines(pixel_positions, max_pixel_distance=0.005, min_line_cutoff=0.5, smoothing=True, line_split_cutoff=5, allow_disjoint=False, spectrum_direction='up')` — Group a set of matched pixel positions into distinct spectral "lines" (traces), by scanning column…
    - `get_dominant_colors(bins=255, color_space='lab', min_counts=2, merge_tolerances=None)` — Find the most common colors in the image by histogramming its pixels in a given color space, then m…
    - `identify_frame_x_boundaries(pixel_positions, min_line_cutoff=0.5, frame_gap_cutoff=0.05)` — Identify the horizontal extent of a plot "frame" (e.g.
    - `identify_frame_boundaries(pixel_positions, min_line_cutoffs=0.5, frame_gap_cutoffs=0.05, identified_components=True)` — Identify both the horizontal and vertical extent of a plot frame within a set of matched pixel posi…
    - `prune_straight_vertical_pixels(pixel_positions, min_line_cutoff=0.5)` — Intended to filter out columns whose matched-pixel run is shorter than `min_line_cutoff` (e.g.
    - `extract_spectra(color=None, use_exact_color=False, image_range=None, max_dominant_percentage=0.2, allow_grayscale=False, color_space='rgb', dominant_color_space='lab', dominant_bins=255, min_dominant_component=50, extract_lines=True, prune_frame_components=True, frame_line_cutoffs=0.5, spectrum_direction='up', x_range=(0, 1), y_range=(0, 1), preserve_x_range=True, preserve_y_range=False, use_entire_pixel_range=True, return_color_code=True, **opts)` — Top-level pipeline for extracting spectral line data from a plot image: resolves the color to searc…