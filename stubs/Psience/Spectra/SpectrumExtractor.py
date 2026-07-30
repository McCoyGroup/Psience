import numpy as np
import McUtils.Devutils as dev
import McUtils.Numputils as nput
from McUtils.ExternalPrograms import PILInterface
from McUtils.Plots import ColorPalette
__all__ = ['SpectrumExtractor']

class SpectrumExtractor:

    def __init__(self, image_data, color_space='rgb', operation_color_space='lab'):
        """
        **LLM Docstring**

        Store image pixel data for spectrum extraction, converting it into the desired working color space if it isn't already in that space.

        :param image_data: the image pixel data, as a `(channels, height, width)` array
        :type image_data: np.ndarray
        :param color_space: the color space `image_data` is currently given in
        :type color_space: str
        :param operation_color_space: the color space to convert the image into (and perform subsequent operations in)
        :type operation_color_space: str
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def from_pil(cls, pil_image, color_space=None, **etc):
        """
        **LLM Docstring**

        Build a `SpectrumExtractor` from a PIL image, inferring the source color space from the image's mode if not given and dropping any alpha channel.

        :param pil_image: the PIL image to extract data from
        :type pil_image: object
        :param color_space: the color space of the image; inferred from `pil_image.mode` if not given (`'l'` maps to `'grayscale'`, `'rgba'` maps to `'rgb'`)
        :type color_space: str | None
        :param etc: extra options forwarded to the constructor
        :type etc: dict
        :return: the constructed extractor
        :rtype: SpectrumExtractor
        """
        ...

    @classmethod
    def from_file(cls, file, color_space=None, **etc):
        """
        **LLM Docstring**

        Build a `SpectrumExtractor` by loading an image from a file, via `PILInterface.from_file` and `from_pil`.

        :param file: the path to the image file
        :type file: str
        :param color_space: the color space of the image, forwarded to `from_pil`
        :type color_space: str | None
        :param etc: extra options forwarded to `from_pil`
        :type etc: dict
        :return: the constructed extractor
        :rtype: SpectrumExtractor
        """
        ...

    @classmethod
    def from_url(cls, file, color_space=None):
        """
        **LLM Docstring**

        Build a `SpectrumExtractor` by downloading an image from a URL, via `PILInterface.from_url` and `from_pil`.

        :param file: the URL to download the image from
        :type file: str
        :param color_space: the color space of the image, forwarded to `from_pil`
        :type color_space: str | None
        :return: the constructed extractor
        :rtype: SpectrumExtractor
        """
        ...
    default_tolerances = [25, 25, 25]

    def find_pixels(self, color, distance_tol=None, tolerances=None, color_space='rgb', search_color_space=None, image_range=None):
        """
        **LLM Docstring**

        Find the pixel coordinates in the image whose color is within a per-channel tolerance of a target color, optionally restricted to a sub-region of the image (given as absolute pixel bounds or fractional bounds of the image dimensions).

        :param color: the target color to search for, as a color string or a tuple of channel values in `color_space`
        :type color: str | tuple
        :param distance_tol: a single Euclidean color-distance tolerance; not currently supported
        :type distance_tol: float | None
        :param tolerances: per-channel tolerance values; defaults to `self.default_tolerances`
        :type tolerances: Iterable[float] | None
        :param color_space: the color space `color` is given in
        :type color_space: str
        :param search_color_space: the color space to perform the comparison in; defaults to `self.color_space`
        :type search_color_space: str | None
        :param image_range: a `(y_range, x_range)` (or single-axis) restriction on the search region, each as `(min, max)` in absolute pixels or as fractions of the corresponding image dimension
        :type image_range: tuple | None
        :return: the `(row_indices, col_indices)` of matching pixels, as returned by `np.where`
        :rtype: tuple[np.ndarray, np.ndarray]
        :raises NotImplementedError: if `distance_tol` is given (Euclidean-distance-based matching isn't supported)
        """
        ...

    def find_spectrum_lines(self, pixel_positions, max_pixel_distance=0.005, min_line_cutoff=0.5, smoothing=True, line_split_cutoff=5, allow_disjoint=False, spectrum_direction='up'):
        """
        **LLM Docstring**

        Group a set of matched pixel positions into distinct spectral "lines" (traces), by scanning column by column, splitting each column's matched rows into contiguous runs, and greedily attaching each run to the nearest existing line within `max_pixel_distance` (or starting a new line otherwise), then filtering out short lines, smoothing/averaging duplicate y-values per column, and (optionally) applying a moving-average smoothing window.

        :param pixel_positions: the `(row_indices, col_indices)` of matched pixels, as returned by `find_pixels`
        :type pixel_positions: tuple[np.ndarray, np.ndarray]
        :param max_pixel_distance: the maximum column+row distance for a new run to be attached to an existing line; if less than 1, treated as a fraction of the image width
        :type max_pixel_distance: float
        :param min_line_cutoff: the minimum number of points required to keep a line; if less than 1, treated as a fraction of the image width
        :type min_line_cutoff: float
        :param smoothing: whether to average duplicate y-values per column; if an integer greater than 1, additionally applies a moving-average smoothing window of that width
        :type smoothing: bool | int
        :param line_split_cutoff: the row-distance gap above which a column's matched rows are split into separate runs
        :type line_split_cutoff: float
        :param allow_disjoint: whether more than one run per column can be attached to the same line
        :type allow_disjoint: bool
        :param spectrum_direction: `'up'` or `'down'`, controlling the order runs within a column are considered (affects which line a given run preferentially attaches to)
        :type spectrum_direction: str
        :return: the extracted lines, each as an `(2, npoints)` array of `(x, y)` coordinates, sorted from longest to shortest
        :rtype: list[np.ndarray]
        """
        ...
    default_merge_tolerances = [10, 10, 10]

    def get_dominant_colors(self, bins=255, color_space='lab', min_counts=2, merge_tolerances=None):
        """
        **LLM Docstring**

        Find the most common colors in the image by histogramming its pixels in a given color space, then merging together histogram bins whose colors are within `merge_tolerances` of each other.

        :param bins: the number of histogram bins per channel
        :type bins: int
        :param color_space: the color space to histogram in
        :type color_space: str
        :param min_counts: the minimum pixel count for a histogram bin to be considered
        :type min_counts: int
        :param merge_tolerances: per-channel tolerance for merging nearby histogram bins together; defaults to `self.default_merge_tolerances`
        :type merge_tolerances: Iterable[float] | None
        :return: a mapping from color (as a channel-value tuple) to pixel count, sorted from most to least common
        :rtype: dict
        """
        ...

    def identify_frame_x_boundaries(self, pixel_positions, min_line_cutoff=0.5, frame_gap_cutoff=0.05):
        """
        **LLM Docstring**

        Identify the horizontal extent of a plot "frame" (e.g. an axis border) within a set of matched pixel positions, by finding the first column with an unusually long run of matched pixels (more than `min_line_cutoff`) and then scanning forward until a gap larger than `frame_gap_cutoff` (or another long-run column) is found.

        :param pixel_positions: the `(column_indices, row_indices)` of matched pixels to scan, grouped by column
        :type pixel_positions: tuple[np.ndarray, np.ndarray]
        :param min_line_cutoff: the minimum run length (in the row/second array) to consider a column part of the frame; if less than 1, treated as a fraction of the image width
        :type min_line_cutoff: float
        :param frame_gap_cutoff: the maximum number of consecutive short-run columns tolerated before considering the frame ended; if less than 1, treated as a fraction of the image width
        :type frame_gap_cutoff: float
        :return: `(min_x, max_x)`, the identified frame boundary column indices
        :rtype: tuple[int, int]
        """
        ...

    def identify_frame_boundaries(self, pixel_positions, min_line_cutoffs=0.5, frame_gap_cutoffs=0.05, identified_components=True):
        """
        **LLM Docstring**

        Identify both the horizontal and vertical extent of a plot frame within a set of matched pixel positions, via `identify_frame_x_boundaries` applied along each axis in turn (or skipped, falling back to the raw min/max, for axes not requested via `identified_components`).

        :param pixel_positions: the `(x_positions, y_positions)` of matched pixels
        :type pixel_positions: tuple[np.ndarray, np.ndarray]
        :param min_line_cutoffs: the minimum run-length cutoff(s) for the x/y frame-boundary detection; a single value is used for both axes
        :type min_line_cutoffs: float | tuple[float, float]
        :param frame_gap_cutoffs: the gap-tolerance cutoff(s) for the x/y frame-boundary detection
        :type frame_gap_cutoffs: float | tuple[float, float]
        :param identified_components: whether to actually detect the frame boundary along the x/y axis (`True` for both, or a per-axis pair); axes not detected fall back to the raw pixel min/max
        :type identified_components: bool | tuple[bool, bool]
        :return: `(x_range, y_range)`, each a `(min, max)` pair
        :rtype: tuple[tuple[int, int], tuple[int, int]]
        """
        ...

    def prune_straight_vertical_pixels(self, pixel_positions, min_line_cutoff=0.5):
        """
        **LLM Docstring**

        Intended to filter out columns whose matched-pixel run is shorter than `min_line_cutoff` (e.g. to discard sparse, likely-noise vertical pixel runs) and return their concatenated `(x, y)` positions. As written, `np.full(x, len(Y))` swaps the intended arguments -- it builds an array of length `x` filled with the value `len(Y)`, rather than an array of length `len(Y)` filled with the value `x` -- so the returned `x`-coordinate array does not actually correspond one-to-one with `new_y`; this looks like a bug.

        :param pixel_positions: the `(column_indices, row_indices)` of matched pixels, grouped by column
        :type pixel_positions: tuple[np.ndarray, np.ndarray]
        :param min_line_cutoff: the run-length cutoff below which a column's pixels are kept; if less than 1, treated as a fraction of the image width
        :type min_line_cutoff: float
        :return: the concatenated (but, due to the bug above, mismatched) `x`/`y` position arrays
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def extract_spectra(self, color=None, use_exact_color=False, image_range=None, max_dominant_percentage=0.2, allow_grayscale=False, color_space='rgb', dominant_color_space='lab', dominant_bins=255, min_dominant_component=50, extract_lines=True, prune_frame_components=True, frame_line_cutoffs=0.5, spectrum_direction='up', x_range=(0, 1), y_range=(0, 1), preserve_x_range=True, preserve_y_range=False, use_entire_pixel_range=True, return_color_code=True, **opts):
        """
        **LLM Docstring**

        Top-level pipeline for extracting spectral line data from a plot image: resolves the color to search for (either the given color, or the most prominent non-background dominant color found via `get_dominant_colors`), finds all matching pixels, optionally prunes out pixels belonging to the plot frame/border, groups the remaining pixels into line traces (via `find_spectrum_lines`), and rescales the resulting line coordinates into the requested `x_range`/`y_range` (preserving the overall span across all lines by default, or normalizing each line independently).

        :param color: the target line color to extract; if `None` (or `use_exact_color` is `False`), a dominant color is chosen automatically
        :type color: str | tuple | None
        :param use_exact_color: whether to search for `color` exactly rather than resolving it to the nearest dominant color
        :type use_exact_color: bool
        :param image_range: a sub-region of the image to restrict the pixel search to, forwarded to `find_pixels`
        :type image_range: tuple | None
        :param max_dominant_percentage: the maximum fraction of image pixels a color can occupy and still be considered a candidate line color (excludes large background regions)
        :type max_dominant_percentage: float
        :param allow_grayscale: whether low-saturation (grayscale-like) dominant colors are allowed as candidates
        :type allow_grayscale: bool
        :param color_space: the color space `color` is given in
        :type color_space: str
        :param dominant_color_space: the color space used for dominant-color detection/matching
        :type dominant_color_space: str
        :param dominant_bins: the number of histogram bins per channel for dominant-color detection
        :type dominant_bins: int
        :param min_dominant_component: the minimum pixel count for a dominant-color candidate
        :type min_dominant_component: int
        :param extract_lines: whether to group matched pixels into line traces (via `find_spectrum_lines`) rather than returning the raw matched pixel positions
        :type extract_lines: bool
        :param prune_frame_components: whether to prune pixels outside the detected plot frame, along one or both axes
        :type prune_frame_components: bool | tuple[bool, bool]
        :param frame_line_cutoffs: the run-length cutoff(s) used for frame-boundary detection, forwarded to `identify_frame_boundaries`
        :type frame_line_cutoffs: float | tuple[float, float]
        :param spectrum_direction: forwarded to `find_spectrum_lines`
        :type spectrum_direction: str
        :param x_range: the target output range to rescale the x-coordinates into
        :type x_range: tuple[float, float] | int | None
        :param y_range: the target output range to rescale the y-coordinates into
        :type y_range: tuple[float, float] | int | None
        :param preserve_x_range: whether to preserve the relative x-alignment across all extracted lines (rescaling them together) rather than normalizing each line independently
        :type preserve_x_range: bool
        :param preserve_y_range: whether to preserve the relative y-alignment across all extracted lines
        :type preserve_y_range: bool
        :param use_entire_pixel_range: whether to use the full matched-pixel range (rather than just the range spanned by the extracted lines) when computing the rescaling span
        :type use_entire_pixel_range: bool | tuple[bool, bool]
        :param return_color_code: whether to convert the resolved color to an RGB hex code before returning it
        :type return_color_code: bool
        :param opts: extra options, split between dominant-color-detection, line-extraction, and pixel-search options
        :type opts: dict
        :return: `(color, lines)` -- the resolved color (or its RGB code) and the extracted/rescaled line data, or `(color, [])` if no matching pixels were found
        :rtype: tuple
        """
        ...