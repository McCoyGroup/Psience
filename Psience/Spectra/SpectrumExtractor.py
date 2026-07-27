import numpy as np
import McUtils.Devutils as dev
import McUtils.Numputils as nput
from McUtils.ExternalPrograms import PILInterface
from McUtils.Plots import ColorPalette

__all__ = [
    "SpectrumExtractor"
]

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
        image_data = np.asanyarray(image_data)
        if color_space != operation_color_space:
            image_data = ColorPalette.color_convert(image_data, color_space, operation_color_space)
        self.img = image_data
        self.color_space = operation_color_space

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
        if color_space is None:
            color_space = pil_image.mode.lower()
            if color_space == 'l':
                color_space = 'grayscale'
            elif color_space == 'rgba':
                color_space = 'rgb'

        image_data = np.moveaxis(np.asarray(pil_image), -1, 0)[:3]
        if color_space == 'rgb':
            image_data = np.clip(image_data, 0, 255)

        return cls(
            image_data,
            color_space=color_space,
            **etc
        )
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
        pil = PILInterface.from_file(file)

        return cls.from_pil(pil, color_space=color_space, **etc)
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
        pil = PILInterface.from_url(file)
        return cls.from_pil(pil, color_space=color_space)

    default_tolerances = [25, 25, 25]
    def find_pixels(
            self,
            color,
            distance_tol=None,
            tolerances=None,
            color_space='rgb',
            search_color_space=None,
            image_range=None
    ):
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
        if isinstance(color, str):
            color = ColorPalette.parse_color_string(color)

        if search_color_space is None:
            search_color_space = self.color_space
        color = ColorPalette.color_convert(color, color_space, search_color_space)

        image_data = self.img
        if search_color_space is not None and search_color_space != self.color_space:
            image_data = ColorPalette.color_convert(image_data, self.color_space, search_color_space)

        if image_range is not None:
            y,x = image_range
            if nput.is_numeric(x):
                y,x = None, [y,x]
            if x is not None:
                xm, xM = x
                if not nput.is_int(xm):
                    xm = int(image_data.shape[-2] * xm)
                if not nput.is_int(xM):
                    xM = int(image_data.shape[-2] * xM)
                image_data = image_data[:, xm:xM, :]
            if y is not None:
                ym, yM = y
                if not nput.is_int(ym):
                    ym = int(image_data.shape[-1] * ym)
                if not nput.is_int(yM):
                    yM = int(image_data.shape[-1] * yM)
                image_data = image_data[:, :, ym:yM]

        if distance_tol is not None:
            raise NotImplementedError("finding pixels from color distance not supported")
        else:
            if tolerances is None:
                tolerances = self.default_tolerances

            mask = np.full(image_data[0].shape, True)
            for channel, center, tol in zip(
                image_data,
                color,
                tolerances
            ):
                mask = np.logical_and(
                    mask,
                    np.logical_and(
                        channel >= center - tol,
                        channel <= center + tol
                    )
                )

        return np.where(mask)

    def find_spectrum_lines(self,
                            pixel_positions,
                            max_pixel_distance=.005,
                            min_line_cutoff=0.5,
                            smoothing=True,
                            line_split_cutoff=5,
                            allow_disjoint=False,
                            spectrum_direction='up'
                            ):
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

        if max_pixel_distance < 1:
            max_pixel_distance = max_pixel_distance * max(self.img.shape[-1:])

        lines = []
        for x,Y in zip(*nput.group_by(pixel_positions[1], pixel_positions[0])[0]):
            Y = np.sort(Y) # just in case
            # if smoothing:
            y_diffs = np.diff(Y)
            split_points = np.where(np.abs(y_diffs) > line_split_cutoff)
            if len(split_points) > 0 and len(split_points[0]) > 0:
                Y = [
                    g
                    for g in np.split(Y, split_points[0]+1)
                ]
            else:
                Y = [Y]

            if spectrum_direction.lower() == 'down':
                Y = reversed(Y)

            if len(lines) == 0:
                lines.extend([[x, y]] for y in Y)
            else:
                add_pos = set()
                for y in Y:
                    for n,l in enumerate(lines):
                        if n in add_pos: continue

                        xx, yy = l[-1]
                        dx = x-xx
                        if abs(dx) > max_pixel_distance: continue

                        dy = min([abs(y[0] - yy[0]), abs(y[-1] - yy[-1])]) # check boundaries first
                        if abs(dx)+abs(dy) > max_pixel_distance:
                            dy = np.min(np.abs(yy[:, np.newaxis] - y[np.newaxis, :]))
                        if abs(dx)+abs(dy) < max_pixel_distance:
                            l.append([x, y])
                            if not allow_disjoint:
                                add_pos.add(n)
                            break
                    else:
                        lines.append([[x, y]])


        if min_line_cutoff < 1:
            min_line_cutoff = self.img.shape[-1] * min_line_cutoff
            lines = [
                ll for ll in lines
                if len(ll) > min_line_cutoff
            ]


        if smoothing:
            lines = [
                [[x, np.average(y)] for x, y in l]
                for l in lines
            ]
            if nput.is_int(smoothing) and smoothing is not True:
                lines = [np.array(l).T for l in lines]
                lines = [
                    np.array([
                        np.convolve(x, np.ones(smoothing), mode='valid') / smoothing,
                        np.convolve(y, np.ones(smoothing), mode='valid') / smoothing
                    ]).T
                    for x,y in lines
                ]
        else:
            lines = [
                sum(
                    ([[x, yy] for yy in y] for x, y in l),
                    []
                )
                for l in lines
            ]

        return [
            np.array(l).T
            for l in sorted(lines, key=lambda x:-len(x))
        ]

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
        image_data = self.img
        if color_space != self.color_space:
            image_data = ColorPalette.color_convert(image_data, self.color_space, color_space)
        dists, (x, y, z) = np.histogramdd(np.moveaxis(image_data, 0, -1).reshape(-1, 3), bins=bins)
        take = np.where(dists > min_counts)
        base_colors = {
            (
                (x[i] + x[i+1])/2,
                (y[j] + y[j+1])/2,
                (z[k] + z[k+1])/2,
            ):int(dists[i,j,k])
            for i,j,k in zip(*take)
        }

        base_colors = {
            k:base_colors[k]
            for k in sorted(base_colors, key=lambda x:-base_colors[x])
        }
        if merge_tolerances is None:
            merge_tolerances = self.default_merge_tolerances

        test_keys = list(base_colors.keys())
        merges = {}
        for n,c in enumerate(base_colors.keys()):
            for k in test_keys[:max(n-1, 0)]:
                if c == k or k in merges: continue
                if all(abs(k[i] - c[i]) <= merge_tolerances[i] for i in range(3)):
                    merges[c] = k
                    break

        for c,k in merges.items():
            base_colors[k] += base_colors[c]
            del base_colors[c]

        # print(
        #     sum(v for v in base_colors.values()) / (image_data.shape[-1] * image_data.shape[-2])
        # )

        return {
            k:base_colors[k]
            for k in sorted(base_colors, key=lambda x:-base_colors[x])
        }


        #TOOD: move this to OpenCV
        # all_pixels = np.moveaxis(self.img, 0, -1).reshape(-1, 3)
        # pix = {}
        # if tolerances is None:
        #     tolerances = self.default_tolerances
        # tolerances = np.asanyarray(tolerances)
        # m = 0
        # for p in all_pixels:
        #     old_pix = None
        #     new_pix = {}
        #     m += 1
        #     for k in pix:
        #         if all(k[i] - p[i] <= tolerances[i] for i in range(3)):
        #             n = pix[k]
        #             p = tuple((n*k[i] + p[i]) / (n+1) for i in range(3))
        #             new_pix[p] = n + 1
        #             old_pix = k
        #             break
        #     else:
        #         p = tuple(p)
        #         new_pix[p] = 1
        #     if old_pix is not None:
        #         del pix[old_pix]
        #     pix.update(new_pix)
        #
        # return pix

    def identify_frame_x_boundaries(self, pixel_positions, min_line_cutoff=.5, frame_gap_cutoff=.05):
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
        min_x = None
        cur_x = None
        max_X = None
        if min_line_cutoff < 1:
            min_line_cutoff = self.img.shape[-1] * min_line_cutoff
        if frame_gap_cutoff < 1:
            frame_gap_cutoff = self.img.shape[-1] * frame_gap_cutoff

        miss_counts = 0
        for x,Y in zip(*nput.group_by(pixel_positions[1], pixel_positions[0])[0]):
            if len(Y) > min_line_cutoff:
                if min_x is not None:
                    max_X = x - 1
                    break
                cur_x = x
            elif cur_x is not None:
                miss_counts += 1
                if miss_counts > frame_gap_cutoff:
                    min_x = cur_x + 1

        if min_x is None:
            if max_X is not None:
                if max_X < .5 * self.img.shape[-1]:
                    min_x = max_X
                    max_X = None
            else:
                min_x = np.min(pixel_positions[0])
        if max_X is None:
            max_X = np.max(pixel_positions[0])

        return min_x, max_X

    def identify_frame_boundaries(self, pixel_positions,
                                  min_line_cutoffs=.5,
                                  frame_gap_cutoffs=.05,
                                  identified_components=True
                                  ):
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
        if nput.is_numeric(min_line_cutoffs): min_line_cutoffs = [min_line_cutoffs, min_line_cutoffs]
        if nput.is_numeric(frame_gap_cutoffs): frame_gap_cutoffs = [frame_gap_cutoffs, frame_gap_cutoffs]
        x_cutoff, y_cutoff = min_line_cutoffs
        x_gap, y_gap = frame_gap_cutoffs

        if x_cutoff < 1:
            x_cutoff = self.img.shape[-2] * x_cutoff
        if y_cutoff < 1:
            y_cutoff = self.img.shape[-1] * y_cutoff

        if x_gap < 1:
            x_gap = self.img.shape[-2] * x_gap
        if y_gap < 1:
            y_gap = self.img.shape[-1] * y_gap

        if identified_components is True: identified_components = [identified_components, identified_components]

        find_x, find_y = identified_components
        if find_x:
            x_range = self.identify_frame_x_boundaries(
                (pixel_positions[0], pixel_positions[1]),
                min_line_cutoff=x_cutoff,
                frame_gap_cutoff=x_gap
            )
        else:
            x_range = [np.min(pixel_positions[0]), np.max(pixel_positions[0])]

        if find_y:
            y_range = self.identify_frame_x_boundaries(
                (pixel_positions[1], pixel_positions[0]),
                min_line_cutoff=y_cutoff,
                frame_gap_cutoff=y_gap
            )
        else:
            y_range = [np.min(pixel_positions[1]), np.max(pixel_positions[1])]

        return x_range, y_range

    def prune_straight_vertical_pixels(self, pixel_positions, min_line_cutoff=.5):
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
        new_x = []
        new_y = []
        if min_line_cutoff < 1:
            min_line_cutoff = self.img.shape[-1] * min_line_cutoff
        for x,Y in zip(*nput.group_by(pixel_positions[1], pixel_positions[0])[0]):
            if len(Y) > min_line_cutoff: continue
            new_x.append(np.full(x, len(Y)))
            new_y.append(Y)
        return np.concatenate(new_x), np.concatenate(new_y)

    def extract_spectra(self,
                        color=None,
                        use_exact_color=False,
                        image_range=None,
                        max_dominant_percentage=.2,
                        allow_grayscale=False,
                        color_space='rgb',
                        dominant_color_space='lab',
                        dominant_bins=255,
                        min_dominant_component=50,
                        extract_lines=True,
                        prune_frame_components=True,
                        frame_line_cutoffs=.5,
                        spectrum_direction='up',
                        x_range=(0, 1),
                        y_range=(0, 1),
                        preserve_x_range=True,
                        preserve_y_range=False,
                        use_entire_pixel_range=True,
                        return_color_code=True,
                        **opts
                        ):
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
        if isinstance(color, str):
            color = ColorPalette.parse_color_string(color)[:3]
            color_space = 'rgb'
        if color is None or not use_exact_color:
            opts = dev.OptionsSet(opts)
            dom_opts, opts = opts.split(None,
                                        {
                                            'bins',
                                            'merge_tolerances'
                                        }
                                        )
            doms = self.get_dominant_colors(dominant_bins,
                                            min_counts=min_dominant_component,
                                            color_space=dominant_color_space,
                                            **dom_opts)
            npixels = self.img.shape[-2] * self.img.shape[-1]
            max_npix = npixels * max_dominant_percentage
            if color is None:
                color = next((
                    d for d in doms
                    if (
                        doms[d] < max_npix
                        and (
                                allow_grayscale
                                or ColorPalette.color_convert(d, dominant_color_space, 'lch')[1] > .1
                        )
                    )
                ))
            else:
                c = ColorPalette.color_convert(color, color_space, dominant_color_space)
                color = next((
                    d for d in
                    sorted(doms,
                           key=lambda k: (np.linalg.norm(c - np.array(k))**2) * (doms[k] / npixels)
                           )
                ))
            color_space = dominant_color_space
        opts = dev.OptionsSet(opts)
        line_opts, pixel_opts = opts.split(None,
                                           {
                                               'max_pixel_distance',
                                               'min_line_cutoff',
                                               'smoothing',
                                               'line_split_cutoff',
                                               'allow_disjoint'
                                           }
                                           )
        y, x = self.find_pixels(color, image_range=image_range, color_space=color_space, **pixel_opts)
        if len(x) == 0 or len(y) == 0:
            return color, []

        if prune_frame_components:
            if prune_frame_components is True:
                prune_frame_components = [prune_frame_components, prune_frame_components]
            prune_x, prune_y = prune_frame_components
            x_bounds, y_bounds = self.identify_frame_boundaries(
                (x, y),
                identified_components=[prune_x, prune_y],
                min_line_cutoffs=frame_line_cutoffs
            )
            # raise Exception(x_bounds, y_bounds, (np.min(x), np.max(x)), (np.min(y), np.max(y)))
            mask = np.logical_and(
                np.logical_and(x_bounds[0] <= x,  x <= x_bounds[1]),
                np.logical_and(y_bounds[0] <= y,  y <= y_bounds[1])
            )

            x = x[mask]
            y = y[mask]

        y = self.img.shape[1] - y
        if extract_lines:
            lines = self.find_spectrum_lines((x, y), spectrum_direction=spectrum_direction, **line_opts)
        else:
            lines = [[x, y]]


        if extract_lines:
            if use_entire_pixel_range is True:
                use_entire_pixel_range = [True, False]
            if x_range is not None:
                if nput.is_int(x_range):
                    x_range = [0, x_range]

                if preserve_x_range:
                    if use_entire_pixel_range[0]:
                        x_span = [
                            np.min(x),
                            np.max(x)
                        ]
                    else:
                        x_mins = [
                            np.min(xx)
                            for xx, yy in lines
                            if len(xx) > 0
                        ]
                        if len(x_mins) == 0:
                            x_mins = [np.min(x)]

                        x_maxes = [
                            np.max(xx)
                            for xx, yy in lines
                            if len(xx) > 0
                        ]
                        if len(x_maxes) == 0:
                            x_maxes = [np.max(x)]

                        x_span = [
                            min(x_mins),
                            max(x_maxes)
                        ]
                    x_scale = (x_range[1] - x_range[0]) / (x_span[1] - x_span[0])
                    lines = [
                        [
                            (x - x_span[0]) * x_scale + x_range[0],
                            y
                        ]
                        for x,y in lines
                    ]
                else:
                    lines = [
                        [
                            (x - np.min(x)) * (x_range[1] - x_range[0]) / (np.max(x) - np.min(x)) + x_range[0],
                            y
                        ]
                        for x, y in lines
                    ]

            if y_range is not None:
                if nput.is_int(y_range):
                    y_range = [0, y_range]

                if preserve_y_range:
                    if use_entire_pixel_range[0]:
                        y_span = [
                            np.min(y),
                            np.max(y)
                        ]
                    else:
                        y_mins = [
                            np.min(yy)
                            for xx, yy in lines
                            if len(yy) > 0
                        ]
                        if len(y_mins) == 0:
                            y_mins = [np.min(y)]

                        y_maxes = [
                            np.max(yy)
                            for xx, yy in lines
                            if len(yy) > 0
                        ]
                        if len(y_maxes) == 0:
                            y_maxes = [np.max(y)]
                        y_span = [
                            min(y_mins),
                            max(y_maxes)
                        ]
                    y_scale = (y_range[1] - y_range[0]) / max([(y_span[1] - y_span[0]), 1])
                    lines = [
                        [
                            x,
                            (y - y_span[0]) * y_scale + y_range[0],
                        ]
                        for x,y in lines
                    ]
                else:
                    lines = [
                        [
                            x,
                            (y - np.min(y)) * (y_range[1] - y_range[0]) / (np.max(y) - np.min(y)) + y_range[0]
                        ]
                        for x, y in lines
                    ]

        if return_color_code:
            color = ColorPalette.rgb_code(
                ColorPalette.color_convert(
                    color, color_space, 'rgb'
                )
            )

        return color, lines

