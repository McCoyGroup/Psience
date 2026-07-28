## <a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor">SpectrumExtractor</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor.py#L11)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor.py#L11?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
default_tolerances: list
default_merge_tolerances: list
```
<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, image_data, color_space='rgb', operation_color_space='lab'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor.py#L13)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor.py#L13?message=Update%20Docs)]
</div>
**LLM Docstring**

Store image pixel data for spectrum extraction, converting it into the desired working color space if it isn't already in that space.
  - `image_data`: `np.ndarray`
    > the image pixel data, as a `(channels, height, width)` array
  - `color_space`: `str`
    > the color space `image_data` is currently given in
  - `operation_color_space`: `str`
    > the color space to convert the image into (and perform subsequent operations in)
  - `:returns`: `None`
    > None


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.from_pil" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_pil(cls, pil_image, color_space=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L34)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L34?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `SpectrumExtractor` from a PIL image, inferring the source color space from the image's mode if not given and dropping any alpha channel.
  - `pil_image`: `object`
    > the PIL image to extract data from
  - `color_space`: `str | None`
    > the color space of the image; inferred from `pil_image.mode` if not given (`'l'` maps to `'grayscale'`, `'rgba'` maps to `'rgb'`)
  - `etc`: `dict`
    > extra options forwarded to the constructor
  - `:returns`: `SpectrumExtractor`
    > the constructed extractor


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.from_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_file(cls, file, color_space=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L66)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L66?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `SpectrumExtractor` by loading an image from a file, via `PILInterface.from_file` and `from_pil`.
  - `file`: `str`
    > the path to the image file
  - `color_space`: `str | None`
    > the color space of the image, forwarded to `from_pil`
  - `etc`: `dict`
    > extra options forwarded to `from_pil`
  - `:returns`: `SpectrumExtractor`
    > the constructed extractor


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.from_url" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_url(cls, file, color_space=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L85)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L85?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `SpectrumExtractor` by downloading an image from a URL, via `PILInterface.from_url` and `from_pil`.
  - `file`: `str`
    > the URL to download the image from
  - `color_space`: `str | None`
    > the color space of the image, forwarded to `from_pil`
  - `:returns`: `SpectrumExtractor`
    > the constructed extractor


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.find_pixels" class="docs-object-method">&nbsp;</a> 
```python
find_pixels(self, color, distance_tol=None, tolerances=None, color_space='rgb', search_color_space=None, image_range=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L103)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L103?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the pixel coordinates in the image whose color is within a per-channel tolerance of a target color, optionally restricted to a sub-region of the image (given as absolute pixel bounds or fractional bounds of the image dimensions).
  - `color`: `str | tuple`
    > the target color to search for, as a color string or a tuple of channel values in `color_space`
  - `distance_tol`: `float | None`
    > a single Euclidean color-distance tolerance; not currently supported
  - `tolerances`: `Iterable[float] | None`
    > per-channel tolerance values; defaults to `self.default_tolerances`
  - `color_space`: `str`
    > the color space `color` is given in
  - `search_color_space`: `str | None`
    > the color space to perform the comparison in; defaults to `self.color_space`
  - `image_range`: `tuple | None`
    > a `(y_range, x_range)` (or single-axis) restriction on the search region, each as `(min, max)` in absolute pixels or as fractions of the corresponding image dimension
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > the `(row_indices, col_indices)` of matching pixels, as returned by `np.where`


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.find_spectrum_lines" class="docs-object-method">&nbsp;</a> 
```python
find_spectrum_lines(self, pixel_positions, max_pixel_distance=0.005, min_line_cutoff=0.5, smoothing=True, line_split_cutoff=5, allow_disjoint=False, spectrum_direction='up'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L185)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L185?message=Update%20Docs)]
</div>
**LLM Docstring**

Group a set of matched pixel positions into distinct spectral "lines" (traces), by scanning column by column, splitting each column's matched rows into contiguous runs, and greedily attaching each run to the nearest existing line within `max_pixel_distance` (or starting a new line otherwise), then filtering out short lines, smoothing/averaging duplicate y-values per column, and (optionally) applying a moving-average smoothing window.
  - `pixel_positions`: `tuple[np.ndarray, np.ndarray]`
    > the `(row_indices, col_indices)` of matched pixels, as returned by `find_pixels`
  - `max_pixel_distance`: `float`
    > the maximum column+row distance for a new run to be attached to an existing line; if less than 1, treated as a fraction of the image width
  - `min_line_cutoff`: `float`
    > the minimum number of points required to keep a line; if less than 1, treated as a fraction of the image width
  - `smoothing`: `bool | int`
    > whether to average duplicate y-values per column; if an integer greater than 1, additionally applies a moving-average smoothing window of that width
  - `line_split_cutoff`: `float`
    > the row-distance gap above which a column's matched rows are split into separate runs
  - `allow_disjoint`: `bool`
    > whether more than one run per column can be attached to the same line
  - `spectrum_direction`: `str`
    > `'up'` or `'down'`, controlling the order runs within a column are considered (affects which line a given run preferentially attaches to)
  - `:returns`: `list[np.ndarray]`
    > the extracted lines, each as an `(2, npoints)` array of `(x, y)` coordinates, sorted from longest to shortest


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.get_dominant_colors" class="docs-object-method">&nbsp;</a> 
```python
get_dominant_colors(self, bins=255, color_space='lab', min_counts=2, merge_tolerances=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L298)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L298?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the most common colors in the image by histogramming its pixels in a given color space, then merging together histogram bins whose colors are within `merge_tolerances` of each other.
  - `bins`: `int`
    > the number of histogram bins per channel
  - `color_space`: `str`
    > the color space to histogram in
  - `min_counts`: `int`
    > the minimum pixel count for a histogram bin to be considered
  - `merge_tolerances`: `Iterable[float] | None`
    > per-channel tolerance for merging nearby histogram bins together; defaults to `self.default_merge_tolerances`
  - `:returns`: `dict`
    > a mapping from color (as a channel-value tuple) to pixel count, sorted from most to least common


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.identify_frame_x_boundaries" class="docs-object-method">&nbsp;</a> 
```python
identify_frame_x_boundaries(self, pixel_positions, min_line_cutoff=0.5, frame_gap_cutoff=0.05): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L386)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L386?message=Update%20Docs)]
</div>
**LLM Docstring**

Identify the horizontal extent of a plot "frame" (e.g. an axis border) within a set of matched pixel positions, by finding the first column with an unusually long run of matched pixels (more than `min_line_cutoff`) and then scanning forward until a gap larger than `frame_gap_cutoff` (or another long-run column) is found.
  - `pixel_positions`: `tuple[np.ndarray, np.ndarray]`
    > the `(column_indices, row_indices)` of matched pixels to scan, grouped by column
  - `min_line_cutoff`: `float`
    > the minimum run length (in the row/second array) to consider a column part of the frame; if less than 1, treated as a fraction of the image width
  - `frame_gap_cutoff`: `float`
    > the maximum number of consecutive short-run columns tolerated before considering the frame ended; if less than 1, treated as a fraction of the image width
  - `:returns`: `tuple[int, int]`
    > `(min_x, max_x)`, the identified frame boundary column indices


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.identify_frame_boundaries" class="docs-object-method">&nbsp;</a> 
```python
identify_frame_boundaries(self, pixel_positions, min_line_cutoffs=0.5, frame_gap_cutoffs=0.05, identified_components=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L433)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L433?message=Update%20Docs)]
</div>
**LLM Docstring**

Identify both the horizontal and vertical extent of a plot frame within a set of matched pixel positions, via `identify_frame_x_boundaries` applied along each axis in turn (or skipped, falling back to the raw min/max, for axes not requested via `identified_components`).
  - `pixel_positions`: `tuple[np.ndarray, np.ndarray]`
    > the `(x_positions, y_positions)` of matched pixels
  - `min_line_cutoffs`: `float | tuple[float, float]`
    > the minimum run-length cutoff(s) for the x/y frame-boundary detection; a single value is used for both axes
  - `frame_gap_cutoffs`: `float | tuple[float, float]`
    > the gap-tolerance cutoff(s) for the x/y frame-boundary detection
  - `identified_components`: `bool | tuple[bool, bool]`
    > whether to actually detect the frame boundary along the x/y axis (`True` for both, or a per-axis pair); axes not detected fall back to the raw pixel min/max
  - `:returns`: `tuple[tuple[int, int], tuple[int, int]]`
    > `(x_range, y_range)`, each a `(min, max)` pair


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.prune_straight_vertical_pixels" class="docs-object-method">&nbsp;</a> 
```python
prune_straight_vertical_pixels(self, pixel_positions, min_line_cutoff=0.5): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L492)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L492?message=Update%20Docs)]
</div>
**LLM Docstring**

Intended to filter out columns whose matched-pixel run is shorter than `min_line_cutoff` (e.g. to discard sparse, likely-noise vertical pixel runs) and return their concatenated `(x, y)` positions. As written, `np.full(x, len(Y))` swaps the intended arguments -- it builds an array of length `x` filled with the value `len(Y)`, rather than an array of length `len(Y)` filled with the value `x` -- so the returned `x`-coordinate array does not actually correspond one-to-one with `new_y`; this looks like a bug.
  - `pixel_positions`: `tuple[np.ndarray, np.ndarray]`
    > the `(column_indices, row_indices)` of matched pixels, grouped by column
  - `min_line_cutoff`: `float`
    > the run-length cutoff below which a column's pixels are kept; if less than 1, treated as a fraction of the image width
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > the concatenated (but, due to the bug above, mismatched) `x`/`y` position arrays


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.extract_spectra" class="docs-object-method">&nbsp;</a> 
```python
extract_spectra(self, color=None, use_exact_color=False, image_range=None, max_dominant_percentage=0.2, allow_grayscale=False, color_space='rgb', dominant_color_space='lab', dominant_bins=255, min_dominant_component=50, extract_lines=True, prune_frame_components=True, frame_line_cutoffs=0.5, spectrum_direction='up', x_range=(0, 1), y_range=(0, 1), preserve_x_range=True, preserve_y_range=False, use_entire_pixel_range=True, return_color_code=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L515)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L515?message=Update%20Docs)]
</div>
**LLM Docstring**

Top-level pipeline for extracting spectral line data from a plot image: resolves the color to search for (either the given color, or the most prominent non-background dominant color found via `get_dominant_colors`), finds all matching pixels, optionally prunes out pixels belonging to the plot frame/border, groups the remaining pixels into line traces (via `find_spectrum_lines`), and rescales the resulting line coordinates into the requested `x_range`/`y_range` (preserving the overall span across all lines by default, or normalizing each line independently).
  - `color`: `str | tuple | None`
    > the target line color to extract; if `None` (or `use_exact_color` is `False`), a dominant color is chosen automatically
  - `use_exact_color`: `bool`
    > whether to search for `color` exactly rather than resolving it to the nearest dominant color
  - `image_range`: `tuple | None`
    > a sub-region of the image to restrict the pixel search to, forwarded to `find_pixels`
  - `max_dominant_percentage`: `float`
    > the maximum fraction of image pixels a color can occupy and still be considered a candidate line color (excludes large background regions)
  - `allow_grayscale`: `bool`
    > whether low-saturation (grayscale-like) dominant colors are allowed as candidates
  - `color_space`: `str`
    > the color space `color` is given in
  - `dominant_color_space`: `str`
    > the color space used for dominant-color detection/matching
  - `dominant_bins`: `int`
    > the number of histogram bins per channel for dominant-color detection
  - `min_dominant_component`: `int`
    > the minimum pixel count for a dominant-color candidate
  - `extract_lines`: `bool`
    > whether to group matched pixels into line traces (via `find_spectrum_lines`) rather than returning the raw matched pixel positions
  - `prune_frame_components`: `bool | tuple[bool, bool]`
    > whether to prune pixels outside the detected plot frame, along one or both axes
  - `frame_line_cutoffs`: `float | tuple[float, float]`
    > the run-length cutoff(s) used for frame-boundary detection, forwarded to `identify_frame_boundaries`
  - `spectrum_direction`: `str`
    > forwarded to `find_spectrum_lines`
  - `x_range`: `tuple[float, float] | int | None`
    > the target output range to rescale the x-coordinates into
  - `y_range`: `tuple[float, float] | int | None`
    > the target output range to rescale the y-coordinates into
  - `preserve_x_range`: `bool`
    > whether to preserve the relative x-alignment across all extracted lines (rescaling them together) rather than normalizing each line independently
  - `preserve_y_range`: `bool`
    > whether to preserve the relative y-alignment across all extracted lines
  - `use_entire_pixel_range`: `bool | tuple[bool, bool]`
    > whether to use the full matched-pixel range (rather than just the range spanned by the extracted lines) when computing the rescaling span
  - `return_color_code`: `bool`
    > whether to convert the resolved color to an RGB hex code before returning it
  - `opts`: `dict`
    > extra options, split between dominant-color-detection, line-extraction, and pixel-search options
  - `:returns`: `tuple`
    > `(color, lines)` -- the resolved color (or its RGB code) and the extracted/rescaled line data, or `(color, [])` if no matching pixels were found
 </div>
</div>












---


<div markdown="1" class="text-secondary">
<div class="container">
  <div class="row">
   <div class="col" markdown="1">
**Feedback**   
</div>
   <div class="col" markdown="1">
**Examples**   
</div>
   <div class="col" markdown="1">
**Templates**   
</div>
   <div class="col" markdown="1">
**Documentation**   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[Bug](https://github.com/McCoyGroup/Psience/issues/new?title=Documentation%20Improvement%20Needed)/[Request](https://github.com/McCoyGroup/Psience/issues/new?title=Example%20Request)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor.py#L11?message=Update%20Docs)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>
</div>