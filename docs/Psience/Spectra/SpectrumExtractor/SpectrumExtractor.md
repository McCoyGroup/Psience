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


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.from_pil" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_pil(cls, pil_image, color_space=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L20)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L20?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.from_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_file(cls, file, color_space=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L38)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L38?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.from_url" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_url(cls, file, color_space=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L43?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.find_pixels" class="docs-object-method">&nbsp;</a> 
```python
find_pixels(self, color, distance_tol=None, tolerances=None, color_space='rgb', search_color_space=None, image_range=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L49)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L49?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.find_spectrum_lines" class="docs-object-method">&nbsp;</a> 
```python
find_spectrum_lines(self, pixel_positions, max_pixel_distance=0.005, min_line_cutoff=0.5, smoothing=True, line_split_cutoff=5, allow_disjoint=False, spectrum_direction='up'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L110)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L110?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.get_dominant_colors" class="docs-object-method">&nbsp;</a> 
```python
get_dominant_colors(self, bins=255, color_space='lab', min_counts=2, merge_tolerances=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L201)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L201?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.identify_frame_x_boundaries" class="docs-object-method">&nbsp;</a> 
```python
identify_frame_x_boundaries(self, pixel_positions, min_line_cutoff=0.5, frame_gap_cutoff=0.05): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L273)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L273?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.identify_frame_boundaries" class="docs-object-method">&nbsp;</a> 
```python
identify_frame_boundaries(self, pixel_positions, min_line_cutoffs=0.5, frame_gap_cutoffs=0.05, identified_components=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L306)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L306?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.prune_straight_vertical_pixels" class="docs-object-method">&nbsp;</a> 
```python
prune_straight_vertical_pixels(self, pixel_positions, min_line_cutoff=0.5): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L349)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L349?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.SpectrumExtractor.SpectrumExtractor.extract_spectra" class="docs-object-method">&nbsp;</a> 
```python
extract_spectra(self, color=None, use_exact_color=False, image_range=None, max_dominant_percentage=0.2, allow_grayscale=False, color_space='rgb', dominant_color_space='lab', dominant_bins=255, min_dominant_component=50, extract_lines=True, prune_frame_components=True, frame_line_cutoffs=0.5, spectrum_direction='up', x_range=(0, 1), y_range=(0, 1), preserve_x_range=True, preserve_y_range=False, use_entire_pixel_range=True, return_color_code=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L360)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/SpectrumExtractor/SpectrumExtractor.py#L360?message=Update%20Docs)]
</div>
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