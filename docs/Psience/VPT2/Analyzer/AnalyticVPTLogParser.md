## <a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser">AnalyticVPTLogParser</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L736)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L736?message=Update%20Docs)]
</div>

Log-file parser for text logs produced by `AnalyticVPTRunner.run_VPT(..., logger=<path>)`,
following the same overall convention as `VPTAnalyzerLogParser` (a `LogParser` subclass
exposing cached, lazily-parsed properties backed by small `StringLineByLineReader` block
parsers) -- but adapted to `AnalyticVPTRunner`'s own, structurally different log format:

- There's no single outer `">>--- Starting Perturbation Theory Runner ---<<"` banner
  wrapping the whole run the way there is for the classic `VPTRunner`; results live
  directly under a top-level `"Running VPT"` block instead (see `tree`).
- Its per-initial-state tables -- the `"Transition Moments:"` block, and the IR-spectrum-
  style table that immediately follows it -- are genuinely `"|"`-delimited (unlike the
  classic parser's fixed-width columns), and use a different state-label convention:
  `"()"` for the ground state, `"k(q)"` for `q` quanta in 1-indexed mode `k`, concatenated
  for combination states (e.g. `"1(1)2(1)"` for one quantum each in modes 1 and 2) --
  rather than the classic parser's space-separated per-mode quanta (`"0 0 0"`, `"1 1 0"`,
  ...). See `parse_state_label`.
- The `"Transition Moments:"` table already prints the fully-summed *total* transition
  moment for each transition as its own leading `x, y, z` columns, ahead of the same
  per-order breakdown the classic log only ever gives in pieces -- so there's no analytic-
  log equivalent of the `load_term_counts` column-slicing bug that had to be fixed for
  `VPTAnalyzerLogParser` (see `claude_drafts/vpt_analyzer_log_parsing_fixes.patch`):
  `transition_moment_corrections` here just reads that leading total straight off the
  table, no order-by-order recombination required.

Only the two tables needed to reconstruct a `transition_dict` for 2D-IR are implemented
here -- `spectra` and `transition_moment_corrections`, mirroring the two properties
`prep_vpt_response_data_from_log` actually calls on `VPTAnalyzerLogParser` for the classic
format. Everything else (`energies` and friends, `deperturbed_spectra`,
`deperturbed_transition_moment_corrections`, ...) raises `NotImplementedError("TBD")` for
now -- mirroring the many properties that are already effectively unsupported (bare
`KeyError`) for `VPTAnalyzerLogParser`+`VPTResultsLoader`'s classic `"log_file"` dispatch,
just made explicit here rather than left as a surprise.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
AnalyticBlockParser: AnalyticBlockParser
```
<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, log_file, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L771)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L771?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a log-file parser for an `AnalyticVPTRunner` run's text log output, accepting
either a file path or an in-memory `io.StringIO` -- mirrors
`VPTAnalyzerLogParser.__init__` exactly, since the underlying `LogParser` needs a real
file either way.
  - `log_file`: `str | io.StringIO`
    > the path to the log file, or an in-memory string buffer containing its contents
  - `opts`: `dict`
    > extra options forwarded to the base `LogParser.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.parse_state_label" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
parse_state_label(cls, label, ndim=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L809)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L809?message=Update%20Docs)]
</div>
**LLM Docstring**

Parse an `AnalyticVPTRunner`-style state label (e.g. `"1(1)2(1)"`, or `"()"` for the
ground state) into an excitation-quanta tuple (e.g. `(0, 1, 1)`), matching the tuple
format the classic parser's labels are converted to elsewhere (e.g.
`NonlinearResponse._parse_vpt_state_label`). Nothing in this class does this conversion
internally -- `spectra`/`transition_moment_corrections` hand back the raw label strings,
same as `VPTAnalyzerLogParser` does for its own differently-formatted labels -- so a
caller reconstructing a `transition_dict` needs to convert them itself, same as it
already has to for the classic format.

The mode-index-to-tuple-position mapping is *not* the naive `label_index - 1`: these
labels are formatted by `Psience.BasisReps.Util.StateMaker.parse_state` in its default
`mode='low-high'` convention, which numbers a state's modes from the *end* of the
excitation tuple -- position 1 is the tuple's last entry, position `ndim` is its first
-- so mode index `k` (1-indexed, as printed) lands at tuple position `ndim - k`
(confirmed against `StateMaker.make_state`'s own inverse, `state[-i] = q`, and verified
against `water_freq_response_analytic.json`: label `"3(1)"` is the fundamental at
1572.7 cm^-1, matching classic-format tuple `(1, 0, 0)`, not `(0, 0, 1)`).
  - `label`: `str`
    > the raw state label, e.g. `"1(1)2(1)"` or `"()"`
  - `ndim`: `int | None`
    > the number of modes; if not given, inferred as the highest mode index that
    appears in `label` (or 1, for the ground state alone) -- which only gives the right
    answer if the *lowest-position* mode (position 1, the tuple's last entry) shows up
    excited somewhere in `label`, so callers that need a reliable `ndim` should pass it
    explicitly instead of relying on the fallback
  - `:returns`: `tuple[int]`
    > excitation-quanta tuple, e.g. `(0, 1, 1)`


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.tree" class="docs-object-method">&nbsp;</a> 
```python
@property
tree(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L854)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L854?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) content of the log's top-level `"Running VPT"` block: a flat list mixing
named sub-block dicts (e.g. `{"Transition Moments:": [...]}`, for content
`AnalyticVPTRunner` wraps in a `logger.block(tag=...)`) with bare, one-string-per-line
entries (for content it prints via a plain `logger.log_print(...)` instead, such as the
energies table and the IR-spectrum-style table -- see `_raw_lines_after`). Unlike
`VPTAnalyzerLogParser.tree`, there's no multi-top-level-block collapsing to do here:
`AnalyticVPTRunner`'s log always has exactly one `"Running VPT"` entry among its
top-level blocks (alongside the separate `"calculating G/potential/dipole derivatives"`
setup blocks, which aren't needed for anything implemented so far), so this just
reaches directly into it.
  - `:returns`: `list`
    > the parsed `"Running VPT"` block content


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.reformat_spectrum_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
reformat_spectrum_block(cls, sb): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1011)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1011?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert one parsed spectrum sub-block (`{initial_label: [(final_label, value_tokens), ...]}`)
into a `{'states':..., 'harmonic':..., 'anharmonic':...}` dict -- the analytic-log
equivalent of `VPTAnalyzerLogParser.reformat_spec_block`, for the same 4-numeric-column
(harmonic freq/intensity, anharmonic freq/intensity) table shape.
  - `sb`: `dict`
    > the parsed block data, as returned by `AnalyticBlockParser`
  - `:returns`: `dict | list[dict]`
    > a single dict (if there was only one initial state) or list of dicts, each with `'states'`, `'harmonic'`, and `'anharmonic'` entries


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.reformat_transition_moment_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
reformat_transition_moment_block(cls, sb): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1048)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1048?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert one parsed transition-moment sub-block (`{initial_label: [(final_label, value_tokens), ...]}`)
into a `{'states':..., 'transition_moment':...}` dict, reading the already fully-summed
total `x, y, z` transition moment straight off the table's own leading 3 columns (see
this class's docstring for why no per-order recombination is needed here, unlike
`VPTAnalyzerLogParser.reformat_tm_block`).
  - `sb`: `dict`
    > the parsed block data, as returned by `AnalyticBlockParser`
  - `:returns`: `dict | list[dict]`
    > a single dict (if there was only one initial state) or list of dicts, each with `'states'` and `'transition_moment'` entries


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.parse_analytic_blocks" class="docs-object-method">&nbsp;</a> 
```python
parse_analytic_blocks(self, spec_str, reformatter): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1078)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1078?message=Update%20Docs)]
</div>
**LLM Docstring**

Parse a raw multi-initial-state table block of log text into a list of reformatted
per-block dicts, via `AnalyticBlockParser` and the given `reformatter`
(`reformat_spectrum_block` or `reformat_transition_moment_block`) -- the analytic-log
equivalent of `VPTAnalyzerLogParser.parse_spectrum_blocks`/`parse_tm_blocks`.
  - `spec_str`: `str`
    > the raw text of the table block
  - `reformatter`: `callable`
    > the classmethod to reformat each parsed sub-block with
  - `:returns`: `list`
    > the list of reformatted blocks, one per initial state


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1098)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1098?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) IR spectrum data parsed from the untagged table trailing the log's
`"Transition Moments:"` block, via `reformat_spectrum_block`. Needed for 2D-IR: gives
frequency + intensity per transition, per initial state -- the analytic-log equivalent
of `VPTAnalyzerLogParser.spectra`. Unlike that method, each returned dict also carries
an `'initial_state'` entry (the raw label, e.g. `"3(1)"` or `"()"`) rather than leaving
the caller to infer it, since there's no diagonal row here to infer it from.
  - `:returns`: `list[dict]`
    > the parsed per-initial-state spectrum data, each entry also keyed by `'initial_state'`


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1118)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1118?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) transition-dipole-moment data parsed from the log's `"Transition Moments:"`
block, via `reformat_transition_moment_block`. Needed for 2D-IR: gives the already
fully-summed `[x, y, z]` transition moment per transition, per initial state -- the
analytic-log equivalent of `VPTAnalyzerLogParser.transition_moment_corrections`, minus
the per-order recombination that one needs (see this class's docstring). Each returned
dict also carries an `'initial_state'` entry, same as `spectra` above.
  - `:returns`: `list[dict]`
    > the parsed per-initial-state transition-moment data, each entry also keyed by `'initial_state'`


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.harmonic_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
harmonic_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1147)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1147?message=Update%20Docs)]
</div>
**LLM Docstring**

Not yet implemented -- see this class's docstring.


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1158)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1158?message=Update%20Docs)]
</div>
**LLM Docstring**

Not yet implemented -- see this class's docstring.


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.zero_order_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1169)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1169?message=Update%20Docs)]
</div>
**LLM Docstring**

Not yet implemented -- see this class's docstring.


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.deperturbed_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1180)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1180?message=Update%20Docs)]
</div>
**LLM Docstring**

Not yet implemented -- see this class's docstring.


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.deperturbed_spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1191)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1191?message=Update%20Docs)]
</div>
**LLM Docstring**

Not yet implemented -- see this class's docstring.


<a id="Psience.VPT2.Analyzer.AnalyticVPTLogParser.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1202)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/AnalyticVPTLogParser.py#L1202?message=Update%20Docs)]
</div>
**LLM Docstring**

Not yet implemented -- see this class's docstring.
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analyzer/AnalyticVPTLogParser.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analyzer/AnalyticVPTLogParser.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analyzer/AnalyticVPTLogParser.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analyzer/AnalyticVPTLogParser.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L736?message=Update%20Docs)   
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