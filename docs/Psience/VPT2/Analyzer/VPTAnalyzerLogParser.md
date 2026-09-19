## <a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser">VPTAnalyzerLogParser</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L109?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
EnergiesBlockParser: EnergiesBlockParser
SpectrumBlockParser: SpectrumBlockParser
TransitionMomentBlockParser: TransitionMomentBlockParser
```
<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, log_file, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L110)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L110?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a log-file parser for a VPT run's text log output, accepting either a file path or an in-memory `io.StringIO` (which is first written out to a temporary file, since the underlying `LogParser` needs a real file).
  - `log_file`: `str | io.StringIO`
    > the path to the log file, or an in-memory string buffer containing its contents
  - `opts`: `dict`
    > extra options forwarded to the base `LogParser.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.tree" class="docs-object-method">&nbsp;</a> 
```python
@property
tree(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L147)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L147?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) parsed block-tree structure of the log file, collapsed down to just the "Computing PT corrections:" subtree (or otherwise condensed) if the raw parse produced multiple top-level blocks.

Every log produced by `VPTRunner.run_simple(..., logger=...)` is wrapped in exactly one outer
named block (`">>--- Starting Perturbation Theory Runner ---<<"`), so `to_tree()` almost always
returns a tree with a *single* top-level key. The `len(self._tree) > 1` branch below never used
to handle that case at all -- it only unwrapped a tree with *more than one* top-level entry --
which meant the named tables underneath that single banner block (`"IR Data"`, `"X/Y/Z Dipole
Contributions"`, etc) were never reachable and any log-based lookup of them (`.spectrum`,
`.transition_moment_corrections`, ...) failed with a bare `IndexError`, even for the one
pre-existing reference log fixture (`methanol_vpt_3.out`) shipped in `ci/tests/TestData`. This
now unwraps that single banner block first, before falling through to the pre-existing
multi-block handling (which still applies to logs that genuinely have more than one top-level
block once unwrapped, or that start directly with a "Computing PT corrections:" block and never
had an outer banner to begin with).
  - `:returns`: `object`
    > the parsed (and condensed) log-file tree


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.reformat_eng_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
reformat_eng_block(cls, sb): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L252)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L252?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert a parsed energies block (a list of `(state_label, value_tokens)` rows) into parallel arrays of state labels, harmonic energies, and anharmonic energies, handling the alternate row format that includes a placeholder `'-'` degeneracy-group marker column.
  - `sb`: `list`
    > the parsed block data, as returned by `EnergiesBlockParser`
  - `:returns`: `tuple[list[str], np.ndarray, np.ndarray]`
    > `(states, harm, anh)` -- the state labels and their harmonic/anharmonic energy arrays


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.parse_energies_blocks" class="docs-object-method">&nbsp;</a> 
```python
parse_energies_blocks(self, spec_str): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L278)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L278?message=Update%20Docs)]
</div>
**LLM Docstring**

Parse a raw energies-table block of log text into `(states, harmonic, anharmonic)` arrays, via `EnergiesBlockParser` and `reformat_eng_block`.
  - `spec_str`: `str`
    > the raw text of the energies block
  - `:returns`: `tuple[list[str], np.ndarray, np.ndarray]`
    > `(states, harm, anh)`


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.harmonic_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
harmonic_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L293)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L293?message=Update%20Docs)]
</div>
**LLM Docstring**

The `(states, harmonic_energies)` pair parsed from the log's "States Energies" (or, for a degenerate run, "Degenerate Energies") table block, cached after the first access.
  - `:returns`: `tuple[list[str], np.ndarray]`
    > `(states, harmonic_energies)`


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L311)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L311?message=Update%20Docs)]
</div>
**LLM Docstring**

The `(states, anharmonic_energies)` pair parsed from the log's "States Energies"/"Degenerate Energies" table block, cached after the first access.
  - `:returns`: `tuple[list[str], np.ndarray]`
    > `(states, anharmonic_energies)`


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.zero_order_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L329)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L329?message=Update%20Docs)]
</div>
**LLM Docstring**

The `(states, harmonic_energies)` pair parsed from the log's "States Energies"/"Degenerate Energies" table block; equivalent to `harmonic_energies`.
  - `:returns`: `tuple[list[str], np.ndarray]`
    > `(states, harmonic_energies)`


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.deperturbed_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L347)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L347?message=Update%20Docs)]
</div>
**LLM Docstring**

The `(states, deperturbed_energies)` pair parsed from the log's "Deperturbed Energies" table block (falling back to "States Energies" if no degenerate treatment was run), cached after the first access.
  - `:returns`: `tuple[list[str], np.ndarray]`
    > `(states, deperturbed_energies)`


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L365)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L365?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) IR spectrum data parsed from the log's "IR Data" block, via `parse_spectrum_blocks`.
  - `:returns`: `list[dict]`
    > the parsed per-initial-state spectrum data


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.deperturbed_spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L379)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L379?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) deperturbed IR spectrum data parsed from the log's "Deperturbed IR Data" block, via `parse_spectrum_blocks`.
  - `:returns`: `list[dict]`
    > the parsed per-initial-state deperturbed spectrum data


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.reformat_spec_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
reformat_spec_block(cls, sb): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L475)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L475?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert a parsed spectrum block (per-initial-state groups of `(state_label, value_tokens)` rows) into a list of per-initial-state dicts holding the final-state labels and their harmonic/anharmonic `(frequency, intensity)` value pairs.
  - `sb`: `dict`
    > the parsed block data, as returned by `SpectrumBlockParser`
  - `:returns`: `dict | list[dict]`
    > a single dict (if there was only one initial state) or list of dicts, each with `'states'`, `'harmonic'`, and `'anharmonic'` entries


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.parse_spectrum_blocks" class="docs-object-method">&nbsp;</a> 
```python
parse_spectrum_blocks(self, spec_str): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L507)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L507?message=Update%20Docs)]
</div>
**LLM Docstring**

Parse a raw multi-initial-state spectrum-table block of log text into a list of reformatted per-block spectrum dicts, via `SpectrumBlockParser` and `reformat_spec_block`.
  - `spec_str`: `str`
    > the raw text of the spectrum block
  - `:returns`: `list`
    > the list of reformatted spectrum blocks


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.load_term_counts" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
load_term_counts(cls, nterms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L591)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L591?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute how many perturbative-order transition-moment correction terms are packed into each
column-chunk of a transition-moment data row, using a `SymmetricGroupGenerator`'s cumulative
term totals (the running total of `(i, j, k)` triples with `i + j + k <= order`, across
increasing `order`) to derive each order's own term count.

`_indexer._cumtotals` is a list of *cumulative* boundaries (e.g. `[0, 1, 4, 10, 20]` for 3
modes: 1 term through order 0, 4 through order 1, 10 through order 2, 20 through order 3).
This used to return those cumulative boundaries themselves, filtered to `< nterms`, and
`reformat_tm_block` then used each one directly as a *chunk width* -- but a cumulative total
is not a per-order width, and the final boundary that actually reaches `nterms` was always
excluded by the strict `<` (a row with exactly `nterms` columns needs the boundary *at*
`nterms` included, not just those strictly below it). For a `nterms=10` row (all 10 order-0
through order-2 correction terms, i.e. `1 + 3 + 6`), this returned `[0, 1, 4]` and got used as
chunk widths `0, 1, 4` -- consuming only the first 5 of the row's 10 printed columns and
silently dropping the rest, which is why `VPTAnalyzerLogParser`-reconstructed transition
moments for combination-band/overtone transitions came out ~5-30% off even once the multi-
initial-state parsing bugs elsewhere in this class were fixed (see
`claude_drafts/vpt_analyzer_log_parsing_fixes.patch`). This now includes the boundary
that reaches `nterms` itself, and converts the cumulative boundaries into genuine per-order
term-count *widths* via consecutive differences, so `reformat_tm_block` consumes every printed
column.
  - `nterms`: `int`
    > the total number of numeric columns present in a transition-moment data row
  - `:returns`: `list[int]`
    > the number of correction terms belonging to each perturbative order, in order (e.g. `[1, 3, 6]` for a 10-column row)


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.reformat_tm_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
reformat_tm_block(cls, sb): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L632)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L632?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert a parsed transition-moment-correction block (per-initial-state groups of `(state_label, value_tokens)` rows) into a list of per-initial-state dicts holding the final-state labels and their per-order correction arrays (split from the flat value-token list using `load_term_counts`).
  - `sb`: `dict`
    > the parsed block data, as returned by `TransitionMomentBlockParser`
  - `:returns`: `dict | list[dict]`
    > a single dict (if there was only one initial state) or list of dicts, each with `'states'` and `'corrections'` (a list of per-order arrays) entries


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.parse_tm_blocks" class="docs-object-method">&nbsp;</a> 
```python
parse_tm_blocks(self, spec_str): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L665)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L665?message=Update%20Docs)]
</div>
**LLM Docstring**

Parse a raw multi-initial-state transition-moment-table block of log text into a list of reformatted per-block correction dicts, via `TransitionMomentBlockParser` and `reformat_tm_block`.
  - `spec_str`: `str`
    > the raw text of the transition-moment block
  - `:returns`: `list`
    > the list of reformatted transition-moment blocks


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L681)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L681?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) transition-dipole-moment corrections parsed from the log's per-axis "X/Y/Z Dipole Contributions" blocks, combined into a single list of per-initial-state dicts each holding the `[x, y, z]` correction data.
  - `:returns`: `list[dict]`
    > the combined per-initial-state transition-moment correction data


<a id="Psience.VPT2.Analyzer.VPTAnalyzerLogParser.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L708)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.py#L708?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) deperturbed transition-dipole-moment corrections parsed from the log's per-axis "X/Y/Z Deperturbed Dipole Contributions" blocks, combined into a single list of per-initial-state dicts each holding the `[x, y, z]` correction data.
  - `:returns`: `list[dict]`
    > the combined per-initial-state deperturbed transition-moment correction data
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analyzer/VPTAnalyzerLogParser.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L109?message=Update%20Docs)   
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