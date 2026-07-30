# Psience -- how to use this directory

This is a **compressed, LLM-oriented reference** for the `Psience` codebase,
generated from its real source. It is NOT the library itself: nothing here
is meant to be imported or run as the real package. Use it to find out what
exists and how to call it correctly, then write code against the real
`Psience` package -- not against this directory.

## Read in this order

1. **`summaries/index.md`** -- one line per top-level package: what it's
   for, how many modules it has, and where its summary/stubs live. Start
   here to pick the right package.
2. **`summaries/<package>.md`** -- every public class/function/method in
   that package, with its exact call signature (parameter names and
   defaults, no type annotations) and a short one-line purpose. This is
   usually enough to write a correct call or import.
3. **`Psience/<package>/<module>.py`** -- the full stub for one module:
   real signatures, real docstrings (see caveats below), nothing
   abbreviated further. Open this when the summary line isn't enough --
   you need the full docstring, an edge case, or something the summary
   truncated.

Do not skip straight to step 3 for everything -- the summaries exist so you
usually don't have to.

## What's real vs. placeholder in the stub files

- **Signatures** (parameter names, defaults, `*args`/`**kwargs`) are exact
  copies of the real source.
- **Docstrings** are the real docstrings. Class docstrings are capped at
  their first paragraph (truncated ones are marked -- see "Compression
  tricks" below); everything else is complete.
- **Function/method bodies are placeholders** (`...`). They carry no
  information about real behavior, return values, or edge cases -- never
  infer runtime behavior from a stub body. The docstring is the only
  source of behavioral information here; if it doesn't say, this
  directory doesn't know either.
- **`__all__`** (however it's constructed in the real source -- plain
  assignment, `+=` accumulation, etc.) is preserved exactly, so each
  stubbed module's real export list is accurate.

## Compression tricks -- how to read them correctly

These exist purely to cut size; none of them drop information, but
misreading the compact form as the real API would be actively wrong:

- **Enum-style / flat constant runs.** A class with many one-line members
  (e.g. an `enum.Enum` with dozens of values) is collapsed into a single
  `_MEMBERS = {...}` dict, immediately preceded by a comment stating the
  REAL access pattern, e.g. `SomeEnum.Primary`, not
  `SomeEnum._MEMBERS['Primary']`. **Always follow that comment's stated
  access pattern** -- the dict is a compact data table, not the calling
  convention.
- **Large literal data** (big lookup tables, constant dictionaries). When
  a top-level literal is large, its value is replaced with a comment
  describing its keys/shape (e.g. `1069 keys: ['Hydrogen1', 'Hydrogen2', ...]`)
  so you know what's queryable without the full payload inflating this
  directory. The real values are NOT included anywhere in this directory at all -- only their keys/shape. If you need an actual value, say so rather than guessing one.
- **Truncated class docstrings** end with the marker
  `*(truncated — see stub for full docstring)*`. There is genuinely more
  text in the real source that isn't reproduced anywhere in this
  directory -- if the missing part matters, say you don't have it rather
  than guessing at what follows.

## Cross-package dependencies

`dependency_graph.json` (at the root of this directory) maps which packages, classes, methods, and functions depend on which others -- both within `Psience` and on external packages (stdlib and a handful of very common third-party packages are excluded as noise). Use it to trace connections between parts of the codebase before writing an example that spans more than one package, or to find everything that uses a particular class or method. Resolution is best-effort static analysis, refined by live introspection where possible -- treat it as a strong hint, not a guarantee, especially for dynamic or conditional imports.

## Freshness

This is a point-in-time snapshot. It can drift out of date relative to the
real `Psience` source. If you need current, authoritative behavior (not
just "does this function/class exist and roughly how do I call it"),
verify against the real source rather than treating this directory as
ground truth.
