# Analytic expression libraries

Use the normal analytic runner to build or reuse the mathematical and prepared
libraries. They contain expressions for the formal requests exercised by the run.

```python
from Psience.VPT2 import AnalyticVPTRunner

results = AnalyticVPTRunner.run_simple(
    system, states,
    ir_file="expressions/library.ir.tar",
    compiled_file="compiled/library.plan.tar",
    storage_backend="pickle",  # use "portable" for non-pickle interchange
    evaluation_folder="evaluations/run-001",
    **your_existing_options
)
```

The first run builds missing expressions. Later runs reuse matching entries,
including across physical systems. Transfer the compiled file and pass only
`compiled_file` to the recipient's `run_simple` call. The IR file is optional;
retaining it permits independent compilation into either backend on another
machine. Use `compilation_options={"reuse_only": True}` to require a complete
prepared library with no expression building or compilation.

Both backends save the same fully prepared plans. New entries default to
`portable`; omitting `storage_backend` preserves existing entries. Selecting a
backend explicitly converts existing entries without recompilation. Libraries
can be extended with new corrections, operator requests or quantum changes.
The `expressions_file` symbolic checkpoint and `expressions_folder` legacy
cache keep their earlier meanings.

Evaluation opens artifacts read-only and writes optional checkpoints only in
`evaluation_folder`, which must be outside the artifact folders. Omitting that
option leaves no persistent evaluation output. File-backed evaluation executes
the stored full-expression programs for either controller evaluation policy,
without runtime symbolic projection or implicit physical-system pruning.

## Explicit stage control

The lower-level interfaces remain available when the stages run separately:

```python
from Psience.VPT2 import AnalyticExpressionStore, CompiledAnalyticExpression

store = AnalyticExpressionStore()
with store.export_ir(wrapper, "expressions/energy.ir.tar", storage_backend="portable"):
    pass
compiled = store.compile_ir("expressions/energy.ir.tar", "compiled/energy.plan.tar",
                            storage_backend="pickle")
compiled = CompiledAnalyticExpression.load("compiled/energy.plan.tar")
values = compiled.evaluate(state_perms, coefficient_expansions, frequencies,
                           evaluation_folder="evaluations/run-001")
```

`compile_ir` also accepts a library produced by `run_simple`; it preserves the
formal request labels and version metadata needed by the recipient runner.
`save(..., storage_backend=...)` converts IR or prepared files without compiling.
`pack_ir` packs named simplified expressions into a single library.

Pickle is for trusted Python use with compatible installed packages. Portable
files contain documented semantic records and primitive NumPy buffers, with no
pickled classes or generated machine code. Wider-than-binary64/complex128
precision and unsupported custom mathematical types currently fail explicitly.

## Wire format

The following schema is experimental and versioned independently from producer
source hashes. It describes the non-pickle model and its Python storage adapter.

```text
Experimental analytic expression IR, revision 1

This draft schema is not yet a stable public interchange commitment. The two
artifact kinds share mathematical field definitions but have separate graph
semantics. Storage is independently selectable: portable (default) or pickle.

Container
  Standard, uncompressed tar using the existing McUtils embedded-index format.
  First member: .__mcutils_index__/locator, payload 1024 bytes. Its two 512-byte
  slots hold an ASCII SHA256 digest, newline, JSON pointer and NUL padding.
  Pointer: version=1, generation, absolute offset, size, end, sha256.
  The pointer addresses .__mcutils_index__/index.jsonl. Its first line carries
  version=2, previous-generation pointer (or null), and replaced prefixes;
  remaining rows are [name, absolute_data_offset,size,serialized_TarInfo_JSON].
  save() consolidates the current live generation. Ordinary tar readers see all
  mathematical members normally. Indexed read_only=True access needs no adjacent
  SQLite database, lock file, extraction or writes. No new IR-specific tar
  container is introduced. Historical first-index prototype output is superseded.

  expression/manifest.json describes a single expression.
  expression/part-NNNNNN-rR.json describes a bounded part of rank R.
  expression/part-NNNNNN-rR-buffer-DDD.npy contains pooled primitive data.
  Directory builds use the same mathematical filenames without expression/.

  A collection instead has root manifest.json with kind
  "analytic-expression-ir-collection" or "compiled-analytic-expression-collection",
  catalog_version:1, complete:true, entries:{label:"entries/000000",...}.
  Each prefix contains an independent manifest, parts and buffers. Labels never
  become member paths. Select an entry before loading; a sole entry auto-selects.

Manifest
  kind: "analytic-expression-ir"; ir_version: 1; complete: true
  parts: file, rank, sha256, bytes and symbolic-count descriptors
  change, baseline_shift, num_fixed, projection, degree/index bounds
  specializations: explicit active_operators and mode_dimension restrictions,
    null by default; no molecule-specific numerical values are captured
  operator_keys: the formal operator universe before optional pruning, allowing
    an evaluator to reject incompatible active-operator specializations
  implementation/signature/options/timings: producer/build provenance; source
    hashes and producer Python/NumPy versions are not reader compatibility gates
  storage_backend: "portable" or "pickle"; absent means portable for old files.
  storage_version: 1 for portable, 2 for current pickle; older pickle version 1
    (pickled generic record plus .npy pools) remains readable/upgradable.
  Each part declares file, physical sha256, bytes, buffers and content_sha256.
    Physical checksums validate the stored bytes. content_sha256 hashes the same
    typed node/array record model across backends; it is not a canonical algebraic
    equivalence hash. Source IR content hashes define compiled logical identity,
    independently of the selected input/output backend. Storage-specific build
    directories use separate keys to avoid collisions.
  Part bytes count metadata plus buffers, not tar/index padding overhead.

Storage modes
  Portable version 1 uses the JSON graph and pooled .npy members described below.
    No pickle data or private class names occur in this interchange representation.
  Pickle version 2 uses a protocol-5 envelope:
    {format:semantic_record_format,version:record_version,
     representation:"python",payload:mathematical_expression_or_prepared_payload}.
    A private Pickler adapter sends only the same whitelisted node fields to the
    shared node restorer. Arbitrary object dictionaries, symbolic/runtime caches
    and generated machine code are excluded. Python memoization preserves sharing.
    Arrays are inline (dtype string,shape,C-order bytes,scalar flag), explicitly
    little endian and restored as immutable views over bytes. Source hash salts
    and table caches are rebuilt by the same constructors as portable loading.
    Exceptionally deep Python graphs use representation:"graph", payload:the flat
    semantic graph record, plus .npy buffers, instead of changing recursion limits.
    This fallback is still pickle storage and applies only to known node types.
  Native pickle is for trusted Python use and requires compatible installed
    Psience/McUtils code. Portable storage is the documented external-library
    format. Both support directory staging and one-file indexed tar publication.
    A collection may mix storage modes; the individual manifests declare them.
  Changing storage decodes/restores fields and rewrites their representation; it
    does not simplify, lower, re-plan or invoke an expression compiler. Single
    artifacts stay single artifacts and collection labels are preserved.
    Equal-backend/version repacking streams declared files without deserialization
    and verifies the copied part checksum before publication.

Semantic graph record
  format: "psience-analytic-expression"; version: 2 (reads earlier revision 1)
  root: encoded value
  nodes: array of typed records, referenced by integer positions
  References may point forward and sharing is explicit. Cycles, unknown kinds
  and out-of-range references are rejected. No private class names are stored.

Values
  null, JSON strings and booleans represent themselves.
  Compact revision 2 uses bare JSON integers for absolute values <= 2**53-1;
    these are exact in ordinary binary64-based readers too.
  Bare JSON arrays mean tuples. Python lists keep their explicit list tag.
  {"integer":"decimal"} preserves arbitrary integer precision.
  {"float":"hex"} preserves a Python binary64 number, signed zero and IEEE
    special values through float.hex/fromhex, rather than JSON NaN extensions.
  {"complex":["real hex","imaginary hex"]}, {"bytes":"hex"}
  {"scalar":[explicit_numpy_dtype, encoded_primitive]} preserves NumPy scalars.
  {"tuple":[values...]}, {"list":[values...]}
    The tuple tag is retained for reading revision 1; new writers use bare arrays.
  {"mapping":[[encoded_key,encoded_value]...]} preserves typed tuple/tag keys.
  {"ref":node_index} refers to a typed mathematical node or array.

Array node
  kind: "array"; file; dtype; shape; offset; count; sha256
  The .npy file is a 1D pool of one explicit primitive dtype. offset/count are
  element positions; slice then reshape to shape in C order. Pools consolidate
  small factors to avoid a tar member for every array. Numeric byte order is
  explicitly little endian. Object/structured arrays and floating precision
  beyond binary64/complex128 are rejected in this initial version. There is no
  pickle payload in portable mode, raw platform-long-double assumption or lossy coercion.

Other nodes have kind plus fields (an encoded tuple/value):
  polynomial-factor: (coefficient_array, integer_shift)
  polynomial-axis: tuple of factor references, multiplied along one mode
  polynomial-term: (tuple_of_axes, steps); divisor is sqrt(2)**steps
  polynomial-path: (weighted_term_items, operation_or_null, ndim, order, reduced)
    weighted_term_items: tuple of (term_reference, scalar_weight)
    operation: (opcode, operands)
    opcodes: add, linear, scale, shift, remap, permutation_sum, project, mul
    add: tuple of child references, summed
    linear: tuple of (child,weight), summed
    scale: (child,weight)
    shift: (child,shift_vector), evaluate child at occupations + shift
    remap: (child,pull,output_ndim), gather occupations by pull
    permutation_sum: (child,permutations,pulls), sum child over the pulls
    project: (child,keep,condensed), fill child occupations with zero except
      selected axes keep, which receive the parent input occupations
    mul: (left,right,left_inds,right_inds,left_pull,right_pull), evaluate the
      children with the respective pulls and multiply; inds retain algebraic
      contraction provenance
    Paths are not flattened into dense multivariate polynomial coefficients.
  product-polynomial: (axis_coefficient_arrays, prefactor, steps)
  polynomial-sum: (polynomial_list, prefactor, reduced)
  energy-sum: (denominator_factor_product_to_numerator_mapping, prefactor,
               reduced, ndim)
  tensor-sum: (coefficient_product_to_polynomial_mapping, prefactor,reduced,ndim)
    A coefficient descriptor is (expansion_order, tensor_identity, mode_indices...)
    tensor_identity uses native IDs V,G,U,Z,u,M,LV,LG.
    Harmonic rows retain their side tag and integer frequency deltas. An empty
    factor product is a unit denominator; tagged factors are never cast to ints.
  sqrt-change: (polynomial, change_array, baseline_shift_array)
  projection: (kind P/Q/R/V, formal_slot)
  change-domain: (retain_boolean, endpoint_0_or_1)
  condition-endpoint: (condition, mask_or_resolvent, retain_boolean, endpoint)
  same-state: ()
  resonance-blocks: (); concrete block membership is supplied at evaluation
  boolean-condition: (and/or/not, tuple_of_conditions)
  symmetric-condition: (condition,)
  energy-window: (tolerance, modes_or_null); frequencies remain runtime inputs
  integer-lattice: (modes, integer_basis_columns); only the mathematical lattice
    is stored here, not the derived evaluation pivot/inverse program
  polynomial-invariant: (functions, modes, freeze_other_modes, tolerance)
  polynomial-relationship: (function,modes,comparison,tolerance,freeze_other_modes)
  elementary-function: (operation, arguments)
    variable:(name,idx); scalar/power/exponent/logarithm:(parameter,idx)
    identity/exp/ln/sin/cos:(idx,)
    elementary-sum/product/composition:(functions,idx)
    sum/product/composition:(functions,indices)
    Original native function evaluation conventions are retained, including
    integer-input Scalar dtype behavior. Known operations only; no code strings.

Unknown custom conditions/functions and wider precision fail export explicitly.
Runtime caches and _free_index_plan_data are omitted from mathematical IR.

Compiled manifest and graph
  kind:"compiled-analytic-expression"; plan_version:1; complete:true
  program_versions:{polynomial:3,function:1,denominator:2,degeneracy:1}
  Readers also accept polynomial revisions 1/2 and denominator revision 1.
  source_ir_signature, semantic change/baseline/projection/specializations,
  parts and checksums, options.free_index_reduction, producer timing/provenance.
  Readers reject unsupported program versions. Source hashes are provenance.
  Each part format:"psience-compiled-expression", version:2 (also reads 1), with the same value
  tags/reference and buffer rules as mathematical IR.

  Root mapping:
    expression: compiled-tensor node
    execution: version, max_order, free-index key groups, multiplicity factorials,
      coefficient_patterns:{arity:(pattern_array,pattern_to_ID_map)}, degeneracy
    plans:{native/reduced:(subexpression_mapping,term_count,index_bindings)}
    free_index: prepared canonical values/bindings/statistics, or null

  Earlier polynomial-program:(instructions,ndim,order,signature,is_zero,occupation_offsets[,graph_program])
    Instruction:(opcode,arguments,extra,release_indices).
    constant: arguments is the numeric constant.
    leaf: arguments is ((weight,divisor,axes),...), axes is
      ((input_axis_or_minus_one,integer_offset,coefficient_array,coefficient_key),...).
      Evaluate each axis polynomial at occupation+offset; -1 uses zero occupation.
      Multiply axes, multiply weight, divide divisor, then sum leaf terms in order.
    alias/scale/mul/sum/linear: arguments are prior value indices. extra is scale
      or ordered weights where applicable. Sums retain native addition order.
    release_indices free values after their final use; output is the last value.
    Coefficient convolution/shifts and structural remaps are already resolved.
    Revision 2 keeps child program references and compiled coordinate bindings
    instead of unrolling every coordinate context; only old graph records use it.

  Current polynomial revision 3:
    table-polynomial:(table,node_index,ndim,order,signature,is_zero,occupation_offsets)
    polynomial-table:(opcodes,dimensions,child_offsets,children,parameter_offsets,
      parameters,terms,axes,coefficients,coefficient_keys,scalars,signature)
    Integer columns are mapped .npy buffers in portable mode (inline arrays
    in native pickle): opcodes int8, the rest int64.
    Node i uses children[child_offsets[i]:child_offsets[i+1]] and the analogous
      parameter slice. Child references are table row indices, possibly forward.
    Opcodes by ordinal: constant,leaf,sum,linear,scale,shift,remap,permutation_sum,
      project,mul. dimensions[i] is the node's formal input dimension.
    Parameters: constant/scale use one scalar ID; linear uses one per child;
      shift/remap store their integer vectors; permutation_sum stores one pull
      per child (width child dimension); project stores (child_dimension,keep...);
      mul concatenates left/right pulls; sum has none; leaf stores (first,last)
      selecting rows of terms. Terms rows are (weight_scalar,divisor_scalar,
      first_axis,last_axis). Axes rows are (input_column,offset,coefficient_ID).
    Coefficients are shared primitive arrays; scalars retain typed numeric values.
    The numeric call stack executes these saved calls/bindings directly, without
      building symbolic graphs, generated code or per-instruction Python objects.
    Array views retain the immutable mapping. Completed call frames are released;
      actual-state numerical caches remain bounded by the existing evaluator.

  function-program:(instructions,signature,elementary,idx)
    Instruction:(opcode,dependency_indices,parameter), output is the last value.
    input, take, identity, constant, sum, product, power, exponent, log, exp,
    ln, sin, cos; composition is lowered to dependencies and selections.
    This is a portable NumPy numerical program, not Python AST or native code.

  condition-program:(name,prepared_fields,flags,key)
    same-state/resonance-blocks:(); boolean:(operation,conditions);
    symmetric:(condition,); energy-window:(tolerance,modes);
    integer-lattice:(modes,basis,pivots,inverse_numerators,denominator);
    invariant:(function_programs,modes,freeze_other_modes,tolerance);
    relationship:(function_program,modes,comparison,tolerance,freeze_other_modes).
    flags:(symmetric,diagonal,equivalence,requires_frequencies).
    key is a stable structural fingerprint; no actual resonance membership.
  compiled-endpoint:(condition_program,role,retain,endpoint,literal_integer_tag)
  compiled-tensor / compiled-energy:(terms,prefactor,reduced,ndim,index_map_or_null)
  subexpression-plan:(compiled_expression,direct_polynomial,prefactor,groups)
  energy-term-plan:(ordinal,energy_changes,polynomial,contraction_shape,projected,
                    constant_numerator,denominator_program)
  energy-program: tuple of named denominator-program fields; a shared graph node
    rather than a repeated inline mapping for every numerator.
  compiled-energy-plan:(ordered_energy_term_references,prefactor,reduced,ndim)
    Reuses the exact stored term keys/numerators to reconstruct the mapping,
    avoiding a duplicate copy of every denominator key in the file.
    Denominator programs store ordinary rows or prepared projection/condition/
    domain instructions, baseline deltas, grouping, slots and projection indices.
  execution.degeneracy stores formal energy side changes, linear/perfect abstract
    requirements, compact projection maps and coefficient-key/energy-key inventory.
    Concrete physical resonance selection remains a runtime binding.

Compiled loading restores prepared fields directly; it must not invoke symbolic
canonicalization, factor convolution, free-index planning, predicate compilation,
request-pattern analysis, degree planning or abstract degeneracy construction.
Evaluation has no default persistent destination and writes only in a separately
configured evaluation_folder. No artifact contains molecular tensor/frequency/
occupation arrays, checkpoint data, or serialized Numba machine code.

Runner libraries and request metadata
  AnalyticVPTRunner.run_simple(..., ir_file=..., compiled_file=...,
    storage_backend='portable'|'pickle', evaluation_folder=...) orchestrates the
    same export/lowering/read interfaces. ir_file is optional; compiled_file
    is required when selecting these artifact options. Legacy expressions_file
    and expressions_folder keep their existing meanings.
  The libraries use the ordinary collection layout. Each entry label is a
    backend-independent fingerprint of request_version:1, the formal operator,
    change/baseline, projection and any explicit specialization. Numerical
    coefficients/frequencies/occupations/resonance membership are not inputs.
    Each entry manifest stores request_version:1 and request_signature:label.
    These fields survive independent compile_ir compilation of an IR library.
    Request metadata is a versioned lookup convention, not an exact producer
    source hash or a canonical algebraic-equivalence claim.
  A compatible compiled entry skips symbolic construction and all structural
    compilers. A retained matching IR entry skips symbolic construction and
    enters independent lowering. Otherwise, the existing symbolic builder
    supplies the source. New entries are staged during the public calculation;
    each completed library is published atomically on normal exit. Exceptions
    before publication leave the previous libraries intact.
  storage_backend=None preserves existing entries; new entries default portable.
    An explicit backend converts existing entries without recompilation.
    compilation_options.reuse_only=True rejects missing/incompatible requests,
    including a requested but absent mathematical source. Adding ir_file later
    builds missing source IR without recompiling existing prepared plans.
  File-backed execution consumes stored full-expression programs, for either
    direct/actions controller policy, without molecule-based symbolic pruning
    or runtime occupation-coefficient projection. Native/reduced saved-plan
    selection and actual numerical routing remain supported. Evaluation output
    is optional and must be outside the artifact folders.
```
