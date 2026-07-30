### `Analytic.py` — Provides a symbolic approach to vibrational perturbation theory based on a Harmonic description
  - **class `DefaultValues`** (enum.Enum)
  - **class `AnalyticPerturbationTheorySolver`**
    > A re-attempt at using the recursive expressions
    > to provide simpler code for getting APT expressions
    - `__init__(hamiltonian_expansion, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None)`
    - `from_order(order, internals=True, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None)`
    - `modify_hamiltonian(hamiltonian_corrections)`
    - `get_correction(key, cls, order, **kw)`
    - `shifted_hamiltonian_correction(order, **kw)`
    - `energy_correction(order, **kw)`
    - `wavefunction_correction(order, **kw)`
    - `overlap_correction(order, degenerate_changes=None, **kw)`
    - `full_wavefunction_correction(order, **kw)`
    - `operator_correction(order, operator_type=None, **kw)`
    - `operator_degenerate_correction(self, order, /, degenerate_changes, operator_type=None, **kw)`
    - `reexpressed_hamiltonian(order, **kw)`
    - `reexpressed_hamiltonian_degenerate_correction(self, order, /, degenerate_changes, **kw)`
    - `operator_expansion_terms(order, logger=None, base_index=None, operator_type=None)`
    - `clear_caches()`
  - **class `PolynomialInterface`**
    > Provides a basic interface to allow for the uniform manipulation
    > of objects that dispatch down to some form of scalar multiplied by
    > a sum of polynomials
    - `format_expr()`
    - `ndim()`
    - `audit(target, ignore_constants=True)`
    - `ensure_dimension(ndim)`
    - `align_dimensions(other)`
    - `shift(shift)`
    - `scale(scaling)`
    - `permute(perm, check_perm=True, allow_padding=False)`
    - `mul_simple(other)`
    - `rmul_simple(other)`
    - `mul_along(other, inds, remainder=None, mapping=None)`
    - `rmul_along(other, inds, remainder=None, mapping=None)`
    - `combine(**kwargs)`
    - `mutate(*args, **kwargs)`
  - **class `PolyPath`**
    > A simple holder class that contains the set of modifications along each dimension
    > used to build a given polynomial
    - `__init__(paths, scaling)`
    - `as_tuple()`
    - `mul_along(other, inds, remainder)`
  - **class `TreeSerializer`**
    - `serialize_iterable(iterable, primitive_test, concat, track_shapes=True)`
    - `build_dict_trees(dict_obj)`
    - `serialize_tree_dict(dict_obj, key_primitive_test=None, vals_primitive_test=None, concat=None)`
    - `default_concat(arrays)`
    - `default_prim(obj)`
    - `deserialize_subiterable(shape, flat, i, pad, depth)`
    - `deserialize_iterable(shape, flat, depth)`
    - `stitch_dict(iterables)`
    - `rebuild_tree_dict(key_trees, tree_shape, val_buffers)`
  - **class `ProductPTPolynomial`** (PolynomialInterface)
    > TODO: include prefactor term so we can divide out energy changes
    - `__init__(coeffs, prefactor=1, idx=None, steps=None)`
    - `lookup(idx)`
    - `mutate(coeffs=default, prefactor=default, idx=None, steps=default)`
    - `to_state(serializer=None)`
    - `ndim()`
    - `order()`
    - `audit(target, ignore_constants=True)`
    - `ensure_dimension(ndim)`
    - `pad(left_right_pads)`
    - `prep_float(c)`
    - `format_simple_poly(coeffs, idx)`
    - `format_expr()`
    - `permute(new_inds, check_perm=True, allow_padding=False)`
    - `constant_rescale()` — rescales so constant term is 1
    - `monic_coeffs()`
    - `monic_coeffs(coeffs)`
    - `combine(other)`
    - `condense(inds=None, return_inds=False, check_inds=True)`
    - `shift(shift)`
    - `scale(scalar)`
    - `mul_simple(other)`
    - `rmul_simple(other)`
    - `fast_ind_remainder(n, diff)`
    - `get_index_mapping(dim_1, dim_2, inds, return_remainder=True)` — Computes the corresponding and indices remainders for multdimensional polynomial multiplications
    - `mul_along(other, inds, remainder=None, mapping=None)`
    - `rmul_along(other, inds, remainder=None, mapping=None)` — We multiply every subpoly along the given indices, transposing the appropriate tensor indices
  - **class `ProductPTPolynomialSum`** (PolynomialInterface)
    - `__init__(polynomials, prefactor=1, ndim=None, order=None, reduced=False)`
    - `prep_serialization_dict()`
    - `from_serialization_dict(big_dict)`
    - `to_state(serializer=None)`
    - `mutate(polynomials=default, prefactor=default, ndim=default, order=default, reduced=default)`
    - `ndim()`
    - `order()`
    - `audit(target=None, ignore_constants=True)`
    - `ensure_dimension(ndim)`
    - `format_expr()`
    - `permute(new_inds, check_perm=True, allow_padding=False)`
    - `combine_polys(poly_set, cache)`
    - `combine()`
    - `condense(inds=None, return_inds=False, check_inds=True)`
    - `shift(shift)`
    - `scale(scalar)`
    - `mul_simple(other)`
    - `rmul_simple(other)`
    - `mul_along(other, inds, remainder=None, mapping=None)`
    - `rmul_along(other, inds, remainder=None, mapping=None)` — We multiply every subpoly along the given indices, transposing the appropriate tensor indices
  - **class `PTEnergyChangeProductSum`** (TensorCoefficientPoly, PolynomialInterface)
    > A representation of a sum of 1/energy * poly sums
    > which is here so we can transpose energy change indices intelligently
    - `__init__(terms, prefactor=1, canonicalize=True, reduced=False)`
    - `mutate(terms=default, prefactor=default, canonicalize=default, reduced=default)`
    - `flip_energy_terms()`
    - `prep_serialization_dict()`
    - `from_serialization_dict(big_dict)`
    - `to_state(serializer=None)`
    - `from_state(state, serializer=None)`
    - `filter(terms, mode='match')`
    - `canonical_key(monomial_tuple)`
    - `side_change_iter(key)`
    - `format_key(key)`
    - `format_energy_prod_key(key)`
    - `format_expr()` — Formats in a Mathematica-ingestible format
    - `shift_key(key, shift)`
    - `shift_energies(shift)`
    - `shift(shift, shift_energies=False)`
    - `scale(scaling)`
    - `get_key_ndim(terms)`
    - `ndim()`
    - `audit(target=None, ignore_constants=True)`
    - `ensure_dimension(ndim)`
    - `sort()`
    - `permute(new_inds, check_perm=True, allow_padding=False)`
    - `find_term_scaling(key)`
    - `combine_energies()`
    - `combine(combine_subterms=True, combine_energies=False)`
    - `mul_along(other, inds, remainder=None, mapping=None)` — We multiply every subpoly along the given indices, transposing the appropriate tensor indices
    - `rmul_along(other, inds, remainder=None, mapping=None)` — We multiply every subpoly along the given indices, transposing the appropriate tensor indices
    - `mul_simple(other)` — We multiply every subpoly along the given indices, transposing the appropriate tensor indices
    - `rmul_simple(other)`
  - **class `PTTensorCoeffProductSum`** (TensorCoefficientPoly, PolynomialInterface)
    > A representation for a sum of tensor coefficients * polynomial sums
    > which is primarily here so we can transpose tensor coefficients intelligently
    - `__init__(terms, prefactor=1, canonicalize=True, ndim=None, inds_map=None, reduced=False)`
    - `mutate(terms=default, *, prefactor=default, ndim=default, inds_map=default, canonicalize=default, reduced=default)`
    - `prep_serialization_dict()`
    - `from_serialization_dict(big_dict)`
    - `to_state(serializer=None)`
    - `from_state(state, serializer=None)`
    - `flip_energy_terms()`
    - `filter(terms, mode='match')`
    - `filter_coefficients(terms, mode='match')`
    - `filter_energies(terms, mode='match')`
    - `format_key(key)`
    - `format_tensor_key(key)`
    - `format_expr()` — Formats in a Mathematica-ingestible format
    - `coeff_product_inds(key, return_counts=False)`
    - `get_inds(key)`
    - `prune_operators(ops)`
    - `ndim()`
    - `audit(target=None, required_dimension=None, ignore_constants=True)` — Checks to ensure that the number of dimensions aligns with
    - `ensure_dimension(ndim)`
    - `ensure_subpoly_dim(terms)`
    - `sort()`
    - `permute(new_inds, check_perm=True, allow_padding=False)`
    - `combine_terms()`
    - `reindex_terms(terms)`
    - `combine(combine_coeffs=False, combine_subterms=True, combine_energies=False)`
    - `free_up_indices(start, stop)`
    - `shift_energies(change)`
    - `shift(shift)`
    - `scale(scaling)`
    - `symmetrizers()`
    - `canonical_key(monomial_tuple, symmetrizers=None)`
    - `mul_along(other, inds, remainder=None, index_classes=None, mapping=None, baseline=None)` — We multiply every subpoly along the given indices, transposing the appropriate tensor indices
    - `rmul_along(other, inds, remainder=None, mapping=None)` — We multiply every subpoly along the given indices, transposing the appropriate tensor indices
    - `mul_simple(other)` — We multiply every subpoly along the given indices, transposing the appropriate tensor indices
    - `rmul_simple(other)`
  - **class `SqrtChangePoly`** (PolynomialInterface)
    - `__init__(poly_obj, change, shift, canonicalize=False)`
    - `mutate(poly_obj=default, poly_change=default, shift_start=default, canonicalize=False)`
    - `to_state(serializer=None)`
    - `from_state(state, serializer=None)`
    - `strip()`
    - `ndim()`
    - `audit(target=None, ignore_constants=True)`
    - `ensure_dimension(ndim)`
    - `short_repr()`
    - `format_sqrt_expr()`
    - `format_expr()` — Formats in a Mathematica-ingestible format
    - `combine(combine_coeffs=False, combine_subterms=True, combine_energies=False)`
    - `sort()`
    - `canonicalize(poly_obj, poly_change, shift_start)`
    - `shift_energies(changes)`
    - `shift(shift, shift_energies=False)`
    - `scale(scaling)`
    - `permute(perm, check_perm=True, allow_padding=False)`
    - `filter_coefficients(terms, mode='match')`
    - `filter_energies(terms, mode='match')`
    - `get_change_poly(initial_changes, extra_changes, initial_shift, extra_shift)`
    - `mul_simple(other)`
    - `rmul_simple(other)`
    - `mul_along(other, inds, remainder=None, mapping=None, baseline=None)`
    - `rmul_along(other, inds, remainder=None, mapping=None)`
  - **class `PerturbationTheoryTerm`**
    > A generic version of one of the three terms in
    > PT that will generate a correction polynomial
    - `__init__(logger=None, checkpoint=None, allowed_terms=None, allowed_energy_changes=None, intermediate_normalization=None, allowed_coefficients=None, disallowed_coefficients=None)`
    - `get_subexpressions()`
    - `expressions()`
    - `change_sort_key(changes)`
    - `change_sort(changes)`
    - `sorted_changes(changes)`
    - `get_changes()`
    - `changes()`
    - `get_serializer_key()`
    - `serializer_key()`
    - `debug_logging()`
    - `default_logger()`
    - `default_filtering()`
    - `default_filters()`
    - `get_core_poly(changes, shift=None)`
    - `get_poly_terms(changes, simplify=None, shift=None)`
  - **class `OperatorExpansionTerm`** (PerturbationTheoryTerm)
    - `__init__(terms, order=None, identities=None, symmetrizers=None, index=None, allowed_terms=None, change_rules=None, **opts)`
    - `change_rules()`
    - `change_rules(rules)`
    - `get_changes()`
    - `get_subexpressions()`
    - `get_core_poly(changes, shift=None)`
  - **class `HamiltonianExpansionTerm`** (OperatorExpansionTerm)
    - `__init__(terms, order=None, identities=None, symmetrizers=None, change_rules=None, **opts)`
  - **class `PerturbationOperator`** (PerturbationTheoryTerm)
    - `__init__(subterm)`
    - `lookup(subterm)`
    - `get_changes()`
    - `get_core_poly(changes, shift=None)`
  - **class `ShiftedEnergyBaseline`** (PerturbationTheoryTerm)
    > Represents a term that will be multipled by on the left rather than the right
    > for evaluating things like Y[1]M[0]Y[1], essentially changing raising operations to lowering
    - `__init__(base_term)`
    - `lookup(base_term)`
    - `get_changes()`
    - `get_poly_terms(changes, shift=None, **opts)`
  - **class `ShiftedHamiltonianCorrection`** (PerturbationTheoryTerm)
    > Adds the wave function correction and the overlap term
    - `__init__(parent, order, allowed_terms=None, **opts)`
    - `get_serializer_key()`
    - `get_changes()`
    - `get_subexpressions()`
  - **class `WavefunctionCorrection`** (PerturbationTheoryTerm)
    - `__init__(parent, order, allowed_terms=None, **opts)`
    - `get_serializer_key()`
    - `get_changes()`
    - `get_subexpressions()`
  - **class `EnergyCorrection`** (PerturbationTheoryTerm)
    - `__init__(parent, order, allowed_terms=None, **opts)`
    - `get_serializer_key()`
    - `get_changes()`
    - `get_subexpressions()`
    - `get_poly_terms(changes, shift=None, **opts)`
  - **class `WavefunctionOverlapCorrection`** (PerturbationTheoryTerm)
    - `__init__(parent, order, allowed_terms=None, degenerate_changes=None, **opts)`
    - `get_serializer_key()`
    - `get_changes()`
    - `get_subexpressions()`
  - **class `FullWavefunctionCorrection`** (PerturbationTheoryTerm)
    > Adds the wave function correction and the overlap term
    - `__init__(parent, order, allowed_terms=None, **opts)`
    - `get_serializer_key()`
    - `get_changes()`
    - `get_subexpressions()`
  - **class `OperatorCorrection`** (PerturbationTheoryTerm)
    - `__init__(parent, order, operator_type=None, allowed_terms=None, wavefunction_generator=None, base_index=None, **opts)`
    - `get_type_key(operator_type)`
    - `get_serializer_key()`
    - `get_repr_key()`
    - `get_changes()`
    - `default_wavefunction_generator(o)`
    - `wavefunction_generator(o)`
    - `get_subexpressions(bra_wavefunction_generator=None, ket_wavefunction_generator=None, bounds=None, const_zeros=None)`
  - **class `OperatorDegenerateCorrection`** (OperatorCorrection)
    - `__init__(parent, order, degenerate_changes=None, operator_type=None, allowed_terms=None, **opts)`
    - `default_wavefunction_generator(o)`
    - **class `Left`** (OperatorCorrection)
      - `__init__(parent, order, real_parent)`
      - `get_subexpressions()`
    - **class `Right`** (OperatorCorrection)
      - `__init__(parent, order, real_parent)`
      - `get_subexpressions()`
    - **class `Both`** (OperatorCorrection)
      - `__init__(parent, order, real_parent)`
      - `get_subexpressions()`
    - `get_subexpressions(bra_wavefunction_generator=None, ket_wavefunction_generator=None)`
    - `get_right_degenerate_expressions(bra_wavefunction_generator=None, ket_wavefunction_generator=None)`
    - `get_left_degenerate_expressions(bra_wavefunction_generator=None, ket_wavefunction_generator=None)`
    - `get_both_degenerate_expressions(bra_wavefunction_generator=None, ket_wavefunction_generator=None)`
  - **class `DiagonalHamiltonian`** (OperatorExpansionTerm)
    - `__init__()`
    - `get_changes()`
  - **class `ReexpressedHamiltonian`** (OperatorCorrection)
    - `__init__(parent, order, allowed_terms=None, degenerate_changes=None, hamiltonian_corrections=None, **opts)`
    - `prep_expansion(base_expansion, corrs, *, logger, allowed_coefficients, disallowed_coefficients, base_index=5)`
    - `get_repr_key()`
  - **class `ReexpressedHamiltonianDegenerateCorrection`** (OperatorDegenerateCorrection)
    - `__init__(parent, order, allowed_terms=None, hamiltonian_corrections=None, **opts)`
    - `get_repr_key()`
  - **class `ScaledPerturbationTheoryTerm`** (PerturbationTheoryTerm)
    - `__init__(base_term, scaling)`
    - `lookup(base_term, scaling)`
    - `get_changes()`
    - `get_core_poly(changes, shift=None)`
  - **class `PerturbationTheoryTermSum`** (PerturbationTheoryTerm)
    - `__init__(*terms)`
    - `get_changes()`
    - `get_subexpressions()`
  - **class `PerturbationTheoryTermProduct`** (PerturbationTheoryTerm)
    - `__init__(post_op, pre_op)`
    - `lookup(post_op, pre_op)`
    - `get_change_classes(changes)`
    - `get_combination_inds(n, r)`
    - `get_combination_comp(n, r)`
    - `get_changes()`
    - `get_expressions()`
    - `get_poly_product_terms(gen1, gen2, change_1, change_2, target_inds, remainder_inds, reorgs, simplify=True)`
    - `get_core_poly(changes, shift=None)`
  - **class `PerturbationTheoryExpressionEvaluator`**
    - `__init__(op, expr, change=None, logger=None, parallelizer=None)`
    - `set_cache_size(max_size)`
    - `get_cache()`
    - `evaluate_polynomial_expression(state_perms, coeffs, freqs, expr, change, baseline_shift, num_fixed, op=None, logger=None, parallelizer=None, degenerate_changes=None, only_degenerate_terms=False, zero_cutoff=None, verbose=False, log_scaled=True)`
    - `evaluate(state_perms, coeffs, freqs, degenerate_changes=None, only_degenerate_terms=False, zero_cutoff=None, parallelizer=None, verbose=False, log_scaled=True)`
  - **class `PerturbationTheoryEvaluator`**
    - `__init__(solver, expansion, freqs=None, parallelizer=None)`
    - `modify_hamiltonian(hamiltonian_corrections)`
    - `get_energy_corrections(states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False, logger=None, parallelizer=None)`
    - `is_single_expansion(expansion, min_order=None)`
    - `get_overlap_corrections(states, order=None, expansions=None, degenerate_states=None, freqs=None, zero_cutoff=None, verbose=False, parallelizer=None)`
    - `get_diff_map(state_map)`
    - `get_finals(initial, change, perms)`
    - `get_degenerate_changes(degenerate_pairs)`
    - `get_state_by_state_corrections(generator, states, order=None, terms=None, epaths=None, expansions=None, freqs=None, verbose=False, allowed_coefficients=None, disallowed_coefficients=None, degenerate_states=None, only_degenerate_terms=False, degenerate_correction_generator=None, include_degenerate_correction_terms=True, log_scaled=False, zero_cutoff=None, return_sorted=False, logger=None, parallelizer=None)`
    - `get_matrix_corrections(states, order=None, expansions=None, freqs=None, zero_cutoff=None, verbose=False)`
    - `get_full_wavefunction_corrections(states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False)`
    - `get_wavefunction_corrections(states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False)`
    - `get_reexpressed_hamiltonian(states, order=None, expansions=None, freqs=None, degenerate_states=None, only_degenerate_terms=False, verbose=False, include_diagonal=False, hamiltonian_corrections=None, **opts)` — :param hamiltonian_corrections:  `[[(order, terms), expansion], ...]`
    - `get_operator_corrections(operator_expansion, states, order=None, expansions=None, freqs=None, degenerate_states=None, operator_type=None, check_single=True, terms=None, min_order=1, verbose=False, **opts)`
    - `evaluate_expressions(states, exprs, expansions=None, operator_expansions=None, degenerate_states=None, zero_cutoff=None, verbose=False)`

### `Analyzer.py` — Provides analyzer class to handle common VPT analyses
  - **class `VPTResultsSource`** (enum.Enum)
    > Enum of sources to load PT results from
- `property_dispatcher(basefn)` — Class decorator factory that turns a method into a source-dependent "virtual property": the decorat…
  - **class `VPTAnalyzerLogParser`** (LogParser)
    - `__init__(log_file, **opts)`
    - `tree()` — The (cached) parsed block-tree structure of the log file, collapsed down to just the "Computing PT…
    - **class `EnergiesBlockParser`** (StringLineByLineReader)
      - `__init__(spec_str, **opts)`
      - `check_tag(line, depth=0, active_tag=None, label=None, history=None)` — Skip blank lines, the table header row (starting with `"State"`), and deeply-indented continuation…
      - `handle_block_line(label, line, depth=0, history=None)` — Split a data row of the energies table into the state label and its numeric columns, inferring the…
    - `reformat_eng_block(sb)` — Convert a parsed energies block (a list of `(state_label, value_tokens)` rows) into parallel arrays…
    - `parse_energies_blocks(spec_str)` — Parse a raw energies-table block of log text into `(states, harmonic, anharmonic)` arrays, via `Ene…
    - `harmonic_energies()` — The `(states, harmonic_energies)` pair parsed from the log's "States Energies" (or, for a degenerat…
    - `energies()` — The `(states, anharmonic_energies)` pair parsed from the log's "States Energies"/"Degenerate Energi…
    - `zero_order_energies()` — The `(states, harmonic_energies)` pair parsed from the log's "States Energies"/"Degenerate Energies…
    - `deperturbed_energies()` — The `(states, deperturbed_energies)` pair parsed from the log's "Deperturbed Energies" table block…
    - `spectra()` — The (cached) IR spectrum data parsed from the log's "IR Data" block, via `parse_spectrum_blocks`.
    - `deperturbed_spectra()` — The (cached) deperturbed IR spectrum data parsed from the log's "Deperturbed IR Data" block, via `p…
    - **class `SpectrumBlockParser`** (StringLineByLineReader)
      - `__init__(spec_str, **opts)`
      - `check_tag(line, depth=0, active_tag=None, label=None, history=None)` — Skip blank lines, header rows, and deeply-indented continuation lines, and recognize `" Initial Sta…
      - `handle_block_line(label, line, depth=0, history=None)` — Split a data row of a spectrum sub-block into the final-state label and its numeric columns, using…
    - `reformat_spec_block(sb)` — Convert a parsed spectrum block (per-initial-state groups of `(state_label, value_tokens)` rows) in…
    - `parse_spectrum_blocks(spec_str)` — Parse a raw multi-initial-state spectrum-table block of log text into a list of reformatted per-blo…
    - **class `TransitionMomentBlockParser`** (StringLineByLineReader)
      - `__init__(spec_str, **opts)`
      - `check_tag(line, depth=0, active_tag=None, label=None, history=None)` — Skip blank lines, header rows, and deeply-indented continuation lines, and recognize `" Initial Sta…
      - `handle_block_line(label, line, depth=0, history=None)` — Split a data row of a transition-moment sub-block into the final-state label and its numeric correc…
    - `load_term_counts(nterms)` — Compute how many perturbative-order transition-moment correction terms exist below a given total co…
    - `reformat_tm_block(sb)` — Convert a parsed transition-moment-correction block (per-initial-state groups of `(state_label, val…
    - `parse_tm_blocks(spec_str)` — Parse a raw multi-initial-state transition-moment-table block of log text into a list of reformatte…
    - `transition_moment_corrections()` — The (cached) transition-dipole-moment corrections parsed from the log's per-axis "X/Y/Z Dipole Cont…
    - `deperturbed_transition_moment_corrections()` — The (cached) deperturbed transition-dipole-moment corrections parsed from the log's per-axis "X/Y/Z…
  - **class `VPTResultsLoader`**
    > Provides tools for loading results into canonical
    > sources from a simulation, both from checkpoint files and from
    > `PerturbationTheoryWavefunctions` objects and potentially more
    - `__init__(res, res_type=None)`
    - `resolve_file_res_type(res)` — Infer whether a given file path refers to a checkpoint file or a text log file, based on its file e…
    - `get_res_type(res)` — :param res:
    - `potential_terms()` — Returns the expansion of the potential
    - `kinetic_terms()` — Returns the expansion of the kinetic energy
    - `dipole_terms()` — Returns the expansion of the dipole moment
    - `basis()` — Returns the basis for the calculation
    - `target_states()` — Returns the target states for the calculation
    - `spectrum()` — Returns the IR spectrum calculated from perturbation theory
    - `zero_order_spectrum()` — Returns the zero-order IR spectrum calculated from perturbation theory
    - `energy_corrections()` — Returns the corrections to the energies
    - `energies()` — :return:
    - `zero_order_energies()` — :return:
    - `deperturbed_energies()` — :return:
    - `wavefunction_corrections()` — Returns the corrections to the wavefunctions
    - `transition_moment_corrections()` — Returns the corrections to the wavefunctions
    - `transition_moments()` — Sum the transition-dipole-moment corrections (over all perturbative orders) into a single `(nstates…
    - `deperturbed_transition_moment_corrections()` — Returns the corrections to the wavefunctions
    - `deperturbed_transition_moments()` — :return:
    - `degenerate_states()` — Returns the deperturbed states used to make the degenerate transform
    - `deperturbed_hamiltonians()` — Returns the deperturbed Hamiltonians used to make the degenerate transform
    - `degenerate_energies()` — Returns the deperturbed states used to make the degenerate transform
    - `degenerate_rotations()` — Returns the deperturbed states used to make the degenerate transform
    - `log_file()` — Returns the log_file for the run
    - `log_parser()` — A `VPTAnalyzerLogParser` for this loader's data: the data itself, if it already is one, otherwise a…
- `loaded_prop(fn)` — Decorator for `VPTAnalyzer` convenience properties that mirror a `VPTResultsLoader` method of the s…
  - **class `VPTAnalyzer`**
    > Provides analysis tools on VPT results
    - `__init__(res)`
    - `run_VPT(*args, logger=None, **kwargs)` — Runs a VPT calculation through `VPTRunner.run_simple` and
    - `potential_terms()` — :return:
    - `kinetic_terms()` — :return:
    - `dipole_terms()` — :return:
    - `basis()` — :return:
    - `target_states()` — :return:
    - `spectrum()` — :return:
    - `energy_corrections()` — :return:
    - `energies()` — :return:
    - `frequencies()` — :return:
    - `zero_order_spectrum()` — :return:
    - `deperturbed_spectrum()` — :return:
    - `deperturbed_frequencies()` — :return:
    - `wavefunction_corrections()` — :return:
    - `transition_moment_corrections()` — :return:
    - `transition_moments()` — :return:
    - `deperturbed_transition_moment_corrections()` — :return:
    - `deperturbed_transition_moments()` — :return:
    - `deperturbed_hamiltonians()` — :return:
    - `zero_order_energies()` — :return:
    - `deperturbed_energies()` — :return:
    - `degenerate_states()` — :return:
    - `degenerate_energies()` — :return:
    - `shift_and_transform_hamiltonian(hams, shifts)` — :param hams:
    - `get_shifted_transformed_transition_moments(deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose')` — :param deg_states:
    - `get_shifted_transformed_spectrum(zpe, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose')` — :param zpe:
    - `shifted_transformed_spectrum(deg_states, hams, shifts, return_transformation=False, handling_mode='transpose')` — :param deg_states:
    - `transition_data(states, keys=['frequency', 'transition_moment', 'intensity'], data='deperturbed')` — :param states:
    - `transition_moment_term_sums(states, terms=None, rotation=None, data='deperturbed')` — :param states:
    - `transition_moment_term_sums_first_order(states, rotation=None, data='deperturbed')` — :param states:
    - `intensity_breakdown(states, terms=None, data='deperturbed')` — :param states:
    - `degenerate_coupling_element(state1, state2)` — :param state1:
    - `format_deperturbed_hamiltonian(which)` — :param which:
    - `log_parser()` — The underlying `VPTAnalyzerLogParser` for this analyzer's loaded results, via `self.loader.log_pars…
    - `print_output_tables(print_energy_corrections=False, print_energies=False, print_transition_moments=False, print_intensities=True, **kwargs)` — Print the standard VPT output tables (energy corrections, energies, transition moments, intensities…

### `Common.py`
  - **class `PerturbationTheoryException`** (Exception)
  - **class `Settings`**

### `Corrections.py`
  - **class `PerturbationTheoryCorrections`**
    > Represents a set of corrections from perturbation theory.
    > Can be used to correct other operators in the basis of the original calculation.
    - `__init__(states, coupled_states, total_basis, energy_corrs, wfn_corrections, all_energy_corrections=None, degenerate_states=None, degenerate_transformation=None, degenerate_energies=None, degenerate_hamiltonians=None, nondeg_hamiltonian_precision=3, logger=None)`
    - `from_dicts(states, corrections, **opts)` — :param states: a dict with the states described by the corrections, the set of states coupled, and…
    - `degenerate()` — :return:
    - `energies()` — :return:
    - `order()` — :return:
    - `take_subspace(space)` — Takes only those elements that are in space
    - `create_coupling_matrix(corrs, states, flat_total_space, nstates, order, filters=None, non_zero_cutoff=1e-14, logger=None)` — :param corrs:
    - `prune(threshold=0.1, in_place=False)` — Returns corrections with couplings less than the given cutoff set to zero
    - `get_transformed_Hamiltonians(hams, deg_group=None)` — :param corrs:
    - `get_degenerate_rotation(deg_group, hams, label=None, zero_point_energy=None, local_coupling_hamiltonian=None, local_coupling_order=None)` — :param deg_group:
    - `get_degenerate_transformation(group, hams, gaussian_resonance_handling=False, label=None, zero_point_energy=None, local_coupling_hamiltonian=None, local_coupling_order=None)` — Compute the degenerate-perturbation-theory rotation for a single group of resonant/degenerate state…
    - `default_state_filter(state, couplings, energy_cutoff=None, energies=None, basis=None, target_modes=None)` — Excludes modes that differ in only one position, prioritizing states with fewer numbers of quanta
    - `find_strong_couplings(threshold=0.1, state_filter=None)` — Finds positions in the expansion matrices where the couplings are too large
    - `format_strong_couplings_report(couplings=None, threshold=0.1, int_fmt='{:>3.0f}', padding='{:<8}', join=True, use_excitations=True)` — Format a human-readable report of the states found by `find_strong_couplings` (or an explicitly sup…
    - `collapse_strong_couplings(sc)` — :param sc:
    - `operator_representation(operator_expansion, order=None, subspace=None, contract=True, logger_symbol='A', logger_conversion=None)` — Generates the representation of the operator in the basis of stored states
    - `get_overlap_matrices()` — Returns the overlap matrices for the set of corrections
    - `savez(file)` — Intended to serialize the corrections to an `npz` file, but currently disabled -- immediately raise…
    - `loadz(file)` — Intended to reconstruct corrections from an `npz` file previously written by `savez`, but currently…
    - `to_state(serializer=None)` — Serialize this object's core data (states, coupled states, total basis, energy/wavefunction correct…
    - `from_state(data, serializer=None)` — Reconstruct a `PerturbationTheoryCorrections` object from a previously serialized state dict, deser…
  - **class `AnalyticPerturbationTheoryCorrections`**
    - `get_zpe_pos()` — Find (and cache) the index of the zero-point-energy (ground) state within `self.states`, falling ba…
    - `energies()` — The final state energies: the sum of the (per-order) energy corrections if there are no degenerate…
    - `deperturbed_energies()` — The (cached) deperturbed state energies -- the sum of the per-order energy corrections, without any…
    - `handle_degenerate_transformation(degenerate_ham)` — Diagonalize a degenerate-block Hamiltonian and reorder its eigenvalues/eigenvectors so that each ou…
    - `get_degenerate_transformations(basis, energies)` — Apply degenerate perturbation theory block by block: for each group of degenerate states, build its…
    - `degenerate_hamiltonians()` — The (cached) per-block degenerate Hamiltonians, computed as a side effect of evaluating `energies`…
    - `degenerate_coefficients()` — The (cached) per-block degenerate-perturbation-theory mixing coefficients, computed as a side effec…
    - `get_freqs()` — Compute the vibrational transition frequencies (final energies) relative to the zero-point energy.
    - `get_deperturbed_freqs()` — Compute the deperturbed vibrational transition frequencies relative to the deperturbed zero-point e…
    - `degenerate_transformation_pairs()` — The (cached) per-block `(row_tf, col_tf)` transformation pairs mapping each `state_lists` block's r…
    - `transition_moments()` — The (cached) final transition-dipole moments: the deperturbed transition moments if there are no de…
    - `harmonic_transition_moments()` — The purely harmonic (zeroth-order) transition-dipole moments, extracted from the first entry of eac…
    - `deperturbed_transition_moments()` — The (cached) deperturbed transition-dipole moments, summing each block's transition-moment correcti…
    - `get_spectra(energies, transition_moments)` — Build a per-block list of discrete IR spectra from a set of state energies and transition moments,…
    - `harmonic_spectra()` — The purely harmonic (zeroth-order) IR spectra, built from the zeroth-order energy corrections and h…
    - `deperturbed_spectra()` — The (cached) deperturbed IR spectra, built from the deperturbed energies and transition moments via…
    - `spectra()` — The (cached) final IR spectra: the deperturbed spectra if there are no degenerate states, otherwise…
    - `deperturbed_operator_values()` — The (cached) deperturbed values of any extra tracked operators, summing each operator's per-block c…
    - `operator_values()` — The (cached) final values of any extra tracked operators: the deperturbed operator values if there…

### `DegeneracySpecs.py`
  - **class `DegenerateSpaceInputFormat`** (enum.Enum)
    > Real access pattern: DegenerateSpaceInputFormat.<MemberName> (this is an enum with 8 members, e.g. DegenerateSpaceInputFormat.Groups == 'groups'). Collapsed into a dict below purely for compactness -- do not index it as a dict in real code:
  - **class `DegeneracySpec`**
    > Provides a container for specifying degeneracies
    > in a way that can be cleanly canonicalized
    - `__init__(application_order=None, group_filter=None, energy_cutoff=0.00225, test_modes=None, max_mode_differences=None, maximize_filtered_groups=True, decoupling_overide=100, extra_groups=None, uncoupled_states=0, inconsistent_polyads=None, wavefunction_corrections=None)`
    - `merge_state_blocks(state_blocks)` — Merge a collection of (possibly overlapping) groups of states into a set of maximal connected group…
    - `get_degenerate_group_filter(solver=None, evaluator=None, corrs=None, frequencies=None, zero_order_energies=None, high_frequency_modes=None, logger=None, low_frequency_mode_cutoff=0.00115, threshold=None)` — Resolve the group-filter specification to actually use when identifying degenerate groups, building…
    - `get_format_mapping()` — The mapping from each `DegenerateSpaceInputFormat` enum value to the concrete `DegeneracySpec` subc…
    - `get_default_spec()` — The default degeneracy-handling spec to use when none is given explicitly (a `StronglyCoupledDegene…
    - `from_spec(spec, format=None, **kwargs)` — Build a concrete `DegeneracySpec` instance from a raw specification, inferring which format it repr…
    - `infer_format(spec)` — Infer which `DegenerateSpaceInputFormat` a raw specification represents: first checking for an expl…
    - `get_groups(input_states, solver=None, **kwargs)` — :param solver:
    - `get_group_polyad_relation(exc)` — Compute the pairwise "polyad" relation vectors for a set of excitation vectors: for every pair of s…
    - `get_polyad_pairs_from_polyad_specs(polyads)` — Flatten a list of polyad groups (each a list of transformation rules) into a flat list of all pairw…
    - `get_polyad_pairs(input_states=None, groups=None, solver=None, **kwargs)` — :param solver:
    - `canonicalize(spec)` — Abstract hook for validating/normalizing a raw specification into the form this subclass expects, u…
    - `get_group_filter(target_modes=None, max_mode_differences=None, maximize_groups=None, decoupling_overide=None, extra_groups=None, uncoupled_states=None, corrections=None, **kwargs)` — Build the group-filter callable used to prune/validate the degenerate groups this spec identifies,…
  - **class `EnergyCutoffDegeneracySpec`** (DegeneracySpec)
    - `__init__(cutoff, **opts)`
    - `canonicalize(spec)` — Check whether a raw specification is a bare numeric energy cutoff, matching this subclass's expecte…
  - **class `MartinTestDegeneracySpec`** (DegeneracySpec)
    - `__init__(threshold=4.6e-06, test_energy_window=0.0046, convert=True, frequencies=None, **opts)`
    - `prep_states(states)` — Identify, for each input state, the nearby (within `self.window`) states reachable via a cubic (`x^…
    - `get_coupled_spaces(input_states, solver=None)` — Evaluate the Martin resonance test for each input state against its candidate coupling partners (co…
    - `get_groups(input_states, solver=None, extra_groups=None, **kwargs)` — Build the final degenerate groups by running the Martin test (via `get_coupled_spaces`) and merging…
    - `canonicalize(spec)` — Check whether a raw specification looks like a `MartinTestDegeneracySpec`, by checking that it has…
    - `get_group_filter(**kwargs)` — Build the group-filter callable for this spec, supplying the Martin-test-specific state/basis/corre…
  - **class `StronglyCoupledDegeneracySpec`** (DegeneracySpec)
    - `__init__(wfc_threshold=None, state_filter=None, extend_spaces=True, iterations=None, evaluator=None, **opts)`
    - `application_order()` — Property getter/setter for when this spec's degeneracy handling is applied: `'pre'` if an `evaluato…
    - `application_order(ord)` — Property getter/setter for when this spec's degeneracy handling is applied: `'pre'` if an `evaluato…
    - `prep_states(input_states)` — Return the input states unchanged; this spec doesn't need to expand the input state space ahead of…
    - `identify_strong_couplings(solver, corrs)` — Find the strongly-coupled states for every state in `corrs`, via `PerturbationTheoryCorrections.fin…
    - `get_input_state_couplings(input_states)` — Compute (and cache) the strongly-coupled partner states for a set of input states, using `self.eval…
    - `get_groups(input_states, couplings=None, solver=None, extra_groups=None, **kwargs)` — :param input_states:
    - `canonicalize(spec)` — Check whether a raw specification looks like a `StronglyCoupledDegeneracySpec`: a bare number (the…
    - `get_strong_coupling_space(states, couplings, extra_groups=None)` — Build the final degenerate groups from a mapping of state-to-coupled-states, merging each state wit…
    - `get_degenerate_group_filter(solver=None, evaluator=None, threshold=None, logger=None, frequencies=None, **kwargs)` — Build the group-filter specification for this spec, supplying the coupling threshold and an appropr…
  - **class `GroupsDegeneracySpec`** (DegeneracySpec)
    - `__init__(groups, **opts)`
    - `canonicalize(spec)` — Check whether a raw specification is a list of valid groups, matching this subclass's expected inpu…
    - `get_groups(input_states, solver=None, **kwargs)` — :param solver:
  - **class `PolyadDegeneracySpec`** (DegeneracySpec)
    - `__init__(polyads, max_quanta=None, iterations=2, require_converged=False, extra_groups=None, extra_polyads=None, full_group_polyads=True, **opts)`
    - `canonicalize(spec)` — Check whether a raw specification is a list of valid polyad rules, matching this subclass's expecte…
    - `get_groups(input_states, solver=None, **kwargs)` — :param solver:
    - `get_degenerate_polyad_space(states, polyadic_pairs, max_quanta=None, iterations=2, require_converged=False, extra_groups=None)` — Gets degenerate spaces by using pairs of transformation rules to
    - `get_polyad_pairs(input_states=None, groups=None, solver=None, full_group_polyads=None, **kwargs)` — Get the pairwise polyad transformation rules used for identifying strong couplings: either derived…
  - **class `TotalQuantaDegeneracySpec`** (PolyadDegeneracySpec)
    - `__init__(n_T_vectors, max_quanta=3, target_modes=None, extra_polyads=None, **opts)`
    - `reduce_rule(a, b)` — Reduce a pair of raising/lowering quantum-number patterns to their simplest equivalent form: dividi…
    - `make_nt_polyad(nt, target_modes=None, max_quanta=3)` — Build the set of pairwise polyad transformation rules connecting states that share the same weighte…
  - **class `CallableDegeneracySpec`** (DegeneracySpec)
    - `__init__(callable, **opts)`
    - `get_groups(input_states, solver=None, **kwargs)` — :param input_states:
  - **class `DegenerateMultiStateSpace`** (BasisMultiStateSpace)
    - `default_group_filter(group, states=None, basis=None, corrections=None, threshold=None, energy_cutoff=None, frequencies=None, energies=None, decoupling_overide=10, maximize_groups=True, target_modes=None, max_mode_differences=None, threshold_step_size=0.1, extra_groups=None, uncoupled_states=None, state_diff_filters=1, logger=None)` — Excludes modes that differ in only one position, prioritizing states with fewer numbers of quanta
    - `construct_filer(spec=None, **kwargs)` — Build the group-filter callable to use when identifying degenerate groups: if a `DegeneracySpec` is…
    - `from_spec(degenerate_states, states=None, solver=None, evaluator=None, full_basis=None, format=None, group_filter=None, log_groups=False, **kwargs)` — Generates a DegenerateMultiStateSpace object from a number

### `Hamiltonian.py` — Provides support for build perturbation theory Hamiltonians
  - **class `PerturbationTheoryHamiltonian`**
    > Represents the main Hamiltonian used in the perturbation theory calculation.
    > Uses a harmonic oscillator basis for representing H0, H1, H2 etc.
    > The PT process is split into a PerturbationTheorySolver and a PerturbationTheoryHamiltonian
    > where the Hamiltonian just implements the split of the perturbation and the Solver manages the equations.
    - `__init__(molecule=None, n_quanta=None, modes=None, mode_selection=None, mode_transformation=None, rephase_modes=None, local_mode_couplings=False, local_mode_coupling_order=None, full_surface_mode_selection=None, potential_derivatives=None, include_potential=True, include_gmatrix=True, include_coriolis_coupling=True, include_pseudopotential=True, include_only_mode_couplings=None, potential_terms=None, allow_higher_potential_terms=False, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, include_dipole=True, dipole_terms=None, selection_rules=None, operator_chunk_size=None, operator_coefficient_threshold=None, matrix_element_threshold=None, logger=None, checkpoint=None, results=None, parallelizer=None, **expansion_options)`
    - `from_fchk(file, internals=None, mode_selection=None, **kw)` — :param file: fchk file to load from
    - **class `TermGetter`**
      - `__init__(base_terms, input_terms, mode_selection=None)`
      - `take_mode_subset(V, sel)` — Restrict a term tensor to a subset of modes along every axis, if a mode selection is configured.
      - `adjust_base_term(V)` — Hook for post-processing a base term before it's returned; on the base `TermGetter` this is a no-op…
    - **class `CoriolisTermGetter`** (TermGetter)
      - `adjust_base_term(Z)` — Collapse the Coriolis term's 3 rotational-axis components down to a single tensor by summing the `x…
    - `dipole_terms()` — The (lazily constructed and cached) `DipoleTerms` object used to expand the dipole surface, or `Non…
    - `prep_local_couplings(local_mode_couplings)` — Normalize the `local_mode_couplings` constructor argument into a `[v0, g0]` pair of local-mode pote…
    - `prep_operator_terms(coeffs, order)` — Build the perturbative expansion terms for an arbitrary operator given as a `[constant, deriv1, der…
    - `get_perturbations(expansion_orders, return_reps=True, order=None)` — Gets the `Representation` objects for the perturbations up through second order
    - `get_Nielsen_xmatrix(freqs=None, v3=None, v4=None, zeta_Be=None)` — Provides Nielsen's X-Matrix when working in Cartesian coordinates
    - `get_Nielsen_energies(states, x_mat=None, freqs=None, v3=None, v4=None, zeta_Be=None, return_split=False, return_X=False)` — :param states:
    - `get_2nd_order_freqs(states, *, freqs=None, V_terms=None, G_terms=None)` — :param states:
    - `get_solver(states, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, **opts)` — Build a `PerturbationTheorySolver` for the given target states: resolves the per-term expansion ord…
    - `get_wavefunctions(states, initial_states=None, degeneracies=None, allow_post_PT_calc=True, ignore_odd_order_energies=True, use_full_basis=True, order=2, expansion_order=None, memory_constrained=None, target_property_rules=None, results=None, degenerate_transformation_layout=None, return_solver=False, **opts)` — Gets a set of `PerturbationTheoryWavefunctions` from the perturbations defined by the Hamiltonian
    - `get_action_expansion(coupled_states=None, degeneracies=None, allow_sakurai_degs=False, allow_post_PT_calc=True, modify_degenerate_perturbations=False, intermediate_normalization=False, ignore_odd_order_energies=True, zero_element_warning=True, state_space_iterations=None, order=2)` — Gets the expansion of the energies in terms of Miller's "classical actions" by
    - `get_breakdown(states, coupled_states=None, degeneracies=None, order=2)` — Intended to compute a term-by-term breakdown of the VPT energies (harmonic-only, +cubic, +quartic,…

### `Runner.py` — A little package of utilities for setting up/running VPT jobs
  - **class `VPTSystem`**
    > Provides a little helper for setting up the input
    > system for a VPT job
    > *(truncated — see stub for full docstring)*
    - `__init__(mol, internals=None, dummy_atoms=None, modes=None, local_modes=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, potential_derivatives=None, potential_function=None, order=2, dipole_derivatives=None, eckart_embed=False, copy_mol=False)`
    - `prep_local_modes(dRdX, dXdR=None, sort_freqs=False)` — Build a set of "local mode" normal-mode data (frequencies, mode matrix, and its inverse) from a set…
    - `nmodes()` — Provides the number of modes in the system
    - `get_potential_derivatives(potential_function, order=2, **fd_opts)` — Computes potential derivatives for the given function through finite difference
    - `from_harmonic_scan(scan_array)` — Intended to build a `VPTSystem` from a harmonic potential-energy scan array.
  - **class `VPTStateSpace`**
    > Provides a helper to make it easier to set up the input
    > state spaces/degenerate spaces to run the perturbation theory
    > *(truncated — see stub for full docstring)*
    - `__init__(states, degeneracy_specs=None, system=None, frequencies=None, evaluator=None)`
    - `from_system_and_spec(system, spec, **opts)` — Build a `VPTStateSpace` from a system and a flexible state specification: passes an already-built `…
    - `from_system_and_quanta(system, quanta, target_modes=None, only_target_modes=False, **opts)` — Takes a system and a number of quanta and constructs a state space
    - `get_state_list_from_quanta(n_quanta, n_modes, target_modes=None, only_target_modes=False)` — Gets states up to `n_quanta` over `n_modes`
    - `build_degenerate_state_spaces(degeneracy_specs, states, system=None, evaluator=None, freqs=None)` — :param degeneracy_specs:
    - `filter_generator(target_property, order=2, initial_states=None, postfilters=None)` — Build a callable that produces a state-space filter for a given target property, by binding `target…
    - `get_filter(target_property, order=2, initial_states=None, postfilters=None)` — Obtains a state space filter for the given target property
    - `get_state_space_filter(states, initial_states=None, n_modes=None, order=2, target='wavefunctions', postfilters=None, **opts)` — Gets `state_space_filters` for the input `states` targeting some property
  - **class `VPTHamiltonianOptions`**
    > Provides a helper to keep track of the levers available for
    > setting up the Hamiltonian
    - `__init__(mode_selection=None, mode_transformation=None, local_mode_couplings=None, local_mode_coupling_order=None, full_surface_mode_selection=None, include_potential=None, include_gmatrix=None, include_coriolis_coupling=None, include_pseudopotential=None, include_only_mode_couplings=None, potential_terms=None, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, include_dipole=None, dipole_terms=None, dipole_derivatives=None, undimensionalize_normal_modes=None, use_numerical_jacobians=None, eckart_embed_derivatives=None, eckart_embed_planar_ref_tolerance=None, strip_dummy_atoms=None, strip_embedding_coordinates=None, mixed_derivative_handling_mode=None, mixed_derivative_warning_threshold=None, mixed_derivative_handle_zeros=None, rephase_modes=None, backpropagate_internals=None, direct_propagate_cartesians=None, zero_mass_term=None, internal_fd_mesh_spacing=None, internal_fd_stencil=None, cartesian_fd_mesh_spacing=None, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=None, cartesian_by_internal_derivative_method=None, internal_by_cartesian_order=None, cartesian_by_internal_order=None, expansion_handling_mode=None, jacobian_warning_threshold=None, check_input_force_constants=None, hessian_tolerance=None, grad_tolerance=None, freq_tolerance=None, g_derivative_threshold=None, gmatrix_tolerance=None, use_internal_modes=None, use_cartesian_kinetic_energy=None, operator_coefficient_threshold=None, imaginary_frequency_handling_mode=None)`
  - **class `VPTRuntimeOptions`**
    > Provides a helper to keep track of the options available
    > for configuring the way the code runs
    - `__init__(operator_chunk_size=None, matrix_element_threshold=None, nondeg_hamiltonian_precision=None, logger=None, verbose=None, checkpoint=None, results=None, parallelizer=None, memory_constrained=None, checkpoint_keys=None, use_cached_representations=None, use_cached_basis=None)`
  - **class `VPTSolverOptions`**
    > Provides a helper to keep track of the options available
    > for configuring the way the perturbation theory is applied
    > *(truncated — see stub for full docstring)*
    - `__init__(order=2, expansion_order=None, coupled_states=None, total_space=None, flat_total_space=None, state_space_iterations=None, state_space_terms=None, state_space_filters=None, extended_state_space_filter_generator=None, extended_state_space_postprocessor=None, allow_post_PT_calc=None, modify_degenerate_perturbations=None, gaussian_resonance_handling=None, ignore_odd_order_energies=None, intermediate_normalization=None, zero_element_warning=None, degenerate_states=None, handle_strong_couplings=None, strong_coupling_test_modes=None, strong_couplings_state_filter=None, strongly_coupled_group_filter=None, extend_strong_coupling_spaces=None, strong_coupling_zero_order_energy_cutoff=None, low_frequency_mode_cutoff=None, zero_order_energy_corrections=None, check_overlap=None)`
    - `get_zero_order_energies(corrected_fundamental_freqs, states)` — :param corrected_fundamental_freqs:
  - **class `VPTRunner`**
    > A helper class to make it easier to run jobs by making the inputs/options
    > clear and making it easier to customize run options
    - `__init__(system, states, initial_states=None, hamiltonian_options=None, solver_options=None, runtime_options=None)`
    - `get_Hamiltonian()` — Build a `PerturbationTheoryHamiltonian` for this runner's system, combining the configured Hamilton…
    - `hamiltonian()` — The (cached) `PerturbationTheoryHamiltonian` for this runner, built lazily via `get_Hamiltonian` th…
    - `get_wavefunctions(**opts)` — Run the full VPT calculation and return the resulting wavefunctions, combining the solver options,…
    - `get_solver()` — Build a `PerturbationTheorySolver` for this runner's target states, without running the full wavefu…
    - `print_output_tables(wfns=None, file=None, print_intensities=True, print_energies=True, print_energy_corrections=True, print_transition_moments=True, operators=None, logger=None, sep_char='=', sep_len=100)` — Prints a bunch of formatted output data from a PT run
    - `print_tables(wfns=None, file=None, print_intensities=True, print_energy_corrections=True, print_transition_moments=True, operators=None, logger=None, sep_char='=', sep_len=100)` — Prints a bunch of formatted output data from a PT run
    - `get_Nielsen_energies(return_split=False, return_X=False, **potential_params)` — Compute harmonic and anharmonic (Nielsen-formula) vibrational energies for the target states, via t…
    - `print_Nielsen_frequencies(logger=None, state_formatting='vector', **potential_params)` — Compute the Nielsen-formula harmonic and anharmonic transition frequencies (relative to the ground…
    - `construct(system, states, target_property=None, extended_space_target_property=None, basis_filters=None, initial_states=None, corrected_fundamental_frequencies=None, **opts)` — Top-level constructor that assembles a fully configured `VPTRunner` (and, if target states are give…
    - `run_simple(system, states, target_property=None, corrected_fundamental_frequencies=None, calculate_intensities=True, plot_spectrum=False, operators=None, **opts)` — Makes a runner using the `construct` method and then calls that
    - **class `helpers`**
      > A stub to be replaced with the AnneInputHelpers interface
      - `run_anne_job(base_dir, states=2, calculate_intensities=None, return_analyzer=False, return_runner=False, modes_file=('nm_int.dat', 'modes.dat'), atoms_file='atom.dat', masses_file='mass.dat', coords_file='cart_ref.dat', zmat_file='z_mat.dat', potential_files=('cub.dat', 'quart.dat', 'quintic.dat', 'sextic.dat'), dipole_files=('lin_dip.dat', 'quad_dip.dat', 'cub_dip.dat', 'quart_dip.dat', 'quintic_dip.dat'), coordinate_transformation=None, coordinate_transformation_file='coordinate_transformation.py', results_file=None, order=None, expansion_order=None, energy_units=None, normalization_type=0, **opts)` — Stub placeholder documenting the intended `AnneInputHelpers.run_anne_job` interface for running a V…
  - **class `AnneInputHelpers`**
    - `get_tensor_idx(line, inds, m, start_at=0)` — Parse one line of a flattened force-constant/tensor data file into its (1-indexed-to-0-indexed) tup…
    - `parse_tensor(block, dims=None)` — Parse a symmetric force-constant tensor from an Anne-format data file (or its raw string content),…
    - `parse_dipole_tensor(block, dims=None)` — Parse a dipole-derivative tensor from an Anne-format data file (or its raw string content), treatin…
    - `parse_freqs_line(line)` — Parse a whitespace-separated line of frequency values into a NumPy array.
    - `parse_modes_line(line, nmodes)` — Parse a whitespace-separated, flattened block of mode-matrix values into a `(nmodes, ncols)` array,…
    - `parse_modes(block, energy_units=None)` — Parse the frequencies/mode matrix/inverse from an Anne-format modes file, or the first file that ex…
    - `parse_coords(block)` — Parse a Cartesian-coordinates data file (or its raw string content) into an `(natoms, 3)` array.
    - `parse_atoms(block)` — Parse an atomic-number data file (or its raw string content) into a list of element symbols.
    - `parse_masses(block)` — Parse a masses data file (or its raw string content) into a list of floating-point mass values.
    - `parse_zmatrix(block)` — Parse a Z-matrix connectivity data file (or its raw string content) into a standard `(natoms, 4)` Z…
    - `standard_sorting(zmat)` — converts from [r1, r2, r3, ..., a1, a2, ..., t1, t2, ...] coords
    - `get_internal_FG(freqs, modes, inv, sorting=None)` — Build the internal-coordinate force-constant (`F`) and kinetic-energy (`G`) matrices from a set of…
    - `renormalize_modes(freqs, modes, inv, sorting=None, type=2)` — Re-derive a set of normal modes (with a chosen dimensionality convention) from the internal-coordin…
    - `rerotate_force_field(old_inv, new_modes, old_field, dim_skips=0, sorting=None)` — Re-express a force-field expansion (a list of derivative tensors) from one mode basis to another, b…
    - `reexpress_normal_modes(base_modes, old_field, dipole, sorting=None, type=2)` — Re-express a potential (and, optionally, dipole) expansion from one normal-mode basis into a renorm…
    - `mass(atom)` — Look up an atom's mass (in atomic units) from its element symbol.
    - `extract_term_lists(checkpoint, terms, skip_dimensions=0, threshold=0, aggregator=None)` — Read a set of stored expansion terms out of a checkpoint file and flatten each one into a list of `…
    - `write_term_lists(terms, file_template=None, int_fmt='{:>3.0f}', float_fmt='{:>16.8e}', index_function=None)` — Write a set of flattened term-list data (as produced by `extract_term_lists`) out to text files (or…
    - `extract_terms(chk, out, terms, default_output='output.hdf5', aggregator=None, index_function=None, skip_dimensions=0)` — Extract and write out a named expansion's terms from a checkpoint (or a directory containing a defa…
    - `extract_potential(chk, out='potential_expansion_{}.dat')` — Extract and write out the potential-energy expansion terms from a checkpoint, via `extract_terms`.
    - `extract_gmatrix(chk, out='gmatrix_expansion_{}.dat')` — Extract and write out the G-matrix expansion terms from a checkpoint, via `extract_terms`.
    - `extract_dipole_expansion(chk, out='dipole_expansion_{}.dat')` — Extract and write out the dipole-derivative expansion terms from a checkpoint, via `extract_terms`,…
    - `run_anne_job(base_dir, states=2, initial_states=0, calculate_intensities=None, return_analyzer=False, return_runner=False, modes_file=('nm_int.dat', 'modes.dat'), atoms_file='atom.dat', masses_file='mass.dat', coords_file='cart_ref.dat', zmat_file='z_mat.dat', potential_files=('cub.dat', 'quart.dat', 'quintic.dat', 'sextic.dat'), dipole_files=('lin_dip.dat', 'quad_dip.dat', 'cub_dip.dat', 'quart_dip.dat', 'quintic_dip.dat'), coordinate_transformation=None, coordinate_transformation_file='coordinate_transformation.py', results_file=None, order=None, expansion_order=None, energy_units=None, normalization_type=0, **opts)` — Run a full VPT calculation from a directory of Anne-format input files (normal modes, atoms, masses…
    - `run_fchk_job(base_dir, states=2, calculate_intensities=None, return_analyzer=False, return_runner=False, zmat_file='z_mat.dat', fchk_file='fchk.fchk', results_file='output.hdf5', **opts)` — Run a VPT calculation using a Gaussian FChk file for the molecule/potential/dipole data, together w…
    - `get_internal_expansion(fchk, internals, states=2, **opts)` — Build a throwaway `VPTRunner` for a given FChk file and internal-coordinate specification, purely t…
    - `run_internal_expansion(expansion_data, calculate_intensities=False, **opts)` — Run a VPT calculation directly from a previously extracted internal-coordinate expansion (as produc…
  - **class `MultiVPTStateSpace`**
    > Generalizes a VPTStateSpace to pairs of initial and final spaces
    - `__init__(state_space_pairs, system=None, degeneracy_specs=None, evaluator=None, **opts)`
    - `state_list_pairs()` — The `[initial_state_list, target_state_list]` pairs for every `(initial, target)` space pair in thi…
  - **class `AnalyticVPTRunner`**
    - `__init__(expansions, order=None, expansion_order=None, freqs=None, internals=True, logger=None, hamiltonian=None, checkpoint=None, dipole_expansion=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None, local_mode_couplings=None, local_mode_coupling_order=None, parallelizer=None)`
    - `from_hamiltonian(ham, order, expansion_order=None, logger=None, checkpoint=None, parallelizer=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, take_diagonal_v4_terms=True, intermediate_normalization=None, corrected_fundamental_frequencies=None, **opts)` — A driver powered by a classic PerturbationTheoryHamiltonian object
    - `construct(system, states=None, *, order=2, expressions_file=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, mixed_derivative_handling_mode='analytical', degeneracy_specs=None, corrected_fundamental_frequencies=None, parallelizer=None, **settings)` — Build an `AnalyticVPTRunner` (and, if target states are given, its resolved `MultiVPTStateSpace`) d…
    - `from_file(file_name, order=2, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, expressions_file=None, **settings)` — Build an `AnalyticVPTRunner` from a molecule file (e.g.
    - `construct_classic_runner(states, system=None, logger=None, corrected_fundamental_frequencies=None, potential_terms=None, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, dipole_terms=None, initial_states=None, **opts)` — Build a classic `VPTRunner` reproducing this analytic evaluator's expansion data (rescaling the sto…
    - `clear_caches()` — Clear the global caches used by the underlying `AnalyticPerturbationTheorySolver` (e.g.
    - `prep_multispace(states, freqs, system=None, degeneracy_specs=None)` — Coerce a raw state specification into a `MultiVPTStateSpace`, passing an already-built one through…
    - `prep_states(states, degeneracy_specs=None)` — Coerce a raw state specification into a `MultiVPTStateSpace` using this evaluator's own frequencies…
    - `evaluate_expressions(states, exprs, zero_cutoff=None, operator_expansions=None, degeneracy_specs=None, verbose=False)` — Evaluate a set of arbitrary symbolic perturbation-theory expressions numerically for the given targ…
    - `get_matrix_corrections(states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False)` — Compute the perturbative matrix-element corrections for the given target states, via the underlying…
    - `get_energy_corrections(states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False)` — Compute the perturbative energy corrections for the given target states, via the underlying `Pertur…
    - `get_overlap_corrections(states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False)` — Compute the perturbative wavefunction-overlap corrections for the given target states, via the unde…
    - `prep_eval_state_pairs(states)` — Flatten a `MultiVPTStateSpace`'s `(initial, final)` state-list pairs into a flat list of `[single_i…
    - `get_full_wavefunction_corrections(states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False)` — Compute the full (all-component) perturbative wavefunction corrections for the given target states,…
    - `get_wavefunction_corrections(states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False)` — Compute the perturbative wavefunction corrections for the given target states, via the underlying `…
    - `unflatten_corr(states, corrs)` — Regroup a flat correction result (expressed over the combined initial/final state spaces) back into…
    - `get_operator_corrections(operator_expansion, states, order=None, terms=None, degeneracy_specs=None, verbose=False, operator_type=None, check_single=True, **opts)` — Compute the perturbative corrections to one or more arbitrary operator expansions for the given tar…
    - `construct_corrections_vectors(states, corrs)` — Assemble a set of flat per-order correction matrices spanning the full flat state space, from one o…
    - `construct_corrections_matrix(group, corrs)` — Assemble a set of square per-order correction matrices restricted to a single group of states, from…
    - `get_transition_moment_corrections(states, dipole_expansion=None, order=None, degeneracy_specs=None, axes=None, **opts)` — Compute the perturbative transition-dipole-moment corrections for the given target states, using th…
    - `get_freqs(states, order=None, degeneracy_specs=None, return_corrections=False, verbose=False)` — Compute the vibrational transition frequencies (in wavenumbers) for the given target states, relati…
    - `get_reexpressed_hamiltonian(states, order=None, degeneracy_specs=None, only_degenerate_terms=True, verbose=False, hamiltonian_corrections=None, **opts)` — Build the deperturbed (degenerate-block) Hamiltonian matrices for each degenerate group of states,…
    - `get_wfc_test_states(input_states, energy_window)` — Identify the candidate states that could plausibly be strongly coupled (via wavefunction-correction…
    - `get_test_wfn_corrs(input_states, energy_window)` — We take the expansions and frequencies that we have and at find the possible terms
    - `format_energies_table(states, energies, energy_corrections, zpe_pos, number_format='.3f')` — Format a table of state energies/frequencies alongside their per-order corrections, converting to w…
    - `format_degenerate_energies_table(states, energies, deperturbed_energies, zpe_pos, number_format='.3f')` — Format a table comparing each state's degenerate-perturbation-theory-corrected energy/frequency aga…
    - `format_transition_moment_table(states, transition_moments, transition_moment_corrections, number_format='.8f')` — Format a table of transition-dipole moments (and their per-order corrections) for each initial-stat…
    - `format_operators_table(states, keys, operator_values, operator_corrections, number_format='.8f')` — Format a table of arbitrary operator expectation values (and their per-order corrections) for each…
    - `format_spectrum_table(states, harmonic_spectra, spectra, deperturbed_spectra=None, number_format='.3f')` — Format a table of harmonic, anharmonic, and (optionally) deperturbed IR spectra for each initial-st…
    - `prep_operators(operator_expansions, operator_terms, order=None)` — Normalize a user-supplied operator specification (raw expansion coefficients, in either list or nam…
    - `format_matrix(ham)` — Format a matrix as a plain-text string using this class's standard print options (fixed precision,…
    - `modify_hamiltonian(hamiltonian_corrections)` — Build a new `AnalyticVPTRunner` whose underlying evaluator has extra Hamiltonian corrections applie…
    - `run_VPT(states, calculate_intensities=True, operator_expansions=None, operator_terms=None, operator_type=None, order=None, verbose=False, degeneracy_specs=None, handle_degeneracies=True, zero_cutoff=None, transition_moment_terms=None, hamiltonian_corrections=None, clear_caches=True, hamiltonian_correction_type=None, only_degenerate_terms=True, force_return_on_crash=True)` — Top-level entry point for running a full analytic VPT calculation on a set of target states: option…
    - `run_simple(system, states, calculate_intensities=True, operator_expansions=None, operator_terms=None, operator_type=None, verbose=False, return_runner=False, degeneracy_specs=None, degeneracy_states=None, handle_degeneracies=True, zero_cutoff=None, clear_caches=True, hamiltonian_correction_type=None, hamiltonian_corrections=None, only_degenerate_terms=True, force_return_on_crash=True, **opts)` — Convenience one-shot entry point: builds an `AnalyticVPTRunner` for the given system (via `construc…

### `Solver.py`
  - **class `PerturbationTheorySolver`**
    > A solver that applies perturbation theory
    > given a series of corrections and population of states.
    > Supports degenerate and non-degenerate PT.
    - `__init__(perturbations, states, coupled_states=None, order=2, total_space=None, flat_total_space=None, state_space_iterations=None, state_space_terms=None, state_space_filters=None, extended_state_space_filter_generator=None, extended_state_space_postprocessor=None, target_property_rules=None, allow_sakurai_degs=False, allow_post_PT_calc=True, modify_degenerate_perturbations=False, gaussian_resonance_handling=False, ignore_odd_order_energies=False, intermediate_normalization=False, reorthogonalize_degenerate_states=None, check_overlap=True, zero_element_warning=True, degenerate_states=None, degeneracy_handlers=None, handle_strong_couplings=False, local_coupling_hamiltonian=None, local_coupling_order=None, low_frequency_mode_cutoff=0.00115, zero_order_energy_corrections=None, nondeg_hamiltonian_precision=3, memory_constrained=False, keep_hamiltonians=None, logger=None, parallelizer=None, checkpointer=None, results=None, checkpoint_keys=None, use_cached_representations=False, use_cached_basis=False)`
    - `coupled_states()` — :return:
    - `total_space_dim()` — :return:
    - `flat_total_space()` — :return:
    - `total_state_space()` — :return:
    - **class `PastIndexableTuple`** (list)
    - `representations()` — :return:
    - `representations(reps)` — The (cached) list of Hamiltonian perturbation representation matrices, computed lazily via `get_VPT…
    - `degenerate_spaces()` — :return:
    - `merge_deg_spaces(new_states)` — Combine several independently-identified sets of degenerate state groups into one consistent set, b…
    - `zero_order_energies()` — :return:
    - `apply_VPT()` — Applies perturbation theory to the held basis of states using the
    - `get_VPT_representations()` — Gets the sparse representations of the passed perturbation inside the basis of coupled states.
    - `extend_VPT_representations(new_flat_space, new_states)` — Extend the cached Hamiltonian perturbation representation matrices to cover a newly added block of…
    - `load_state_spaces()` — :return:
    - `extend_state_spaces(new_targets, degenerate_states=None)` — Extend the solver's state spaces to additionally include a set of new target states: (re)computes t…
    - `load_coupled_spaces(degenerate_spaces=None, spaces=None, wavefunction_terms=None, property_filter=None, filter_spaces=None)` — Determines which states need to be coupled at which levels of correction
    - **class `StateSpaceWrapper`**
      > Wraps a state space so that it can define stuff like __add__, __mul__, and __neg__
      - `__init__(space)`
      - `simple_union(other)` — Combine this wrapped state space with another (treating `None`/`0` as "nothing to add"), unwrapping…
    - **class `ProjectionOperatorWrapper`**
      > Generates a symbolic form of a perturbation operator that
      > either projects onto a degenerate space or projects it out
      - `__init__(space, complement=False)`
      - `get_transformed_space(other)` — :param other:
    - **class `ProjectedOperator`**
      > Generates a symbolic form of an operator where
      > an operator can first be applied, then unused terms projected
      > out, before returning the state space
      - `__init__(projector, operator)`
      - `get_transformed_space(other)` — :param other:
    - `get_coupled_space(input_state_space, degenerate_space, use_second_deg, allow_PT_degs=True, wavefunction_terms=None, spaces=None, property_filter=None, filter_spaces=None)` — Applies the VPT equations semi-symbolically, dispatching based on how many
    - `get_nondeg_coupled_space(input_state_space, degenerate_space=None, spaces=None, wavefunction_terms=None, property_filter=None, filter_spaces=None)` — Applies the non-degenerate equations in semi-symbolic form to determine
    - `get_corrections(non_zero_cutoff=None, handle_strong_couplings=None, check_overlap=None)` — Applies the perturbation theory equations to obtain
    - `high_frequency_modes()` — The indices of the vibrational modes whose fundamental transition frequency exceeds `self.low_frequ…
    - `identify_strong_couplings(corrs)` — Find the strongly-coupled state pairs among the solved corrections, delegating to whichever configu…
    - `construct_strong_coupling_spaces(spec, sc, corrs, states, threshold)` — Build the degenerate state groups implied by a set of identified strong couplings, and, if the spec…
    - `apply_VPT_equations(state_index, degenerate_space_indices, degenerate_energies, zero_order_state, degenerate_subspace, degenerate_subsubspace, perturbations=None, allow_PT_degs=None, ignore_odd_orders=None, intermediate_normalization=None, non_zero_cutoff=None)` — Applies VPT equations, dispatching based on how many
    - `apply_VPT_nondeg_equations(state_index, deg_group, perturbations=None, non_zero_cutoff=None, check_overlap=None, intermediate_normalization=False, ignore_odd_orders=False)` — Does the dirty work of doing the VPT iterative equations.
    - `apply_VPT_2k1_rules(existing_corrs, perturbations=None)` — Apply expressions allowing for obtaining higher-order
    - `apply_post_PT_variational_calc(degenerate_states, corrs)` — Applies degenerate perturbation theory by building a representation
    - `drop_deg_pert_els(perts, deg_groups)` — :param perts:

### `Terms.py` — Stores all of the terms used inside the VPT2 representations
  - **class `DumbTensor`**
    > A wrapper to make tensor algebra suck less
    - `__init__(tensor)`
    - `shape()` — **LLM Docstring**
    - `dot(b, *args, **kwargs)` — Contract this tensor with `b` using `DumbTensor._dot`, unwrapping `b` first if it is itself a `Dumb…
    - `shift(*args, **kwargs)` — Apply `_shift` to this tensor's data and wrap the result in a new `DumbTensor`.
    - `transpose(*perm)` — Transpose the wrapped tensor according to `perm` and wrap the result in a new `DumbTensor`.
    - `contract_dim(targ_dim)` — Apply `_contract_dim` to this tensor's data and wrap the result in a new `DumbTensor`.
  - **class `MixedDerivativeHandlingModes`** (enum.Enum)
  - **class `ImaginaryFrequencyHandlingMode`** (enum.Enum)
  - **class `JacobianKeys`** (enum.Enum)
    > Real access pattern: JacobianKeys.<MemberName> (this is an enum with 12 members, e.g. JacobianKeys.CartesiansByInternals == 'CartesiansByInternals'). Collapsed into a dict below purely for compactness -- do not index it as a dict in real code:
  - **class `ExpansionTerms`**
    > Base class for kinetic, potential, and dipole derivative terms
    - `__init__(molecule, modes=None, mode_selection=None, mode_transformation=None, use_internal_modes=None, logger=None, parallelizer=None, checkpointer=None, undimensionalize=None, numerical_jacobians=True, eckart_embed_derivatives=True, eckart_embed_planar_ref_tolerance=None, strip_dummies=False, strip_embedding=True, mixed_derivative_handling_mode=None, mixed_derivative_warning_threshold=0.00025, mixed_derivative_handle_zeros=False, backpropagate_internals=False, direct_propagate_cartesians=False, zero_mass_term=10000000.0, expansion_handling_mode='old', internal_fd_mesh_spacing=0.01, internal_fd_stencil=None, cartesian_fd_mesh_spacing=0.01, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=None, cartesian_by_internal_derivative_method=None, internal_by_cartesian_order=3, cartesian_by_internal_order=4, jacobian_warning_threshold=10000.0, coordinate_transformations=None, coordinate_derivatives=None, imaginary_frequency_handling_mode='abs')`
    - `num_atoms()` — Gets the number of atoms (excluding dummies if `strip_dummies` is `True`)
    - `modes()` — The stored mode object (in whatever basis -- Cartesian or internal -- it was constructed with).
    - `get_terms(order=None)` — Gets the terms up to the given order
    - `get_term(t)` — Provides the term at order `t`
    - `terms()` — The (cached) full set of expansion terms, computed lazily via `get_terms()` the first time they're…
    - `get_int_jacobs(jacs)` — Gets the specified Internal->Cartesian Jacobians
    - `get_cart_jacobs(jacs)` — Gets the specified Cartesian->Internal Jacobians
    - `inertial_frame()` — Provides the inertial axis frame
    - `inertial_frame_derivatives()` — Compute the first and second derivatives of the (mass-weighted) inertia tensor with respect to mass…
    - `moment_of_inertia_derivs(order)` — Compute the Taylor-series derivatives of the inverse inertia tensor with respect to the normal-mode…
    - `get_coordinate_transforms(internal_by_cartesian_order=None, cartesian_by_internal_order=None, current_cache=None)` — Compute (and cache, per-molecule, both in memory and via the checkpointer) the full set of Jacobian…
    - `cartesian_L_matrix()` — :return: the leading term of `get_cartesians_by_cartesian_modes(1)`
    - `get_cartesians_by_cartesian_modes(order=None)` — Fetch the Cartesians-by-Cartesian-normal-modes Jacobians up to the requested order, computing them…
    - `cartesian_L_inverse()` — :return: the leading term of `get_cartesian_modes_by_cartesians(1)`
    - `get_cartesian_modes_by_cartesians(order=None)` — Fetch the Cartesian-normal-modes-by-Cartesians Jacobians up to the requested order, computing them…
    - `internal_L_matrix()` — :return: the leading term of `get_internal_modes_by_internals(1)`
    - `get_internal_modes_by_internals(order=None, strip_embedding=True)` — Fetch the internal-normal-modes-by-internals Jacobians up to the requested order, optionally stripp…
    - `internal_L_inverse()` — :return: the leading term of `get_internals_by_internal_modes(1)`
    - `get_internals_by_internal_modes(order=None, strip_embedding=True)` — Fetch the internals-by-internal-normal-modes Jacobians up to the requested order, optionally stripp…
    - `cartesians_by_modes()` — All cached Cartesians-by-internal-modes Jacobians, computing the default set if not already cached.
    - `get_cartesians_by_modes(order=None)` — Fetch the Cartesians-by-internal-normal-modes Jacobians up to the requested order, computing them i…
    - `modes_by_cartesians()` — All cached internal-normal-modes-by-Cartesians Jacobians, computing the default set if not already…
    - `get_modes_by_cartesians(order=None, strip_embedding=True)` — Fetch the internal-normal-modes-by-Cartesians Jacobians up to the requested order.
    - `cartesians_by_internals()` — All cached Cartesians-by-internals Jacobians, computing the default set if not already cached.
    - `get_cartesians_by_internals(order=None, strip_embedding=False)` — Fetch the Cartesians-by-internals Jacobians up to the requested order, optionally stripping embeddi…
    - `internals_by_cartesians()` — All cached internals-by-Cartesians Jacobians, computing the default set if not already cached.
    - `get_internals_by_cartesians(order=None, strip_embedding=False)` — Fetch the internals-by-Cartesians Jacobians up to the requested order, optionally stripping embeddi…
    - `cartesian_modes_by_internal_modes()` — All cached Cartesian-normal-modes-by-internal-normal-modes Jacobians, computing the default set if…
    - `get_cartesian_modes_by_internal_modes(order=None)` — Fetch the Cartesian-normal-modes-by-internal-normal-modes Jacobians up to the requested order.
    - `internal_modes_by_cartesian_modes()` — All cached internal-normal-modes-by-Cartesian-normal-modes Jacobians, computing the default set if…
    - `get_internal_modes_by_cartesian_modes(order=None)` — Fetch the internal-normal-modes-by-Cartesian-normal-modes Jacobians up to the requested order.
  - **class `PotentialTerms`** (ExpansionTerms)
    > A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates
    - `__init__(molecule, mixed_derivs=None, modes=None, potential_derivatives=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, logger=None, parallelizer=None, checkpointer=None, check_input_force_constants=True, allow_higher_potential_terms=False, hessian_tolerance=0.0001, grad_tolerance=0.0001, freq_tolerance=0.002, **opts)`
    - `v_derivs()` — Property getter/setter for the canonicalized potential-energy derivative tensors.
    - `v_derivs(v)` — Property getter/setter for the canonicalized potential-energy derivative tensors.
    - `get_terms(order=None, logger=None)` — Compute the potential-energy expansion terms in the molecule's normal-mode coordinates, zeroing out…
    - `get_potential_optimized_coordinates(V_expansion, order=2)` — Find, order by order, the coordinate transformation that eliminates as much of the potential-energy…
    - `optimize_coordinates(order=2)` — Build the potential-optimized coordinate transformation for this potential expansion and re-express…
  - **class `KineticTerms`** (ExpansionTerms)
    > Represents the KE coefficients
    - `__init__(molecule, g_derivative_threshold=0.001, gmatrix_tolerance=1e-06, use_cartesian_kinetic_energy=False, check_input_gmatrix=True, freq_tolerance=0.002, **opts)`
    - `get_terms(order=None, logger=None, return_expressions=False)` — Compute the G-matrix (kinetic-energy coefficient) Taylor-expansion terms in the molecule's normal-m…
    - `reexpress_G(G_expansion, forward_derivs, reverse_derivs=None, order=2)` — Apply a coordinate transformation to the G-matrix
    - `reexpress(forward_derivs, reverse_derivs=None, order=2)` — Finds a coordinate transformation the give 0 contribution to the G-matrix
    - `get_kinetic_optimized_coordinates(G_expansion, order=2)` — Intended to iteratively find the coordinate transformation that eliminates as much of the G-matrix…
    - `optimize_coordinates(order=2)` — Build the kinetic-energy-optimized coordinate transformation for this G-matrix expansion and re-exp…
  - **class `CoriolisTerm`** (ExpansionTerms)
    > Calculates the Coriolis coupling term
    - `get_zetas_and_momi()` — Compute the Coriolis zeta constants and the inertial-frame moments of inertia: mass-weights and fre…
    - `get_zetas()` — The Coriolis zeta-constant tensor alone, via `get_zetas_and_momi`.
    - `get_terms(order=None, J=0)` — Compute the Coriolis rotational-coupling operator's Taylor-expansion terms by combining the frequen…
  - **class `PotentialLikeTerm`** (KineticTerms)
    > This accounts for the potential-like term.
    > In Cartesian diplacement modes this is the Watson U.
    > In proper internals, this is the V' term.
    - `get_terms(order=None, logger=None)` — Compute the potential-like (Watson `U` / internal-coordinate `V'`) correction term's Taylor-expansi…
  - **class `DipoleTerms`** (ExpansionTerms)
    - `__init__(molecule, dipole_derivatives=None, mixed_derivs=None, modes=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, logger=None, parallelizer=None, checkpointer=None, **opts)`
    - `get_terms(order=None, logger=None)` — Compute the dipole-moment Taylor-expansion terms (one component at a time) in the molecule's normal…
  - **class `OperatorTerms`** (ExpansionTerms)
    > Literally as simple as it comes for an operator expansion.
    > One dimensional, no mixed derivative stuff.
    - `__init__(molecule, operator_derivatives=None, modes=None, mode_selection=None, logger=None, parallelizer=None, checkpointer=None, **opts)`
    - `get_terms(order=None, logger=None)` — Compute the operator's Taylor-expansion terms in the molecule's normal-mode coordinates, returning…

### `Wavefunctions.py` — Provides classes to support wave functions coming from VPT calculations
  - **class `PerturbationTheoryWavefunctions`** (Wavefunctions)
    > These things are fed the first and second order corrections
    - `__init__(mol, basis, corrections, initial_states=None, modes=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, logger=None, checkpoint=None, results=None, operator_settings=None, expansion_options=None, degenerate_transformation_layout=None)`
    - `get_dimension()`
    - `energies()`
    - `energies(engs)`
    - `to_state(serializer=None)`
    - `from_state(data, serializer=None)`
    - `degenerate_transformation()`
    - `initial_state_indices()`
    - `energies_to_order(order=None)`
    - `deperturbed_energies()`
    - `deperturbed_frequencies(order=None)`
    - `order()`
    - `expectation(operator, other=None)`
    - `zero_order_energies()`
    - `get_M0(mu_0)`
    - `get_M1(mu_1)`
    - `get_M2(mu_2)`
    - `get_M3(mu_3)`
    - `get_Mi(i, mu, base_sym='M')`
    - **class `TermHolder`** (tuple)
      > symbolic wrapper on Tuple so we can know that we've canonicalized some term
    - `dipole_terms()`
    - `dipole_terms(v)`
    - **class `DipolePartitioningMethod`** (enum.Enum)
    - `dipole_partitioning()`
    - `dipole_partitioning(p)`
    - `transition_moments()` — Computes the transition moments between wavefunctions stored in the object
    - `deperturbed_transition_moments()` — Computes the transition moments between wavefunctions stored in the object
    - `transition_moment_corrections()` — Computes the transition moment corrections between wavefunctions stored in the object
    - `deperturbed_transition_moment_corrections()` — Computes the transition moment corrections between wavefunctions stored in the object
    - `oscillator_strengths()` — Computes the oscillator strengths for transitions from the ground state to the other states
    - `deperturbed_oscillator_strengths()` — Computes the oscillator strengths for transitions from the ground state to the other states
    - `oscillator_strengths_to_order(order)` — :param tms:
    - `deperturbed_oscillator_strengths_to_order(order)` — :param tms:
    - `intensities()` — Computes the intensities for transitions from the ground state to the other states
    - `deperturbed_intensities()` — Computes the intensities for transitions from the ground state to the other states
    - `intensities_to_order(order, return_freqs=False)` — Computes the intensities for transitions from the ground state to the other states
    - `deperturbed_intensities_to_order(order, return_freqs=False)` — Computes the intensities for transitions from the ground state to the other states
    - `zero_order_intensities()` — Computes the harmonic intensities for transitions from the ground state to the other states
    - `prep_operator_terms(coeffs)`
    - `operator_correction_data(operator_coeffs, order=None)`
    - `generate_intensity_breakdown(include_wavefunctions=True)` — Generates a breakdown of the terms that contribute to the intensity
    - `write_CSV_breakdown(file, intensity_breakdown, padding=None)` — Writes an intensity breakdown to a CSV by annoyingly flattening all the arrays
    - `format_energies_table(states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}')`
    - `format_deperturbed_energies_table(states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}')`
    - `format_property_matrices(states, prop_corrs, real_fmt='{:>.8e}', padding_fmt='{:>16}')`
    - `format_dipole_contribs_tables()`
    - `format_deperturbed_dipole_contribs_tables()`
    - `format_energy_corrections_table(real_fmt='{:>12.5f}')`
    - `format_intensities_table(real_fmt='{:>12.5f}')`
    - `format_deperturbed_intensities_table(real_fmt='{:>12.5f}')`
    - `get_spectrum()`
    - `get_deperturbed_spectrum()`
    - `format_operator_table(operators, real_fmt='{:>12.5f}')`