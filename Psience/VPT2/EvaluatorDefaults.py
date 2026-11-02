"""Named runner defaults. Explicit settings always take precedence.

The legacy profile selects supported legacy algorithms, not historical bugs.
Derivative sub-options are resolved by ExpansionTerms from expansion_handling_mode.
"""

__all__ = ['resolve_evaluator_defaults']

_PROFILES = {
    'new': dict(
        expansion_handling_mode='new',
        use_block_recurrence=True,
        operator_contraction_method='direct',
        factorized_coriolis=True,
        state_space_storage='auto',
        strong_coupling_property_filter=True
    ),
    'old': dict(
        expansion_handling_mode='old',
        use_block_recurrence=False,
        operator_contraction_method='tensor',
        factorized_coriolis=False,
        state_space_storage='memory',
        strong_coupling_property_filter=False
    )
}


def resolve_evaluator_defaults(evaluator_defaults='new', options=None):
    """Return a fresh option map; None retains the selected profile's default."""
    if not isinstance(evaluator_defaults, str) or evaluator_defaults not in _PROFILES:
        raise ValueError("evaluator_defaults must be 'new' or 'old'")
    result = dict(_PROFILES[evaluator_defaults])
    if options is not None:
        result.update((key, value) for key, value in options.items()
                      if value is not None or key not in result)
    return result
