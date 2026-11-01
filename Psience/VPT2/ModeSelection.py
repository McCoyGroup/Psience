"""Exact selection of interaction mode subsets in analytic contractions."""
import dataclasses,itertools,math

@dataclasses.dataclass(frozen=True)
class RequiredModeSubsets:
    """Retain subsets containing every specified original mode position.

    This selects contributions, not states, coefficients, or denominators.
    No perturbation order is assumed. The evaluator currently restricts this
    interface to diagonal requests in the original mode ordering.
    """
    required: tuple

    def __post_init__(self):
        values=tuple(int(i) for i in self.required)
        if any(i<0 or i!=j for i,j in zip(values,self.required)) or len(set(values))!=len(values):
            raise ValueError('required mode positions must be distinct nonnegative integers')
        object.__setattr__(self,'required',tuple(sorted(values)))

    def _split(self,ndim,num_fixed,free_count):
        if any(i>=ndim for i in self.required):raise ValueError('required mode position outside basis')
        required=tuple(i for i in self.required if i>=num_fixed)
        available=tuple(i for i in range(num_fixed,ndim) if i not in required)
        choose=free_count-len(required)
        return required,available,choose

    def count(self,ndim,num_fixed,free_count):
        required,available,choose=self._split(ndim,num_fixed,free_count)
        return 0 if choose<0 or choose>len(available) else math.comb(len(available),choose)

    def iter_subsets(self,ndim,num_fixed,free_count):
        required,available,choose=self._split(ndim,num_fixed,free_count)
        if choose<0 or choose>len(available):return
        for extra in itertools.combinations(available,choose):
            yield tuple(sorted(required+extra))
