"""
An implementation of vibrational perturbation theory (VPT) that uses sparse matrix methods to obtain
corrections.
Makes heavy use of the `BasisReps` package as well as `McUtils.Coordinerds` to obtain representations
of the corrections to the vibrational Hamiltonian.
Is technically not restricted to VPT in a harmonic basis, but no other forms of PT are likely to
be implemented in the near future.
For the purposes of papers, we've been calling this implementation `PyVibPTn`.

The easiest way to run jobs is through the `VPTRunner.run_simple` interface.
The options for jobs along with short descriptions are detailed in
[`VPTSystem`](VPT2/Runner/VPTSystem.md) for molecule/system-related options,
[`VPTStateSpace`](VPT2/Runner/VPTStateSpace.md) for state-space & degeneracy-related options,
[`VPTHamiltonianOptions`](VPT2/Runner/VPTHamiltonianOptions.md) for expansion-related options,
[`VPTHamiltonianOptions`](VPT2/Runner/VPTHamiltonianOptions.md) for expansion-related options,
[`VPTSolverOptions`](VPT2/Runner/VPTSolverOptions.md) for options related to constructing representations and applying VPT,
and [`VPTRuntimeOptions`](VPT2/Runner/VPTRuntimeOptions.md) for options related to the runtime/logging

**A basic tutorial to provide more extensive hand-holding can be found [here](VPT2/tutorial.md).**

Finally, the general code flow is detailed below

![pt design](/Psience/img/PyVibPTnDesign.png){:width="100%"}

:long_description: The implementation of vibrational perturbation theory provided here uses a kernel/config/driver type of design.
The kernel that actually solves the perturbation theory equations is the [`PerturbationTheorySolver`](PerturbationTheorySolver.md) object.
The config comes the [`PerturbationTheoryHamiltonian`](PerturbationTheoryHamiltonian.md), which implements the expansion of the Hamiltonian
with respect to normal modes.
Finally, the driver is the [`VPTRunner`](VPTRunner.md) which through its `run_simple` method aggregates all of the possible options needed
for VPT and sends them to the right parts of the architecture.

Because of this, while it can be helpful to know how the `PerturbationTheorySolver` and `PerturbationTheoryHamiltonian` work, all that one
usually needs to do to run VPT is to call `VPTRunner.run_simple`.
The most basic call looks like
```python
VPTRunner.run_simple(system_spec, states)
```
where the `system_spec` is often an `fchk` file from an electronic structure calculation that provides a restricted quartic potential
and the `states` is an `int` that specifies the max number of quanta of excitation in the states that will be corrected.

The system spec can also be more explicit, being either a `Molecule` object or a list like `[atoms, coords, opts]`
where `atoms` is a list of strings of the atoms, `coords` are the accompanying Cartesian coordinates, and `opts` is a `dict`
of options for the molecule, such as the `masses`.
It is also possible (and sometimes necessary) to supply custom normal modes, which can be done through the `modes` option (examples below).
"""

__all__ = []
from .Runner import *; from .Runner import __all__ as exposed
__all__ += exposed
from .Analyzer import *; from .Analyzer import __all__ as exposed
__all__ += exposed
from .Hamiltonian import *; from .Hamiltonian import __all__ as exposed
__all__ += exposed
from .Solver import *; from .Solver import __all__ as exposed
__all__ += exposed
from .Corrections import *; from .Corrections import __all__ as exposed
__all__ += exposed
from .Wavefunctions import *; from .Wavefunctions import __all__ as exposed
__all__ += exposed
from .Terms import *; from .Terms import __all__ as exposed
__all__ += exposed