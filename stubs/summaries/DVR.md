### `BaseDVR.py` — Redoes what was originally PyDVR but in the _right_ way using proper subclassing and abstract prope…
  - **class `BaseDVR`**
    > Provides the abstract interface for creating a
    > convenient runnable DVR that can be cleanly subclassed to provide
    > extensions
    - `__init__(domain=None, divs=None, potential_function=None, logger=None, **base_opts)`
    - `get_grid(domain=None, divs=None, **kwargs)` — Abstract hook for building this DVR's 1D grid over the given domain/division count.
    - `grid(domain=None, divs=None, **kwargs)` — Build (or retrieve) this DVR's grid, falling back to the stored `domain`/`divs` if not given explic…
    - `get_kinetic_energy(grid=None, mass=None, hb=1, **kwargs)` — Abstract hook for building this DVR's 1D kinetic-energy operator matrix on the given grid.
    - `handle_kinetic_coupling(grid, ke_1D, g, g_deriv, hb=1, logger=None, **kwargs)` — Apply a (possibly coordinate-dependent) kinetic-coupling correction to a 1D kinetic-energy matrix:…
    - `kinetic_energy(grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs)` — Build the full kinetic-energy operator, computing the base 1D operator (via `get_kinetic_energy`) a…
    - `real_momentum(grid=None, mass=None, hb=1, **kwargs)` — Abstract hook for the real part of the momentum-operator matrix on the given grid.
    - `potential_energy(grid=None, potential_function=None, potential_values=None, potential_grid=None, logger=None, **pars)` — Calculates the potential energy at the grid points based
    - `hamiltonian(kinetic_energy=None, potential_energy=None, potential_threshold=None, **pars)` — Calculates the total Hamiltonian from the kinetic and potential matrices
    - `wavefunctions(hamiltonian=None, num_wfns=25, nodeless_ground_state=False, diag_mode=None, logger=None, **pars)` — Calculates the wavefunctions for the given Hamiltonian.
    - `run(result='wavefunctions', logger=None, grid=None, potential_energy=None, kinetic_energy=None, hamiltonian=None, **opts)` — :return:
  - **class `DVRResults`**
    > A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension
    - `__init__(grid=None, kinetic_energy=None, potential_energy=None, hamiltonian=None, wavefunctions=None, parent=None, **opts)`
    - `dimension()` — The number of spatial dimensions of the underlying DVR grid.
    - `plot_potential(plot_class=None, figure=None, plot_units=None, energy_threshold=None, zero_shift=False, **opts)` — Simple plotting function for the potential.
  - **class `DVRException`** (Exception)
    > Base exception class for working with DVRs

### `ColbertMiller.py` — Provides implementations of the Colbert-Miller DVR types defined in
  - **class `CartesianDVR`** (BaseDVR)
    > Provides the Colbert Miller DVR on the Cartesian [-inf, inf] range
    - `get_grid(domain=None, divs=None, **kw)` — Provides the Colbert-Miller DVR grid for the [-inf, inf] range
    - `get_kinetic_energy(grid=None, mass=None, hb=1, **kwargs)` — Build the Colbert-Miller kinetic-energy matrix for the free-particle (Cartesian, `[-inf, inf]`) DVR…
    - `real_momentum(grid=None, mass=None, hb=1, **kwargs)` — Provides the real part of the momentum for the [-inf, inf] range
  - **class `RingDVR`** (BaseDVR)
    > Provides a DVR for working on the (0, 2Pi) range with periodicity from Colbert and Miller
    - `__init__(domain=None, **opts)`
    - `get_grid(domain=None, divs=None, **kw)` — Provides the Colbert-Miller 1D grid for the [0, 2Pi] range
    - `get_kinetic_energy(grid=None, mass=1, hb=1, **kw)` — Colbert-Miller kinetic energy for the [0, 2pi] range
    - `real_momentum(grid=None, hb=1, **kw)` — Provides the real part of the momentum for the [0, 2pi] range
  - **class `PolarDVR`** (BaseDVR)
    > Provides a DVR for working on the (0, pi) range from Colbert and Miller
    - `__init__(domain=None, **opts)`
    - `get_grid(domain=(0, np.pi), divs=None, **kwargs)` — Provides the grid appropriate for the Colbert-Miller (0, Pi) range
    - `get_kinetic_energy(grid=None, mass=None, hb=1, **kwargs)` — Colbert-Miller kinetic energy for the [0, pi] range
  - **class `RadialDVR`** (BaseDVR)
    > Provides a DVR for working on the (0, inf) range from Colbert and Miller
    - `get_grid(domain=(0, np.pi), divs=None, **kwargs)` — Provides the grid appropriate for the Colbert-Miller (0, Pi) range
    - `get_kinetic_energy(grid=None, mass=None, hb=1, **kwargs)` — Colbert-Miller kinetic energy for the (0, inf) range

### `DVR.py`
  - **class `DVRConstructor`**
    - `load_domain_map()` — Build the default mapping from special-cased coordinate domains to the appropriate specialized 1D D…
    - `infer_DVR_type(domain)` — Infer which specialized 1D DVR class to use for a given coordinate domain, matching it (via `np.all…
    - `construct(domain=None, divs=None, potential_function=None, g=None, g_deriv=None, mass=None, po_divs=25, classes=None, scf=False, potential_optimize=False, logger=True, **base_opts)` — Dispatch to build an appropriate DVR object from a domain/division specification, either a single 1…
- `DVR(domain=None, divs=None, classes=None, potential_function=None, g=None, g_deriv=None, scf=False, potential_optimize=False, **base_opts)` — Constructs a DVR object

### `DirectProduct.py` — Provides a direct-product extension to a system of 1D DVRs
  - **class `DirectProductDVR`** (BaseDVR)
    - `__init__(dvrs_1D, zero_threshold=1e-14, **base_opts)`
    - `get_grid(domain=None, divs=None, **kwargs)` — Build the full multi-dimensional DVR grid as the Cartesian product of each component 1D DVR's own g…
    - `grid(domain=None, divs=None, **kwargs)` — Build (or retrieve) the multi-dimensional DVR grid, falling back to this DVR's stored `domain`/`div…
    - `get_kinetic_energy(grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, include_kinetic_coupling=True, **kwargs)` — Assemble the full multi-dimensional kinetic-energy operator as a sparse matrix, either as a simple…
    - `kinetic_energy(grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, **kwargs)` — Computes the N-dimensional kinetic energy
  - **class `CartesianNDDVR`** (DirectProductDVR)
    > Provides an ND-DVR over different domains
    - `__init__(domains, **base_opts)`
  - **class `RingNDDVR`** (DirectProductDVR)
    > Provides an ND-DVR for products of periodic (0, 2Pi) ranges
    - `__init__(divs, **base_opts)`
  - **class `SphericalDVR`** (DirectProductDVR)
    - `__init__(r_max, divs, **base_opts)`

### `Extensions.py`
  - **class `SCFWavefunctionGenerator`**
    - `__init__(dvr_1D)`
  - **class `SelfConsistentDVR`** (GridSCF)
    - `__init__(base_dvr, **opts)`
    - `create_grid_vals()` — Compute the full grid and potential-energy values for the underlying multi-dimensional DVR, used to…
    - `create_solvers(grid, pe)` — Build the per-dimension SCF wavefunction generators, first rebinding each 1D DVR's `g`/`g_deriv` ki…
  - **class `PotentialOptimizedDVR`** (DirectProductDVR)
    - `__init__(wfns_1D, **base_opts)`
    - `from_minimum(base_dvr, **opts)` — Build a `PotentialOptimizedDVR` using the SCF wavefunctions computed from an initial (unconverged,…
    - `from_scf(scf_dvr, wfns=None, **opts)` — Build a `PotentialOptimizedDVR` using the wavefunctions from a fully-converged SCF run (or an expli…

### `FiniteBasisDVR.py`
  - **class `InitialBasis`** (typing.Protocol)
    - `dimensions()` — Protocol property for the number of basis functions along each dimension of the underlying finite b…
    - `x(n)` — Protocol stub for the position-operator matrix representation in the first `n` basis functions.
    - `p2(n)` — Protocol stub for the momentum-squared operator matrix representation in the first `n` basis functi…
    - `p(n)` — Protocol stub for the momentum-operator matrix representation in the first `n` basis functions.
  - **class `FiniteBasisDVR`** (BaseDVR)
    - `__init__(basis, domain=None, divs=None, **opts)`
    - `get_grid(domain=None, divs=None, **kwargs)` — Build the DVR grid by diagonalizing the position-operator matrix of the underlying finite basis.
    - `real_momentum(grid=None, mass=None, hb=1, **kwargs)` — Build the momentum-operator matrix in the DVR grid basis, by transforming the initial basis's momen…
    - `get_kinetic_energy(grid=None, mass=None, hb=1, **kwargs)` — Build the kinetic-energy matrix in the DVR grid basis, using the initial basis's momentum-squared o…
    - `potential_energy(grid=None, potential_function=None, potential_values=None, potential_grid=None, logger=None, **pars)` — Compute the potential-energy matrix using just the grid points (rather than the full `(grid_points,…
  - **class `HarmonicDVR`** (FiniteBasisDVR)
    - `__init__(divs=None, **opts)`
  - **class `WavefunctionBasisDVR`** (FiniteBasisDVR)
    - `__init__(wfns=None, **opts)`

### `Wavefunctions.py` — Provides a DVRWavefunction class that inherits from the base Psience wavefunction
  - **class `DVRWavefunction`** (Wavefunction)
    - `__init__(energy, data, parent=None, grid=None, index=None, **opts)`
    - `get_dimension()` — The number of degrees of freedom this wavefunction is defined over, inferred from the trailing dime…
    - `plot(figure=None, grid=None, **opts)` — Plot the wavefunction using its own DVR grid and values, delegating to the base `Wavefunction.plot`.
    - `expectation(op, other=None)` — Computes the expectation value of operator op over the wavefunction other and self
    - `interp()` — A (lazily built and cached) continuous interpolant of the wavefunction's grid values, used by `eval…
    - `evaluate(points)` — Evaluates the functions at the given points
    - `marginalize_out(dofs)` — Computes the projection of the current wavefunction onto a set of degrees
  - **class `DVRWavefunctions`** (Wavefunctions)
    - `__init__(energies=None, wavefunctions=None, grid=None, results=None, **opts)`
    - `plot(figure=None, **opts)` — Plots the held wavefunctions
    - `expectation(op, other=None, multiplicative=True)` — Computes the expectation value of operator op over the wavefunction other and self
    - `transform_operator(M)` — Transform an operator matrix given in the DVR grid-point basis into the basis of these wavefunction…
    - `coordinate()` — The position-operator matrix in the wavefunction basis, computed as the expectation value of the DV…
    - `momentum()` — The real part of the momentum-operator matrix in the wavefunction basis, computed from the underlyi…
    - `laplacian()` — The Laplacian operator matrix in the wavefunction basis, derived from a fresh (unit-mass, zero-pote…
    - `kinetic_energy()` — The kinetic-energy operator matrix in the wavefunction basis, transformed from the stored `DVRResul…
    - `potential_energy()` — The potential-energy operator matrix in the wavefunction basis, transformed from the stored `DVRRes…