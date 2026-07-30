### `AnalyticJacobianDotCalculator.py` — Temporary storage for an attempt to get G matrix elements purely symbolically
  - **class `OrientationVector`**
    - `__init__(x, y, z, basis)`
    - `dot(other)`
    - `rotation_matrix(angle, axis)`
  - **class `RotationMatrix`**
    - `__init__(v1, v2, v3)`
    - `dot(v)`
  - **class `BondVector`**
    - `__init__(i, j, norm=1, embedding=None)`
    - `one_common_atom(bond, mode='shared')`
    - `shared_atom(bond)`
    - `unshared_atom(bond)`
    - `direction()`
    - `as_expr()`
    - `get_embedding_angle_cos(i, j, k, l)`
    - `get_positioning_cos(i, j, k, l)`
    - `angle_cos(other, require_int=False)`
    - `dot(other)`
    - `cross(other)`
    - `polar_components(embedding_atoms, vector)`
  - **class `BondNormal`**
    - `__init__(a, b, norm=1, embedding=None)`
    - `direction()`
    - `embedding_basis(embedding)`
    - `polar_components(embedding, self)`
    - `angle_cos(other, require_int=True)`
    - `as_expr()`
    - `cross(other)`
    - `dot(other)`
  - **class `BondVectorSum`**
    - `__init__(*terms)`
    - `simplify()`
    - `as_expr()`
    - `distribute(f)`
    - `dot(other)`
  - **class `InternalJacobianDisplacements`**
    - `symbolic_m(i)`
    - `symbolic_r(i, j)`
    - `symbolic_a(i, j, k)`
    - `symbolic_t(i, j, k, l)`
    - `symbolic_y(i, j, k, l)`
    - `lam(i, j, k)`
    - `dr(i, j)`
    - `da(i, j, k)`
    - `dt(i, j, k, l)`
    - `dy(i, j, k, l)`
    - `displacement_vectors(inds, coord_type=None)`

### `AnalyticModelConstructors.py`
  - **class `AnalyticPotentialConstructor`** (AnalyticModelBase)
    > Provides a set of symbolic potentials for use in models
    > *(truncated — see stub for full docstring)*
    - `morse(i, j, De=None, a=None, re=None, eq=None, w=None, wx=None)` — Returns a fully symbolic form of a Morse potential
    - `calc_morse(De, a, r, re)`
    - `harm(k, x, x_e)`
    - `harmonic(*args, k=None, eq=None, qe=None)` — Returns a fully symbolic form of a Morse potential
    - `lin(k, x, x_e)`
    - `linear(*args, k=1, eq=None, xe=None)` — Returns a fully symbolic form of a linear function
    - `pow(k, x, x_e, n)`
    - `power(*args, k=1, eq=None, n=None, xe=None)` — Returns a fully symbolic form of a linear function
    - `cos(*args, eq=None, qe=None)` — Returns a fully symbolic form of cos
    - `sin(*args, eq=None, qe=None)` — Returns a fully symbolic form of sin
    - `multiwell(*args, turning_points=None, origin=None, eq=None, minimum=0, depth=None)` — :param args:
  - **class `GeometricFunction`**
    - `__init__(coords, coord_funcs, val_func)`
    - `position_function(i)`
    - `normal_function(i, j, k)`
    - `mass_function(i)`
    - `distance_function(i, j)`
    - `angle_function(i, j, k)`
    - `dihedral_function(i, j, k, l)`
    - `book_function(i, j, k, l)`
    - `create_symbol_function(sym)`
    - `sorted_vars(vars)`
    - `from_expr(expr)`
  - **class `AnalyticKineticEnergyConstructor`** (AnalyticModelBase)
    > Provides G and V' elements from Frederick and Woywood
    > *(truncated — see stub for full docstring)*
    - `kinetic_exprs(inds1, inds2, coord_types=None, target_symbols=None)`
    - `kinetic_exprs_direct(inds1, inds2, coord_types=None, return_vp=True, target_symbols=None)` — Evaluated using the simple expressions in Table 1 from Frederick and Woywood
    - `g(inds1, inds2, coord_types=None, target_symbols=None, return_function=False, method='lookup')`
    - `vp(inds1, inds2, coord_types=None, target_symbols=None, return_function=False, method='lookup')`
    - `g_matrix(coord_specs, return_function=False, return_matrix=True, method='lookup')`
  - **class `AnalyticModel`**
    > Provides a symbolic representation of an analytically evaluatable Hamiltonian
    > which can be used to get derived expressions to evaluate.
    > Supplies methods to automatically run DVR and VPT calculations from the model
    > specifications as well.
    > *(truncated — see stub for full docstring)*
    - `__init__(coordinates, potential, dipole=None, values=None, rotation=None)`
    - `from_potential(potential, dipole=None, values=None, rotation=None)`
    - `internal_coordinates()`
    - `constants()`
    - `normal_modes(dimensionless=True)`
    - `to_normal_modes(dimensionless=True)`
    - `get_VPT_expansions(order=2, expansion_order=None, include_potential=None, include_gmatrix=None, include_pseudopotential=None, evaluate=True)`
    - `run_VPT(order=2, states=2, return_analyzer=True, expansion_order=None, include_potential=None, include_gmatrix=None, include_pseudopotential=None, atoms=None, coords=None, **kwargs)`
    - **class `SympyExpr`**
      - `__init__(expr, core, ndim)`
      - `lambdify(vars, expr, ndim, lambdify_arrays=False)`
    - `prep_lambda_expr(base_coords, expr, dummify=False, rewrite_trig=True)`
    - `lambdify(coord_vec, expr, coordinate_transform=None, mode=None, dummify=False, rewrite_trig=True, lambdify_arrays=False)`
    - `wrap_function(expr, transform_coordinates=True, mode=None)`
    - `expand_potential(order, lambdify=True, evaluate=True, contract=True)`
    - `get_DVR_parameters(expansion_order=None, lambdify=True, evaluate='constants')`
    - `setup_DVR(domain=None, divs=None, use_normal_modes=False, expansion_order=None, potential_function=None, **params)`
    - `evaluate(expr, mode='all', numericize=False)`
    - `parse_symbol(sym)`
    - `jacobian(order=0, evaluate=False, lambdify=False)`
    - `jacobian_inverse(order=0, evaluate=False)`
    - `g(order=0, evaluate=False, lambdify=False)`
    - `v(order=2, evaluate=False, lambdify=False)`
    - `vp(order=0, evaluate=False, lambdify=False)`
    - `mu(order=1, evaluate=False, lambdify=False)`
    - **class `NamespaceContext`**
      - `__init__(context=None)`
      - `insert_vars()`
      - `prune_vars()`
    - `sym(base, *args)`
    - `m(i)`
    - `r(i, j)`
    - `a(i, j, k)`
    - `t(i, j, k, l)`
    - `y(i, j, k, l)`
    - `mass(atom)`
    - `molecular_potential(mol)`
    - `molecular_dipole(mol)`
    - `molecular_gmatrix(mol)`
  - **class `MolecularModel`** (AnalyticModel)
    - `__init__(mol, coords, potential, dipole=None, values=None, rotation=None)`
    - `potential()`
    - `gmatrix()`
    - `vprime()`
    - `dipole()`
    - `setup_AIMD(timestep=0.5, initial_energies=None, initial_displacements=None, displaced_coords=None, track_kinetic_energy=False, track_velocities=False, **aimd_opts)`
    - `setup_DGB(centers, *, masses=None, modes='normal', transformations=None, alphas='auto', cartesians=None, potential_function=None, dipole_function=None, optimize_centers=None, quadrature_degree=None, expansion_degree=None, pairwise_potential_functions=None, internals=False, logger=True, **opts)`
  - **class `MolecularModelFunction`**
    - `__init__(deriv_evaluator, mol)`
    - `evaluate(carts, deriv_order=None, internals=False, which=None, sel=None, axes=None, derivs=None)`
  - **class `MolecularModelPotentialFunction`** (MolecularModelFunction)
    - `__init__(model, mol)`
  - **class `MolecularModelDipoleFunction`** (MolecularModelFunction)
    - `__init__(model, mol)`
    - `evaluate(carts, deriv_order=None, internals=False, which=None, sel=None, axes=None)` — This has the added complication of needing to dispatch over the axes...
  - **class `MolecularModelGMatrixFunction`** (MolecularModelFunction)
    - `__init__(model, mol)`
  - **class `MolecularModelVPrimeFunction`** (MolecularModelFunction)
    - `__init__(model, mol)`

### `Helpers.py`
  - **class `SympyShim`**
    > Provides a loader that can either load sympy when requested
    > or throw an error if it can't be loaded
  - **class `SymbolicCaller`**
    > Delegates to `__call__` through `__getitem__` for symbolics
    - `__init__(sym)`
  - **class `AnalyticModelBase`**
    > Provides a base class for analytic models
    - `get_numeric_types()`
    - `take_derivs(expr, vars)` — Takes derivatives of `expr` with respect to `vars` even if `expr` is an array
    - `eval_exprs(expr, subs)` — Evaluates `expr` with the given substitutions
    - `symbol_list(names, instance=None)` — Gets a list of symbols for `names` with a given instance number
    - `symbolic_x(i)` — Provides a symbolic representation of a position
    - `symbolic_n(i, j, k)` — Provides a symbolic representation of a normal to a plane
    - `symbolic_m(i)` — Provides a symbolic representation of a mass
    - `symbol(base, *args, **kwargs)`
    - `symbolic_r(i, j)` — Provides a symbolic representation of a bond length
    - `symbolic_a(i, j, k)` — Provides a symbolic representation of a bond angle
    - `symbolic_t(i, j, k, l)` — Provides a symbolic representation of a dihedral angle
    - `symbolic_y(i, j, k, l)` — Provides a symbolic representation of a book angle
    - `infer_coord_type(inds)`
    - `var(*args, coord_type=None)`
    - `reindex_symbol(symbol, mapping, target_symbols=None)` — Changes the indices on symbols using the given mapping
    - `lam(i, j, k)` — Provides the `lambda` expression from Frederick and Woywood
    - `is_identity(A)`
    - `transpose(A)`
    - `dot(a, b, axes=None)`
    - `contract(a, axes)`
    - `transform_coordinates(rotation, coord_vec=None, coord_name_fmt='q{id}[{num}]')`