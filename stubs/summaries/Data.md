### `KEData.py` — Provides a class for handling a compiled set of atomic data
  - **class `KEDataHandler`** (DataHandler)
    > A DataHandler that's built for use with the G-matrix and V' terms
    > from the 1999 Frederick and Woywood paper
    > Usually used through the `KETermData` object.
    - `__init__()`
    - `find_expressions(k, return_permutation=False)`

### `PotentialRegistry.py`
  - **class `PotentialReleaseZIPManager`** (ReleaseZIPManager)
  - **class `PotentialRegistryAPI`** (GitHubReleaseManager)
    - `__init__(token=None, request_delay_time=None, release_manager=None, **opts)`
    - `list_repos(owner=None)`
    - `list_potentials()`
    - `list_releases(repo_or_owner, repo=None)`
    - `latest_release(repo_or_owner, repo=None)`
    - `get_release_list(repo_or_owner, repo=None, update=None)`
    - `get_potential(name, update=None, version=None)`

### `ScanManager.py` — `ScanManager` is the single entry point for both halves of a scan:
- `shape_scan_iterator(base_shape, zigzag=False)`
- `scan_iterator(domains, expand_domains=True, index_iterator=None, zigzag=False)`
- `molecule_scan_iterator(mol, geometries, index_iterator=None, zigzag=False, return_molecules=True)`
- `molecule_scan_geometries_iterator(mol, domains, which, return_values=False, return_molecules=True, index_iterator=None, zigzag=False, coordinate_generator=None, **etc)`
- `molecule_displaced_geometries_iterator(mol, displacement_positions, which, return_molecules=True, return_values=False, index_iterator=None, zigzag=False, coordinate_generator=None, **etc)`
- `molecule_atom_position_scan_iterator(mol, atom_indices, domains, which=None, embedding=None, return_molecules=True, **iterator_options)`
  - **class `ScanManager`**
    - `__init__(output_directory, scan_id=None, job_prefix='scan', index_format='03d')`
    - `scan_dir()` — Directory jobs are written to / read from.
    - `scan_info_file()`
    - `default_job_builder(mol, *, job_type, commands=None, **etc)`
    - `generate(scan_iterator, job_builder=None, coord_labels=None, extra_info=None, overwrite=False, append=False, job_prefix=None, job_file_ext=None, job_type=None, **job_kwargs)`
    - `load_scan_info()` — Loads this scan's `scan_info.json` manifest.
    - `default_output_file_generator(input_file)` — Default `output_file_generator`: swaps the input job file's extension
    - `load_molecules(output_file_generator=None, molecule_loader=None, scan_info=None, skip_missing=True)` — Rebuilds a `Molecule` for every completed step of the scan.
    - `parse(molecular_property_extractor, output_file_generator=None, molecule_loader=None, scan_info=None, skip_missing=True, fill_value=np.nan)` — Rebuilds the `Molecule` for every completed scan step, runs

### `Surfaces.py` — Provides concrete tools for dealing with two of the most useful types of surfaces we have
  - **class `DipoleSurface`** (MultiSurface)
    > Provides a unified interface to working with dipole surfaces.
    > Currently basically no fancier than a regular surface (although with convenient loading functions), but dipole-specific
    > stuff could come
    - `__init__(mu_x, mu_y, mu_z)`
    - `center()`
    - `ref()`
    - `expansion_tensors()`
    - `get_log_values(log_file, keys=('StandardCartesianCoordinates', 'DipoleMoments'))`
    - `from_log_file(log_file, coord_transf, keys=('StandardCartesianCoordinates', 'DipoleMoments'), tol=0.001, **opts)` — Loads dipoles from a Gaussian log file and builds a dipole surface by interpolating.
    - `from_fchk_file(fchk_file, **opts)` — Loads dipoles from a Gaussian formatted checkpoint file and builds a dipole surface via a linear ap…
    - `from_derivatives(expansion, center=None, **opts)`
    - `from_mol(mol, expansion=None, center=None, transforms=None, use_internals=True, **opts)`
  - **class `PotentialSurface`** (Surface)
    > A potential surface structure to go along with the DipoleSurface.
    > Provides convenient access to potential data + a unified interface to things like energy minimization
    - `get_log_values(log_file, keys=('StandardCartesianCoordinates', 'ScanEnergies'))`
    - `from_log_file(log_file, coord_transf, keys=('StandardCartesianCoordinates', 'ScanEnergies'), tol=0.001, **opts)` — Loads dipoles from a Gaussian log file and builds a potential surface by interpolating.
    - `from_fchk_file(fchk_file, ref=None, **opts)` — Loads potential from a Gaussian formatted checkpoint file and builds a potential surface via a quar…
    - `from_mol(mol, expansion=None, center=None, transforms=None, transformed_derivatives=False, use_internals=True, **opts)`
    - `from_derivatives(expansion, center=None, ref=None, transforms=None, transformed_derivatives=False, **opts)`