# Changelog

## 1.6.0

Changes since 1.5.0. This release rebuilds the electronic-structure backend and
VPT evaluation engine. Use with McUtils 1.8.0.

### Added

- Reusable electronic-structure backends for ML potentials, xTB and PySCF; charge, dipole, polarizability, electrostatic-potential and charge/spin-density descriptors.
- Lazy DAG and indexed analytic VPT evaluation for higher-order calculations, with bounded caches, parallel evaluation and improved degeneracy handling.
- VPT-to-2D-IR response preparation and caching, including imports from classic and analytic VPT logs.
- Molecular readouts and exportable normal-mode animations with mode selection and amplitude controls.
- Pairing-based Franck–Condon overlap evaluation, optional DVR pseudopotentials, Weyl eigenvalue bounds and additional DGB generalized-eigenproblem solvers.
- Experimental conformer alignment, reaction-coordinate searches, NEB/string path optimization and transition-state refinement.

### Changed

- DGB defaults to `lobpcg`; select `mode="similarity"` to retain the previous solver.
- Electronic-structure evaluators delegate to `Psience.ElectronicStructure` backends. Existing evaluator names remain; custom model setup and calculation integrations should use the evaluator's `.backend`.
- `MolecularHamiltonian.get_VPT_expansions(n)` interprets integer `n` as perturbation order, using potential derivatives through `n + 2`; dictionaries still specify derivative orders explicitly.
- RMSD bias defaults to unweighted per-atom RMSD and deposits hills after completed integration steps.

### Fixed

- Franck–Condon rotations, high-quanta overlap counting and ezFCF mass weighting.
- DGB integrals and kinetic terms for rotated and momentum-bearing Gaussian basis functions, and molecular Hamiltonian derivative terms.
- RMSD metadynamics forces and walker histories, Langevin noise scaling, state-space quantum bounds, InChI round trips and ring bond rendering.

[Full diff](https://github.com/McCoyGroup/Psience/compare/v1.5.0...d727b62ce87d7cc549e3a80172243f74f2bfd046)
