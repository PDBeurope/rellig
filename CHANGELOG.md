# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-07-27

### Changed

- **Breaking:** `Cofactors`, `Reactants`, and `Drugs` pipeline classes now take
  explicit constructor arguments (`ligand_cif`, `ligand_type`, `out_dir`, `logger`,
  etc.) instead of an `argparse.Namespace` object, so they can be used
  programmatically without going through the CLI.
- Cofactor reference data (RDKit templates, cofactor details, and EC numbers) is
  now loaded via cached module-level functions (`init_rdkit_templates`,
  `get_cofactor_details`, `get_cofactor_ec` in `pdberellig.helpers.utils`) instead
  of being loaded once per `Cofactors` instance.
- ChEBI structure file is now parsed using its native/perceived data types
  instead of forcing every column to `str`, fixing comparisons against
  `status_id` and `default_structure`.

### Added

- Pytest suite covering the cofactors, drugs, and reactants pipelines.
- GitHub Actions workflow to run the test suite on pushes/PRs.

### Fixed

- `Drugs.get_drugbank_targets` now returns an empty `DataFrame` (instead of
  `None`) when no DrugBank targets category is present in the CIF file,
  preventing downstream failures.
