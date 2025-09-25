# Repository Guidelines

## Project Structure & Module Organization
- `src/`: Fortran 90 modules for ERSEM (e.g., `pelagic_base.F90`, `benthic_*`, `oxygen.F90`). Built as library `fabm_models_ersem` via CMake.
- `docs/`: Sphinx documentation (`make html`), with developer setup under `docs/source/developers/`.
- `github-actions/`: CI build and system-test scripts for PyFABM, GOTM+FABM, and FABM0d.
- `testcases/`: Example FABM/ERSEM YAML setups and utilities.
- `conda.recipe/`: Recipe to build ERSEM together with FABM, GOTM, and PyFABM.

## Build, Test, and Development Commands
- Configure + build (standalone/developer layout):
  - Clone dependencies (see docs) then: `mkdir build && cd build && cmake .. -DFABM_ERSEM_BASE=.. [-DFABM_BASE=../fabm -DGOTM_BASE=../gotm] && make -j`.
- Conda build (bundles deps): `conda build conda.recipe`.
- Docs: `cd docs && pip install -r requirements.txt && make html`.
- Tests (run what CI runs):
  - PyFABM: `pytest github-actions/pyfabm-ersem`
  - GOTM+FABM: `pytest github-actions/gotm-fabm-ersem`
  - FABM0d: `pytest github-actions/fabm0d-gotm-ersem`

## Coding Style & Naming Conventions
- Language: Fortran 90 free-form (`.F90`). Prefer lower_snake_case for modules, files, and variables; types often `type_ersem_*`; IDs `id_*`.
- Indentation: consistent small indent (2–3 spaces); align continuations thoughtfully.
- Lines: keep under ~100 chars; comment with `!` and write clear procedure headers.
- CMake options: `-DERSEM_USE_IRON=ON` to include iron; use `-DCMAKE_Fortran_FLAGS="-fcheck=all"` for debug (as in CI).

## Testing Guidelines
- Framework: pytest-based system tests in `github-actions/*`. Ensure tests pass locally before PRs.
- Add or update expected outputs where applicable (see `regen_expected_results.py`).
- Name new tests to match the area they cover and keep inputs under version control.

## Commit & Pull Request Guidelines
- Commits: short, imperative subject; reference issues/PRs (e.g., `Fix oxygen parameter (#135)`). Group related changes.
- PRs: include description, rationale, linked issues, and any result plots or logs if behavior changes. Update docs/testcases when interfaces or parameters change. CI must pass.

## Configuration Tips
- Building requires a Fortran compiler, CMake, and NetCDF-Fortran. The `conda.recipe/build.sh` shows exact steps used in automation and is a good reference.
