# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

ERSEM (European Regional Seas Ecosystem Model) is a marine biogeochemical and ecosystem model describing carbon, nitrogen, phosphorus, silicon, oxygen, and iron cycling through pelagic and benthic ecosystems. The codebase is Fortran 90, integrating with the FABM (Framework for Aquatic Biogeochemical Models) framework.

This is a UTokyo development version optimized for estuarine and coastal systems. Upstream: https://github.com/pmlmodelling/ersem

## Build Commands

### PyFABM-ERSEM (Python interface)
```bash
git clone https://github.com/fabm-model/fabm.git
cd fabm && python -m pip install .
```

### GOTM-FABM-ERSEM (1D water column)

**Generic build:**
```bash
mkdir build && cd build
cmake ../gotm -DFABM_BASE=../fabm -DFABM_ERSEM_BASE=../ersem
make install -j $(nproc)
```

**Tokyo Bay Test Case (TB-GOTM) Build:**
```bash
# Clean previous build and configure
rm -rf ~/build && mkdir ~/build && cd ~/build

cmake ~/Github/gotm/code \
  -DFABM_BASE=~/Github/fabm \
  -DFABM_ERSEM_BASE=~/Github/ersem \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_INSTALL_PREFIX=~/local/fabm-ifx/gotm

# Build and install (use 8 cores)
make install -j8

# The gotm executable will be at: ~/local/fabm-ifx/gotm/bin/gotm
# Symlink already exists at: ~/Github/TB-GOTM/ersem/run/gotm
```

**Run Tokyo Bay simulation:**
```bash
cd ~/Github/TB-GOTM/ersem/run
./gotm --ignore_unknown_config
```

**Plot results:**
```bash
conda activate xgotm
python plot_wq_contourf.py temp salt O2_o pel_sulfur_H2S
```

### FABM0d (0D box model)
```bash
cmake ../fabm/src/drivers/0d -DGOTM_BASE=../gotm/code -DFABM_ERSEM_BASE=../ersem
make install -j $(nproc)
```

### CMake Options
- `-DERSEM_USE_IRON=ON` - Enable iron cycling (adds `-DIRON` preprocessor flag)
- `-DCMAKE_Fortran_COMPILER=ifx` or `gfortran`
- `-DCMAKE_Fortran_FLAGS="-fcheck=all"` - Enable runtime checks for debugging

## Testing

Tests run via pytest in the `github-actions/` directory:
```bash
pytest github-actions/pyfabm-ersem           # PyFABM tests
pytest github-actions/gotm-fabm-ersem        # GOTM tests
pytest github-actions/fabm0d-gotm-ersem      # FABM0d tests
```

## Code Architecture

### FABM Integration Pattern

All ERSEM modules follow the FABM model pattern:
1. Each module defines a type extending a base class (e.g., `type_ersem_primary_producer` extends `type_ersem_pelagic_base`)
2. Types implement `initialize` (parameter/variable registration) and `do` (rate calculations) procedures
3. The model factory in `ersem_model_library.F90` registers all models and maps YAML names to Fortran types

### Module Hierarchy

```
ersem_model_library.F90  -- Factory: maps "ersem/xyz" names to types
    |
    +-- shared.F90           -- Constants, flags (use_iron), standard variables
    |
    +-- pelagic_base.F90     -- Base for pelagic state variables (C, N, P, Si, Fe, Chl)
    |   +-- primary_producer.F90   -- P1-P4 phytoplankton
    |   +-- microzooplankton.F90   -- Heterotrophic nanoflagellates
    |   +-- mesozooplankton.F90    -- Copepods, etc.
    |   +-- bacteria.F90           -- Standard bacterial remineralization
    |   +-- bacteria_docdyn.F90    -- Dynamic DOC bacteria
    |
    +-- benthic_base.F90     -- Base for benthic state variables
        +-- benthic_column.F90              -- Sediment layering/diffusion
        +-- benthic_bacteria.F90            -- Sediment bacteria
        +-- benthic_fauna.F90               -- Deposit/filter feeders
        +-- benthic_nitrogen_cycle.F90      -- Nitrification/denitrification
```

### Key Files

- `src/ersem_model_library.F90` - Model factory; add new models to the `select case` block
- `src/shared.F90` - Shared constants (`CMass`, `ZeroX`), `use_iron` flag, aggregate standard variables
- `src/primary_producer.F90` - Primary production with optional silicate, calcification, iron
- `testcases/*.yaml` - Example configurations showing model coupling

### Adding a New Model

1. Create `src/new_model.F90` following existing module patterns
2. Add `use ersem_new_model` to `ersem_model_library.F90`
3. Add case to the factory: `case ('new_model'); allocate(type_ersem_new_model::model)`
4. Add source file to `src/CMakeLists.txt`
5. Create YAML configuration entry in testcases

## Coding Conventions

- Fortran 90 free-form (`.F90`), lower_snake_case for modules/variables
- Types: `type_ersem_*`; variable IDs: `id_*`
- Keep lines under ~100 chars
- Use `!` for comments with clear procedure headers

## Configuration

Model setup via YAML files (examples in `testcases/`):
- Instance names map to `ersem/<model_name>`
- Parameters, initial conditions, and variable coupling defined per instance
- Template: `fabm-ersem.yaml.template`

## Documentation

- Build docs: `cd docs && pip install -r requirements.txt && make html`
- Online: https://ersem.readthedocs.io/
- Reference: Butenschön et al. (2016), GMD, doi:10.5194/gmd-9-1293-2016

## Language

- **Chat (conversation)**: Always respond in Japanese (日本語で応答すること)
- **Documentation files** (README.md, etc.): English
- **Code comments**: English
- **Commit messages**: English
