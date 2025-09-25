# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

ERSEM (European Regional Seas Ecosystem Model) is a marine biogeochemical and ecosystem model describing carbon, nitrogen, phosphorus, silicon, oxygen, and iron cycling through lower trophic level pelagic and benthic ecosystems. The codebase is primarily written in Fortran 90 and integrates with the FABM (Framework for Aquatic Biogeochemical Models) framework.

## Build System

ERSEM uses CMake for building and requires integration with FABM. The main build configurations:

### PyFABM-ERSEM (Python interface)
```bash
# Clone dependencies
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/pmlmodelling/ersem.git

# Build and install
cd fabm
python -m pip install .
```

### GOTM-FABM-ERSEM (with General Ocean Turbulence Model)
```bash
# Clone all components
git clone https://github.com/pmlmodelling/ersem.git
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/gotm-model/code.git gotm

# Build with CMake
mkdir build && cd build
cmake ../gotm -DFABM_BASE=../fabm -DFABM_ERSEM_BASE=../ersem
make install -j $(nproc)
```

### FABM0d (0-dimensional aquarium setup)
Similar to GOTM build but uses the 0d configuration for box model simulations.

## Testing

Tests are run via GitHub Actions (`.github/workflows/ersem.yml`):

- **PyFABM tests**: `pytest github-actions/pyfabm-ersem`
- **GOTM tests**: `pytest github-actions/gotm-fabm-ersem/test_state_variables.py github-actions/gotm-fabm-ersem/test_gotm.py`
- **FABM0d tests**: `pytest github-actions/fabm0d-gotm-ersem`

Test configurations are in `testcases/` directory with various YAML setups.

## Code Architecture

### Core Modules (`src/`)

The model is organized into functional modules:

**Pelagic Components:**
- `primary_producer.F90` - Phytoplankton dynamics
- `microzooplankton.F90`, `mesozooplankton.F90` - Zooplankton
- `bacteria.F90`, `bacteria_docdyn.F90` - Bacterial processes

**Benthic Components:**
- `benthic_base.F90` - Base benthic functionality
- `benthic_bacteria.F90`, `benthic_fauna.F90` - Benthic organisms
- `benthic_column_*.F90` - Sediment-water column interactions

**Biogeochemical Processes:**
- `carbonate.F90`, `calcification.F90` - Carbon system
- `oxygen.F90` - Oxygen dynamics
- `nitrification.F90`, `nitrous_oxide.F90` - Nitrogen cycle
- `light.F90`, `light_iop*.F90` - Light attenuation

**Framework Integration:**
- `ersem_model_library.F90` - Main FABM interface
- `shared.F90` - Shared utilities
- Iron support optional via `-DIRON` flag

### Configuration

Model setup through YAML files (examples in `testcases/`):
- Define state variables, parameters, and coupling
- Different configurations for various complexity levels
- Template: `fabm-ersem.yaml.template`

## Documentation

- User docs: https://ersem.readthedocs.io/
- Build documentation: `cd docs && make html`
- Main publication: Butenschön et al. (2016), GMD, doi:10.5194/gmd-9-1293-2016