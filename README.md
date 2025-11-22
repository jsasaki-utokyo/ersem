![Tests](https://img.shields.io/github/actions/workflow/status/pmlmodelling/ersem/ersem.yml?label=tests&style=flat-square)
[![Documentation Status](https://readthedocs.org/projects/ersem/badge/?version=latest)](https://ersem.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/302390544.svg)](https://zenodo.org/badge/latestdoi/302390544)

# ERSEM

![ERSEM diagram](docs/images/ERSEM.png)

## Overview

[ERSEM](http://ersem.com) (European Regional Seas Ecosystem Model) is a marine biogeochemical and ecosystem model. It describes the cycling of carbon, nitrogen, phosphorus, silicon, oxygen and iron through the lower trophic level pelagic and benthic ecosystems.

**Key Features:**
- Multi-element biogeochemical cycling (C, N, P, Si, O, Fe)
- Pelagic and benthic ecosystem components
- Integration with FABM (Framework for Aquatic Biogeochemical Models)
- Multiple deployment options (0D, 1D, 3D)

## Quick Start

### Conda Installation (Recommended for Users)

```bash
conda install -c conda-forge ersem
```

For detailed tutorials, see: https://ersem.readthedocs.io/en/latest/tutorials/index.html

### Python Interface (PyFABM-ERSEM)

```bash
# Clone dependencies
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/pmlmodelling/ersem.git

# Build and install
cd fabm
python -m pip install .
```

## Installation from Source

For developers and advanced users who need to build ERSEM from source.

### Prerequisites

- Fortran compiler (Intel ifort or gfortran)
- CMake (>= 3.0)
- NetCDF libraries
- Git

### 1. Clone Repositories

The following assumes repositories are cloned in `~/Github/`:

```bash
# Create directory structure
mkdir -p ~/Github && cd ~/Github

# Clone GOTM with submodules (v6.0 branch)
git clone --recurse-submodules -b v6.0 https://github.com/gotm-model/code.git gotm

# Clone FABM (main branch)
git clone https://github.com/fabm-model/fabm.git

# Clone ERSEM (main branch)
git clone https://github.com/pmlmodelling/ersem.git

# Optional: Create topic branches for development
cd fabm && git checkout -b topic && cd ..
cd ersem && git checkout -b topic && cd ..
cd gotm && git checkout -b topic && cd ..
```

### 2. Build Options

Choose the configuration that matches your needs:

#### FABM0D (Box Model / Aquarium Setup)

For 0-dimensional box model simulations. Build files are created in `~/build/`, executables installed in `~/local/fabm-{compiler}/0d/bin/`.

**Installation Script** (`install_ersem_fabm0d.sh`):

```bash
#!/bin/bash

# Compiler selection: ifort or gfortran
COMPILER=ifort
# COMPILER=gfortran

# Number of cores for compilation
NCORES=$(nproc)

# Set paths
GITHUB_DIR=~/Github
BUILD_DIR=~/build
INSTALL_BASE=~/local/fabm-${COMPILER}

# Create build directory
mkdir -p ${BUILD_DIR}/fabm-${COMPILER}/0d
cd ${BUILD_DIR}/fabm-${COMPILER}/0d

# Configure with CMake
cmake ${GITHUB_DIR}/fabm/src/drivers/0d \
  -DFABM_HOST=0d \
  -DFABM_BASE=${GITHUB_DIR}/fabm \
  -DFABM_ERSEM_BASE=${GITHUB_DIR}/ersem \
  -DCMAKE_Fortran_COMPILER=${COMPILER} \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_BASE}/0d

# Build and install
make -j${NCORES}
make install
```

Make executable and run:
```bash
chmod +x install_ersem_fabm0d.sh
./install_ersem_fabm0d.sh
```

#### GOTM-ERSEM (1D Water Column Model)

For 1-dimensional water column simulations with turbulence. Build files are created in `~/build/`, executables installed in `~/local/fabm-{compiler}/gotm/bin/`.

**Installation Script** (`install_ersem_gotm.sh`):

```bash
#!/bin/bash

# Compiler selection: ifort or gfortran
COMPILER=ifort
# COMPILER=gfortran

# Number of cores for compilation
NCORES=$(nproc)

# Set paths
GITHUB_DIR=~/Github
BUILD_DIR=~/build
INSTALL_BASE=~/local/fabm-${COMPILER}

# Create build directory
mkdir -p ${BUILD_DIR}/fabm-${COMPILER}/gotm
cd ${BUILD_DIR}/fabm-${COMPILER}/gotm

# Configure with CMake
cmake ${GITHUB_DIR}/gotm \
  -DFABM_BASE=${GITHUB_DIR}/fabm \
  -DFABM_ERSEM_BASE=${GITHUB_DIR}/ersem \
  -DCMAKE_Fortran_COMPILER=${COMPILER} \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_BASE}/gotm

# Build and install
make -j${NCORES}
make install
```

Make executable and run:
```bash
chmod +x install_ersem_gotm.sh
./install_ersem_gotm.sh
```

#### Build Notes

- Add `-DIRON` CMake flag to enable iron cycling
- For custom NetCDF locations, add: `-DNETCDF_CONFIG=/path/to/nc-config`
- Default installation creates executables: `fabm0d` and `gotm` respectively

## Testing

ERSEM includes comprehensive automated tests run via GitHub Actions:

```bash
# PyFABM tests
pytest github-actions/pyfabm-ersem

# GOTM-FABM-ERSEM tests
pytest github-actions/gotm-fabm-ersem/test_state_variables.py
pytest github-actions/gotm-fabm-ersem/test_gotm.py

# FABM0d tests
pytest github-actions/fabm0d-gotm-ersem
```

Test configurations are available in the `testcases/` directory.

**View test results:** https://github.com/pmlmodelling/ersem/actions

## Documentation & Resources

- **User Documentation:** https://ersem.readthedocs.io/en/latest/
- **Tutorials:** https://ersem.readthedocs.io/en/latest/tutorials/index.html
- **Developer Guide:** https://ersem.readthedocs.io/en/latest/developers/index.html
- **Test Cases:** `testcases/` directory contains various YAML configurations
- **Model Configuration:** `fabm-ersem.yaml.template` for setup examples

## Support

We strongly encourage everyone using the ERSEM code to register as a user by filling a [short registration form](https://forms.office.com/r/X0iXv8AvTC).

This helps us:
- Understand who is using ERSEM and for what applications
- Provide better support to the user community
- Keep you informed about latest model developments and news

## How to Cite

To refer to ERSEM in publications, please cite:

Butenschön, M., Clark, J., Aldridge, J.N., Allen, J.I., Artioli, Y., Blackford, J., Bruggeman, J., Cazenave, P., Ciavatta, S., Kay, S., Lessin, G., van Leeuwen, S., van der Molen, J., de Mora, L., Polimene, L., Sailley, S., Stephens, N., Torres, R. (2016). ERSEM 15.06: a generic model for marine biogeochemistry and the ecosystem dynamics of the lower trophic levels. Geoscientific Model Development, 9(4), 1293–1339. doi: [10.5194/gmd-9-1293-2016](https://doi.org/10.5194/gmd-9-1293-2016).

To refer specifically to the ERSEM source code, you may use its Zenodo DOI (badge at top of page).

## Acknowledgements & License

- **Acknowledgements:** https://ersem.readthedocs.io/en/latest/acknowledgements.html
- **License:** https://ersem.readthedocs.io/en/latest/license.html
