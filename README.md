[![DOI](https://zenodo.org/badge/302390544.svg)](https://zenodo.org/badge/latestdoi/302390544)

# ERSEM - Japanese Waters Configuration

![ERSEM diagram](docs/images/ERSEM.png)

## Overview

This repository contains a development version of ERSEM (European Regional Seas Ecosystem Model) configured for Japanese waters. ERSEM is a marine biogeochemical and ecosystem model describing the cycling of carbon, nitrogen, phosphorus, silicon, oxygen and iron through the lower trophic level pelagic and benthic ecosystems.

**Key Features:**
- Multi-element biogeochemical cycling (C, N, P, Si, O, Fe)
- Pelagic and benthic ecosystem components
- Integration with FABM (Framework for Aquatic Biogeochemical Models)
- Multiple deployment options (0D box model, 1D water column)

**Note:** This is a development branch for Japanese waters. The upstream version is maintained at https://github.com/pmlmodelling/ersem (master branch).

## Installation from Source

### Prerequisites

- Fortran compiler (Intel ifx, ifort, or gfortran)
- CMake (>= 3.0)
- NetCDF libraries
- Git

### 1. Clone Repositories

The following assumes repositories are cloned in `~/Github/`:

```bash
# Create directory structure
mkdir -p ~/Github && cd ~/Github

# Clone GOTM with submodules
git clone --recurse-submodules https://github.com/gotm-model/code.git gotm/code

# Clone FABM
git clone https://github.com/fabm-model/fabm.git

# Clone this ERSEM repository
git clone https://github.com/jsasaki-utokyo/ersem.git

# Optional: Create topic branches for development
cd fabm && git checkout -b topic && cd ..
cd ersem && git checkout -b topic && cd ..
cd gotm/code && git checkout -b topic && cd ../..
```

### 2. Build Options

Choose the configuration that matches your needs:

#### FABM0D (Box Model / Aquarium Setup)

For 0-dimensional box model simulations. Executables installed in `~/local/fabm-{compiler}/0d/bin/`.

**Installation Script** (`install_ersem_fabm0d.sh`):

```bash
#!/bin/bash
# install_ersem_fabm0d.sh

CMAKE_FLAG=""
while getopts "f:" flag; do
    case "${flag}" in
        f) CMAKE_FLAG=${OPTARG};;
    esac
done

CPU="$(nproc)"
FABM_HOST=0d

# Specify ifx, ifort, or gfortran
CMAKE_Fortran_COMPILER=ifx

# You may want to use -O3 for better performance but check whether
# the change in results are acceptable.
# Default of -O2 will be used with commenting out the next line.
# FFLAGS=-O3

GOTM_BASE=~/Github/gotm/code
FABM_ERSEM_BASE=~/Github/ersem
FABM_BASE=~/Github/fabm/src/drivers/$FABM_HOST
# CMAKE_BUILD_TYPE=Debug

CMAKE_INSTALL_PREFIX=~/local/fabm-${CMAKE_Fortran_COMPILER}/${FABM_HOST}

echo $CMAKE_FLAG

rm -rf ~/build  # Delete if exits
mkdir ~/build && cd ~/build
cmake $FABM_BASE \
  -DGOTM_BASE=$GOTM_BASE \
  -DFABM_ERSEM_BASE=$FABM_ERSEM_BASE \
  -DFABM_HOST=$FABM_HOST \
  -DCMAKE_Fortran_COMPILER=$CMAKE_Fortran_COMPILER \
  -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
  -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
  $CMAKE_FLAG
make install -j $CPU
```

Make executable and run:
```bash
chmod +x install_ersem_fabm0d.sh
./install_ersem_fabm0d.sh
# Or with additional CMake flags:
./install_ersem_fabm0d.sh -f "-DIRON=ON"
```

#### GOTM-ERSEM (1D Water Column Model)

For 1-dimensional water column simulations with turbulence. Executables installed in `~/local/fabm-{compiler}/gotm/bin/`.

**Installation Script** (`install_ersem_gotm.sh`):

```bash
#!/bin/bash
# install_ersem_gotm.sh

CMAKE_FLAG=""
while getopts ":f:" flag; do
    case "${flag}" in
        f) CMAKE_FLAG=${OPTARG};;
    esac
done

CPU="$(nproc)"
FABM_HOST=gotm

# Select ifx, ifort, or gfortran
CMAKE_Fortran_COMPILER=ifx

# You may want to use -O3 for better performance but check whether
# the change in results are acceptable.
# Default of -O2 will be used with commenting out the next line.
# FFLAGS=-O3

GOTM_BASE=~/Github/gotm/code
FABM_ERSEM_BASE=~/Github/ersem
FABM_BASE=~/Github/fabm/
CMAKE_INSTALL_PREFIX=~/local/fabm-${CMAKE_Fortran_COMPILER}/${FABM_HOST}

echo "Building GOTM-FABM-ERSEM"
rm -rf ~/build  # Delete if exits
mkdir ~/build && cd ~/build
cmake $GOTM_BASE \
  -DFABM_BASE=$FABM_BASE \
  -DFABM_ERSEM_BASE=$FABM_ERSEM_BASE \
  -DCMAKE_Fortran_COMPILER=$CMAKE_Fortran_COMPILER \
  -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
  -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
  $CMAKE_FLAG
make install -j $CPU
```

Make executable and run:
```bash
chmod +x install_ersem_gotm.sh
./install_ersem_gotm.sh
# Or with additional CMake flags:
./install_ersem_gotm.sh -f "-DIRON=ON"
```

### Compiler Options

**Supported Compilers:**
- `ifx` - Intel Fortran Compiler (oneAPI, recommended)
- `ifort` - Intel Fortran Compiler (classic)
- `gfortran` - GNU Fortran Compiler

**Optimization Levels:**
- Default: `-O2` (balanced optimization)
- High performance: `-O3` (uncomment `FFLAGS=-O3` in scripts, verify results)

**Additional CMake Options:**
- Enable iron cycling: `-f "-DIRON=ON"`
- Debug build: Uncomment `CMAKE_BUILD_TYPE=Debug` in scripts
- Custom NetCDF: `-f "-DNETCDF_CONFIG=/path/to/nc-config"`

### Python Interface (PyFABM-ERSEM)

```bash
# From fabm directory
cd ~/Github/fabm
python -m pip install .
```

## Testing

Test configurations are available in the `testcases/` directory with various YAML setups.

## Configuration

Model setup through YAML files (examples in `testcases/`):
- Define state variables, parameters, and coupling
- Different configurations for various complexity levels
- Template: `fabm-ersem.yaml.template`

## Documentation

**General ERSEM Documentation:**
- Official ERSEM website: http://ersem.com
- Upstream repository: https://github.com/pmlmodelling/ersem
- Main publication: Butenschön et al. (2016), GMD, doi:10.5194/gmd-9-1293-2016

## How to Cite

To refer to ERSEM in publications, please cite:

Butenschön, M., Clark, J., Aldridge, J.N., Allen, J.I., Artioli, Y., Blackford, J., Bruggeman, J., Cazenave, P., Ciavatta, S., Kay, S., Lessin, G., van Leeuwen, S., van der Molen, J., de Mora, L., Polimene, L., Sailley, S., Stephens, N., Torres, R. (2016). ERSEM 15.06: a generic model for marine biogeochemistry and the ecosystem dynamics of the lower trophic levels. Geoscientific Model Development, 9(4), 1293–1339. doi: [10.5194/gmd-9-1293-2016](https://doi.org/10.5194/gmd-9-1293-2016).

To refer specifically to the ERSEM source code, you may use its Zenodo DOI (badge at top of page).

## License

See upstream repository for license information: https://ersem.readthedocs.io/en/latest/license.html
