[![DOI](https://zenodo.org/badge/302390544.svg)](https://zenodo.org/badge/latestdoi/302390544)

# ERSEM - Japanese Waters Configuration

![ERSEM diagram](docs/images/ERSEM.png)

## Overview

This repository contains a development version of ERSEM (European Regional Seas Ecosystem Model) configured for Japanese waters. ERSEM is a marine biogeochemical and ecosystem model describing the cycling of carbon, nitrogen, phosphorus, silicon, oxygen and iron through the lower trophic level pelagic and benthic ecosystems.

**Key Features:**
- Multi-element biogeochemical cycling (C, N, P, Si, O, Fe)
- Pelagic and benthic ecosystem components
- Integration with FABM (Framework for Aquatic Biogeochemical Models)
- Multiple deployment options:
  - **0D**: Box models (FABM0D)
  - **1D**: Water column models (GOTM-ERSEM)
  - **3D**: Coupled with FVCOM, ROMS, NEMO, and other ocean models

**Note:** This is a development branch for Japanese waters. The upstream version is maintained at https://github.com/pmlmodelling/ersem (master branch).

## Installation from Source

### Prerequisites

**Required for all configurations:**
- Fortran compiler (Intel ifx or gfortran)
- CMake (>= 3.10)
- NetCDF libraries
- Git

**Additional requirements depend on your use case:**
- **FABM0D / GOTM-ERSEM**: GOTM (General Ocean Turbulence Model)
- **FVCOM-ERSEM**: FVCOM (Finite Volume Community Ocean Model) with FABM coupler
- **Other 3D models**: ROMS, NEMO, MOM, HYCOM, SCHISM (not covered here)

### 1. Clone Repositories

#### Core Components (Required for All)

```bash
# Create directory structure
mkdir -p ~/Github && cd ~/Github

# Clone FABM
git clone git@github.com:jsasaki-utokyo/fabm.git

# Clone this ERSEM repository
git clone git@github.com:jsasaki-utokyo/ersem.git
```

#### GOTM (For 0D and 1D Water Column Models)

Only needed for FABM0D and GOTM-ERSEM configurations (GOTM v6 required):

```bash
cd ~/Github
mkdir gotm
cd gotm
git clone git@github.com:jsasaki-utokyo/code.git
cd code
git submodule update --init --recursive
```

#### FVCOM (For 3D Unstructured Grid Applications)

For coastal and estuary applications with 3D hydrodynamics:

```bash
cd ~/Github

# Clone FVCOM with FABM support (estuarine-utokyo version)
git clone -b uk-fabm-v5.1.0-dev git@github.com:estuarine-utokyo/FVCOM.git
```

**Note:** This uses a customized FVCOM v5.1.0-dev with FABM integration developed for Japanese estuarine applications. See `FVCOM/branches/uk-fabm-v5.1.0-dev/README-FABM.md` for detailed documentation.

### 2. Build Options

ERSEM can be coupled with different hydrodynamic models through FABM. Choose the configuration that matches your application:

| Configuration | Description | Use Case |
|--------------|-------------|----------|
| **FABM0D** | 0-dimensional box model | Aquarium setups, mesocosm experiments, testing |
| **GOTM-ERSEM** | 1-dimensional water column | Vertical profiling, L4 station, single-point simulations |
| **FVCOM-ERSEM** | 3D unstructured grid | Coastal regions, estuaries, complex bathymetry |

#### Build Instructions by Configuration

---

#### Option 1: FABM0D (Box Model / Aquarium Setup)

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

# Specify ifx or gfortran
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

Save the script to any directory (e.g., your home directory or a scripts folder), make it executable, and run:

```bash
chmod +x install_ersem_fabm0d.sh
./install_ersem_fabm0d.sh
# Or with additional CMake flags:
./install_ersem_fabm0d.sh -f "-DIRON=ON"
```

**Note:** The script can be run from any directory. It uses absolute paths (`~/Github/`, `~/build/`, `~/local/`) and will create build directories as needed.

---

#### Option 2: GOTM-ERSEM (1D Water Column Model)

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

# Select ifx or gfortran
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

Save the script to any directory, make it executable, and run:

```bash
chmod +x install_ersem_gotm.sh
./install_ersem_gotm.sh
# Or with additional CMake flags:
./install_ersem_gotm.sh -f "-DIRON=ON"
```

**Note:** The script can be run from any directory. It uses absolute paths (`~/Github/`, `~/build/`, `~/local/`) and will create build directories as needed.

---

#### Option 3: FVCOM-ERSEM (3D Unstructured Grid Model)

For 3-dimensional coastal and estuarine simulations with unstructured grids.

**Prerequisites:**
- FVCOM with FABM integration (see [FVCOM repository cloning](#fvcom-for-3d-unstructured-grid-applications) above)
- All core components (FABM, ERSEM)
- NetCDF libraries

**Step 1: Build FABM with FVCOM Driver**

```bash
#!/bin/bash
# Build FABM with FVCOM host driver

FABM_SRC="${HOME}/Github/fabm"
ERSEM_SRC="${HOME}/Github/ersem"
INSTALL_DIR="${HOME}/local/fabm-ifx-ersem"
BUILD_DIR="${HOME}/build/fabm-fvcom-ersem"

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

cmake ${FABM_SRC} \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DFABM_HOST=fvcom \
    -DFABM_ERSEM_BASE=${ERSEM_SRC} \
    -DCMAKE_Fortran_COMPILER=ifx \
    -DCMAKE_BUILD_TYPE=Release \
    -DFABM_EMBED_VERSION=ON

make -j $(nproc)
make install
```

**Important:** Use `-DFABM_HOST=fvcom` to ensure proper 2D driver configuration.

**Step 2: Configure FVCOM Build**

Edit `FVCOM/src/make.inc`:

```fortran
# Enable FABM (around line 540)
FLAG_48 = -DFABM

# Set FABM installation paths
FABMDIR = ${HOME}/local/fabm-ifx-ersem
FABMLIB = -L$(FABMDIR)/lib64 -lfabm
FABMINC = -I$(FABMDIR)/include -I$(FABMDIR)/include/yaml
```

**Step 3: Build FVCOM**

```bash
cd ~/Github/FVCOM/src
make clean
make -j $(nproc)
```

**Step 4: Runtime Configuration**

Create runtime namelist with FABM section (`CASENAME_run.nml`):

```fortran
&NML_FABM
  FABM_MODEL = T                    ! Enable FABM
  FABM_YAML = 'fabm.yaml'           ! FABM configuration file
  NC_FABM = 'fabm_output.nc'        ! FABM output file
  FABM_DIAG_OUT = T                 ! Output diagnostic variables
  OBC_FABM_NUDGING = T              ! Enable open boundary nudging
  OBC_FABM_FILE = 'ersem_obc.nc'    ! Boundary conditions (if needed)
  STARTUP_FABM_TYPE = 'set_values'  ! or 'constant'
/
```

Create `fabm.yaml` configuration file in your case directory (see `testcases/` for examples).

**Key Features:**
- NetCDF boundary conditions for open boundaries
- River source input with biogeochemical variables
- Hot start from restart files
- Parallel MPI execution support

**Detailed Documentation:**
- See `~/Github/FVCOM/branches/uk-fabm-v5.1.0-dev/README-FABM.md` for comprehensive guide
- Includes OBC setup, river forcing, restart files, and troubleshooting

---

### Compiler Options

**Supported Compilers:**
- `ifx` - Intel Fortran Compiler (oneAPI, recommended)
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
