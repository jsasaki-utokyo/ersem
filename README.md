[![DOI](https://zenodo.org/badge/302390544.svg)](https://zenodo.org/badge/latestdoi/302390544)

# ERSEM - UTokyo Estuarine Configuration

![ERSEM diagram](docs/images/ERSEM.png)

## Overview

This repository contains a development version of ERSEM (European Regional Seas Ecosystem Model) configured for estuarine and coastal applications at UTokyo. ERSEM is a marine biogeochemical and ecosystem model describing the cycling of carbon, nitrogen, phosphorus, silicon, oxygen and iron through the lower trophic level pelagic and benthic ecosystems.

**Key Features:**
- Multi-element biogeochemical cycling (C, N, P, Si, O, Fe)
- Pelagic and benthic ecosystem components
- Integration with FABM (Framework for Aquatic Biogeochemical Models)
- Multiple deployment options:
  - **0D**: Box models (FABM-ERSEM-0D)
  - **1D**: Water column models (FABM-ERSEM-GOTM)
  - **3D**: Coupled with FVCOM, ROMS, NEMO, and other ocean models (FABM-ERSEM-FVCOM, etc.)

**Note:** This is a development version optimized for estuarine and coastal systems. The upstream version is maintained at https://github.com/pmlmodelling/ersem (master branch).

## Installation from Source

### Prerequisites

**Required for all configurations:**
- Fortran compiler (Intel ifx or gfortran)
- CMake (>= 3.10)
- NetCDF libraries
- Git

**Additional requirements depend on your use case:**
- **FABM-ERSEM-0D / FABM-ERSEM-GOTM**: GOTM (General Ocean Turbulence Model)
- **FABM-ERSEM-FVCOM**: FVCOM (Finite Volume Community Ocean Model) with FABM coupler
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

Only needed for FABM-ERSEM-0D and FABM-ERSEM-GOTM configurations (GOTM v6 required):

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

**Note:** This uses a customized FVCOM v5.1.0-dev with FABM integration developed for estuarine applications at UTokyo. See `FVCOM/branches/uk-fabm-v5.1.0-dev/README-FABM.md` for detailed documentation.

### 2. Build Options

ERSEM can be coupled with different hydrodynamic models through FABM. Choose the configuration that matches your application:

| Configuration | Description | Use Case |
|--------------|-------------|----------|
| **FABM-ERSEM-0D** | 0-dimensional box model | Aquarium setups, mesocosm experiments, testing |
| **FABM-ERSEM-GOTM** | 1-dimensional water column | Vertical profiling, L4 station, single-point simulations |
| **FABM-ERSEM-FVCOM** | 3D unstructured grid | Coastal regions, estuaries, complex bathymetry |

#### Build Instructions by Configuration

---

#### Option 1: FABM-ERSEM-0D (Box Model / Aquarium Setup)

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

# Use minimum of 8 cores or nproc/4, with floor of 4
NPROC=$(nproc)
CPU=$(( NPROC < 32 ? NPROC : NPROC / 4 ))
CPU=$(( CPU > 8 ? 8 : CPU ))
CPU=$(( CPU < 4 ? 4 : CPU ))

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

#### Option 2: FABM-ERSEM-GOTM (1D Water Column Model)

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

# Use minimum of 8 cores or nproc/4, with floor of 4
NPROC=$(nproc)
CPU=$(( NPROC < 32 ? NPROC : NPROC / 4 ))
CPU=$(( CPU > 8 ? 8 : CPU ))
CPU=$(( CPU < 4 ? 4 : CPU ))

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

echo "Building FABM-ERSEM-GOTM"
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

#### Option 3: FABM-ERSEM-FVCOM (3D Unstructured Grid Model)

For 3-dimensional coastal and estuarine simulations with unstructured grids.

**Prerequisites:**
- FVCOM with FABM integration (see [FVCOM repository cloning](#fvcom-for-3d-unstructured-grid-applications) above)
- All core components (FABM, ERSEM)
- NetCDF libraries

**Step 1: Build FABM with FVCOM Driver**

```bash
#!/bin/bash
# Build FABM with FVCOM host driver

# Use minimum of 8 cores or nproc/4, with floor of 4
NPROC=$(nproc)
CPU=$(( NPROC < 32 ? NPROC : NPROC / 4 ))
CPU=$(( CPU > 8 ? 8 : CPU ))
CPU=$(( CPU < 4 ? 4 : CPU ))

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

make -j ${CPU}
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
# Use minimum of 8 cores or nproc/4, with floor of 4
NPROC=$(nproc)
CPU=$(( NPROC < 32 ? NPROC : NPROC / 4 ))
CPU=$(( CPU > 8 ? 8 : CPU ))
CPU=$(( CPU < 4 ? 4 : CPU ))

cd ~/Github/FVCOM/src
make clean
make -j ${CPU}
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

For carbonate engine testing, see `tests/carbonate_engine/` directory which contains:
- Reference test cases comparing Engine 0 and Engine 1 outputs
- Validation against PyCO2SYS calculations
- Test driver scripts

## Carbonate Chemistry Engine

### Overview

ERSEM includes a configurable carbonate chemistry solver with two engine options:

| Engine | Description | pH/pCO2 Calculation | Default Constants |
|--------|-------------|---------------------|-------------------|
| **0** (Legacy) | Original ERSEM solver | Follows et al. (2006) iterative | Millero (1995) K1/K2 |
| **1** (New) | PyCO2SYS-compatible solver | Newton-Raphson iteration | Lueker et al. (2000) K1/K2 |

**Engine 1** was developed to reproduce PyCO2SYS calculations and provides more accurate carbonate chemistry, particularly for estuarine applications with variable salinity.

### Building with the New Engine

The carbonate engine is included in the standard ERSEM build. No special build flags are required.

```bash
# Standard build includes both engines
./install_ersem_gotm.sh

# The engine is selected at runtime via fabm.yaml configuration
```

### Configuration Parameters

Configure the carbonate system in your `fabm.yaml` file under the `O3` (carbonate) instance:

```yaml
instances:
  O3:
    model: ersem/carbonate
    parameters:
      # Engine selection
      engine: 1                    # 0: Legacy ERSEM, 1: PyCO2SYS-style

      # Constants for Engine 1 (ignored when engine=0)
      opt_k_carbonic: 1            # 1: Lueker 2000, 2: Millero 2010
      opt_total_borate: 1          # 1: Uppstrom 1974, 2: Lee 2010

      # Alkalinity formulation
      iswtalk: 5                   # 1-4: from S/T, 5: dynamic, 6: custom TA(S)

      # Custom TA(S) regression (for iswtalk=6)
      ta_slope: 43.626             # TA = ta_intercept + ta_slope × S
      ta_intercept: 846.48         # (µmol/kg)

      # pH scale
      pHscale: 1                   # 1: total, 0: SWS, -1: SWS backward compatible

      # DIC relaxation (optional, for 1D models)
      relax_c: 0.0                 # Relaxation timescale (days), 0=off
      c_ta_ratio: 0.0              # Target DIC/TA ratio, 0=use fixed target
      c_relax_target: 2100.0       # Fixed DIC target (mmol/m³)
```

### Engine Parameter Reference

#### `engine` - Carbonate Chemistry Solver

| Value | Description | Use Case |
|-------|-------------|----------|
| 0 | Legacy ERSEM solver (Follows et al. 2006) | Backward compatibility, open ocean |
| 1 | PyCO2SYS-style solver | Estuaries, validation studies |

#### `opt_k_carbonic` - K1/K2 Constants (Engine 1 only)

| Value | Reference | Valid Range |
|-------|-----------|-------------|
| 1 | Lueker et al. (2000) | S: 19-43, T: 2-35°C |
| 2 | Millero (2010) | S: 1-50, T: 0-50°C |

#### `opt_total_borate` - Total Boron (Engine 1 only)

| Value | Reference | Notes |
|-------|-----------|-------|
| 1 | Uppstrom (1974) | PyCO2SYS default |
| 2 | Lee et al. (2010) | Higher accuracy |

#### `iswtalk` - Alkalinity Formulation

| Value | Description |
|-------|-------------|
| 1 | Millero et al. (1998) - Atlantic |
| 2 | Millero et al. (1998) - S < 35 |
| 3 | Linear TA(T) - North Sea |
| 4 | Linear TA(T) - high accuracy |
| 5 | Dynamic alkalinity (state variable) |
| 6 | Custom TA(S) regression |

#### `pHscale` - pH Scale

| Value | Description |
|-------|-------------|
| 1 | Total scale (recommended) |
| 0 | Seawater scale (SWS) |
| -1 | SWS with backward-compatible offset |

### Backward Compatibility

The default configuration maintains full backward compatibility:

```yaml
O3:
  parameters:
    engine: 0          # Legacy solver (default)
    pHscale: -1        # SWS with backward-compatible offset
    iswtalk: 5         # Dynamic alkalinity
```

**Key points for backward compatibility:**

1. **Engine 0 is default**: Existing configurations without `engine` parameter will use the legacy solver

2. **pH scale offset**: When `pHscale: -1`, the legacy pH offset (0.009) is applied for SWS scale consistency with older ERSEM versions

3. **Alkalinity options**: All existing `iswtalk` options (1-5) work identically in both engines

4. **No rebuild required**: Engine selection is a runtime parameter in `fabm.yaml`

### Switching Between Engines

To switch from legacy to new engine:

```yaml
# Before (legacy)
O3:
  parameters:
    engine: 0
    pHscale: -1

# After (PyCO2SYS-compatible)
O3:
  parameters:
    engine: 1
    opt_k_carbonic: 1      # Lueker 2000
    opt_total_borate: 1    # Uppstrom 1974
    pHscale: 1             # Total scale
```

**Expected differences:**
- pH values differ by ~0.01-0.02 due to different K1/K2 constants
- pCO2 values may differ by ~5-20 µatm
- Largest differences occur at low salinity and extreme temperatures

### Validation

Engine 1 has been validated against PyCO2SYS v1.8.0:

| Parameter | Agreement |
|-----------|-----------|
| pH | < 0.001 units |
| pCO2 | < 0.1 µatm |
| HCO3⁻ | < 0.01 mmol/m³ |
| CO3²⁻ | < 0.01 mmol/m³ |

See `tests/carbonate_engine/reference_cases.csv` for detailed validation data.

### DIC Relaxation (1D Models)

For 1D water column models, DIC may drift due to lack of lateral transport. The relaxation feature nudges DIC toward a target value:

```yaml
O3:
  parameters:
    relax_c: 2.0           # 2-day relaxation timescale
    c_ta_ratio: 0.88       # Target DIC = 0.88 × TA
    # OR use fixed target:
    # c_ta_ratio: 0.0
    # c_relax_target: 2000.0
```

**Recommendation:** Use `c_ta_ratio` instead of `c_relax_target` to prevent carbonate singularity when salinity (and thus TA) varies.

### References

- Follows, M.J., et al. (2006). On the solution of the carbonate chemistry system. Ocean Modelling, 12, 290-301.
- Lueker, T.J., et al. (2000). Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2. Marine Chemistry, 70, 105-119.
- Millero, F.J. (1995). Thermodynamics of the carbon dioxide system in the oceans. Geochimica et Cosmochimica Acta, 59, 661-677.
- Millero, F.J. (2010). Carbonate constants for estuarine waters. Marine and Freshwater Research, 61, 139-142.
- Uppstrom, L.R. (1974). The boron/chlorinity ratio of deep-sea water. Deep Sea Research, 21, 161-162.

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
