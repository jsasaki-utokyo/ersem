# Developer Notes: Carbonate Engine

## Overview

The carbonate-engine is an alternative carbonate chemistry solver for ERSEM that
reproduces PyCO2SYS results for the DIC+TA -> pH+pCO2 calculation without requiring
Python at runtime.

## Engine Switch

The `engine` parameter in the carbonate model controls which solver is used:

| engine | Description |
|--------|-------------|
| 0 | Legacy ERSEM carbonate solver (CO2DYN) - **default** |
| 1 | New carbonate-engine solver (PyCO2SYS-style) |

The default (`engine=0`) ensures full backward compatibility with existing simulations.

## What engine=1 Does

When `engine=1` is selected, the carbonate model uses `carbonate_engine_solve()`
instead of `CO2DYN()` for computing:

- pH (on the requested scale: Total, SWS, or backward-compatible SWS)
- pCO2 (partial pressure of CO2 in atm)
- Carbonate species: H2CO3, HCO3, CO3 (all in mol/kg)
- K0 (Henry's law constant)

The solver:
1. Calculates equilibrium constants (K0, K1, K2, KB, KW, KS, KF) consistent with PyCO2SYS
2. Uses Brent's method for robust pH root-finding
3. Returns species concentrations and pCO2 from the solved pH

## What Remains Legacy

Even with `engine=1`, the following still use existing ERSEM code:

- Air-sea CO2 flux calculations (Schmidt number, gas transfer velocity)
- CaCO3 saturation state calculations (Om_cal, Om_arg)
- Alkalinity approximation formulas (when iswtalk != 5)
- Unit conversions (mmol/m³ <-> mol/kg)

## Configuration

### YAML Example

```yaml
ersem:
  O3:                           # carbonate model
    model: ersem/carbonate
    parameters:
      iswCO2: 1                 # carbonate diagnostics on
      iswASFLUX: 6              # Wanninkhof 2014 gas exchange
      iswtalk: 5                # dynamic alkalinity (state variable)
      pHscale: 1                # Total pH scale
      engine: 1                 # carbonate-engine solver
      ta_slope: 43.626          # TA(S) regression slope (for iswtalk=6)
      ta_intercept: 846.48      # TA(S) regression intercept (for iswtalk=6)
```

### TA(S) Regression Parameters

When `iswtalk=6`, the model uses a linear TA(S) regression:

```
TA = ta_intercept + ta_slope * S
```

Default values (Endo et al., 2023, Tokyo Bay):
- `ta_slope`: 43.626 µmol/kg/PSU
- `ta_intercept`: 846.48 µmol/kg

These can be adjusted for different regional settings.

## Running the Fortran Test

```bash
cd tests/carbonate_engine

# Compile (standalone, no FABM needed)
gfortran -O2 -o test_driver test_driver.f90

# Run tests
./test_driver reference_cases.csv
```

The test driver validates the solver against PyCO2SYS reference values.

## Regenerating Reference Cases

```bash
cd tests/carbonate_engine

# Requires PyCO2SYS
pip install PyCO2SYS

# Generate fresh reference values
python generate_reference.py
```

## Implementation Details

### Equilibrium Constants

The carbonate-engine uses the following constant formulations:

| Constant | Source |
|----------|--------|
| K0 | Weiss (1974) |
| K1, K2 | Millero (2010) - Total scale (opt 10 in PyCO2SYS) |
| KB | Dickson (1990) via Millero (1995) |
| KW | Millero (1995) |
| KS | Dickson (1990) + Perez & Fraga (1987) |
| KF | Perez & Fraga (1987) |
| Total boron | Lee et al. (2010) - BT/S = 0.0004157 |

All constants include pressure corrections from Millero (1995).

### Solver Algorithm

1. Calculate equilibrium constants at given T, S, P
2. Calculate total boron from salinity
3. Define F(H) = TA_model(DIC, H, K1, K2, KB, KW, BT) - TA_input
4. Use Brent's method to find H where F(H) = 0
5. Calculate carbonate species from solved H
6. Return pH = -log10(H) on requested scale

### pH Scale Handling

The solver always works internally on the Total pH scale. For other scales:

- `phscale=1` (Total): pH returned directly
- `phscale=0` (SWS): pH converted using total2sws factor
- `phscale=-1` (backward compat): Uses old K1/K2 formulas, SWS output

### Fallback Behavior

If `carbonate_engine_solve()` returns `success=.false.`, the carbonate model
automatically falls back to using the previous timestep's diagnostic values,
identical to the legacy solver's behavior.

## File Locations

- **Solver implementation**: `src/carbonate_engine.F90`
- **Integration point**: `src/carbonate.F90` (do() and do_surface() subroutines)
- **Test harness**: `tests/carbonate_engine/`
- **Build integration**: `src/CMakeLists.txt`

## Known Limitations

1. The carbonate-engine currently only supports surface conditions (Pbar=0).
   Pressure correction code is present but should be validated for deep water.

2. Minor numerical differences (~0.001 pH) from PyCO2SYS are expected due to
   floating-point precision and minor constant formulation differences.

3. The solver uses fixed pH bounds [2, 12]. Extremely unusual DIC/TA combinations
   outside typical marine ranges may fail to converge.

## References

- Weiss, R. F. (1974). CO2 in water and seawater: the solubility of a non-ideal gas.
- Millero, F. J. (2010). Carbonate constants for estuarine waters.
- Dickson, A. G. (1990). Standard potential of the reaction: AgCl(s) + 1/2H2(g) = Ag(s) + HCl(aq).
- Perez, F. F., & Fraga, F. (1987). Association constant of fluoride and hydrogen ions in seawater.
- Lee, K., et al. (2010). The universal ratio of boron to chlorinity for the North Pacific and North Atlantic oceans.
- Endo, T., et al. (2023). Carbonate chemistry in Tokyo Bay. Frontiers in Marine Science.
