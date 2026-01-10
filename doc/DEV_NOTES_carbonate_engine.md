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
      opt_k_carbonic: 1         # K1/K2 formulation (1: Lueker 2000, 2: Millero 2010)
      opt_total_borate: 2       # Total boron (1: Uppstrom 1974, 2: Lee 2010)
      ta_slope: 43.626          # TA(S) regression slope (for iswtalk=6)
      ta_intercept: 846.48      # TA(S) regression intercept (for iswtalk=6)
```

### Engine Parameters

| Parameter | Values | Description |
|-----------|--------|-------------|
| `opt_k_carbonic` | 1 (default) | Lueker et al. (2000) - valid for S=19-43 |
|                  | 2 | Millero (2010) - valid for S=1-50 |
| `opt_total_borate` | 1 | Uppstrom (1974) - BT = 416 × S/35 µmol/kg |
|                    | 2 (default) | Lee et al. (2010) - BT = 432.6 × S/35 µmol/kg |

**Note for estuarine applications (e.g., Tokyo Bay):** For low salinity waters (S < 19),
use `opt_k_carbonic: 2` (Millero 2010) which is valid for S=1-50.

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

| Constant | Source | Notes |
|----------|--------|-------|
| K0 | Weiss (1974) | Pressure correction only at depth (Pbar > 0) |
| K1, K2 | Lueker et al. (2000) | Default (`opt_k_carbonic=1`), valid S=19-43 |
|        | Millero (2010) | `opt_k_carbonic=2`, valid S=1-50 |
| KB | Dickson (1990) via Millero (1995) | Total scale |
| KW | Millero (1995) | Converted from SWS to Total scale |
| KS | Dickson (1990) | Free scale |
| KF | Perez & Fraga (1987) | Free scale |
| Total sulfate (ST) | Morris & Riley (1966) | ST = 0.14 × Cl / 96.062 |
| Total fluoride (FT) | Riley (1965) | FT = 0.000067 × Cl / 18.998 |
| Total boron (BT) | Uppstrom (1974) | `opt_total_borate=1`: 416 × S/35 µmol/kg |
|                  | Lee et al. (2010) | `opt_total_borate=2`: 432.6 × S/35 µmol/kg |

All constants include pressure corrections from Millero (1995).

### Low Salinity Support

The solver handles S=0 (freshwater) gracefully:
- At S=0, total sulfate, fluoride, and boron are zero
- pH scale conversion factors become unity
- K1/K2 formulas are extrapolated (use Millero 2010 for best results at low S)

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

1. **Salinity range for K1/K2:**
   - Lueker 2000 (`opt_k_carbonic=1`): Valid for S=19-43
   - Millero 2010 (`opt_k_carbonic=2`): Valid for S=1-50
   - For estuarine applications with S < 19, use Millero 2010

2. **Deep water validation:** Pressure correction code is present but should be
   validated for deep water applications (Pbar > 0).

3. **pH bounds:** The solver uses fixed pH bounds [2, 12]. Extremely unusual
   DIC/TA combinations outside typical marine ranges may fail to converge.

## Accuracy

The carbonate-engine achieves excellent agreement with PyCO2SYS:
- pH accuracy: ~6.5 × 10⁻⁷ (sub-µpH level)
- fCO2 accuracy: ~0.001 µatm

Test suite includes 20 cases covering:
- Temperature: 5-30°C
- Salinity: 5-38 PSU (including estuarine conditions)
- All tests pass with tolerances of 10⁻⁶ pH and 0.01 µatm

## References

- Dickson, A. G. (1990). Standard potential of the reaction: AgCl(s) + 1/2H2(g) = Ag(s) + HCl(aq). J. Chem. Thermodyn., 22, 113-127.
- Endo, T., et al. (2023). Carbonate chemistry in Tokyo Bay. Frontiers in Marine Science.
- Lee, K., et al. (2010). The universal ratio of boron to chlorinity for the North Pacific and North Atlantic oceans. Geochim. Cosmochim. Acta, 74, 1801-1811.
- Lueker, T. J., Dickson, A. G., & Keeling, C. D. (2000). Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. Mar. Chem., 70, 105-119.
- Millero, F. J. (1995). Thermodynamics of the carbon dioxide system in the oceans. Geochim. Cosmochim. Acta, 59, 661-677.
- Millero, F. J. (2010). Carbonate constants for estuarine waters. Mar. Freshwater Res., 61, 139-142.
- Morris, A. W., & Riley, J. P. (1966). The bromide/chlorinity and sulphate/chlorinity ratio in sea water. Deep-Sea Res., 13, 699-705.
- Perez, F. F., & Fraga, F. (1987). Association constant of fluoride and hydrogen ions in seawater. Mar. Chem., 21, 161-168.
- Riley, J. P. (1965). The occurrence of anomalously high fluoride concentrations in the North Atlantic. Deep-Sea Res., 12, 219-220.
- Uppstrom, L. R. (1974). The boron/chlorinity ratio of deep-sea water from the Pacific Ocean. Deep-Sea Res., 21, 161-162.
- Weiss, R. F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. Mar. Chem., 2, 203-215.
