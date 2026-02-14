# Phytoplankton Temperature Response Enhancement: CTMI Model

## Overview

This document proposes enhancing the ERSEM primary producer temperature response by adding the **Cardinal Temperature Model with Inflection (CTMI)** as an alternative to the current Q10-based formulation. The CTMI allows species-specific optimal temperatures (Topt), enabling representation of seasonal phytoplankton succession driven by thermal niche partitioning.

The primary motivation is improving dissolved oxygen (DO) reproduction in highly eutrophic systems such as Tokyo Bay, where the current Q10 formulation cannot represent winter diatom blooms or the diatom-to-dinoflagellate succession.

## Problem Statement

### Current Temperature Response in ERSEM

The temperature response in `primary_producer.F90` (line 356) is:

```fortran
et = max(0.0_rk, self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk))
```

Mathematically:

```
f(T) = max(0, Q10^((T - 10)/10) - Q10^((T - 32)/3))
```

- First term: exponential increase referenced to 10 degC
- Second term: steep suppression above ~32 degC

With the standard Q10 = 2.0 (shared by P1-P4), the effective temperature response is:

| Temperature (degC) | f(T)  | Interpretation |
|---------------------|-------|----------------|
| 5                   | 0.71  | Winter minimum |
| 10                  | 1.00  | Reference      |
| 15                  | 1.41  |                |
| 20                  | 1.99  |                |
| 25                  | 2.78  |                |
| 28                  | 3.36  | Near peak      |
| 30                  | 3.56  | Effective peak |
| 32                  | 2.83  | Rapid decline  |
| 35                  | 0     | Zero growth    |

**Key limitations:**

1. **Monotonically increasing up to ~30 degC**: No mechanism for species that grow best at lower temperatures. A diatom adapted to 10-15 degC cannot outcompete a dinoflagellate at those temperatures if both share the same curve.

2. **All functional groups share identical temperature responses**: P1-P4 all use Q10 = 2.0 with hardcoded reference (10 degC) and suppression (32 degC) temperatures. The groups can only be differentiated by their `sum` (maximum productivity) parameter, not by their thermal niche.

3. **Cannot represent seasonal succession**: The observed pattern in eutrophic bays -- winter diatom blooms, spring-summer transition to dinoflagellates, autumn diatom resurgence -- is fundamentally driven by differential temperature optima across functional groups. A single monotonic curve cannot capture this crossover.

### Impact on Model Performance

In the TB-GOTM optimization study (ERSEM_GOTM_opt_v3, 1200 trials), surface DO correlation was only 0.43 despite optimizing P1-P4 growth rates. The optimizer pushed P1 (diatoms) to the lower bound and P4 (microphytoplankton) to the upper bound, indicating the model cannot simultaneously reproduce winter-high and summer-high primary production with a shared temperature function.

## Proposed Solution: CTMI

### The CTMI Formula

The implementation uses a polynomial form of the Cardinal Temperature Model with Inflection,
satisfying the same four defining conditions as the original Rosso et al. (1993) formulation:
f(Tmin)=0, f(Tmax)=0, f(Topt)=1, f'(Topt)=0.

The original Rosso rational form (cubic/linear) has a singularity when Topt < (Tmin+Tmax)/2,
which occurs for cold-adapted species like diatoms. The polynomial form eliminates this issue.

```
f(T) = 0                                              if T <= Tmin or T >= Tmax

f(T) = (T - Tmin)(T - Tmax)(ctmi_a * T + ctmi_b)      if Tmin < T < Tmax

where:
  a = Topt - Tmin
  b = Topt - Tmax
  ctmi_a = -(a + b) / (a * b)^2
  ctmi_b = (a * b + (a + b) * Topt) / (a * b)^2
```

The coefficients ctmi_a and ctmi_b are precomputed from Tmin, Topt, Tmax during initialization.
A safety clamp max(0, min(1, f)) is applied at runtime.

**Properties:**
- f(Topt) = 1.0 exactly (by construction)
- f(Tmin) = f(Tmax) = 0.0
- Unimodal cubic polynomial, singularity-free for all valid Tmin < Topt < Tmax
- All three parameters have direct biological meaning:
  - **Tmin**: minimum temperature for growth (degC)
  - **Topt**: temperature at which growth rate is maximum (degC)
  - **Tmax**: maximum temperature for growth (degC)

### Comparison of Temperature Response Models

| Model | Formula | Captures Topt? | Per-Group? | Parameters |
|-------|---------|---------------|------------|------------|
| ERSEM Q10 (current) | Exponential + high-T cutoff | No (monotonic) | No (same Q10) | Q10 |
| Eppley (1972) | Exponential envelope | No (monotonic) | No | a, b |
| Eppley-Norberg (2004) | Envelope x quadratic niche | Yes | Yes | a, b, z, w |
| **CTMI (proposed)** | Cardinal temperature | **Yes** | **Yes** | **Tmin, Topt, Tmax** |
| Sharpe-Schoolfield | Boltzmann-Arrhenius | Yes | Yes | Ea, Ed, Topt |

The CTMI is preferred over the Eppley-Norberg model because:
1. Parameters are directly measurable from growth experiments
2. No coupling to the Eppley envelope (which may overestimate growth; Bissinger et al. 2008)
3. Simplest parameterization (3 vs 4 parameters)
4. Validated across 15 microalgal species (Bernard & Remond 2012)

### Note on `sum` Semantics

With the current Q10 formulation, `sum` is defined as "maximum specific productivity at reference temperature" (10 degC). With CTMI, `sum` becomes "maximum specific productivity at optimal temperature" (Topt). Since f(Topt) = 1.0 for both formulations (Q10 at T=10 and CTMI at T=Topt), the `sum` parameter retains consistent units and meaning. However, numerical values may need adjustment when switching between formulations.

## Implementation Design

### New Parameters

Added to the `type_ersem_primary_producer` type definition:

| Parameter | Type | YAML Key | Units | Default | Description |
|-----------|------|----------|-------|---------|-------------|
| `iswTemp` | integer | `iswTemp` | - | 1 | Temperature response switch (1: Q10, 2: CTMI) |
| `Tmin` | real | `Tmin` | degC | 0.0 | Minimum temperature for growth (CTMI only) |
| `Topt` | real | `Topt` | degC | 20.0 | Optimal temperature for growth (CTMI only) |
| `Tmax` | real | `Tmax` | degC | 35.0 | Maximum temperature for growth (CTMI only) |

When `iswTemp = 1` (default), the original Q10 formulation is used and Tmin/Topt/Tmax are ignored. Full backward compatibility is preserved.

### Code Changes

Only **one file** requires modification: `src/primary_producer.F90`

**Type definition** (after line 69):
```fortran
integer :: iswTemp
real(rk) :: Tmin, Topt, Tmax
```

**Parameter registration** (after line 106, in `initialize` subroutine):
```fortran
call self%get_parameter(self%iswTemp, 'iswTemp', '', &
     'temperature response (1: Q10 with high-T suppression, 2: CTMI cardinal temperature model)', &
     default=1, minimum=1, maximum=2)
call self%get_parameter(self%Tmin, 'Tmin', 'degrees_Celsius', &
     'minimum temperature for growth (CTMI only)', default=0.0_rk)
call self%get_parameter(self%Topt, 'Topt', 'degrees_Celsius', &
     'optimal temperature for growth (CTMI only)', default=20.0_rk)
call self%get_parameter(self%Tmax, 'Tmax', 'degrees_Celsius', &
     'maximum temperature for growth (CTMI only)', default=35.0_rk)
```

**Temperature response** (in `do` subroutine):
```fortran
! Temperature response
if (self%iswTemp == 1) then
   ! Original Q10 formulation with high-temperature suppression
   et = max(0.0_rk, self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk))
else
   ! Polynomial CTMI (Cardinal Temperature Model with Inflection)
   ! Singularity-free cubic: f(Tmin)=0, f(Tmax)=0, f(Topt)=1, f'(Topt)=0
   ! Coefficients ctmi_a, ctmi_b precomputed in initialize.
   if (ETW <= self%Tmin .or. ETW >= self%Tmax) then
      et = 0.0_rk
   else
      et = (ETW - self%Tmin) * (ETW - self%Tmax) &
         * (self%ctmi_a * ETW + self%ctmi_b)
      et = max(0.0_rk, min(1.0_rk, et))
   end if
end if
```

### No Changes Required

- `ersem_model_library.F90`: No new model class (modifying existing `primary_producer`)
- `CMakeLists.txt`: No new source files
- Build script: Same compilation process

## Tokyo Bay Application

### Phytoplankton Seasonal Succession in Tokyo Bay

Based on long-term monitoring (Nakada et al. 2021; Nagai et al. 2022), Tokyo Bay phytoplankton can be categorized by thermal niche:

| Season | Temperature | Dominant Taxa | Functional Group |
|--------|------------|---------------|------------------|
| Winter (Dec-Feb) | 8-15 degC | *S. japonicum*, *S. marinoi-dohrnii* | Cold-adapted diatoms |
| Spring (Mar-May) | 15-20 degC | *Chaetoceros*, *Thalassiosira* | Spring diatoms |
| Summer (Jun-Aug) | 25-30 degC | *Prorocentrum*, *Heterosigma* | Dinoflagellates/raphidophytes |
| Autumn (Sep-Nov) | 15-25 degC | Mixed diatoms | Transitional community |

Key observations:
- *Skeletonema japonicum* blooms at the coldest temperatures (8-10 degC) and disappears above 15 degC (Nagai et al. 2022)
- Summer red tides driven by *Heterosigma akashiwo* and dinoflagellates require warm water (>22 degC)
- The transition period (May-June) sees a sharp shift from diatom to dinoflagellate dominance

### Proposed CTMI Parameters for Tokyo Bay

Based on literature values and thermal ecology:

| Group | ERSEM Instance | Represents | Tmin (degC) | Topt (degC) | Tmax (degC) |
|-------|---------------|------------|------|------|------|
| P1 | Diatoms | Cold-water diatoms (*S. japonicum*, *Chaetoceros*) | 2 | 15 | 30 |
| P2 | Nanophytoplankton | Small flagellates, Cryptophytes | 5 | 20 | 33 |
| P3 | Picophytoplankton | Small autotrophs, cyanobacteria-like | 8 | 25 | 35 |
| P4 | Microphytoplankton | Dinoflagellates (*Prorocentrum*, *Heterocapsa*) | 10 | 25 | 35 |

**Rationale for P1 (Diatoms, Topt = 15 degC)**:
- *S. japonicum* is absent above 25 degC, peaks at 8-15 degC (Nagai et al. 2022)
- Chesapeake Bay diatom model uses Topt ~12-18 degC (Cerco & Cole 1993)
- This enables February blooming and suppresses diatoms in summer, matching observations

**Rationale for P4 (Dinoflagellates, Topt = 25 degC)**:
- *P. minimum* bloom model uses Topt ~20-25 degC (Li et al. 2021)
- Anderson & Barton (2021) show dinoflagellates have the weakest Q10 but highest Topt among functional types
- Summer red tides in Tokyo Bay require warm-water specialization

### Expected Temperature Response Curves

With the proposed CTMI parameters:

| Temperature | P1 (diatoms) | P2 (nano) | P3 (pico) | P4 (dino) |
|------------|-------------|-----------|-----------|-----------|
| 5 degC | 0.42 | 0 | 0 | 0 |
| 10 degC | 0.86 | 0.53 | 0.11 | 0 |
| 15 degC | **1.00** | 0.88 | 0.48 | 0.44 |
| 20 degC | 0.88 | **1.00** | 0.84 | 0.83 |
| 25 degC | 0.53 | 0.86 | **1.00** | **1.00** |
| 28 degC | 0.23 | 0.64 | 0.93 | 0.92 |
| 30 degC | 0 | 0.42 | 0.78 | 0.78 |

This produces the desired succession: diatoms dominate in winter-spring (5-15 degC), nanophytoplankton peak in spring (20 degC), and pico/dinoflagellates dominate in summer (25-30 degC).

## Optimization Considerations

### CTMI Parameters as Tuning Targets

The Topt values for each group are prime candidates for Optuna optimization:

```python
tuning_targets = [
    # Maximum productivity
    ("P1", "sum", "parameters", 0.5, 3.0),
    ("P2", "sum", "parameters", 0.5, 3.0),
    ("P3", "sum", "parameters", 0.5, 3.0),
    ("P4", "sum", "parameters", 0.5, 3.0),
    # Optimal temperatures
    ("P1", "Topt", "parameters", 10.0, 20.0),  # diatoms: cold
    ("P4", "Topt", "parameters", 20.0, 30.0),  # dinoflagellates: warm
]
```

Note: Tmin and Tmax are less sensitive and can be fixed at literature values, while Topt directly controls the timing of seasonal peaks.

### Interaction with Other Parameters

The CTMI interacts strongly with:
- **a0w** (light attenuation): Controls euphotic zone depth, affecting which species can bloom
- **B1 sum** (bacterial decomposition): Determines O2 consumption relative to production
- **R6 rm** (POM sinking): Determines whether organic matter is remineralized in the water column or sediment

A comprehensive optimization (v4) should include CTMI Topt values alongside these parameters.

## Build and Test Procedure

1. Modify `~/Github/ersem/src/primary_producer.F90` as described above
2. Clean build: `rm -rf ~/build`
3. Compile: `~/Github/fabm/src/drivers/gotm/install_ersem_gotm.sh`
4. Copy binary to TB-GOTM directories (not symlink)
5. Backward compatibility test: Run with `iswTemp: 1` (default), compare output to pre-modification run
6. CTMI test: Update `fabm.yaml` with proposed CTMI parameters, run simulation

## References

### CTMI Model
- Rosso, L., Lobry, J.R., & Flandrois, J.P. (1993). An unexpected correlation between cardinal temperatures of microbial growth highlighted by a new model. *J. Theor. Biol.*, 162(4), 447-463.
- Bernard, O. & Remond, B. (2012). Validation of a simple model accounting for light and temperature effect on microalgal growth. *Bioresource Technology*, 123, 520-527.

### Temperature-Growth Relationships
- Eppley, R.W. (1972). Temperature and phytoplankton growth in the sea. *Fishery Bulletin*, 70(4), 1063-1085.
- Norberg, J. (2004). Biodiversity and ecosystem functioning: A complex adaptive systems approach. *Limnol. Oceanogr.*, 49(4 part 2), 1269-1277.
- Bissinger, J.E. et al. (2008). Predicting marine phytoplankton maximum growth rates from temperature: Improving on the Eppley curve. *Limnol. Oceanogr.*, 53(2), 487-493.
- Anderson, S.I. & Barton, A.D. (2021). Marine phytoplankton functional types exhibit diverse responses to thermal change. *Nature Communications*, 12, 6413.
- Thomas, M.K. et al. (2012). A global pattern of thermal adaptation in marine phytoplankton. *Science*, 338, 1085-1088.

### Tokyo Bay Phytoplankton Ecology
- Nakada, S. et al. (2021). Phytoplankton species abundance in Tokyo Bay (Japan) from 1998 to 2019. *Ecological Research*, 36.
- Nagai, S. et al. (2022). Temporal niche partitioning of *Skeletonema*: seasonal succession in Tokyo Bay. *Aquatic Microbial Ecology*, 89, ame02000.
- Ogawa, Y. & Ichimura, S. (1997). Changes in Red Tide Events and Phytoplankton Community Composition in Tokyo Bay 1907-1997. *Umi/La Mer*, 7(3), 159.

### Eutrophic Bay Models
- Li, M. et al. (2021). A three-dimensional mechanistic model of *Prorocentrum minimum* blooms in eutrophic Chesapeake Bay. *Sci. Total Environ.*, 769, 144528.
- Cerco, C.F. & Cole, T. (1993). Three-dimensional eutrophication model of Chesapeake Bay. *J. Environ. Eng.*, 119(6), 1006-1025.
- Sohma, A. et al. (2008). A benthic-pelagic coupled ecosystem model for Tokyo Bay. *Ecol. Model.*, 215(1-3).

### ERSEM Model
- Butenschon, M. et al. (2016). ERSEM 15.06: a generic model for marine biogeochemistry. *Geosci. Model Dev.*, 9, 1293-1339.
- Blackford, J.C. et al. (2004). Ecosystem dynamics at six contrasting sites. *J. Mar. Syst.*, 52, 191-215.

---

*Document created: 2026-02-14*
*Target file: `~/Github/ersem/src/primary_producer.F90` line 356*
