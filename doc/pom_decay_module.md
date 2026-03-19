# POM Decay Module (`ersem/pom_decay`)

## Overview

The `pom_decay` module adds first-order particulate organic matter (POM)
decomposition in the water column with direct O2 consumption. This process
operates in parallel with the existing bacteria-mediated POM uptake (B1 sRP
pathway) and represents particle-attached microbial communities that are
ecologically distinct from free-living pelagic bacteria.

**Added:** 2026-03-19 (jsasaki)

## Motivation

Standard ERSEM lacks direct POM decomposition in the water column. All POM
processing requires free-living bacteria (B1) via the sRP pathway, giving
effective decomposition rates of ~0.001 /d — an order of magnitude slower
than comparable models (BROM: 0.01–0.16 /d, ERGOM: 0.003–0.02 /d,
CE-QUAL-ICM: 0.02–0.075 /d).

This limitation prevents reproduction of mid-water hypoxia in eutrophic
systems where rapid decomposition of sinking organic matter (phytoplankton,
fecal pellets, marine snow) is the primary mechanism for O2 depletion below
the pycnocline. See `TB-FVCOM/ersem/docs/pom_decomposition_model_comparison.md`
for a detailed model comparison.

## Physical Basis

The first-order POM decay term parameterizes:

1. **Extracellular enzymatic hydrolysis** by particle-attached bacteria
   that colonize sinking particles within hours of formation.
2. **Respiration by particle-attached microbial communities** that consume
   the solubilized organic matter and directly consume O2.
3. **Cell autolysis** of dead phytoplankton through intracellular enzymes.

Because particle-attached bacteria scale proportionally with POM concentration
(rapid colonization), the net decomposition rate is well approximated by
first-order kinetics: dPOM/dt = −k × POM.

This process is NOT redundant with B1 (free-living bacteria). B1 primarily
consumes dissolved organic matter (DOM: R1, R2) and has low affinity for POM.
Particle-attached bacteria have 10–100× higher cell-specific activity and
direct access to POM substrate.

## Equations

For each POM class (R4, R6, R8), the module computes:

```
k_eff = k_decomp × Q10^((T − Tref)/10) × Q10^((T − 32)/3) [high-T cutoff]
      × O2/(O2 + K_O2)                                      [O2 limitation]

dPOM_c/dt  = −k_eff × POM_c                    [carbon loss]
dPOM_n/dt  = −k_eff × POM_n                    [nitrogen loss]
dPOM_p/dt  = −k_eff × POM_p                    [phosphorus loss]
dPOM_s/dt  = −k_eff × POM_s                    [silicate loss, if present]

dO2/dt     = −k_eff × POM_c × ur_O2            [oxygen consumption]
dDIC/dt    = +k_eff × POM_c / 12               [CO2 production]
dNH4/dt    = +k_eff × POM_n                    [ammonium release]
dPO4/dt    = +k_eff × POM_p                    [phosphate release]
dSi/dt     = +k_eff × POM_s                    [silicate release]
dTA/dt     = +k_eff × POM_n − k_eff × POM_p   [alkalinity: +1/NH4, −1/PO4]
```

The temperature function follows ERSEM's standard Q10 formulation with a
high-temperature cutoff at 32°C, consistent with all other ERSEM modules.

## Parameters

| Parameter | Unit | Default | Description |
|-----------|------|:-------:|-------------|
| `k_decomp` | 1/d | 0.02 | Specific decomposition rate at Tref |
| `q10` | - | 2.0 | Q10 temperature coefficient |
| `Tref` | °C | 20.0 | Reference temperature |
| `K_O2` | mmol O2/m³ | 5.0 | Half-saturation O2 for decomposition |
| `ur_O2` | mmol O2/mg C | 0.1 | O2 consumed per carbon decomposed |

### Recommended Values by POM Class

| POM class | k_decomp | Rationale |
|-----------|:--------:|-----------|
| R4 (small, labile) | 0.04–0.10 | Fresh phytoplankton, fecal pellets |
| R6 (medium) | 0.02–0.05 | Marine snow, aged detritus |
| R8 (large, refractory) | 0.004–0.01 | Slow-sinking refractory material |

Values are based on BROM (Yakushev et al. 2017), ERGOM (Neumann et al. 2022),
and CE-QUAL-ICM (Cerco & Cole 1993). See the model comparison document for
details.

## FABM Configuration

Each POM class requires a separate instance in `fabm.yaml`. Couplings to
POM state variables and inorganic pools must be specified explicitly.
Nutrient couplings (N, P, Si, Fe) are optional (`required=.false.`).

### Example: R6 (medium POM with Si)

```yaml
R6_decay:
  model: ersem/pom_decay
  parameters:
    k_decomp: 0.03
    q10: 2.0
    Tref: 20.0
    K_O2: 5.0
    ur_O2: 0.1
  coupling:
    RPc: R6/c           # POM carbon (required)
    RPn: R6/n           # POM nitrogen (optional)
    RPp: R6/p           # POM phosphorus (optional)
    RPs: R6/s           # POM silicate (optional, only if present)
    O2o: O2/o           # dissolved oxygen (required)
    O3c: O3/c           # DIC (required)
    N4n: N4/n           # ammonium (optional)
    N1p: N1/p           # phosphate (optional)
    N5s: N5/s           # silicate (optional)
```

### Example: R4 (small POM without Si)

```yaml
R4_decay:
  model: ersem/pom_decay
  parameters:
    k_decomp: 0.06
  coupling:
    RPc: R4/c
    RPn: R4/n
    RPp: R4/p
    O2o: O2/o
    O3c: O3/c
    N4n: N4/n
    N1p: N1/p
```

Note: Do NOT couple `RPf` or `N7f` (iron) unless `ERSEM_USE_IRON=ON` in
CMakeLists.txt. Similarly, do not couple `RPs`/`N5s` for POM classes that
lack silicate (e.g., R4 with `composition: cnpf`).

## Source Files

| File | Change |
|------|--------|
| `src/pom_decay.F90` | New module (jsasaki 2026-03-19) |
| `src/ersem_model_library.F90` | Register `pom_decay` in factory |
| `src/CMakeLists.txt` | Add `pom_decay.F90` to build |

## Interaction with Existing Processes

The `pom_decay` module operates **in parallel** with:

- **B1 sRP pathway**: Free-living bacteria continue to process POM via
  sRP1R1/sRP2R1/sRP3R1. This pathway remains active and is not modified.
- **POM sinking**: Sinking (via `rm` parameter in pelagic_base) continues
  to transport POM vertically. With pom_decay active, R6 sinking speed
  has negligible effect on mid-water O2 because decomposition dominates.
- **Benthic processing**: H1/H2 benthic bacteria continue to process Q6
  (POM that reaches the bottom). With active water-column decomposition,
  less POM reaches the bottom, potentially reducing SOD.

## Validation (Tokyo Bay)

Tested in Phase 36 (26 cases) of the TB-FVCOM Tokyo Bay simulation:

- k_R6=0.03: Mid-water (8–14 m) August bias improved from +1.83 to +0.73 mg/L
- k_R6=0.05: Mid-water bias −0.51 (first reproduction of mid-water hypoxia)
- R6 sinking speed (5–30 m/d) has no effect when pom_decay is active
- Phase 37 (40 cases) is re-evaluating all ERSEM parameters with pom_decay

## References

- Yakushev, E.V. et al. (2017). Bottom RedOx Model (BROM v.1.1). *GMD* 10:453–482.
- Neumann, T. et al. (2022). Non-Redfieldian carbon model for the Baltic Sea (ERGOM v1.2). *GMD* 15:8473–8540.
- Cerco, C.F. & Cole, T. (1993). Three-dimensional eutrophication model of Chesapeake Bay. *J. Environ. Eng.* 119:1006–1025.
- Butenschön, M. et al. (2016). ERSEM 15.06. *GMD* 9:1293–1339.
