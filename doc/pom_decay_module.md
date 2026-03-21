# Organic Matter Decomposition Redesign: POM Hydrolysis and DOM Oxidation

## Overview

This document describes the redesign of particulate organic matter (POM)
decomposition and dissolved organic matter (DOM) oxidation in ERSEM, motivated
by the inability to reproduce mid-water hypoxia in eutrophic estuaries.

The redesign introduces two new modules and one structural change:

1. **`pom_decay`** (v3): POM hydrolysis (POM → DOM) and direct oxidation (POM → CO2 + O2)
   with **independent rate constants** for each pathway
2. **`dom_decay`**: Bacteria-independent DOM oxidation (DOM → CO2 + O2)
3. **B1 sRP removal**: B1 no longer consumes POM directly; B1 consumes only DOM

This follows BROM's hybrid philosophy: particle-attached microbial processes
are parameterized as first-order reactions (implicit bacteria), while
free-living bacteria (B1) remain explicit state variables for DOM consumption.

**History:**
- 2026-03-19 (jsasaki): Initial implementation (direct oxidation only, v1)
- 2026-03-20: Redesign with autolysis pathway and G-model DOM split (v2)
- 2026-03-21: Separate independent rate constants k_autolysis and k_oxidation (v3)

## Design Philosophy

### BROM's Hybrid Approach

BROM (Yakushev et al. 2017) uses independent first-order rate constants:

| Pathway | Rate constant | Half-life | Products |
|---------|:------------:|:---------:|----------|
| Autolysis (POM → DOM) | K_PON_DON = 0.15 /d | 4.6 days | DOM |
| Direct oxidation (POM → CO2) | K_PON_ox = 0.01 /d | 69 days | CO2 + O2 consumed |
| DOM oxidation (DOM → CO2) | K_DON_ox = 0.10 /d | 7 days | CO2 + O2 consumed |

These are **independent** rates, not fractions of a single rate. The autolysis
and direct oxidation have fundamentally different timescales.

### Application to ERSEM

- **`pom_decay`**: Two independent rate constants (k_autolysis, k_oxidation)
- **`dom_decay`**: Background DOM oxidation (equivalent to BROM's K_DON_ox)
- **B1**: Consumes DOM only (sRP = 0). POM hydrolysis handled by pom_decay.

### POM Class-Specific Parameters

ERSEMs three POM classes (R4, R6, R8) represent different particle types
with different physical and biochemical properties:

| POM | Type | Sinking | k_autolysis | Rationale |
|:---:|------|:-------:|:-----------:|-----------|
| R4 | Small detritus | Slow | Higher (2× R6) | High surface/volume ratio |
| R6 | Fecal pellets, aggregates | **Fast** | Base rate | Dense, mineral-ballasted |
| R8 | Diatom shells, refractory | Slow | Lower (R6/3) | Resistant to hydrolysis |

Sinking speeds follow Stokes' law relationships:
- R4 ≈ R6 × 0.1 (small, low density)
- R8 ≈ R6 × 0.2 (large but porous, low effective density)
- R6 is fastest due to high particle density (fecal pellets)

## Equations

### POM Decomposition (`pom_decay` v3)

For each POM class, two independent pathways:

```
f(T)  = max(0, Q10^((T-Tref)/10) - Q10^((T-32)/3))
f(O2) = O2 / (O2 + K_O2)

Autolysis (POM → DOM):
  rate_aut = k_autolysis × f(T) × f(O2) × POM
  dPOM/dt -= rate_aut
  dR1/dt  += rate_aut × f_G1       [65% → R1 labile DOM]
  dR2/dt  += rate_aut × f_G2       [25% → R2 semi-labile]
  dR3/dt  += rate_aut × f_G3       [10% → R3 semi-refractory]
  (NO O2 consumed)

Direct oxidation (POM → CO2 + O2):
  rate_ox = k_oxidation × f(T) × f(O2) × POM
  dPOM/dt -= rate_ox
  dO2/dt  -= rate_ox × ur_O2
  dDIC/dt += rate_ox / CMass
```

G-model fractions (DiToro 2001): f_G1=0.65, f_G2=0.25, f_G3=0.10.

### DOM Oxidation (`dom_decay`)

```
dDOM/dt = -k_ox × f(T) × f(O2) × DOM
dO2/dt  = -k_ox × f(T) × f(O2) × DOM × ur_O2
```

Both dom_decay and B1 act on the same R1/R2/R3 pools simultaneously.

## Parameters

### POM Decomposition (`pom_decay` v3)

| Parameter | Unit | Default | Description |
|-----------|------|:-------:|-------------|
| `k_autolysis` | 1/d | 0.15 | Hydrolysis rate: POM → DOM (BROM: 0.15) |
| `k_oxidation` | 1/d | 0.01 | Direct oxidation: POM → CO2+O2 (BROM: 0.01) |
| `f_G1` | - | 0.65 | Labile fraction → R1 (DiToro G1) |
| `f_G2` | - | 0.25 | Semi-labile fraction → R2 (DiToro G2) |
| `f_G3` | - | 0.10 | Semi-refractory fraction → R3 (DiToro G3) |
| `q10` | - | 2.0 | Q10 temperature coefficient |
| `Tref` | °C | 20.0 | Reference temperature |
| `K_O2` | mmol O2/m³ | 5.0 | O2 half-saturation |
| `ur_O2` | mmol O2/mg C | 0.1 | O2 per C (direct oxidation only) |

**Note on defaults:** BROM defaults (k_autolysis=0.15) are too high for
Tokyo Bay. Phase 39b found k_aut_R6=0.02-0.03 is optimal (see Validation).

### Recommended Values per POM Class

| POM | k_autolysis | k_oxidation | Sinking rm |
|:---:|:-----------:|:-----------:|:----------:|
| R4 | R6 × 2 | 0.01 | R6 × 0.1 (min 0.5) |
| R6 | 0.02–0.03 | 0.01 | 10 m/d (base) |
| R8 | R6 / 3 | 0.01 | R6 × 0.2 (min 0.5) |

### DOM Oxidation (`dom_decay`)

| DOM class | k_ox (/d) | Rationale |
|-----------|:---------:|-----------|
| R1 (labile) | 0.005–0.01 | Lower than BROM (0.10) to avoid excess consumption |
| R2 (semi-labile) | k_R1/5 | |
| R3 (semi-refractory) | k_R1/50 | Preserves long-term accumulation |

### B1 Modification

| Parameter | Original | Redesign | Reason |
|-----------|:--------:|:--------:|--------|
| sRP1R1 | 0.011 | **0.0** | POM hydrolysis handled by pom_decay |
| sRP2R1 | 0.020 | **0.0** | Same |
| sRP3R1 | 0.002 | **0.0** | Same |
| sR1, rR2, rR3 | unchanged | unchanged | B1 continues DOM consumption |

## Source Files

| File | Version | Description |
|------|---------|-------------|
| `src/pom_decay.F90` | v3 (2026-03-21) | Independent k_autolysis + k_oxidation |
| `src/dom_decay.F90` | v1 (2026-03-20) | Bacteria-independent DOM oxidation |
| `src/ersem_model_library.F90` | Modified | Factory registration |
| `src/CMakeLists.txt` | Modified | Build list |

## Validation (Tokyo Bay)

### Phase 36 (26 cases): Initial pom_decay v1 (direct oxidation only)

- k_R6=0.03: Mid-water bias improved +1.83 → +0.73 mg/L
- k_R6=0.05: First reproduction of mid-water hypoxia (bias −0.51)
- **Issue**: POM_bot=0, Q6 depleted (1/5 of original)

### Phase 37-38 (114 cases): pom_decay v2 with autolysis + dom_decay

- P5 Tmax fix remains essential
- BROM defaults (k_aut=0.15) cause extreme over-consumption
- k_dom_R1=0.05 (BROM-like) too high; k_dom_R1=0.01 better
- f_aut fraction approach was abandoned (v2 → v3) because autolysis
  and direct oxidation have fundamentally different timescales

### Phase 39/39b/39c (166 cases): pom_decay v3 with independent rates

**Critical finding: POM sinking speed has NO effect on dissolved oxygen.**

R4=1-50, R6=1-70 m/d (including Stokes-consistent R4/R8 linkage) all
produce identical results when k_autolysis is the same. This is because
pom_decay acts on the standing stock of POM at every grid point; the
total DOM production (and thus O2 consumption) depends only on the total
POM present in the domain, not on where it sinks.

**Optimal k_autolysis for R6 (dominant POM class):**

| k_aut_R6 | Mid-water Aug | Bottom Aug | Assessment |
|:--------:|:-------------:|:----------:|------------|
| 0.005 | +0.95 | +0.44 | Too weak |
| 0.01 | +0.88 | +0.35 | Weak |
| **0.02** | **+0.70** | **+0.13** | **Good balance** |
| **0.03** | **+0.31** | **−0.20** | **Best mid-water** |
| 0.05 | −1.08 | −1.40 | Over-consumption |
| 0.15 (BROM) | −2.86 | −2.75 | Severe over-consumption |

k_aut_R6 = 0.02–0.03 (half-life 23–35 days) is optimal for Tokyo Bay.
This is ~5–8× slower than BROM's default (0.15 /d, half-life 4.6 days).

### Remaining Issues (as of Phase 39c)

1. **Shallow station bottom DO**: Chi-1LH and Urayasu bottom DO bias not
   yet evaluated with full postprocessing for the redesigned system.
2. **Q6 depletion**: With pom_decay active, less POM reaches the benthos.
   The impact on SOD and benthic nutrient recycling needs assessment.
3. **POM sinking speed insensitivity**: The finding that sinking speed has
   zero effect is physically unexpected and may indicate a model behavior
   that warrants further investigation (e.g., DOM-mediated feedback loops).

### Next Steps (to be resumed)

1. Run postprocess (Stage 1) for P39b (64 cases) and P39c (30 cases)
   to obtain full validation statistics (KGE, RMSE, bias per station)
2. Select best cases for Stage 2 (full pipeline with plots)
3. Evaluate all 4 stations (Kawasaki, Chi-1LH, Kemigawa, Urayasu)
4. Investigate why POM sinking speed has no effect
5. Test with adjusted POM production ratios (R4:R6:R8 split in
   zooplankton/phytoplankton death terms)

## References

- Yakushev, E.V. et al. (2017). Bottom RedOx Model (BROM v.1.1). *GMD* 10:453–482.
- DiToro, D.M. (2001). *Sediment Flux Modeling*. Wiley-Interscience.
- Fasham, M.J.R. et al. (1990). A nitrogen-based model of plankton dynamics. *J. Mar. Res.* 48:591–639.
- Anderson, T.R. (2005). Plankton functional type modelling. *J. Plankton Res.* 27:1073–1081.
- Ogura, N. et al. (2003). Decomposition of phytoplankton in seawater. *J. Oceanogr.* 59:417–426.
- Westrich, J.T. & Berner, R.A. (1984). The role of sedimentary organic matter. *Limnol. Oceanogr.* 29:236–249.
- Cerco, C.F. & Cole, T. (1993). Chesapeake Bay eutrophication model. *J. Environ. Eng.* 119:1006–1025.
- Butenschön, M. et al. (2016). ERSEM 15.06. *GMD* 9:1293–1339.
