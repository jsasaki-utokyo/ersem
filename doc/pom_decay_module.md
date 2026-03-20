# Organic Matter Decomposition Redesign: POM Hydrolysis and DOM Oxidation

## Overview

This document describes the redesign of particulate organic matter (POM)
decomposition and dissolved organic matter (DOM) oxidation in ERSEM, motivated
by the inability to reproduce mid-water hypoxia in eutrophic estuaries.

The redesign introduces two new modules and one structural change:

1. **`pom_decay`**: POM hydrolysis (POM → DOM) and direct oxidation (POM → CO2 + O2)
2. **`dom_decay`** (planned): Bacteria-independent DOM oxidation (DOM → CO2 + O2)
3. **B1 sRP removal**: B1 no longer consumes POM directly; B1 consumes only DOM

This follows BROM's hybrid philosophy: particle-attached microbial processes
are parameterized as first-order reactions (implicit bacteria), while
free-living bacteria (B1) remain explicit state variables for DOM consumption.

**Initial implementation:** 2026-03-19 (jsasaki)
**Redesign documented:** 2026-03-20

## Background

### Problem: Mid-Water Hypoxia

Standard ERSEM cannot reproduce mid-water (8–14 m) hypoxia in Tokyo Bay
despite extensive parameter tuning (680+ cases over Phases 26–35).
Layer-by-layer O2 budget diagnosis showed:

- Mid-water O2 consumption: 6.6 mmol O2/m³/d
- Lateral O2 supply: ~6.4 mmol O2/m³/d (96% of consumption)
- Net dO2/dt: −0.008 mg/L/d (model) vs −0.077 mg/L/d (observed)

### Root Cause: Missing POM Decomposition Pathway

Standard ERSEM has no mechanism for POM decomposition in the water column
independent of free-living bacteria (B1). B1's POM uptake via the sRP
pathway gives effective decomposition rates of ~0.001 /d, which is 10–75×
slower than comparable models.

Critically, the sRP pathway in `bacteria_docdyn.F90` does not produce DOM
as an intermediate. Despite the parameter name "remineralisation of substrate
to DOM", the code transfers POM carbon directly into B1 biomass:

```fortran
! bacteria_docdyn.F90 L309: POM uptake goes directly to B1, not to R1/R2/R3
fRPB1c = sugB1 * RPcP * sRPR1 / sutB1
rugB1 = sugB1*(R1cP + R2cP*rR2 + R3cP*rR3) + sum(fRPB1c)
```

This is physically incorrect: bacteria cannot directly ingest particles.
All POM decomposition must proceed through extracellular enzymatic hydrolysis,
producing DOM as an intermediate.

### Comparison with Other Models

| Model | POM decomposition | DOM decomposition | Explicit bacteria |
|-------|-------------------|-------------------|:-----------------:|
| **BROM** | First-order autolysis + oxidation | First-order + bacteria | Yes (4 types) |
| **ERGOM** | First-order mineralization | First-order | No |
| **CE-QUAL-ICM** | First-order (3 G-classes) | First-order | No |
| **Fasham (1990)** | Via bacteria | Via bacteria | Yes |
| **ERSEM (original)** | B1 direct uptake only | B1 only | Yes (B1) |
| **ERSEM (redesign)** | First-order + DOM intermediate | First-order + B1 | Yes (B1) |

The redesigned ERSEM adopts BROM's hybrid approach while retaining ERSEM's
multi-class DOM structure (R1/R2/R3).

See `TB-FVCOM/ersem/docs/pom_decomposition_model_comparison.md` for the
full model comparison with equations and parameter values.

## Design Philosophy

### BROM's Hybrid Approach

BROM (Yakushev et al. 2017) uses two parallel mechanisms for organic matter
decomposition:

1. **Implicit (first-order)**: Represents background microbial activity that
   is always present proportional to substrate concentration. Based on the
   assumption that microbial communities rapidly colonize substrates
   (steady-state assumption; Anderson 2005).
2. **Explicit (bacteria state variable)**: Captures dynamic microbial responses
   to environmental changes (redox transitions, bloom events).

Both mechanisms act on the **same substrate pools**. They are not redundant:
the implicit pathway ensures decomposition occurs even when explicit bacteria
are scarce, while explicit bacteria add dynamic variability.

### Application to ERSEM

ERSEM's B1 (free-living bacteria) serves as the explicit component.
The new modules add the implicit component:

- **`pom_decay`**: Implicit POM decomposition (particle-attached bacteria)
- **`dom_decay`**: Implicit DOM decomposition (background microbial activity)

B1's direct POM uptake (sRP pathway) is removed because:
1. Bacteria cannot directly ingest particles (must hydrolyze to DOM first)
2. POM hydrolysis is now handled by `pom_decay`
3. Keeping both would double-count POM decomposition

## Redesigned Organic Matter Flow

```
                    pom_decay module
                   ┌─────────────────────────────────────────┐
                   │                                         │
POM (R4/R6/R8) ───┤  Autolysis (f_aut × k × POM):           │
                   │    → R1 (G1: 65%)  ─┐                  │
                   │    → R2 (G2: 25%)  ─┤  DOM production   │
                   │    → R3 (G3: 10%)  ─┘  (no O2 consumed) │
                   │                                         │
                   │  Direct oxidation ((1-f_aut) × k × POM):│
                   │    → CO2 + O2 consumed                  │
                   │    → NH4, PO4, Si released              │
                   └─────────────────────────────────────────┘

                    dom_decay module (planned)
                   ┌─────────────────────────────────────────┐
DOM (R1/R2/R3) ───┤  First-order oxidation (k_ox × DOM):     │
                   │    → CO2 + O2 consumed                  │
                   │    → NH4, PO4 released                  │
                   └─────────────────────────────────────────┘

                    B1 bacteria (existing, modified)
                   ┌─────────────────────────────────────────┐
DOM (R1/R2/R3) ───┤  B1 uptake (sR1, rR2, rR3):              │
                   │    → B1 biomass → respiration → O2      │
                   │  (sRP pathway REMOVED: B1 no longer      │
                   │   consumes POM directly)                 │
                   └─────────────────────────────────────────┘

                    Benthic system (existing, unchanged)
                   ┌─────────────────────────────────────────┐
POM settling ─────┤  → Q6 (benthic POM)                      │
                   │  → H1/H2 (benthic bacteria) → O2/SOD    │
                   └─────────────────────────────────────────┘
```

### Role Assignment

| Process | Module | Substrate | Products |
|---------|--------|-----------|----------|
| POM hydrolysis (autolysis) | `pom_decay` | POM | R1 + R2 + R3 |
| POM surface oxidation | `pom_decay` | POM | CO2 + O2 consumed |
| DOM background oxidation | `dom_decay` | R1, R2, R3 | CO2 + O2 consumed |
| DOM uptake by free-living bacteria | B1 (existing) | R1, R2, R3 | B1 biomass → respiration |
| POM settling to benthos | Sinking (existing) | POM | Q6 |
| Benthic decomposition | H1/H2 (existing) | Q6 | SOD |

### Key Design Decisions

1. **B1 sRP pathway is removed** (sRP1R1 = sRP2R1 = sRP3R1 = 0).
   B1 consumes DOM only. POM hydrolysis is handled by `pom_decay`.

2. **POM autolysis produces R1/R2/R3** following DiToro's G-model fractions
   (G1=65%, G2=25%, G3=10%). This ensures the reactivity continuum is
   preserved in the DOM pool.

3. **dom_decay and B1 act on the same DOM pools simultaneously**.
   This is identical to BROM's design where K_DON_ox and explicit bacteria
   consume the same DOM. The two pathways are not redundant: when B1 is
   abundant (surface), B1 dominates; when B1 is scarce (mid-water, deep),
   dom_decay provides the baseline decomposition.

4. **R3 (semi-refractory DOM) is preserved** for COD environmental assessment.
   Both dom_decay (with low k_R3_ox) and B1 (with rR3=0.0025) process R3
   slowly, allowing long-term accumulation.

## Equations

### POM Decomposition (`pom_decay`)

For each POM class (R4, R6, R8):

```
k_eff = k_decomp × f(T) × f(O2)
f(T)  = max(0, Q10^((T-Tref)/10) - Q10^((T-32)/3))
f(O2) = O2 / (O2 + K_O2)

Autolysis pathway (fraction f_aut of total decomposition):
  dPOM_c/dt -= k_eff × f_aut × POM_c
  dR1_c/dt  += k_eff × f_aut × POM_c × f_G1        [G1 fraction → R1]
  dR2_c/dt  += k_eff × f_aut × POM_c × f_G2        [G2 fraction → R2]
  dR3_c/dt  += k_eff × f_aut × POM_c × f_G3        [G3 fraction → R3]
  dNH4/dt   += k_eff × f_aut × POM_n               [nutrient release]
  dPO4/dt   += k_eff × f_aut × POM_p
  dTA/dt    += k_eff × f_aut × (POM_n - POM_p)     [alkalinity]
  (NO O2 consumption — O2 is consumed later when DOM is oxidized)

Direct oxidation pathway (fraction 1-f_aut):
  dPOM_c/dt -= k_eff × (1-f_aut) × POM_c
  dO2/dt    -= k_eff × (1-f_aut) × POM_c × ur_O2   [O2 consumption]
  dDIC/dt   += k_eff × (1-f_aut) × POM_c / CMass   [CO2 production]
  dNH4/dt   += k_eff × (1-f_aut) × POM_n
  dPO4/dt   += k_eff × (1-f_aut) × POM_p
  dTA/dt    += k_eff × (1-f_aut) × (POM_n - POM_p)
```

### DOM Oxidation (`dom_decay`, planned)

For each DOM class (R1, R2, R3):

```
k_ox_eff = k_ox × f(T) × f(O2)

dDOM_c/dt -= k_ox_eff × DOM_c
dO2/dt    -= k_ox_eff × DOM_c × ur_O2
dDIC/dt   += k_ox_eff × DOM_c / CMass
dTA/dt    += k_ox_eff × (DOM_n - DOM_p)    [if N,P tracked]
```

Note: R1 has C, N, P composition; R2 and R3 have C only. For R2/R3,
nutrient release follows Redfield stoichiometry applied to the C flux.

## Parameters

### POM Decomposition (`pom_decay`)

| Parameter | Unit | Default | Description |
|-----------|------|:-------:|-------------|
| `k_decomp` | 1/d | 0.03 | Total POM decomposition rate at Tref |
| `f_aut` | - | 0.9 | Fraction routed to DOM (autolysis) |
| `f_G1` | - | 0.65 | Labile fraction of autolysis products (→R1) |
| `f_G2` | - | 0.25 | Semi-labile fraction (→R2) |
| `f_G3` | - | 0.10 | Semi-refractory fraction (→R3) |
| `q10` | - | 2.0 | Q10 temperature coefficient |
| `Tref` | °C | 20.0 | Reference temperature |
| `K_O2` | mmol O2/m³ | 5.0 | O2 half-saturation for decomposition |
| `ur_O2` | mmol O2/mg C | 0.1 | O2 consumed per carbon (direct oxidation only) |

f_aut=0.9 is based on BROM where autolysis(0.15/d)/total(0.16/d) ≈ 0.94.

f_G1/f_G2/f_G3 values are based on DiToro (2001) G-model fractions for
estuarine organic matter.

### DOM Oxidation (`dom_decay`, planned)

| Parameter | Unit | Default | Description |
|-----------|------|:-------:|-------------|
| `k_ox` | 1/d | varies | First-order oxidation rate at Tref |

Suggested values per DOM class:

| DOM class | k_ox (/d) | Rationale |
|-----------|:---------:|-----------|
| R1 (labile) | 0.05–0.10 | BROM K_DON_ox=0.10; fast turnover |
| R2 (semi-labile) | 0.005–0.02 | Slower; R2 = semi-labile by definition |
| R3 (semi-refractory) | 0.0005–0.002 | Very slow; preserves long-term accumulation |

### B1 Modification

| Parameter | Original | Redesign | Reason |
|-----------|:--------:|:--------:|--------|
| sRP1R1 (R4 uptake) | 0.011 | **0.0** | POM hydrolysis handled by pom_decay |
| sRP2R1 (R6 uptake) | 0.020 | **0.0** | Same |
| sRP3R1 (R8 uptake) | 0.002 | **0.0** | Same |
| sR1, rR2, rR3 | unchanged | unchanged | B1 continues DOM consumption |

## Validation Status

### Phase 36 (initial pom_decay, direct oxidation only)

The initial implementation (POM → CO2 + O2, no DOM intermediate) was tested
in Phase 36 (26 cases) and Phase 37 (48 cases):

- k_R6=0.03: Mid-water (8–14 m) August bias improved from +1.83 to +0.73 mg/L
- k_R6=0.05: Mid-water bias −0.51 (first reproduction of mid-water hypoxia)
- R6 sinking speed (5–30 m/d) had no effect when pom_decay was active
- Phase 37 showed P5 Tmax fix remains essential even with pom_decay
- Best configuration: P37_G1 (k=0.035, OBC100%, P5fix, defaults, KGE=0.500)

### Issue Identified: Q6 Depletion

With direct-oxidation-only pom_decay, bottom POM dropped to zero and Q6
(benthic organic matter) decreased from 33,000 to 6,200 mg C/m².
This caused shallow-station bottom DO to increase (worsened reproduction)
because sediment oxygen demand (SOD) collapsed.

The redesign with DOM intermediate (autolysis pathway) and adjusted POM
sinking speeds is expected to address this by:
1. Routing most POM → DOM (which doesn't directly affect bottom POM flux)
2. Only 10% (direct oxidation) consumes POM with O2 in the water column
3. POM continues sinking while DOM is produced; balance between
   decomposition and sinking determines bottom POM flux

### Pending Validation

The full redesign (pom_decay with autolysis + dom_decay + B1 sRP removal)
has not yet been implemented or tested. Phase 38 will validate the
complete redesigned system.

## Source Files

### Current (Phase 36–37)

| File | Status | Description |
|------|--------|-------------|
| `src/pom_decay.F90` | Implemented | Direct oxidation only (to be updated) |
| `src/ersem_model_library.F90` | Modified | Factory registration |
| `src/CMakeLists.txt` | Modified | Build list |

### Planned (Redesign)

| File | Status | Description |
|------|--------|-------------|
| `src/pom_decay.F90` | To update | Add autolysis pathway + f_G1/G2/G3 |
| `src/dom_decay.F90` | To create | Bacteria-independent DOM oxidation |
| `src/ersem_model_library.F90` | To update | Register dom_decay |
| `src/bacteria_docdyn.F90` | Config only | Set sRP1R1=sRP2R1=sRP3R1=0 in YAML |

Note: B1's sRP removal does not require source code changes. Setting
sRP values to 0.0 in `fabm.yaml` is sufficient.

## References

- Yakushev, E.V. et al. (2017). Bottom RedOx Model (BROM v.1.1). *GMD* 10:453–482.
- DiToro, D.M. (2001). *Sediment Flux Modeling*. Wiley-Interscience.
- Fasham, M.J.R. et al. (1990). A nitrogen-based model of plankton dynamics
  in the oceanic mixed layer. *J. Mar. Res.* 48:591–639.
- Anderson, T.R. (2005). Plankton functional type modelling: running before
  we can walk? *J. Plankton Res.* 27:1073–1081.
- Ogura, N. et al. (2003). Decomposition of phytoplankton in seawater.
  *J. Oceanogr.* 59:417–426.
- Westrich, J.T. & Berner, R.A. (1984). The role of sedimentary organic matter
  in bacterial sulfate reduction. *Limnol. Oceanogr.* 29:236–249.
- Neumann, T. et al. (2022). Non-Redfieldian carbon model (ERGOM v1.2). *GMD* 15:8473–8540.
- Cerco, C.F. & Cole, T. (1993). Three-dimensional eutrophication model of
  Chesapeake Bay. *J. Environ. Eng.* 119:1006–1025.
- Butenschön, M. et al. (2016). ERSEM 15.06. *GMD* 9:1293–1339.
