# ERSEM Nutrient Limitation Model: Analysis and Improvement Proposals for Tokyo Bay

## 1. Current ERSEM Nutrient Model Architecture

### 1.1 Hybrid Droop Quota + Affinity System

ERSEM uses a **hybrid system** that combines the Droop cell-quota model with affinity-based nutrient uptake. This is fundamentally different from the simple Monod (external concentration) approach.

**Growth limitation** (internal quota, Droop model):
```
iNp = MIN(1, MAX(0, (qpc - qplc) / (xqcp * qpRPIcX - qplc)))
iNn = MIN(1, MAX(0, (qnc - qnlc) / (xqcn * qnRPIcX - qnlc)))
```
where `qpc` and `qnc` are the current internal P:C and N:C ratios, `qplc` and `qnlc` are the minimum (subsistence) quotas, and `xqcp * qpRPIcX` is the threshold quota for full growth capacity. Growth depends on **internal** nutrient stores, not directly on external concentrations.

**Nutrient uptake** (affinity-based, linear model):
```
rump  = qurp * N1p * c        (max achievable P uptake)
rumn3 = qun3 * N3n * c        (max achievable NO3 uptake)
rumn4 = qun4 * N4n * c        (max achievable NH4 uptake)
```
where `qurp`, `qun3`, `qun4` are affinity parameters (m^3/mg C/d). This is a **first-order linear uptake model** without explicit saturation.

**Uptake regulation** (demand-capped):
```
fN1PIp = MIN(rump, runp + misp)
```
Actual uptake is the minimum of the affinity-based maximum and the demand-driven quota replenishment term. This provides implicit saturation when internal quotas are full.

**Silicate limitation** (only Monod-type in the model):
```
iNs = MIN(1, N5s / (N5s + chs))
```
Only silicate uses the classical Monod half-saturation formulation with external concentrations.

### 1.2 Key Parameters and Their Roles

| Parameter | Units | Role | Ecological interpretation |
|-----------|-------|------|--------------------------|
| `sum` | 1/d | Maximum specific productivity | Maximum growth rate at optimal conditions |
| `qun3`, `qun4` | m^3/mg C/d | N affinity | Uptake efficiency at low external [N] |
| `qurp` | m^3/mg C/d | P affinity | Uptake efficiency at low external [P] |
| `qnlc` | mmol N/mg C | Minimum N:C quota | Subsistence threshold; growth ceases below this |
| `qplc` | mmol P/mg C | Minimum P:C quota | Subsistence threshold for P |
| `xqcn`, `xqcp` | - | Limitation threshold | Quota level at which limitation is relieved |
| `xqn`, `xqp` | - | Maximum quota multiplier | Luxury uptake capacity (relative to Redfield) |
| `snplux` | 1/d | Luxury uptake rate | How fast internal stores fill beyond minimum |
| `chs` | mmol Si/m^3 | Si half-saturation | Monod Ks for silicate (diatoms only) |
| `phim` | mg Chl/mg C | Max Chl:C ratio | Controls chlorophyll synthesis per unit biomass |

### 1.3 Relationship Between ERSEM Affinity and Monod's Ks

The Michaelis-Menten uptake equation at low substrate concentration approximates to:

```
V ≈ (Vmax / Ks) × S
```

This initial slope is the **nutrient affinity** (Healey, 1980; Aksnes & Egge, 1991). Therefore:

```
ERSEM affinity ≈ Vmax / Ks   (per unit biomass carbon)
```

ERSEM's linear uptake model represents the **low-concentration regime** of Monod kinetics. In eutrophic environments like Tokyo Bay (DIN 10-50 mmol N/m^3), the lack of explicit saturation means uptake can be unrealistically high, but the Droop quota system provides implicit saturation as internal stores fill.

## 2. Phytoplankton Competition Theory: r-K Selection and Nutrient Strategies

### 2.1 Functional Trait Framework

Phytoplankton competition for nutrients follows a functional trait framework (Litchman et al., 2007; Litchman & Klausmeier, 2008):

**r-strategists** (bloom-formers, e.g., diatoms like *Skeletonema*):
- High maximum growth rate (Vmax)
- High half-saturation constant (Ks) -- need high nutrients for maximum uptake
- Low to moderate Vmax/Ks ratio → **lower biomass-specific affinity**
- Thrive in pulsed, high-nutrient environments (blooms)

**K-strategists** (steady-state specialists, e.g., picophytoplankton like *Synechococcus*):
- Low maximum growth rate
- Very low Ks -- efficient uptake even at trace nutrient levels
- High Vmax/Ks ratio → **higher biomass-specific affinity**
- Thrive in stable, low-nutrient environments

**Dinoflagellates** (mixotrophic strategists):
- Low growth rate AND high Ks → **lowest affinity**
- Compensate via vertical migration, mixotrophy, toxin production
- Success depends on behavioral strategies, not nutrient uptake efficiency

### 2.2 Competition Outcome in the Droop Model

In the Droop quota framework, the steady-state competitive outcome depends on R* (Tilman, 1982), which for the variable-internal-stores model depends on the interaction of multiple parameters (Grover, 1991):

```
Effective R* ~ (mortality × q_min) / affinity
```

The **competitive index** can be approximated as `q_min / affinity`. A species with low minimum quota AND high affinity has the lowest R* and wins under nutrient-depleted conditions.

**Counter-intuitive aspect**: Neither high affinity alone nor low minimum quota alone guarantees competitive success. Both must be considered together.

### 2.3 Cell Size Scaling Laws

Edwards et al. (2012) established allometric scaling relationships:

| Trait | Scaling with cell volume | Consequence |
|-------|-------------------------|-------------|
| Vmax (absolute) | Increases with size | Larger cells take up more per cell |
| Ks | Increases with size | Larger cells need more ambient nutrient |
| Vmax/Ks (biomass-specific affinity) | **Decreases** with size | Smaller cells are better low-nutrient competitors |
| Maximum growth rate | **Decreases** with size | Smaller cells grow faster per unit biomass |
| Minimum quota (q_min) | Increases with size | Larger cells need more stored nutrient |

Expected affinity ordering:
```
Picophytoplankton > Nanoflagellates > Diatoms > Dinoflagellates
```

## 3. Current Parameter Assessment for Tokyo Bay

### 3.1 Current P1-P4 Configuration

| Parameter | P1 Diatom | P2 Nano | P3 Pico | P4 Micro | Ecological expectation |
|-----------|-----------|---------|---------|----------|----------------------|
| sum (1/d) | 1.375 | 1.625 | 2.0 | 1.1125 | P1 could be higher for eutrophic diatoms |
| qun3 | 0.0025 | 0.004 | 0.006 | 0.002 | Ordering CORRECT: pico > nano > diatom > micro |
| qun4 | 0.0025 | 0.004 | 0.007 | 0.002 | Ordering CORRECT |
| qurp | 0.003 | 0.004 | 0.006 | 0.002 | Ordering CORRECT |
| qnlc | 0.0042 | 0.005 | 0.006 | 0.0042 | **INCORRECT**: P3 highest, should be lowest |
| qplc | 0.0001 | 0.000225 | 0.00035 | 0.0001 | **INCORRECT**: P3 highest, should be lowest |
| phim | 0.06 | 0.025 | 0.015 | 0.045 | P2, P3 too low for eutrophic conditions |
| chs (Si) | 0.2 | - | - | - | Too low; literature 0.5-5.0 mmol/m^3 |

### 3.2 Competitive Index Analysis (q_min / affinity)

| Group | qnlc/qun3 (N) | qplc/qurp (P) | Interpretation |
|-------|----------------|----------------|----------------|
| P1 (diatom) | 1.68 | 0.033 | Moderate N-competitor, strong P-competitor |
| P2 (nano) | 1.25 | 0.056 | Moderate in both |
| P3 (pico) | 1.00 | 0.058 | Best N-competitor, but advantage is weak |
| P4 (micro) | 2.10 | 0.050 | Worst N-competitor |

**Problem**: P3's competitive advantage in nitrogen is too small (1.00 vs P2's 1.25). With corrected minimum quotas, P3 should have a much lower ratio (~0.5), creating stronger differentiation.

### 3.3 Key Issues Identified

1. **Minimum quotas inversely ordered**: P3 (picophytoplankton) has the HIGHEST minimum quotas (qnlc=0.006, qplc=0.00035), contradicting theory. K-strategists should have the LOWEST minimum quotas.

2. **P1 Topt too low**: Current Topt=15°C suppresses diatom growth during summer (25-28°C), but Tokyo Bay diatoms (*Skeletonema*) bloom throughout the year. Topt should be raised to 18-22°C.

3. **phim too low for P2, P3**: Maximum Chl:C ratios of 0.025 and 0.015 are too low for eutrophic conditions. Literature values for nutrient-replete conditions: 0.03-0.06 mg Chl/mg C.

4. **Grazing pressure**: ERSEM defaults were calibrated for the North Sea. Tokyo Bay's hypereutrophic conditions likely require reduced grazing rates.

## 4. Diagnosis: 15x Chlorophyll-a Underestimation

### 4.1 Observed vs Modeled

- **Observed**: ~20 mg Chl/m^3 surface annual mean (inner Tokyo Bay, Kawasaki station)
- **Model**: ~1.35 mg Chl/m^3
- **Factor**: ~15x underestimation

### 4.2 Contributing Factors (Ranked by Estimated Impact)

| Factor | Estimated Impact | Mechanism |
|--------|-----------------|-----------|
| **Grazing too strong** | 2-5x | Z4/Z5/Z6 calibrated for North Sea suppress eutrophic blooms |
| **P1 Topt too low** | 2-3x | Diatom production suppressed in summer when T > 20°C |
| **phim too low for P2/P3** | 1.3-1.5x | Insufficient Chl per unit phytoplankton carbon |
| **sum could be higher** | 1.3-1.5x | Especially P1 in eutrophic conditions |
| **1D model limitations** | 1.5-2x | Missing lateral nutrient inputs (rivers, sewage) |
| **Minimum quota ordering** | 1.1-1.3x | Weak species differentiation reduces total community productivity |

Combined effect: 2.5 × 2 × 1.4 × 1.3 × 1.5 × 1.2 ≈ **16x**, consistent with observed discrepancy.

### 4.3 Typical Chl:C Ratios for Eutrophic Coastal Waters

| Condition | C:Chl (g:g) | Chl:C (g:g) |
|-----------|-------------|-------------|
| Eutrophic estuary, mean | 20-40 | 0.025-0.050 |
| Summer high light | 50-100 | 0.010-0.020 |
| Spring bloom | 20-30 | 0.033-0.050 |
| Deep/low light | 15-25 | 0.040-0.067 |

Sources: Jakobsen & Markager (2016), Geider et al. (1997, 1998).

## 5. Improvement Proposals

### 5.1 Proposal A: Fix Minimum Quota Ordering (No Code Change)

Correct the ecologically inconsistent minimum quotas in `fabm.yaml`:

| Parameter | P1 Diatom | P2 Nano | P3 Pico | P4 Micro |
|-----------|-----------|---------|---------|----------|
| qnlc (current) | 0.0042 | 0.005 | 0.006 | 0.0042 |
| qnlc (proposed) | 0.0042 | 0.0035 | 0.0025 | 0.0042 |
| qplc (current) | 0.0001 | 0.000225 | 0.00035 | 0.0001 |
| qplc (proposed) | 0.0001 | 0.00012 | 0.00008 | 0.0001 |

**Rationale**: Picophytoplankton have smaller cells with lower metabolic requirements per unit biomass. They can subsist on lower internal nutrient levels (Litchman et al., 2007; Edwards et al., 2012).

### 5.2 Proposal B: Increase phim for Non-Diatom Groups (No Code Change)

| Parameter | P1 Diatom | P2 Nano | P3 Pico | P4 Micro |
|-----------|-----------|---------|---------|----------|
| phim (current) | 0.06 | 0.025 | 0.015 | 0.045 |
| phim (proposed) | 0.06 | 0.04 | 0.03 | 0.05 |

**Rationale**: In eutrophic waters, nutrient-replete cells maintain higher Chl:C ratios. The current P2/P3 values are more appropriate for oligotrophic conditions (Geider et al., 1998).

### 5.3 Proposal C: Adjust Grazing Parameters (No Code Change)

Reduce top-down control for Tokyo Bay's hypereutrophic conditions:

| Parameter | Current | Proposed | Group |
|-----------|---------|----------|-------|
| Z5 sum | 1.25 | 0.8-1.0 | Microzooplankton max uptake rate |
| Z5 pu | 0.5 | 0.3-0.4 | Microzooplankton assimilation efficiency |
| Z4 sum | 1.0 | 0.8 | Mesozooplankton max uptake rate |
| Z4 pu | 0.6 | 0.4-0.5 | Mesozooplankton assimilation efficiency |

**Rationale**: In hypereutrophic bays, phytoplankton growth regularly outpaces zooplankton grazing during bloom events. North Sea-calibrated grazing rates are too strong for Tokyo Bay conditions.

### 5.4 Proposal D: Widen CTMI Temperature Ranges for P1 (No Code Change)

| Parameter | Current | Proposed | Rationale |
|-----------|---------|----------|-----------|
| P1 Topt | 15.0 | 18-22 | Tokyo Bay diatoms bloom year-round; 15°C too restrictive |
| P1 Tmax | 30.0 | 32-33 | Allow some growth at summer peak temperatures |
| P4 Tmin | 10.0 | 12-14 | Narrow range forces stronger seasonality |

### 5.5 Proposal E: Increase Silicate Half-Saturation for P1 (No Code Change)

| Parameter | Current | Proposed | Rationale |
|-----------|---------|----------|-----------|
| P1 chs | 0.2 | 1.0-3.0 | Literature Ks(Si) for *Skeletonema*: 1-3 mmol/m^3 |

**Rationale**: Current value (0.2 mmol/m^3) means diatoms essentially never experience Si limitation. Raising to 1.0-3.0 creates a realistic silicate-depletion mechanism for bloom termination (Paasche, 1973; Nelson & Treguer, 1992).

### 5.6 Proposal F: Add Optional Monod Saturation to Nutrient Uptake (Code Change)

**Motivation**: ERSEM's linear affinity model has no uptake saturation. In Tokyo Bay's eutrophic conditions (DIN 10-50 mmol/m^3), this can produce unrealistic uptake rates. Adding an optional Michaelis-Menten saturation would:
- Cap uptake rates at high nutrient concentrations
- Allow half-saturation constant (Ks) to be a tunable parameter per species
- Enable explicit r-K strategy differentiation

**Proposed implementation** in `primary_producer.F90`:

Add parameters:
```fortran
integer :: iswNutUptake    ! 1: linear affinity (default), 2: Michaelis-Menten
real(rk) :: VmaxN, KsN     ! Max N uptake rate and half-saturation (Monod)
real(rk) :: VmaxP, KsP     ! Max P uptake rate and half-saturation (Monod)
```

Replace uptake equations:
```fortran
if (self%iswNutUptake == 1) then
   ! Original linear affinity model
   rumn3 = self%qun3 * N3nP * c
else
   ! Michaelis-Menten saturation
   rumn3 = self%VmaxN * N3nP / (self%KsN + N3nP) * c
end if
```

**Assessment**: This is a larger code change and may not be the highest priority. The linear affinity model works adequately when combined with the Droop quota cap. However, explicit Ks would be more intuitive for parameter tuning and directly comparable to literature values.

**Recommendation**: Defer this code change. Focus first on Proposals A-E (parameter tuning only) and evaluate results before modifying the Fortran code.

## 6. Recommended v5 Optuna Tuning Strategy

### 6.1 Parameter Selection

Based on the analysis, the following 14 parameters are recommended for v5:

**Growth rates** (4 parameters):
| Parameter | Instance | Range | Rationale |
|-----------|----------|-------|-----------|
| sum | P1 | 1.0-3.0 | Diatoms: high growth in eutrophic conditions |
| sum | P2 | 1.0-2.5 | Nanoflagellates |
| sum | P3 | 1.0-3.0 | Picophytoplankton |
| sum | P4 | 0.3-1.0 | Dinoflagellates: slow growers |

**Temperature** (2 parameters):
| Parameter | Instance | Range | Rationale |
|-----------|----------|-------|-----------|
| Topt | P1 | 15-25 | Diatom optimal temperature |
| Topt | P4 | 22-30 | Dinoflagellate optimal temperature |

**Nutrient uptake affinities** (4 parameters):
| Parameter | Instance | Range | Rationale |
|-----------|----------|-------|-----------|
| qun3 | P1 | 0.001-0.005 | Diatom nitrate affinity |
| qun3 | P4 | 0.0005-0.003 | Dinoflagellate nitrate affinity |
| qurp | P1 | 0.001-0.006 | Diatom phosphate affinity |
| chs | P1 | 0.2-5.0 | Diatom silicate half-saturation |

**Grazing** (2 parameters):
| Parameter | Instance | Range | Rationale |
|-----------|----------|-------|-----------|
| sum | Z5 | 0.5-1.5 | Microzooplankton grazing rate |
| pu | Z5 | 0.2-0.6 | Microzooplankton assimilation efficiency |

**Ecosystem processes** (2 parameters):
| Parameter | Instance | Range | Rationale |
|-----------|----------|-------|-----------|
| sN4N3 | pel_nit | 0.05-1.0 | Nitrification rate |
| sum | B1 | 1.5-5.0 | Bacterial decomposition rate |

### 6.2 Pre-tuning Fixed Changes

Before Optuna optimization, apply Proposal A (fix minimum quotas) and Proposal B (increase phim for P2/P3) as fixed parameter changes. These are based on clear theoretical grounds and should not be optimized.

### 6.3 Objective Function Design

The v4 Chl-a objective function had issues (RMSE_rel ~1.5 dominated optimization). Recommended changes:

**Option 1: Log-transform for Chl-a**
```python
score_chl = 0.5 * (1 - corr(log(obs+1), log(mod+1))) + 0.5 * RMSE(log(obs+1), log(mod+1)) / std(log(obs+1))
```
This handles the order-of-magnitude scale difference naturally.

**Option 2: Correlation-only for Chl-a**
```python
score_chl = 1 - corr(obs, mod)
```
Focus on temporal pattern rather than absolute magnitude.

**Recommended**: Option 1 (log-transform). This rewards both correct timing and correct magnitude without being dominated by the absolute scale mismatch.

## 7. Literature References

- Aksnes, D.L. & Egge, J.K. (1991). A theoretical model for nutrient uptake in phytoplankton. *Mar. Ecol. Prog. Ser.*, 70, 65-72.
- Butenschon, M. et al. (2016). ERSEM 15.06: a generic model for marine biogeochemistry. *Geosci. Model Dev.*, 9, 1293-1339.
- Droop, M.R. (1968). Vitamin B12 and marine ecology. IV. *J. Mar. Biol. Assoc. UK*, 48, 689-733.
- Droop, M.R. (1973). Some thoughts on nutrient limitation in algae. *J. Phycol.*, 9, 264-272.
- Edwards, K.F. et al. (2012). Allometric scaling and taxonomic variation in nutrient utilization traits. *Limnol. Oceanogr.*, 57(2), 554-566.
- Eppley, R.W. et al. (1969). Half-saturation constants for uptake of nitrate and ammonium by marine phytoplankton. *Limnol. Oceanogr.*, 14, 912-920.
- Flynn, K.J. (2003). Modelling multi-nutrient interactions in phytoplankton. *Prog. Oceanogr.*, 56, 249-279.
- Geider, R.J. et al. (1997). Dynamic model of phytoplankton growth and acclimation. *Mar. Ecol. Prog. Ser.*, 148, 187-200.
- Geider, R.J. et al. (1998). A dynamic regulatory model of phytoplanktonic acclimation. *Limnol. Oceanogr.*, 43, 679-694.
- Grover, J.P. (1991). Resource competition in a variable environment. *Am. Nat.*, 138, 811-835.
- Healey, F.P. (1980). Slope of the Monod equation as an indicator of advantage in nutrient competition. *Microb. Ecol.*, 5, 281-286.
- Jakobsen, H.H. & Markager, S. (2016). Carbon-to-chlorophyll ratio for phytoplankton in temperate coastal waters. *Limnol. Oceanogr.*, 61, 1853-1866.
- Litchman, E. et al. (2007). The role of functional traits and trade-offs in structuring phytoplankton communities. *Ecol. Lett.*, 10, 1170-1181.
- Litchman, E. & Klausmeier, C.A. (2008). Trait-based community ecology of phytoplankton. *Annu. Rev. Ecol. Evol. Syst.*, 39, 615-639.
- Nakada, S. (2021). Phytoplankton species abundance in Tokyo Bay. *Ecol. Res.*
- Nelson, D.M. & Treguer, P. (1992). Role of silicon as a limiting nutrient to Antarctic diatoms. *Mar. Ecol. Prog. Ser.*, 85, 263-272.
- Paasche, E. (1973). Silicon and the ecology of marine plankton diatoms. *Mar. Biol.*, 19, 117-126.
- Smayda, T.J. (1997). Harmful algal blooms: their ecophysiology and general relevance. *Limnol. Oceanogr.*, 42, 1137-1153.
- Tilman, D. (1982). *Resource Competition and Community Structure*. Princeton University Press.
