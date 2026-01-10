# Sulfur Cycle Implementation Guide for ERSEM

## Overview

This document provides implementation instructions for adding an explicit sulfur cycle to ERSEM, enabling simulation of hydrogen sulfide (H2S) dynamics and "blue tide" (aoshio) phenomena caused by elemental sulfur (S0) in coastal waters.

## Background

Under anoxic conditions in bottom waters and sediments, sulfate-reducing bacteria convert sulfate ions (SO4^2-) to hydrogen sulfide (H2S). When H2S is transported to oxygenated waters, it is oxidized first to elemental sulfur (S0) and then to sulfate. The intermediate S0 is responsible for the milky blue-green coloration known as "blue tide" observed in eutrophic coastal areas such as Tokyo Bay.

Currently, ERSEM handles this process implicitly through "oxygen debt" (K6), but this approach:
- Does not explicitly track sulfur species
- Cannot model S0 accumulation (blue tide indicator)
- Does not properly connect benthic H2S to the water column

## Design Decisions

1. **Phased K6 transition**:
   - **Phase 1**: Keep K6 in parallel with new H2S for backward compatibility
   - **Phase 2**: Deprecate K6, make H2S the primary mechanism
   - **Phase 3**: Remove K6 completely after validation
2. **Sulfur cycle only** - no iron-sulfur interactions (FeS, FeS2) in Phase 1
3. **Explicit sulfur species** in both water column and sediments: SO4^2-, H2S, S0
4. **Minimal implementation first** (Phase 1), with expansion possible in later phases
5. **Sulfate (SO4) treatment**: While SO4 is nearly conservative in seawater (~28 mM), it is explicitly tracked in the model to maintain mass balance. Changes due to sulfate reduction are small relative to the background concentration but are computed for completeness

## Best Practice: Simplified Implementation Strategy

Given ERSEM's existing 3-layer benthic architecture, the recommended approach is to **leverage existing infrastructure** rather than introducing a continuous BROM-style model.

### Core Principle

```
ERSEM's 3-layer model already handles:
- Dynamic layer boundaries based on O2/NO3 penetration
- Equilibrium concentration profiles within layers
- Diffusive transport between layers and to water column
- Layer collapse under anoxic conditions

→ Add sulfur species to this existing framework
→ Add reaction terms to appropriate layers
→ Let existing transport handle the rest
```

### Essential Components (Must Have)

| Component | Location | Configuration | Purpose |
|-----------|----------|---------------|---------|
| H2S state variable | Benthic + Pelagic | Default (`last_layer: 3`) | Track sulfide in all layers |
| Sulfate reduction | Benthic Layer 3 | - | Primary H2S source |
| H2S oxidation | Benthic Layer 1 + Pelagic | **K_H2S_ox: 5.0/10.0** | Remove H2S when O2 present (increased rates) |
| **Oxic barrier** | Benthic Layer 1 → Pelagic | **K_barrier: 1000, K_barrier_rate: 100** | Block H2S flux when oxic layer exists |
| Transport | Via `benthic_column_dissolved_matter` | - | Move H2S between layers |

### Important Components (Recommended)

| Component | Location | Purpose |
|-----------|----------|---------|
| S0 state variable | Benthic + Pelagic | Blue tide indicator |
| Two-step oxidation | H2S → S0 → SO4 | Capture partial oxidation |
| SO4 state variable | Benthic + Pelagic | Mass balance (optional: can assume constant) |

### Layer 2 H2S-NO3 Oxidation (Implemented - Critical for Winter)

| Component | Location | Purpose |
|-----------|----------|---------|
| **H2S-NO3 oxidation** | Benthic Layer 2 | **Prevents H2S from reaching pelagic when NO3 zone exists** |

**Physical basis**: H2S diffusing from Layer 3 through Layer 2 (denitrification zone) is oxidized by nitrate through chemolithotrophic sulfur oxidation coupled to denitrification:
```
5 H2S + 2 NO3⁻ → 5 S0 + N2 + 4 H2O

Stoichiometry: 0.4 mol NO3 consumed per mol H2S oxidized
```

**Why this is critical**: In winter, Layer 2 (NO3 zone) is thick (~20-26 cm). H2S and NO3 cannot coexist - the presence of NO3 indicates that H2S has been oxidized. Without this mechanism, H2S can pass through Layer 2 and appear in the pelagic even during winter when it should be zero.

**Validation results** (Tokyo Bay test case):
- Before implementation: H2S-NO3 coexistence in 99.8% of timesteps
- After implementation: H2S-NO3 coexistence reduced to 3.1%
- Winter H2S: Reduced from 0.001-0.002 mmol/m³ to 0.000000 mmol/m³

**Rate formulation**:
```fortran
! Convert depth-integrated NO3 to concentration
NO3_conc_2 = NO3_2 / max(D2m - D1m, 0.0001_rk)

! NO3 limitation (Michaelis-Menten)
f_NO3 = NO3_conc_2 / (NO3_conc_2 + K_NO3_half)

! H2S oxidation by NO3 in Layer 2
R_H2S_NO3_ox = K_H2S_NO3_ox * H2S_2 * f_NO3

! ODEs
dH2S_2/dt = -R_H2S_NO3_ox  ! H2S consumed
dNO3_2/dt = -0.4 * R_H2S_NO3_ox  ! NO3 consumed (stoichiometry)
```

**Parameters**:
| Parameter | Value | Units | Description |
|-----------|-------|-------|-------------|
| `K_H2S_NO3_ox` | 50.0 | 1/d | H2S oxidation rate by NO3 (fast reaction) |
| `K_NO3_half` | 10.0 | mmol/m³ | Half-saturation NO3 for H2S-NO3 oxidation |

**YAML configuration**:
```yaml
ben_sulfur:
  parameters:
    K_H2S_NO3_ox: 50.0    # H2S oxidation by NO3 in Layer 2 (1/d)
    K_NO3_half: 10.0      # half-saturation NO3 (mmol/m^3)
  coupling:
    H2S_2: K_H2S/per_layer/h2  # H2S in layer 2
    NO3_2: K4/per_layer/n2     # NO3 in layer 2 (ammonium pool has NO3 in layer 2)
```

**Note**: This mechanism is now implemented as part of Phase 1, not Phase 2 as originally planned. It was found to be essential for correct winter behavior.

### Optional Components (Phase 2+)

| Component | Purpose | Complexity |
|-----------|---------|------------|
| Thiosulfate (S2O3) | BROM-compatible pathway | Medium |
| Iron-sulfur (FeS, FeS2) | Full redox chemistry | High |

### Recommended Implementation Order

```
Step 1: Add H2S to benthic_column_dissolved_matter (composition 'h')
        - Use default last_layer: 3 (do NOT use last_layer: 1)
        - Enables transport through all layers
        - Enables pelagic-benthic exchange

Step 2: Add sulfate reduction in Layer 3
        - Couple to H2 bacteria respiration rate
        - Simple: H2S production = stoich_S_C × remin_rate

Step 3: Add H2S oxidation in Layer 1 (CRITICAL)
        - Use increased rate: K_H2S_ox = 5.0 (10x default)
        - This consumes H2S in Layer 1 before it escapes
        - Rate limited by O2 availability - natural barrier when O2 present

Step 4: Add pelagic H2S oxidation (CRITICAL)
        - Use increased rate: K_H2S_ox = 10.0 (20x default)
        - Rapidly oxidizes any H2S that escapes to water column
        - Ensures H2S < 0.01 mmol/m³ under oxic conditions

Step 5: (Optional) Add S0 intermediate
        - Split oxidation: H2S → S0 → SO4
        - Enables blue tide simulation
```

### Key Simplification: Layer 3 is Always Anoxic

In ERSEM's 3-layer model, Layer 3 is **by definition** below the denitrification zone. Therefore:
- No need for electron acceptor cascade check in Layer 3
- Sulfate reduction is always enabled in Layer 3
- The rate is controlled by organic matter availability (H2 bacteria respiration)

This is simpler than BROM's continuous approach where cascade logic is needed everywhere.

### Implementation Options

Choose based on your modeling needs:

#### Option A: Minimal (H2S only)
**Use case**: Basic sulfur dynamics, oxygen consumption by sulfide oxidation

```
Species: H2S only (SO4 assumed constant at 28 mM)
Reactions:
  - Benthic Layer 3: OM → H2S (sulfate reduction)
  - Benthic Layer 1: H2S + O2 → (removed)
  - Pelagic: H2S + O2 → (removed)

Pros: Simplest implementation, captures H2S toxicity
Cons: No S0 (cannot simulate blue tide color), no sulfur mass balance
Files: 1 new module + modifications to benthic_column_dissolved_matter
```

#### Option B: H2S + S0 (Recommended for Blue Tide)
**Use case**: Blue tide simulation, partial oxidation dynamics

```
Species: H2S, S0 (SO4 assumed constant)
Reactions:
  - Benthic Layer 3: OM → H2S
  - Benthic Layer 1: H2S + O2 → S0; S0 + O2 → (removed)
  - Pelagic: same as Layer 1

Pros: Captures blue tide (S0 accumulation), moderate complexity
Cons: No explicit sulfur mass balance
Files: 2 new modules + modifications to benthic_column_dissolved_matter
```

#### Option C: Full Sulfur Cycle (Current Guide)
**Use case**: Complete sulfur mass balance, future BROM integration

```
Species: SO4, H2S, S0 (all explicit)
Reactions:
  - Benthic Layer 3: SO4 + OM → H2S
  - Benthic Layer 1: H2S + O2 → S0; S0 + O2 → SO4
  - Pelagic: same as Layer 1

Pros: Full mass balance, ready for Phase 2 expansion
Cons: Most complex, requires tracking SO4 (mostly constant)
Files: 2 new modules + modifications to benthic_column_dissolved_matter
```

**Recommendation**: Start with **Option B** (H2S + S0) for blue tide simulation. SO4 can be added later if mass balance is important.

## New State Variables

### Pelagic (Water Column) Variables

| Variable | Units | Description | Notes |
|----------|-------|-------------|-------|
| `N8s_SO4` | mmol S/m^3 | Sulfate ion | Conservative in seawater, ~28 mM |
| `N8s_H2S` | mmol S/m^3 | Hydrogen sulfide | Replaces current N6 |
| `N8s_S0` | mmol S/m^3 | Elemental sulfur | Blue tide indicator |

### Benthic (Sediment) Variables

| Variable | Units | Description | Layers |
|----------|-------|-------------|--------|
| `G2s_SO4` | mmol S/m^2 | Benthic sulfate | All layers (1-3) |
| `G2s_H2S` | mmol S/m^2 | Benthic hydrogen sulfide | All layers (1-3) |
| `G2s_S0` | mmol S/m^2 | Benthic elemental sulfur | All layers (1-3) |

### Variables to Transition (Phased Removal)

| Variable | Phase 1 | Phase 2 | Phase 3 | Reason |
|----------|---------|---------|---------|--------|
| `K6` | Keep (parallel) | Deprecate | Remove | Replaced by explicit H2S |
| `N6` | Keep (parallel) | Deprecate | Remove | Replaced by N8s_H2S |
| `o_deep` (negative values) | Keep | Constrain to ≥0 | Remove | Oxygen should be non-negative; use H2S instead |

**Note**: During Phase 1, both K6/N6 and the new H2S variables coexist. This allows validation against existing behavior before transitioning.

## Comparison with BROM Model

This implementation is based on the BROM (Bottom RedOx Model) sulfur chemistry but adapted for ERSEM's 3-layer benthic structure. Key differences:

### BROM Characteristics
- **Water column model**: BROM operates in a single vertically-resolved domain without explicit benthic layers
- **4 sulfur species**: SO4, S2O3 (thiosulfate), S0, H2S
- **Iron-sulfur coupling**: FeS, FeS2, MnS precipitation/dissolution
- **Thiodenitrification**: H2S oxidation coupled to NO3 reduction

### Phase 1 Simplifications

| Aspect | BROM | ERSEM Phase 1 | Rationale |
|--------|------|---------------|-----------|
| Sulfur species | SO4, S2O3, S0, H2S | SO4, S0, H2S | S2O3 omitted for simplicity |
| Oxidation pathway | H2S → S0 → S2O3 → SO4 | H2S → S0 → SO4 | Direct S0 → SO4 in one step |
| H2S oxidation kinetics | Linear: K × O2 × H2S | Michaelis-Menten for O2 | Prevents excessive rates at high O2 |
| Benthic structure | None (water column) | 3-layer sediment | Adapted to ERSEM architecture |
| Iron coupling | Full Fe-S system | None | Deferred to Phase 3 |
| Thiodenitrification | H2S + NO3 → S0 + N2 | Implemented in Layer 2 | Critical for winter behavior |

### Parameter Correspondence

The following parameters are taken directly from BROM (confirmed from FABM testcases):
- `K_SO4_rd = 0.000005` (1/d) - Sulfate reduction rate
- `K_H2S_ox = 0.5` (1/d) - H2S oxidation rate
- `K_S0_ox = 0.02` (1/d) - S0 oxidation rate
- `stoich_S_C = 0.5` (mol S/mol C) - From BROM stoichiometry: 53 SO4 per 106 C

### Adaptation to ERSEM Benthic Layers

BROM's redox zonation (oxic → suboxic → anoxic) maps to ERSEM's benthic layers:
- **Layer 1** (oxic): H2S and S0 oxidation active
- **Layer 2** (denitrification): Transport only in Phase 1; thiodenitrification in Phase 2
- **Layer 3** (anoxic): Sulfate reduction active

## Chemical Reactions

### Reaction Equations

```
R1: Sulfate Reduction (Layer 3 only in Phase 1, anaerobic)
    SO4^2- + 2CH2O -> H2S + 2HCO3^-

    Note: Thermodynamically, sulfate reduction can occur in Layer 2 when NO3
    is depleted. However, in Phase 1 implementation, R1 is restricted to
    Layer 3 for simplicity. Layer 2 sulfate reduction may be added in Phase 2.

R2: Sulfide Oxidation (Layer 1 and water column, aerobic)
    H2S + 0.5 O2 -> S0 + H2O

R3: Elemental Sulfur Oxidation (Layer 1 and water column, aerobic)
    S0 + 1.5 O2 + H2O -> SO4^2- + 2H+
```

### Rate Formulations

#### Electron Acceptor Cascade Control Functions

Based on BROM model approach using smooth tanh transitions:

```fortran
! Anoxia function (sulfate reduction inhibited by oxygen)
! O2_scale controls transition sharpness (smaller = sharper)
f_anox(O2) = 0.5 * (1.0 - tanh((O2 - O2_threshold) / O2_scale))

! Nitrate inhibition (sulfate reduction inhibited by nitrate)
! NO3_scale controls transition sharpness (smaller = sharper)
f_no_NO3(NO3) = 0.5 * (1.0 - tanh((NO3 - NO3_threshold) / NO3_scale))

! Combined sulfate reduction factor
! f_sulfate_reduction ≈ 1 when both O2 and NO3 are below thresholds
! f_sulfate_reduction ≈ 0 when either O2 or NO3 is above threshold
f_sulfate_reduction = f_anox * f_no_NO3
```

**Behavior of tanh transition**:
- When O2 >> O2_threshold: f_anox → 0 (sulfate reduction suppressed)
- When O2 << O2_threshold: f_anox → 1 (sulfate reduction enabled)
- O2_scale determines the width of the transition zone

#### R1: Sulfate Reduction Rate

```fortran
R_sulfate_red = K_SO4_rd * SO4 / (SO4 + K_SO4_half) * remin_rate * f_sulfate_reduction
```

#### R2: H2S Oxidation Rate

```fortran
R_H2S_ox = K_H2S_ox * H2S * O2 / (O2 + K_O2_half)
```

#### R3: S0 Oxidation Rate

```fortran
R_S0_ox = K_S0_ox * S0 * O2 / (O2 + K_O2_half)
```

### Default Parameter Values (from BROM)

**Important**: All concentrations use **mmol/m³** as the standard unit in ERSEM.
1 mmol/m³ = 1 µmol/L = 1000 µmol/m³

| Parameter | Value | Units | Description |
|-----------|-------|-------|-------------|
| `K_SO4_rd` | 0.000005 | 1/d | Sulfate reduction rate constant |
| `K_SO4_half` | 1.6 | mmol/m³ | Half-saturation for sulfate |
| `K_H2S_ox` | 0.5 | 1/d | H2S oxidation rate constant |
| `K_S0_ox` | 0.02 | 1/d | S0 oxidation rate constant |
| `K_O2_half` | 0.002 | mmol/m³ | Half-saturation for oxygen |
| `stoich_S_C` | 0.5 | mol S/mol C | Stoichiometry of sulfate reduction (53 SO4 : 106 C) |

**Phase 2 parameters** (for Layer 2 sulfate reduction and thiodenitrification):
| Parameter | Value | Units | Description |
|-----------|-------|-------|-------------|
| `O2_threshold` | 0.01 | mmol/m³ | O2 threshold for sulfate reduction |
| `NO3_threshold` | 0.005 | mmol/m³ | NO3 threshold for sulfate reduction |
| `K_hs_no3` | 0.8 | 1/d | Thiodenitrification rate (H2S + NO3) |

**Note on scale parameters**: Setting scale equal to threshold gives a smooth transition where f ≈ 0.12 at concentration = 0 and f ≈ 0.88 at concentration = 2×threshold. Smaller scale values create sharper transitions.

### Stoichiometry

| Reaction | O2 Consumption |
|----------|----------------|
| R1 | 0 (organic matter substitutes) |
| R2 | 0.5 mol O2 / mol H2S |
| R3 | 1.5 mol O2 / mol S0 |

## Integration with 3-Layer Benthic Model

### Dynamic Layer Depth Mechanism

ERSEM's benthic layer depths (D1m, D2m) are **state variables** that dynamically adjust based on chemical equilibrium. This is a key architectural feature that must be understood for proper sulfur cycle integration.

#### How Layer Depths are Controlled

The `benthic_column_dissolved_matter.F90` module calculates equilibrium concentration profiles and relaxes layer depths toward equilibrium:

```fortran
! Equilibrium depth where concentration drops to zero (line 595)
D = -2*sigma*C0/(P+2*P_deep)

! Layer depth relaxes toward equilibrium (line 349)
_SET_BOTTOM_ODE_(self%id_layer, (d_top - Dm(self%last_layer)) / self%relax)
```

Where:
- `sigma` = diffusivity (m²/d)
- `C0` = concentration at layer top (from pelagic)
- `P` = source/sink in current layer (negative for consumption)
- `P_deep` = source/sink in deeper layers

**Key insight**: When pelagic concentration (C0) is HIGH, the equilibrium depth D increases proportionally.

#### Layer Control Mapping

| Module | `last_layer` | Controls | Behavior |
|--------|--------------|----------|----------|
| G2 (oxygen) | 1 | D1m | D1m expands when bottom water O2 is high |
| K3 (nitrate) | 2 | D2m | D2m expands when NO3 supply is high |

#### Observed Correlations (Tokyo Bay Test Case)

Analysis of model output confirms the dynamic layer mechanism is working:

| Correlation | Value | Interpretation |
|-------------|-------|----------------|
| O2_bottom vs D1m | **0.39** | Moderate positive - D1m responds to O2 |
| O2_bottom vs D2m | **0.39** | Moderate positive - D2m also responds |
| D1m vs nitrification | **0.97** | Very strong - larger D1m → more nitrification |

**Layer depths at different O2 conditions**:
- High O2 (>200 mmol/m³): D1m = 0.0028 m (**40% larger** than mean)
- Overall mean: D1m = 0.002 m

This confirms that oxygen penetration depth increases when bottom water O2 is high, leading to enhanced nitrification and NO3 production in Layer 1.

### Sulfur Cycle Architectural Issue

**Problem Identified**: The current implementation allows H2S to diffuse directly from Layer 3 to the pelagic without being fully oxidized at the oxic/anoxic interface.

#### Root Cause Analysis

The nitrogen cycle handles this correctly through the `last_layer` mechanism:
- **G2 (oxygen)**: `last_layer: 1` → O2 drops to zero at Layer 1 bottom
- **K3 (nitrate)**: `last_layer: 2` → NO3 drops to zero at Layer 2 bottom

However, the original sulfur implementation used:
- **K_H2S**: `last_layer: 3` (default) → H2S has non-zero profile throughout all layers

This means:
1. H2S produced in Layer 3 diffuses upward through all layers
2. The oxidation in Layer 1 only acts on the H2S **pool** (H2S_1), not the diffusing **flux**
3. When diffusion rate exceeds oxidation rate, H2S escapes to pelagic

#### Quantitative Evidence

From Tokyo Bay simulation with ben_sulfur module:
- H2S production rate: 3.0 mmol S/m²/d (in Layer 3)
- H2S oxidation rate: 0.05 mmol S/m²/d (in Layer 1)
- **Imbalance ratio: ~60x** (production >> oxidation)
- Result: H2S appears in pelagic even when O2 > 200 mmol/m³

#### Solution Approaches

**Approach 1: `last_layer: 1` for H2S (NOT RECOMMENDED)**

Initially, it seemed logical to configure K_H2S with `last_layer: 1` (like G2/oxygen) to create a chemical barrier. However, **this approach causes a conflict**:

- Both G2 (oxygen) and K_H2S would try to control D1m via their ODE contributions
- These contributions are summed, causing D1m to grow to unrealistic values (~150 mm instead of ~2 mm)
- The layer depth mechanism is designed for **one constituent per layer boundary**:
  - G2 controls D1m
  - K3 controls D2m
  - K_H2S should NOT also control D1m

**Approach 2: Increased Oxidation Rates (PARTIALLY EFFECTIVE)**

Increasing oxidation rates helps consume H2S but may not be sufficient in all conditions:

```yaml
# Benthic sulfur cycle - moderate benthic oxidation
ben_sulfur:
  parameters:
    K_H2S_ox: 5.0          # 10x increase from default 0.5

# Pelagic sulfur cycle - rapid water column oxidation
pel_sulfur:
  parameters:
    K_H2S_ox: 10.0         # 20x increase from default 0.5
```

**Limitation**: This approach relies on oxidation of H2S *after* it enters the bottom water. During fall/winter when H2S flux is high, some H2S can still escape to the pelagic even with high oxidation rates.

**Approach 3: Oxic Barrier Mechanism (RECOMMENDED)**

The most physically correct solution is to implement an **oxic barrier** that prevents H2S from passing through the oxic layer. This is implemented in `benthic_sulfur_cycle.F90`:

```yaml
ben_sulfur:
  parameters:
    K_barrier: 1000.0      # Oxic barrier effectiveness (1/m)
    K_barrier_rate: 100.0  # Rate of barrier H2S oxidation (1/d)
```

**Physical basis**: When H2S diffuses upward through an oxic layer, it is rapidly oxidized by O2 with metal catalysis (Fe, Mn). The oxidation rate is so fast (half-life of minutes to hours) that effectively no H2S can pass through even a thin oxic layer.

**Mathematical formulation**:
```
f_barrier = 1 - exp(-K_barrier * D1m * f_O2_pel)
R_barrier_ox = K_barrier_rate * H2S_pel * f_barrier
```

Where:
- `K_barrier` controls how effective the oxic layer is at blocking H2S
- `D1m` is the oxic layer depth (m)
- `f_O2_pel` is oxygen limitation factor based on bottom water O2
- `R_barrier_ox` is the rate at which bottom-water H2S is oxidized by the barrier

**Example**: With K_barrier = 1000 (1/m):
- 2mm oxic layer: f_barrier = 0.86 (86% blocked)
- 5mm oxic layer: f_barrier = 0.993 (99.3% blocked)

**Effect of oxic barrier**:
1. When oxic layer exists (D1m > ~1mm) and O2 is present, H2S flux is effectively blocked
2. During anoxic events when D1m collapses, barrier weakens and H2S can escape
3. Provides physically correct behavior: H2S only enters pelagic under truly anoxic conditions
4. Combined with oxidation rates, achieves H2S ≈ 0 under normal oxic conditions

### Layer-Specific Process Assignment

```
Water Column:
  Variables: N8s_SO4, N8s_H2S, N8s_S0, O2
  Reactions: R2 (H2S oxidation), R3 (S0 oxidation)

============ Sediment-Water Interface ============

Layer 1 (Oxic): 0 to D1m
  Variables: G2s_SO4_1, G2s_H2S_1, G2s_S0_1, G2o_1
  Reactions: R2 (H2S oxidation), R3 (S0 oxidation)
  Note: H2S mostly oxidized here under normal conditions

Layer 2 (Denitrification): D1m to D2m
  Variables: G2s_SO4_2, G2s_H2S_2, G2s_S0_2, G2o_2(~0), K3n_2
  Reactions: H2S-NO3 oxidation (chemolithotrophic denitrification)
  Note: H2S diffusing from Layer 3 is oxidized by NO3, preventing H2S-NO3 coexistence

Layer 3 (Anoxic): D2m to d_tot
  Variables: G2s_SO4_3, G2s_H2S_3, G2s_S0_3
  Reactions: R1 (sulfate reduction) - PRIMARY H2S SOURCE
  Note: Linked to H2 (anaerobic bacteria) respiration rate
```

**Phase 1 Simplification**: Sulfate reduction is implemented only in Layer 3, where conditions are guaranteed to be anoxic with depleted nitrate. This avoids complexity in handling the transition zone in Layer 2.

**ERSEM constraint (Layer 2)**: Layer 2 exists only while NO3 penetrates (D2m > D1m). When NO3 is depleted, D2m → D1m and Layer 2 collapses, leaving no volume for Layer-2 reactions. Therefore sulfate reduction is strictly a Layer-3 process in ERSEM and should not be implemented in Layer 2.

### Transport Mechanism

Use the existing `benthic_column_dissolved_matter` transport scheme:

- Equilibrium concentration profile calculation through all layers
- Diffusive transport based on layer-specific diffusivities (EDZ_1, EDZ_2, EDZ_3)
- Pelagic-benthic exchange at sediment-water interface

**Important**: H2S produced in Layer 3 diffuses through Layer 2 and Layer 1 before reaching the water column. It does NOT bypass intermediate layers. The existing transport code handles this correctly.

### S0 Settling and Surface Benthic Fate (Current Implementation)

- **Diffusive exchange**: Dissolved/colloidal S0 exchanges with the sediment via `benthic_column_dissolved_matter` (pelagic top boundary condition).
- **Settling term**: The pelagic S0 settling term (`K_S0_sink` in `sulfur_cycle.F90`) is implemented as a sink in the water column; it does **not** add S0 to the benthic pools.
- **Benthic surface S0**: In `benthic_sulfur_cycle.F90`, S0 in Layer 1 is produced locally from H2S oxidation and removed by S0 oxidation and burial. It can accumulate transiently if production exceeds removal, but there is no explicit long-term storage term.
- **If deposition is desired**: Add an explicit bottom-exchange/source term that routes some or all of the pelagic S0 settling flux into `G2s_S0_1` (or an equivalent benthic S0 state).

### Oxygen Budget Integration

The sulfur cycle modules modify oxygen through ODE contributions. This integrates with ERSEM's existing oxygen model as follows:

**FABM's ODE accumulation**: Multiple modules can contribute to the same state variable's rate of change. FABM accumulates all `_SET_ODE_` contributions before time integration.

```
Total dO2/dt = (existing ERSEM O2 terms)
             + (sulfur_cycle contribution: -0.5*R_H2S_ox - 1.5*R_S0_ox)
             + (benthic_sulfur_cycle contribution to G2o)
```

**Oxygen non-negativity**:
- Phase 1: Rely on existing ERSEM oxygen handling (may allow negative values as "oxygen debt")
- Phase 2+: Enforce O2 ≥ 0 constraint; negative demand is implicitly handled by H2S accumulation

**Consistency check**: The sulfur module should not consume more O2 than available. The Michaelis-Menten formulation `O2/(O2+K_O2_half)` naturally reduces reaction rates as O2 approaches zero, preventing excessive consumption.

### Behavior Under Bottom Water Anoxia

#### BROM vs ERSEM: Key Architectural Difference

**BROM (Continuous)**: Designed for vertically-resolved water+sediment continuum. Redox boundaries shift smoothly based on O2/NO3 concentrations. No discrete layers - reactions occur at any depth based on local conditions.

**ERSEM (3-Layer Discrete)**: Layer boundaries (D1m, D2m) are **state variables** that adjust dynamically based on oxygen/nitrate penetration depths calculated in `benthic_column_dissolved_matter.F90`. The key mechanism:

```fortran
! From benthic_column_dissolved_matter.F90:
! Equilibrium concentration at sediment surface
c_top = c_pel + cmix * sms   ! c_pel = bottom water conc, sms = net benthic production

! Layer depth relaxes towards equilibrium (where constituent drops to zero)
_SET_BOTTOM_ODE_(self%id_layer, (d_top + max(self%minD, H_eq) - D1m) * self%relax)
```

When bottom water O2 decreases:
1. The equilibrium oxygen profile shifts upward
2. D1m (oxic layer depth) relaxes toward minD (typically 0.0001 m)
3. Layer 1 effectively collapses, removing the oxidation barrier

#### Schematic: Normal vs Anoxic Conditions (with `last_layer: 1` for H2S)

```
NORMAL (Oxic bottom water)        ANOXIC (O2 → 0 in bottom water)

Water: O2 = 200 mmol/m³           Water: O2 ≈ 0 (anoxic event)
       ↓ O2 diffuses down                ↓ No O2 supply to sediment
       H2S ≈ 0 (blocked by barrier)      H2S escapes to water!
═══════════════════════════       ═══════════════════════════
Layer 1: D1m ≈ 2.8 mm             Layer 1: D1m → minD (0.1 mm)
  - O2 present in pore water        - O2 ≈ 0 (no supply from above)
  - H2S → 0 at bottom (barrier!)    - Barrier collapses with layer
  - H2S flux BLOCKED                - H2S passes through!
───────────────────────────       ───────────────────────────
Layer 2: D2m - D1m ≈ 43 mm        Layer 2: Expands
  - Denitrification zone            - NO3 may still be present initially
  - H2S = 0 (below barrier)         - H2S gradient forms
  - Transport only (Phase 1)        - Thiodenitrification possible (Phase 2)
───────────────────────────       ───────────────────────────
Layer 3: d_tot - D2m ≈ 250 mm     Layer 3: Dominates column
  - Sulfate reduction active        - Enhanced H2S production
  - H2S produced and trapped        - H2S diffuses directly to water
  - No escape when barrier exists
═══════════════════════════       ═══════════════════════════

KEY: With last_layer: 1 for H2S, the equilibrium profile forces H2S → 0
at the Layer 1 bottom, creating an effective oxidation barrier. H2S can
only escape to pelagic when D1m collapses during anoxia (D1m → minD).
```

#### H2S Flux Behavior with `last_layer: 1`

With `last_layer: 1` configured for H2S, the `benthic_column_dissolved_matter` module calculates equilibrium profiles where H2S concentration drops to zero at the bottom of Layer 1. The benthic-pelagic flux is determined by this equilibrium profile, not simple bulk diffusion.

**Under normal (oxic) conditions**:
- D1m expands (e.g., 2.8 mm at high O2)
- H2S equilibrium profile → 0 at Layer 1 bottom
- Benthic-pelagic H2S flux ≈ 0 (chemical barrier effective)
- H2S produced in Layer 3 is "trapped" below the barrier

**Under anoxic conditions**:
- D1m collapses toward minD (0.1 mm)
- Chemical barrier becomes ineffective (layer too thin)
- H2S equilibrium profile extends to sediment surface
- Flux to water column increases dramatically
- In water column: if any O2 present → partial oxidation → S0 (blue tide)

**Key mechanism**: The `last_layer` parameter creates an implicit oxidation barrier by forcing the equilibrium profile to zero at the specified layer boundary. This is more robust than relying on explicit oxidation reactions, which may not keep pace with diffusive flux.

#### Limitation of 3-Layer Approach

The discrete 3-layer model cannot represent:
1. **Gradual redox transition**: BROM allows reactions at any depth; ERSEM confines reactions to specific layers
2. **Sub-layer O2 gradients**: Within Layer 1, ERSEM assumes uniform conditions
3. **Dynamic micro-zonation**: Rapid O2 fluctuations create transient micro-zones not captured

**Mitigation in Phase 1**: Use the existing `benthic_column_dissolved_matter` transport mechanism which calculates equilibrium profiles and diffusive fluxes correctly even when layers collapse. The Michaelis-Menten O2 limitation naturally reduces oxidation rates as O2 approaches zero.

#### Critical Implementation Note

The benthic sulfur cycle should NOT use Layer 2 O2/NO3 for cascade control when checking sulfate reduction conditions. Layer 3 is by definition anoxic. Instead:
- Sulfate reduction in Layer 3: Always active (Layer 3 = anoxic by definition)
- Oxidation in Layer 1: Rate-limited by actual O2 concentration in Layer 1

```fortran
! Correct: Use Layer 1 O2 for oxidation rate (naturally zero when D1m → minD)
O2_conc_1 = max(0.0_rk, O2_1) / max(D1m, 0.0001_rk)
R_H2S_ox_1 = K_H2S_ox * H2S_1 * O2_conc_1 / (O2_conc_1 + K_O2_half)
! When O2_1 → 0 and D1m → minD, R_H2S_ox_1 → 0 automatically
```

## New Module Structure

### Files to Create

```
src/
├── sulfur_cycle.F90              # Pelagic sulfur cycle (NEW)
├── benthic_sulfur_cycle.F90      # Benthic sulfur cycle (NEW)
```

### Files to Modify

**Phase 1 (minimal changes to existing code):**
```
src/
├── CMakeLists.txt                       # Add new source files to build
├── benthic_column_dissolved_matter.F90  # Add s_so4, s_h2s, s_s0 composition types
├── ersem_model_library.F90              # Register new modules
```

**Phase 2+ (K6 deprecation - NOT in Phase 1):**
```
src/
├── benthic_nitrogen_cycle.F90           # Remove K6 dependency (Phase 2)
├── benthic_bacteria.F90                 # Remove K6 reference from H2 (Phase 2)
├── oxygen.F90                           # Enforce non-negative oxygen (Phase 2)
```

## Module Implementation

### sulfur_cycle.F90 (Pelagic)

```fortran
#include "fabm_driver.h"

module ersem_sulfur_cycle

   use fabm_types
   use ersem_shared

   implicit none
   private

   type, extends(type_base_model), public :: type_ersem_sulfur_cycle
      ! State variable IDs
      type(type_state_variable_id) :: id_SO4, id_H2S, id_S0
      type(type_state_variable_id) :: id_O2

      ! Diagnostic variable IDs
      type(type_diagnostic_variable_id) :: id_R_H2S_ox, id_R_S0_ox

      ! Parameters
      real(rk) :: K_H2S_ox      ! H2S oxidation rate constant (1/d)
      real(rk) :: K_S0_ox       ! S0 oxidation rate constant (1/d)
      real(rk) :: K_O2_half     ! Half-saturation O2 for oxidation (mmol/m3)

   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class(type_ersem_sulfur_cycle), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! Set time unit to d-1 (ERSEM convention: rates are per day)
      self%dt = 86400._rk

      ! Get parameters
      call self%get_parameter(self%K_H2S_ox, 'K_H2S_ox', '1/d', &
           'H2S oxidation rate constant', default=0.5_rk)
      call self%get_parameter(self%K_S0_ox, 'K_S0_ox', '1/d', &
           'S0 oxidation rate constant', default=0.02_rk)
      call self%get_parameter(self%K_O2_half, 'K_O2_half', 'mmol/m^3', &
           'half-saturation O2 for oxidation', default=0.002_rk)

      ! Register state variables
      call self%register_state_variable(self%id_SO4, 'SO4', 'mmol S/m^3', &
           'sulfate', minimum=0.0_rk)
      call self%register_state_variable(self%id_H2S, 'H2S', 'mmol S/m^3', &
           'hydrogen sulfide', minimum=0.0_rk)
      call self%register_state_variable(self%id_S0, 'S0', 'mmol S/m^3', &
           'elemental sulfur', minimum=0.0_rk)

      ! Register dependency on oxygen (link to external O2 state variable)
      call self%register_state_dependency(self%id_O2, 'O2', 'mmol O_2/m^3', 'oxygen')

      ! Register diagnostic variables for output
      call self%register_diagnostic_variable(self%id_R_H2S_ox, 'R_H2S_ox', &
           'mmol S/m^3/d', 'H2S oxidation rate', source=source_do)
      call self%register_diagnostic_variable(self%id_R_S0_ox, 'R_S0_ox', &
           'mmol S/m^3/d', 'S0 oxidation rate', source=source_do)

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_sulfur_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: SO4, H2S, S0, O2
      real(rk) :: R_H2S_ox, R_S0_ox, f_O2

      _LOOP_BEGIN_
         _GET_(self%id_SO4, SO4)
         _GET_(self%id_H2S, H2S)
         _GET_(self%id_S0, S0)
         _GET_(self%id_O2, O2)

         ! Ensure non-negative oxygen for rate calculation
         O2 = max(0.0_rk, O2)

         ! Oxygen limitation function (Michaelis-Menten)
         f_O2 = O2 / (O2 + self%K_O2_half)

         ! R2: H2S + 0.5*O2 -> S0 + H2O
         R_H2S_ox = self%K_H2S_ox * H2S * f_O2

         ! R3: S0 + 1.5*O2 + H2O -> SO4 + 2H+
         R_S0_ox = self%K_S0_ox * S0 * f_O2

         ! Set ODEs (rates in per day, matching self%dt)
         _SET_ODE_(self%id_H2S, -R_H2S_ox)
         _SET_ODE_(self%id_S0,   R_H2S_ox - R_S0_ox)
         _SET_ODE_(self%id_SO4,  R_S0_ox)
         _SET_ODE_(self%id_O2,  -0.5_rk * R_H2S_ox - 1.5_rk * R_S0_ox)

         ! Set diagnostics
         _SET_DIAGNOSTIC_(self%id_R_H2S_ox, R_H2S_ox)
         _SET_DIAGNOSTIC_(self%id_R_S0_ox, R_S0_ox)

      _LOOP_END_
   end subroutine do

end module ersem_sulfur_cycle
```

### benthic_sulfur_cycle.F90 (Benthic)

```fortran
#include "fabm_driver.h"

module ersem_benthic_sulfur_cycle

   use fabm_types
   use ersem_shared

   implicit none
   private

   type, extends(type_base_model), public :: type_ersem_benthic_sulfur_cycle
      ! Layer-specific state variable dependencies
      ! Phase 1: Only Layer 1 (oxidation) and Layer 3 (sulfate reduction) are active
      type(type_bottom_state_variable_id) :: id_SO4_1, id_SO4_3
      type(type_bottom_state_variable_id) :: id_H2S_1, id_H2S_3
      type(type_bottom_state_variable_id) :: id_S0_1
      type(type_bottom_state_variable_id) :: id_O2_1
      type(type_horizontal_dependency_id) :: id_D1m, id_D2m, id_Dtot

      ! Link to organic matter decomposition rate from H2 bacteria
      type(type_horizontal_dependency_id) :: id_remin_rate

      ! Diagnostic variables
      type(type_horizontal_diagnostic_variable_id) :: id_R_sulfate_red
      type(type_horizontal_diagnostic_variable_id) :: id_R_H2S_ox_ben
      type(type_horizontal_diagnostic_variable_id) :: id_R_S0_ox_ben

      ! Parameters
      real(rk) :: K_SO4_rd       ! Sulfate reduction rate constant (1/d)
      real(rk) :: K_SO4_half     ! Half-saturation for sulfate (mmol/m3)
      real(rk) :: K_H2S_ox       ! H2S oxidation rate constant (1/d)
      real(rk) :: K_S0_ox        ! S0 oxidation rate constant (1/d)
      real(rk) :: K_O2_half      ! Half-saturation for oxygen (mmol/m3)
      real(rk) :: stoich_S_C     ! Stoichiometry: mol S per mol C oxidized

   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self, configunit)
      class(type_ersem_benthic_sulfur_cycle), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! Set time unit to d-1 (ERSEM convention: rates are per day)
      self%dt = 86400._rk

      ! Get parameters (values from BROM)
      call self%get_parameter(self%K_SO4_rd, 'K_SO4_rd', '1/d', &
           'sulfate reduction rate constant', default=0.000005_rk)
      call self%get_parameter(self%K_SO4_half, 'K_SO4_half', 'mmol/m^3', &
           'half-saturation sulfate', default=1.6_rk)
      call self%get_parameter(self%K_H2S_ox, 'K_H2S_ox', '1/d', &
           'H2S oxidation rate constant', default=0.5_rk)
      call self%get_parameter(self%K_S0_ox, 'K_S0_ox', '1/d', &
           'S0 oxidation rate constant', default=0.02_rk)
      call self%get_parameter(self%K_O2_half, 'K_O2_half', 'mmol/m^3', &
           'half-saturation O2', default=0.002_rk)
      call self%get_parameter(self%stoich_S_C, 'stoich_S_C', 'mol S/mol C', &
           'stoichiometry of sulfate reduction', default=0.5_rk)

      ! Register dependencies for layer-specific sulfur variables
      ! Phase 1: Only Layer 1 (oxidation) and Layer 3 (reduction) active
      ! Layer 2 handles transport only via benthic_column_dissolved_matter
      call self%register_state_dependency(self%id_SO4_1, 'SO4_1', 'mmol S/m^2', &
           'sulfate in layer 1')
      call self%register_state_dependency(self%id_SO4_3, 'SO4_3', 'mmol S/m^2', &
           'sulfate in layer 3')

      call self%register_state_dependency(self%id_H2S_1, 'H2S_1', 'mmol S/m^2', &
           'hydrogen sulfide in layer 1')
      call self%register_state_dependency(self%id_H2S_3, 'H2S_3', 'mmol S/m^2', &
           'hydrogen sulfide in layer 3')

      call self%register_state_dependency(self%id_S0_1, 'S0_1', 'mmol S/m^2', &
           'elemental sulfur in layer 1')

      ! Oxygen dependency for Layer 1 oxidation
      call self%register_state_dependency(self%id_O2_1, 'O2_1', 'mmol O_2/m^2', &
           'oxygen in layer 1')

      ! Layer depths (using standard variables from ersem_shared)
      call self%register_dependency(self%id_D1m, depth_of_bottom_interface_of_layer_1)
      call self%register_dependency(self%id_D2m, depth_of_bottom_interface_of_layer_2)
      call self%register_dependency(self%id_Dtot, depth_of_sediment_column)

      ! Organic matter remineralization rate (link to H2 bacteria)
      call self%register_dependency(self%id_remin_rate, 'remin_rate', 'mmol C/m^2/d', &
           'organic matter remineralization rate in layer 3')

      ! Diagnostic variables (domain=domain_bottom required for benthic diagnostics)
      call self%register_diagnostic_variable(self%id_R_sulfate_red, 'R_sulfate_red', &
           'mmol S/m^2/d', 'sulfate reduction rate', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_H2S_ox_ben, 'R_H2S_ox_ben', &
           'mmol S/m^2/d', 'benthic H2S oxidation rate', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_S0_ox_ben, 'R_S0_ox_ben', &
           'mmol S/m^2/d', 'benthic S0 oxidation rate', &
           domain=domain_bottom, source=source_do_bottom)

   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_ersem_benthic_sulfur_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: SO4_3, H2S_1, H2S_3, S0_1
      real(rk) :: O2_1
      real(rk) :: D1m, D2m, Dtot
      real(rk) :: D1m_safe, O2_conc_1   ! Guarded values for safe division
      real(rk) :: remin_rate
      real(rk) :: R_sulfate_red, R_H2S_ox_1, R_S0_ox_1
      real(rk) :: SO4_conc_3, layer3_thickness

      _HORIZONTAL_LOOP_BEGIN_

         ! Get sulfur variables needed for reactions
         ! Layer 3: sulfate reduction
         _GET_HORIZONTAL_(self%id_SO4_3, SO4_3)
         _GET_HORIZONTAL_(self%id_H2S_3, H2S_3)
         ! Layer 1: oxidation reactions
         _GET_HORIZONTAL_(self%id_H2S_1, H2S_1)
         _GET_HORIZONTAL_(self%id_S0_1, S0_1)
         _GET_HORIZONTAL_(self%id_O2_1, O2_1)

         ! Get layer depths
         _GET_HORIZONTAL_(self%id_D1m, D1m)
         _GET_HORIZONTAL_(self%id_D2m, D2m)
         _GET_HORIZONTAL_(self%id_Dtot, Dtot)

         ! Get remineralization rate from H2 bacteria
         _GET_HORIZONTAL_(self%id_remin_rate, remin_rate)

         ! Calculate layer 3 thickness and SO4 concentration
         ! Guard against division by zero when layer collapses
         layer3_thickness = max(Dtot - D2m, 0.0001_rk)
         SO4_conc_3 = SO4_3 / layer3_thickness

         ! Guard D1m against extremely small values to prevent numerical issues
         ! minD = 0.0001 m is the minimum layer thickness in ERSEM
         D1m_safe = max(D1m, 0.0001_rk)

         ! ============================================================
         ! PHASE 1 SIMPLIFICATION: No electron acceptor cascade check
         ! ============================================================
         ! In ERSEM's 3-layer model, Layer 3 is BY DEFINITION anoxic
         ! (below the denitrification layer). The layer boundary D2m
         ! is calculated based on nitrate penetration depth.
         ! Therefore, sulfate reduction in Layer 3 is always enabled.
         !
         ! Note: In Phase 2, Layer 2 sulfate reduction may be added
         ! when NO3 is depleted, requiring cascade logic there.
         ! ============================================================

         ! R1: Sulfate reduction in Layer 3 (always active - Layer 3 is anoxic)
         ! Rate proportional to organic matter remineralization by H2 bacteria
         R_sulfate_red = self%K_SO4_rd * SO4_conc_3 / (SO4_conc_3 + self%K_SO4_half) &
                       * remin_rate * self%stoich_S_C

         ! R2: H2S oxidation in layer 1
         ! Convert depth-integrated O2 to concentration for rate calculation
         ! Use D1m_safe to prevent division by very small numbers
         O2_conc_1 = max(0.0_rk, O2_1) / D1m_safe
         R_H2S_ox_1 = self%K_H2S_ox * H2S_1 * O2_conc_1 / (O2_conc_1 + self%K_O2_half)

         ! R3: S0 oxidation in layer 1
         R_S0_ox_1 = self%K_S0_ox * S0_1 * O2_conc_1 / (O2_conc_1 + self%K_O2_half)

         ! Set ODEs for layer 3 (sulfate reduction)
         _SET_BOTTOM_ODE_(self%id_SO4_3, -R_sulfate_red)
         _SET_BOTTOM_ODE_(self%id_H2S_3,  R_sulfate_red)

         ! Set ODEs for layer 1 (oxidation reactions)
         _SET_BOTTOM_ODE_(self%id_H2S_1, -R_H2S_ox_1)
         _SET_BOTTOM_ODE_(self%id_S0_1,   R_H2S_ox_1 - R_S0_ox_1)
         _SET_BOTTOM_ODE_(self%id_SO4_1,  R_S0_ox_1)
         _SET_BOTTOM_ODE_(self%id_O2_1,  -0.5_rk * R_H2S_ox_1 - 1.5_rk * R_S0_ox_1)

         ! Set diagnostics
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_sulfate_red, R_sulfate_red)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_H2S_ox_ben, R_H2S_ox_1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_S0_ox_ben, R_S0_ox_1)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module ersem_benthic_sulfur_cycle
```

### Modifications to benthic_column_dissolved_matter.F90

Add new composition types in the `initialize` subroutine:

```fortran
! In the select case block for composition:
case ('s_so4')
   call initialize_constituent(self, self%constituents(iconstituent), &
        profile, profile%constituents(iconstituent), &
        's_so4', 'mmol S/m^2', 'sulfate')
case ('s_h2s')
   call initialize_constituent(self, self%constituents(iconstituent), &
        profile, profile%constituents(iconstituent), &
        's_h2s', 'mmol S/m^2', 'hydrogen sulfide')
case ('s_s0')
   call initialize_constituent(self, self%constituents(iconstituent), &
        profile, profile%constituents(iconstituent), &
        's_s0', 'mmol S/m^2', 'elemental sulfur')
```

### Modifications to ersem_model_library.F90

Register the new modules. ERSEM uses a factory pattern with `select case` to map YAML model names to Fortran types:

```fortran
! Add use statements at the top of the module (after existing use statements)
use ersem_sulfur_cycle
use ersem_benthic_sulfur_cycle

! In the 'create' subroutine, add cases to the select case block:
! (add before the 'case default' line)
case ('sulfur_cycle');         allocate(type_ersem_sulfur_cycle::model)
case ('benthic_sulfur_cycle'); allocate(type_ersem_benthic_sulfur_cycle::model)
```

After these changes, the models can be referenced in YAML as `ersem/sulfur_cycle` and `ersem/benthic_sulfur_cycle`.

### Modifications to src/CMakeLists.txt

Add the new source files to the `add_library` command:

```cmake
add_library(fabm_models_ersem OBJECT
            # ... existing files ...
            sulfur_cycle.F90
            benthic_sulfur_cycle.F90
           )
```

### Modifications to oxygen.F90 (Phase 2+)

**Note**: This modification is NOT required for Phase 1. K6 continues to handle oxygen debt.

```fortran
! Phase 2+: Change the comment and behavior
! Note: negative oxygen concentrations are NO LONGER permitted.
! Oxygen debt is now represented explicitly by H2S.

! Ensure nonnegative constraint is applied by setting minimum=0.0_rk
! in register_state_variable call
```

### Modifications to benthic_nitrogen_cycle.F90 (Phase 2+)

**Note**: This modification is NOT required for Phase 1. K6 coexists with H2S.

```fortran
! Phase 2+: Remove K6 dependency
! - Remove id_K6_sms dependency
! - Remove K6_calculator child model
! - Remove K6-related calculations in do_bottom

! The electron acceptor cascade will then be handled by:
! 1. Oxygen consumption (primary) - existing
! 2. Denitrification (when O2 depleted) - existing
! 3. Sulfate reduction (when O2 and NO3 depleted) - new benthic_sulfur_cycle
```

## YAML Configuration Example

```yaml
# ==============================================================
# Sulfur Cycle Configuration (Phase 1: Minimal Implementation)
# ==============================================================

# --- Pelagic Sulfur Species ---
N8_sulfur:
  long_name: pelagic sulfur cycle
  model: ersem/sulfur_cycle
  parameters:
    K_H2S_ox: 0.5              # H2S oxidation rate constant (1/d)
    K_S0_ox: 0.02              # S0 oxidation rate constant (1/d)
    K_O2_half: 0.002           # Half-saturation O2 (mmol/m3)
  initialization:
    SO4: 28000.0               # Seawater sulfate ~28 mM
    H2S: 0.0                   # Initial H2S
    S0: 0.0                    # Initial elemental sulfur
  coupling:
    O2: O2/o                   # Link to oxygen

# --- Benthic Sulfate ---
G2s_SO4:
  long_name: benthic sulfate
  model: ersem/benthic_column_dissolved_matter
  parameters:
    composition: s_so4
  initialization:
    s_so4: 8400.0              # ~28mM * 0.3m depth
  coupling:
    s_so4_pel: N8_sulfur/SO4   # Link to pelagic sulfate

# --- Benthic Hydrogen Sulfide ---
# NOTE: Do NOT use last_layer: 1 here - it conflicts with G2 (oxygen) which controls D1m
# Instead, rely on increased oxidation rates in ben_sulfur and pel_sulfur modules
G2s_H2S:
  long_name: benthic hydrogen sulfide
  model: ersem/benthic_column_dissolved_matter
  parameters:
    composition: s_h2s
    # last_layer: 3 (default) - H2S exists in all layers, oxidation controls flux
  initialization:
    s_h2s: 0.0
  coupling:
    s_h2s_pel: N8_sulfur/H2S   # Link to pelagic H2S - CRITICAL!

# --- Benthic Elemental Sulfur ---
G2s_S0:
  long_name: benthic elemental sulfur
  model: ersem/benthic_column_dissolved_matter
  parameters:
    composition: s_s0
    last_layer: 3
  initialization:
    s_s0: 0.0
  coupling:
    s_s0_pel: N8_sulfur/S0     # Link to pelagic S0

# --- Benthic Sulfur Reactions ---
ben_sulfur:
  long_name: benthic sulfur cycle
  model: ersem/benthic_sulfur_cycle
  parameters:
    K_SO4_rd: 0.000005         # Sulfate reduction rate (1/d)
    K_SO4_half: 1.6            # Half-saturation sulfate (mmol/m3)
    K_H2S_ox: 0.5              # H2S oxidation rate (1/d)
    K_S0_ox: 0.02              # S0 oxidation rate (1/d)
    K_O2_half: 0.002           # Half-saturation O2 (mmol/m3)
    O2_threshold: 0.01         # O2 threshold for sulfate reduction (mmol/m3)
    NO3_threshold: 0.005       # NO3 threshold for sulfate reduction (mmol/m3)
    O2_scale: 0.01             # O2 transition width for tanh (mmol/m3)
    NO3_scale: 0.005           # NO3 transition width for tanh (mmol/m3)
    stoich_S_C: 0.5            # mol S per mol C
  coupling:
    # Layer-specific sulfur species
    SO4_1: G2s_SO4/per_layer/s_so4_1
    SO4_2: G2s_SO4/per_layer/s_so4_2
    SO4_3: G2s_SO4/per_layer/s_so4_3
    H2S_1: G2s_H2S/per_layer/s_h2s_1
    H2S_2: G2s_H2S/per_layer/s_h2s_2
    H2S_3: G2s_H2S/per_layer/s_h2s_3
    S0_1: G2s_S0/per_layer/s_s0_1
    S0_2: G2s_S0/per_layer/s_s0_2
    S0_3: G2s_S0/per_layer/s_s0_3
    # Oxygen and nitrate
    O2_1: G2/per_layer/o1
    O2_2: G2/per_layer/o2
    NO3_2: K3/per_layer/n2
    # Layer depths
    D1m: ben_col/D1m
    D2m: ben_col/D2m
    Dtot: ben_col/Dtot
    # Organic matter remineralization (link to H2 bacteria)
    remin_rate: H2/respiration_rate
```

## Testing Strategy

### Unit Tests

Create `github-actions/pyfabm-ersem/test_sulfur_cycle.py`:

```python
import pytest
import pyfabm
import numpy as np

class TestSulfurCycle:
    """Test suite for sulfur cycle implementation"""

    # Test conditions
    O2_OXIC = 0.2          # mmol/m3 - oxic conditions
    O2_ANOXIC = 0.001      # mmol/m3 - anoxic conditions
    H2S_INIT = 0.1         # mmol/m3 - initial H2S
    S0_INIT = 0.05         # mmol/m3 - initial S0
    SO4_INIT = 28000.0     # mmol/m3 - seawater sulfate

    # Tolerances
    RATE_TOLERANCE = 1e-10  # For "zero" rate checks
    MASS_TOLERANCE = 1e-6   # For mass conservation (relative)

    def test_sulfate_reduction_requires_anoxia(self, model):
        """Sulfate reduction should only occur under anoxic conditions"""
        # Test Case 1: Oxic conditions (O2 = 0.2 mmol/m3)
        # Expected: f_anox ≈ 0, R_sulfate_red ≈ 0
        model.state['O2'] = self.O2_OXIC
        model.state['NO3'] = 0.001  # Below threshold
        rate_oxic = model.get_diagnostic('R_sulfate_red')
        assert rate_oxic < self.RATE_TOLERANCE, \
            f"Sulfate reduction should be ~0 under oxic conditions, got {rate_oxic}"

        # Test Case 2: Anoxic conditions (O2 = 0.001 mmol/m3)
        # Expected: f_anox ≈ 1, R_sulfate_red > 0
        model.state['O2'] = self.O2_ANOXIC
        rate_anoxic = model.get_diagnostic('R_sulfate_red')
        assert rate_anoxic > self.RATE_TOLERANCE, \
            f"Sulfate reduction should be > 0 under anoxic conditions, got {rate_anoxic}"

    def test_h2s_oxidation_requires_oxygen(self, model):
        """H2S oxidation should require oxygen"""
        model.state['H2S'] = self.H2S_INIT

        # Test Case 1: With oxygen
        model.state['O2'] = self.O2_OXIC
        rate_with_O2 = model.get_diagnostic('R_H2S_ox')
        assert rate_with_O2 > 0, "H2S oxidation rate should be > 0 with oxygen"

        # Test Case 2: Without oxygen
        model.state['O2'] = 0.0
        rate_no_O2 = model.get_diagnostic('R_H2S_ox')
        assert rate_no_O2 < self.RATE_TOLERANCE, \
            f"H2S oxidation rate should be ~0 without oxygen, got {rate_no_O2}"

    def test_sulfur_mass_conservation(self, model):
        """Total sulfur should be conserved in closed system"""
        # Initial total sulfur
        S_total_init = (model.state['SO4'] + model.state['H2S'] + model.state['S0'])

        # Run for 1 day
        model.integrate(dt=86400.0)

        # Final total sulfur
        S_total_final = (model.state['SO4'] + model.state['H2S'] + model.state['S0'])

        # Check conservation (relative error)
        rel_error = abs(S_total_final - S_total_init) / S_total_init
        assert rel_error < self.MASS_TOLERANCE, \
            f"Sulfur not conserved: initial={S_total_init}, final={S_total_final}, error={rel_error}"

    def test_oxygen_consumption_stoichiometry(self, model):
        """Verify O2 consumption matches stoichiometry"""
        model.state['H2S'] = self.H2S_INIT
        model.state['S0'] = self.S0_INIT
        model.state['O2'] = self.O2_OXIC

        R_H2S_ox = model.get_diagnostic('R_H2S_ox')
        R_S0_ox = model.get_diagnostic('R_S0_ox')
        dO2_dt = model.get_rate('O2')

        # Expected O2 consumption: -0.5*R_H2S_ox - 1.5*R_S0_ox
        expected_dO2 = -0.5 * R_H2S_ox - 1.5 * R_S0_ox

        assert abs(dO2_dt - expected_dO2) < self.RATE_TOLERANCE, \
            f"O2 stoichiometry mismatch: got {dO2_dt}, expected {expected_dO2}"

    def test_blue_tide_scenario(self, model):
        """Test blue tide conditions: anoxic bottom water -> S0 accumulation"""
        # Initial conditions: anoxic water, H2S flux from sediment
        model.state['O2'] = self.O2_ANOXIC
        model.state['H2S'] = 0.5  # mmol/m3 - elevated H2S from sediment
        model.state['S0'] = 0.0

        # Run for 1 hour
        model.integrate(dt=3600.0)

        # Under anoxic conditions:
        # - H2S oxidation should be minimal
        # - S0 should remain low (no significant oxidation of H2S)
        # But if slight O2 present, some S0 forms

        # Now simulate water column mixing bringing some O2
        model.state['O2'] = 0.05  # mmol/m3 - slight oxygenation
        model.integrate(dt=3600.0)

        # S0 should increase as H2S is partially oxidized
        assert model.state['S0'] > 0, "S0 should increase during partial H2S oxidation"
```

### Specific Test Conditions

| Test | Initial Conditions | Expected Result | Tolerance |
|------|-------------------|-----------------|-----------|
| Anoxic suppression | O2=0.2, NO3=0.001 mmol/m³ | R_sulfate_red < 1e-10 | 1e-10 mmol/m²/d |
| Sulfate reduction | O2=0.001, NO3=0.001 mmol/m³, SO4=28 mM | R_sulfate_red > 0 | - |
| H2S oxidation | H2S=0.1, O2=0.2 mmol/m³ | R_H2S_ox > 0 | - |
| Mass conservation | SO4+H2S+S0 initial | SO4+H2S+S0 after 1d same | <0.0001% |
| O2 stoichiometry | H2S=0.1, S0=0.05, O2=0.2 | dO2/dt = -0.5*R2 - 1.5*R3 | 1e-10 |

### Integration Tests

Create a test configuration in `testcases/fabm-ersem-sulfur-test.yaml` that:
1. Simulates a shallow coastal system (10m depth, 0.3m sediment)
2. Initial conditions: O2=200 mmol/m³, SO4=28000 mmol/m³, H2S=0, S0=0
3. Forces bottom water anoxia by reducing surface O2 flux
4. Verifies:
   - H2S accumulation in sediment Layer 3 after 10 days
   - H2S flux to water column when Layer 1 O2 depletes
   - S0 peak in water column during partial re-oxygenation

## Implementation Phases

### Phase 1 (Current Proposal)

**Scope**:
- New sulfur species: SO4, H2S, S0 in pelagic and benthic compartments
- New reactions: R1-R3 (sulfate reduction in Layer 3, H2S/S0 oxidation)
- K6/N6 coexistence: Existing oxygen debt mechanism remains unchanged
- Minimal changes to existing modules (only add composition types and register new modules)

**Files to create**: `sulfur_cycle.F90`, `benthic_sulfur_cycle.F90`
**Files to modify**: `CMakeLists.txt`, `benthic_column_dissolved_matter.F90`, `ersem_model_library.F90`

**Validation**: Compare model behavior with and without sulfur cycle enabled

### Phase 2 (Future)

**Scope**:
- Deprecate K6: Remove oxygen debt mechanism, rely on explicit H2S
- Enforce O2 ≥ 0 constraint
- Add Layer 2 sulfate reduction (when NO3 is depleted)
- Add thiosulfate (S2O3) as intermediate species (completing BROM pathway: H2S → S0 → S2O3 → SO4)
- ~~Add thiodenitrification reactions~~ (H2S-NO3 oxidation now implemented in Phase 1)
- Add additional thiodenitrification pathways from BROM:
  - ~~R4: 5 H2S + 2 NO3⁻ → 5 S0 + N2 + 4 H2O~~ (implemented in Phase 1 as H2S-NO3 oxidation)
  - R5: 5 S0 + 6 NO3⁻ → 5 SO4²⁻ + 3 N2 (K_s0_no3 = 0.9 1/d in BROM)
  - R6: 5 S2O3²⁻ + 8 NO3⁻ → 10 SO4²⁻ + 4 N2 (K_s2o3_no3 = 0.01 1/d in BROM)
- Add S0 disproportionation: 4 S0 + 3 H2O → 2 H2S + S2O3²⁻ + 2 H⁺ (K_s0_disp = 0.001 1/d)

**Files to modify**: `sulfur_cycle.F90`, `benthic_sulfur_cycle.F90`, `benthic_nitrogen_cycle.F90`, `benthic_bacteria.F90`, `oxygen.F90`

### Phase 3 (Future)

**Scope**:
- Remove K6 completely
- Add iron-sulfur species from BROM:
  - FeS (iron monosulfide): H2S + Fe²⁺ → FeS + 2H⁺
  - FeS2 (pyrite): FeS + S0 → FeS2
  - MnS (manganese sulfide): H2S + Mn²⁺ → MnS + 2H⁺
- Integrate with ERSEM iron module (-DERSEM_USE_IRON flag)
- Add temperature dependencies (Q10 formulations as in BROM)
- Full BROM-style redox chemistry including Mn cycling

**New parameters from BROM**:
- K_fes_form = 1.0e-6 (1/d) - FeS formation rate
- K_fes_ox = 1.0e-5 (1/d) - FeS oxidation rate
- K_fes2_form = 1.0e-7 (1/d) - Pyrite formation rate
- K_fes2_ox = 1.0e-6 (1/d) - Pyrite oxidation rate

## References

1. BROM model source: https://github.com/fabm-model/fabm/tree/master/src/models/niva/brom
2. BROM test configuration: https://github.com/fabm-model/fabm/blob/master/testcases/fabm-niva-brom.yaml
3. Butenschön et al. (2016). ERSEM 15.06. Geosci. Model Dev., 9, 1293-1339
4. Jørgensen, B. B. (1982). Mineralization of organic matter in the sea bed. Nature, 296, 643-645
5. Canfield, D. E. et al. (1993). Pathways of organic carbon oxidation in three continental margin sediments. Mar. Geol., 113, 27-40
6. Yakushev, E. V. et al. (2017). Bottom RedOx Model (BROM v.1.1): a coupled benthic–pelagic model for simulation of water and sediment biogeochemistry. Geosci. Model Dev., 10, 453-482

## Notes for Implementers

1. **Critical: Do NOT use `last_layer: 1` for H2S**:
   - Using `last_layer: 1` for K_H2S conflicts with G2 (oxygen) which also uses `last_layer: 1`
   - Both modules would try to control D1m, causing it to grow to unrealistic values (~150 mm)
   - The `last_layer` mechanism is designed for ONE constituent per layer boundary
   - **Correct approach**: Use the **oxic barrier mechanism** (Approach 3):
     - Set `K_barrier: 1000.0` (1/m) - barrier effectiveness
     - Set `K_barrier_rate: 100.0` (1/d) - oxidation rate in bottom water
   - Combined with increased oxidation rates:
     - Benthic: K_H2S_ox = 5.0 (10x default)
     - Pelagic: K_H2S_ox = 10.0 (20x default)
   - This achieves H2S ≈ 0 under oxic conditions (D1m > ~1mm and O2 present)

2. **Oxic barrier coupling requirements**:
   - The `ben_sulfur` module requires coupling to pelagic variables:
     - `H2S_pel`: pelagic H2S for barrier oxidation
     - `S0_pel`: pelagic S0 (product of barrier oxidation)
     - `O2_pel`: pelagic O2 (consumed by barrier oxidation)
   - Add these couplings to the YAML configuration (see example below)

3. **Phase 1: K6 coexistence**:
   - Keep K6 and N6 in parallel with new H2S variables
   - Do NOT modify benthic_nitrogen_cycle.F90 or benthic_bacteria.F90 in Phase 1
   - This allows validation of new sulfur cycle against existing behavior
   - K6 removal is deferred to Phase 2/3

2. **Test incrementally**:
   - First implement pelagic sulfur_cycle.F90 alone
   - Verify with simple box model (0D)
   - Then add benthic_sulfur_cycle.F90
   - Test benthic-pelagic coupling
   - Finally, in Phase 2+, modify existing modules to remove K6

3. **Watch for unit conversions**:
   - Pelagic: mmol/m³ (concentration)
   - Benthic: mmol/m² (depth-integrated)
   - Layer concentrations need division by layer thickness
   - All parameters in this document use **mmol/m³** (not µmol/m³)

4. **Guard against numerical issues**:
   - Use `max(D1m, 0.0001)` before dividing by layer thickness
   - Layer thickness can become very small under anoxic conditions
   - The 0.0001 m minimum matches ERSEM's `minD` parameter

5. **Verify mass conservation**: Total sulfur (SO4 + H2S + S0) should be conserved in closed systems. Run conservation tests before integration.

6. **Oxygen budget**: The sulfur modules contribute to O2 rates via `_SET_ODE_`. FABM accumulates all contributions. Verify that total O2 consumption is physically reasonable.

7. **Consult BROM source**: https://github.com/fabm-model/fabm/blob/master/src/models/niva/brom/brom_redox.F90 for reference implementation details. Note that BROM is a water column model; this guide adapts its chemistry to ERSEM's 3-layer benthic structure.
