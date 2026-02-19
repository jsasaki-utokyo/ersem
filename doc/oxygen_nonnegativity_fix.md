# Fix: Oxygen Non-Negativity Constraint and Benthic Exchange Instability

## Summary

A critical numerical instability in ERSEM allowed dissolved oxygen (O2) to reach
physically impossible values — both extremely negative (-2158 mmol/m3 = -69 mg/L)
and extremely positive (+4542 mmol/m3 = +145 mg/L) — in the bottom layer of 1D
GOTM-ERSEM simulations. The root cause was a missing non-negativity constraint on
the O2 state variable, combined with an unstable feedback loop in the benthic
column dissolved matter module.

## Affected Files

| File | Change |
|------|--------|
| `src/oxygen.F90` | Added `minimum=0` to O2 state variable registration |
| `src/benthic_column_dissolved_matter.F90` | Changed O2 constituent from `nonnegative=legacy_ersem_compatibility` (.false.) to `nonnegative=.true.` |

A corresponding change to the host model configuration is also recommended:

| File | Change |
|------|--------|
| `gotm.yaml` | Changed `ode_method: fe` (Forward Euler) to `ode_method: emp2` (Extended Modified Patankar 2nd order) |

## Bug Description

### Symptom

In a 1D GOTM-ERSEM simulation of Tokyo Bay (Kawasaki station, 20m depth, 30
vertical layers), dissolved oxygen in the bottom layer exhibited catastrophic
oscillations during summer stratification:

- O2 dropped to -2158 mmol/m3 (physically: this would mean -69 mg/L)
- O2 spiked to +4542 mmol/m3 (physically: this would mean +145 mg/L, ~15x saturation)
- Oscillations persisted for weeks with ~4-day period

### Root Cause Analysis

The instability arose from three interacting factors:

#### 1. Missing `minimum=0` on O2 state variable (`oxygen.F90`)

The O2 state variable was registered without a minimum bound:

```fortran
! BEFORE (buggy)
call self%register_state_variable(self%id_O2o,'o','mmol O_2/m^3','oxygen',300._rk)
```

This had two consequences in GOTM:
- **`repair_state`** (which clips variables to [minimum, maximum]) had no minimum to enforce
- **`posconc`** flag was set to 0 (only set to 1 when `minimum >= 0`), so the vertical
  diffusion solver did not apply positivity constraints

As a result, O2 could freely go negative via Forward Euler integration.

#### 2. Benthic O2 with `nonnegative=.false.` (`benthic_column_dissolved_matter.F90`)

In the benthic column dissolved matter module, the O2 constituent was registered with:

```fortran
! BEFORE (buggy)
call initialize_constituent(..., 'o', ..., nonnegative=legacy_ersem_compatibility)
```

Since `legacy_ersem_compatibility = .false.` (in `shared.F90`), the O2 constituent
had `nonnegative = .false.`. This disabled the Patankar trick for benthic O2 and
enabled a secondary state variable `c_int_deep` (oxygen mass below the zero
isocline in the sediment profile).

#### 3. Equilibrium Profile Amplification Mechanism

The benthic column dissolved matter module computes an equilibrium O2 profile in
the sediment using a parabolic (quadratic) function. The curvature of this profile
is proportional to `P / sigma`, where `P` is the production/consumption rate and
`sigma` is the sediment diffusivity.

With typical ERSEM parameters:
- Sediment diffusivities: EDZ = 5 x 10^-5 m2/d (extremely small)
- This creates enormous curvature (~50,000 m^-2) in the equilibrium profile

The depth-integrated equilibrium O2 (`c_int_eq`) becomes extremely negative
(~-1000 to -1300 mmol/m2) when pelagic O2 goes even slightly negative.

The benthic-pelagic exchange flux is:

```
exchange = sms + (c_int + c_int_deep - c_int_eq) / relax
```

When pelagic O2 first goes slightly negative (~-20 mmol/m3):

1. `c_int_eq` jumps to ~-1000 mmol/m2 (due to the parabolic profile with small diffusivity)
2. `c_int_deep` is only ~-26 mmol/m2 (hasn't caught up yet)
3. The mismatch `(c_int + c_int_deep - c_int_eq)` is ~+974 mmol/m2
4. Divided by `relax` (5 days), this produces ~+195 mmol/m2/d of O2 flux into the pelagic
5. Divided by the bottom layer height (~0.667 m), this is ~+292 mmol/m3/d — a massive
   source of O2 in the bottom cell

This positive flux **pumps O2 back into the pelagic** at an enormous rate. Over
the next ~4 days, `c_int_deep` relaxes toward `c_int_eq`, eventually overshooting.
When the exchange reverses sign, O2 is driven negative again, and the cycle repeats.

The diagnostic evidence was definitive:
- `G2_o_pb_flux` reached +10,577 mmol/m2/d (= +15,865 mmol/m3/d in pelagic)
- `G2_o_deep` decreased from -26 to -16,620 mmol/m2
- `D1m` (oxygenated layer depth) was stuck at `minD` = 0.0001 m throughout

### Timeline of a Typical Oscillation Cycle

```
Jul 15:  O2 = -20 mmol/m3 (slightly negative from respiration)
         c_int_eq jumps to ~-1000 mmol/m2
         G2_o_deep = -26 mmol/m2 (lags behind)
         → exchange = +195 mmol/m2/d (positive = O2 injected into pelagic)

Jul 16:  O2 spikes to +1000-4500 mmol/m3 (massive injection)
         G2_o_deep rapidly decreasing toward c_int_eq

Jul 17-19: O2 remains elevated (1000-4000 mmol/m3)
           G2_o_deep still relaxing, exchange flux decreasing

Jul 20:  G2_o_deep overshoots c_int_eq
         Exchange reverses → O2 pulled back to benthic
         O2 crashes to -2000 mmol/m3
         Cycle may repeat
```

## Fix

### Fix 1: Add `minimum=0` to O2 (`oxygen.F90`, line 51)

```fortran
! AFTER (fixed)
call self%register_state_variable(self%id_O2o,'o','mmol O_2/m^3','oxygen',300._rk,minimum=0._rk)
```

This enables:
- `repair_state` clips O2 to [0, inf) after each timestep
- `posconc = 1` enables positivity constraints in the vertical diffusion solver
- Prevents O2 from going negative in the first place

**Rationale**: When sulfate reduction is modeled explicitly (e.g., via `pel_sulfur`),
oxygen debt is properly tracked through H2S concentration rather than negative O2.
The original ERSEM design of allowing negative O2 as an "oxygen debt proxy" is no
longer needed.

### Fix 2: Set `nonnegative=.true.` for benthic O2 (`benthic_column_dissolved_matter.F90`, line 127)

```fortran
! AFTER (fixed)
call initialize_constituent(..., 'o', ..., nonnegative=.true.)
```

This:
- Enables the Patankar trick for benthic O2 exchange calculations
- Disables `c_int_deep` (no longer tracks negative mass below zero isocline)
- Prevents the equilibrium profile amplification mechanism entirely

### Fix 3: Switch ODE solver from Forward Euler to EMP2 (`gotm.yaml`)

```yaml
# AFTER (fixed)
ode_method: emp2  # Extended Modified Patankar, 2nd order
```

The EMP2 solver:
- Guarantees non-negativity of all state variables by construction
- Is 2nd-order accurate (vs 1st-order for Forward Euler)
- Prevents numerical undershoot/overshoot that triggers the amplification cascade

## Impact

This fix is relevant for any ERSEM configuration where:
1. Bottom O2 can approach zero (e.g., eutrophic environments, deep stratified waters)
2. The benthic column dissolved matter module is used with O2
3. `legacy_ersem_compatibility = .false.` (the current default)

Without this fix, any configuration that develops bottom-layer hypoxia is at risk
of triggering the benthic exchange instability, leading to unphysical O2 values
and potential simulation failure.
