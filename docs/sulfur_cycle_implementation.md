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

1. **Replace K6 (oxygen debt)** with explicit H2S
2. **Sulfur cycle only** - no iron-sulfur interactions (FeS, FeS2) in Phase 1
3. **Explicit sulfur species** in both water column and sediments: SO4^2-, H2S, S0
4. **Minimal implementation first** (Phase 1), with expansion possible in later phases

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

### Variables to Remove

| Variable | Reason |
|----------|--------|
| `K6` | Replaced by explicit H2S |
| `N6` | Replaced by N8s_H2S |
| `o_deep` (negative values) | Oxygen should be non-negative; use H2S instead |

## Chemical Reactions

### Reaction Equations

```
R1: Sulfate Reduction (Layer 2-3, anaerobic)
    SO4^2- + 2CH2O -> H2S + 2HCO3^-

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
f_anox(O2) = 0.5 * (1.0 - tanh((O2 - O2_threshold) / O2_scale))

! Nitrate inhibition (sulfate reduction inhibited by nitrate)
f_no_NO3(NO3) = 0.5 * (1.0 - tanh((NO3 - NO3_threshold) / NO3_scale))

! Combined sulfate reduction factor
f_sulfate_reduction = f_anox * f_no_NO3
```

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

| Parameter | Value | Units | Description |
|-----------|-------|-------|-------------|
| `K_SO4_rd` | 0.000005 | 1/d | Sulfate reduction rate constant |
| `K_SO4_half` | 1600 | umol/m^3 | Half-saturation for sulfate |
| `K_H2S_ox` | 0.5 | 1/d | H2S oxidation rate constant |
| `K_S0_ox` | 0.02 | 1/d | S0 oxidation rate constant |
| `K_O2_half` | 2 | umol/m^3 | Half-saturation for oxygen |
| `O2_threshold` | 10 | umol/m^3 | O2 threshold for sulfate reduction |
| `NO3_threshold` | 5 | umol/m^3 | NO3 threshold for sulfate reduction |

### Stoichiometry

| Reaction | O2 Consumption |
|----------|----------------|
| R1 | 0 (organic matter substitutes) |
| R2 | 0.5 mol O2 / mol H2S |
| R3 | 1.5 mol O2 / mol S0 |

## Integration with 3-Layer Benthic Model

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
  Reactions: Limited sulfate reduction (NO3 takes priority)

Layer 3 (Anoxic): D2m to d_tot
  Variables: G2s_SO4_3, G2s_H2S_3, G2s_S0_3
  Reactions: R1 (sulfate reduction) - PRIMARY H2S SOURCE
  Note: Linked to H2 (anaerobic bacteria) respiration
```

### Transport Mechanism

Use the existing `benthic_column_dissolved_matter` transport scheme:

- Equilibrium concentration profile calculation through all layers
- Diffusive transport based on layer-specific diffusivities (EDZ_1, EDZ_2, EDZ_3)
- Pelagic-benthic exchange at sediment-water interface

**Important**: H2S produced in Layer 3 diffuses through Layer 2 and Layer 1 before reaching the water column. It does NOT bypass intermediate layers. The existing transport code handles this correctly.

### Behavior Under Bottom Water Anoxia

When bottom water becomes anoxic:

```
Normal Conditions:          Anoxic Conditions:

Water: O2 > 0              Water: O2 ~ 0
  | H2S oxidized             | H2S NOT oxidized -> accumulates
-----------                -----------
Layer 1: D1m ~ 9mm         Layer 1: D1m -> minD (nearly zero)
  | H2S mostly oxidized      | No oxidation barrier
-----------                -----------
Layer 2: ~40mm             Layer 2: very thin
-----------                -----------
Layer 3: ~250mm            Layer 3: nearly entire column
  H2S production             Massive H2S production -> direct flux to water
```

## New Module Structure

### Files to Create

```
src/
├── sulfur_cycle.F90              # Pelagic sulfur cycle (NEW)
├── benthic_sulfur_cycle.F90      # Benthic sulfur cycle (NEW)
```

### Files to Modify

```
src/
├── benthic_column_dissolved_matter.F90  # Add s_so4, s_h2s, s_s0 composition
├── benthic_nitrogen_cycle.F90           # Remove K6 dependency
├── benthic_bacteria.F90                 # Remove K6 reference from H2
├── oxygen.F90                           # Enforce non-negative oxygen
├── ersem_model_library.F90              # Register new modules
```

## Module Implementation

### sulfur_cycle.F90 (Pelagic)

```fortran
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

      ! Set time unit to d-1
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

      ! Register dependency on oxygen
      call self%register_state_dependency(self%id_O2, 'O2', 'mmol O2/m^3', 'oxygen')

      ! Register diagnostic variables
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

         ! Oxygen limitation function
         f_O2 = O2 / (O2 + self%K_O2_half)

         ! R2: H2S + 0.5*O2 -> S0 + H2O
         R_H2S_ox = self%K_H2S_ox * H2S * f_O2

         ! R3: S0 + 1.5*O2 + H2O -> SO4 + 2H+
         R_S0_ox = self%K_S0_ox * S0 * f_O2

         ! Set ODEs
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
module ersem_benthic_sulfur_cycle

   use fabm_types
   use ersem_shared

   implicit none
   private

   type, extends(type_base_model), public :: type_ersem_benthic_sulfur_cycle
      ! Layer-specific state variable dependencies
      type(type_bottom_state_variable_id) :: id_SO4_1, id_SO4_2, id_SO4_3
      type(type_bottom_state_variable_id) :: id_H2S_1, id_H2S_2, id_H2S_3
      type(type_bottom_state_variable_id) :: id_S0_1, id_S0_2, id_S0_3
      type(type_bottom_state_variable_id) :: id_O2_1
      type(type_horizontal_dependency_id) :: id_O2_2, id_NO3_2
      type(type_horizontal_dependency_id) :: id_D1m, id_D2m, id_Dtot

      ! Link to organic matter decomposition rate
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
      real(rk) :: O2_threshold   ! O2 threshold for sulfate reduction (mmol/m3)
      real(rk) :: NO3_threshold  ! NO3 threshold for sulfate reduction (mmol/m3)
      real(rk) :: stoich_S_C     ! Stoichiometry: mol S per mol C oxidized

   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self, configunit)
      class(type_ersem_benthic_sulfur_cycle), intent(inout), target :: self
      integer, intent(in) :: configunit

      self%dt = 86400._rk

      ! Get parameters
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
      call self%get_parameter(self%O2_threshold, 'O2_threshold', 'mmol/m^3', &
           'O2 threshold for sulfate reduction', default=0.01_rk)
      call self%get_parameter(self%NO3_threshold, 'NO3_threshold', 'mmol/m^3', &
           'NO3 threshold for sulfate reduction', default=0.005_rk)
      call self%get_parameter(self%stoich_S_C, 'stoich_S_C', 'mol S/mol C', &
           'stoichiometry of sulfate reduction', default=0.5_rk)

      ! Register dependencies for layer-specific sulfur variables
      call self%register_state_dependency(self%id_SO4_1, 'SO4_1', 'mmol S/m^2', &
           'sulfate in layer 1')
      call self%register_state_dependency(self%id_SO4_2, 'SO4_2', 'mmol S/m^2', &
           'sulfate in layer 2')
      call self%register_state_dependency(self%id_SO4_3, 'SO4_3', 'mmol S/m^2', &
           'sulfate in layer 3')

      call self%register_state_dependency(self%id_H2S_1, 'H2S_1', 'mmol S/m^2', &
           'hydrogen sulfide in layer 1')
      call self%register_state_dependency(self%id_H2S_2, 'H2S_2', 'mmol S/m^2', &
           'hydrogen sulfide in layer 2')
      call self%register_state_dependency(self%id_H2S_3, 'H2S_3', 'mmol S/m^2', &
           'hydrogen sulfide in layer 3')

      call self%register_state_dependency(self%id_S0_1, 'S0_1', 'mmol S/m^2', &
           'elemental sulfur in layer 1')
      call self%register_state_dependency(self%id_S0_2, 'S0_2', 'mmol S/m^2', &
           'elemental sulfur in layer 2')
      call self%register_state_dependency(self%id_S0_3, 'S0_3', 'mmol S/m^2', &
           'elemental sulfur in layer 3')

      ! Oxygen and nitrate dependencies
      call self%register_state_dependency(self%id_O2_1, 'O2_1', 'mmol O2/m^2', &
           'oxygen in layer 1')
      call self%register_dependency(self%id_O2_2, 'O2_2', 'mmol O2/m^2', &
           'oxygen in layer 2')
      call self%register_dependency(self%id_NO3_2, 'NO3_2', 'mmol N/m^2', &
           'nitrate in layer 2')

      ! Layer depths
      call self%register_dependency(self%id_D1m, depth_of_bottom_interface_of_layer_1)
      call self%register_dependency(self%id_D2m, depth_of_bottom_interface_of_layer_2)
      call self%register_dependency(self%id_Dtot, depth_of_sediment_column)

      ! Organic matter remineralization rate (link to H2 bacteria)
      call self%register_dependency(self%id_remin_rate, 'remin_rate', 'mmol C/m^2/d', &
           'organic matter remineralization rate in layer 3')

      ! Diagnostic variables
      call self%register_diagnostic_variable(self%id_R_sulfate_red, 'R_sulfate_red', &
           'mmol S/m^2/d', 'sulfate reduction rate', source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_H2S_ox_ben, 'R_H2S_ox_ben', &
           'mmol S/m^2/d', 'benthic H2S oxidation rate', source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_S0_ox_ben, 'R_S0_ox_ben', &
           'mmol S/m^2/d', 'benthic S0 oxidation rate', source=source_do_bottom)

   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_ersem_benthic_sulfur_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: SO4_1, SO4_2, SO4_3
      real(rk) :: H2S_1, H2S_2, H2S_3
      real(rk) :: S0_1, S0_2, S0_3
      real(rk) :: O2_1, O2_2, NO3_2
      real(rk) :: D1m, D2m, Dtot
      real(rk) :: remin_rate
      real(rk) :: f_anox, f_no_NO3, f_sulfate_red
      real(rk) :: R_sulfate_red, R_H2S_ox_1, R_S0_ox_1
      real(rk) :: SO4_conc_3, layer3_thickness

      _HORIZONTAL_LOOP_BEGIN_

         ! Get sulfur variables
         _GET_HORIZONTAL_(self%id_SO4_1, SO4_1)
         _GET_HORIZONTAL_(self%id_SO4_2, SO4_2)
         _GET_HORIZONTAL_(self%id_SO4_3, SO4_3)
         _GET_HORIZONTAL_(self%id_H2S_1, H2S_1)
         _GET_HORIZONTAL_(self%id_H2S_2, H2S_2)
         _GET_HORIZONTAL_(self%id_H2S_3, H2S_3)
         _GET_HORIZONTAL_(self%id_S0_1, S0_1)
         _GET_HORIZONTAL_(self%id_S0_2, S0_2)
         _GET_HORIZONTAL_(self%id_S0_3, S0_3)

         ! Get oxygen and nitrate
         _GET_HORIZONTAL_(self%id_O2_1, O2_1)
         _GET_HORIZONTAL_(self%id_O2_2, O2_2)
         _GET_HORIZONTAL_(self%id_NO3_2, NO3_2)

         ! Get layer depths
         _GET_HORIZONTAL_(self%id_D1m, D1m)
         _GET_HORIZONTAL_(self%id_D2m, D2m)
         _GET_HORIZONTAL_(self%id_Dtot, Dtot)

         ! Get remineralization rate
         _GET_HORIZONTAL_(self%id_remin_rate, remin_rate)

         ! Calculate layer 3 thickness and concentration
         layer3_thickness = max(Dtot - D2m, 0.0001_rk)
         SO4_conc_3 = SO4_3 / layer3_thickness

         ! Electron acceptor cascade control
         ! Layer 3 is by definition anoxic (O2=0, NO3~0), so f_sulfate_red ~ 1
         f_anox = 0.5_rk * (1.0_rk - tanh((O2_2 - self%O2_threshold) * 100._rk))
         f_no_NO3 = 0.5_rk * (1.0_rk - tanh((NO3_2 - self%NO3_threshold) * 100._rk))
         f_sulfate_red = f_anox * f_no_NO3

         ! R1: Sulfate reduction in layer 3
         ! Rate proportional to organic matter remineralization
         R_sulfate_red = self%K_SO4_rd * SO4_conc_3 / (SO4_conc_3 + self%K_SO4_half) &
                       * remin_rate * self%stoich_S_C * f_sulfate_red

         ! R2: H2S oxidation in layer 1
         ! Convert depth-integrated O2 to concentration for rate calculation
         R_H2S_ox_1 = self%K_H2S_ox * H2S_1 * (O2_1/D1m) / ((O2_1/D1m) + self%K_O2_half)

         ! R3: S0 oxidation in layer 1
         R_S0_ox_1 = self%K_S0_ox * S0_1 * (O2_1/D1m) / ((O2_1/D1m) + self%K_O2_half)

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

Register the new modules:

```fortran
! Add use statements
use ersem_sulfur_cycle
use ersem_benthic_sulfur_cycle

! In the registration section, add:
call add(factory, 'sulfur_cycle', type_ersem_sulfur_cycle)
call add(factory, 'benthic_sulfur_cycle', type_ersem_benthic_sulfur_cycle)
```

### Modifications to oxygen.F90

Enforce non-negative oxygen:

```fortran
! Change the comment and behavior
! Note: negative oxygen concentrations are NO LONGER permitted.
! Oxygen debt is now represented explicitly by H2S.

! Ensure nonnegative constraint is applied
```

### Modifications to benthic_nitrogen_cycle.F90

Remove K6 dependency:

```fortran
! Remove:
! - id_K6_sms dependency
! - K6_calculator child model
! - K6-related calculations in do_bottom

! The electron acceptor cascade is now handled by:
! 1. Oxygen consumption (primary)
! 2. Denitrification (when O2 depleted)
! 3. Sulfate reduction (when O2 and NO3 depleted) - handled by benthic_sulfur_cycle
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
G2s_H2S:
  long_name: benthic hydrogen sulfide
  model: ersem/benthic_column_dissolved_matter
  parameters:
    composition: s_h2s
    last_layer: 3              # Present in all layers
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
    O2_threshold: 0.01         # O2 threshold for sulfate reduction
    NO3_threshold: 0.005       # NO3 threshold for sulfate reduction
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

def test_sulfate_reduction_requires_anoxia():
    """Sulfate reduction should only occur under anoxic conditions"""
    # Test with O2 > threshold: R_sulfate_red should be ~0
    # Test with O2 < threshold: R_sulfate_red should be > 0
    pass

def test_h2s_oxidation_requires_oxygen():
    """H2S oxidation should require oxygen"""
    # Test with O2 > 0: H2S should decrease
    # Test with O2 = 0: H2S should remain constant
    pass

def test_sulfur_mass_conservation():
    """Total sulfur should be conserved"""
    # SO4 + H2S + S0 = constant (in closed system)
    pass

def test_oxygen_consumption_stoichiometry():
    """Verify O2 consumption matches stoichiometry"""
    # R2: 0.5 mol O2 per mol H2S
    # R3: 1.5 mol O2 per mol S0
    pass

def test_blue_tide_scenario():
    """Test blue tide conditions"""
    # Under bottom water anoxia:
    # 1. H2S should accumulate in sediments
    # 2. H2S should flux to water column
    # 3. S0 should increase in water column
    pass
```

### Integration Tests

Create a test configuration in `testcases/fabm-ersem-sulfur-test.yaml` that:
1. Simulates a shallow coastal system
2. Induces bottom water hypoxia/anoxia
3. Verifies H2S production and S0 accumulation

## Implementation Phases

### Phase 1 (Current Proposal)

- Sulfur species: SO4, H2S, S0 only
- Reactions: R1-R3 (sulfate reduction, H2S oxidation, S0 oxidation)
- No iron interactions
- Estimated effort: 1-2 weeks

### Phase 2 (Future)

- Add thiosulfate (S2O3) as intermediate
- Add thiodenitrification (R4: H2S + NO3 -> S0 + N2)
- Temperature dependence
- pH effects on rates

### Phase 3 (Future)

- Add iron-sulfur species: FeS, FeS2
- Integrate with ERSEM iron module (-DIRON)
- Full BROM-style redox chemistry

## References

1. BROM model: https://github.com/fabm-model/fabm/tree/master/src/models/ihamocc/brom
2. Butenschön et al. (2016). ERSEM 15.06. Geosci. Model Dev., 9, 1293-1339
3. Jørgensen, B. B. (1982). Mineralization of organic matter in the sea bed. Nature, 296, 643-645
4. Canfield, D. E. et al. (1993). Pathways of organic carbon oxidation in three continental margin sediments. Mar. Geol., 113, 27-40

## Notes for Implementers

1. **Start with K6 compatibility**: Initially keep K6 in parallel with H2S to avoid breaking existing configurations. Remove K6 in a later cleanup phase.

2. **Test incrementally**:
   - First implement pelagic sulfur_cycle.F90 alone
   - Verify with simple box model
   - Then add benthic_sulfur_cycle.F90
   - Finally, modify existing modules

3. **Watch for unit conversions**:
   - Pelagic: mmol/m³ (concentration)
   - Benthic: mmol/m² (depth-integrated)
   - Layer concentrations need division by layer thickness

4. **Verify mass conservation**: Total sulfur (SO4 + H2S + S0) should be conserved in closed systems.

5. **Consult BROM source**: https://github.com/fabm-model/fabm/blob/master/src/models/ihamocc/brom/brom_redox.F90 for reference implementation details.
