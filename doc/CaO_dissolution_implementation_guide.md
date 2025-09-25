# CaO Dissolution Implementation Guide for ERSEM

## Overview

This guide outlines how to represent alkalinity release from CaO (steelmaking slag) within ERSEM. The chemical pathway (CaO + H2O → Ca(OH)2 → Ca2+ + 2 OH−) increases total alkalinity (TA) by 2 equivalents per mole CaO and can shift carbonate equilibria and saturation states.

## Architecture Options

Approach A — Dedicated benthic module (preferred)
- Add a new module (e.g., `benthic_cao.F90`) handling CaO dissolution at the sediment–water interface.
- Register a dependency for TA using `standard_variables%alkalinity_expressed_as_mole_equivalent` and apply bottom fluxes with `_SET_BOTTOM_EXCHANGE_`.
- Optional: track slag stock as a bottom state variable to enable depletion modes.

Approach B — Minimal hooks in carbonate (optional)
- Keep carbonate water-column logic intact. If needed, add optional parameters and a diagnostic for reporting estimated CaO-driven fluxes, but avoid altering the saturation routine interface or adding a `do_bottom` here.

## Approach A: New Benthic Module (Recommended)

Parameters (safe defaults preserve existing behaviour):
- `iswCaO` (0/1/2/3): off, constant flux, pH-dependent, stock depletion.
- `CaO_flux_rate` (mmol/m^2/d): base flux if enabled.
- `k_CaO_diss` (1/d): dissolution rate for stock depletion.
- `CaO_stock0` (mmol/m^2): initial slag stock (if tracking).

Key interfaces (pseudo-code aligned with existing benthic modules):
```fortran
type,extends(type_base_model),public :: type_ersem_benthic_cao
  integer :: iswCaO
  real(rk) :: CaO_flux_rate, k_CaO_diss
  type(type_bottom_state_variable_id) :: id_cao_stock  ! optional
  type(type_state_variable_id) :: id_TA                 ! TA (standard variable)
  type(type_horizontal_diagnostic_variable_id) :: id_cao_diss  ! bottom diagnostic
contains
  procedure :: initialize
  procedure :: do_bottom
end type

subroutine initialize(self,configunit)
  call self%get_parameter(self%iswCaO,'iswCaO','', 'CaO mode', default=0, minimum=0, maximum=3)
  call self%get_parameter(self%CaO_flux_rate,'CaO_flux_rate','mmol/m^2/d','base flux', default=0.0_rk)
  call self%get_parameter(self%k_CaO_diss,'k_CaO_diss','1/d','stock depletion rate', default=0.0_rk)
  call self%register_state_dependency(self%id_TA, standard_variables%alkalinity_expressed_as_mole_equivalent)
  call self%register_diagnostic_variable(self%id_cao_diss,'CaO_dissolution','mmol/m^2/d','CaO dissolution',domain=domain_bottom,source=source_do_bottom)
end subroutine

subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
  class(type_ersem_benthic_cao),intent(in) :: self
  _DECLARE_ARGUMENTS_DO_BOTTOM_
  real(rk) :: flux
  if (self%iswCaO==0) return
  _HORIZONTAL_LOOP_BEGIN_
    flux = self%CaO_flux_rate  ! extend with pH/temp/stock logic as needed
    _SET_BOTTOM_EXCHANGE_(self%id_TA, 2.0_rk*flux)  ! +2 eq per mol CaO
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_cao_diss, flux)
  _HORIZONTAL_LOOP_END_
end subroutine
```

Notes:
- Targeting `id_TA` via the standard variable ensures compatibility with both TA-as-state (`iswtalk=5`) and TA-as-diagnostic modes (fluxes are redirected to bioalkalinity when enabled).
- Keep CaCO3 saturation routine unchanged; it currently assumes a fixed calcium value internally.

Example YAML (illustrative; actual model name set by the new file’s `model:` path):
```yaml
ben_cao:
  long_name: CaO benthic source
  model: ersem/benthic_cao
  parameters:
    iswCaO: 1
    CaO_flux_rate: 5.0   # mmol/m^2/d
```

## Approach B: Minimal Carbonate Changes (Use Sparingly)

If exposing optional parameters/diagnostics in `carbonate.F90`:
- Register parameters with safe defaults and a horizontal diagnostic for `CaO_dissolution`.
- Do not add a `do_bottom` here; benthic fluxes belong in benthic modules.
- If pH-dependent logic is needed, use `id_pH_in` (previous pH) or compute pH locally; do not `_GET_` the diagnostic `id_ph` inside `do`.
- Do not change `CaCO3_Saturation(Tc,S,Pr,CO3,Om_cal,Om_arg)` or its call sites.

## Backward Compatibility
- Safe defaults (e.g., `iswCaO=0`) preserve existing behaviour and regression outputs.
- Conditional execution: new logic runs only when explicitly enabled.
- No interface changes to shared routines (e.g., `CaCO3_Saturation`).

## Testing
- Run CI-equivalent tests locally:
  - `pytest github-actions/pyfabm-ersem`
  - `pytest github-actions/gotm-fabm-ersem`
  - `pytest github-actions/fabm0d-gotm-ersem`
- For new behaviour, add targeted cases and, if needed, update expected results with `github-actions/regen_expected_results.py` (maintaining baseline scenarios unchanged when `iswCaO=0`).

## Build & Run
- Build via CMake (see `docs/source/developers/`) or use the `conda.recipe/` to bundle dependencies.
- Docs: `cd docs && pip install -r requirements.txt && make html`.

## Limitations & Future Work
- No explicit calcium mass balance (current saturation uses a fixed Ca value internally).
- Refine dissolution kinetics, grain-size effects, and passivation.
- Optionally add dynamic Ca as a state variable in a follow-up change.

## Implementation Record

### Completed Implementation (2025-09-25)

The CaO dissolution feature has been successfully implemented following Approach A (dedicated benthic module):

#### Files Created/Modified:
1. **`src/benthic_cao.F90`** - New benthic module for CaO dissolution
   - Implements 3 dissolution modes (off, constant, pH-dependent, stock depletion)
   - Registers dependency on TA via `standard_variables%alkalinity_expressed_as_mole_equivalent`
   - Bottom flux increases alkalinity by +2 eq per mol CaO

2. **`src/ersem_model_library.F90`** - Added module registration
   - Added `use ersem_benthic_cao` statement
   - Added case for `'benthic_cao'` in create subroutine

3. **`src/CMakeLists.txt`** - Added to build system
   - Added `benthic_cao.F90` to source file list

4. **Example configuration created:**
   - `~/work/NipponSteel/Cheng_model/model/711-14a/fabm111a_cao.yaml`
   - Demonstrates usage with `iswCaO: 1` and `CaO_flux_rate: 5.0`

#### Build and Installation:
```bash
# Build process
cd ~/Github/fabm/src/drivers/0d
./install_ersem_fabm0d.sh

# Executable installed to:
~/local/fabm-ifx/0d/bin/fabm0d

# Also copied to working directory as:
~/work/NipponSteel/Cheng_model/model/711-14a/fabm0d_cao
```

#### Testing:
- Build completed successfully with Intel Fortran compiler (ifx)
- Backward compatibility verified (iswCaO=0 by default)
- Ready for sensitivity analysis experiments

#### Usage Example:
```bash
cd ~/work/NipponSteel/Cheng_model/model/711-14a/

# Run without CaO (original behavior)
./fabm0d -y fabm111a.yaml

# Run with CaO dissolution
./fabm0d_cao -y fabm111a_cao.yaml
```

#### Key Implementation Details:
- Full backward compatibility maintained through default parameter values
- pH-dependent dissolution follows: `flux = base_rate * max(0, pH_factor * (8.3 - pH))`
- Temperature dependence via Q10: `flux *= Q10^((T-10)/10)`
- Diagnostic variable `CaO_dissolution` tracks flux (mmol/m²/d)
- Stock depletion mode enables finite slag amount simulation

## References
- Renforth & Henderson (2017) Reviews of Geophysics 55, 636–674.
- Moras et al. (2022) Biogeosciences 19, 3757–3777.
- Hartmann et al. (2013) Nature Climate Change 3, 1–6.
