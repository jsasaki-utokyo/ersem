#include "fabm_driver.h"

!-----------------------------------------------------------------------
! Pelagic sulfur cycle module for ERSEM
!
! Handles H2S and S0 dynamics in the water column.
! Simplified implementation (Option B): H2S + S0 only, no explicit SO4.
!
! Reactions:
!   R1: H2S + 0.5 O2 -> S0 + H2O           (H2S oxidation)
!   R2: S0 + 1.5 O2 + H2O -> SO4 + 2H+     (S0 oxidation, SO4 not tracked)
!   R3: S0 -> settling/sinking             (particle settling)
!   R4: 4S0 + 4H2O -> 3H2S + SO4 + 2H+     (disproportionation under anoxia)
!
! S0 Removal Mechanisms:
!   Elemental sulfur (S0) can accumulate if only oxidation is considered.
!   Additional removal pathways are essential:
!   - Settling: S0 particles sink out of the water column
!   - Disproportionation: Under anoxic conditions, S0 reacts with water
!     to produce H2S and sulfate (microbially mediated)
!
! References:
!   - BROM model (Yakushev et al. 2017)
!   - Rate constants from fabm-niva-brom
!   - Thamdrup et al. (1994) - S0 disproportionation
!-----------------------------------------------------------------------

module ersem_sulfur_cycle

   use fabm_types
   use ersem_shared

   implicit none
   private

   type, extends(type_base_model), public :: type_ersem_sulfur_cycle
      ! State variable IDs
      type(type_state_variable_id) :: id_H2S, id_S0
      type(type_state_variable_id) :: id_O2
      type(type_state_variable_id) :: id_TA   ! alkalinity (isw_ta = 1 only)

      ! Diagnostic variable IDs
      type(type_diagnostic_variable_id) :: id_R_H2S_ox, id_R_S0_ox
      type(type_diagnostic_variable_id) :: id_R_S0_sink, id_R_S0_disp

      ! Parameters
      real(rk) :: K_H2S_ox      ! H2S oxidation rate constant (1/d)
      real(rk) :: K_S0_ox       ! S0 oxidation rate constant (1/d)
      real(rk) :: K_O2_half     ! Half-saturation O2 for oxidation (mmol/m3)
      real(rk) :: K_S0_sink     ! S0 settling/sinking rate (1/d)
      real(rk) :: K_S0_disp     ! S0 disproportionation rate (1/d)
      integer  :: isw_ta        ! 1 = charge the alkalinity of S0 -> SO4 (docs/114, 2026-09-11); 0 = legacy

      ! ledger counters (isw_ledger = 1 only; nippon-steel docs/119 FC1)
      integer  :: isw_ledger
      type(type_diagnostic_variable_id) :: id_ledger_sulf_H2S, id_ledger_sulf_S0
      type(type_diagnostic_variable_id) :: id_ledger_sulf_O2, id_ledger_sulf_TA

   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class(type_ersem_sulfur_cycle), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! Set time unit to d-1 (ERSEM convention)
      self%dt = 86400._rk

      ! Get parameters (default values from BROM)
      call self%get_parameter(self%K_H2S_ox, 'K_H2S_ox', '1/d', &
           'H2S oxidation rate constant', default=0.5_rk)
      call self%get_parameter(self%K_S0_ox, 'K_S0_ox', '1/d', &
           'S0 oxidation rate constant', default=0.02_rk)
      call self%get_parameter(self%K_O2_half, 'K_O2_half', 'mmol/m^3', &
           'half-saturation O2 for oxidation', default=1.0_rk)

      ! S0 removal parameters
      ! K_S0_sink: Elemental sulfur particles settle out of the water column.
      ! S0 forms colloidal or particulate aggregates that sink.
      ! Typical settling rate: 0.1-0.5 /d (equivalent to ~1-5 m/d for 10m depth)
      call self%get_parameter(self%K_S0_sink, 'K_S0_sink', '1/d', &
           'S0 settling/sinking rate', default=0.2_rk)

      ! K_S0_disp: Disproportionation under anoxic conditions
      ! 4S0 + 4H2O -> 3H2S + SO4 + 2H+ (microbially mediated)
      ! This reaction occurs when O2 is depleted and recycles S0 back to H2S.
      ! Typical rate: 0.05-0.2 /d under fully anoxic conditions
      call self%get_parameter(self%K_S0_disp, 'K_S0_disp', '1/d', &
           'S0 disproportionation rate under anoxia', default=0.1_rk)

      ! ALKALINITY of the sulfur chain (2026-09-11, nippon-steel docs/114,
      ! Fable review of EM-ES). The benthic module charges -2 TA per S when S0
      ! is oxidised to sulfate (returning the +2 TA of sulfate reduction), but
      ! this pelagic module oxidised S0 with oxygen and NO alkalinity term, so
      ! S0 leaving the bed was reoxidised in the water without returning the
      ! sulfate-reduction alkalinity: a charge-balance defect. isw_ta = 1 adds
      !   S0 + 1.5 O2 + H2O -> SO4(2-) + 2H+        : -2   TA per S0 oxidised
      !   4 S0 + 4 H2O -> 3 H2S + SO4(2-) + 2H+     : -0.5 TA per S0 disproportionated
      ! (H2S <-> S0 is TA-neutral in the sulfate-reduction convention used
      ! throughout). isw_ta = 0 (default) is the legacy behaviour, bit-identical.
      call self%get_parameter(self%isw_ta, 'isw_ta', '', &
           'charge the alkalinity of S0 oxidation and disproportionation (0: legacy, 1: on)', &
           default=0, minimum=0, maximum=1)
      if (self%isw_ta == 1) &
         call self%register_state_dependency(self%id_TA, standard_variables%alkalinity_expressed_as_mole_equivalent)

      ! LEDGER COUNTERS (2026-09-16, nippon-steel docs/119 §4.1 and §4.3 FC1).
      ! isw_ledger = 1 registers one diagnostic per source term this module
      ! writes to a ledger target (DIC, TA, O2, NO3, NH4, N2, PO4, Si, H2S, S0,
      ! CaCO3; organic matter only where it crosses between the water and the
      ! bed). Each value is the expression passed to the _SET_ODE_ /
      ! _SET_BOTTOM_ODE_ / _SET_BOTTOM_EXCHANGE_ call beside it, in this
      ! module's time unit (per day), with the sign applied to the target, so
      ! post-run gates can compare what the code applies against independently
      ! specified stoichiometry and against the state changes. Diagnostics only:
      ! never read by any reaction. isw_ledger = 0 (default) registers and
      ! computes nothing, bit-identical to the legacy build.
      call self%get_parameter(self%isw_ledger, 'isw_ledger', '', &
           'ledger counters: diagnostics of the source terms applied (0: off, 1: on)', &
           default=0, minimum=0, maximum=1)
      if (self%isw_ledger == 1) then
         call self%register_diagnostic_variable(self%id_ledger_sulf_H2S, 'ledger_sulf_H2S', &
              'mmol S/m^3/d', 'ledger: H2S oxidation and S0 disproportionation -> hydrogen sulfide', source=source_do)
         call self%register_diagnostic_variable(self%id_ledger_sulf_S0, 'ledger_sulf_S0', &
              'mmol S/m^3/d', 'ledger: H2S oxidation, S0 oxidation, settling and disproportionation -> elemental sulfur', &
              source=source_do)
         call self%register_diagnostic_variable(self%id_ledger_sulf_O2, 'ledger_sulf_O2', &
              'mmol O_2/m^3/d', 'ledger: H2S and S0 oxidation -> oxygen', source=source_do)
         if (self%isw_ta == 1) &
            call self%register_diagnostic_variable(self%id_ledger_sulf_TA, 'ledger_sulf_TA', &
                 'mmol/m^3/d', 'ledger: S0 oxidation and disproportionation -> total alkalinity', source=source_do)
      end if

      ! Register state variables
      call self%register_state_variable(self%id_H2S, 'H2S', 'mmol S/m^3', &
           'hydrogen sulfide', minimum=0.0_rk)
      call self%register_state_variable(self%id_S0, 'S0', 'mmol S/m^3', &
           'elemental sulfur', minimum=0.0_rk)

      ! Register dependency on oxygen
      call self%register_state_dependency(self%id_O2, 'O2', 'mmol O_2/m^3', 'oxygen')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_R_H2S_ox, 'R_H2S_ox', &
           'mmol S/m^3/d', 'H2S oxidation rate', source=source_do)
      call self%register_diagnostic_variable(self%id_R_S0_ox, 'R_S0_ox', &
           'mmol S/m^3/d', 'S0 oxidation rate', source=source_do)
      call self%register_diagnostic_variable(self%id_R_S0_sink, 'R_S0_sink', &
           'mmol S/m^3/d', 'S0 settling rate', source=source_do)
      call self%register_diagnostic_variable(self%id_R_S0_disp, 'R_S0_disp', &
           'mmol S/m^3/d', 'S0 disproportionation rate', source=source_do)

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_sulfur_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: H2S, S0, O2
      real(rk) :: R_H2S_ox, R_S0_ox, R_S0_sink, R_S0_disp
      real(rk) :: f_O2, f_anox

      _LOOP_BEGIN_
         _GET_(self%id_H2S, H2S)
         _GET_(self%id_S0, S0)
         _GET_(self%id_O2, O2)

         ! Ensure non-negative values
         O2 = max(0.0_rk, O2)
         H2S = max(0.0_rk, H2S)
         S0 = max(0.0_rk, S0)

         ! Oxygen limitation (Michaelis-Menten)
         f_O2 = O2 / (O2 + self%K_O2_half)

         ! Anoxia factor (inverse of f_O2, high when O2 is low)
         f_anox = self%K_O2_half / (O2 + self%K_O2_half)

         ! ============================================================
         ! H2S oxidation
         ! ============================================================
         ! R1: H2S + 0.5 O2 -> S0 + H2O
         R_H2S_ox = self%K_H2S_ox * H2S * f_O2

         ! ============================================================
         ! S0 removal mechanisms
         ! ============================================================
         ! R2: S0 + 1.5 O2 + H2O -> SO4 + 2H+ (oxidation when O2 present)
         R_S0_ox = self%K_S0_ox * S0 * f_O2

         ! R3: S0 settling/sinking (always active, particle settling)
         ! Elemental sulfur forms colloidal/particulate aggregates that sink
         R_S0_sink = self%K_S0_sink * S0

         ! R4: S0 disproportionation under anoxic conditions
         ! 4S0 + 4H2O -> 3H2S + SO4 + 2H+ (microbially mediated)
         ! Only active when O2 is low, produces H2S from S0
         R_S0_disp = self%K_S0_disp * S0 * f_anox

         ! ============================================================
         ! Set ODEs
         ! ============================================================
         ! H2S: consumed by oxidation, produced by disproportionation (3/4 of S0)
         _SET_ODE_(self%id_H2S, -R_H2S_ox + 0.75_rk * R_S0_disp)

         ! S0: produced by H2S oxidation, consumed by oxidation/settling/disproportionation
         _SET_ODE_(self%id_S0, R_H2S_ox - R_S0_ox - R_S0_sink - R_S0_disp)

         ! O2: consumed by H2S oxidation and S0 oxidation
         _SET_ODE_(self%id_O2, -0.5_rk * R_H2S_ox - 1.5_rk * R_S0_ox)

         ! Alkalinity of sulfate regeneration (isw_ta = 1; see initialize)
         if (self%isw_ta == 1) _SET_ODE_(self%id_TA, -2.0_rk * R_S0_ox - 0.5_rk * R_S0_disp)

         ! Ledger counters (isw_ledger = 1; see initialize): the terms passed above, verbatim
         if (self%isw_ledger == 1) then
            _SET_DIAGNOSTIC_(self%id_ledger_sulf_H2S, -R_H2S_ox + 0.75_rk * R_S0_disp)
            _SET_DIAGNOSTIC_(self%id_ledger_sulf_S0, R_H2S_ox - R_S0_ox - R_S0_sink - R_S0_disp)
            _SET_DIAGNOSTIC_(self%id_ledger_sulf_O2, -0.5_rk * R_H2S_ox - 1.5_rk * R_S0_ox)
            if (self%isw_ta == 1) _SET_DIAGNOSTIC_(self%id_ledger_sulf_TA, -2.0_rk * R_S0_ox - 0.5_rk * R_S0_disp)
         end if

         ! Set diagnostics
         _SET_DIAGNOSTIC_(self%id_R_H2S_ox, R_H2S_ox)
         _SET_DIAGNOSTIC_(self%id_R_S0_ox, R_S0_ox)
         _SET_DIAGNOSTIC_(self%id_R_S0_sink, R_S0_sink)
         _SET_DIAGNOSTIC_(self%id_R_S0_disp, R_S0_disp)

      _LOOP_END_
   end subroutine do

end module ersem_sulfur_cycle
