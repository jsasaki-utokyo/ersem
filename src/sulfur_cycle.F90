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
   use ersem_sulfur_laws, only: srk, theta_q10, sulfur_extents, sulfur_nu   ! jsasaki 2026-11-23: isw_law = 1 (shared sulfur laws)

   implicit none
   private

   ! jsasaki 2026-11-23: isw_law = 1: reactive sulfur in the water (H2S + S0 + SO4a), a conserved quantity
   type(type_universal_standard_variable), parameter :: total_reactive_sulfur = type_universal_standard_variable( &
      name='total_reactive_sulfur', units='mmol S/m^3', aggregate_variable=.true., conserved=.true.)

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

      !---> jsasaki 2026-11-23: shared water/bed sulfur laws (TB-FVCOM box-ersem/docs/sulfur_unification.md, isw_law = 1)
      integer  :: isw_law
      real(rk) :: k_hs_tot, f_s0, k_s0ox, k_disp, h2s_stop, q10_hs, q10_s0, q10_disp, fw_hs, fw_s0
      integer  :: n_sub
      type(type_state_variable_id)      :: id_SO4a
      type(type_diagnostic_variable_id) :: id_SO4
      type(type_dependency_id)          :: id_temp, id_salt
      type(type_global_dependency_id)   :: id_hstep
      !<--- jsasaki 2026-11-23

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

      !---> jsasaki 2026-11-23: isw_law (TB-FVCOM box-ersem/docs/sulfur_unification.md section 4): 0 = the legacy
      ! Monod laws below, untouched; 1 = the shared water/bed laws of MUSE (mass action in both reactants, per-reaction
      ! Q10, MUSE parameter names and values), integrated over the host step by exact extents (ersem_sulfur_laws); the
      ! legacy rate keys must then be absent (FABM rejects unused keys)
      call self%get_parameter(self%isw_law, 'isw_law', '', &
           'sulfur laws (0: legacy Monod, 1: shared water/bed mass-action laws)', default=0, minimum=0, maximum=1)
      if (self%isw_law == 1) then
         call self%get_parameter(self%k_hs_tot, 'k_hs_tot', '(mmol/m^3)^-1/d', 'total aerobic sulfide oxidation (MUSE C2)')
         call self%get_parameter(self%f_s0, 'f_s0', '-', 'share of sulfide oxidation to S0', minimum=0.0_rk, maximum=1.0_rk)
         call self%get_parameter(self%k_s0ox, 'k_s0ox', '(mmol/m^3)^-1/d', 'S0 oxidation (MUSE C2)')
         call self%get_parameter(self%k_disp, 'k_disp', '1/d', 'S0 disproportionation (MUSE C6)', minimum=0.0_rk)
         call self%get_parameter(self%h2s_stop, 'h2s_stop', 'mmol/m^3', 'sulfide stop of disproportionation')
         call self%get_parameter(self%q10_hs, 'q10_hs', '-', 'Q10 of sulfide oxidation (MUSE q10(10) = q10(11))')
         call self%get_parameter(self%q10_s0, 'q10_s0', '-', 'Q10 of S0 oxidation (MUSE q10(12))')
         call self%get_parameter(self%q10_disp, 'q10_disp', '-', 'Q10 of disproportionation (MUSE q10(19))')
         call self%get_parameter(self%fw_hs, 'fw_hs', '-', 'water activity factor of sulfide oxidation', default=1.0_rk, &
              minimum=0.0_rk)
         call self%get_parameter(self%fw_s0, 'fw_s0', '-', 'water activity factor of S0 oxidation', default=1.0_rk, &
              minimum=0.0_rk)
         call self%get_parameter(self%n_sub, 'n_sub', '', 'Strang sub-passes per host step', default=32, minimum=1)
         ! preflight (section 4.1): no O2 inhibition in the disproportionation law, so it must stay slow
         if (self%k_disp > 0.01_rk) call self%fatal_error('initialize', 'isw_law = 1 needs k_disp <= 0.01 d-1')
         if (self%h2s_stop <= 0.0_rk) call self%fatal_error('initialize', 'isw_law = 1 needs h2s_stop > 0')
         call self%register_dependency(self%id_temp, standard_variables%temperature)
         call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
         ! the host step length (s; 0 outside a trial = the instantaneous law), supplied by the box driver
         call self%register_global_dependency(self%id_hstep, type_global_standard_variable(name='box_step_length', units='s'))
      end if
      !<--- jsasaki 2026-11-23

      ! Get parameters (default values from BROM)
      if (self%isw_law == 0) then   ! jsasaki 2026-11-23: legacy keys only with isw_law = 0
      call self%get_parameter(self%K_H2S_ox, 'K_H2S_ox', '1/d', &
           'H2S oxidation rate constant', default=0.5_rk)
      call self%get_parameter(self%K_S0_ox, 'K_S0_ox', '1/d', &
           'S0 oxidation rate constant', default=0.02_rk)
      call self%get_parameter(self%K_O2_half, 'K_O2_half', 'mmol/m^3', &
           'half-saturation O2 for oxidation', default=1.0_rk)
      end if   ! jsasaki 2026-11-23

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
      if (self%isw_law == 0) &   ! jsasaki 2026-11-23: legacy key only with isw_law = 0
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

      !---> jsasaki 2026-11-23: isw_law = 1 (sulfur_unification.md sections 4.1, 5): sulfate anomaly SO4a (departure of
      ! water sulfate from the salinity reference; may be negative), diagnostic SO4 = SO4_ref(S) + SO4a, and the
      ! conserved quantity total_reactive_sulfur = H2S + S0 + SO4a (closed with the MUSE bed's element 5)
      if (self%isw_law == 1) then
         if (self%isw_ta /= 1) call self%fatal_error('initialize', 'isw_law = 1 needs isw_ta = 1')
         if (self%K_S0_sink /= 0.0_rk) call self%fatal_error('initialize', 'isw_law = 1 needs K_S0_sink = 0')
         call self%register_state_variable(self%id_SO4a, 'SO4a', 'mmol S/m^3', &
              'sulfate anomaly from the salinity reference', minimum=-1.e20_rk)
         call self%register_diagnostic_variable(self%id_SO4, 'SO4', 'mmol S/m^3', &
              'sulfate (salinity reference 28240 S/35 x 1.025 + SO4a)', source=source_do)
         call self%add_to_aggregate_variable(total_reactive_sulfur, self%id_H2S)
         call self%add_to_aggregate_variable(total_reactive_sulfur, self%id_S0)
         call self%add_to_aggregate_variable(total_reactive_sulfur, self%id_SO4a)
      end if
      !<--- jsasaki 2026-11-23

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_sulfur_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: H2S, S0, O2
      real(rk) :: R_H2S_ox, R_S0_ox, R_S0_sink, R_S0_disp
      real(rk) :: f_O2, f_anox
      ! jsasaki 2026-11-23: isw_law = 1 locals
      real(rk) :: T, Sal, SO4a, hstep, hd, k_hs, k_s0, k_dis, nu, e_hs, e_s0, e_dis
      real(rk) :: dH2S, dS0, dO2, dSO4, dTA
      real(rk), parameter :: eps4 = 4.0_rk * epsilon(1.0_rk)

      _LOOP_BEGIN_
         _GET_(self%id_H2S, H2S)
         _GET_(self%id_S0, S0)
         _GET_(self%id_O2, O2)

         !---> jsasaki 2026-11-23: isw_law = 1 (sulfur_unification.md sections 3, 4.2): the shared laws, integrated
         ! over the host step h by n_sub Strang sub-passes of exact extents; the rate passed to the host is the summed
         ! extent / h. h = 0 (outside a trial): the instantaneous laws. A per-stock clamp keeps the host's Euler update
         ! of H2S, S0 and O2 from rounding below zero (changes the closure by at most 4 eps of the stock)
         if (self%isw_law == 1) then
            _GET_(self%id_temp, T)
            _GET_(self%id_salt, Sal)
            _GET_(self%id_SO4a, SO4a)
            _GET_GLOBAL_(self%id_hstep, hstep)
            if (hstep < 0.0_rk) call self%fatal_error('do', 'box_step_length < 0')
            O2 = max(0.0_rk, O2)
            H2S = max(0.0_rk, H2S)
            S0 = max(0.0_rk, S0)
            k_hs = theta_q10(self%q10_hs, T) * self%k_hs_tot * self%fw_hs
            k_s0 = theta_q10(self%q10_s0, T) * self%k_s0ox * self%fw_s0
            k_dis = theta_q10(self%q10_disp, T) * self%k_disp
            nu = sulfur_nu(self%f_s0)
            if (hstep > 0.0_rk) then
               hd = hstep / 86400.0_rk
               call sulfur_extents(H2S, S0, O2, k_hs, self%f_s0, k_s0, k_dis, self%h2s_stop, hd, self%n_sub, &
                                   e_hs, e_s0, e_dis)
               R_H2S_ox = e_hs / hd
               R_S0_ox = e_s0 / hd
               R_S0_disp = e_dis / hd
            else
               hd = 0.0_rk
               R_H2S_ox = k_hs * H2S * O2
               R_S0_ox = k_s0 * S0 * O2
               R_S0_disp = k_dis * S0 * max(0.0_rk, 1.0_rk - H2S / self%h2s_stop)
            end if
            R_S0_sink = 0.0_rk
            dH2S = -R_H2S_ox + 0.75_rk * R_S0_disp
            dS0 = self%f_s0 * R_H2S_ox - R_S0_ox - R_S0_disp
            dO2 = -nu * R_H2S_ox - 1.5_rk * R_S0_ox
            dSO4 = (1.0_rk - self%f_s0) * R_H2S_ox + R_S0_ox + 0.25_rk * R_S0_disp
            dTA = -2.0_rk * (1.0_rk - self%f_s0) * R_H2S_ox - 2.0_rk * R_S0_ox - 0.5_rk * R_S0_disp
            if (hd > 0.0_rk) then
               dH2S = max(dH2S, -H2S * (1.0_rk - eps4) / hd)
               dS0 = max(dS0, -S0 * (1.0_rk - eps4) / hd)
               dO2 = max(dO2, -O2 * (1.0_rk - eps4) / hd)
            end if
            _SET_ODE_(self%id_H2S, dH2S)
            _SET_ODE_(self%id_S0, dS0)
            _SET_ODE_(self%id_O2, dO2)
            _SET_ODE_(self%id_SO4a, dSO4)
            _SET_ODE_(self%id_TA, dTA)
            if (self%isw_ledger == 1) then
               _SET_DIAGNOSTIC_(self%id_ledger_sulf_H2S, dH2S)
               _SET_DIAGNOSTIC_(self%id_ledger_sulf_S0, dS0)
               _SET_DIAGNOSTIC_(self%id_ledger_sulf_O2, dO2)
               _SET_DIAGNOSTIC_(self%id_ledger_sulf_TA, dTA)
            end if
            _SET_DIAGNOSTIC_(self%id_SO4, 28240.0_rk * Sal / 35.0_rk * 1.025_rk + SO4a)
         else
         !<--- jsasaki 2026-11-23

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
         end if   ! jsasaki 2026-11-23: isw_law

         ! Set diagnostics
         _SET_DIAGNOSTIC_(self%id_R_H2S_ox, R_H2S_ox)
         _SET_DIAGNOSTIC_(self%id_R_S0_ox, R_S0_ox)
         _SET_DIAGNOSTIC_(self%id_R_S0_sink, R_S0_sink)
         _SET_DIAGNOSTIC_(self%id_R_S0_disp, R_S0_disp)

      _LOOP_END_
   end subroutine do

end module ersem_sulfur_cycle
