#include "fabm_driver.h"

! =========================================================================
! dom_decay.F90 — Bacteria-independent DOM oxidation module for ERSEM
!
! First-order dissolved organic matter oxidation with direct O2 consumption.
! Represents background microbial activity that is always present,
! independent of the explicit B1 bacteria state variable.
!
! Both dom_decay and B1 act on the same DOM pools simultaneously:
!   - When B1 is abundant (surface), B1 dominates DOM consumption
!   - When B1 is scarce (mid-water, deep), dom_decay provides baseline
!
! This follows BROM's K_DON_ox design (Yakushev et al. 2017).
! Each FABM instance targets one DOM class (R1, R2, or R3).
!
! jsasaki 2026-03-20: Initial implementation
! =========================================================================

module ersem_dom_decay

   use fabm_types
   use ersem_shared
   use benthos_shared_laws, only: th_q10   ! jsasaki 2026-11-24: isw_law = 1: the pure Q10 law shared with the bed DOC decay

   implicit none

   private

   type, extends(type_base_model), public :: type_ersem_dom_decay
      ! Target DOM state variable
      type (type_state_variable_id) :: id_RDc   ! DOM carbon
      type (type_state_variable_id) :: id_RDn   ! DOM nitrogen (optional, R1 only)
      type (type_state_variable_id) :: id_RDp   ! DOM phosphorus (optional, R1 only)

      ! Inorganic pools
      type (type_state_variable_id) :: id_O2o   ! dissolved oxygen
      type (type_state_variable_id) :: id_O3c   ! dissolved inorganic carbon
      type (type_state_variable_id) :: id_N4n   ! ammonium
      type (type_state_variable_id) :: id_N1p   ! phosphate
      type (type_state_variable_id) :: id_TA    ! total alkalinity

      ! Environmental dependencies
      type (type_dependency_id) :: id_ETW       ! temperature

      ! Diagnostic
      type (type_diagnostic_variable_id) :: id_rate  ! oxidation rate

      ! Parameters
      real(rk) :: k_ox        ! specific oxidation rate at Tref (1/d)
      real(rk) :: q10         ! Q10 temperature coefficient
      real(rk) :: Tref        ! reference temperature (deg C)
      real(rk) :: K_O2        ! O2 half-saturation (mmol O2/m3)
      real(rk) :: ur_O2       ! O2 per C oxidized (mmol O2/mg C)
      ! jsasaki 2026-11-05: switch so that the rates declared in 1/d act per day (FABM rates are per second unless dt is set)
      logical  :: rates_per_day
      !---> jsasaki 2026-11-24: isw_law (muse docs/carb_s10m/rev52/design.md T01, T02): 0 (default) = the former behaviour, bit identical
      ! (modified Q10 with heat cut; DIC produced with 12.0 mg C per mmol C); 1 = unified law (pure Q10 th_q10, as the bed DOC decay
      ! of MUSE; DIC produced with ersem_shared's CMass = 12.011, which conserves carbon)
      integer  :: isw_law
      !<--- jsasaki 2026-11-24
   contains
      procedure :: initialize
      procedure :: do
   end type

   ! jsasaki 2026-11-24: renamed from CMass (it shadowed ersem_shared's CMass = 12.011); used only with isw_law = 0 (former behaviour)
   real(rk), parameter :: CMass_legacy = 12.0_rk  ! mg C per mmol C

contains

   subroutine initialize(self, configunit)
      class(type_ersem_dom_decay), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! --- Parameters ---
      call self%get_parameter(self%k_ox, 'k_ox', '1/d', &
         'specific oxidation rate at reference temperature', default=0.01_rk)
      call self%get_parameter(self%q10, 'q10', '-', &
         'Q10 temperature coefficient', default=2.0_rk)
      call self%get_parameter(self%Tref, 'Tref', 'degrees_Celsius', &
         'reference temperature', default=20.0_rk)
      call self%get_parameter(self%K_O2, 'K_O2', 'mmol O_2/m^3', &
         'half-saturation oxygen concentration', default=5.0_rk)
      call self%get_parameter(self%ur_O2, 'ur_O2', 'mmol O_2/mg C', &
         'oxygen consumed per carbon oxidized', default=0.1_rk)
      !---> jsasaki 2026-11-05: k_ox is declared in 1/d but self%dt was never set, so FABM applied it per second (86400x too fast;
      ! found by the S47 pfAc2 audit: lit-arm water DOC oxidised within ~1 s).
      ! jsasaki 2026-11-05: default switched to .true. (owner decision 2026-11-05): the declared 1/d now acts per day; set
      ! rates_per_day: false only to reproduce runs made before this date.
      call self%get_parameter(self%rates_per_day, 'rates_per_day', '', &
         'apply k_ox per day (true) instead of per second (false, former behaviour before 2026-11-05)', default=.true.)
      if (self%rates_per_day) self%dt = 86400._rk
      !<--- jsasaki 2026-11-05
      !---> jsasaki 2026-11-24: isw_law: temperature law and carbon mass (default 0 reproduces every earlier run of every host)
      call self%get_parameter(self%isw_law, 'isw_law', '', &
         'temperature law and carbon mass (0: modified Q10 with heat cut, 12.0 mg C/mmol C (former); 1: pure Q10, 12.011)', default=0)
      if (self%isw_law /= 0 .and. self%isw_law /= 1) call self%fatal_error('initialize', 'isw_law must be 0 or 1')
      !<--- jsasaki 2026-11-24

      ! --- Target DOM ---
      call self%register_state_dependency(self%id_RDc, 'RDc', 'mg C/m^3', &
         'dissolved organic carbon')
      call self%register_state_dependency(self%id_RDn, 'RDn', 'mmol N/m^3', &
         'dissolved organic nitrogen', required=.false.)
      call self%register_state_dependency(self%id_RDp, 'RDp', 'mmol P/m^3', &
         'dissolved organic phosphorus', required=.false.)

      ! --- Inorganic pools ---
      call self%register_state_dependency(self%id_O2o, 'O2o', 'mmol O_2/m^3', &
         'dissolved oxygen')
      call self%register_state_dependency(self%id_O3c, 'O3c', 'mmol C/m^3', &
         'dissolved inorganic carbon')
      call self%register_state_dependency(self%id_N4n, 'N4n', 'mmol N/m^3', &
         'ammonium', required=.false.)
      call self%register_state_dependency(self%id_N1p, 'N1p', 'mmol P/m^3', &
         'phosphate', required=.false.)
      call self%register_state_dependency(self%id_TA, &
         standard_variables%alkalinity_expressed_as_mole_equivalent)

      ! --- Environment ---
      call self%register_dependency(self%id_ETW, standard_variables%temperature)

      ! --- Diagnostic ---
      call self%register_diagnostic_variable(self%id_rate, 'rate', 'mg C/m^3/d', &
         'DOM oxidation rate', output=output_time_step_averaged)

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_dom_decay), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ETW, O2o
      real(rk) :: RDc, RDn, RDp
      real(rk) :: et, fO2, k_eff
      real(rk) :: fOx_c, fOx_n, fOx_p

      _LOOP_BEGIN_

         ! --- Temperature dependence ---
         _GET_(self%id_ETW, ETW)
         !---> jsasaki 2026-11-24: isw_law = 1: pure Q10 (law P of the unification table), exactly 1 at Tref
         if (self%isw_law == 1) then
            et = th_q10(ETW, self%q10, self%Tref)
         else
            et = max(0.0_rk, self%q10**((ETW - self%Tref) / 10.0_rk) &
                 - self%q10**((ETW - 32.0_rk) / 3.0_rk))
         end if
         !<--- jsasaki 2026-11-24

         ! --- Oxygen limitation ---
         _GET_(self%id_O2o, O2o)
         O2o = max(0.0_rk, O2o)
         fO2 = O2o / (O2o + self%K_O2)

         ! --- Effective rate ---
         k_eff = self%k_ox * et * fO2

         ! --- DOM carbon ---
         _GET_(self%id_RDc, RDc)
         RDc = max(0.0_rk, RDc)
         fOx_c = k_eff * RDc

         ! DOM loss
         _SET_ODE_(self%id_RDc, -fOx_c)

         ! O2 consumption
         _SET_ODE_(self%id_O2o, -fOx_c * self%ur_O2)

         ! DIC production
         !---> jsasaki 2026-11-24: isw_law = 1: ersem_shared's CMass (12.011); 0: the former local 12.0
         if (self%isw_law == 1) then
            _SET_ODE_(self%id_O3c, +fOx_c / CMass)
         else
            _SET_ODE_(self%id_O3c, +fOx_c / CMass_legacy)
         end if
         !<--- jsasaki 2026-11-24

         ! --- Nitrogen (if tracked by this DOM class) ---
         if (_AVAILABLE_(self%id_RDn)) then
            _GET_(self%id_RDn, RDn)
            fOx_n = k_eff * max(0.0_rk, RDn)
            _SET_ODE_(self%id_RDn, -fOx_n)
            if (_AVAILABLE_(self%id_N4n)) _SET_ODE_(self%id_N4n, +fOx_n)
         else
            fOx_n = 0.0_rk
         end if

         ! --- Phosphorus (if tracked) ---
         if (_AVAILABLE_(self%id_RDp)) then
            _GET_(self%id_RDp, RDp)
            fOx_p = k_eff * max(0.0_rk, RDp)
            _SET_ODE_(self%id_RDp, -fOx_p)
            if (_AVAILABLE_(self%id_N1p)) _SET_ODE_(self%id_N1p, +fOx_p)
         else
            fOx_p = 0.0_rk
         end if

         ! --- Total Alkalinity: +1/NH4, -1/PO4 ---
         _SET_ODE_(self%id_TA, +fOx_n - fOx_p)

         ! --- Diagnostic ---
         _SET_DIAGNOSTIC_(self%id_rate, fOx_c)

      _LOOP_END_

   end subroutine do

end module ersem_dom_decay
