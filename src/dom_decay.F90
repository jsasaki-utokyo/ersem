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
   contains
      procedure :: initialize
      procedure :: do
   end type

   real(rk), parameter :: CMass = 12.0_rk  ! mg C per mmol C

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
         et = max(0.0_rk, self%q10**((ETW - self%Tref) / 10.0_rk) &
              - self%q10**((ETW - 32.0_rk) / 3.0_rk))

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
         _SET_ODE_(self%id_O3c, +fOx_c / CMass)

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
