#include "fabm_driver.h"

! =========================================================================
! pom_decay.F90 — POM hydrolysis and direct oxidation module for ERSEM
!
! Two independent pathways with separate rate constants:
!
!   1. Autolysis (k_autolysis): POM → DOM (R1/R2/R3)
!      Extracellular enzymatic hydrolysis by particle-attached bacteria.
!      Products distributed by DiToro G-model: G1(65%)→R1, G2(25%)→R2, G3(10%)→R3.
!      No O2 consumed. O2 is consumed later when DOM is oxidized.
!
!   2. Direct oxidation (k_oxidation): POM → CO2 + O2 consumed
!      Respiration by particle-attached bacteria at particle surface.
!
! The two rates are INDEPENDENT (not a fraction of a single rate):
!   BROM reference: k_autolysis=0.15/d, k_oxidation=0.01/d
!
! B1's sRP pathway should be set to zero when using this module.
!
! References:
!   BROM (Yakushev et al. 2017): autolysis=0.15/d, oxidation=0.01/d
!   DiToro (2001): G1=65%, G2=25%, G3=10%
!
! jsasaki 2026-03-19: Initial implementation (direct oxidation only)
! jsasaki 2026-03-20: Redesign with autolysis + G-model DOM split
! jsasaki 2026-03-21: Separate independent rate constants for autolysis
!                     and direct oxidation (no longer coupled via f_aut)
! =========================================================================

module ersem_pom_decay

   use fabm_types
   use ersem_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ersem_pom_decay
      ! Target POM state variables
      type (type_state_variable_id) :: id_RPc, id_RPn, id_RPp, id_RPs, id_RPf

      ! DOM targets for autolysis products
      type (type_state_variable_id) :: id_R1c, id_R2c, id_R3c
      type (type_state_variable_id) :: id_R1n, id_R1p

      ! Inorganic pools (for direct oxidation pathway)
      type (type_state_variable_id) :: id_O2o, id_O3c
      type (type_state_variable_id) :: id_N4n, id_N1p, id_N5s, id_N7f
      type (type_state_variable_id) :: id_TA

      ! Environmental dependencies
      type (type_dependency_id) :: id_ETW

      ! Diagnostics
      type (type_diagnostic_variable_id) :: id_rate
      type (type_diagnostic_variable_id) :: id_rate_aut
      type (type_diagnostic_variable_id) :: id_rate_ox

      ! Parameters — two independent rate constants
      real(rk) :: k_autolysis  ! autolysis rate at Tref (1/d): POM → DOM
      real(rk) :: k_oxidation  ! direct oxidation rate at Tref (1/d): POM → CO2 + O2
      real(rk) :: f_G1, f_G2, f_G3  ! G-model fractions for autolysis products
      real(rk) :: q10, Tref, K_O2, ur_O2
   contains
      procedure :: initialize
      procedure :: do
   end type

   real(rk), parameter :: CMass = 12.0_rk

contains

   subroutine initialize(self, configunit)
      class(type_ersem_pom_decay), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! --- Two independent rate constants (jsasaki 2026-03-21) ---
      call self%get_parameter(self%k_autolysis, 'k_autolysis', '1/d', &
         'autolysis rate: POM to DOM (BROM: 0.15)', default=0.15_rk)
      call self%get_parameter(self%k_oxidation, 'k_oxidation', '1/d', &
         'direct oxidation rate: POM to CO2+O2 (BROM: 0.01)', default=0.01_rk)

      ! --- G-model fractions (DiToro 2001) ---
      call self%get_parameter(self%f_G1, 'f_G1', '-', &
         'labile fraction of autolysis products (to R1)', default=0.65_rk)
      call self%get_parameter(self%f_G2, 'f_G2', '-', &
         'semi-labile fraction (to R2)', default=0.25_rk)
      call self%get_parameter(self%f_G3, 'f_G3', '-', &
         'semi-refractory fraction (to R3)', default=0.10_rk)

      ! --- Common parameters ---
      call self%get_parameter(self%q10, 'q10', '-', &
         'Q10 temperature coefficient', default=2.0_rk)
      call self%get_parameter(self%Tref, 'Tref', 'degrees_Celsius', &
         'reference temperature', default=20.0_rk)
      call self%get_parameter(self%K_O2, 'K_O2', 'mmol O_2/m^3', &
         'half-saturation oxygen concentration', default=5.0_rk)
      call self%get_parameter(self%ur_O2, 'ur_O2', 'mmol O_2/mg C', &
         'oxygen consumed per carbon in direct oxidation', default=0.1_rk)

      ! --- Target POM ---
      call self%register_state_dependency(self%id_RPc, 'RPc', 'mg C/m^3', &
         'particulate organic carbon')
      call self%register_state_dependency(self%id_RPn, 'RPn', 'mmol N/m^3', &
         'particulate organic nitrogen', required=.false.)
      call self%register_state_dependency(self%id_RPp, 'RPp', 'mmol P/m^3', &
         'particulate organic phosphorus', required=.false.)
      call self%register_state_dependency(self%id_RPs, 'RPs', 'mmol Si/m^3', &
         'particulate organic silicate', required=.false.)
      call self%register_state_dependency(self%id_RPf, 'RPf', 'umol Fe/m^3', &
         'particulate organic iron', required=.false.)

      ! --- DOM targets for autolysis ---
      call self%register_state_dependency(self%id_R1c, 'R1c', 'mg C/m^3', &
         'labile dissolved organic carbon')
      call self%register_state_dependency(self%id_R2c, 'R2c', 'mg C/m^3', &
         'semi-labile dissolved organic carbon')
      call self%register_state_dependency(self%id_R3c, 'R3c', 'mg C/m^3', &
         'semi-refractory dissolved organic carbon')
      call self%register_state_dependency(self%id_R1n, 'R1n', 'mmol N/m^3', &
         'labile dissolved organic nitrogen', required=.false.)
      call self%register_state_dependency(self%id_R1p, 'R1p', 'mmol P/m^3', &
         'labile dissolved organic phosphorus', required=.false.)

      ! --- Inorganic pools ---
      call self%register_state_dependency(self%id_O2o, 'O2o', 'mmol O_2/m^3', &
         'dissolved oxygen')
      call self%register_state_dependency(self%id_O3c, 'O3c', 'mmol C/m^3', &
         'dissolved inorganic carbon')
      call self%register_state_dependency(self%id_N4n, 'N4n', 'mmol N/m^3', &
         'ammonium', required=.false.)
      call self%register_state_dependency(self%id_N1p, 'N1p', 'mmol P/m^3', &
         'phosphate', required=.false.)
      call self%register_state_dependency(self%id_N5s, 'N5s', 'mmol Si/m^3', &
         'silicate', required=.false.)
      call self%register_state_dependency(self%id_N7f, 'N7f', 'umol Fe/m^3', &
         'dissolved iron', required=.false.)
      call self%register_state_dependency(self%id_TA, &
         standard_variables%alkalinity_expressed_as_mole_equivalent)

      ! --- Environment ---
      call self%register_dependency(self%id_ETW, standard_variables%temperature)

      ! --- Diagnostics ---
      call self%register_diagnostic_variable(self%id_rate, 'rate', 'mg C/m^3/d', &
         'total POM decomposition rate', output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_rate_aut, 'rate_aut', 'mg C/m^3/d', &
         'autolysis rate (POM to DOM)', output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_rate_ox, 'rate_ox', 'mg C/m^3/d', &
         'direct oxidation rate (POM to CO2)', output=output_time_step_averaged)

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_pom_decay), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ETW, O2o, RPc, RPn, RPp, RPs, RPf
      real(rk) :: et, fO2
      real(rk) :: k_aut_eff, k_ox_eff   ! effective rates after T and O2 correction
      real(rk) :: fAut_c, fOx_c, fTotal_c
      real(rk) :: fAut_n, fAut_p, fOx_n, fOx_p
      real(rk) :: fDecomp_n, fDecomp_p

      _LOOP_BEGIN_

         ! --- Temperature dependence (ERSEM standard Q10) ---
         _GET_(self%id_ETW, ETW)
         et = max(0.0_rk, self%q10**((ETW - self%Tref) / 10.0_rk) &
              - self%q10**((ETW - 32.0_rk) / 3.0_rk))

         ! --- Oxygen limitation (Michaelis-Menten) ---
         _GET_(self%id_O2o, O2o)
         O2o = max(0.0_rk, O2o)
         fO2 = O2o / (O2o + self%K_O2)

         ! --- Independent effective rates ---
         k_aut_eff = self%k_autolysis * et * fO2
         k_ox_eff  = self%k_oxidation * et * fO2

         ! --- POM carbon ---
         _GET_(self%id_RPc, RPc)
         RPc = max(0.0_rk, RPc)

         ! --- Decomposition rates (independent) ---
         fAut_c   = k_aut_eff * RPc       ! POM → DOM
         fOx_c    = k_ox_eff  * RPc       ! POM → CO2 + O2
         fTotal_c = fAut_c + fOx_c        ! total POM loss

         ! ================================================================
         ! POM carbon loss (both pathways)
         ! ================================================================
         _SET_ODE_(self%id_RPc, -fTotal_c)

         ! ================================================================
         ! Pathway 1: AUTOLYSIS (POM → DOM)
         ! G-model distribution to R1/R2/R3. No O2 consumed.
         ! ================================================================
         _SET_ODE_(self%id_R1c, +fAut_c * self%f_G1)
         _SET_ODE_(self%id_R2c, +fAut_c * self%f_G2)
         _SET_ODE_(self%id_R3c, +fAut_c * self%f_G3)

         ! ================================================================
         ! Pathway 2: DIRECT OXIDATION (POM → CO2 + O2)
         ! ================================================================
         _SET_ODE_(self%id_O2o, -fOx_c * self%ur_O2)
         _SET_ODE_(self%id_O3c, +fOx_c / CMass)

         ! ================================================================
         ! Nutrient release
         ! ================================================================
         ! --- Nitrogen ---
         if (_AVAILABLE_(self%id_RPn)) then
            _GET_(self%id_RPn, RPn)
            RPn = max(0.0_rk, RPn)
            fAut_n = k_aut_eff * RPn
            fOx_n  = k_ox_eff  * RPn
            fDecomp_n = fAut_n + fOx_n
            _SET_ODE_(self%id_RPn, -fDecomp_n)
            ! Autolysis N → R1 (DOM-N) if coupled, else → NH4
            if (_AVAILABLE_(self%id_R1n)) then
               _SET_ODE_(self%id_R1n, +fAut_n)
            else
               if (_AVAILABLE_(self%id_N4n)) _SET_ODE_(self%id_N4n, +fAut_n)
            end if
            ! Direct oxidation N → NH4
            if (_AVAILABLE_(self%id_N4n)) _SET_ODE_(self%id_N4n, +fOx_n)
         else
            fDecomp_n = 0.0_rk; fAut_n = 0.0_rk; fOx_n = 0.0_rk
         end if

         ! --- Phosphorus ---
         if (_AVAILABLE_(self%id_RPp)) then
            _GET_(self%id_RPp, RPp)
            RPp = max(0.0_rk, RPp)
            fAut_p = k_aut_eff * RPp
            fOx_p  = k_ox_eff  * RPp
            _SET_ODE_(self%id_RPp, -(fAut_p + fOx_p))
            if (_AVAILABLE_(self%id_R1p)) then
               _SET_ODE_(self%id_R1p, +fAut_p)
            else
               if (_AVAILABLE_(self%id_N1p)) _SET_ODE_(self%id_N1p, +fAut_p)
            end if
            if (_AVAILABLE_(self%id_N1p)) _SET_ODE_(self%id_N1p, +fOx_p)
         else
            fAut_p = 0.0_rk; fOx_p = 0.0_rk
         end if

         ! --- Silicate ---
         if (_AVAILABLE_(self%id_RPs)) then
            _GET_(self%id_RPs, RPs)
            _SET_ODE_(self%id_RPs, -(k_aut_eff + k_ox_eff) * max(0.0_rk, RPs))
            if (_AVAILABLE_(self%id_N5s)) &
               _SET_ODE_(self%id_N5s, +(k_aut_eff + k_ox_eff) * max(0.0_rk, RPs))
         end if

         ! --- Iron ---
         if (_AVAILABLE_(self%id_RPf)) then
            _GET_(self%id_RPf, RPf)
            _SET_ODE_(self%id_RPf, -(k_aut_eff + k_ox_eff) * max(0.0_rk, RPf))
            if (_AVAILABLE_(self%id_N7f)) &
               _SET_ODE_(self%id_N7f, +(k_aut_eff + k_ox_eff) * max(0.0_rk, RPf))
         end if

         ! ================================================================
         ! Total Alkalinity: +1/NH4, -1/PO4
         ! ================================================================
         if (_AVAILABLE_(self%id_R1n)) then
            _SET_ODE_(self%id_TA, +fOx_n - fOx_p)
         else
            _SET_ODE_(self%id_TA, +fDecomp_n - (fAut_p + fOx_p))
         end if

         ! --- Diagnostics ---
         _SET_DIAGNOSTIC_(self%id_rate, fTotal_c)
         _SET_DIAGNOSTIC_(self%id_rate_aut, fAut_c)
         _SET_DIAGNOSTIC_(self%id_rate_ox, fOx_c)

      _LOOP_END_

   end subroutine do

end module ersem_pom_decay
