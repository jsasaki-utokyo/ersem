#include "fabm_driver.h"

! =========================================================================
! pom_decay.F90 — POM hydrolysis and direct oxidation module for ERSEM
!
! Parameterizes particle-attached microbial decomposition of POM in the
! water column via two pathways:
!
!   1. Autolysis (POM → DOM): Extracellular enzymatic hydrolysis by
!      particle-attached bacteria. Products are distributed to R1/R2/R3
!      following DiToro G-model fractions (G1=65%, G2=25%, G3=10%).
!      No O2 consumed — O2 is consumed when DOM is later oxidized.
!
!   2. Direct oxidation (POM → CO2 + O2): Respiration by particle-
!      attached bacteria at the particle surface. O2 consumed directly.
!
! The split between pathways is controlled by f_aut (~0.9, BROM-like).
! B1's sRP pathway should be set to zero when using this module.
!
! References:
!   BROM (Yakushev et al. 2017 GMD): autolysis=0.15/d, oxidation=0.01/d
!   DiToro (2001): G1=65%, G2=25%, G3=10% reactivity classes
!   CE-QUAL-ICM (Cerco & Cole 1993): LPOC=0.075/d, RPOC=0.02/d
!
! jsasaki 2026-03-19: Initial implementation (direct oxidation only)
! jsasaki 2026-03-20: Redesign with autolysis pathway and G-model DOM split
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
      type (type_state_variable_id) :: id_R1c   ! labile DOM carbon (G1)
      type (type_state_variable_id) :: id_R2c   ! semi-labile DOM carbon (G2)
      type (type_state_variable_id) :: id_R3c   ! semi-refractory DOM carbon (G3)
      type (type_state_variable_id) :: id_R1n, id_R1p  ! R1 nutrients (if available)

      ! Inorganic pools (for direct oxidation pathway)
      type (type_state_variable_id) :: id_O2o   ! dissolved oxygen
      type (type_state_variable_id) :: id_O3c   ! dissolved inorganic carbon
      type (type_state_variable_id) :: id_N4n   ! ammonium
      type (type_state_variable_id) :: id_N1p   ! phosphate
      type (type_state_variable_id) :: id_N5s   ! silicate
      type (type_state_variable_id) :: id_N7f   ! dissolved iron
      type (type_state_variable_id) :: id_TA    ! total alkalinity

      ! Environmental dependencies
      type (type_dependency_id) :: id_ETW       ! temperature

      ! Diagnostics
      type (type_diagnostic_variable_id) :: id_rate       ! total decomposition rate
      type (type_diagnostic_variable_id) :: id_rate_aut   ! autolysis rate
      type (type_diagnostic_variable_id) :: id_rate_ox    ! direct oxidation rate

      ! Parameters
      real(rk) :: k_decomp    ! total decomposition rate at Tref (1/d)
      real(rk) :: f_aut       ! fraction routed to DOM (autolysis)
      real(rk) :: f_G1        ! G1 fraction of autolysis products (→R1)
      real(rk) :: f_G2        ! G2 fraction (→R2)
      real(rk) :: f_G3        ! G3 fraction (→R3)
      real(rk) :: q10         ! Q10 temperature coefficient
      real(rk) :: Tref        ! reference temperature (deg C)
      real(rk) :: K_O2        ! O2 half-saturation (mmol O2/m3)
      real(rk) :: ur_O2       ! O2 per C for direct oxidation (mmol O2/mg C)
   contains
      procedure :: initialize
      procedure :: do
   end type

   real(rk), parameter :: CMass = 12.0_rk  ! mg C per mmol C

contains

   subroutine initialize(self, configunit)
      class(type_ersem_pom_decay), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! --- Parameters ---
      call self%get_parameter(self%k_decomp, 'k_decomp', '1/d', &
         'total decomposition rate at reference temperature', default=0.03_rk)
      call self%get_parameter(self%f_aut, 'f_aut', '-', &
         'fraction of decomposition routed to DOM (autolysis)', default=0.9_rk)
      call self%get_parameter(self%f_G1, 'f_G1', '-', &
         'labile fraction of autolysis products (to R1)', default=0.65_rk)
      call self%get_parameter(self%f_G2, 'f_G2', '-', &
         'semi-labile fraction of autolysis products (to R2)', default=0.25_rk)
      call self%get_parameter(self%f_G3, 'f_G3', '-', &
         'semi-refractory fraction of autolysis products (to R3)', default=0.10_rk)
      call self%get_parameter(self%q10, 'q10', '-', &
         'Q10 temperature coefficient', default=2.0_rk)
      call self%get_parameter(self%Tref, 'Tref', 'degrees_Celsius', &
         'reference temperature', default=20.0_rk)
      call self%get_parameter(self%K_O2, 'K_O2', 'mmol O_2/m^3', &
         'half-saturation oxygen concentration', default=5.0_rk)
      call self%get_parameter(self%ur_O2, 'ur_O2', 'mmol O_2/mg C', &
         'oxygen consumed per carbon in direct oxidation', default=0.1_rk)

      ! --- Target POM (carbon required, others optional) ---
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

      ! --- DOM targets for autolysis (all required when f_aut > 0) ---
      call self%register_state_dependency(self%id_R1c, 'R1c', 'mg C/m^3', &
         'labile dissolved organic carbon')
      call self%register_state_dependency(self%id_R2c, 'R2c', 'mg C/m^3', &
         'semi-labile dissolved organic carbon')
      call self%register_state_dependency(self%id_R3c, 'R3c', 'mg C/m^3', &
         'semi-refractory dissolved organic carbon')
      ! R1 nutrients (R2/R3 are C-only in ERSEM)
      call self%register_state_dependency(self%id_R1n, 'R1n', 'mmol N/m^3', &
         'labile dissolved organic nitrogen', required=.false.)
      call self%register_state_dependency(self%id_R1p, 'R1p', 'mmol P/m^3', &
         'labile dissolved organic phosphorus', required=.false.)

      ! --- Inorganic pools (for direct oxidation + nutrient release) ---
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

      real(rk) :: ETW, O2o
      real(rk) :: RPc, RPn, RPp, RPs, RPf
      real(rk) :: et, fO2, k_eff
      real(rk) :: fTotal_c         ! total decomposition (mg C/m3/d)
      real(rk) :: fAut_c           ! autolysis fraction (mg C/m3/d)
      real(rk) :: fOx_c            ! direct oxidation fraction (mg C/m3/d)
      real(rk) :: fDecomp_n, fDecomp_p, fDecomp_s, fDecomp_f
      real(rk) :: fAut_n, fAut_p   ! autolysis nutrient fluxes
      real(rk) :: fOx_n, fOx_p     ! direct oxidation nutrient fluxes

      _LOOP_BEGIN_

         ! --- Temperature dependence (ERSEM standard Q10 with high-T cutoff) ---
         _GET_(self%id_ETW, ETW)
         et = max(0.0_rk, self%q10**((ETW - self%Tref) / 10.0_rk) &
              - self%q10**((ETW - 32.0_rk) / 3.0_rk))

         ! --- Oxygen limitation (Michaelis-Menten) ---
         _GET_(self%id_O2o, O2o)
         O2o = max(0.0_rk, O2o)
         fO2 = O2o / (O2o + self%K_O2)

         ! --- Effective rate and total decomposition ---
         k_eff = self%k_decomp * et * fO2

         _GET_(self%id_RPc, RPc)
         RPc = max(0.0_rk, RPc)

         fTotal_c = k_eff * RPc
         fAut_c   = fTotal_c * self%f_aut        ! → DOM
         fOx_c    = fTotal_c * (1.0_rk - self%f_aut)  ! → CO2 + O2

         ! ================================================================
         ! POM carbon loss (both pathways)
         ! ================================================================
         _SET_ODE_(self%id_RPc, -fTotal_c)

         ! ================================================================
         ! Pathway 1: AUTOLYSIS (POM → DOM)
         ! Carbon distributed to R1/R2/R3 by G-model fractions.
         ! No O2 consumed. O2 is consumed later by B1 or dom_decay.
         ! ================================================================
         _SET_ODE_(self%id_R1c, +fAut_c * self%f_G1)
         _SET_ODE_(self%id_R2c, +fAut_c * self%f_G2)
         _SET_ODE_(self%id_R3c, +fAut_c * self%f_G3)

         ! DIC from autolysis: none (carbon goes to DOM, not CO2)

         ! ================================================================
         ! Pathway 2: DIRECT OXIDATION (POM → CO2 + O2)
         ! Particle-surface respiration by attached bacteria.
         ! ================================================================
         _SET_ODE_(self%id_O2o, -fOx_c * self%ur_O2)
         _SET_ODE_(self%id_O3c, +fOx_c / CMass)

         ! ================================================================
         ! Nutrient release (both pathways release nutrients)
         ! Autolysis N → R1 (labile DOM N) if available, else → NH4
         ! Direct oxidation N/P → inorganic (NH4, PO4)
         ! ================================================================

         ! --- Nitrogen ---
         if (_AVAILABLE_(self%id_RPn)) then
            _GET_(self%id_RPn, RPn)
            fDecomp_n = k_eff * max(0.0_rk, RPn)
            fAut_n = fDecomp_n * self%f_aut
            fOx_n  = fDecomp_n * (1.0_rk - self%f_aut)

            _SET_ODE_(self%id_RPn, -fDecomp_n)

            ! Autolysis N: to R1 (labile DOM has N tracking) if coupled
            if (_AVAILABLE_(self%id_R1n)) then
               _SET_ODE_(self%id_R1n, +fAut_n)
            else
               ! Fallback: release as NH4
               if (_AVAILABLE_(self%id_N4n)) _SET_ODE_(self%id_N4n, +fAut_n)
            end if

            ! Direct oxidation N: to NH4
            if (_AVAILABLE_(self%id_N4n)) _SET_ODE_(self%id_N4n, +fOx_n)
         else
            fDecomp_n = 0.0_rk
            fAut_n = 0.0_rk
            fOx_n = 0.0_rk
         end if

         ! --- Phosphorus ---
         if (_AVAILABLE_(self%id_RPp)) then
            _GET_(self%id_RPp, RPp)
            fDecomp_p = k_eff * max(0.0_rk, RPp)
            fAut_p = fDecomp_p * self%f_aut
            fOx_p  = fDecomp_p * (1.0_rk - self%f_aut)

            _SET_ODE_(self%id_RPp, -fDecomp_p)

            ! Autolysis P: to R1 (labile DOM has P tracking) if coupled
            if (_AVAILABLE_(self%id_R1p)) then
               _SET_ODE_(self%id_R1p, +fAut_p)
            else
               if (_AVAILABLE_(self%id_N1p)) _SET_ODE_(self%id_N1p, +fAut_p)
            end if

            ! Direct oxidation P: to PO4
            if (_AVAILABLE_(self%id_N1p)) _SET_ODE_(self%id_N1p, +fOx_p)
         else
            fDecomp_p = 0.0_rk
            fAut_p = 0.0_rk
            fOx_p = 0.0_rk
         end if

         ! --- Silicate → dissolved Si (both pathways) ---
         if (_AVAILABLE_(self%id_RPs)) then
            _GET_(self%id_RPs, RPs)
            fDecomp_s = k_eff * max(0.0_rk, RPs)
            _SET_ODE_(self%id_RPs, -fDecomp_s)
            if (_AVAILABLE_(self%id_N5s)) _SET_ODE_(self%id_N5s, +fDecomp_s)
         end if

         ! --- Iron → dissolved Fe (both pathways) ---
         if (_AVAILABLE_(self%id_RPf)) then
            _GET_(self%id_RPf, RPf)
            fDecomp_f = k_eff * max(0.0_rk, RPf)
            _SET_ODE_(self%id_RPf, -fDecomp_f)
            if (_AVAILABLE_(self%id_N7f)) _SET_ODE_(self%id_N7f, +fDecomp_f)
         end if

         ! ================================================================
         ! Total Alkalinity
         ! Autolysis: N to R1 (no TA effect if DOM-N tracked),
         !   or N to NH4 (+1) if fallback. P same logic.
         ! Direct oxidation: NH4 (+1), PO4 (-1)
         ! Simplified: TA += total_N_to_NH4 - total_P_to_PO4
         ! When R1n/R1p coupled: only direct oxidation contributes
         ! When not coupled: both pathways contribute
         ! ================================================================
         if (_AVAILABLE_(self%id_R1n)) then
            ! Autolysis N goes to R1n (no TA effect), only direct ox to NH4
            _SET_ODE_(self%id_TA, +fOx_n - fOx_p)
         else
            ! All N goes to NH4
            _SET_ODE_(self%id_TA, +fDecomp_n - fDecomp_p)
         end if

         ! --- Diagnostics ---
         _SET_DIAGNOSTIC_(self%id_rate, fTotal_c)
         _SET_DIAGNOSTIC_(self%id_rate_aut, fAut_c)
         _SET_DIAGNOSTIC_(self%id_rate_ox, fOx_c)

      _LOOP_END_

   end subroutine do

end module ersem_pom_decay
