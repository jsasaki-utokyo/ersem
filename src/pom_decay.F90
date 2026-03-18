#include "fabm_driver.h"

! =========================================================================
! pom_decay.F90 — Direct first-order POM decomposition module for ERSEM
!
! Parameterizes particle-attached microbial decomposition of POM in the
! water column. This process is independent of free-living bacteria (B1)
! and represents:
!   1. Extracellular enzymatic hydrolysis by particle-attached bacteria
!   2. Respiration by particle-attached microbial communities
!   3. Cell autolysis of dead phytoplankton
!
! The rate is first-order in POM concentration (proportional to particle
! surface area and rapid colonization of fresh particles).
!
! References:
!   BROM (Yakushev et al. 2017 GMD): K_PON_ox=0.01/d, autolysis=0.15/d
!   ERGOM (Neumann et al. 2022 GMD): detritus mineralization 0.003-0.02/d
!   CE-QUAL-ICM (Cerco & Cole 1993): LPOC=0.075/d, RPOC=0.02/d
!
! Each FABM instance targets one POM class (R4, R6, or R8).
! Multiple instances can be defined in fabm.yaml with different rates.
!
! jsasaki 2026-03-19: Initial implementation
! =========================================================================

module ersem_pom_decay

   use fabm_types
   use ersem_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ersem_pom_decay
      ! Target POM state variables (carbon + optional N, P, Si, Fe)
      type (type_state_variable_id) :: id_RPc, id_RPn, id_RPp, id_RPs, id_RPf

      ! Inorganic pools receiving decomposition products
      type (type_state_variable_id) :: id_O2o     ! dissolved oxygen
      type (type_state_variable_id) :: id_O3c     ! dissolved inorganic carbon
      type (type_state_variable_id) :: id_N4n     ! ammonium
      type (type_state_variable_id) :: id_N1p     ! phosphate
      type (type_state_variable_id) :: id_N5s     ! silicate
      type (type_state_variable_id) :: id_N7f     ! dissolved iron
      type (type_state_variable_id) :: id_TA      ! total alkalinity

      ! Environmental dependencies
      type (type_dependency_id) :: id_ETW         ! temperature

      ! Diagnostic
      type (type_diagnostic_variable_id) :: id_rate  ! decomposition rate

      ! Parameters
      real(rk) :: k_decomp    ! specific decomposition rate at Tref (1/d)
      real(rk) :: q10         ! Q10 temperature coefficient
      real(rk) :: Tref        ! reference temperature (deg C)
      real(rk) :: K_O2        ! O2 half-saturation for decomposition (mmol O2/m3)
      real(rk) :: ur_O2       ! O2 consumed per C respired (mmol O2/mg C)

      ! Flags for optional constituents
      logical :: has_n, has_p, has_s, has_f
   contains
      procedure :: initialize
      procedure :: do
   end type

   ! Molar mass of carbon (mg C/mmol C)
   real(rk), parameter :: CMass = 12.0_rk

contains

   subroutine initialize(self, configunit)
      class(type_ersem_pom_decay), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! Read parameters
      call self%get_parameter(self%k_decomp, 'k_decomp', '1/d', &
         'specific decomposition rate at reference temperature', default=0.02_rk)
      call self%get_parameter(self%q10, 'q10', '-', &
         'Q10 temperature coefficient', default=2.0_rk)
      call self%get_parameter(self%Tref, 'Tref', 'degrees_Celsius', &
         'reference temperature', default=20.0_rk)
      call self%get_parameter(self%K_O2, 'K_O2', 'mmol O_2/m^3', &
         'half-saturation oxygen concentration for decomposition', default=5.0_rk)
      call self%get_parameter(self%ur_O2, 'ur_O2', 'mmol O_2/mg C', &
         'oxygen consumed per carbon decomposed', default=0.1_rk)

      ! Register links to target POM (carbon is required)
      call self%register_state_dependency(self%id_RPc, 'RPc', 'mg C/m^3', &
         'particulate organic carbon')

      ! Optional POM constituents (required=.false. for composition-dependent)
      call self%register_state_dependency(self%id_RPn, 'RPn', 'mmol N/m^3', &
         'particulate organic nitrogen', required=.false.)
      call self%register_state_dependency(self%id_RPp, 'RPp', 'mmol P/m^3', &
         'particulate organic phosphorus', required=.false.)
      call self%register_state_dependency(self%id_RPs, 'RPs', 'mmol Si/m^3', &
         'particulate organic silicate', required=.false.)
      call self%register_state_dependency(self%id_RPf, 'RPf', 'umol Fe/m^3', &
         'particulate organic iron', required=.false.)

      ! Register links to inorganic pools
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

      ! Environmental dependencies
      call self%register_dependency(self%id_ETW, standard_variables%temperature)

      ! Diagnostic output
      call self%register_diagnostic_variable(self%id_rate, 'rate', 'mg C/m^3/d', &
         'POM decomposition rate', output=output_time_step_averaged)

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_pom_decay), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ETW, O2o
      real(rk) :: RPc, RPn, RPp, RPs, RPf
      real(rk) :: et, fO2, k_eff
      real(rk) :: fDecomp_c    ! decomposition rate (mg C/m3/d)
      real(rk) :: fDecomp_n    ! N release (mmol N/m3/d)
      real(rk) :: fDecomp_p    ! P release (mmol P/m3/d)
      real(rk) :: fDecomp_s    ! Si release (mmol Si/m3/d)
      real(rk) :: fDecomp_f    ! Fe release (umol Fe/m3/d)

      _LOOP_BEGIN_

         ! Get temperature
         _GET_(self%id_ETW, ETW)

         ! Temperature dependence (same Q10 formulation as ERSEM bacteria)
         et = max(0.0_rk, self%q10**((ETW - self%Tref) / 10.0_rk) &
              - self%q10**((ETW - 32.0_rk) / 3.0_rk))

         ! Oxygen limitation (Michaelis-Menten, like BROM/ERGOM)
         _GET_(self%id_O2o, O2o)
         O2o = max(0.0_rk, O2o)
         fO2 = O2o / (O2o + self%K_O2)

         ! Effective decomposition rate
         k_eff = self%k_decomp * et * fO2

         ! Get POM carbon
         _GET_(self%id_RPc, RPc)
         RPc = max(0.0_rk, RPc)

         ! Carbon decomposition rate (mg C/m3/d)
         fDecomp_c = k_eff * RPc

         ! ---- POM loss ----
         _SET_ODE_(self%id_RPc, -fDecomp_c)

         ! ---- O2 consumption (direct, like BROM/ERGOM/ICM) ----
         _SET_ODE_(self%id_O2o, -fDecomp_c * self%ur_O2)

         ! ---- DIC production ----
         _SET_ODE_(self%id_O3c, +fDecomp_c / CMass)

         ! ---- Nutrient release (proportional to POM stoichiometry) ----
         ! Nitrogen → NH4
         if (_AVAILABLE_(self%id_RPn)) then
            _GET_(self%id_RPn, RPn)
            if (RPc > 0.0_rk) then
               fDecomp_n = k_eff * max(0.0_rk, RPn)
            else
               fDecomp_n = 0.0_rk
            end if
            _SET_ODE_(self%id_RPn, -fDecomp_n)
            _SET_ODE_(self%id_N4n, +fDecomp_n)
         else
            fDecomp_n = 0.0_rk
         end if

         ! Phosphorus → PO4
         if (_AVAILABLE_(self%id_RPp)) then
            _GET_(self%id_RPp, RPp)
            if (RPc > 0.0_rk) then
               fDecomp_p = k_eff * max(0.0_rk, RPp)
            else
               fDecomp_p = 0.0_rk
            end if
            _SET_ODE_(self%id_RPp, -fDecomp_p)
            _SET_ODE_(self%id_N1p, +fDecomp_p)
         else
            fDecomp_p = 0.0_rk
         end if

         ! Silicate → dissolved Si
         if (_AVAILABLE_(self%id_RPs)) then
            _GET_(self%id_RPs, RPs)
            if (RPc > 0.0_rk) then
               fDecomp_s = k_eff * max(0.0_rk, RPs)
            else
               fDecomp_s = 0.0_rk
            end if
            _SET_ODE_(self%id_RPs, -fDecomp_s)
            _SET_ODE_(self%id_N5s, +fDecomp_s)
         end if

         ! Iron → dissolved Fe
         if (_AVAILABLE_(self%id_RPf)) then
            _GET_(self%id_RPf, RPf)
            if (RPc > 0.0_rk) then
               fDecomp_f = k_eff * max(0.0_rk, RPf)
            else
               fDecomp_f = 0.0_rk
            end if
            _SET_ODE_(self%id_RPf, -fDecomp_f)
            _SET_ODE_(self%id_N7f, +fDecomp_f)
         end if

         ! ---- Total Alkalinity ----
         ! Aerobic OM decomposition: NH4 release increases TA (+1/mmol N),
         ! PO4 release decreases TA (-1/mmol P).
         ! This follows the ERSEM convention used in bacteria_docdyn.F90
         ! (lines 418, 454) and mesozooplankton.F90 (line 487).
         _SET_ODE_(self%id_TA, +fDecomp_n - fDecomp_p)

         ! Diagnostic
         _SET_DIAGNOSTIC_(self%id_rate, fDecomp_c)

      _LOOP_END_

   end subroutine do

end module ersem_pom_decay
