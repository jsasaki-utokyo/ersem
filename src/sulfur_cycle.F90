#include "fabm_driver.h"

!-----------------------------------------------------------------------
! Pelagic sulfur cycle module for ERSEM
!
! Handles H2S and S0 oxidation in the water column.
! Simplified implementation (Option B): H2S + S0 only, no explicit SO4.
!
! Reactions:
!   R1: H2S + 0.5 O2 -> S0 + H2O
!   R2: S0 + 1.5 O2 + H2O -> SO4 + 2H+ (SO4 not tracked, just O2 consumed)
!
! References:
!   - BROM model (Yakushev et al. 2017)
!   - Rate constants from fabm-niva-brom
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

      ! Set time unit to d-1 (ERSEM convention)
      self%dt = 86400._rk

      ! Get parameters (default values from BROM)
      call self%get_parameter(self%K_H2S_ox, 'K_H2S_ox', '1/d', &
           'H2S oxidation rate constant', default=0.5_rk)
      call self%get_parameter(self%K_S0_ox, 'K_S0_ox', '1/d', &
           'S0 oxidation rate constant', default=0.02_rk)
      call self%get_parameter(self%K_O2_half, 'K_O2_half', 'mmol/m^3', &
           'half-saturation O2 for oxidation', default=1.0_rk)

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

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_sulfur_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: H2S, S0, O2
      real(rk) :: R_H2S_ox, R_S0_ox, f_O2

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

         ! R1: H2S + 0.5 O2 -> S0 + H2O
         R_H2S_ox = self%K_H2S_ox * H2S * f_O2

         ! R2: S0 + 1.5 O2 + H2O -> SO4 + 2H+
         R_S0_ox = self%K_S0_ox * S0 * f_O2

         ! Set ODEs
         _SET_ODE_(self%id_H2S, -R_H2S_ox)
         _SET_ODE_(self%id_S0,   R_H2S_ox - R_S0_ox)
         _SET_ODE_(self%id_O2,  -0.5_rk * R_H2S_ox - 1.5_rk * R_S0_ox)

         ! Set diagnostics
         _SET_DIAGNOSTIC_(self%id_R_H2S_ox, R_H2S_ox)
         _SET_DIAGNOSTIC_(self%id_R_S0_ox, R_S0_ox)

      _LOOP_END_
   end subroutine do

end module ersem_sulfur_cycle
