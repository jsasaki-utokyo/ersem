#include "fabm_driver.h"

! -----------------------------------------------------------------------------
! Simple nutrient relaxation module for GOTM-FABM-ERSEM.
!
! Nudges a pelagic state variable toward a prescribed target concentration
! with a configurable timescale (tau, in days). This mimics lateral nutrient
! supply from rivers and ocean exchange in a 1D model that lacks advection.
!
! Usage in fabm.yaml:
!   relax_N3:
!     model: ersem/nutrient_relaxation
!     parameters:
!       target: 20.0    # target concentration (same units as variable)
!       tau: 20.0       # relaxation timescale (days)
!     coupling:
!       variable: N3/n  # state variable to relax
! -----------------------------------------------------------------------------

module ersem_nutrient_relaxation

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_ersem_nutrient_relaxation
      type(type_state_variable_id) :: id_variable
      real(rk) :: target_value
      real(rk) :: tau
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class(type_ersem_nutrient_relaxation), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! Set time unit to days (consistent with ERSEM convention)
      self%dt = 86400.0_rk

      call self%get_parameter(self%target_value, 'target', '', &
         'target concentration for relaxation')
      call self%get_parameter(self%tau, 'tau', 'd', &
         'relaxation timescale', default=20.0_rk, minimum=0.1_rk)

      call self%register_state_dependency(self%id_variable, &
         'variable', '', 'state variable to relax toward target')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_nutrient_relaxation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value

      _LOOP_BEGIN_
         _GET_(self%id_variable, value)
         ! Nudging: d(variable)/dt += (target - value) / tau
         ! Rate in [variable_units / day] since self%dt = 86400
         _ADD_SOURCE_(self%id_variable, (self%target_value - value) / self%tau)
      _LOOP_END_

   end subroutine do

end module ersem_nutrient_relaxation
