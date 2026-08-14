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
      type(type_global_dependency_id) :: id_yday
      real(rk) :: target_value
      real(rk) :: tau
      logical  :: seasonal
      real(rk) :: target_month(12)
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class(type_ersem_nutrient_relaxation), intent(inout), target :: self
      integer, intent(in) :: configunit

      integer :: imonth
      character(len=10) :: pname

      ! Set time unit to days (consistent with ERSEM convention)
      self%dt = 86400.0_rk

      call self%get_parameter(self%target_value, 'target', '', &
         'target concentration for relaxation')
      call self%get_parameter(self%tau, 'tau', 'd', &
         'relaxation timescale', default=20.0_rk, minimum=0.1_rk)

      ! Optional seasonal mode (2026-08-14): with seasonal=true, twelve
      ! monthly targets (target_m01..target_m12, default = target) are
      ! linearly interpolated between month centres on the day of year.
      ! Default false preserves the original constant-target behaviour and
      ! registers nothing extra.
      call self%get_parameter(self%seasonal, 'seasonal', '', &
         'interpolate monthly targets target_m01..target_m12', default=.false.)
      if (self%seasonal) then
         do imonth = 1, 12
            write (pname, '(a,i2.2)') 'target_m', imonth
            call self%get_parameter(self%target_month(imonth), trim(pname), '', &
               'monthly target, month '//pname(9:10), default=self%target_value)
         end do
         call self%register_global_dependency(self%id_yday, &
            standard_variables%number_of_days_since_start_of_the_year)
      end if

      call self%register_state_dependency(self%id_variable, &
         'variable', '', 'state variable to relax toward target')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_nutrient_relaxation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value, target, yday, pos, w
      integer :: m1, m2

      target = self%target_value
      if (self%seasonal) then
         ! Piecewise-linear interpolation between month centres (day 15 of a
         ! 30.4-day nominal month), wrapping December -> January.
         _GET_GLOBAL_(self%id_yday, yday)
         pos = modulo(yday - 15.2_rk, 365.0_rk) / 30.4167_rk
         m1 = min(11, int(pos))
         w = pos - m1
         m2 = modulo(m1 + 1, 12) + 1
         m1 = m1 + 1
         target = (1.0_rk - w) * self%target_month(m1) + w * self%target_month(m2)
      end if

      _LOOP_BEGIN_
         _GET_(self%id_variable, value)
         ! Nudging: d(variable)/dt += (target - value) / tau
         ! Rate in [variable_units / day] since self%dt = 86400
         _ADD_SOURCE_(self%id_variable, (target - value) / self%tau)
      _LOOP_END_

   end subroutine do

end module ersem_nutrient_relaxation
