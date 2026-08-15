#include "fabm_driver.h"

! Benthic variable that supports resuspension and remineralization.
! Both processes return material to the pelagic.

module ersem_benthic_calcite

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_calcite
      type (type_horizontal_diagnostic_variable_id) :: id_dissolution
      type (type_horizontal_diagnostic_variable_id) :: id_precipitation
      type (type_dependency_id)                     :: id_Om_Cal
      type (type_dependency_id)                     :: id_par

      ! Parameters
      real(rk) :: fdissmax, fdissmin, ndiss, KcalomX
      real(rk) :: K_prec, K_par_prec, n_prec
      integer  :: iswcal
      
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_ersem_benthic_calcite), intent(inout), target :: self
   integer,                         intent(in)            :: configunit

   real(rk) :: c0
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%initialize_ersem_benthic_base()

      call self%get_parameter(self%iswcal,'iswcal','','dissolution dependence on calcite saturation (0: none, 1: power law, 2: hyperbolic)', minimum=0, maximum=2)
      select case (self%iswcal)
         case (1)
            call self%get_parameter(self%ndiss,'ndiss','-','power of the dissolution law (Keir 1980)')
         case (2)
            call self%get_parameter(self%KcalomX,'KcalomX','-','half-saturation constant for calcification limitation from saturation state')
      end select
      if (self%iswcal == 0) then
         call self%get_parameter(self%fdissmin, 'fdiss', '1/d','specific dissolution rate', default=0.0_rk)
         self%fdissmax = 0.0_rk
      else
         call self%get_parameter(self%fdissmax, 'fdissmax', '1/d','maximum specific dissolution rate', minimum=0._rk, default=0.0_rk)
         call self%get_parameter(self%fdissmin, 'fdissmin', '1/d','minimum specific dissolution rate', minimum=0._rk, default=0.001_rk * self%fdissmax)
         call self%register_dependency(self%id_Om_Cal,'Om_Cal','-','calcite saturation')
      end if
      ! Light-driven mat calcification (jsasaki 2026-08-15; design:
      ! nippon-steel/docs/15-mat-calcification.md). Benthic photosynthesis
      ! elevates the mat microenvironment's pH/saturation in the light and
      ! precipitates CaCO3 (TA -2, DIC -1 per C); the existing dissolution
      ! law provides the dark return path. The light gate is the
      ! mat-photosynthesis proxy (same convention as the sulfur module's
      ! K_par_ox); water-column Om_cal stands in for the porewater state.
      ! Default K_prec = 0 keeps the module bit-identical.
      call self%get_parameter(self%K_prec, 'K_prec', 'mg C/m^2/d', &
         'maximum light-driven calcification rate (0: off)', &
         default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%K_par_prec, 'K_par_prec', 'W/m^2', &
         'PAR half-saturation of the calcification light gate', &
         default=4.0_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%n_prec, 'n_prec', '-', &
         'power of the supersaturation drive (Om_cal - 1)', &
         default=1.0_rk, minimum=0.0_rk)
      if (self%K_prec > 0.0_rk .and. self%iswcal == 0) &
         call self%register_dependency(self%id_Om_Cal,'Om_Cal','-','calcite saturation')
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_diagnostic_variable(self%id_precipitation,'precipitation','mg C/m^2/d','light-driven calcification',source=source_do_bottom)
      call self%get_parameter(c0,'c0','mg C/m^2','background calcite concentration',default=0.0_rk)

      call self%add_constituent('c',0.0_rk,c0)
      call self%register_state_dependency(self%id_O3c,'O3c','mmol/m^3','dissolved inorganic carbon')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_diagnostic_variable(self%id_dissolution,'dissolution','mg C/m^2/d','dissolution',source=source_do_bottom)

   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_calcite),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: bL2c
      real(rk) :: Om_Cal
      real(rk) :: fdiss
      real(rk) :: par, F_prec

      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_c, bL2c)
         if (self%iswcal>0 .or. self%K_prec>0._rk) then
            _GET_(self%id_om_cal, om_cal)
            om_cal=max(om_cal,0._rk)
         end if

         if (self%iswcal==0) then  ! NB select case would be cleaner but makes vectorization impossible for ifort 14
            fdiss = 0._rk
         elseif (self%iswcal==1) then
            fdiss = (max(1._rk-om_cal,0._rk))**self%ndiss
         else
            fdiss = max(0._rk,(1._rk-om_cal)/(1._rk-om_cal+self%KcalomX))
         end if

         fdiss = max(fdiss * self%fdissmax, self%fdissmin)

         _SET_BOTTOM_ODE_(self%id_c, -fdiss*bL2c)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dissolution, -fdiss*bL2c)
         _SET_BOTTOM_EXCHANGE_(self%id_O3c, fdiss*bL2c/CMass)
         _SET_BOTTOM_EXCHANGE_(self%id_TA, 2*fdiss*bL2c/CMass)  ! Dissolution of CaCO3 increases alkalinity by 2 units

         ! Light-driven mat calcification (docs/15; inert when K_prec = 0)
         F_prec = 0.0_rk
         if (self%K_prec > 0.0_rk) then
            _GET_(self%id_par, par)
            F_prec = self%K_prec &
                   * max(0.0_rk, par) / (max(0.0_rk, par) + self%K_par_prec) &
                   * (max(0.0_rk, om_cal - 1.0_rk))**self%n_prec
            _SET_BOTTOM_ODE_(self%id_c, F_prec)
            _SET_BOTTOM_EXCHANGE_(self%id_O3c, -F_prec/CMass)
            _SET_BOTTOM_EXCHANGE_(self%id_TA, -2*F_prec/CMass)  ! Precipitation of CaCO3 removes 2 units of alkalinity
         end if
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_precipitation, F_prec)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
