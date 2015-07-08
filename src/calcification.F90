#include "fabm_driver.h"
module ersem_calcification

   use fabm_types

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_calcification
!     Variable identifiers
      type (type_diagnostic_variable_id)   :: id_RainR
      type (type_state_variable_id)        :: id_O3c,id_TA
      type (type_bottom_state_variable_id) :: id_BL2c
      type (type_dependency_id)            :: id_om_cal

      integer  :: iswcal
      real(rk) :: Rain0
      real(rk) :: ncalc,ndiss,KcalomX
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
   end type
      
contains
      
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_ersem_calcification), intent(inout), target :: self
      integer,                              intent(in)            :: configunit
!
      real(rk) :: sedL2,c0
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%iswcal,'iswcal','','calcification/dissolution dependence on calcite saturation (1: power law, 2: hyperbolic)')
      if (self%iswcal<0.or.self%iswcal>2) then
         call self%log_message('ISWcal set to 1')
         self%iswcal = 1
      end if
      select case (self%iswcal)
         case (1)
            call self%get_parameter(self%ncalc,'ncalc','-','power of the calcification law (Ridgwell et al. 2007, mineral calcite)')
            call self%get_parameter(self%ndiss,'ndiss','-','power of the dissolution law (Keir 1980)')
         case (2)
            call self%get_parameter(self%KcalomX,'KcalomX','-','half-saturation constant for calcification limitation from saturation state')
      end select
      call self%get_parameter(self%Rain0,'Rain0','-','maximum rain ratio from PISCES')
      call self%get_parameter(sedL2,'sedL2','m/d','sinking velocity')
      call self%get_parameter(c0,'c0','mg C/m^3','background concentration',default=0.0_rk)

      call self%initialize_ersem_base(rm=sedL2,sedimentation=.false.)
      call self%get_parameter(self%sedimentation,'sedimentation','','enable sedimentation',default=.true.)
      call self%add_constituent('c',0.0_rk,c0)

      call self%register_diagnostic_variable(self%id_RainR,'RainR','1','rain ratio')

      call self%register_dependency(self%id_om_cal,'om_cal','-','calcite saturation')
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','total dissolved inorganic carbon')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)    

      if (self%sedimentation) then
          call self%register_bottom_state_dependency(self%id_bL2c,'bL2c','mg C/m^2','benthic calcite')
          call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
          call self%register_dependency(self%id_dens,     standard_variables%density)
      end if
   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_ersem_calcification), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: om_cal,L2c
      real(rk) :: fcalc,fdiss

      _LOOP_BEGIN_
         _GET_WITH_BACKGROUND_(self%id_c,L2c)  ! Jorn: legacy ersem includes background, but excluding background seems more appropriate
         _GET_(self%id_om_cal,om_cal)

         if (self%iswcal==0) then  ! NB select case would be cleaner but makes vectorization impossible for ifort 14
            fcalc = 0._rk
            fdiss = 0._rk
         elseif (self%iswcal==1) then
            fcalc = max(om_cal-1._rk,0._rk)**self%ncalc
            fdiss = max(1._rk-om_cal,0._rk)**self%ndiss
         else
            fcalc = max(0._rk,(om_cal-1._rk)/(om_cal-1._rk+self%KcalomX))
            fdiss = max(0._rk,(1._rk-om_cal)/(1._rk-om_cal+self%KcalomX))
         end if

         _SET_ODE_(self%id_c,  -fdiss*L2c)
         _SET_ODE_(self%id_O3c, fdiss*L2c/CMass)
         _SET_ODE_(self%id_TA,2*fdiss*L2c/CMass)  ! Dissolution of CaCO3 increases alkalinity by 2 units
         _SET_DIAGNOSTIC_(self%id_RainR,fcalc * self%Rain0)
      _LOOP_END_

   end subroutine do

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_calcification), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: tbed,density
      real(rk) :: fac,sdrate,fsd
      real(rk) :: L2c

      ! Bed characteristics - from Puls and Sundermann 1990
      ! Critical bed shear velocity = 0.01 m/s
      real(rk) :: tdep = 0.01_rk**2

      if (.not.self%sedimentation) return

      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_bedstress,tbed)
         _GET_(self%id_dens,density)

         _GET_(self%id_c,L2c)

         fsd = self%get_sinking_rate(_ARGUMENTS_LOCAL_)

         ! Divide actual stress (Pa) by density (kg/m^3) to obtain square of bed shear velocity.
         tbed = tbed/density

         if(tbed<tdep) then
            fac=1._rk-tbed/tdep
         else
            fac=0._rk
         endif

         !sdrate = min(fsd*fac,pdepth(I)/timestep) ! Jorn: CFL criterion disabled because FABM does not provide timestep
         sdrate = fsd*fac

         _SET_BOTTOM_ODE_(self%id_bL2c, sdrate*L2c)
         _SET_BOTTOM_EXCHANGE_(self%id_c,-sdrate*L2c)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module