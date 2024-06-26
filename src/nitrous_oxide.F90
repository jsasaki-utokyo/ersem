#include "fabm_driver.h"

!------------------------------------------------------------------------
! This module calculates saturation concentrations and air-sea
! exchange of nitrous oxide. Production of nitrous oxide is implemented
! in nitrification module (nitrification.F90).
! For implementation and validation of nitrous oxide on the North-West 
! European Shelf see Lessin et al. (2020): doi.org/10.1029/2019JG005613
!------------------------------------------------------------------------

module ersem_nitrous_oxide

   use fabm_types
   use fabm_particle

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_nitrous_oxide
      ! Variables
      type (type_dependency_id)                     :: id_ETW,id_X1X
      type (type_horizontal_dependency_id)          :: id_wnd,id_pN2Oa
      type (type_horizontal_diagnostic_variable_id) :: id_airsea
      type (type_diagnostic_variable_id)            :: id_satp

      integer :: iswN2O

   contains
!     Model procedures
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_ersem_nitrous_oxide),intent(inout),target :: self
      integer,                         intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%initialize_ersem_base(sedimentation=.false.)

      call self%add_constituent('n',0.0_rk)

      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_X1X,standard_variables%practical_salinity)
      call self%register_dependency(self%id_wnd,standard_variables%wind_speed)
      call self%register_dependency(self%id_pN2Oa,partial_pressure_of_n2o)
      call self%register_diagnostic_variable(self%id_satp,'satp','%','nitrous oxide % saturation')
      call self%register_diagnostic_variable(self%id_airsea,'airsea','mmol N2O/m^2/d','airsea flux of N2O',source=source_do_surface)
      call self%get_parameter(self%iswN2O,'iswN2O','','air-sea flux switch (0: off, 1: on)',default=1)

   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_ersem_nitrous_oxide), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: O5n,ETW,X1X
      real(rk) :: koN2O,pN2Oa

       _LOOP_BEGIN_
          _GET_(self%id_n,O5n)
          _GET_(self%id_ETW,ETW)
          _GET_(self%id_X1X,X1X)
          _GET_HORIZONTAL_(self%id_pN2Oa,pN2Oa)

         koN2O = n2o_transfer_coefficient(self,ETW,X1X)

       _SET_DIAGNOSTIC_(self%id_satp,(100._rk*O5n*0.5_rk/(KoN2O*pN2Oa*1.e-9*1.e6)))
      _LOOP_END_
   end subroutine

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_ersem_nitrous_oxide), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: O5n,ETW,X1X,wnd, fwind
      real(rk) :: koN2O,sc,n2oflux,pN2Oa

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_n,O5n)
         _GET_(self%id_ETW,ETW)
         _GET_(self%id_X1X,X1X)
         _GET_HORIZONTAL_(self%id_wnd,wnd)
         _GET_HORIZONTAL_(self%id_pN2Oa,pN2Oa)

! Schmidt number for N2O (Wanninkhof, 1992)
 
        sc=2301._rk-151._rk*ETW+4.74_rk*ETW**2._rk-0.06_rk*ETW**3
        fwind =  0.39_rk * wnd**2 *(sc/660._rk)**(-0.5_rk)
        fwind=fwind*24._rk/100._rk   ! convert to m/day

        koN2O = n2o_transfer_coefficient(self,ETW,X1X)

!        PNOatm is given in natm, so it is here converted to atm  by means of a 1.-9 factor. As KoNO is
!        given in mol L-1 atm -1 (Weiss and Price, 1980), a conversion factor of
!        1.e6 is used to have the final units of mmol m-3. N0 is given in mmol
!        of N so it is multiplied by 0.5 to have mmol of N2O. (Luca, July 2016)

        if (self%iswN2O .eq. 1) then
         n2oflux = fwind*(KoN2O*pN2Oa*1.e-9*1.e6-O5n*0.5_rk)
        else
         n2oflux=0._rk
        endif

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_airsea,n2oflux)

         _SET_SURFACE_EXCHANGE_(self%id_n,n2oflux*2._rk)

      _HORIZONTAL_LOOP_END_
   end subroutine do_surface

   function n2o_transfer_coefficient(self,ETW,X1X) result(koN2O)
       class (type_ersem_nitrous_oxide), intent(in) :: self
       real(rk),                         intent(in) :: ETW,X1X
       real(rk)                                     :: koN2O
       real(rk)           :: tk,tk100
       real(rk),parameter :: A1 = -62.7062_rk
       real(rk),parameter :: A2 = 97.3066_rk
       real(rk),parameter :: A3 = 24.1406_rk
       real(rk),parameter :: B1 = -0.05842_rk
       real(rk),parameter :: B2 = 0.033193_rk
       real(rk),parameter :: B3 = -0.0051_rk

       TK=ETW+273.15_rk
       TK100=TK/100._rk

        koN2O = exp(A1+A2/tk100 + A3 * log(tk100) + &
        &       X1X * (B1 + B2 * tk100 + B3 * tk100 ** 2._rk))
    end function
   
end module
