#include "fabm_driver.h"
module ersem_calcification

   use fabm_types

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_ersem_calcification
!     Variable identifiers
      type (type_diagnostic_variable_id) :: id_RainR,id_L2O3c
      type (type_state_variable_id)      :: id_O3c,id_TA
      type (type_dependency_id)          :: id_om_cal

      integer  :: iswcal
      real(rk) :: Rain0,sL2O3X
      real(rk) :: ncalc,ndiss,KcalomX

      ! ledger counters (isw_ledger = 1 only; nippon-steel docs/119 FC1)
      integer :: isw_ledger
      type (type_diagnostic_variable_id) :: id_ledger_diss_c, id_ledger_diss_O3c, id_ledger_diss_TA
   contains
      procedure :: initialize
      procedure :: do
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
      call self%get_parameter(self%iswcal,'iswcal','','calcification/dissolution dependence on calcite saturation (1: power law, 2: hyperbolic)',minimum=1,maximum=2)
      select case (self%iswcal)
         case (1)
            call self%get_parameter(self%ncalc,'ncalc','-','power of the calcification law (Ridgwell et al. 2007, mineral calcite)')
            call self%get_parameter(self%ndiss,'ndiss','-','power of the dissolution law (Keir 1980)')
         case (2)
            call self%get_parameter(self%KcalomX,'KcalomX','-','half-saturation constant for calcification limitation from saturation state')
      end select
      call self%get_parameter(self%Rain0,'Rain0','-','maximum rain ratio from PISCES')
      call self%get_parameter(sedL2,'sedL2','m/d','sinking velocity')
      call self%get_parameter(self%sL2O3X,'sL2O3','1/d','maximum specific dissolution rate', default=1.0_rk)
      call self%get_parameter(c0,'c0','mg C/m^3','background concentration',default=0.0_rk)

      call self%initialize_ersem_base(rm=sedL2,sedimentation=.true.)
      call self%add_constituent('c',0.0_rk,c0)

      call self%register_diagnostic_variable(self%id_RainR,'RainR','1','rain ratio')
      call self%register_diagnostic_variable(self%id_L2O3c,'L2O3c','mg C/m^3/d','calcite dissolution rate')
      call self%register_dependency(self%id_om_cal,'om_cal','-','calcite saturation')
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','total dissolved inorganic carbon')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)

      ! LEDGER COUNTERS (2026-09-16, nippon-steel docs/119 §4.1 and §4.3 FC1).
      ! isw_ledger = 1 registers one diagnostic per source term this module
      ! writes to a ledger target (DIC, TA, O2, NO3, NH4, N2, PO4, Si, H2S, S0,
      ! CaCO3; organic matter only where it crosses between the water and the
      ! bed). Each value is the expression passed to the _SET_ODE_ /
      ! _SET_BOTTOM_ODE_ / _SET_BOTTOM_EXCHANGE_ call beside it, in this
      ! module's time unit (per day), with the sign applied to the target, so
      ! post-run gates can compare what the code applies against independently
      ! specified stoichiometry and against the state changes. Diagnostics only:
      ! never read by any reaction. isw_ledger = 0 (default) registers and
      ! computes nothing, bit-identical to the legacy build.
      call self%get_parameter(self%isw_ledger,'isw_ledger','', &
           'ledger counters: diagnostics of the source terms applied (0: off, 1: on)', &
           default=0,minimum=0,maximum=1)
      if (self%isw_ledger == 1) then
         call self%register_diagnostic_variable(self%id_ledger_diss_c,'ledger_diss_c','mg C/m^3/d', &
              'ledger: calcite dissolution -> free calcite',source=source_do)
         call self%register_diagnostic_variable(self%id_ledger_diss_O3c,'ledger_diss_O3c','mmol C/m^3/d', &
              'ledger: calcite dissolution -> dissolved inorganic carbon',source=source_do)
         call self%register_diagnostic_variable(self%id_ledger_diss_TA,'ledger_diss_TA','mmol/m^3/d', &
              'ledger: calcite dissolution -> total alkalinity',source=source_do)
      end if
   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_ersem_calcification), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: om_cal,L2c
      real(rk) :: fcalc,fdiss

      _LOOP_BEGIN_
         if (legacy_ersem_compatibility) then
            ! Legacy ERSEM includes background value, but this is inappropriate as it is used in a sink term.
            _GET_WITH_BACKGROUND_(self%id_c,L2c)
         else
            _GET_(self%id_c,L2c)
         end if
         _GET_(self%id_om_cal,om_cal)

         om_cal=max(om_cal,0._rk)

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

         fdiss = fdiss * self%sL2O3X

         _SET_ODE_(self%id_c,  -fdiss*L2c)
         _SET_DIAGNOSTIC_(self%id_L2O3c,-fdiss*L2c)
         _SET_ODE_(self%id_O3c, fdiss*L2c/CMass)
         _SET_ODE_(self%id_TA,2*fdiss*L2c/CMass)  ! Dissolution of CaCO3 increases alkalinity by 2 units
         if (self%isw_ledger == 1) then
            _SET_DIAGNOSTIC_(self%id_ledger_diss_c,  -fdiss*L2c)
            _SET_DIAGNOSTIC_(self%id_ledger_diss_O3c, fdiss*L2c/CMass)
            _SET_DIAGNOSTIC_(self%id_ledger_diss_TA,2*fdiss*L2c/CMass)
         end if
         _SET_DIAGNOSTIC_(self%id_RainR,fcalc * self%Rain0)
      _LOOP_END_

   end subroutine do

end module
