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
      ! ledger counters (isw_ledger = 1 only; nippon-steel docs/119 FC1)
      integer  :: isw_ledger
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_diss_c,id_ledger_diss_O3c,id_ledger_diss_TA
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_prec_c,id_ledger_prec_O3c,id_ledger_prec_TA

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
      ! The precipitation counters are registered whether or not K_prec > 0 and
      ! read exactly zero while the precipitation block is inactive.
      call self%get_parameter(self%isw_ledger, 'isw_ledger', '', &
           'ledger counters: diagnostics of the source terms applied (0: off, 1: on)', &
           default=0, minimum=0, maximum=1)
      if (self%isw_ledger == 1) then
         call self%register_diagnostic_variable(self%id_ledger_diss_c, 'ledger_diss_c', 'mg C/m^2/d', &
              'ledger: calcite dissolution -> benthic calcite', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_diss_O3c, 'ledger_diss_O3c', 'mmol C/m^2/d', &
              'ledger: calcite dissolution -> pelagic dissolved inorganic carbon (bottom flux)', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_diss_TA, 'ledger_diss_TA', 'mmol eq/m^2/d', &
              'ledger: calcite dissolution -> pelagic alkalinity (bottom flux)', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_prec_c, 'ledger_prec_c', 'mg C/m^2/d', &
              'ledger: light-driven calcification -> benthic calcite', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_prec_O3c, 'ledger_prec_O3c', 'mmol C/m^2/d', &
              'ledger: light-driven calcification -> pelagic dissolved inorganic carbon (bottom flux)', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_prec_TA, 'ledger_prec_TA', 'mmol eq/m^2/d', &
              'ledger: light-driven calcification -> pelagic alkalinity (bottom flux)', source=source_do_bottom)
      end if

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
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_c, -fdiss*bL2c)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_O3c, fdiss*bL2c/CMass)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_TA, 2*fdiss*bL2c/CMass)
         end if

         ! Light-driven mat calcification (docs/15; inert when K_prec = 0)
         F_prec = 0.0_rk
         if (self%isw_ledger == 1) then
            ! zero unless the precipitation block below applies its terms
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_prec_c, 0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_prec_O3c, 0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_prec_TA, 0.0_rk)
         end if
         if (self%K_prec > 0.0_rk) then
            _GET_(self%id_par, par)
            F_prec = self%K_prec &
                   * max(0.0_rk, par) / (max(0.0_rk, par) + self%K_par_prec) &
                   * (max(0.0_rk, om_cal - 1.0_rk))**self%n_prec
            _SET_BOTTOM_ODE_(self%id_c, F_prec)
            _SET_BOTTOM_EXCHANGE_(self%id_O3c, -F_prec/CMass)
            _SET_BOTTOM_EXCHANGE_(self%id_TA, -2*F_prec/CMass)  ! Precipitation of CaCO3 removes 2 units of alkalinity
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_prec_c, F_prec)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_prec_O3c, -F_prec/CMass)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_prec_TA, -2*F_prec/CMass)
            end if
         end if
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_precipitation, F_prec)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
