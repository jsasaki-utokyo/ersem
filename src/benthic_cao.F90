#include "fabm_driver.h"

! Benthic CaO dissolution module for ocean alkalinity enhancement
! Simulates alkalinity release from steelmaking slag containing CaO
! CaO + H2O → Ca(OH)2 → Ca2+ + 2OH-
! Each mole of CaO increases alkalinity by 2 equivalents
!
! SLAG CARBONATION (2026-09-11, nippon-steel docs/114 EG; mechanism M1 of
! docs/reviews/2026-09-11-fable-late-calcia-carbon.md). Residual portlandite
! in the slag takes up dissolved CO2 in the solid phase:
!    Ca(OH)2(s) + CO2(aq) -> CaCO3(s) + H2O
! The water loses one mole of DIC per mole reacted and NO alkalinity (neither
! the reactant solid nor the product solid is in the water), and no oxygen.
! Rate = k_carb * [CO2*], k_carb a transfer velocity (m/d), [CO2*] the bottom
! cell's carbonic acid (carbonate module diagnostic CarbA, mmol/m^3); first
! order in CO2*, so it is strongest where pCO2 is highest (late in a dark
! closure). The carbon removed is kept in the bottom state CaCO3_slag so total
! carbon is conserved. Independent of iswCaO; k_carb = 0 (default) is OFF and
! bit-identical. Not represented: exhaustion or passivation of the portlandite
! (negligible over a 72-h closure at the rates of interest), temperature.

module ersem_benthic_cao

   use fabm_types
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_benthic_cao
      ! Parameters
      integer  :: iswCaO           ! CaO dissolution mode (0=off, 1=constant, 2=pH-dependent, 3=stock depletion, 4=surface passivation)
      real(rk) :: CaO_flux_rate    ! Base dissolution flux (mmol/m^2/d)
      real(rk) :: k_CaO_diss       ! Dissolution rate constant for stock depletion (1/d)
      real(rk) :: pH_factor        ! pH sensitivity factor for pH-dependent mode
      real(rk) :: temp_Q10         ! Q10 temperature factor
      real(rk) :: CaO_half_sat     ! Half-saturation constant for stock limitation (mmol/m^2)
      real(rk) :: CaO_stock0       ! Initial CaO stock (mmol/m^2)
      real(rk) :: k_carb           ! slag carbonation transfer velocity (m/d); 0 = off

      ! State variables and dependencies
      type (type_bottom_state_variable_id)          :: id_cao_stock  ! CaO stock (if tracking)
      type (type_state_variable_id)                 :: id_TA         ! Total alkalinity
      type (type_dependency_id)                     :: id_pH         ! pH (for pH-dependent mode)
      type (type_dependency_id)                     :: id_temp       ! Temperature
      type (type_state_variable_id)                 :: id_O3c        ! DIC (carbonation sink)
      type (type_dependency_id)                     :: id_CO2aq      ! CO2* (carbonation driver)
      type (type_bottom_state_variable_id)          :: id_carb_c     ! carbonated slag carbon
      type (type_horizontal_diagnostic_variable_id) :: id_carb_flux  ! carbonation flux

      ! Diagnostics
      type (type_horizontal_diagnostic_variable_id) :: id_cao_diss   ! CaO dissolution flux

   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self, configunit)
      class (type_ersem_benthic_cao), intent(inout), target :: self
      integer,                         intent(in)           :: configunit

      ! Set time unit to d-1
      self%dt = 86400._rk

      ! Register parameters with safe defaults that preserve existing behavior
      call self%get_parameter(self%iswCaO, 'iswCaO', '', &
         'CaO dissolution mode (0: off, 1: constant flux, 2: pH-dependent, 3: stock depletion, 4: surface passivation)', &
         default=0, minimum=0, maximum=4)

      call self%get_parameter(self%CaO_flux_rate, 'CaO_flux_rate', 'mmol/m^2/d', &
         'base CaO dissolution flux', &
         default=0.0_rk, minimum=0.0_rk)

      call self%get_parameter(self%k_CaO_diss, 'k_CaO_diss', '1/d', &
         'first-order dissolution rate constant (mode 4: F = k * accessible stock)', &
         default=0.01_rk, minimum=0.0_rk)

      call self%get_parameter(self%pH_factor, 'pH_factor', '-', &
         'pH sensitivity factor (enhanced dissolution at low pH)', &
         default=2.0_rk, minimum=0.0_rk)

      call self%get_parameter(self%temp_Q10, 'temp_Q10', '-', &
         'Q10 temperature factor for dissolution', &
         default=2.0_rk, minimum=1.0_rk)

      call self%get_parameter(self%CaO_stock0, 'CaO_stock0', 'mmol/m^2', &
         'initial CaO stock at bottom (only used in mode 3)', &
         default=1000.0_rk, minimum=0.0_rk)

      ! Slag carbonation (docs/114 EG). Registered unconditionally so the same
      ! coupling block works in every arm; k_carb = 0 leaves it inert.
      call self%get_parameter(self%k_carb, 'k_carb', 'm/d', &
         'slag carbonation transfer velocity: CO2(aq) taken up by residual portlandite (0 = off)', &
         default=0.0_rk, minimum=0.0_rk)
      call self%register_state_dependency(self%id_O3c, 'O3c', 'mmol C/m^3', &
         'dissolved inorganic carbon (carbonation sink)')
      call self%register_dependency(self%id_CO2aq, 'CO2aq', 'mmol/m^3', &
         'carbonic acid concentration CO2* (carbonation driver)')
      call self%register_bottom_state_variable(self%id_carb_c, 'CaCO3_slag', 'mmol C/m^2', &
         'carbon fixed by slag carbonation', 0.0_rk, minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_carb_c)
      call self%register_diagnostic_variable(self%id_carb_flux, 'carbonation', 'mmol C/m^2/d', &
         'slag carbonation flux (DIC sink, alkalinity-neutral)', &
         domain=domain_bottom, source=source_do_bottom)

      ! Only proceed with registration if CaO dissolution is enabled
      if (self%iswCaO > 0) then
         ! Register dependency on total alkalinity (standard variable)
         call self%register_state_dependency(self%id_TA, &
            standard_variables%alkalinity_expressed_as_mole_equivalent)

         ! Register diagnostic for CaO dissolution flux
         call self%register_diagnostic_variable(self%id_cao_diss, 'CaO_dissolution', 'mmol/m^2/d', &
            'CaO dissolution flux from bottom slag', &
            domain=domain_bottom, source=source_do_bottom)

         ! Register dependencies for pH-dependent mode
         if (self%iswCaO == 2) then
            call self%register_dependency(self%id_pH, standard_variables%ph_reported_on_total_scale)
            call self%register_dependency(self%id_temp, standard_variables%temperature)
         end if

         ! Register stock state variable for the stock-carrying modes (3, 4)
         if (self%iswCaO == 3 .or. self%iswCaO == 4) then
            call self%register_bottom_state_variable(self%id_cao_stock, 'CaO_stock', 'mmol/m^2', &
               'CaO stock at bottom', &
               self%CaO_stock0, minimum=0.0_rk)
         end if
         if (self%iswCaO == 3) then
            call self%get_parameter(self%CaO_half_sat, 'CaO_half_sat', 'mmol/m^2', &
               'half-saturation constant for stock limitation', &
               default=100.0_rk, minimum=0.0_rk)

            call self%register_dependency(self%id_pH, standard_variables%ph_reported_on_total_scale)
            call self%register_dependency(self%id_temp, standard_variables%temperature)
         end if
      end if

   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_cao), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: cao_flux, pH, temp, stock
      real(rk) :: f_pH, f_temp, f_stock
      real(rk) :: co2aq, f_carb

      ! Slag carbonation first: it is independent of the dissolution mode.
      _HORIZONTAL_LOOP_BEGIN_
         f_carb = 0.0_rk
         if (self%k_carb > 0.0_rk) then
            _GET_(self%id_CO2aq, co2aq)
            f_carb = self%k_carb * max(co2aq, 0.0_rk)
            _SET_BOTTOM_EXCHANGE_(self%id_O3c, -f_carb)
            _SET_BOTTOM_ODE_(self%id_carb_c, f_carb)
         end if
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_carb_flux, f_carb)
      _HORIZONTAL_LOOP_END_

      ! Exit immediately if CaO dissolution is disabled
      if (self%iswCaO == 0) return

      _HORIZONTAL_LOOP_BEGIN_

         select case (self%iswCaO)
            case (1)  ! Constant flux mode
               cao_flux = self%CaO_flux_rate

            case (2)  ! pH-dependent mode
               _GET_(self%id_pH, pH)
               _GET_(self%id_temp, temp)

               ! pH function: enhanced dissolution at low pH
               ! Dissolution increases as pH drops below 8.3 (typical seawater)
               f_pH = max(0.0_rk, self%pH_factor * (8.3_rk - pH))

               ! Temperature function using Q10
               f_temp = self%temp_Q10 ** ((temp - 10.0_rk) / 10.0_rk)

               cao_flux = self%CaO_flux_rate * f_pH * f_temp

            case (3)  ! Stock depletion mode
               _GET_HORIZONTAL_(self%id_cao_stock, stock)

               if (stock > 0.0_rk) then
                  _GET_(self%id_pH, pH)
                  _GET_(self%id_temp, temp)

                  ! pH and temperature functions as above
                  f_pH = max(0.0_rk, self%pH_factor * (8.3_rk - pH))
                  f_temp = self%temp_Q10 ** ((temp - 10.0_rk) / 10.0_rk)
                  ! Flux depends on remaining stock
                  ! bug: cao_flux = self%k_CaO_diss * stock * f_pH * f_temp
                  ! Stock limitation using Michaelis-Menten kinetics
                  f_stock = stock / (stock + self%CaO_half_sat)

                  cao_flux = self%CaO_flux_rate * f_pH * f_temp * f_stock

                  ! Deplete the stock
                  _SET_BOTTOM_ODE_(self%id_cao_stock, -cao_flux)
               else
                  cao_flux = 0.0_rk
               end if

            case (4)  ! Surface-passivation mode: first-order in the ACCESSIBLE stock
               _GET_HORIZONTAL_(self%id_cao_stock, stock)

               ! F = k * S with S the accessible SURFACE stock, not the bulk
               ! inventory: its depletion IS the passivation, so the flux
               ! declines exponentially with e-folding time 1/k. CaO_stock0
               ! is therefore a small fitted number (order 1e2 mmol/m^2),
               ! orders below the bulk Ca inventory. No pH or temperature
               ! factor: the 2024-10 dissolution series shows Q10 ~ 1 and
               ! tank pH cannot constrain the pH law.
               cao_flux = self%k_CaO_diss * max(stock, 0.0_rk)

               ! Deplete the accessible stock
               _SET_BOTTOM_ODE_(self%id_cao_stock, -cao_flux)
         end select

         ! Apply alkalinity flux to water column
         ! CaO dissolution increases alkalinity by 2 equivalents per mole
         _SET_BOTTOM_EXCHANGE_(self%id_TA, 2.0_rk * cao_flux)

         ! Set diagnostic
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_cao_diss, cao_flux)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module ersem_benthic_cao