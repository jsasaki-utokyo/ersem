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
!
! E5 FORMULATIONS (2026-09-22, nippon-steel docs/125 s21, s21.9). All default
! OFF and skipped by branch, so every existing configuration is bit-identical.
!  * S0/S1 calcium-silicate dissolution SURROGATE, C2S as the stoichiometric
!    endmember: Ca2SiO4 + 4H2O -> 2Ca2+ + H4SiO4 + 4OH-. GROSS water-side
!    stoichiometry per mol Si: TA +r_ta_si (4 for C2S), DIC 0, Si +1 (to N5),
!    O2 0; the OH- converts CO2 into bicarbonate/carbonate inside DIC.
!    Rate R = k_sil + v_sil*[CO2*] (mmol Si/m^2/d); k_sil zero order (S0),
!    v_sil a transfer velocity on the bottom cell's CO2* (S1, a declared
!    phenomenological proxy for the acidity of the bed's respiration).
!    Constant accessibility is an ASSUMPTION: the donor state CaSi_stock
!    (initial CaSi_stock0, the declared accessible inventory) is debited for
!    conservation only and never limits the rate. Couplings N5s and the donor
!    are registered only when k_sil > 0 or v_sil > 0.
!  * L1 acid-promoted free lime: mode 4's flux times ([CO2*]/co2_ref)**n_co2;
!    n_co2 = 0 skips the factor (exactly mode 4). Mode 4 only.
!  * isw_ledger = 1 registers counters of every exchange this module applies
!    (CaO TA, carbonation DIC and solid, silicate TA, N5 and donor debit), for
!    the reaction-ledger gates (docs/119 FC1); diagnostics only.

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
      real(rk) :: k_sil            ! E5 S0: zero-order silicate dissolution (mmol Si/m^2/d); 0 = off
      real(rk) :: v_sil            ! E5 S1: CO2*-promoted silicate dissolution velocity (m/d); 0 = off
      real(rk) :: r_ta_si          ! E5: alkalinity per mol Si dissolved (eq/mol; 4 = C2S)
      real(rk) :: CaSi_stock0      ! E5: declared accessible silicate inventory (mmol Si/m^2), bookkeeping only
      real(rk) :: n_co2            ! E5 L1: exponent of the CO2* factor on mode 4; 0 = off
      real(rk) :: co2_ref          ! E5 L1: pivot CO2* (mmol/m^3) at which L1 equals mode 4
      logical  :: sil_on
      integer  :: isw_ledger

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

      ! E5 silicate surrogate
      type (type_state_variable_id)                 :: id_N5s        ! dissolved silicate (receives Si)
      type (type_bottom_state_variable_id)          :: id_sil_stock  ! donor: accessible silicate inventory
      type (type_horizontal_diagnostic_variable_id) :: id_sil_diss   ! silicate dissolution flux
      ! ledger counters (isw_ledger = 1 only)
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_cao_TA, id_ledger_carb_O3c, id_ledger_carb_c
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_sil_TA, id_ledger_sil_N5s, id_ledger_sil_stock

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

      ! E5 (docs/125 s21.9): parameters registered unconditionally with inert
      ! defaults; states, couplings and diagnostics only when a term is active.
      call self%get_parameter(self%k_sil, 'k_sil', 'mmol Si/m^2/d', &
         'E5 S0: zero-order calcium-silicate (C2S surrogate) dissolution rate (0 = off)', &
         default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%v_sil, 'v_sil', 'm/d', &
         'E5 S1: CO2*-promoted calcium-silicate dissolution velocity (0 = off)', &
         default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%r_ta_si, 'r_ta_si', 'eq/mol', &
         'E5: alkalinity released per mol Si dissolved (4 = C2S endmember)', &
         default=4.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%CaSi_stock0, 'CaSi_stock0', 'mmol Si/m^2', &
         'E5: declared accessible silicate inventory (donor bookkeeping; never rate-limiting)', &
         default=7000.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%n_co2, 'n_co2', '-', &
         'E5 L1: exponent of the ([CO2*]/co2_ref) factor on mode 4 (0 = off)', &
         default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%co2_ref, 'co2_ref', 'mmol/m^3', &
         'E5 L1: pivot CO2* at which the acid-promoted law equals mode 4', &
         default=15.0_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%isw_ledger, 'isw_ledger', '', &
         'ledger counters: diagnostics of the source terms applied (0: off, 1: on)', &
         default=0, minimum=0, maximum=1)
      if (self%n_co2 > 0.0_rk .and. self%iswCaO /= 4) &
         call self%fatal_error('initialize', 'n_co2 (E5 L1) applies to iswCaO = 4 only')
      self%sil_on = (self%k_sil > 0.0_rk .or. self%v_sil > 0.0_rk)
      if (self%sil_on) then
         call self%register_state_dependency(self%id_N5s, 'N5s', 'mmol Si/m^3', &
            'dissolved silicate (receives the silicate surrogate''s Si)')
         call self%register_bottom_state_variable(self%id_sil_stock, 'CaSi_stock', 'mmol Si/m^2', &
            'accessible calcium-silicate inventory (donor; bookkeeping only)', self%CaSi_stock0, minimum=0.0_rk)
         call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_sil_stock)
         call self%register_diagnostic_variable(self%id_sil_diss, 'silicate_dissolution', 'mmol Si/m^2/d', &
            'calcium-silicate surrogate dissolution flux', domain=domain_bottom, source=source_do_bottom)
      end if
      if (self%isw_ledger == 1) then
         call self%register_diagnostic_variable(self%id_ledger_cao_TA, 'ledger_cao_TA', 'mmol eq/m^2/d', &
              'ledger: CaO dissolution -> pelagic alkalinity (bottom flux)', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_carb_O3c, 'ledger_carb_O3c', 'mmol C/m^2/d', &
              'ledger: slag carbonation -> pelagic dissolved inorganic carbon (bottom flux)', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_carb_c, 'ledger_carb_c', 'mmol C/m^2/d', &
              'ledger: slag carbonation -> carbonated slag carbon', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_sil_TA, 'ledger_sil_TA', 'mmol eq/m^2/d', &
              'ledger: silicate surrogate -> pelagic alkalinity (bottom flux)', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_sil_N5s, 'ledger_sil_N5s', 'mmol Si/m^2/d', &
              'ledger: silicate surrogate -> pelagic silicate (bottom flux)', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_sil_stock, 'ledger_sil_stock', 'mmol Si/m^2/d', &
              'ledger: silicate surrogate -> donor inventory', source=source_do_bottom)
      end if

      ! Only proceed with registration if CaO dissolution is enabled (or the
      ! E5 silicate surrogate needs the alkalinity)
      if (self%iswCaO > 0 .or. self%sil_on) then
         ! Register dependency on total alkalinity (standard variable)
         call self%register_state_dependency(self%id_TA, &
            standard_variables%alkalinity_expressed_as_mole_equivalent)
      end if
      if (self%iswCaO > 0) then

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
      real(rk) :: co2aq, f_carb, f_sil

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
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_carb_O3c, -f_carb)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_carb_c, f_carb)
         end if
      _HORIZONTAL_LOOP_END_

      ! E5 S0/S1 calcium-silicate surrogate (docs/125 s21): independent of the
      ! dissolution mode; skipped entirely when inactive.
      if (self%sil_on) then
         _HORIZONTAL_LOOP_BEGIN_
            f_sil = self%k_sil
            if (self%v_sil > 0.0_rk) then
               _GET_(self%id_CO2aq, co2aq)
               f_sil = f_sil + self%v_sil * max(co2aq, 0.0_rk)
            end if
            _SET_BOTTOM_EXCHANGE_(self%id_TA, self%r_ta_si * f_sil)
            _SET_BOTTOM_EXCHANGE_(self%id_N5s, f_sil)
            _SET_BOTTOM_ODE_(self%id_sil_stock, -f_sil)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sil_diss, f_sil)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_sil_TA, self%r_ta_si * f_sil)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_sil_N5s, f_sil)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_sil_stock, -f_sil)
            end if
         _HORIZONTAL_LOOP_END_
      elseif (self%isw_ledger == 1) then
         _HORIZONTAL_LOOP_BEGIN_
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_sil_TA, 0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_sil_N5s, 0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_sil_stock, 0.0_rk)
         _HORIZONTAL_LOOP_END_
      end if

      ! Exit immediately if CaO dissolution is disabled
      if (self%iswCaO == 0) then
         if (self%isw_ledger == 1) then
            _HORIZONTAL_LOOP_BEGIN_
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_cao_TA, 0.0_rk)
            _HORIZONTAL_LOOP_END_
         end if
         return
      end if

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
               ! E5 L1 (docs/125 s21): acid-promoted free lime; n_co2 = 0 skips it
               if (self%n_co2 > 0.0_rk) then
                  _GET_(self%id_CO2aq, co2aq)
                  cao_flux = cao_flux * (max(co2aq, 0.0_rk) / self%co2_ref)**self%n_co2
               end if

               ! Deplete the accessible stock
               _SET_BOTTOM_ODE_(self%id_cao_stock, -cao_flux)
         end select

         ! Apply alkalinity flux to water column
         ! CaO dissolution increases alkalinity by 2 equivalents per mole
         _SET_BOTTOM_EXCHANGE_(self%id_TA, 2.0_rk * cao_flux)

         ! Set diagnostic
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_cao_diss, cao_flux)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_cao_TA, 2.0_rk * cao_flux)
         end if

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module ersem_benthic_cao