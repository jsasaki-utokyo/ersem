#include "fabm_driver.h"

!-----------------------------------------------------------------------
! Benthic sulfur cycle module for ERSEM
!
! Handles sulfate reduction in Layer 3 and sulfide oxidation in Layer 1.
! Simplified implementation (Option B): H2S + S0 only, no explicit SO4.
!
! Reactions:
!   Layer 3 (anoxic): OM -> H2S (sulfate reduction, coupled to H2 bacteria)
!   Layer 1 (oxic):   H2S + 0.5 O2 -> S0
!                     S0 + 1.5 O2 -> (removed from system)
!
! Key simplification: Layer 3 is by definition anoxic in ERSEM's 3-layer
! model, so no electron acceptor cascade check is needed.
!
! References:
!   - BROM model (Yakushev et al. 2017)
!   - ERSEM benthic_nitrogen_cycle.F90 for patterns
!-----------------------------------------------------------------------

module ersem_benthic_sulfur_cycle

   use fabm_types
   use ersem_shared

   implicit none
   private

   type, extends(type_base_model), public :: type_ersem_benthic_sulfur_cycle
      ! State variable dependencies (layer-specific via benthic_column_dissolved_matter)
      type(type_bottom_state_variable_id) :: id_H2S_1, id_H2S_3
      type(type_bottom_state_variable_id) :: id_S0_1
      type(type_bottom_state_variable_id) :: id_G2o  ! Oxygen in Layer 1

      ! Layer depth dependencies
      type(type_horizontal_dependency_id) :: id_D1m, id_D2m, id_Dtot

      ! Organic matter remineralization rate from H2 bacteria
      type(type_horizontal_dependency_id) :: id_remin_rate

      ! Diagnostic variables
      type(type_horizontal_diagnostic_variable_id) :: id_R_sulfate_red
      type(type_horizontal_diagnostic_variable_id) :: id_R_H2S_ox_ben
      type(type_horizontal_diagnostic_variable_id) :: id_R_S0_ox_ben

      ! Parameters
      real(rk) :: K_H2S_prod     ! H2S production rate per unit remineralization (mol S/mol C)
      real(rk) :: K_H2S_ox       ! H2S oxidation rate constant (1/d)
      real(rk) :: K_S0_ox        ! S0 oxidation rate constant (1/d)
      real(rk) :: K_O2_half      ! Half-saturation for oxygen (mmol/m3)

   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self, configunit)
      class(type_ersem_benthic_sulfur_cycle), intent(inout), target :: self
      integer, intent(in) :: configunit

      ! Set time unit to d-1 (ERSEM convention)
      self%dt = 86400._rk

      ! Get parameters
      call self%get_parameter(self%K_H2S_prod, 'K_H2S_prod', 'mol S/mol C', &
           'H2S production per C remineralized (stoichiometry)', default=0.5_rk)
      call self%get_parameter(self%K_H2S_ox, 'K_H2S_ox', '1/d', &
           'H2S oxidation rate constant', default=0.5_rk)
      call self%get_parameter(self%K_S0_ox, 'K_S0_ox', '1/d', &
           'S0 oxidation rate constant', default=0.02_rk)
      call self%get_parameter(self%K_O2_half, 'K_O2_half', 'mmol/m^3', &
           'half-saturation O2 for oxidation', default=1.0_rk)

      ! Register dependencies for layer-specific sulfur variables
      ! These link to variables created by benthic_column_dissolved_matter with composition 'h' and 'e'
      call self%register_state_dependency(self%id_H2S_1, 'H2S_1', 'mmol S/m^2', &
           'hydrogen sulfide in layer 1')
      call self%register_state_dependency(self%id_H2S_3, 'H2S_3', 'mmol S/m^2', &
           'hydrogen sulfide in layer 3')
      call self%register_state_dependency(self%id_S0_1, 'S0_1', 'mmol S/m^2', &
           'elemental sulfur in layer 1')

      ! Oxygen in Layer 1 for oxidation reactions
      call self%register_state_dependency(self%id_G2o, 'G2o', 'mmol O_2/m^2', &
           'oxygen in layer 1')

      ! Layer depths
      call self%register_dependency(self%id_D1m, depth_of_bottom_interface_of_layer_1)
      call self%register_dependency(self%id_D2m, depth_of_bottom_interface_of_layer_2)
      call self%register_dependency(self%id_Dtot, depth_of_sediment_column)

      ! Organic matter remineralization rate from H2 bacteria
      ! This should be coupled to the H2 bacteria respiration output
      call self%register_dependency(self%id_remin_rate, 'remin_rate', 'mmol C/m^2/d', &
           'anaerobic remineralization rate in layer 3')

      ! Diagnostic variables
      call self%register_diagnostic_variable(self%id_R_sulfate_red, 'R_sulfate_red', &
           'mmol S/m^2/d', 'sulfate reduction rate (H2S production)', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_H2S_ox_ben, 'R_H2S_ox_ben', &
           'mmol S/m^2/d', 'benthic H2S oxidation rate', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_S0_ox_ben, 'R_S0_ox_ben', &
           'mmol S/m^2/d', 'benthic S0 oxidation rate', &
           domain=domain_bottom, source=source_do_bottom)

   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_ersem_benthic_sulfur_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: H2S_1, H2S_3, S0_1, G2o
      real(rk) :: D1m, remin_rate
      real(rk) :: O2_conc_1, f_O2
      real(rk) :: R_sulfate_red, R_H2S_ox_1, R_S0_ox_1

      _HORIZONTAL_LOOP_BEGIN_

         ! Get state variables
         _GET_HORIZONTAL_(self%id_H2S_1, H2S_1)
         _GET_HORIZONTAL_(self%id_H2S_3, H2S_3)
         _GET_HORIZONTAL_(self%id_S0_1, S0_1)
         _GET_HORIZONTAL_(self%id_G2o, G2o)

         ! Get layer depth
         _GET_HORIZONTAL_(self%id_D1m, D1m)

         ! Get remineralization rate from H2 bacteria
         _GET_HORIZONTAL_(self%id_remin_rate, remin_rate)

         ! Ensure non-negative
         H2S_1 = max(0.0_rk, H2S_1)
         H2S_3 = max(0.0_rk, H2S_3)
         S0_1 = max(0.0_rk, S0_1)
         G2o = max(0.0_rk, G2o)
         remin_rate = max(0.0_rk, remin_rate)

         ! ============================================================
         ! LAYER 3: Sulfate reduction (always active - Layer 3 is anoxic)
         ! ============================================================
         ! H2S production is proportional to anaerobic OM remineralization
         ! Stoichiometry: 53 SO4 per 106 C -> 0.5 mol S per mol C
         R_sulfate_red = self%K_H2S_prod * remin_rate

         ! ============================================================
         ! LAYER 1: H2S and S0 oxidation (limited by O2 availability)
         ! ============================================================
         ! Convert depth-integrated O2 to concentration
         ! Guard against division by very small D1m (layer collapse)
         O2_conc_1 = G2o / max(D1m, 0.0001_rk)

         ! Oxygen limitation (Michaelis-Menten)
         f_O2 = O2_conc_1 / (O2_conc_1 + self%K_O2_half)

         ! H2S + 0.5 O2 -> S0
         R_H2S_ox_1 = self%K_H2S_ox * H2S_1 * f_O2

         ! S0 + 1.5 O2 -> (removed)
         R_S0_ox_1 = self%K_S0_ox * S0_1 * f_O2

         ! ============================================================
         ! Set ODEs
         ! ============================================================
         ! Layer 3: H2S production from sulfate reduction
         _SET_BOTTOM_ODE_(self%id_H2S_3, R_sulfate_red)

         ! Layer 1: H2S consumption, S0 production/consumption, O2 consumption
         _SET_BOTTOM_ODE_(self%id_H2S_1, -R_H2S_ox_1)
         _SET_BOTTOM_ODE_(self%id_S0_1,   R_H2S_ox_1 - R_S0_ox_1)
         _SET_BOTTOM_ODE_(self%id_G2o,   -0.5_rk * R_H2S_ox_1 - 1.5_rk * R_S0_ox_1)

         ! Set diagnostics
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_sulfate_red, R_sulfate_red)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_H2S_ox_ben, R_H2S_ox_1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_S0_ox_ben, R_S0_ox_1)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module ersem_benthic_sulfur_cycle
