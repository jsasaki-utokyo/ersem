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
!   All layers:       H2S + Fe(II) -> FeS (burial, irreversible)
!
! Key simplification: Layer 3 is by definition anoxic in ERSEM's 3-layer
! model, so no electron acceptor cascade check is needed.
!
! FeS Precipitation (Iron Sulfide Burial):
!   Free sulfide is buffered by reactive iron, forming FeS which is buried.
!   This is the primary irreversible sink for H2S in marine sediments.
!   Without this sink, H2S accumulates unrealistically.
!   Implemented as first-order removal: R_FeS = k_FeS * H2S
!   where k_FeS represents effective iron availability and reaction rate.
!
! Oxic Barrier Mechanism:
!   When an oxic layer (Layer 1) exists, H2S diffusing upward from deeper
!   layers is chemically oxidized before reaching the pelagic. This is
!   modeled as a "barrier" that suppresses the H2S flux to the pelagic:
!     flux_suppression = 1 - exp(-K_barrier * D1m * f_O2)
!   where K_barrier controls barrier effectiveness, D1m is oxic layer depth,
!   and f_O2 is oxygen availability. Any H2S that would have escaped is
!   oxidized to S0 within the benthic system.
!
! References:
!   - BROM model (Yakushev et al. 2017)
!   - ERSEM benthic_nitrogen_cycle.F90 for patterns
!   - Luther et al. (2011) Thermodynamics and Kinetics of Sulfide Oxidation
!   - Rickard & Luther (2007) Chemistry of Iron Sulfides
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

      ! Pelagic H2S at bottom for oxic barrier mechanism
      type(type_state_variable_id) :: id_H2S_pel
      type(type_state_variable_id) :: id_S0_pel
      type(type_state_variable_id) :: id_O2_pel

      ! Layer depth dependencies
      type(type_horizontal_dependency_id) :: id_D1m, id_D2m, id_Dtot

      ! Organic matter remineralization rate from H2 bacteria
      type(type_horizontal_dependency_id) :: id_remin_rate

      ! Diagnostic variables
      type(type_horizontal_diagnostic_variable_id) :: id_R_sulfate_red
      type(type_horizontal_diagnostic_variable_id) :: id_R_H2S_ox_ben
      type(type_horizontal_diagnostic_variable_id) :: id_R_S0_ox_ben
      type(type_horizontal_diagnostic_variable_id) :: id_R_barrier_ox
      type(type_horizontal_diagnostic_variable_id) :: id_R_FeS_ben
      type(type_horizontal_diagnostic_variable_id) :: id_R_FeS_pel
      type(type_horizontal_diagnostic_variable_id) :: id_f_barrier

      ! Parameters
      real(rk) :: K_H2S_prod     ! H2S production rate per unit remineralization (mol S/mol C)
      real(rk) :: K_H2S_ox       ! H2S oxidation rate constant (1/d)
      real(rk) :: K_S0_ox        ! S0 oxidation rate constant (1/d)
      real(rk) :: K_O2_half      ! Half-saturation for oxygen (mmol/m3)
      real(rk) :: K_barrier      ! Oxic barrier effectiveness (1/m)
      real(rk) :: K_barrier_rate ! Rate at which barrier oxidizes H2S (1/d)
      real(rk) :: K_FeS_ben      ! FeS precipitation rate in benthic layers (1/d)
      real(rk) :: K_FeS_pel      ! FeS precipitation rate in pelagic (1/d)

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

      ! Oxic barrier parameters
      ! K_barrier: controls how effective the oxic layer is at blocking H2S
      ! Physical basis: exp(-K_barrier * D1m) represents the fraction of H2S
      ! that can pass through an oxic layer of depth D1m without being oxidized.
      ! Typical value: 500-2000 1/m (very effective barrier even for thin layers)
      call self%get_parameter(self%K_barrier, 'K_barrier', '1/m', &
           'oxic barrier effectiveness (higher = stronger barrier)', default=1000.0_rk)
      call self%get_parameter(self%K_barrier_rate, 'K_barrier_rate', '1/d', &
           'rate of barrier H2S oxidation in bottom water', default=100.0_rk)

      ! FeS precipitation parameters (iron sulfide burial)
      ! This is the key irreversible sink that prevents runaway H2S accumulation.
      ! In real sediments, reactive iron (Fe(II) from Fe(III) reduction) rapidly
      ! precipitates with H2S to form FeS, which is then buried or further
      ! converted to pyrite (FeS2). This effectively removes sulfide from the
      ! active biogeochemical cycle.
      ! K_FeS_ben: rate constant for FeS precipitation in benthic layers
      ! Typical values: 0.1-1.0 1/d (represents effective Fe availability)
      call self%get_parameter(self%K_FeS_ben, 'K_FeS_ben', '1/d', &
           'FeS precipitation rate in benthic layers (iron sulfide burial)', default=0.5_rk)
      ! K_FeS_pel: rate constant for FeS precipitation in pelagic bottom water
      ! Usually lower than benthic because less reactive Fe available in water column
      call self%get_parameter(self%K_FeS_pel, 'K_FeS_pel', '1/d', &
           'FeS precipitation rate in pelagic (scavenging)', default=0.1_rk)

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

      ! Pelagic variables at bottom for oxic barrier mechanism
      call self%register_state_dependency(self%id_H2S_pel, 'H2S_pel', 'mmol S/m^3', &
           'pelagic hydrogen sulfide')
      call self%register_state_dependency(self%id_S0_pel, 'S0_pel', 'mmol S/m^3', &
           'pelagic elemental sulfur')
      call self%register_state_dependency(self%id_O2_pel, 'O2_pel', 'mmol O_2/m^3', &
           'pelagic oxygen')

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
      call self%register_diagnostic_variable(self%id_R_barrier_ox, 'R_barrier_ox', &
           'mmol S/m^2/d', 'H2S oxidation by oxic barrier', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_FeS_ben, 'R_FeS_ben', &
           'mmol S/m^2/d', 'FeS precipitation in benthic (H2S removal)', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_FeS_pel, 'R_FeS_pel', &
           'mmol S/m^2/d', 'FeS precipitation in pelagic (H2S scavenging)', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_f_barrier, 'f_barrier', &
           '-', 'oxic barrier suppression factor (0=no barrier, 1=complete)', &
           domain=domain_bottom, source=source_do_bottom)

   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_ersem_benthic_sulfur_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: H2S_1, H2S_3, S0_1, G2o
      real(rk) :: H2S_pel, S0_pel, O2_pel
      real(rk) :: D1m, remin_rate
      real(rk) :: O2_conc_1, f_O2, f_O2_pel
      real(rk) :: R_sulfate_red, R_H2S_ox_1, R_S0_ox_1
      real(rk) :: f_barrier, R_barrier_ox
      real(rk) :: R_FeS_1, R_FeS_3, R_FeS_ben, R_FeS_pel

      _HORIZONTAL_LOOP_BEGIN_

         ! Get state variables
         _GET_HORIZONTAL_(self%id_H2S_1, H2S_1)
         _GET_HORIZONTAL_(self%id_H2S_3, H2S_3)
         _GET_HORIZONTAL_(self%id_S0_1, S0_1)
         _GET_HORIZONTAL_(self%id_G2o, G2o)

         ! Get pelagic variables at bottom
         _GET_(self%id_H2S_pel, H2S_pel)
         _GET_(self%id_S0_pel, S0_pel)
         _GET_(self%id_O2_pel, O2_pel)

         ! Get layer depth
         _GET_HORIZONTAL_(self%id_D1m, D1m)

         ! Get remineralization rate from H2 bacteria
         _GET_HORIZONTAL_(self%id_remin_rate, remin_rate)

         ! Ensure non-negative
         H2S_1 = max(0.0_rk, H2S_1)
         H2S_3 = max(0.0_rk, H2S_3)
         S0_1 = max(0.0_rk, S0_1)
         G2o = max(0.0_rk, G2o)
         H2S_pel = max(0.0_rk, H2S_pel)
         S0_pel = max(0.0_rk, S0_pel)
         O2_pel = max(0.0_rk, O2_pel)
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
         ! OXIC BARRIER MECHANISM
         ! ============================================================
         ! When an oxic layer exists (D1m > 0) and O2 is available, the
         ! oxic layer acts as a chemical barrier that oxidizes H2S before
         ! it can escape to the pelagic. This is implemented as an
         ! additional oxidation term for bottom-water H2S.
         !
         ! The barrier factor represents what fraction of H2S would be
         ! oxidized while diffusing through the oxic layer:
         !   f_barrier = 1 - exp(-K_barrier * D1m * f_O2_pel)
         !
         ! When f_barrier ~ 1 (thick oxic layer with high O2), nearly all
         ! H2S entering the bottom water is immediately oxidized.
         ! When f_barrier ~ 0 (no oxic layer or no O2), H2S passes freely.
         !
         ! Physical basis: The oxidation rate of H2S with metal catalysis
         ! (Fe, Mn) is very fast (half-life of minutes to hours). With
         ! K_barrier = 1000 1/m, even a 2mm oxic layer gives:
         !   f_barrier = 1 - exp(-1000 * 0.002 * 1) = 0.86 (86% blocked)
         ! and a 5mm layer gives 99.3% blocking.

         ! Oxygen limitation for barrier (based on bottom water O2)
         f_O2_pel = O2_pel / (O2_pel + self%K_O2_half)

         ! Barrier suppression factor
         ! Use effective D1m with minimum of 1mm (0.001m) to account for:
         ! 1. Diffusive boundary layer at sediment-water interface
         ! 2. The fact that H2S oxidation occurs at the interface even when
         !    sediment oxic layer is thin, as long as bottom water O2 is present
         ! This prevents barrier failure when D1m collapses during anoxic events
         ! but bottom water O2 has recovered.
         f_barrier = 1.0_rk - exp(-self%K_barrier * max(D1m, 0.001_rk) * f_O2_pel)

         ! H2S oxidation rate by barrier (removes H2S from bottom water)
         ! This is proportional to H2S concentration and barrier strength
         R_barrier_ox = self%K_barrier_rate * H2S_pel * f_barrier

         ! ============================================================
         ! FeS PRECIPITATION (IRON SULFIDE BURIAL)
         ! ============================================================
         ! This is the key irreversible sink for H2S that prevents runaway
         ! accumulation. In real sediments, reactive iron (from Fe(III)
         ! reduction in anoxic zones) rapidly precipitates with H2S to form
         ! iron sulfides (FeS, eventually pyrite FeS2).
         !
         ! The reaction: H2S + Fe(II) -> FeS(s) + 2H+
         !
         ! Since we don't explicitly track benthic Fe(II), we use a first-order
         ! parameterization where the rate constant K_FeS represents effective
         ! reactive iron availability times the intrinsic reaction rate.
         !
         ! FeS precipitation in benthic Layer 1 (oxic layer)
         ! Lower rate because Fe is mostly in oxidized form (Fe(III))
         R_FeS_1 = self%K_FeS_ben * 0.1_rk * H2S_1  ! Reduced rate in oxic layer

         ! FeS precipitation in benthic Layer 3 (anoxic layer)
         ! Higher rate because Fe(II) is abundant from Fe(III) reduction
         R_FeS_3 = self%K_FeS_ben * H2S_3

         ! Total benthic FeS precipitation
         R_FeS_ben = R_FeS_1 + R_FeS_3

         ! FeS precipitation (scavenging) in pelagic bottom water
         ! Removes H2S through reaction with particulate Fe or settling FeS
         R_FeS_pel = self%K_FeS_pel * H2S_pel

         ! ============================================================
         ! Set ODEs
         ! ============================================================
         ! Layer 3: H2S production from sulfate reduction, loss from FeS burial
         _SET_BOTTOM_ODE_(self%id_H2S_3, R_sulfate_red - R_FeS_3)

         ! Layer 1: H2S consumption by oxidation and FeS precipitation
         _SET_BOTTOM_ODE_(self%id_H2S_1, -R_H2S_ox_1 - R_FeS_1)
         _SET_BOTTOM_ODE_(self%id_S0_1,   R_H2S_ox_1 - R_S0_ox_1)
         _SET_BOTTOM_ODE_(self%id_G2o,   -0.5_rk * R_H2S_ox_1 - 1.5_rk * R_S0_ox_1)

         ! Pelagic: H2S removal by oxic barrier oxidation and FeS scavenging
         ! Barrier oxidation produces S0, FeS scavenging is irreversible removal
         ! Note: Using _SET_BOTTOM_EXCHANGE_ applies flux to bottom cell only
         _SET_BOTTOM_EXCHANGE_(self%id_H2S_pel, -R_barrier_ox - R_FeS_pel)
         _SET_BOTTOM_EXCHANGE_(self%id_S0_pel,   R_barrier_ox)
         _SET_BOTTOM_EXCHANGE_(self%id_O2_pel,  -0.5_rk * R_barrier_ox)

         ! Set diagnostics
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_sulfate_red, R_sulfate_red)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_H2S_ox_ben, R_H2S_ox_1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_S0_ox_ben, R_S0_ox_1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_barrier_ox, R_barrier_ox)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_FeS_ben, R_FeS_ben)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_FeS_pel, R_FeS_pel)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_barrier, f_barrier)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module ersem_benthic_sulfur_cycle
