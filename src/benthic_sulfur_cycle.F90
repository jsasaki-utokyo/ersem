#include "fabm_driver.h"

!-----------------------------------------------------------------------
! Benthic sulfur cycle module for ERSEM
!
! Handles sulfate reduction in Layer 3 and sulfide oxidation in Layers 1-2.
! Simplified implementation (Option B): H2S + S0 only, no explicit SO4.
!
! Reactions:
!   Layer 3 (anoxic): OM -> H2S (sulfate reduction, coupled to H2 bacteria)
!   Layer 2 (suboxic): H2S + 0.4 NO3 -> S0 + 0.2 N2 (chemolithotrophic oxidation)
!   Layer 1 (oxic):   H2S + 0.5 O2 -> S0
!                     S0 + 1.5 O2 -> SO4 (removed from system)
!                     S0 -> burial (settling into deeper sediment)
!   All layers:       H2S + Fe(II) -> FeS (burial, irreversible)
!
! Key simplification: Layer 3 is by definition anoxic in ERSEM's 3-layer
! model, so no electron acceptor cascade check is needed.
!
! H2S-NO3 Oxidation (Layer 2 - Suboxic Zone):
!   H2S diffusing from Layer 3 through Layer 2 is oxidized by nitrate
!   through chemolithotrophic sulfur oxidation coupled to nitrate reduction.
!   The nitrogen product is partitioned by f_DNRA (0-1):
!     Denitrification: 5 H2S + 2 NO3- -> 5 S0 + N2 + 4 H2O
!     DNRA:            4 H2S + NO3- + 2H+ -> 4 S0 + NH4+ + 3 H2O
!     Mixed:           r_NO3 = 2/(5+3*f_DNRA) mol NO3 per mol H2S
!   This prevents H2S from reaching Layer 1 and the pelagic when
!   a thick suboxic zone exists (as in winter).
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
!   - Brunet & Garcia-Gil (1996) Sulfide-induced dissimilatory nitrate reduction
!-----------------------------------------------------------------------

module ersem_benthic_sulfur_cycle

   use fabm_types
   use ersem_shared

   implicit none
   private

   ! SR3 fixed grid (jsasaki 2026-09-23; nippon-steel docs/127 s9, review round 40): the bed sulfide as dynamic
   ! pore-water inventories on a grid FIXED in absolute depth (isw_h2s_layers = 2), replacing the homogeneous G2_H2S
   ! column. Cells: top width h2s_grid_top growing by h2s_grid_growth up to h2s_grid_cap, then equal cells no wider
   ! than the cap down to the column base (tools/sr3_grid_reference.py mesh(), the frozen instrument). A cell cut by an
   ! ERSEM interface takes each layer's first-order sink constant and sulfate-reduction placement in proportion to
   ! the length of the cell inside that layer; the conductance between cell centres is 1/int(dz/diff) (ERSEM's
   ! convention: areal flux = diff * dc/dz with pore-water c = m/(poro*dz)); the surface exchange is
   ! J = (c_1 - c_pel)/(cmix + int_0^{zc_1} dz/diff) with the uptake limiter h_supply_h2s. Three callbacks keep FABM's
   ! graph acyclic: the parent's reactions read states and interface depths only; a summary child reports the bed
   ! total (read by the H2 bacteria's H2S_col) and the ERSEM-layer totals; a transport child computes the fluxes.
   ! The moving-sub-box variant (isw_h2s_layers = 1, SR3-B) failed its registered transport gate (docs/127 s8.3)
   ! and is retired.
   type, extends(type_base_model) :: type_h2s_bed_summary
      integer :: n = 0
      real(rk), allocatable :: z(:)
      type(type_bottom_state_variable_id), allocatable :: id_m(:)
      type(type_horizontal_dependency_id) :: id_D1m, id_D2m, id_Dtot
      type(type_horizontal_diagnostic_variable_id) :: id_total
      type(type_horizontal_diagnostic_variable_id) :: id_layer(3)
   contains
      procedure :: do_bottom => summary_do_bottom
   end type

   type, extends(type_base_model) :: type_h2s_bed_transport
      integer  :: n = 0
      real(rk), allocatable :: z(:)
      real(rk) :: h_supply = 0.0_rk, minD = 0.0_rk, dtot = 0.0_rk
      type(type_bottom_state_variable_id), allocatable :: id_m(:)
      type(type_state_variable_id) :: id_H2S_pel
      type(type_horizontal_dependency_id) :: id_D1m, id_D2m, id_Dtot, id_poro, id_cmix
      type(type_horizontal_dependency_id) :: id_diff(3)
      type(type_horizontal_diagnostic_variable_id) :: id_J_req, id_J_app
      type(type_bottom_state_variable_id) :: id_cumJ, id_guard      ! accepted-step counters (docs/127 s10)
      logical  :: test = .false.
      real(rk) :: t_D1 = 0.0_rk, t_D2 = 0.0_rk, t_diff(3) = 0.0_rk, t_cpel = 0.0_rk
   contains
      procedure :: do_bottom => transport_do_bottom
   end type

   type, extends(type_base_model), public :: type_ersem_benthic_sulfur_cycle
      integer  :: isw_h2s_layers = 0   ! 0: legacy homogeneous G2_H2S column; 2: fixed grid (1, SR3-B, retired)
      integer  :: n_h2s = 0            ! fixed-grid cells
      real(rk), allocatable :: z_h2s(:)   ! fixed-grid edges (m), 0 = sediment surface
      real(rk) :: minD_h2s = 0.0_rk, dtot_h2s = 0.0_rk
      type(type_bottom_state_variable_id), allocatable :: id_h2s(:)
      type(type_bottom_state_variable_id) :: id_cumP, id_cumS, id_guard   ! accepted-step counters (docs/127 s10)
      logical  :: h2s_test = .false.                                      ! frozen-coefficient operator test
      real(rk) :: t_D1 = 0.0_rk, t_D2 = 0.0_rk, t_k(3) = 0.0_rk, t_P1 = 0.0_rk, t_P3 = 0.0_rk
      type(type_horizontal_diagnostic_variable_id) :: id_k_h2s(3), id_P_h2s_1, id_P_h2s_3
      ! State variable dependencies (layer-specific via benthic_column_dissolved_matter)
      type(type_bottom_state_variable_id) :: id_H2S_1, id_H2S_2, id_H2S_3
      type(type_bottom_state_variable_id) :: id_S0_1
      type(type_bottom_state_variable_id) :: id_S0s    ! layer-1 S0 as a SOLID (isw_S0_solid = 1)
      type(type_bottom_state_variable_id) :: id_S0_2   ! S0 in Layer 2 (H2S-NO3 product)
      type(type_bottom_state_variable_id) :: id_G2o  ! Oxygen in Layer 1
      type(type_bottom_state_variable_id) :: id_NO3_2  ! NO3 in Layer 2 for H2S-NO3 oxidation
      type(type_bottom_state_variable_id) :: id_G4n    ! Dinitrogen gas (N2 product)
      type(type_bottom_state_variable_id) :: id_K4n2   ! Ammonium in Layer 2 (DNRA product)

      ! Pelagic H2S at bottom for oxic barrier mechanism
      type(type_state_variable_id) :: id_H2S_pel
      type(type_state_variable_id) :: id_S0_pel
      type(type_state_variable_id) :: id_O2_pel

      ! Bottom cell thickness for pelagic exchange dimension conversion
      type(type_dependency_id) :: id_h_bottom

      ! Bottom PAR for light-driven interface oxidation (jsasaki 2026-08-15)
      type(type_dependency_id) :: id_par

      ! Alkalinity coupling
      type(type_bottom_state_variable_id) :: id_benTA   ! Alkalinity in Layer 1
      type(type_bottom_state_variable_id) :: id_benTA2  ! Alkalinity in Layer 2
      type(type_bottom_state_variable_id) :: id_benTA3  ! Alkalinity in Layer 3

      ! Layer depth dependencies
      type(type_horizontal_dependency_id) :: id_D1m, id_D2m, id_Dtot

      ! Organic matter remineralization rate from H2 bacteria
      type(type_horizontal_dependency_id) :: id_remin_rate

      ! Diagnostic variables
      type(type_horizontal_diagnostic_variable_id) :: id_R_sulfate_red
      type(type_horizontal_diagnostic_variable_id) :: id_R_H2S_ox_ben
      type(type_horizontal_diagnostic_variable_id) :: id_R_H2S_NO3_ox  ! H2S oxidation by NO3 in Layer 2
      type(type_horizontal_diagnostic_variable_id) :: id_R_S0_NO3_ox   ! S0 oxidation to SULFATE by NO3 in Layer 2
      type(type_horizontal_diagnostic_variable_id) :: id_R_S0_ox_ben
      type(type_horizontal_diagnostic_variable_id) :: id_R_S0_burial
      type(type_horizontal_diagnostic_variable_id) :: id_R_barrier_ox
      type(type_horizontal_diagnostic_variable_id) :: id_R_FeS_ben
      type(type_horizontal_diagnostic_variable_id) :: id_R_FeS_pel
      type(type_horizontal_diagnostic_variable_id) :: id_f_barrier

      ! ledger counters (isw_ledger = 1 only; nippon-steel docs/119 FC1)
      integer  :: isw_ledger
      type(type_horizontal_diagnostic_variable_id) :: id_ledger_srfes_H2S_3, id_ledger_no3fes_H2S_2
      type(type_horizontal_diagnostic_variable_id) :: id_ledger_thio_NO3_2, id_ledger_thio_S0_2
      type(type_horizontal_diagnostic_variable_id) :: id_ledger_thio_G4n, id_ledger_thio_K4n2, id_ledger_thio_benTA2
      type(type_horizontal_diagnostic_variable_id) :: id_ledger_srox_H2S_1, id_ledger_oxbur_S0s, id_ledger_oxbur_S0_1
      type(type_horizontal_diagnostic_variable_id) :: id_ledger_ox1_G2o, id_ledger_srox_benTA, id_ledger_sr3_benTA3
      type(type_horizontal_diagnostic_variable_id) :: id_ledger_barfes_H2S_pel, id_ledger_barrier_S0s
      type(type_horizontal_diagnostic_variable_id) :: id_ledger_barrier_S0_1, id_ledger_barrier_S0_pel
      type(type_horizontal_diagnostic_variable_id) :: id_ledger_barrier_O2_pel

      ! Parameters
      real(rk) :: K_H2S_prod     ! H2S production rate per unit remineralization (mol S/mol C)
      real(rk) :: K_H2S_ox       ! H2S oxidation rate constant (1/d)
      real(rk) :: f_ox_direct    ! fraction of interface H2S oxidation completed to SO4 in place (0 = legacy, all via S0)
      real(rk) :: K_H2S_NO3_ox   ! H2S oxidation by NO3 rate constant (1/d)
      real(rk) :: K_S0_NO3_ox    ! S0 -> SO4 oxidation by NO3 rate constant (1/d); 0 = off (legacy)
      real(rk) :: K_NO3_half     ! Half-saturation for NO3 (mmol/m3)
      real(rk) :: K_S0_ox        ! S0 oxidation rate constant (1/d)
      real(rk) :: K_S0_burial    ! S0 burial rate (1/d)
      real(rk) :: K_O2_half      ! Half-saturation for oxygen (mmol/m3)
      real(rk) :: K_barrier      ! Oxic barrier effectiveness (1/m)
      real(rk) :: K_barrier_rate ! Rate at which barrier oxidizes H2S (1/d)
      real(rk) :: K_FeS_ben      ! FeS precipitation rate in benthic layers (1/d)
      integer  :: isw_barrier_dest  ! 0: barrier S0 to the water (legacy), 1: to bed layer 1
      integer  :: isw_S0_solid      ! 0: layer-1 S0 in the dissolved column (legacy), 1: own solid pool
      real(rk) :: K_FeS_pel      ! FeS precipitation rate in pelagic (1/d)
      real(rk) :: p_sr_1         ! fraction of sulfate reduction delivered at the interface / layer 1 (jsasaki 2026-08-15, docs/14)
      real(rk) :: K_O2_half_pel  ! bottom-water O2 half-saturation (cubic Hill) for interface oxidation; 0 = legacy layer-1 Monod (jsasaki 2026-08-15, docs/14)
      real(rk) :: K_par_ox       ! PAR half-saturation for light-driven interface oxidation (mat photosynthesis O2); 0 = off (jsasaki 2026-08-15, docs/14)
      ! SR2 (nippon-steel docs/127 s7.2, s8 item 5): interface sulfate reduction excluded where an oxidant is available
      real(rk) :: O2_thr_sr      ! O2 threshold (bottom water AND layer-1 pore water, mmol O2/m3); 0 = factor off
      real(rk) :: NO3_thr_sr     ! layer-1 pore-water NO3 threshold (mmol N/m3); 0 = factor off
      real(rk) :: PAR_thr_sr     ! light threshold (W/m2); 0 = factor off
      logical  :: sr2_on         ! any threshold > 0
      type(type_bottom_state_variable_id) :: id_NO3_1
      type(type_horizontal_dependency_id) :: id_poro
      type(type_horizontal_diagnostic_variable_id) :: id_sr_share_1
      real(rk) :: f_DNRA         ! Fraction of H2S-NO3 N going to NH4 (0-1)

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
      ! COMPLETE interface oxidation (2026-09-11, nippon-steel docs/114 EP).
      ! The column's dissolved pools are spread UNIFORMLY over the column
      ! (benthic_column_dissolved_matter per_layer), so S0 made at the oxic
      ! interface is diluted at once to the whole column depth: with a 5 mm
      ! oxic layer over an 11 cm column only ~4 % of it stays where it can be
      ! oxidised, and the S0 -> SO4 step that returns the sulfate-reduction
      ! alkalinity barely runs (0.01-0.03 mmol S/m2/d against 2-3 produced).
      ! Coastal sediments reoxidise 70-92 % of their sulfide, mostly to
      ! sulfate at the interface. f_ox_direct completes that fraction of the
      ! layer-1 H2S oxidation in place: H2S + 2 O2 -> SO4(2-), -2 TA per S,
      ! no S0 released. 0 = legacy two-step chain, bit-identical.
      call self%get_parameter(self%f_ox_direct, 'f_ox_direct', '-', &
           'fraction of interface H2S oxidation completed to SO4 in place (0 = legacy)', &
           default=0.0_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%K_S0_ox, 'K_S0_ox', '1/d', &
           'S0 oxidation rate constant', default=0.02_rk)
      call self%get_parameter(self%K_S0_burial, 'K_S0_burial', '1/d', &
           'S0 burial rate (settling into deeper sediment)', default=0.5_rk)
      call self%get_parameter(self%K_O2_half, 'K_O2_half', 'mmol/m^3', &
           'half-saturation O2 for oxidation', default=1.0_rk)

      ! H2S-NO3 oxidation parameters (Layer 2 - denitrification zone)
      ! Chemolithotrophic sulfur oxidation coupled to denitrification:
      !   5 H2S + 2 NO3- -> 5 S0 + N2 + 4 H2O
      ! This is the key mechanism that prevents H2S from reaching Layer 1
      ! when a thick denitrification zone (Layer 2) exists.
      ! Typical rate: 10-100 1/d (fast reaction when both substrates present)
      call self%get_parameter(self%K_H2S_NO3_ox, 'K_H2S_NO3_ox', '1/d', &
           'H2S oxidation rate by NO3 in Layer 2', default=50.0_rk)
      ! SECOND STEP of thiodenitrification (jsasaki 2026-08-19, docs/18 §17).
      ! The H2S -> S0 step above is real and is what the large sulfur bacteria
      ! (Beggiatoa, Thioploca, Thiomargarita) do: they store the S0 in
      ! intracellular globules. What was missing is its FATE -- the same
      ! organisms oxidise the stored S0 on to SULFATE with nitrate when
      ! sulfide runs short. Without it S0_2 has production and no consumption
      ! at all (doc/sulfur_process_review.md item 7), and, because only
      ! oxidation BACK TO SULFATE regenerates the SO4^2- that sulfate
      ! reduction removed, the alkalinity layer 2 carries can never be given
      ! back. The model therefore retained 99.98 % of the sulfate-reduction
      ! alkalinity where coastal sediments reoxidise 70-92 % of their sulfide.
      ! 0 = off, so every run before this reproduces unchanged.
      call self%get_parameter(self%K_S0_NO3_ox, 'K_S0_NO3_ox', '1/d', &
           'S0 oxidation to sulfate by NO3 in Layer 2 (0 = off)', &
           default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%K_NO3_half, 'K_NO3_half', 'mmol/m^3', &
           'half-saturation NO3 for H2S-NO3 oxidation', default=10.0_rk)

      ! DNRA (Dissimilatory Nitrate Reduction to Ammonium) parameter
      ! Fraction of nitrogen from H2S-NO3 oxidation going to NH4 instead of N2.
      ! Default 0.0 preserves current behavior (all N to N2).
      call self%get_parameter(self%f_DNRA, 'f_DNRA', '-', &
           'fraction of H2S-NO3 nitrogen routed to NH4 (DNRA)', default=0.0_rk, &
           minimum=0.0_rk, maximum=1.0_rk)

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
      ! Where the oxic barrier's product goes (nippon-steel docs/114, 2026-09-12).
      ! The barrier oxidises sulfide AT THE SEDIMENT INTERFACE, but the legacy
      ! code puts the S0 it makes into the WATER (id_S0_pel), so elemental
      ! sulphur accumulates in the water column and its later oxidation charges
      ! the water's alkalinity. isw_barrier_dest = 1 puts the product in the
      ! bed's layer-1 S0 pool instead, where it is subject to the bed's own
      ! oxidation, burial and nitrate pathways. Default 0 = the former model.
      call self%get_parameter(self%isw_barrier_dest, 'isw_barrier_dest', '', &
           'barrier S0 destination (0: pelagic, 1: benthic layer 1)', default=0)

      ! Elemental sulphur is a SOLID (nippon-steel docs/114 CY, owner decision
      ! 2026-09-12). Layer-1 S0 has been a constituent of
      ! benthic_column_dissolved_matter, whose equilibrium stock is 36-55x the
      ! actual one for this species, so the bed returns to the water whatever is
      ! put in it -- measured. At isw_S0_solid = 1 layer-1 S0 lives in this
      ! module's own bottom state variable: same production, oxidation, burial
      ! and alkalinity terms, but it stays where it is made. The variable is
      ! registered only when the switch is on, so the default build is unchanged
      ! and gains no output variable.
      call self%get_parameter(self%isw_S0_solid, 'isw_S0_solid', '', &
           'layer-1 elemental sulphur (0: dissolved column, 1: own solid pool)', default=0)
      if (self%isw_S0_solid == 1) &
         call self%register_state_variable(self%id_S0s, 'S0s', 'mmol S/m^2', &
              'elemental sulfur, solid, layer 1', minimum=0.0_rk)
      ! K_FeS_pel: rate constant for FeS precipitation in pelagic bottom water
      ! Usually lower than benthic because less reactive Fe available in water column
      call self%get_parameter(self%K_FeS_pel, 'K_FeS_pel', '1/d', &
           'FeS precipitation rate in pelagic (scavenging)', default=0.1_rk)
      ! Vertical placement of sulfate reduction (jsasaki 2026-08-15; design:
      ! nippon-steel/docs/14-anaerobic-pathway.md section 4). The original
      ! wiring put ALL SR products (H2S + alkalinity) in Layer 3, whose
      ! exchange timescale is far longer than the tank's 72-h closures: the
      ! SR alkalinity accumulated at depth while the observed water-column
      ! dark TA rise (+79..+112 umol/kg/72 h) never appeared. In a thin,
      ! organic-rich sediment SR runs just below the interface, so a
      ! fraction p_sr_1 of the production (H2S and its +2 TA/S) is delivered
      ! to Layer 1, where the existing O2-dependent oxidation chain also
      ! yields the observed dark-light TA sign flip (lit: prompt reoxidation
      ! retracts the TA; dark: it survives). Default 0 = legacy layer-3-only
      ! behaviour (bit-identical).
      call self%get_parameter(self%p_sr_1, 'p_sr_1', '-', &
           'fraction of sulfate reduction delivered at the interface (Layer 1)', &
           default=0.0_rk, minimum=0.0_rk, maximum=1.0_rk)
      ! O2 control of the interface oxidation (jsasaki 2026-08-15, docs/14
      ! section 5). The legacy Monod uses the layer-1 CONCENTRATION
      ! G2o/D1m, which in the tank runs sits at a supersaturated
      ! 2900-3300 mmol/m^3 in lit AND dark alike (a pool/thickness
      ! bookkeeping artefact) - it carries no light/dark information, so
      ! sulfide oxidation cannot discriminate the regimes and the observed
      ! dark-light TA sign flip cannot be expressed. Physically the
      ! near-interface oxidation is supplied by O2 from the overlying
      ! water; with K_O2_half_pel > 0 the layer-1 H2S/S0 oxidation uses a
      ! cubic-Hill response to BOTTOM-WATER O2 instead (lit ~200-230 vs
      ! dark ~80-200 mmol/m^3 gives the contrast). Default 0 = legacy
      ! (bit-identical).
      call self%get_parameter(self%K_O2_half_pel, 'K_O2_half_pel', 'mmol/m^3', &
           'bottom-water O2 half-saturation (cubic Hill) for interface oxidation (0: legacy layer-1 Monod)', &
           default=0.0_rk, minimum=0.0_rk)
      ! Light-driven interface oxidation (jsasaki 2026-08-15, docs/14
      ! section 6). Even with the bottom-water O2 response, lit and dark
      ! regimes barely separate: the gated tank (like the real one) keeps
      ! dark bottom water at 110-200 mmol/m^3. The controlling physics is
      ! the MAT MICROENVIRONMENT: benthic-producer photosynthesis
      ! super-oxygenates the top millimetres in the light, while in
      ! darkness the mat goes anoxic within minutes regardless of the
      ! overlying water (classic microsensor observation). With
      ! K_par_ox > 0 the interface H2S/S0 oxidation factor becomes
      ! max(f_O2, PAR/(PAR + K_par_ox)) - light guarantees oxidation, and
      ! in the dark the O2 term takes over. Default 0 = off (bit-identical).
      call self%get_parameter(self%K_par_ox, 'K_par_ox', 'W/m^2', &
           'PAR half-saturation for light-driven interface oxidation (0: off)', &
           default=0.0_rk, minimum=0.0_rk)

      ! SR2 (jsasaki 2026-09-22; nippon-steel docs/127 s7.2 and s8 item 5, review rounds 38-39). The layer-1 share of
      ! sulfate reduction becomes p_sr_1 * product of linear exclusion ramps, each exactly zero at or above its
      ! threshold: bottom-water O2, layer-1 PORE-WATER O2 (G2o/(poro*D1m)), layer-1 PORE-WATER NO3 (NO3_1/(poro*D1m))
      ! and light. A threshold of 0 switches its factor off; all 0 = legacy (share = p_sr_1; bit-identical, no new
      ! coupling, no new output). The same share is used for H2S, its +2 TA/S and every ledger counter; the remainder
      ! goes to layer 3, so total sulfate reduction is unchanged.
      call self%get_parameter(self%O2_thr_sr, 'O2_thr_sr', 'mmol O_2/m^3', &
           'SR2: O2 availability threshold for interface sulfate reduction (0: off)', default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%NO3_thr_sr, 'NO3_thr_sr', 'mmol N/m^3', &
           'SR2: layer-1 pore-water NO3 availability threshold (0: off)', default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%PAR_thr_sr, 'PAR_thr_sr', 'W/m^2', &
           'SR2: light threshold excluding interface sulfate reduction (0: off)', default=0.0_rk, minimum=0.0_rk)
      self%sr2_on = self%O2_thr_sr > 0.0_rk .or. self%NO3_thr_sr > 0.0_rk .or. self%PAR_thr_sr > 0.0_rk
      if (self%NO3_thr_sr > 0.0_rk) call self%register_state_dependency(self%id_NO3_1, 'NO3_1', 'mmol N/m^2', &
           'nitrate in layer 1 (SR2 exclusion)')
      if (self%O2_thr_sr > 0.0_rk .or. self%NO3_thr_sr > 0.0_rk) &
           call self%register_dependency(self%id_poro, sediment_porosity)
      if (self%sr2_on) call self%register_diagnostic_variable(self%id_sr_share_1, 'sr_share_1', '-', &
           'SR2: share of sulfate reduction placed in layer 1', domain=domain_bottom, source=source_do_bottom)

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
      ! Counters whose SET call depends on isw_S0_solid, isw_barrier_dest or
      ! legacy_ersem_compatibility are registered under the same condition.
      call self%get_parameter(self%isw_ledger, 'isw_ledger', '', &
           'ledger counters: diagnostics of the source terms applied (0: off, 1: on)', &
           default=0, minimum=0, maximum=1)
      if (self%isw_ledger == 1) then
         call self%register_diagnostic_variable(self%id_ledger_srfes_H2S_3, 'ledger_srfes_H2S_3', 'mmol S/m^2/d', &
              'ledger: sulfate reduction (layer-3 share) and FeS precipitation -> H2S layer 3', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_no3fes_H2S_2, 'ledger_no3fes_H2S_2', 'mmol S/m^2/d', &
              'ledger: H2S oxidation by NO3 and FeS precipitation -> H2S layer 2', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_thio_NO3_2, 'ledger_thio_NO3_2', 'mmol N/m^2/d', &
              'ledger: thiodenitrification (H2S -> S0 and S0 -> SO4 by NO3) -> nitrate layer 2', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_thio_S0_2, 'ledger_thio_S0_2', 'mmol S/m^2/d', &
              'ledger: thiodenitrification (H2S -> S0 and S0 -> SO4 by NO3) -> S0 layer 2', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_thio_G4n, 'ledger_thio_G4n', 'mmol N/m^2/d', &
              'ledger: thiodenitrification (H2S -> S0 and S0 -> SO4 by NO3) -> dinitrogen gas', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_thio_K4n2, 'ledger_thio_K4n2', 'mmol N/m^2/d', &
              'ledger: thiodenitrification DNRA share -> ammonium layer 2', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_srox_H2S_1, 'ledger_srox_H2S_1', 'mmol S/m^2/d', &
              'ledger: sulfate reduction (layer-1 share), H2S oxidation and FeS precipitation -> H2S layer 1', &
              domain=domain_bottom, source=source_do_bottom)
         if (self%isw_S0_solid == 1) then
            call self%register_diagnostic_variable(self%id_ledger_oxbur_S0s, 'ledger_oxbur_S0s', 'mmol S/m^2/d', &
                 'ledger: H2S oxidation to S0, S0 oxidation and S0 burial -> solid S0 layer 1', &
                 domain=domain_bottom, source=source_do_bottom)
         else
            call self%register_diagnostic_variable(self%id_ledger_oxbur_S0_1, 'ledger_oxbur_S0_1', 'mmol S/m^2/d', &
                 'ledger: H2S oxidation to S0, S0 oxidation and S0 burial -> S0 layer 1', &
                 domain=domain_bottom, source=source_do_bottom)
         end if
         call self%register_diagnostic_variable(self%id_ledger_ox1_G2o, 'ledger_ox1_G2o', 'mmol O_2/m^2/d', &
              'ledger: H2S oxidation (to S0 and direct to SO4) and S0 oxidation -> benthic oxygen layer 1', &
              domain=domain_bottom, source=source_do_bottom)
         if (.not.legacy_ersem_compatibility) then
            call self%register_diagnostic_variable(self%id_ledger_thio_benTA2, 'ledger_thio_benTA2', 'mmol eq/m^2/d', &
                 'ledger: thiodenitrification (H2S -> S0 and S0 -> SO4 by NO3) -> benthic alkalinity layer 2', &
                 domain=domain_bottom, source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_srox_benTA, 'ledger_srox_benTA', 'mmol eq/m^2/d', &
                 'ledger: sulfate reduction (layer-1 share), S0 oxidation and direct H2S oxidation -> benthic alkalinity layer 1', &
                 domain=domain_bottom, source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_sr3_benTA3, 'ledger_sr3_benTA3', 'mmol eq/m^2/d', &
                 'ledger: sulfate reduction (layer-3 share) -> benthic alkalinity layer 3', &
                 domain=domain_bottom, source=source_do_bottom)
         end if
         call self%register_diagnostic_variable(self%id_ledger_barfes_H2S_pel, 'ledger_barfes_H2S_pel', 'mmol S/m^2/d', &
              'ledger: oxic barrier oxidation and pelagic FeS scavenging -> pelagic H2S (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         if (self%isw_barrier_dest == 1) then
            if (self%isw_S0_solid == 1) then
               call self%register_diagnostic_variable(self%id_ledger_barrier_S0s, 'ledger_barrier_S0s', 'mmol S/m^2/d', &
                    'ledger: oxic barrier oxidation -> solid S0 layer 1', &
                    domain=domain_bottom, source=source_do_bottom)
            else
               call self%register_diagnostic_variable(self%id_ledger_barrier_S0_1, 'ledger_barrier_S0_1', 'mmol S/m^2/d', &
                    'ledger: oxic barrier oxidation -> S0 layer 1', &
                    domain=domain_bottom, source=source_do_bottom)
            end if
         else
            call self%register_diagnostic_variable(self%id_ledger_barrier_S0_pel, 'ledger_barrier_S0_pel', 'mmol S/m^2/d', &
                 'ledger: oxic barrier oxidation -> pelagic S0 (exchange)', &
                 domain=domain_bottom, source=source_do_bottom)
         end if
         call self%register_diagnostic_variable(self%id_ledger_barrier_O2_pel, 'ledger_barrier_O2_pel', 'mmol O_2/m^2/d', &
              'ledger: oxic barrier oxidation -> pelagic oxygen (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
      end if

      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

      ! Register dependencies for layer-specific sulfur variables
      ! These link to variables created by benthic_column_dissolved_matter with composition 'h' and 'e'
      call self%get_parameter(self%isw_h2s_layers, 'isw_h2s_layers', '', &
           'bed sulfide: 0 homogeneous G2_H2S column (legacy), 2 fixed grid (1: retired SR3-B)', default=0, minimum=0, maximum=2)
      if (self%isw_h2s_layers == 1) call self%fatal_error('initialize', &
           'isw_h2s_layers = 1 (SR3-B moving sub-boxes) failed its transport gate and is retired (nippon-steel docs/127 s8.3); use 2')
      if (self%isw_h2s_layers == 2) then
         call register_h2s_grid(self)
      else
         call self%register_state_dependency(self%id_H2S_1, 'H2S_1', 'mmol S/m^2', &
              'hydrogen sulfide in layer 1')
         call self%register_state_dependency(self%id_H2S_2, 'H2S_2', 'mmol S/m^2', &
              'hydrogen sulfide in layer 2')
         call self%register_state_dependency(self%id_H2S_3, 'H2S_3', 'mmol S/m^2', &
              'hydrogen sulfide in layer 3')
      end if
      call self%register_state_dependency(self%id_S0_1, 'S0_1', 'mmol S/m^2', &
           'elemental sulfur in layer 1')
      call self%register_state_dependency(self%id_S0_2, 'S0_2', 'mmol S/m^2', &
           'elemental sulfur in layer 2')

      ! Oxygen in Layer 1 for oxidation reactions
      call self%register_state_dependency(self%id_G2o, 'G2o', 'mmol O_2/m^2', &
           'oxygen in layer 1')

      ! NO3 in Layer 2 for H2S-NO3 oxidation (chemolithotrophic denitrification)
      call self%register_state_dependency(self%id_NO3_2, 'NO3_2', 'mmol N/m^2', &
           'nitrate in layer 2')

      ! Dinitrogen gas - couples to benthic_nitrogen_cycle's G4n
      call self%register_state_dependency(self%id_G4n, 'G4n', 'mmol N/m^2', &
           'dinitrogen gas')

      ! Ammonium in Layer 2 for DNRA pathway
      call self%register_state_dependency(self%id_K4n2, 'K4n2', 'mmol N/m^2', &
           'ammonium in layer 2')

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

      ! Bottom cell thickness for converting volumetric rates to areal fluxes
      call self%register_dependency(self%id_h_bottom, standard_variables%cell_thickness)

      ! Alkalinity coupling for sulfur redox reactions
      ! Following benthic_nitrogen_cycle.F90 pattern
      if (.not.legacy_ersem_compatibility) then
         call self%register_state_dependency(self%id_benTA, 'benTA', 'mEq/m^2', &
              'benthic alkalinity in aerobic layer')
         call self%register_state_dependency(self%id_benTA2, 'benTA2', 'mEq/m^2', &
              'benthic alkalinity in anaerobic layer')
         call self%register_state_dependency(self%id_benTA3, 'benTA3', 'mEq/m^2', &
              'benthic alkalinity in layer 3')
      end if

      ! Organic matter remineralization rate from H2 bacteria
      ! This should be coupled to H2/fHG3c (benthic_bacteria respiration diagnostic).
      ! fHG3c is in mg C/m^2/d; conversion to mmol C via CMass is done in do_bottom.
      call self%register_dependency(self%id_remin_rate, 'remin_rate', 'mg C/m^2/d', &
           'anaerobic remineralization rate in layer 3')

      ! Diagnostic variables
      call self%register_diagnostic_variable(self%id_R_sulfate_red, 'R_sulfate_red', &
           'mmol S/m^2/d', 'sulfate reduction rate (H2S production)', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_H2S_ox_ben, 'R_H2S_ox_ben', &
           'mmol S/m^2/d', 'benthic H2S oxidation rate (O2)', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_H2S_NO3_ox, 'R_H2S_NO3_ox', &
           'mmol S/m^2/d', 'H2S oxidation by NO3 in Layer 2', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_S0_NO3_ox, 'R_S0_NO3_ox', &
           'mmol S/m^2/d', 'S0 oxidation to sulfate by NO3 in Layer 2', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_S0_ox_ben, 'R_S0_ox_ben', &
           'mmol S/m^2/d', 'benthic S0 oxidation rate', &
           domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_S0_burial, 'R_S0_burial', &
           'mmol S/m^2/d', 'benthic S0 burial rate', &
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

      real(rk) :: H2S_1, H2S_2, H2S_3, S0_1, S0_2, G2o, NO3_2
      real(rk) :: H2S_pel, S0_pel, O2_pel
      real(rk) :: D1m, D2m, remin_rate, h_bottom
      real(rk) :: O2_conc_1, NO3_conc_2, f_O2, f_O2_pel, f_NO3
      real(rk) :: R_sulfate_red, R_H2S_ox_1, R_H2S_NO3_ox, R_S0_NO3_ox, R_S0_ox_1, R_S0_burial
      real(rk) :: R_ox_direct, R_ox_to_S0
      real(rk) :: r_no3_S0, r_ta_S0
      real(rk) :: r_no3, r_ta
      real(rk) :: f_barrier, R_barrier_ox
      real(rk) :: par, f_par
      real(rk) :: share_1, poro, NO3_1, f_ex
      real(rk), allocatable :: mh(:), fr(:, :), ov(:, :)
      integer  :: kbox
      real(rk) :: Dtot, hl(3), kl(3), Pl1, Pl3, Sj, sumS
      real(rk) :: R_FeS_1, R_FeS_2, R_FeS_3, R_FeS_ben, R_FeS_pel

      _HORIZONTAL_LOOP_BEGIN_

         ! Get state variables
         if (self%isw_h2s_layers == 2) then
            if (.not. allocated(mh)) allocate(mh(self%n_h2s), fr(3, self%n_h2s), ov(3, self%n_h2s))
            _GET_HORIZONTAL_(self%id_D1m, D1m)
            _GET_HORIZONTAL_(self%id_D2m, D2m)
            _GET_HORIZONTAL_(self%id_Dtot, Dtot)
            call h2s_check_geometry(self, D1m, D2m, Dtot, self%minD_h2s, self%dtot_h2s)
            hl = (/ D1m, D2m - D1m, Dtot - D2m /)
            call h2s_coverage(self%z_h2s, D1m, D2m, Dtot, ov)
            do kbox = 1, self%n_h2s
               _GET_HORIZONTAL_(self%id_h2s(kbox), mh(kbox))
               fr(:, kbox) = ov(:, kbox) / (self%z_h2s(kbox + 1) - self%z_h2s(kbox))
            end do
            ! the ERSEM-layer inventories the kinetics act on (each cell weighted by its length in the layer), clipped
            ! per CELL so that the cell sinks below and the layer rates driving the products are the same sums
            H2S_1 = sum(fr(1, :) * max(0.0_rk, mh))
            H2S_2 = sum(fr(2, :) * max(0.0_rk, mh))
            H2S_3 = sum(fr(3, :) * max(0.0_rk, mh))
         else
            _GET_HORIZONTAL_(self%id_H2S_1, H2S_1)
            _GET_HORIZONTAL_(self%id_H2S_2, H2S_2)
            _GET_HORIZONTAL_(self%id_H2S_3, H2S_3)
         end if
         _GET_HORIZONTAL_(self%id_S0_2, S0_2)   ! only WRITTEN before K_S0_NO3_ox
         if (self%isw_S0_solid == 1) then
            _GET_HORIZONTAL_(self%id_S0s, S0_1)
         else
            _GET_HORIZONTAL_(self%id_S0_1, S0_1)
         end if
         _GET_HORIZONTAL_(self%id_G2o, G2o)
         _GET_HORIZONTAL_(self%id_NO3_2, NO3_2)

         ! Get pelagic variables at bottom
         _GET_(self%id_H2S_pel, H2S_pel)
         _GET_(self%id_S0_pel, S0_pel)
         _GET_(self%id_O2_pel, O2_pel)

         ! Get layer depths
         _GET_HORIZONTAL_(self%id_D1m, D1m)
         _GET_HORIZONTAL_(self%id_D2m, D2m)

         ! Get remineralization rate from H2 bacteria
         _GET_HORIZONTAL_(self%id_remin_rate, remin_rate)

         ! Get bottom cell thickness for pelagic exchange conversion (m)
         _GET_(self%id_h_bottom, h_bottom)

         ! Ensure non-negative
         H2S_1 = max(0.0_rk, H2S_1)
         H2S_2 = max(0.0_rk, H2S_2)
         S0_2 = max(0.0_rk, S0_2)
         H2S_3 = max(0.0_rk, H2S_3)
         S0_1 = max(0.0_rk, S0_1)
         G2o = max(0.0_rk, G2o)
         NO3_2 = max(0.0_rk, NO3_2)
         H2S_pel = max(0.0_rk, H2S_pel)
         S0_pel = max(0.0_rk, S0_pel)
         O2_pel = max(0.0_rk, O2_pel)
         remin_rate = max(0.0_rk, remin_rate)

         ! ============================================================
         ! LAYER 3: Sulfate reduction (always active - Layer 3 is anoxic)
         ! ============================================================
         ! H2S production is proportional to anaerobic OM remineralization
         ! Stoichiometry: 53 SO4 per 106 C -> 0.5 mol S per mol C
         ! remin_rate is in mg C/m^2/d (from H2/fHG3c); convert to mmol C via CMass
         R_sulfate_red = self%K_H2S_prod * remin_rate / CMass

         ! SR2: the layer-1 share (legacy: the constant p_sr_1)
         share_1 = self%p_sr_1
         if (self%sr2_on) then
            f_ex = 1.0_rk
            if (self%O2_thr_sr > 0.0_rk .or. self%NO3_thr_sr > 0.0_rk) then
               _GET_HORIZONTAL_(self%id_poro, poro)
               if (.not. (poro > 0.0_rk) .or. .not. (D1m > 0.0_rk)) &
                  call self%fatal_error('do_bottom', 'SR2: porosity and D1m must be positive')
            end if
            if (self%O2_thr_sr > 0.0_rk) f_ex = f_ex * max(0.0_rk, 1.0_rk - max(0.0_rk, O2_pel) / self%O2_thr_sr) &
                                                    * max(0.0_rk, 1.0_rk - max(0.0_rk, G2o) / (poro * D1m) / self%O2_thr_sr)
            if (self%NO3_thr_sr > 0.0_rk) then
               _GET_HORIZONTAL_(self%id_NO3_1, NO3_1)
               f_ex = f_ex * max(0.0_rk, 1.0_rk - max(0.0_rk, NO3_1) / (poro * D1m) / self%NO3_thr_sr)
            end if
            if (self%PAR_thr_sr > 0.0_rk) then
               _GET_(self%id_par, par)
               f_ex = f_ex * max(0.0_rk, 1.0_rk - max(0.0_rk, par) / self%PAR_thr_sr)
            end if
            share_1 = self%p_sr_1 * min(1.0_rk, max(0.0_rk, f_ex))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sr_share_1, share_1)
         end if

         ! ============================================================
         ! LAYER 1: H2S and S0 oxidation (limited by O2 availability)
         ! ============================================================
         ! Convert depth-integrated O2 to concentration
         ! Guard against division by very small D1m (layer collapse)
         O2_conc_1 = G2o / max(D1m, 0.0001_rk)

         ! Oxygen limitation: cubic-Hill on bottom-water O2 when
         ! K_O2_half_pel > 0 (docs/14 section 5), else legacy Monod on the
         ! layer-1 concentration.
         if (self%K_O2_half_pel > 0.0_rk) then
            f_O2 = max(0.0_rk, O2_pel)**3 &
                 / (max(0.0_rk, O2_pel)**3 + self%K_O2_half_pel**3)
         else
            f_O2 = O2_conc_1 / (O2_conc_1 + self%K_O2_half)
         end if
         if (self%K_par_ox > 0.0_rk) then
            ! Mat-photosynthesis oxygenation: light guarantees interface
            ! oxidation regardless of the water-column O2 (docs/14 sec. 6)
            _GET_(self%id_par, par)
            f_par = max(0.0_rk, par) / (max(0.0_rk, par) + self%K_par_ox)
            f_O2 = max(f_O2, f_par)
         end if

         ! H2S + 0.5 O2 -> S0
         R_H2S_ox_1 = self%K_H2S_ox * H2S_1 * f_O2
         R_ox_direct = self%f_ox_direct * R_H2S_ox_1      ! completed to SO4 in place
         R_ox_to_S0  = R_H2S_ox_1 - R_ox_direct            ! legacy route via S0

         ! S0 + 1.5 O2 -> SO4 (removed from system)
         R_S0_ox_1 = self%K_S0_ox * S0_1 * f_O2

         ! S0 burial (settling into deeper sediment, irreversible removal)
         ! Elemental sulfur settles/buries into anoxic layers where it may
         ! undergo disproportionation or further reactions
         R_S0_burial = self%K_S0_burial * S0_1

         ! ============================================================
         ! LAYER 2: H2S oxidation by NO3 (chemolithotrophic nitrate reduction)
         ! ============================================================
         ! H2S diffusing from Layer 3 through Layer 2 is oxidized by nitrate.
         ! This is the key mechanism that prevents H2S from reaching Layer 1
         ! (and the pelagic) when a thick suboxic zone exists.
         !
         ! N product is partitioned by f_DNRA between N2 and NH4.
         ! r_NO3 = 2/(5+3*f_DNRA) adjusts for the different electron demands.
         !
         ! The presence of NO3 in Layer 2 indicates H2S cannot coexist.
         !
         ! Convert depth-integrated NO3 to concentration
         ! Layer 2 thickness = D2m - D1m
         NO3_conc_2 = NO3_2 / max(D2m - D1m, 0.0001_rk)

         ! NO3 limitation (Michaelis-Menten)
         f_NO3 = NO3_conc_2 / (NO3_conc_2 + self%K_NO3_half)

         ! H2S oxidation by NO3 in Layer 2
         ! This is a fast reaction when both substrates are present
         R_H2S_NO3_ox = self%K_H2S_NO3_ox * H2S_2 * f_NO3

         ! S0 oxidation to SULFATE by NO3, the second step (see K_S0_NO3_ox)
         R_S0_NO3_ox = self%K_S0_NO3_ox * S0_2 * f_NO3

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
         R_barrier_ox = self%K_barrier_rate * H2S_pel * f_barrier * h_bottom

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

         ! FeS precipitation in benthic Layer 2 (suboxic layer)
         ! Intermediate rate - some Fe(II) available from reduction
         R_FeS_2 = self%K_FeS_ben * 0.5_rk * H2S_2

         ! FeS precipitation in benthic Layer 3 (anoxic layer)
         ! Higher rate because Fe(II) is abundant from Fe(III) reduction
         R_FeS_3 = self%K_FeS_ben * H2S_3

         ! Total benthic FeS precipitation
         R_FeS_ben = R_FeS_1 + R_FeS_2 + R_FeS_3

         ! FeS precipitation (scavenging) in pelagic bottom water
         ! Removes H2S through reaction with particulate Fe or settling FeS
         R_FeS_pel = self%K_FeS_pel * H2S_pel * h_bottom

         ! ============================================================
         ! Set ODEs
         ! ============================================================
         ! Layer 3: H2S production from sulfate reduction, loss from FeS burial.
         ! A p_sr_1 fraction of the production is delivered at the interface
         ! (Layer 1) instead - see the p_sr_1 parameter note (docs/14).
         if (self%isw_h2s_layers == 2) then
            ! fixed grid: sulfate reduction placed by each cell's length in layers 1 and 3; the layers' first-order sink
            ! constants weighted by the cell's fractions (their cell sums are R_H2S_ox_1 + R_FeS_1, R_H2S_NO3_ox + R_FeS_2
            ! and R_FeS_3 above)
            kl = (/ self%K_H2S_ox * f_O2 + self%K_FeS_ben * 0.1_rk, self%K_H2S_NO3_ox * f_NO3 + self%K_FeS_ben * 0.5_rk, &
                    self%K_FeS_ben /)
            Pl1 = share_1 * R_sulfate_red
            Pl3 = (1.0_rk - share_1) * R_sulfate_red
            if (self%h2s_test) then          ! frozen-coefficient operator test: the cell ODE only
               kl = self%t_k
               Pl1 = self%t_P1
               Pl3 = self%t_P3
               hl = (/ self%t_D1, self%t_D2 - self%t_D1, Dtot - self%t_D2 /)
               call h2s_coverage(self%z_h2s, self%t_D1, self%t_D2, Dtot, ov)
               do kbox = 1, self%n_h2s
                  fr(:, kbox) = ov(:, kbox) / (self%z_h2s(kbox + 1) - self%z_h2s(kbox))
               end do
            end if
            sumS = 0.0_rk
            do kbox = 1, self%n_h2s
               Sj = sum(kl * fr(:, kbox)) * max(0.0_rk, mh(kbox))
               sumS = sumS + Sj
               _SET_BOTTOM_ODE_(self%id_h2s(kbox), Pl1 * ov(1, kbox) / hl(1) + Pl3 * ov(3, kbox) / hl(3) - Sj)
            end do
            _SET_BOTTOM_ODE_(self%id_cumP, Pl1 + Pl3)
            _SET_BOTTOM_ODE_(self%id_cumS, sumS)
            _SET_BOTTOM_ODE_(self%id_guard, real(count(mh < 0.0_rk), rk))
            ! first-order constants and production, for the per-closure steady initialisation (docs/127 s9)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_k_h2s(1), kl(1))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_k_h2s(2), kl(2))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_k_h2s(3), kl(3))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_P_h2s_1, share_1 * R_sulfate_red)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_P_h2s_3, (1.0_rk - share_1) * R_sulfate_red)
         else
            _SET_BOTTOM_ODE_(self%id_H2S_3, (1.0_rk - share_1) * R_sulfate_red - R_FeS_3)
         end if
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_srfes_H2S_3, (1.0_rk - share_1) * R_sulfate_red - R_FeS_3)
         end if

         ! Layer 2: H2S consumption by NO3 oxidation and FeS precipitation
         ! Electron-balanced stoichiometry for H2S-NO3 coupling:
         !   H2S -> S0 + 2e- (2 electrons per H2S)
         !   Denitrification: NO3 + 5e- -> 0.5 N2 (5 e-/mol NO3)
         !   DNRA:            NO3 + 8e- -> NH4     (8 e-/mol NO3)
         !   Mixed: (5+3*f_DNRA) e-/mol NO3 -> r_NO3 = 2/(5+3*f_DNRA)
         r_no3 = 2.0_rk / (5.0_rk + 3.0_rk * self%f_DNRA)

         ! S0 -> SO4 needs 6 e-/mol S against H2S -> S0's 2, so the same
         ! electron balance gives exactly THREE times the nitrate per sulfur:
         !   N2 route:   5 S0 + 6 NO3- + 2 H2O -> 5 SO4^2- + 3 N2 + 4 H+
         !   DNRA route: 4 S0 + 3 NO3- + 7 H2O -> 4 SO4^2- + 3 NH4+ + 2 H+
         r_no3_S0 = 6.0_rk / (5.0_rk + 3.0_rk * self%f_DNRA)

         if (self%isw_h2s_layers == 0) _SET_BOTTOM_ODE_(self%id_H2S_2, -R_H2S_NO3_ox - R_FeS_2)
         _SET_BOTTOM_ODE_(self%id_NO3_2, -r_no3 * R_H2S_NO3_ox &
                                         - r_no3_S0 * R_S0_NO3_ox)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_no3fes_H2S_2, -R_H2S_NO3_ox - R_FeS_2)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_thio_NO3_2, -r_no3 * R_H2S_NO3_ox &
                                         - r_no3_S0 * R_S0_NO3_ox)
         end if
         ! S0 produced 1:1 with H2S consumed, and consumed by the second step
         _SET_BOTTOM_ODE_(self%id_S0_2,   R_H2S_NO3_ox - R_S0_NO3_ox)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_thio_S0_2, R_H2S_NO3_ox - R_S0_NO3_ox)
         end if

         ! Partition N between N2 (denitrification) and NH4 (DNRA)
         _SET_BOTTOM_ODE_(self%id_G4n,  (1.0_rk - self%f_DNRA) &
                                        * (r_no3 * R_H2S_NO3_ox &
                                           + r_no3_S0 * R_S0_NO3_ox))
         _SET_BOTTOM_ODE_(self%id_K4n2, self%f_DNRA &
                                        * (r_no3 * R_H2S_NO3_ox &
                                           + r_no3_S0 * R_S0_NO3_ox))
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_thio_G4n, (1.0_rk - self%f_DNRA) &
                                        * (r_no3 * R_H2S_NO3_ox &
                                           + r_no3_S0 * R_S0_NO3_ox))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_thio_K4n2, self%f_DNRA &
                                        * (r_no3 * R_H2S_NO3_ox &
                                           + r_no3_S0 * R_S0_NO3_ox))
         end if

         ! Alkalinity: denitrification +1 TA/mol NO3, DNRA +2 TA/mol NO3
         ! r_TA = (1+f_DNRA) * r_NO3 = 2*(1+f_DNRA)/(5+3*f_DNRA)
         r_ta = (1.0_rk + self%f_DNRA) * r_no3
         ! The second step ALSO regenerates the SO4^2- that sulfate reduction
         ! removed, which is -2 eq per S and is the whole point of it: at
         ! f_DNRA = 0 the net is -2 + 1.2 = -0.8 eq per S0 (a RETRACTION),
         ! and at f_DNRA = 1 it is -2 + 2 x 0.75 = -0.5.
         r_ta_S0 = -2.0_rk + (1.0_rk + self%f_DNRA) * r_no3_S0
         if (.not.legacy_ersem_compatibility) &
            _SET_BOTTOM_ODE_(self%id_benTA2, r_ta * R_H2S_NO3_ox &
                                             + r_ta_S0 * R_S0_NO3_ox)
         if (self%isw_ledger == 1 .and. .not.legacy_ersem_compatibility) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_thio_benTA2, r_ta * R_H2S_NO3_ox &
                                             + r_ta_S0 * R_S0_NO3_ox)
         end if

         ! Layer 1: H2S delivery from interface sulfate reduction (p_sr_1),
         !          consumption by oxidation and FeS precipitation,
         !          S0 production from H2S oxidation, loss from oxidation and burial
         if (self%isw_h2s_layers == 0) _SET_BOTTOM_ODE_(self%id_H2S_1, share_1 * R_sulfate_red - R_H2S_ox_1 - R_FeS_1)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_srox_H2S_1, share_1 * R_sulfate_red - R_H2S_ox_1 - R_FeS_1)
         end if
         if (self%isw_S0_solid == 1) then
            _SET_BOTTOM_ODE_(self%id_S0s,  R_ox_to_S0 - R_S0_ox_1 - R_S0_burial)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_oxbur_S0s, R_ox_to_S0 - R_S0_ox_1 - R_S0_burial)
            end if
         else
            _SET_BOTTOM_ODE_(self%id_S0_1, R_ox_to_S0 - R_S0_ox_1 - R_S0_burial)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_oxbur_S0_1, R_ox_to_S0 - R_S0_ox_1 - R_S0_burial)
            end if
         end if
         _SET_BOTTOM_ODE_(self%id_G2o,   -0.5_rk * R_ox_to_S0 - 2.0_rk * R_ox_direct - 1.5_rk * R_S0_ox_1)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_ox1_G2o, -0.5_rk * R_ox_to_S0 - 2.0_rk * R_ox_direct - 1.5_rk * R_S0_ox_1)
         end if

         ! Alkalinity, layer 1: interface sulfate reduction +2 TA per mol H2S;
         ! S0 + 1.5 O2 + H2O -> SO4^2- + 2H+ => -2 TA per mol S0
         if (.not.legacy_ersem_compatibility) &
            _SET_BOTTOM_ODE_(self%id_benTA, 2.0_rk * share_1 * R_sulfate_red - 2.0_rk * R_S0_ox_1 - 2.0_rk * R_ox_direct)
         if (self%isw_ledger == 1 .and. .not.legacy_ersem_compatibility) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_srox_benTA, 2.0_rk * share_1 * R_sulfate_red - 2.0_rk * R_S0_ox_1 - 2.0_rk * R_ox_direct)
         end if

         ! Layer 3: sulfate reduction produces +2 TA per mol H2S
         ! SO4^2- + 2C_org -> H2S + 2HCO3- (net +2 mEq per mol H2S)
         if (.not.legacy_ersem_compatibility) &
            _SET_BOTTOM_ODE_(self%id_benTA3, 2.0_rk * (1.0_rk - share_1) * R_sulfate_red)
         if (self%isw_ledger == 1 .and. .not.legacy_ersem_compatibility) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_sr3_benTA3, 2.0_rk * (1.0_rk - share_1) * R_sulfate_red)
         end if

         ! Pelagic: H2S removal by oxic barrier oxidation and FeS scavenging
         ! Barrier oxidation produces S0, FeS scavenging is irreversible removal
         ! Note: Using _SET_BOTTOM_EXCHANGE_ applies flux to bottom cell only
         _SET_BOTTOM_EXCHANGE_(self%id_H2S_pel, -R_barrier_ox - R_FeS_pel)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_barfes_H2S_pel, -R_barrier_ox - R_FeS_pel)
         end if
         if (self%isw_barrier_dest == 1) then
            ! the interface keeps the S0 it makes
            if (self%isw_S0_solid == 1) then
               _SET_BOTTOM_ODE_(self%id_S0s,  R_barrier_ox)
               if (self%isw_ledger == 1) then
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_barrier_S0s, R_barrier_ox)
               end if
            else
               _SET_BOTTOM_ODE_(self%id_S0_1, R_barrier_ox)
               if (self%isw_ledger == 1) then
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_barrier_S0_1, R_barrier_ox)
               end if
            end if
         else
            _SET_BOTTOM_EXCHANGE_(self%id_S0_pel, R_barrier_ox)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_barrier_S0_pel, R_barrier_ox)
            end if
         end if
         _SET_BOTTOM_EXCHANGE_(self%id_O2_pel,  -0.5_rk * R_barrier_ox)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_barrier_O2_pel, -0.5_rk * R_barrier_ox)
         end if

         ! Set diagnostics
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_sulfate_red, R_sulfate_red)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_H2S_ox_ben, R_H2S_ox_1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_H2S_NO3_ox, R_H2S_NO3_ox)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_S0_NO3_ox, R_S0_NO3_ox)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_S0_ox_ben, R_S0_ox_1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_S0_burial, R_S0_burial)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_barrier_ox, R_barrier_ox)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_FeS_ben, R_FeS_ben)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_FeS_pel, R_FeS_pel)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_barrier, f_barrier)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

   subroutine register_h2s_grid(self)
      class(type_ersem_benthic_sulfur_cycle), intent(inout), target :: self
      class(type_h2s_bed_summary),   pointer :: summ
      class(type_h2s_bed_transport), pointer :: tran
      integer :: ilay, k, n
      character(len=16) :: nm
      real(rk) :: h_supply, top, cap, growth

      call self%get_parameter(top, 'h2s_grid_top', 'm', 'fixed grid: width of the top cell', default=1.0e-3_rk, &
           minimum=1.0e-5_rk)
      call self%get_parameter(cap, 'h2s_grid_cap', 'm', 'fixed grid: largest cell width', default=5.0e-3_rk, &
           minimum=1.0e-5_rk)
      call self%get_parameter(growth, 'h2s_grid_growth', '-', 'fixed grid: width ratio of successive cells up to the cap', &
           default=1.25_rk, minimum=1.0_rk)
      call self%get_parameter(self%dtot_h2s, 'h2s_grid_dtot', 'm', &
           'fixed grid: sediment column depth (must equal the column''s depth_of_sediment_column)', default=0.0_rk)
      call self%get_parameter(self%minD_h2s, 'minD_h2s', 'm', &
           'fixed grid: smallest admissible ERSEM layer thickness (fatal below; > 0)', default=1.0e-4_rk, minimum=0.0_rk)
      if (.not. (self%minD_h2s > 0.0_rk)) call self%fatal_error('register_h2s_grid', 'minD_h2s must be > 0')
      call self%get_parameter(h_supply, 'h_supply_h2s', 'mmol S/m^3', &
           'uptake limiter c_pel/(c_pel + h) on water-to-bed sulfide exchange (0: off)', default=0.0_rk, minimum=0.0_rk)
      if (.not. (self%dtot_h2s > 0.0_rk)) call self%fatal_error('register_h2s_grid', 'h2s_grid_dtot must be set (> 0)')
      if (cap < top) call self%fatal_error('register_h2s_grid', 'h2s_grid_cap must be >= h2s_grid_top')
      call h2s_mesh(top, cap, growth, self%dtot_h2s, self%z_h2s)
      n = size(self%z_h2s) - 1
      self%n_h2s = n
      allocate(self%id_h2s(n))
      do k = 1, n
         write(nm, '(a,i0)') 'h2s_c', k
         call self%register_state_variable(self%id_h2s(k), trim(nm), 'mmol S/m^2', &
              'bed sulfide, fixed-grid cell (pore water)', 0.0_rk, minimum=0.0_rk)
      end do
      do ilay = 1, 3
         write(nm, '(a,i1)') 'k_h2s_', ilay
         call self%register_diagnostic_variable(self%id_k_h2s(ilay), trim(nm), '1/d', &
              'first-order sulfide sink constant of the ERSEM layer', domain=domain_bottom, source=source_do_bottom)
      end do
      call self%register_diagnostic_variable(self%id_P_h2s_1, 'P_h2s_1', 'mmol S/m^2/d', &
           'sulfate reduction placed in layer 1', domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_P_h2s_3, 'P_h2s_3', 'mmol S/m^2/d', &
           'sulfate reduction placed in layer 3', domain=domain_bottom, source=source_do_bottom)
      ! accepted-step counters: integrated by the driver with the same stages as the cells, restored with them on a
      ! rejected attempt, so the cell balance and any guard activation are exact records (docs/127 s10, round 41)
      call self%register_state_variable(self%id_cumP, 'h2s_cum_P', 'mmol S/m^2', &
           'fixed grid: cumulative sulfate reduction placed in the cells', 0.0_rk)
      call self%register_state_variable(self%id_cumS, 'h2s_cum_S', 'mmol S/m^2', &
           'fixed grid: cumulative first-order removal from the cells', 0.0_rk)
      call self%register_state_variable(self%id_guard, 'h2s_guard', 'd', &
           'fixed grid: cell-days with a negative cell at a stage evaluation (0 = the max(0, m) guard never acted)', 0.0_rk)
      ! frozen-coefficient operator test (default off): the cell ODE uses these constants instead of the model's
      call self%get_parameter(self%h2s_test, 'h2s_test', '', 'fixed grid: frozen-coefficient operator test', default=.false.)
      if (self%h2s_test) then
         call self%get_parameter(self%t_D1, 'h2s_test_D1', 'm', 'test: depth of interface 1')
         call self%get_parameter(self%t_D2, 'h2s_test_D2', 'm', 'test: depth of interface 2')
         call self%get_parameter(self%t_k(1), 'h2s_test_k1', '1/d', 'test: sink constant, layer 1')
         call self%get_parameter(self%t_k(2), 'h2s_test_k2', '1/d', 'test: sink constant, layer 2')
         call self%get_parameter(self%t_k(3), 'h2s_test_k3', '1/d', 'test: sink constant, layer 3')
         call self%get_parameter(self%t_P1, 'h2s_test_P1', 'mmol S/m^2/d', 'test: production in layer 1')
         call self%get_parameter(self%t_P3, 'h2s_test_P3', 'mmol S/m^2/d', 'test: production in layer 3')
      end if

      allocate(summ)
      summ%dt = 86400._rk
      summ%n = n
      summ%z = self%z_h2s
      call self%add_child(summ, 'h2s_summary', configunit=-1)
      allocate(summ%id_m(n))
      do k = 1, n
         write(nm, '(a,i0)') 'm', k
         call summ%register_state_dependency(summ%id_m(k), trim(nm), 'mmol S/m^2', 'bed sulfide cell')
      end do
      call summ%register_dependency(summ%id_D1m, depth_of_bottom_interface_of_layer_1)
      call summ%register_dependency(summ%id_D2m, depth_of_bottom_interface_of_layer_2)
      call summ%register_dependency(summ%id_Dtot, depth_of_sediment_column)
      call summ%register_diagnostic_variable(summ%id_total, 'total', 'mmol S/m^2', 'bed sulfide, total', &
           domain=domain_bottom, source=source_do_bottom)
      do ilay = 1, 3
         write(nm, '(a,i1)') 'layer', ilay
         call summ%register_diagnostic_variable(summ%id_layer(ilay), trim(nm), 'mmol S/m^2', 'bed sulfide, ERSEM layer', &
              domain=domain_bottom, source=source_do_bottom)
      end do

      allocate(tran)
      tran%dt = 86400._rk
      tran%n = n
      tran%z = self%z_h2s
      tran%h_supply = h_supply
      tran%minD = self%minD_h2s
      tran%dtot = self%dtot_h2s
      tran%test = self%h2s_test
      if (self%h2s_test) then
         tran%t_D1 = self%t_D1
         tran%t_D2 = self%t_D2
         call self%get_parameter(tran%t_diff(1), 'h2s_test_diff1', 'm^2/d', 'test: diffusivity, layer 1')
         call self%get_parameter(tran%t_diff(2), 'h2s_test_diff2', 'm^2/d', 'test: diffusivity, layer 2')
         call self%get_parameter(tran%t_diff(3), 'h2s_test_diff3', 'm^2/d', 'test: diffusivity, layer 3')
         call self%get_parameter(tran%t_cpel, 'h2s_test_cpel', 'mmol S/m^3', 'test: water sulfide seen by the exchange')
      end if
      call self%add_child(tran, 'h2s_transport', configunit=-1)
      allocate(tran%id_m(n))
      do k = 1, n
         write(nm, '(a,i0)') 'm', k
         call tran%register_state_dependency(tran%id_m(k), trim(nm), 'mmol S/m^2', 'bed sulfide cell')
      end do
      call tran%register_state_dependency(tran%id_H2S_pel, 'H2S_pel', 'mmol S/m^3', 'pelagic hydrogen sulfide')
      call tran%register_dependency(tran%id_D1m, depth_of_bottom_interface_of_layer_1)
      call tran%register_dependency(tran%id_D2m, depth_of_bottom_interface_of_layer_2)
      call tran%register_dependency(tran%id_Dtot, depth_of_sediment_column)
      call tran%register_dependency(tran%id_poro, sediment_porosity)
      call tran%register_dependency(tran%id_cmix, pelagic_benthic_transfer_constant)
      call tran%register_dependency(tran%id_diff(1), diffusivity_in_sediment_layer_1)
      call tran%register_dependency(tran%id_diff(2), diffusivity_in_sediment_layer_2)
      call tran%register_dependency(tran%id_diff(3), diffusivity_in_sediment_layer_3)
      call tran%register_diagnostic_variable(tran%id_J_req, 'J_requested', 'mmol S/m^2/d', &
           'bed-to-water sulfide exchange before the uptake limiter', domain=domain_bottom, source=source_do_bottom)
      call tran%register_diagnostic_variable(tran%id_J_app, 'J_applied', 'mmol S/m^2/d', &
           'bed-to-water sulfide exchange applied', domain=domain_bottom, source=source_do_bottom)
      call tran%register_state_variable(tran%id_cumJ, 'cum_J', 'mmol S/m^2', &
           'cumulative bed-to-water sulfide exchange applied (accepted steps)', 0.0_rk)
      call tran%register_state_variable(tran%id_guard, 'guard', 'd', &
           'cell-days with a negative cell or negative water sulfide at a stage evaluation', 0.0_rk)
      do k = 1, n
         write(nm, '(a,i0)') 'h2s_c', k
         call summ%request_coupling(summ%id_m(k), '../'//trim(nm))
         call tran%request_coupling(tran%id_m(k), '../'//trim(nm))
      end do
      call tran%request_coupling(tran%id_H2S_pel, '../H2S_pel')
   end subroutine register_h2s_grid

   subroutine h2s_mesh(top, cap, growth, dtot, z)
      ! the frozen mesh rule of tools/sr3_grid_reference.py mesh(): widths top, top*growth, ... while below the cap
      ! and inside the column, then the remaining depth in equal cells no wider than the cap
      real(rk), intent(in) :: top, cap, growth, dtot
      real(rk), allocatable, intent(out) :: z(:)
      real(rk) :: w, e, step
      integer :: ng, nr, i
      ng = 0
      w = top
      e = 0.0_rk
      do while (w < cap .and. e + w < dtot)
         e = e + w
         w = w * growth
         ng = ng + 1
      end do
      nr = ceiling((dtot - e) / cap - 1.0e-12_rk)
      allocate(z(ng + nr + 1))
      z(1) = 0.0_rk
      w = top
      do i = 1, ng
         z(i + 1) = z(i) + w
         w = w * growth
      end do
      step = (dtot - z(ng + 1)) / nr
      do i = 1, nr - 1
         z(ng + 1 + i) = z(ng + 1) + i * step
      end do
      z(ng + nr + 1) = dtot
   end subroutine h2s_mesh

   pure real(rk) function h2s_overlap(a, b, lo, hi)
      real(rk), intent(in) :: a, b, lo, hi
      h2s_overlap = max(0.0_rk, min(b, hi) - max(a, lo))
   end function h2s_overlap

   pure subroutine h2s_coverage(z, D1m, D2m, Dtot, ov)
      ! ov(l, j): length of cell j inside ERSEM layer l
      real(rk), intent(in) :: z(:), D1m, D2m, Dtot
      real(rk), intent(out) :: ov(:, :)
      real(rk) :: cuts(4)
      integer :: j, l
      cuts = (/ 0.0_rk, D1m, D2m, Dtot /)
      do j = 1, size(z) - 1
         do l = 1, 3
            ov(l, j) = h2s_overlap(z(j), z(j + 1), cuts(l), cuts(l + 1))
         end do
      end do
   end subroutine h2s_coverage

   subroutine h2s_check_geometry(model, D1m, D2m, Dtot, minD, dtot_grid)
      ! fail fast (docs/127 s8 item 6, s9): ordered interfaces, every ERSEM layer at least minD thick, and the grid
      ! built for this column depth
      class(type_base_model), intent(in) :: model
      real(rk), intent(in) :: D1m, D2m, Dtot, minD, dtot_grid
      if (abs(Dtot - dtot_grid) > 1.0e-9_rk * dtot_grid) &
         call model%fatal_error('h2s_check_geometry', 'fixed grid: h2s_grid_dtot differs from depth_of_sediment_column')
      if (.not. (D1m >= minD .and. D2m - D1m >= minD .and. Dtot - D2m >= minD)) &
         call model%fatal_error('h2s_check_geometry', 'fixed grid: an ERSEM layer is thinner than minD_h2s or unordered')
   end subroutine h2s_check_geometry

   subroutine summary_do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_h2s_bed_summary), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      integer :: k
      real(rk) :: m, tot, lay(3), D1m, D2m, Dtot, ov(3, self%n)

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_D1m, D1m)
         _GET_HORIZONTAL_(self%id_D2m, D2m)
         _GET_HORIZONTAL_(self%id_Dtot, Dtot)
         call h2s_coverage(self%z, D1m, D2m, Dtot, ov)
         tot = 0.0_rk
         lay = 0.0_rk
         do k = 1, self%n
            _GET_HORIZONTAL_(self%id_m(k), m)
            lay = lay + ov(:, k) / (self%z(k + 1) - self%z(k)) * m
            tot = tot + m
         end do
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_total, tot)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layer(1), lay(1))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layer(2), lay(2))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_layer(3), lay(3))
      _HORIZONTAL_LOOP_END_
   end subroutine summary_do_bottom

   subroutine transport_do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_h2s_bed_transport), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      integer :: k, l
      real(rk) :: D1m, D2m, Dtot, poro, cmix, H2S_pel, diff(3), cuts(4)
      real(rk) :: m(self%n), c(self%n), r(self%n), zc(self%n)
      real(rk) :: F, J, Jreq, cp, Rc, Rs
      integer  :: nneg

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_D1m, D1m)
         _GET_HORIZONTAL_(self%id_D2m, D2m)
         _GET_HORIZONTAL_(self%id_Dtot, Dtot)
         _GET_HORIZONTAL_(self%id_poro, poro)
         _GET_HORIZONTAL_(self%id_cmix, cmix)
         _GET_HORIZONTAL_(self%id_diff(1), diff(1))
         _GET_HORIZONTAL_(self%id_diff(2), diff(2))
         _GET_HORIZONTAL_(self%id_diff(3), diff(3))
         _GET_(self%id_H2S_pel, H2S_pel)
         nneg = 0
         if (H2S_pel < 0.0_rk) nneg = 1
         if (self%test) then                 ! frozen-coefficient operator test
            D1m = self%t_D1
            D2m = self%t_D2
            diff = self%t_diff
            H2S_pel = self%t_cpel
         end if
         call h2s_check_geometry(self, D1m, D2m, Dtot, self%minD, self%dtot)
         if (.not. (poro > 0.0_rk .and. minval(diff) > 0.0_rk .and. cmix >= 0.0_rk)) &
            call self%fatal_error('transport_do_bottom', 'fixed grid: porosity, diffusivities or cmix invalid')
         cuts = (/ 0.0_rk, D1m, D2m, Dtot /)
         do k = 1, self%n
            _GET_HORIZONTAL_(self%id_m(k), m(k))
            if (m(k) < 0.0_rk) nneg = nneg + 1
            c(k) = max(0.0_rk, m(k)) / (poro * (self%z(k + 1) - self%z(k)))
            zc(k) = 0.5_rk * (self%z(k) + self%z(k + 1))
         end do
         r = 0.0_rk
         do k = 1, self%n - 1
            Rc = 0.0_rk             ! resistance between the two centres, int dz/diff over the layers it crosses
            do l = 1, 3
               Rc = Rc + h2s_overlap(zc(k), zc(k + 1), cuts(l), cuts(l + 1)) / diff(l)
            end do
            F = (c(k + 1) - c(k)) / Rc                    ! upward flux from cell k+1 into cell k
            r(k) = r(k) + F
            r(k + 1) = r(k + 1) - F
         end do
         Rs = cmix
         do l = 1, 3
            Rs = Rs + h2s_overlap(0.0_rk, zc(1), cuts(l), cuts(l + 1)) / diff(l)
         end do
         cp = max(0.0_rk, H2S_pel)
         Jreq = (c(1) - cp) / Rs
         J = Jreq
         if (J < 0.0_rk .and. self%h_supply > 0.0_rk) J = J * cp / (cp + self%h_supply)
         r(1) = r(1) - J
         do k = 1, self%n
            _SET_BOTTOM_ODE_(self%id_m(k), r(k))
         end do
         _SET_BOTTOM_EXCHANGE_(self%id_H2S_pel, J)
         _SET_BOTTOM_ODE_(self%id_cumJ, J)
         _SET_BOTTOM_ODE_(self%id_guard, real(nneg, rk))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_J_req, Jreq)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_J_app, J)
      _HORIZONTAL_LOOP_END_
   end subroutine transport_do_bottom

end module ersem_benthic_sulfur_cycle
